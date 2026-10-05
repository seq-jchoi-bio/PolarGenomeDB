from __future__ import annotations

import datetime as dt
import fcntl
import json
import logging
import os
import shutil
import threading
import uuid
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Callable

from ..api_limits import get_api_limits
from ..config import JOB_DIR


JobCallable = Callable[..., None]
logger = logging.getLogger(__name__)
INTERRUPTED_MESSAGE = "This job was interrupted by a server restart. Please submit the analysis again."


class JobManager:
    def __init__(self) -> None:
        limits = get_api_limits().jobs
        self.max_workers = limits.max_concurrent_jobs
        self.retention_days = limits.job_retention_days
        self._executor = ThreadPoolExecutor(max_workers=self.max_workers, thread_name_prefix="apdb-job")
        self._cleanup_lock = threading.Lock()
        self._cleanup_done = False
        self._owner_lock = None
        self._accepting = False
        self._closed = False

    def start(self) -> None:
        if self._owner_lock is not None:
            raise RuntimeError("Job manager is already started")
        JOB_DIR.mkdir(parents=True, exist_ok=True)
        lock = (JOB_DIR / ".worker.lock").open("a+")
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            lock.close()
            raise RuntimeError("Job directory is owned by another API process. Use one worker per job directory.") from None
        self._owner_lock = lock
        try:
            # The exclusive directory lock makes all old active records orphaned.
            for path in JOB_DIR.glob("*/*/status.json"):
                self._fail_active(path, INTERRUPTED_MESSAGE, "server_restart")
            if self._closed:
                self._executor = ThreadPoolExecutor(max_workers=self.max_workers, thread_name_prefix="apdb-job")
                self._closed = False
            self._accepting = True
        except BaseException:
            self._release_lock()
            raise

    def _release_lock(self) -> None:
        if self._owner_lock is not None:
            fcntl.flock(self._owner_lock.fileno(), fcntl.LOCK_UN)
            self._owner_lock.close()
            self._owner_lock = None

    def shutdown(self) -> None:
        self._accepting = False
        # Keep ownership until every accepted task has stopped writing results.
        self._executor.shutdown(wait=True)
        self._closed = True
        self._release_lock()

    def _fail_active(self, path: Path, message: str, reason: str) -> None:
        try:
            status = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            logger.warning("Cannot read job status during recovery: %s", path)
            return
        if not isinstance(status, dict) or status.get("status") not in {"queued", "running"}:
            return
        now = dt.datetime.now(dt.timezone.utc).isoformat()
        status.update(status="failed", error=message, user_error=message,
                      failure_reason=reason, failed_at=now, finished_at=now, updated_at=now)
        temporary = path.with_name(f".status-{uuid.uuid4().hex}.tmp")
        try:
            temporary.write_text(json.dumps(status, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
            os.replace(temporary, path)
        finally:
            temporary.unlink(missing_ok=True)
        logger.warning("Marked interrupted job failed: %s (%s)", path.parent.name, reason)

    def _run(self, job_type: str, job_id: str, target: JobCallable, args: tuple) -> None:
        try:
            target(*args)
        except Exception:
            logger.exception("Unhandled worker exception: %s/%s", job_type, job_id)
            self._fail_active(JOB_DIR / job_type / job_id / "status.json",
                              "The analysis stopped unexpectedly. Please submit the analysis again.", "worker_exception")

    def submit(self, job_type: str, job_id: str, target: JobCallable, *args: object) -> None:
        if not self._accepting:
            raise RuntimeError("Job manager is not accepting jobs")
        self.cleanup_once()
        self._executor.submit(self._run, job_type, job_id, target, args)

    def cleanup_once(self) -> None:
        with self._cleanup_lock:
            if self._cleanup_done:
                return
            self.cleanup_expired_jobs()
            self._cleanup_done = True

    def cleanup_expired_jobs(self) -> int:
        if self.retention_days <= 0 or not JOB_DIR.exists():
            return 0
        cutoff = dt.datetime.now(dt.timezone.utc) - dt.timedelta(days=self.retention_days)
        removed = 0
        for status_path in JOB_DIR.glob("*/*/status.json"):
            try:
                job_dir = status_path.parent
                mtime = dt.datetime.fromtimestamp(status_path.stat().st_mtime, tz=dt.timezone.utc)
                if mtime >= cutoff:
                    continue
                # Never remove a live or queued job while a worker may use it.
                try:
                    if json.loads(status_path.read_text()).get("status") in {"queued", "running"}:
                        continue
                except (ValueError, AttributeError):
                    continue
                shutil.rmtree(job_dir)
                removed += 1
            except OSError:
                continue
        return removed

    def status(self) -> dict[str, int]:
        return {
            "max_concurrent_jobs": self.max_workers,
            "job_retention_days": self.retention_days,
            "job_count": count_job_dirs(JOB_DIR),
        }


def count_job_dirs(root: Path) -> int:
    if not root.exists():
        return 0
    return sum(1 for path in root.glob("*/*") if path.is_dir())


job_manager = JobManager()
