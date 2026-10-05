import json
import tempfile
import threading
import unittest
from pathlib import Path
from unittest.mock import patch

from backend.services.job_manager import JobManager


class JobRecoveryTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.patch = patch('backend.services.job_manager.JOB_DIR', self.root)
        self.patch.start()
        self.manager = JobManager()

    def tearDown(self):
        self.manager.shutdown()
        self.patch.stop()
        self.tmp.cleanup()

    def record(self, kind, job_id, status):
        path = self.root / kind / job_id / 'status.json'
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps({'job_id': job_id, 'status': status, 'created_at': 'original'}))
        return path

    def test_recovery_all_job_types_preserves_terminal_results(self):
        for kind in ('primer', 'search', 'similar', 'similar_export', 'similar_pairwise', 'phylogeny'):
            for state in ('queued', 'running', 'completed', 'failed'):
                path = self.record(kind, state, state)
                (path.parent / 'results.json').write_text('original results')
        self.manager.start()
        for path in self.root.glob('*/*/status.json'):
            data = json.loads(path.read_text())
            self.assertEqual(data['created_at'], 'original')
            self.assertEqual((path.parent / 'results.json').read_text(), 'original results')
            if path.parent.name in ('queued', 'running'):
                self.assertEqual(data['status'], 'failed')
                self.assertEqual(data['failure_reason'], 'server_restart')
            else:
                self.assertEqual(data, {'job_id': path.parent.name, 'status': path.parent.name, 'created_at': 'original'})

    def test_second_owner_cannot_change_live_jobs(self):
        self.manager.start()
        path = self.record('primer', 'live', 'running')
        second = JobManager()
        try:
            with self.assertRaisesRegex(RuntimeError, 'another API process'):
                second.start()
            self.assertEqual(json.loads(path.read_text())['status'], 'running')
        finally:
            second.shutdown()

    def test_shutdown_drains_queue_and_releases_ownership(self):
        self.manager.start()
        gate = threading.Event()
        entered = threading.Event()
        finished = []

        def work(number):
            entered.set()
            gate.wait(5)
            finished.append(number)

        for number in range(3):
            self.manager.submit('primer', str(number), work, number)
        self.assertTrue(entered.wait(2))
        stopped = threading.Event()
        thread = threading.Thread(target=lambda: (self.manager.shutdown(), stopped.set()))
        thread.start()
        self.assertFalse(stopped.wait(.05))
        gate.set()
        thread.join(5)
        self.assertTrue(stopped.is_set())
        self.assertEqual(sorted(finished), [0, 1, 2])
        self.manager.start()

    def test_unhandled_worker_exception_marks_failed(self):
        self.manager.start()
        path = self.record('primer', 'broken', 'queued')

        def fail():
            raise ValueError('outside service try block')

        self.manager.submit('primer', 'broken', fail)
        self.manager.shutdown()
        self.assertEqual(json.loads(path.read_text())['failure_reason'], 'worker_exception')

    def test_corrupt_status_does_not_prevent_other_recovery(self):
        bad = self.record('search', 'bad', 'running')
        bad.write_text('{')
        good = self.record('search', 'good', 'queued')
        self.manager.start()
        self.assertEqual(bad.read_text(), '{')
        self.assertEqual(json.loads(good.read_text())['status'], 'failed')

    def test_recovery_is_idempotent(self):
        path = self.record('primer', 'old', 'running')
        self.manager.start()
        first = path.read_bytes()
        self.manager.shutdown()
        self.manager.start()
        self.assertEqual(first, path.read_bytes())


if __name__ == '__main__':
    unittest.main()
