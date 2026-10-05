from __future__ import annotations

import csv
import datetime as dt
import json
import re
import subprocess
import uuid
from io import StringIO
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional

from pydantic import BaseModel, Field
from Bio import SeqIO

from ..api_limits import get_api_limits
from ..config import JOB_DIR, ORTHOLOG_BIN, ROOT
from ..db import connect
from ..services.gene_service import annotate_search_hit
from ..runtime import load_runtime_json
from .job_files import write_json


SearchProgram = Literal["blastn", "blastp", "diamond"]

SEARCH_LIMITS = get_api_limits().search
BLAST_DISPLAY_EVALUE = 1e-5


class SearchRequest(BaseModel):
    species_id: str = Field(min_length=1)
    program: SearchProgram
    query: str = Field(min_length=1)
    max_hits: int = Field(default=20, ge=1, le=SEARCH_LIMITS.max_hits)


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def job_root() -> Path:
    path = JOB_DIR / "search"
    path.mkdir(parents=True, exist_ok=True)
    return path


def job_path(job_id: str) -> Path:
    if not re.fullmatch(r"[a-f0-9-]{36}", job_id):
        raise ValueError("invalid job_id")
    return job_root() / job_id


def read_json(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def max_sequences_for_program(program: SearchProgram) -> int:
    if program == "diamond":
        return SEARCH_LIMITS.max_diamond_sequences_per_request
    return SEARCH_LIMITS.max_blast_sequences_per_request


def max_hits_for_program(program: SearchProgram) -> int:
    return 1 if program == "diamond" else 10


def query_sequence_count(query: str) -> int:
    normalized = query.strip()
    if not normalized:
        return 0
    if not normalized.startswith(">"):
        return 1
    return sum(1 for line in normalized.splitlines() if line.startswith(">"))


def normalize_query(query: str, *, max_sequences: int, program: SearchProgram) -> str:
    query = query.strip()
    if not query:
        raise ValueError("Enter a sequence or upload a FASTA file.")
    if len(query) > SEARCH_LIMITS.max_query_characters:
        raise ValueError(f"query is too long; max {SEARCH_LIMITS.max_query_characters} characters")
    if not query.startswith(">"):
        if ">" in query:
            raise ValueError("FASTA input must begin with a > header before the first sequence.")
        query = f">query\n{query}"
    if query_sequence_count(query) > max_sequences:
        raise ValueError(f"too many FASTA query sequences; max {max_sequences}")

    records = list(SeqIO.parse(StringIO(query), "fasta"))
    seen = set()
    normalized = []
    for index, record in enumerate(records, 1):
        if not record.id:
            raise ValueError(f"FASTA record {index} needs a name after >.")
        if record.id in seen:
            raise ValueError(f"FASTA record {index} repeats a sequence ID. Use a unique ID for each query.")
        seen.add(record.id)
        sequence = re.sub(r"\s+", "", str(record.seq)).upper()
        if not sequence:
            raise ValueError(f"FASTA record {index} has no sequence.")
        alphabet = "ACGTRYSWKMBDHVN" if program == "blastn" else "ACDEFGHIKLMNPQRSTVWYBXZJUO*"
        if set(sequence) - set(alphabet):
            kind = "DNA letters (A, C, G, T or IUPAC ambiguity codes)" if program == "blastn" else "one-letter amino-acid codes (including X and *)"
            raise ValueError(f"Sequence {index} must contain only {kind}. Remove numbers, gaps and other symbols.")
        normalized.append(f">{record.description}\n{sequence}\n")
    return "".join(normalized)


def search_failure_message(error: Exception) -> str:
    if isinstance(error, subprocess.TimeoutExpired):
        return "Search exceeded the time limit. Try fewer or shorter sequences."
    if isinstance(error, FileNotFoundError):
        return "The search service is unavailable. Please contact the platform administrator."
    return "Search could not be completed. Please retry; if it persists, contact the platform administrator."


def asset_for_program(species_id: str, program: SearchProgram) -> Dict[str, Any]:
    assets = load_runtime_json(f"{species_id}/search_assets.json")["assets"]
    asset_type = {
        "blastn": "blast_genome",
        "blastp": "blast_protein",
        "diamond": "diamond_protein",
    }[program]
    for asset in assets:
        if asset["asset_type"] == asset_type and asset["status"] == "ready":
            return asset
    raise ValueError(f"search asset not available: {species_id}/{asset_type}")


def resolve_project_path(path: str) -> str:
    value = Path(path)
    if value.is_absolute():
        return str(value)
    return str((ROOT / value).resolve())


def create_search_job(request: SearchRequest) -> Dict[str, Any]:
    normalized_query = normalize_query(
        request.query,
        max_sequences=max_sequences_for_program(request.program),
        program=request.program,
    )
    asset = asset_for_program(request.species_id, request.program)
    job_id = str(uuid.uuid4())
    path = job_path(job_id)
    path.mkdir(parents=True, exist_ok=False)

    query_path = path / "query.fa"
    result_path = path / "raw_results.tsv"
    log_path = path / "run.log"
    status_path = path / "status.json"
    request_path = path / "request.json"

    query_path.write_text(normalized_query, encoding="utf-8")
    write_json(request_path, request.model_dump())
    write_json(
        status_path,
        {
            "job_id": job_id,
            "status": "queued",
            "created_at": utc_now(),
            "updated_at": utc_now(),
            "species_id": request.species_id,
            "program": request.program,
            "result_count": 0,
        },
    )
    return {
        "job_id": job_id,
        "status": "queued",
        "species_id": request.species_id,
        "program": request.program,
        "created_at": read_json(status_path)["created_at"],
        "asset": asset,
    }


def execute_search_job(job_id: str) -> None:
    path = job_path(job_id)
    request = SearchRequest(**read_json(path / "request.json"))
    asset = asset_for_program(request.species_id, request.program)
    run_search_job(
        job_id,
        asset,
        path / "query.fa",
        path / "raw_results.tsv",
        path / "run.log",
        path / "status.json",
        request,
    )


def command_for_request(
    request: SearchRequest,
    asset: Dict[str, Any],
    query_path: Path,
    result_path: Path,
    *,
    raw_export: bool = False,
) -> List[str]:
    outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    max_hits = min(request.max_hits, max_hits_for_program(request.program))
    if request.program == "diamond":
        command = [
            str(ORTHOLOG_BIN / "diamond"),
            "blastp",
            "--query",
            str(query_path),
            "--db",
            resolve_project_path(asset["index_prefix"]),
            "--out",
            str(result_path),
            "--outfmt",
            "6",
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "mismatch",
            "gapopen",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore",
            "--threads",
            "2",
        ]
        return command if raw_export else [*command, "--max-target-seqs", str(max_hits)]
    binary = str(ORTHOLOG_BIN / request.program)
    command = [
        binary,
        "-query",
        str(query_path),
        "-db",
        resolve_project_path(asset["index_prefix"]),
        "-out",
        str(result_path),
        "-outfmt",
        outfmt,
    ]
    if raw_export:
        # Use BLAST defaults and retain every returned target. The UI applies its
        # stricter cutoff only after this raw file has been stored.
        return command
    return [*command, "-max_target_seqs", str(max_hits)]


def run_search_job(
    job_id: str,
    asset: Dict[str, Any],
    query_path: Path,
    result_path: Path,
    log_path: Path,
    status_path: Path,
    request: SearchRequest,
) -> None:
    started_at = utc_now()
    write_json(
        status_path,
        {
            "job_id": job_id,
            "status": "running",
            "created_at": read_json(status_path)["created_at"],
            "started_at": started_at,
            "updated_at": started_at,
            "species_id": request.species_id,
            "program": request.program,
            "result_count": 0,
        },
    )
    cmd = command_for_request(request, asset, query_path, result_path, raw_export=True)
    try:
        with log_path.open("w", encoding="utf-8") as log:
            subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT, timeout=SEARCH_LIMITS.timeout_seconds)
        hits = parse_tabular_results(result_path)
        finished_at = utc_now()
        write_json(
            status_path,
            {
                "job_id": job_id,
                "status": "completed",
                "created_at": read_json(status_path)["created_at"],
                "started_at": started_at,
                "finished_at": finished_at,
                "updated_at": finished_at,
                "species_id": request.species_id,
                "program": request.program,
                "result_count": len(hits),
            },
        )
    except Exception as exc:
        failed_at = utc_now()
        write_json(
            status_path,
            {
                "job_id": job_id,
                "status": "failed",
                "created_at": read_json(status_path)["created_at"],
                "started_at": started_at,
                "failed_at": failed_at,
                "updated_at": failed_at,
                "species_id": request.species_id,
                "program": request.program,
                "result_count": 0,
                "error": str(exc),
                "user_error": search_failure_message(exc),
            },
        )


def parse_tabular_results(path: Path) -> List[Dict[str, Any]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    fields = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    hits: List[Dict[str, Any]] = []
    with path.open("rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, fieldnames=fields, delimiter="\t")
        for row in reader:
            for key in ["pident", "evalue", "bitscore"]:
                row[key] = float(row[key])
            for key in ["length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send"]:
                row[key] = int(row[key])
            hits.append(row)
    return hits


def get_search_job(job_id: str) -> Dict[str, Any]:
    path = job_path(job_id)
    status_path = path / "status.json"
    if not status_path.exists():
        raise FileNotFoundError(job_id)
    status = read_json(status_path)
    raw_result_path = path / "raw_results.tsv"
    # Jobs created before raw export used results.tsv directly.
    result_path = raw_result_path if raw_result_path.exists() else path / "results.tsv"
    hits = parse_tabular_results(result_path) if status["status"] == "completed" else []
    if status.get("program") in {"blastn", "blastp"}:
        grouped: Dict[str, List[Dict[str, Any]]] = {}
        for hit in hits:
            if float(hit["evalue"]) > BLAST_DISPLAY_EVALUE:
                continue
            query_hits = grouped.setdefault(str(hit["qseqid"]), [])
            if len(query_hits) < max_hits_for_program(status["program"]):
                query_hits.append(hit)
        hits = [hit for query_hits in grouped.values() for hit in query_hits]
    elif status.get("program") == "diamond":
        grouped: Dict[str, List[Dict[str, Any]]] = {}
        for hit in hits:
            query_hits = grouped.setdefault(str(hit["qseqid"]), [])
            if len(query_hits) < max_hits_for_program("diamond"):
                query_hits.append(hit)
        hits = [hit for query_hits in grouped.values() for hit in query_hits]
    if hits:
        with connect() as conn:
            hits = [annotate_search_hit(conn, status["species_id"], status["program"], hit) for hit in hits]
    status["results"] = hits
    return status


def raw_search_results(job_id: str, query_id: Optional[str] = None) -> List[str]:
    path = job_path(job_id)
    status = read_json(path / "status.json")
    if status.get("status") != "completed":
        raise ValueError("search job has not completed")
    result_path = path / "raw_results.tsv"
    if not result_path.exists():
        raise FileNotFoundError(job_id)
    lines: List[str] = []
    with result_path.open("rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if query_id is None or line.split("\t", 1)[0] == query_id:
                lines.append(line)
    return lines
