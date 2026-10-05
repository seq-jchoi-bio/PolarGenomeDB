from __future__ import annotations

import datetime as dt
import json
import re
import subprocess
import uuid
from pathlib import Path
from typing import Any, Dict, Literal

from pydantic import BaseModel, Field

from ..config import JOB_DIR, ORTHOLOG_BIN
from ..db import connect, row, rows
from .job_files import write_json
from .protein_similarity import protein_similarity
from .similar_jobs import get_similar_job, normalized_protein_ids, selected_reference_proteins


TreeModel = Literal["LG+G", "WAG+G", "JTT+G", "Auto"]
TreeAnalysis = Literal["isoform", "reference", "similar"]


class PhylogenyRequest(BaseModel):
    species_id: str = Field(min_length=1)
    gene_id: str = Field(min_length=1)
    model: TreeModel = "LG+G"
    analysis: TreeAnalysis = "reference"
    protein_id: str | None = None
    similar_job_id: str | None = None
    reference_protein_ids: list[str] | None = None


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def job_root() -> Path:
    path = JOB_DIR / "phylogeny"
    path.mkdir(parents=True, exist_ok=True)
    return path


def job_path(job_id: str) -> Path:
    if not re.fullmatch(r"[a-f0-9-]{36}", job_id):
        raise ValueError("invalid job_id")
    return job_root() / job_id


def read_json(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def clean_protein_sequence(value: object) -> str:
    return re.sub(r"[^A-Za-z]", "", str(value or "")).upper().rstrip("*")


def prepare_tree_input(request: PhylogenyRequest) -> Dict[str, Any]:
    if request.analysis == "similar":
        return prepare_similar_tree_input(request)
    with connect() as conn:
        query_proteins = rows(conn, "SELECT protein_id, transcript_id, sequence FROM proteins WHERE species_id = ? AND gene_id = ? ORDER BY transcript_id, protein_id", (request.species_id, request.gene_id))
    valid_proteins = [protein for protein in query_proteins if len(clean_protein_sequence(protein["sequence"])) >= 20]
    if request.analysis == "reference":
        selected = next((protein for protein in valid_proteins if protein["protein_id"] == request.protein_id), None)
        if selected is None:
            selected = max(valid_proteins, key=lambda protein: len(clean_protein_sequence(protein["sequence"])), default=None)
        valid_proteins = [selected] if selected else []
    sequences = []
    all_hits = []
    for protein in valid_proteins:
        sequence = clean_protein_sequence(protein["sequence"])
        if len(sequence) < 20:
            continue
        query_id = f"Q{len(sequences) + 1}"
        sequences.append({"id": query_id, "label": protein["transcript_id"], "kind": "query", "sequence": sequence})
        if request.analysis == "reference":
            all_hits.extend(protein_similarity(request.species_id, request.gene_id, protein["protein_id"], 5)["hits"])
    if not sequences:
        raise ValueError("no query protein sequence is available for this gene")
    if request.analysis == "isoform":
        if len(sequences) < 2:
            raise ValueError("at least two translated isoforms are required for isoform comparison")
        return {"species_id": request.species_id, "gene_id": request.gene_id, "model": request.model, "analysis": request.analysis, "query_protein_count": len(sequences), "sequences": sequences, "reference_hits": []}
    # Reference mode is transcript-specific: use only the selected query
    # transcript and its own five DIAMOND hits as tree inputs.
    hits = sorted(all_hits, key=lambda hit: (int(hit["rank"]), -float(hit["bitscore"] or 0), float(hit["evalue"] or 1)))[:5]
    if not hits:
        raise ValueError("no reference similarity hits are available for the selected transcript")
    with connect() as conn:
        for index, hit in enumerate(hits, start=1):
            reference = row(conn, "SELECT sequence FROM reference_proteins WHERE reference_id = ? AND protein_id = ?", (hit["reference_id"], hit["reference_protein_id"]))
            if reference and len(clean_protein_sequence(reference["sequence"])) >= 20:
                hit["reference_sequence"] = reference["sequence"]
    reference_hits = []
    for index, hit in enumerate(hits, start=1):
        sequence = clean_protein_sequence(hit.get("reference_sequence"))
        if len(sequence) < 20:
            continue
        sequence_id = f"R{index}"
        label = str(hit.get("reference_transcript_id") or hit["reference_protein_id"] or hit["reference_gene_id"])
        sequences.append({"id": sequence_id, "label": label, "kind": "reference", "sequence": sequence})
        reference_hits.append(
            {
                "id": sequence_id,
                "label": label,
                "reference_id": hit["reference_id"],
                "reference_display_name": hit["reference_display_name"],
                "reference_protein_id": hit["reference_protein_id"],
                "rank": hit["rank"],
                "pident": hit["pident"],
                "bitscore": hit["bitscore"],
                "evalue": hit["evalue"],
            }
        )

    if len(sequences) < 3:
        raise ValueError("at least two reference protein sequences are required for a tree")
    return {
        "species_id": request.species_id,
        "gene_id": request.gene_id,
        "model": request.model,
        "analysis": request.analysis,
        "query_protein_id": selected["protein_id"] if selected else None,
        "query_transcript_id": selected["transcript_id"] if selected else None,
        "query_protein_count": len(sequences) - len(reference_hits),
        "sequences": sequences,
        "reference_hits": reference_hits,
    }


def prepare_similar_tree_input(request: PhylogenyRequest) -> Dict[str, Any]:
    if not request.similar_job_id:
        raise ValueError("similar_job_id is required for a similar-gene tree")
    protein_ids = normalized_protein_ids(request.reference_protein_ids or [])
    if len(protein_ids) > 8:
        raise ValueError("a maximum of 8 input isoforms can be included in one phylogenetic tree")
    similar_job = get_similar_job(request.similar_job_id)
    if similar_job.get("status") != "completed" or not similar_job.get("results"):
        raise ValueError("the similar-gene result must be completed before building a tree")
    if similar_job.get("species_id") != request.species_id:
        raise ValueError("similar-gene job species does not match the requested species")

    result = similar_job["results"]
    groups = result.get("groups") if isinstance(result, dict) else []
    groups = groups if isinstance(groups, list) else []
    available_ids = {
        str(group.get("query_protein", {}).get("protein_id"))
        for group in groups
        if isinstance(group, dict) and isinstance(group.get("query_protein"), dict)
    }
    invalid = [protein_id for protein_id in protein_ids if protein_id not in available_ids]
    if invalid:
        raise ValueError("tree isoforms must originate from the selected similar-gene result")

    reference, references = selected_reference_proteins(request.species_id, protein_ids)
    candidate_ids: list[str] = []
    for group in groups:
        if not isinstance(group, dict):
            continue
        query = group.get("query_protein")
        if not isinstance(query, dict) or str(query.get("protein_id")) not in protein_ids:
            continue
        candidates = group.get("candidates")
        if not isinstance(candidates, list):
            continue
        candidate_ids.extend(str(item.get("protein_id")) for item in candidates[:5] if isinstance(item, dict) and item.get("protein_id"))

    target_by_id: dict[str, dict[str, Any]] = {}
    unique_candidate_ids = list(dict.fromkeys(candidate_ids))
    if unique_candidate_ids:
        placeholders = ",".join("?" for _ in unique_candidate_ids)
        with connect() as conn:
            target_rows = rows(
                conn,
                f"SELECT protein_id, gene_id, transcript_id, sequence FROM proteins WHERE species_id = ? AND protein_id IN ({placeholders})",
                [request.species_id, *unique_candidate_ids],
            )
        target_by_id = {str(item["protein_id"]): item for item in target_rows}

    sequences: list[dict[str, Any]] = []
    for item in references:
        sequence = clean_protein_sequence(item.get("sequence"))
        if len(sequence) < 20:
            continue
        symbol = str(item.get("gene_symbol") or item.get("gene_id") or item["protein_id"])
        transcript = str(item.get("transcript_id") or item["protein_id"])
        sequences.append({"id": f"Q{len(sequences) + 1}", "label": f"{symbol} · {transcript}", "kind": "query", "sequence": sequence})
    query_count = len(sequences)
    for protein_id in unique_candidate_ids:
        item = target_by_id.get(protein_id)
        if not item:
            continue
        sequence = clean_protein_sequence(item.get("sequence"))
        if len(sequence) < 20:
            continue
        gene_id = str(item.get("gene_id") or protein_id)
        transcript = str(item.get("transcript_id") or protein_id)
        sequences.append({"id": f"T{len(sequences) - query_count + 1}", "label": f"{gene_id} · {transcript}", "kind": "candidate", "sequence": sequence})

    if len(sequences) < 3:
        raise ValueError("at least three valid protein sequences are required for a similar-gene tree")
    return {
        "species_id": request.species_id,
        "gene_id": request.gene_id,
        "model": request.model,
        "analysis": request.analysis,
        "query_protein_count": query_count,
        "sequences": sequences,
        "reference_hits": [],
        "reference": reference,
    }


def create_similar_phylogeny_job(
    species_id: str,
    similar_job_id: str,
    reference_protein_ids: list[str],
    model: TreeModel = "LG+G",
) -> Dict[str, Any]:
    return create_phylogeny_job(
        PhylogenyRequest(
            species_id=species_id,
            gene_id="Selected similar-gene inputs",
            model=model,
            analysis="similar",
            similar_job_id=similar_job_id,
            reference_protein_ids=reference_protein_ids,
        )
    )


def create_phylogeny_job(request: PhylogenyRequest) -> Dict[str, Any]:
    job_id = str(uuid.uuid4())
    path = job_path(job_id)
    path.mkdir(parents=True, exist_ok=False)
    created_at = utc_now()
    write_json(path / "request.json", request.model_dump())
    write_json(
        path / "status.json",
        {
            "job_id": job_id,
            "status": "queued",
            "created_at": created_at,
            "updated_at": created_at,
            "species_id": request.species_id,
            "gene_id": request.gene_id,
            "model": request.model,
            "sequence_count": 0,
        },
    )
    return {
        "job_id": job_id,
        "status": "queued",
        "created_at": created_at,
        "species_id": request.species_id,
        "gene_id": request.gene_id,
        "model": request.model,
        "sequence_count": 0,
    }


def write_fasta(path: Path, sequences: list[Dict[str, Any]]) -> None:
    with path.open("w", encoding="ascii") as handle:
        for item in sequences:
            handle.write(f">{item['id']}\n{item['sequence']}\n")


def selected_model(report_path: Path, requested_model: str) -> str:
    if requested_model != "Auto" or not report_path.exists():
        return requested_model
    report = report_path.read_text(encoding="utf-8", errors="replace")
    match = re.search(r"Best-fit model according to BIC:\s*([^\s]+)", report)
    return match.group(1) if match else "Auto"


def alignment_length(path: Path) -> int:
    sequences = []
    current = []
    for line in path.read_text(encoding="ascii", errors="replace").splitlines():
        if line.startswith(">"):
            if current:
                sequences.append("".join(current))
                current = []
        else:
            current.append(line.strip())
    if current:
        sequences.append("".join(current))
    return len(sequences[0]) if sequences else 0


def execute_phylogeny_job(job_id: str) -> None:
    path = job_path(job_id)
    request = PhylogenyRequest(**read_json(path / "request.json"))
    status_path = path / "status.json"
    started_at = utc_now()
    write_json(
        status_path,
        {
            **read_json(status_path),
            "status": "running",
            "started_at": started_at,
            "updated_at": started_at,
        },
    )
    try:
        prepared = prepare_tree_input(request)
        write_json(path / "input.json", prepared)
    except Exception as exc:
        write_json(status_path, {**read_json(status_path), "status": "failed", "error": str(exc), "finished_at": utc_now(), "updated_at": utc_now()})
        return
    write_json(status_path, {**read_json(status_path), "sequence_count": len(prepared["sequences"]), "updated_at": utc_now()})
    input_fasta = path / "input.fa"
    alignment_fasta = path / "alignment.fa"
    prefix = path / "tree"
    log_path = path / "run.log"
    try:
        write_fasta(input_fasta, prepared["sequences"])
        with alignment_fasta.open("w", encoding="ascii") as alignment, log_path.open("w", encoding="utf-8") as log:
            subprocess.run(
                [str(ORTHOLOG_BIN / "mafft"), "--auto", "--quiet", str(input_fasta)],
                check=True,
                stdout=alignment,
                stderr=log,
                timeout=90,
            )
            model_args = ["-m", request.model] if request.model != "Auto" else ["-m", "MFP", "-mset", "LG,WAG,JTT"]
            subprocess.run(
                [
                    str(ORTHOLOG_BIN / "iqtree3"),
                    "-s", str(alignment_fasta),
                    *model_args,
                    "-nt", "1",
                    "-seed", "1",
                    "-quiet",
                    "-redo",
                    "-pre", str(prefix),
                ],
                check=True,
                stdout=log,
                stderr=subprocess.STDOUT,
                timeout=120,
            )

        tree_path = Path(f"{prefix}.treefile")
        report_path = Path(f"{prefix}.iqtree")
        if not tree_path.exists():
            raise RuntimeError("IQ-TREE did not create a Newick tree file")
        result = {
            "species_id": request.species_id,
            "gene_id": request.gene_id,
            "analysis": prepared["analysis"],
            "query_protein_id": prepared.get("query_protein_id"),
            "query_transcript_id": prepared.get("query_transcript_id"),
            "query_protein_count": prepared["query_protein_count"],
            "requested_model": request.model,
            "selected_model": selected_model(report_path, request.model),
            "alignment_length": alignment_length(alignment_fasta),
            "newick": tree_path.read_text(encoding="utf-8").strip(),
            "leaves": [{key: item[key] for key in ("id", "label", "kind")} for item in prepared["sequences"]],
            "reference_hits": prepared["reference_hits"],
        }
        write_json(path / "results.json", result)
        finished_at = utc_now()
        write_json(
            status_path,
            {
                **read_json(status_path),
                "status": "completed",
                "finished_at": finished_at,
                "updated_at": finished_at,
                "result_count": len(prepared["sequences"]),
            },
        )
    except Exception as exc:
        failed_at = utc_now()
        write_json(
            status_path,
            {
                **read_json(status_path),
                "status": "failed",
                "failed_at": failed_at,
                "updated_at": failed_at,
                "error": str(exc),
            },
        )


def get_phylogeny_job(job_id: str) -> Dict[str, Any]:
    path = job_path(job_id)
    status_path = path / "status.json"
    if not status_path.exists():
        raise FileNotFoundError(job_id)
    status = read_json(status_path)
    result_path = path / "results.json"
    status["results"] = read_json(result_path) if status["status"] == "completed" and result_path.exists() else None
    return status
