from __future__ import annotations

import datetime as dt
import json
import re
import subprocess
import tempfile
import uuid
from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field

from ..api_limits import get_api_limits
from ..config import JOB_DIR, ORTHOLOG_BIN
from ..db import connect, rows
from .job_files import write_json
from .protein_similarity import annotation_symbols, reference_gene_symbol
from .reference_lookup import reference_for_species
from .search_jobs import asset_for_program, parse_tabular_results, resolve_project_path


SIMILAR_LIMITS = get_api_limits().similar


class SimilarGeneRequest(BaseModel):
    species_id: str = Field(min_length=1)
    reference_protein_ids: list[str] = Field(min_length=1)


class PairwiseDomainRequest(BaseModel):
    species_id: str = Field(min_length=1)
    reference_protein_id: str = Field(min_length=1)
    target_protein_id: str = Field(min_length=1)


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def job_root() -> Path:
    path = JOB_DIR / "similar"
    path.mkdir(parents=True, exist_ok=True)
    return path


def pairwise_job_root() -> Path:
    path = JOB_DIR / "similar_pairwise"
    path.mkdir(parents=True, exist_ok=True)
    return path


def job_path(job_id: str) -> Path:
    if not re.fullmatch(r"[a-f0-9-]{36}", job_id):
        raise ValueError("invalid job_id")
    return job_root() / job_id


def pairwise_job_path(job_id: str) -> Path:
    if not re.fullmatch(r"[a-f0-9-]{36}", job_id):
        raise ValueError("invalid job_id")
    return pairwise_job_root() / job_id


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def normalized_protein_ids(values: list[str]) -> list[str]:
    unique: list[str] = []
    seen: set[str] = set()
    for value in values:
        protein_id = str(value).strip()
        if protein_id and protein_id not in seen:
            seen.add(protein_id)
            unique.append(protein_id)
    if not unique:
        raise ValueError("select at least one reference protein isoform")
    if len(unique) > SIMILAR_LIMITS.max_selected_isoforms:
        raise ValueError(f"a maximum of {SIMILAR_LIMITS.max_selected_isoforms} reference isoforms is allowed")
    return unique


def clean_protein_sequence(value: object) -> str:
    return re.sub(r"[^A-Za-z]", "", str(value or "")).upper().rstrip("*")


def selected_reference_proteins(species_id: str, protein_ids: list[str]) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    with connect() as conn:
        reference = reference_for_species(conn, species_id)
        placeholders = ",".join("?" for _ in protein_ids)
        records = rows(
            conn,
            f"""
            SELECT protein_id, gene_id, transcript_id, description, length, sequence
            FROM reference_proteins
            WHERE reference_id = ? AND protein_id IN ({placeholders})
            """,
            [reference["reference_id"], *protein_ids],
        )
        by_id = {str(item["protein_id"]): item for item in records}
        missing = [protein_id for protein_id in protein_ids if protein_id not in by_id]
        if missing:
            raise ValueError(f"reference protein not found: {', '.join(missing[:3])}")
        symbols = annotation_symbols(conn, str(reference["reference_id"]), sorted({str(item["gene_id"]) for item in records if item.get("gene_id")}))
    selected = []
    for protein_id in protein_ids:
        item = by_id[protein_id]
        item["gene_symbol"] = symbols.get(str(item.get("gene_id") or "")) or reference_gene_symbol(item.get("description"))
        selected.append(item)
    return reference, selected


def create_similar_job(request: SimilarGeneRequest) -> dict[str, Any]:
    protein_ids = normalized_protein_ids(request.reference_protein_ids)
    reference, proteins = selected_reference_proteins(request.species_id, protein_ids)
    asset = asset_for_program(request.species_id, "diamond")
    job_id = str(uuid.uuid4())
    path = job_path(job_id)
    path.mkdir(parents=True, exist_ok=False)
    created_at = utc_now()
    (path / "query.fa").write_text(
        "".join(f">{item['protein_id']}\n{item['sequence']}\n" for item in proteins),
        encoding="ascii",
    )
    write_json(
        path / "request.json",
        {"species_id": request.species_id, "reference_protein_ids": protein_ids, "reference_id": reference["reference_id"]},
    )
    write_json(
        path / "status.json",
        {
            "job_id": job_id,
            "status": "queued",
            "created_at": created_at,
            "updated_at": created_at,
            "species_id": request.species_id,
            "reference_id": reference["reference_id"],
            "selected_isoforms": len(proteins),
            "result_count": 0,
        },
    )
    return {"job_id": job_id, "status": "queued", "created_at": created_at, "species_id": request.species_id, "reference": reference, "selected_isoforms": len(proteins), "asset": asset}


def execute_similar_job(job_id: str) -> None:
    path = job_path(job_id)
    request = read_json(path / "request.json")
    status_path = path / "status.json"
    created_at = read_json(status_path)["created_at"]
    started_at = utc_now()
    status = {
        "job_id": job_id,
        "status": "running",
        "created_at": created_at,
        "started_at": started_at,
        "updated_at": started_at,
        "species_id": request["species_id"],
        "reference_id": request["reference_id"],
        "selected_isoforms": len(request["reference_protein_ids"]),
        "result_count": 0,
    }
    write_json(status_path, status)
    try:
        asset = asset_for_program(request["species_id"], "diamond")
        command = [
            str(ORTHOLOG_BIN / "diamond"), "blastp",
            "--query", str(path / "query.fa"),
            "--db", resolve_project_path(asset["index_prefix"]),
            "--out", str(path / "results.tsv"),
            "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore",
            "--max-target-seqs", str(SIMILAR_LIMITS.raw_hits_per_isoform),
            "--evalue", "1e-5", "--threads", "2", "--quiet",
        ]
        with (path / "run.log").open("w", encoding="utf-8") as log:
            subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT, timeout=SIMILAR_LIMITS.timeout_seconds)
        payload = build_result_payload(request["species_id"], request["reference_id"], request["reference_protein_ids"], parse_tabular_results(path / "results.tsv"))
        write_json(path / "result.json", payload)
        finished_at = utc_now()
        status.update({"status": "completed", "finished_at": finished_at, "updated_at": finished_at, "result_count": payload["candidate_count"]})
        write_json(status_path, status)
    except Exception as exc:
        failed_at = utc_now()
        status.update({"status": "failed", "failed_at": failed_at, "updated_at": failed_at, "error": str(exc)})
        write_json(status_path, status)


def build_result_payload(species_id: str, reference_id: str, protein_ids: list[str], hits: list[dict[str, Any]]) -> dict[str, Any]:
    reference, references = selected_reference_proteins(species_id, protein_ids)
    target_ids = sorted({str(hit["sseqid"]) for hit in hits})
    target_by_id: dict[str, dict[str, Any]] = {}
    if target_ids:
        placeholders = ",".join("?" for _ in target_ids)
        with connect() as conn:
            records = rows(
                conn,
                f"""
                SELECT p.protein_id, p.gene_id, p.transcript_id, p.length AS subject_length,
                       g.seqid, g.start, g.end
                FROM proteins p
                LEFT JOIN genes g ON g.species_id = p.species_id AND g.gene_id = p.gene_id
                WHERE p.species_id = ? AND p.protein_id IN ({placeholders})
                """,
                [species_id, *target_ids],
            )
        target_by_id = {str(item["protein_id"]): item for item in records}
    hits_by_query: dict[str, list[dict[str, Any]]] = {}
    for hit in hits:
        hits_by_query.setdefault(str(hit["qseqid"]), []).append(hit)
    groups = []
    for reference_protein in references:
        query_id = str(reference_protein["protein_id"])
        candidates = []
        seen_genes: set[str] = set()
        for hit in hits_by_query.get(query_id, []):
            target = target_by_id.get(str(hit["sseqid"]))
            if not target or not target.get("gene_id") or str(target["gene_id"]) in seen_genes:
                continue
            seen_genes.add(str(target["gene_id"]))
            candidates.append({
                "rank": len(candidates) + 1,
                "gene_id": target["gene_id"],
                "transcript_id": target.get("transcript_id") or target["protein_id"],
                "protein_id": target["protein_id"],
                "pident": hit["pident"],
                "alignment_length": hit["length"],
                "query_length": reference_protein["length"],
                "subject_length": target.get("subject_length"),
                "evalue": hit["evalue"],
                "bitscore": hit["bitscore"],
                "seqid": target.get("seqid"),
                "start": target.get("start"),
                "end": target.get("end"),
            })
            if len(candidates) >= SIMILAR_LIMITS.candidates_per_isoform:
                break
        groups.append({
            "query_protein": {
                "protein_id": query_id,
                "gene_id": reference_protein.get("gene_id"),
                "transcript_id": reference_protein.get("transcript_id") or query_id,
                "gene_symbol": reference_protein.get("gene_symbol"),
                "length": reference_protein["length"],
            },
            "candidates": candidates,
        })
    candidate_gene_ids = [str(candidate["gene_id"]) for group in groups for candidate in group["candidates"]]
    return {
        "species_id": species_id,
        "reference": reference,
        "reference_id": reference_id,
        "method": "Reference protein isoforms -> DIAMOND blastp -> top distinct target genes",
        "groups": groups,
        "candidate_count": sum(len(group["candidates"]) for group in groups),
        "expression": expression_matrix(species_id, candidate_gene_ids),
    }


def get_similar_job(job_id: str) -> dict[str, Any]:
    path = job_path(job_id)
    status_path = path / "status.json"
    if not status_path.exists():
        raise FileNotFoundError(job_id)
    status = read_json(status_path)
    if status.get("status") == "completed":
        status["results"] = read_json(path / "result.json") if (path / "result.json").exists() else None
    return status


def expression_matrix(species_id: str, gene_ids: list[str]) -> dict[str, Any]:
    unique_gene_ids = list(dict.fromkeys(gene_ids))
    if not unique_gene_ids:
        return {"available": False, "samples": [], "values": {}}
    placeholders = ",".join("?" for _ in unique_gene_ids)
    with connect() as conn:
        records = rows(
            conn,
            f"""
            SELECT c.gene_id, c.sample_id, c.value,
                   s.display_name, s.condition, s.treatment, s.timepoint, s.replicate,
                   r.metadata_json, COALESCE(r.row_order, 999999) AS row_order
            FROM expression_counts c
            LEFT JOIN samples s ON s.species_id = c.species_id AND s.sample_id = c.sample_id
            LEFT JOIN raw_sample_metadata r ON r.species_id = c.species_id AND r.sample_id = c.sample_id
            WHERE c.species_id = ?
              AND c.value_type = 'cpm'
              AND c.gene_id IN ({placeholders})
            ORDER BY row_order, c.sample_id
            """,
            [species_id, *unique_gene_ids],
        )
    samples: dict[str, dict[str, Any]] = {}
    values: dict[str, dict[str, float]] = {gene_id: {} for gene_id in unique_gene_ids}
    for item in records:
        sample_id = str(item["sample_id"])
        samples.setdefault(sample_id, {
            "sample_id": sample_id,
            "display_name": item.get("display_name") or sample_id,
            "condition": item.get("condition"),
            "treatment": item.get("treatment"),
            "timepoint": item.get("timepoint"),
            "replicate": item.get("replicate"),
            "metadata": json.loads(item["metadata_json"]) if item.get("metadata_json") else {},
            "row_order": item.get("row_order"),
        })
        values.setdefault(str(item["gene_id"]), {})[sample_id] = float(item["value"])
    return {"available": bool(records), "value_type": "cpm", "transform": "log2(CPM + 1)", "samples": list(samples.values()), "values": values}


def create_pairwise_domain_job(request: PairwiseDomainRequest) -> dict[str, Any]:
    with connect() as conn:
        reference = reference_for_species(conn, request.species_id)
        query = rows(
            conn,
            "SELECT protein_id FROM reference_proteins WHERE reference_id = ? AND protein_id = ?",
            (reference["reference_id"], request.reference_protein_id),
        )
        target = rows(
            conn,
            "SELECT protein_id FROM proteins WHERE species_id = ? AND protein_id = ?",
            (request.species_id, request.target_protein_id),
        )
    if not query or not target:
        raise ValueError("reference or target protein was not found")
    job_id = str(uuid.uuid4())
    path = pairwise_job_path(job_id)
    path.mkdir(parents=True, exist_ok=False)
    created_at = utc_now()
    payload = request.model_dump() | {"reference_id": reference["reference_id"]}
    write_json(path / "request.json", payload)
    write_json(path / "status.json", {
        "job_id": job_id, "status": "queued", "created_at": created_at, "updated_at": created_at,
        "species_id": request.species_id, "reference_id": reference["reference_id"], "result_count": 0,
    })
    return {"job_id": job_id, "status": "queued", "created_at": created_at, "species_id": request.species_id}


def execute_pairwise_domain_job(job_id: str) -> None:
    path = pairwise_job_path(job_id)
    request = read_json(path / "request.json")
    status_path = path / "status.json"
    created_at = read_json(status_path)["created_at"]
    started_at = utc_now()
    status = {"job_id": job_id, "status": "running", "created_at": created_at, "started_at": started_at, "updated_at": started_at, "species_id": request["species_id"], "reference_id": request["reference_id"], "result_count": 0}
    write_json(status_path, status)
    try:
        result = pairwise_domain_payload(request["species_id"], request["reference_id"], request["reference_protein_id"], request["target_protein_id"])
        write_json(path / "result.json", result)
        finished_at = utc_now()
        status.update({"status": "completed", "finished_at": finished_at, "updated_at": finished_at, "result_count": len(result.get("shared_domains", []))})
        write_json(status_path, status)
    except Exception as exc:
        failed_at = utc_now()
        status.update({"status": "failed", "failed_at": failed_at, "updated_at": failed_at, "error": str(exc)})
        write_json(status_path, status)


def pairwise_domain_payload(species_id: str, reference_id: str, reference_protein_id: str, target_protein_id: str) -> dict[str, Any]:
    with connect() as conn:
        query = rows(conn, "SELECT protein_id, gene_id, transcript_id, sequence FROM reference_proteins WHERE reference_id = ? AND protein_id = ?", (reference_id, reference_protein_id))
        target = rows(conn, "SELECT protein_id, gene_id, transcript_id, sequence FROM proteins WHERE species_id = ? AND protein_id = ?", (species_id, target_protein_id))
        query_domains = rows(conn, "SELECT pfam_id, domain_name, start, end, evalue FROM reference_pfam_domains WHERE reference_id = ? AND protein_id = ? ORDER BY start, evalue", (reference_id, reference_protein_id))
        target_domains = rows(conn, "SELECT pfam_id, domain_name, start, end, i_evalue AS evalue FROM pfam_domains WHERE species_id = ? AND protein_id = ? ORDER BY start, i_evalue", (species_id, target_protein_id))
    if not query or not target:
        raise ValueError("reference or target protein was not found")
    query_length = len(clean_protein_sequence(query[0]["sequence"]))
    target_length = len(clean_protein_sequence(target[0]["sequence"]))
    if not query_domains:
        return {
            "available": False,
            "reason": "Reference Pfam annotations are not built for this reference proteome.",
            "query_protein_id": reference_protein_id,
            "target_protein_id": target_protein_id,
            "query_transcript_id": query[0].get("transcript_id") or reference_protein_id,
            "target_transcript_id": target[0].get("transcript_id") or target_protein_id,
            "query_length": query_length,
            "target_length": target_length,
            "shared_domains": [],
        }
    target_by_accession: dict[str, list[dict[str, Any]]] = {}
    for domain in target_domains:
        target_by_accession.setdefault(str(domain["pfam_id"]).split(".", 1)[0], []).append(domain)
    shared = []
    for query_domain in query_domains:
        accession = str(query_domain["pfam_id"]).split(".", 1)[0]
        candidates = target_by_accession.get(accession, [])
        if not candidates:
            continue
        candidate = min(candidates, key=lambda item: float(item.get("evalue") or float("inf")))
        identity, aligned_positions, alignment_error = pairwise_identity(
            str(query[0]["sequence"])[int(query_domain["start"]) - 1:int(query_domain["end"])],
            str(target[0]["sequence"])[int(candidate["start"]) - 1:int(candidate["end"])],
        )
        shared.append({
            "pfam_id": query_domain["pfam_id"],
            "domain_name": query_domain.get("domain_name") or candidate.get("domain_name"),
            "query_range": f"{query_domain['start']}-{query_domain['end']}",
            "target_range": f"{candidate['start']}-{candidate['end']}",
            "identity": identity,
            "aligned_positions": aligned_positions,
            "alignment_error": alignment_error,
        })
    return {
        "available": True,
        "query_protein_id": reference_protein_id,
        "target_protein_id": target_protein_id,
        "query_transcript_id": query[0].get("transcript_id") or reference_protein_id,
        "target_transcript_id": target[0].get("transcript_id") or target_protein_id,
        "query_length": query_length,
        "target_length": target_length,
        "shared_domains": shared,
    }


def pairwise_identity(query_sequence: str, target_sequence: str) -> tuple[float | None, int | None, str | None]:
    if not query_sequence or not target_sequence:
        return None, None, "empty domain sequence"
    try:
        with tempfile.TemporaryDirectory(prefix="pgdb_pairwise_") as directory:
            input_path = Path(directory) / "pair.fa"
            input_path.write_text(f">query\n{query_sequence}\n>target\n{target_sequence}\n", encoding="ascii")
            proc = subprocess.run([str(ORTHOLOG_BIN / "mafft"), "--auto", "--quiet", str(input_path)], check=True, capture_output=True, text=True, timeout=60)
        sequences: dict[str, str] = {}
        identifier = ""
        for line in proc.stdout.splitlines():
            if line.startswith(">"):
                identifier = line[1:].strip()
                sequences[identifier] = ""
            elif identifier:
                sequences[identifier] += line.strip()
        aligned = [(left, right) for left, right in zip(sequences.get("query", ""), sequences.get("target", "")) if left != "-" and right != "-"]
        if not aligned:
            return None, 0, None
        return round(sum(left == right for left, right in aligned) * 100 / len(aligned), 2), len(aligned), None
    except Exception as exc:
        return None, None, str(exc)


def get_pairwise_domain_job(job_id: str) -> dict[str, Any]:
    path = pairwise_job_path(job_id)
    status_path = path / "status.json"
    if not status_path.exists():
        raise FileNotFoundError(job_id)
    status = read_json(status_path)
    if status.get("status") == "completed":
        status["results"] = read_json(path / "result.json") if (path / "result.json").exists() else None
    return status
