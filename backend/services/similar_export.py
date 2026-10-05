from __future__ import annotations

import csv
import uuid
from pathlib import Path

from ..config import JOB_DIR
from ..db import connect, rows
from .search_jobs import parse_tabular_results
from .similar_jobs import get_similar_job, job_path, pairwise_domain_payload, read_json, utc_now, write_json


EXPORT_COLUMNS = [
    "record_type", "species_id", "reference_id", "query_gene_id", "query_transcript_id", "query_protein_id",
    "rank", "target_gene_id", "target_transcript_id", "target_protein_id",
    "pident", "alignment_length", "query_length", "subject_length", "evalue", "bitscore",
    "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "seqid", "start", "end",
    "pfam_id", "domain_name", "query_range", "target_range", "identity", "aligned_positions",
    "pairwise_status", "alignment_error",
]


def export_path(job_id: str) -> Path:
    if str(uuid.UUID(job_id)) != job_id:
        raise ValueError("invalid job_id")
    return JOB_DIR / "similar_export" / job_id


def create_export(species_id: str, source_job_id: str, protein_id: str) -> dict:
    source = get_similar_job(source_job_id)
    if source.get("species_id") != species_id:
        raise FileNotFoundError(source_job_id)
    if source.get("status") != "completed" or not source.get("results"):
        raise ValueError("Candidate search must complete before exporting.")
    result = source["results"]
    group = next((g for g in result["groups"] if g["query_protein"]["protein_id"] == protein_id), None)
    if group is None:
        raise ValueError("Query protein does not belong to this search.")
    raw_hits = [h for h in parse_tabular_results(job_path(source_job_id) / "results.tsv") if h["qseqid"] == protein_id]
    job_id = str(uuid.uuid4())
    path = export_path(job_id)
    path.mkdir(parents=True)
    write_json(path / "request.json", {"species_id": species_id, "reference_id": result["reference_id"],
                                      "group": group, "raw_hits": raw_hits})
    status = {"job_id": job_id, "species_id": species_id, "status": "queued", "created_at": utc_now(),
              "query_protein_id": protein_id, "source_job_id": source_job_id}
    write_json(path / "status.json", status)
    return status


def get_export(species_id: str, job_id: str) -> dict:
    status = read_json(export_path(job_id) / "status.json")
    if status.get("species_id") != species_id:
        raise FileNotFoundError(job_id)
    return status


def export_records(request: dict):
    species_id, reference_id = request["species_id"], request["reference_id"]
    group, raw_hits = request["group"], request["raw_hits"]
    query = group["query_protein"]
    base = {"species_id": species_id, "reference_id": reference_id, "query_gene_id": query.get("gene_id"),
            "query_transcript_id": query.get("transcript_id"), "query_protein_id": query["protein_id"],
            "query_length": query.get("length")}
    targets = dict.fromkeys([c["protein_id"] for c in group["candidates"]] + [h["sseqid"] for h in raw_hits])
    with connect() as conn:
        target_info = {p["protein_id"]: p for p in rows(conn, "SELECT protein_id, gene_id, transcript_id, length FROM proteins WHERE species_id = ? AND protein_id IN (" + ",".join("?" for _ in targets) + ")", [species_id, *targets])} if targets else {}
    for candidate in group["candidates"]:
        values = {k: v for k, v in candidate.items() if k in EXPORT_COLUMNS}
        yield {**base, **values, "record_type": "candidate_rank", "target_gene_id": candidate["gene_id"],
               "target_transcript_id": candidate.get("transcript_id"), "target_protein_id": candidate["protein_id"]}
    for index, hit in enumerate(raw_hits, 1):
        target = target_info.get(hit["sseqid"], {})
        values = {k: v for k, v in hit.items() if k in EXPORT_COLUMNS}
        yield {**base, **values, "record_type": "raw_hit", "rank": index, "target_protein_id": hit["sseqid"],
               "target_gene_id": target.get("gene_id"), "target_transcript_id": target.get("transcript_id"),
               "subject_length": target.get("length"), "alignment_length": hit["length"]}
    for protein_id in targets:
        target = target_info.get(protein_id, {})
        context = {**base, "record_type": "pairwise_domain", "target_protein_id": protein_id,
                   "target_gene_id": target.get("gene_id"), "target_transcript_id": target.get("transcript_id"),
                   "subject_length": target.get("length")}
        comparison = pairwise_domain_payload(species_id, reference_id, query["protein_id"], protein_id)
        if not comparison["shared_domains"]:
            yield {**context, "pairwise_status": "no_shared_domains" if comparison["available"] else "no_reference_domains"}
        for domain in comparison["shared_domains"]:
            if domain.get("alignment_error"):
                raise RuntimeError(f"Domain alignment failed for {protein_id}: {domain['alignment_error']}")
            yield {**context, **domain, "pairwise_status": "completed"}


def execute_export(job_id: str) -> None:
    path = export_path(job_id)
    status = read_json(path / "status.json")
    started_at = utc_now()
    status.update(status="running", started_at=started_at, updated_at=started_at)
    write_json(path / "status.json", status)
    try:
        count = 0
        with (path / "results.tmp").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=EXPORT_COLUMNS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for record in export_records(read_json(path / "request.json")):
                writer.writerow(record)
                count += 1
        (path / "results.tmp").replace(path / "results.tsv")
        status.update(status="completed", result_count=count)
    except Exception as exc:
        status.update(status="failed", error=str(exc))
    finished_at = utc_now()
    status.update(finished_at=finished_at, updated_at=finished_at)
    if status["status"] == "failed":
        status["failed_at"] = finished_at
    write_json(path / "status.json", status)
