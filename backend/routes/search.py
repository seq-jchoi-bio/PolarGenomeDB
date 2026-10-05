from __future__ import annotations

from fastapi import APIRouter, HTTPException, Query
from fastapi.responses import PlainTextResponse

from ..db import connect, rows
from .genes import fetch_fasta_interval, reverse_complement
from ..services.search_jobs import SearchRequest, create_search_job, execute_search_job, get_search_job, raw_search_results
from ..services.job_manager import job_manager
from ..runtime import load_runtime_json


router = APIRouter()
SEARCH_RAW_TSV_HEADER = "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n"


@router.get("/species/{species_id}/search-examples")
def get_search_examples(
    species_id: str,
    program: str = Query(pattern="^(blastn|blastp|diamond)$"),
) -> dict:
    record_count = 40 if program == "diamond" else 3
    if program == "blastn":
        with connect() as conn:
            genes = rows(
                conn,
                """
                SELECT gene_id, seqid, start, end, strand
                FROM genes
                WHERE species_id = ? AND end - start + 1 >= 400
                ORDER BY gene_id
                LIMIT ?
                """,
                (species_id, record_count),
            )
        records = []
        for gene in genes:
            sequence, _ = fetch_fasta_interval(species_id, gene["seqid"], int(gene["start"]), int(gene["start"]) + 399)
            if gene.get("strand") == "-":
                sequence = reverse_complement(sequence)
            records.append({"source_id": gene["gene_id"], "sequence": sequence[:400]})
    else:
        with connect() as conn:
            protein_rows = rows(
                conn,
                """
                SELECT protein_id AS source_id, gene_id, sequence
                FROM proteins
                WHERE species_id = ? AND length(sequence) >= 100
                ORDER BY gene_id, protein_id
                LIMIT ?
                """,
                (species_id, record_count * 20),
            )
        seen_genes = set()
        records = []
        for record in protein_rows:
            if record["gene_id"] in seen_genes:
                continue
            seen_genes.add(record["gene_id"])
            records.append({"source_id": record["source_id"], "sequence": str(record["sequence"])[:100]})
            if len(records) == record_count:
                break

    if len(records) < record_count:
        raise HTTPException(status_code=404, detail=f"not enough example sequences for {species_id}/{program}")
    return {"species_id": species_id, "program": program, "records": records}


@router.get("/species/{species_id}/search-assets")
def get_search_assets(species_id: str) -> dict:
    try:
        return load_runtime_json(f"{species_id}/search_assets.json")
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"species not found: {species_id}") from None


@router.post("/search/blast")
def blast_search(request: SearchRequest) -> dict:
    if request.program not in {"blastn", "blastp"}:
        raise HTTPException(status_code=400, detail="BLAST endpoint accepts program blastn or blastp")
    try:
        job = create_search_job(request)
        job_manager.submit("search", job["job_id"], execute_search_job, job["job_id"])
        return {key: value for key, value in job.items() if key != "asset"}
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"species not found: {request.species_id}") from None
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from None


@router.post("/search/diamond")
def diamond_search(request: SearchRequest) -> dict:
    if request.program != "diamond":
        raise HTTPException(status_code=400, detail="DIAMOND endpoint accepts program diamond")
    try:
        job = create_search_job(request)
        job_manager.submit("search", job["job_id"], execute_search_job, job["job_id"])
        return {key: value for key, value in job.items() if key != "asset"}
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"species not found: {request.species_id}") from None
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from None


@router.get("/search/jobs/{job_id}")
def read_search_job(job_id: str) -> dict:
    try:
        return get_search_job(job_id)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"job not found: {job_id}") from None
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from None


@router.get("/search/jobs/{job_id}/raw-results")
def download_raw_search_results(job_id: str, query_id: str | None = Query(default=None, min_length=1)) -> PlainTextResponse:
    try:
        lines = raw_search_results(job_id, query_id)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail="raw BLAST results are not available for this job") from None
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from None
    filename = f"polar-genome-db-{query_id or 'all-queries'}-raw-hits.tsv"
    return PlainTextResponse(
        SEARCH_RAW_TSV_HEADER + "".join(lines),
        media_type="text/tab-separated-values",
        headers={"Content-Disposition": f'attachment; filename="{filename}"'},
    )
