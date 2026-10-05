from __future__ import annotations

import subprocess
import csv
from io import StringIO
from urllib.parse import quote

from fastapi import APIRouter, HTTPException, Query
from fastapi.responses import FileResponse, PlainTextResponse
from pydantic import BaseModel

from ..db import connect, parse_json, row, rows
from ..services.protein_similarity import protein_similarity, reference_id_for_protein, result_rows
from ..services.reference_lookup import lookup_similar_genes, resolve_reference_queries
from ..services.similar_jobs import (
    PairwiseDomainRequest,
    SimilarGeneRequest,
    create_pairwise_domain_job,
    create_similar_job,
    execute_pairwise_domain_job,
    execute_similar_job,
    get_pairwise_domain_job,
    get_similar_job,
)
from ..services.phylogeny_jobs import TreeModel, create_similar_phylogeny_job, execute_phylogeny_job
from ..services.job_manager import job_manager
from ..services.similar_export import create_export, execute_export, export_path, get_export


router = APIRouter()


class ReferenceLookupRequest(BaseModel):
    queries: list[str]


class SimilarExportRequest(BaseModel):
    reference_protein_id: str


@router.post("/species/{species_id}/similar/jobs/{source_job_id}/export")
def start_similar_export(species_id: str, source_job_id: str, request: SimilarExportRequest) -> dict:
    try:
        job = create_export(species_id, source_job_id, request.reference_protein_id)
        job_manager.submit("similar_export", job["job_id"], execute_export, job["job_id"])
        return job
    except FileNotFoundError:
        raise HTTPException(404, "Search result not found.") from None
    except ValueError as exc:
        raise HTTPException(400, str(exc)) from None


@router.get("/species/{species_id}/similar/exports/{job_id}")
def similar_export_status(species_id: str, job_id: str) -> dict:
    try:
        return get_export(species_id, job_id)
    except (FileNotFoundError, ValueError):
        raise HTTPException(404, "Export job not found.") from None


@router.get("/species/{species_id}/similar/exports/{job_id}/download")
def download_similar_export(species_id: str, job_id: str) -> FileResponse:
    job = similar_export_status(species_id, job_id)
    if job["status"] != "completed":
        raise HTTPException(409, "Export is not complete.")
    return FileResponse(export_path(job_id) / "results.tsv", media_type="text/tab-separated-values",
                        filename=f"{species_id}_{job['query_protein_id']}_candidate_pairwise.tsv")


class SimilarTreeRequest(BaseModel):
    similar_job_id: str
    reference_protein_ids: list[str]
    model: TreeModel = "LG+G"


@router.post("/species/{species_id}/similar/jobs")
def create_similar_gene_job(species_id: str, request: SimilarGeneRequest) -> dict:
    if request.species_id != species_id:
        raise HTTPException(status_code=400, detail="species_id must match the URL")
    try:
        job = create_similar_job(request)
        job_manager.submit("similar", job["job_id"], execute_similar_job, job["job_id"])
        return {key: value for key, value in job.items() if key != "asset"}
    except (FileNotFoundError, ValueError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from None


@router.get("/species/{species_id}/similar/jobs/{job_id}")
def read_similar_gene_job(species_id: str, job_id: str) -> dict:
    try:
        job = get_similar_job(job_id)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"job not found: {job_id}") from None
    if job.get("species_id") != species_id:
        raise HTTPException(status_code=404, detail=f"job not found: {job_id}")
    return job


@router.post("/species/{species_id}/similar/tree")
def create_similar_gene_tree(species_id: str, request: SimilarTreeRequest) -> dict:
    try:
        job = create_similar_phylogeny_job(
            species_id=species_id,
            similar_job_id=request.similar_job_id,
            reference_protein_ids=request.reference_protein_ids,
            model=request.model,
        )
        job_manager.submit("phylogeny", job["job_id"], execute_phylogeny_job, job["job_id"])
        return job
    except (FileNotFoundError, ValueError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from None


@router.post("/species/{species_id}/similar/pairwise")
def create_pairwise_domain_comparison(species_id: str, request: PairwiseDomainRequest) -> dict:
    if request.species_id != species_id:
        raise HTTPException(status_code=400, detail="species_id must match the URL")
    try:
        job = create_pairwise_domain_job(request)
        job_manager.submit("similar_pairwise", job["job_id"], execute_pairwise_domain_job, job["job_id"])
        return job
    except (FileNotFoundError, ValueError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from None


@router.get("/species/{species_id}/similar/pairwise/{job_id}")
def read_pairwise_domain_comparison(species_id: str, job_id: str) -> dict:
    try:
        job = get_pairwise_domain_job(job_id)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"job not found: {job_id}") from None
    if job.get("species_id") != species_id:
        raise HTTPException(status_code=404, detail=f"job not found: {job_id}")
    return job


@router.post("/species/{species_id}/similar/resolve")
def resolve_reference_genes(species_id: str, request: ReferenceLookupRequest) -> dict:
    try:
        return resolve_reference_queries(species_id, request.queries)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from None


@router.post("/species/{species_id}/similar/lookup")
def lookup_reference_genes(species_id: str, request: ReferenceLookupRequest) -> dict:
    try:
        return lookup_similar_genes(species_id, request.queries)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from None


@router.get("/species/{species_id}/similar/status")
def get_similarity_status(species_id: str) -> dict:
    with connect() as conn:
        species = row(
            conn,
            "SELECT species_id, scientific_name, display_name, taxon_group FROM species WHERE species_id = ?",
            (species_id,),
        )
        if not species:
            raise HTTPException(status_code=404, detail=f"species not found: {species_id}")
        data = rows(
            conn,
            """
            SELECT
              h.reference_id,
              r.display_name AS reference_display_name,
              r.taxon_group AS reference_taxon_group,
              h.method,
              COUNT(*) AS hit_count,
              COUNT(DISTINCT h.gene_id) AS genes_with_hits,
              MIN(h.evalue) AS best_evalue,
              MAX(h.bitscore) AS best_bitscore,
              MAX(h.created_at) AS updated_at
            FROM similarity_hits h
            JOIN reference_species r ON r.reference_id = h.reference_id
            WHERE h.species_id = ?
            GROUP BY h.reference_id, h.method
            ORDER BY h.reference_id, h.method
            """,
            (species_id,),
        )
    return {"species_id": species_id, "species": species, "similarity": data}


@router.get("/species/{species_id}/similar/{gene_id}/export")
def export_protein_similarity(species_id: str, gene_id: str, protein_id: str) -> PlainTextResponse:
    with connect() as conn:
        protein = row(conn, "SELECT transcript_id FROM proteins WHERE species_id = ? AND gene_id = ? AND protein_id = ?",
                      (species_id, gene_id, protein_id))
        if not protein:
            raise HTTPException(status_code=404, detail="Protein not found for this gene.")
        reference_id = reference_id_for_protein(conn, species_id)
        hits = result_rows(conn, species_id, gene_id, protein_id, reference_id, None)
    columns = ["species_id", "gene_id", "query_transcript_id", "query_protein_id", "reference_id",
               "reference_display_name", "rank", "reference_gene_id", "reference_transcript_id",
               "reference_protein_id", "reference_gene_symbol", "pident", "alignment_length", "evalue",
               "bitscore", "reference_seqid", "reference_start", "reference_end", "reference_strand",
               "reference_protein_description", "method", "created_at"]
    output = StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=columns, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for hit in hits:
        writer.writerow({**hit, "species_id": species_id, "gene_id": gene_id,
                         "query_transcript_id": protein["transcript_id"]})
    filename = quote(f"{species_id}_{gene_id}_{protein_id}_reference_similarity.tsv", safe="")
    return PlainTextResponse(output.getvalue(), media_type="text/tab-separated-values",
                             headers={"Content-Disposition": f"attachment; filename*=UTF-8''{filename}"})


@router.get("/species/{species_id}/similar/{gene_id}")
def get_similar_genes(
    species_id: str,
    gene_id: str,
    protein_id: str | None = None,
    reference_id: str | None = None,
    limit: int = Query(default=10, ge=1, le=50),
) -> dict:
    if protein_id:
        try:
            return protein_similarity(species_id, gene_id, protein_id, min(limit, 5))
        except (FileNotFoundError, subprocess.CalledProcessError) as exc:
            raise HTTPException(status_code=503, detail=f"DIAMOND similarity search failed: {exc}") from None
        except ValueError as exc:
            raise HTTPException(status_code=400, detail=str(exc)) from None
    clauses = ["h.species_id = ?", "h.gene_id = ?"]
    params: list[object] = [species_id, gene_id]
    if reference_id:
        clauses.append("h.reference_id = ?")
        params.append(reference_id)
    params.append(limit)
    with connect() as conn:
        if not row(conn, "SELECT 1 FROM genes WHERE species_id = ? AND gene_id = ?", (species_id, gene_id)):
            raise HTTPException(status_code=404, detail=f"gene not found: {species_id}/{gene_id}")
        hits = rows(
            conn,
            f"""
            SELECT
              h.reference_id,
              rs.display_name AS reference_display_name,
              h.reference_gene_id,
              rg.seqid AS reference_seqid,
              rg.start AS reference_start,
              rg.end AS reference_end,
              rg.strand AS reference_strand,
              h.reference_protein_id,
              rp.description AS reference_protein_description,
              h.protein_id AS query_protein_id,
              h.pident,
              h.alignment_length,
              h.evalue,
              h.bitscore,
              h.rank,
              h.method,
              h.created_at,
              rg.attributes_json AS reference_gene_attributes
            FROM similarity_hits h
            JOIN reference_species rs ON rs.reference_id = h.reference_id
            LEFT JOIN reference_proteins rp
              ON rp.reference_id = h.reference_id
             AND rp.protein_id = h.reference_protein_id
            LEFT JOIN reference_genes rg
              ON rg.reference_id = h.reference_id
             AND rg.gene_id = h.reference_gene_id
            WHERE {' AND '.join(clauses)}
            ORDER BY h.rank, h.bitscore DESC, h.evalue ASC
            LIMIT ?
            """,
            params,
        )
    for hit in hits:
        hit["reference_gene_attributes"] = parse_json(hit.pop("reference_gene_attributes", None))
    return {"species_id": species_id, "gene_id": gene_id, "hit_count": len(hits), "hits": hits}
