from __future__ import annotations

import csv
from io import StringIO
from urllib.parse import quote

from fastapi import APIRouter, HTTPException, Query
from fastapi.responses import PlainTextResponse

from ..db import connect, row, rows
from ..services.pfam import ensure_pfam_schema, insert_pfam_domains, mark_pfam_scanned, run_hmmscan

router = APIRouter()


@router.get("/species/{species_id}/domains/{gene_id}/export")
def export_gene_domains(species_id: str, gene_id: str, protein_id: str) -> PlainTextResponse:
    columns = [
        "species_id", "gene_id", "transcript_id", "protein_id", "pfam_id",
        "domain_name", "domain_length", "start", "end", "hmm_start", "hmm_end",
        "i_evalue", "c_evalue", "full_evalue", "score", "accuracy", "description", "source",
    ]
    with connect() as conn:
        if not row(conn, "SELECT 1 FROM proteins WHERE species_id = ? AND gene_id = ? AND protein_id = ?",
                   (species_id, gene_id, protein_id)):
            raise HTTPException(status_code=404, detail="Protein not found for this gene.")
        records = rows(conn, f"SELECT {', '.join(columns)} FROM pfam_domains "
                       "WHERE species_id = ? AND protein_id = ? ORDER BY start, i_evalue, score DESC",
                       (species_id, protein_id))
    output = StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=columns, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(records)
    filename = quote(f"{species_id}_{gene_id}_{protein_id}_pfam.tsv", safe="")
    return PlainTextResponse(output.getvalue(), media_type="text/tab-separated-values",
                             headers={"Content-Disposition": f"attachment; filename*=UTF-8''{filename}"})


def normalize_domain(domain: dict) -> dict:
    return {
        "domain": domain.get("domain_name"),
        "domain_name": domain.get("domain_name"),
        "pfam_id": domain.get("pfam_id"),
        "accession": domain.get("pfam_id"),
        "domain_length": domain.get("domain_length"),
        "ali_from": domain.get("start"),
        "ali_to": domain.get("end"),
        "start": domain.get("start"),
        "end": domain.get("end"),
        "hmm_from": domain.get("hmm_start"),
        "hmm_to": domain.get("hmm_end"),
        "i_evalue": domain.get("i_evalue"),
        "c_evalue": domain.get("c_evalue"),
        "full_evalue": domain.get("full_evalue"),
        "domain_score": domain.get("score"),
        "score": domain.get("score"),
        "accuracy": domain.get("accuracy"),
        "description": domain.get("description"),
    }


@router.get("/species/{species_id}/domains/{gene_id}")
def get_gene_domains(
    species_id: str,
    gene_id: str,
    protein_id: str | None = Query(default=None),
    refresh: bool = Query(default=False),
    cached_only: bool = Query(default=False),
    cpu: int = Query(default=4, ge=1, le=16),
) -> dict:
    with connect() as conn:
        protein = row(
            conn,
            """
            SELECT protein_id, gene_id, transcript_id, length, sequence
            FROM proteins
            WHERE species_id = ? AND gene_id = ?
              AND (? IS NULL OR protein_id = ?)
            ORDER BY transcript_id, protein_id
            LIMIT 1
            """,
            (species_id, gene_id, protein_id, protein_id),
        )
        if not protein:
            if not row(conn, "SELECT 1 FROM genes WHERE species_id = ? AND gene_id = ?", (species_id, gene_id)):
                raise HTTPException(status_code=404, detail=f"gene not found: {species_id}/{gene_id}")
            if protein_id:
                raise HTTPException(status_code=404, detail=f"protein does not belong to gene: {protein_id}")
            raise HTTPException(status_code=404, detail=f"protein not found: {species_id}/{gene_id}")

        ensure_pfam_schema(conn)
        cached = [] if refresh else rows(
            conn,
            """
            SELECT pfam_id, domain_name, domain_length, start, end, hmm_start, hmm_end,
                   i_evalue, c_evalue, full_evalue, score, accuracy, description
            FROM pfam_domains
            WHERE species_id = ? AND protein_id = ?
            ORDER BY start, i_evalue, score DESC
            """,
            (species_id, protein["protein_id"]),
        )
        if cached:
            return {
                "species_id": species_id,
                "gene_id": gene_id,
                "query_transcript": protein.get("transcript_id"),
                "protein_id": protein["protein_id"],
                "protein_length": protein["length"],
                "method": "hmmscan against Pfam-A; domains filtered by i-Evalue <= 1e-3",
                "cached": True,
                "domains": [normalize_domain(item) for item in cached],
            }

        scanned = None if refresh else row(
            conn,
            """
            SELECT status, domain_count, updated_at
            FROM pfam_protein_status
            WHERE species_id = ? AND protein_id = ?
            """,
            (species_id, protein["protein_id"]),
        )
        if scanned:
            return {
                "species_id": species_id,
                "gene_id": gene_id,
                "query_transcript": protein.get("transcript_id"),
                "protein_id": protein["protein_id"],
                "protein_length": protein["length"],
                "method": "cached Pfam domain lookup; protein was scanned but no domains passed the cutoff",
                "cached": True,
                "domains": [],
            }

        if cached_only:
            return {
                "species_id": species_id,
                "gene_id": gene_id,
                "query_transcript": protein.get("transcript_id"),
                "protein_id": protein["protein_id"],
                "protein_length": protein["length"],
                "method": "cached Pfam domain lookup only",
                "cached": True,
                "domains": [],
            }

        try:
            domains = run_hmmscan([(protein["protein_id"], protein["sequence"])], cpu=cpu, timeout=180)
        except FileNotFoundError as exc:
            raise HTTPException(status_code=503, detail=str(exc)) from exc
        except Exception as exc:
            raise HTTPException(status_code=500, detail=f"hmmscan failed: {exc}") from exc

        protein_rows = {protein["protein_id"]: protein}
        insert_pfam_domains(conn, species_id, protein_rows, domains)
        mark_pfam_scanned(conn, species_id, protein_rows, domains)
        conn.commit()
        return {
            "species_id": species_id,
            "gene_id": gene_id,
            "query_transcript": protein.get("transcript_id"),
            "protein_id": protein["protein_id"],
            "protein_length": protein["length"],
            "method": "hmmscan against Pfam-A; domains filtered by i-Evalue <= 1e-3",
            "cached": False,
            "domains": [normalize_domain(item) for item in domains],
        }
