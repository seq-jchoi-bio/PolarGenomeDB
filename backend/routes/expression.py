from __future__ import annotations

from typing import Optional

from fastapi import APIRouter, HTTPException, Query

from ..api_limits import get_api_limits
from ..db import connect, parse_json, rows


router = APIRouter()
EXPRESSION_LIMITS = get_api_limits().expression


@router.get("/species/{species_id}/expression/matrix")
def get_expression_matrix(
    species_id: str,
    genes: str = Query(min_length=1, description="Comma-separated gene IDs."),
    value_type: str = Query(default="cpm", min_length=1),
    samples: Optional[str] = Query(default=None, description="Optional comma-separated sample IDs."),
) -> dict:
    gene_ids = [item.strip() for item in genes.split(",") if item.strip()]
    sample_ids = [item.strip() for item in samples.split(",") if item.strip()] if samples else []
    if not gene_ids:
        raise HTTPException(status_code=400, detail="at least one gene ID is required")
    if len(gene_ids) > EXPRESSION_LIMITS.max_genes_per_matrix_request:
        raise HTTPException(
            status_code=400,
            detail=f"too many genes; max {EXPRESSION_LIMITS.max_genes_per_matrix_request}",
        )

    gene_placeholders = ",".join("?" for _ in gene_ids)
    sample_clause = ""
    params: list[object] = [species_id, *gene_ids, value_type]
    if sample_ids:
        sample_placeholders = ",".join("?" for _ in sample_ids)
        sample_clause = f" AND c.sample_id IN ({sample_placeholders})"
        params.extend(sample_ids)

    with connect() as conn:
        species_exists = conn.execute("SELECT 1 FROM species WHERE species_id = ?", (species_id,)).fetchone()
        if not species_exists:
            raise HTTPException(status_code=404, detail=f"species not found: {species_id}")

        gene_rows = rows(
            conn,
            f"""
            SELECT gene_id, seqid, start, end, strand, source, attributes_json
            FROM genes
            WHERE species_id = ? AND gene_id IN ({gene_placeholders})
            ORDER BY seqid, start, end, gene_id
            """,
            [species_id, *gene_ids],
        )
        found_gene_ids = {item["gene_id"] for item in gene_rows}
        missing_gene_ids = [gene_id for gene_id in gene_ids if gene_id not in found_gene_ids]

        matrix_rows = rows(
            conn,
            f"""
            SELECT c.gene_id, c.sample_id, s.display_name, s.tissue, s.condition, s.treatment,
                   s.timepoint, s.replicate, c.value, c.source_file
            FROM expression_counts c
            LEFT JOIN samples s
              ON s.species_id = c.species_id AND s.sample_id = c.sample_id
            WHERE c.species_id = ?
              AND c.gene_id IN ({gene_placeholders})
              AND c.value_type = ?
              {sample_clause}
            ORDER BY c.gene_id, c.sample_id
            """,
            params,
        )

    for item in gene_rows:
        item["attributes"] = parse_json(item.pop("attributes_json", None))

    seen_samples: dict[str, dict] = {}
    values: dict[str, dict[str, float]] = {gene_id: {} for gene_id in gene_ids}
    for item in matrix_rows:
        sample_id = item["sample_id"]
        values.setdefault(item["gene_id"], {})[sample_id] = item["value"]
        seen_samples.setdefault(
            sample_id,
            {
                "sample_id": sample_id,
                "display_name": item["display_name"],
                "tissue": item["tissue"],
                "condition": item["condition"],
                "treatment": item["treatment"],
                "timepoint": item["timepoint"],
                "replicate": item["replicate"],
            },
        )

    return {
        "species_id": species_id,
        "value_type": value_type,
        "requested_genes": gene_ids,
        "missing_genes": missing_gene_ids,
        "gene_count": len(gene_rows),
        "sample_count": len(seen_samples),
        "genes": gene_rows,
        "samples": list(seen_samples.values()),
        "values": values,
    }


@router.get("/species/{species_id}/expression/{gene_id}")
def get_gene_expression(species_id: str, gene_id: str) -> dict:
    with connect() as conn:
        gene_exists = conn.execute(
            "SELECT 1 FROM genes WHERE species_id = ? AND gene_id = ?",
            (species_id, gene_id),
        ).fetchone()
        if not gene_exists:
            raise HTTPException(status_code=404, detail=f"gene not found: {species_id}/{gene_id}")
        data = rows(
            conn,
            """
            SELECT c.sample_id, s.display_name, s.tissue, s.condition, s.treatment, s.timepoint,
                   s.replicate, c.value_type, c.value, c.source_file
            FROM expression_counts c
            LEFT JOIN samples s
              ON s.species_id = c.species_id AND s.sample_id = c.sample_id
            WHERE c.species_id = ? AND c.gene_id = ?
            ORDER BY c.value_type, c.sample_id
            """,
            (species_id, gene_id),
        )
    grouped: dict[str, list[dict]] = {}
    for item in data:
        grouped.setdefault(item.pop("value_type"), []).append(item)
    return {"species_id": species_id, "gene_id": gene_id, "expression": grouped}
