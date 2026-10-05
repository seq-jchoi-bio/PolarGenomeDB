from __future__ import annotations

from typing import Optional

from fastapi import APIRouter, HTTPException, Query

from ..api_limits import get_api_limits
from ..db import connect, parse_json, rows


router = APIRouter()
FEATURE_LIMITS = get_api_limits().features
GENE_LIMITS = get_api_limits().genes


def parse_feature_types(value: Optional[str]) -> list[str]:
    if not value:
        return []
    return [item.strip() for item in value.split(",") if item.strip()]


@router.get("/species/{species_id}/features/region")
def get_region_features(
    species_id: str,
    seqid: str = Query(min_length=1),
    start: int = Query(ge=1),
    end: int = Query(ge=1),
    feature_type: Optional[str] = Query(default=None, description="Comma-separated feature types, e.g. exon,CDS"),
    limit: int = Query(default=FEATURE_LIMITS.default_page_size, ge=1, le=FEATURE_LIMITS.max_page_size),
    offset: int = Query(default=0, ge=0),
    gene_limit: int = Query(default=GENE_LIMITS.default_page_size, ge=1, le=GENE_LIMITS.max_page_size),
) -> dict:
    if start > end:
        raise HTTPException(status_code=400, detail="start must be less than or equal to end")

    feature_types = parse_feature_types(feature_type)
    feature_clauses = ["species_id = ?", "seqid = ?", "end >= ?", "start <= ?"]
    feature_params: list[object] = [species_id, seqid, start, end]
    if feature_types:
        placeholders = ",".join("?" for _ in feature_types)
        feature_clauses.append(f"feature_type IN ({placeholders})")
        feature_params.extend(feature_types)

    gene_sql = """
        SELECT gene_id, seqid, start, end, strand, source, attributes_json
        FROM genes
        WHERE species_id = ? AND seqid = ? AND end >= ? AND start <= ?
        ORDER BY start, end, gene_id
        LIMIT ?
    """
    feature_sql = f"""
        SELECT feature_id, feature_type, parent_id, gene_id, transcript_id, seqid, start, end, strand, phase, attributes_json
        FROM features
        WHERE {' AND '.join(feature_clauses)}
        ORDER BY start, end, feature_type, feature_id
        LIMIT ? OFFSET ?
    """
    count_sql = f"""
        SELECT feature_type, COUNT(*) AS count
        FROM features
        WHERE {' AND '.join(feature_clauses)}
        GROUP BY feature_type
        ORDER BY count DESC, feature_type
    """

    with connect() as conn:
        species_exists = conn.execute("SELECT 1 FROM species WHERE species_id = ?", (species_id,)).fetchone()
        if not species_exists:
            raise HTTPException(status_code=404, detail=f"species not found: {species_id}")
        chromosome_exists = conn.execute(
            "SELECT 1 FROM chromosomes WHERE species_id = ? AND seqid = ?",
            (species_id, seqid),
        ).fetchone()
        if not chromosome_exists:
            raise HTTPException(status_code=404, detail=f"chromosome not found: {species_id}/{seqid}")

        gene_rows = rows(conn, gene_sql, (species_id, seqid, start, end, gene_limit))
        feature_rows = rows(conn, feature_sql, [*feature_params, limit, offset])
        feature_counts = rows(conn, count_sql, feature_params)

    for item in gene_rows:
        item["attributes"] = parse_json(item.pop("attributes_json", None))
    for item in feature_rows:
        item["attributes"] = parse_json(item.pop("attributes_json", None))

    return {
        "species_id": species_id,
        "region": {"seqid": seqid, "start": start, "end": end},
        "feature_types": feature_types,
        "limit": limit,
        "offset": offset,
        "gene_limit": gene_limit,
        "gene_count": len(gene_rows),
        "feature_count": len(feature_rows),
        "feature_type_counts": feature_counts,
        "genes": gene_rows,
        "features": feature_rows,
    }
