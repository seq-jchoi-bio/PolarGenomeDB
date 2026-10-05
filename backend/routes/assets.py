from __future__ import annotations

from pathlib import Path
from typing import Any

from fastapi import APIRouter, HTTPException

from ..config import ROOT
from ..db import connect, parse_json, rows


router = APIRouter()


def resolve_asset_path(path: str) -> Path:
    value = Path(path)
    return value if value.is_absolute() else ROOT / value


def asset_record(path: str, extra: dict[str, Any] | None = None) -> dict[str, Any]:
    resolved = resolve_asset_path(path)
    stat = resolved.stat() if resolved.is_file() else None
    payload: dict[str, Any] = {
        "path": path,
        "exists": bool(stat and stat.st_size > 0),
        "size_bytes": stat.st_size if stat else None,
    }
    if extra:
        payload.update(extra)
    return payload


def prefix_asset_record(path: str, extra: dict[str, Any] | None = None) -> dict[str, Any]:
    resolved = resolve_asset_path(path)
    parent = resolved.parent
    pattern = f"{resolved.name}*"
    files = [item for item in sorted(parent.glob(pattern)) if item.is_file() and item.stat().st_size > 0] if parent.exists() else []
    kind = (extra or {}).get("asset_type")
    present = {str(item) for item in files}
    complete = False
    if kind == "bowtie2_genome":
        complete = any(all(f"{resolved}.{part}.{extension}" in present
                           for part in ("1", "2", "3", "4", "rev.1", "rev.2"))
                       for extension in ("bt2", "bt2l"))
    elif kind == "diamond_protein":
        complete = str(resolved if resolved.suffix == ".dmnd" else Path(str(resolved) + ".dmnd")) in present
    elif kind in ("blast_genome", "blast_protein"):
        letter = "n" if kind == "blast_genome" else "p"
        volumes = [str(item)[:-4] for item in files if item.name.endswith(f".{letter}sq")]
        complete = bool(volumes) and all(f"{volume}.{letter}{suffix}" in present
                                         for volume in volumes for suffix in ("hr", "in", "sq"))
    total_size = sum(item.stat().st_size for item in files if item.is_file())
    payload: dict[str, Any] = {
        "path": path,
        "exists": complete,
        "size_bytes": total_size if files else None,
        "file_count": len(files),
    }
    if extra:
        payload.update(extra)
    return payload


@router.get("/species/{species_id}/assets/status")
def get_asset_status(species_id: str) -> dict:
    with connect() as conn:
        species = conn.execute(
            "SELECT species_id, scientific_name, display_name FROM species WHERE species_id = ?",
            (species_id,),
        ).fetchone()
        if not species:
            raise HTTPException(status_code=404, detail=f"species not found: {species_id}")

        genome_rows = rows(conn, "SELECT asset_type, path FROM genome_assets WHERE species_id = ? ORDER BY asset_type", (species_id,))
        annotation_rows = rows(conn, "SELECT format, path FROM annotation_assets WHERE species_id = ? ORDER BY format, path", (species_id,))
        search_rows = rows(
            conn,
            """
            SELECT asset_type, source_fasta, index_prefix, tool, tool_version, created_at, status, metadata_json
            FROM search_assets
            WHERE species_id = ?
            ORDER BY asset_type
            """,
            (species_id,),
        )
        track_rows = rows(
            conn,
            """
            SELECT sample_id, track_type, path, source_bam, bin_size, normalization, created_at, status, metadata_json
            FROM track_assets
            WHERE species_id = ?
            ORDER BY sample_id, track_type, bin_size, normalization
            """,
            (species_id,),
        )
        expression_summary = rows(
            conn,
            """
            SELECT value_type, COUNT(*) AS value_count, COUNT(DISTINCT gene_id) AS genes, COUNT(DISTINCT sample_id) AS samples
            FROM expression_counts
            WHERE species_id = ?
            GROUP BY value_type
            ORDER BY value_type
            """,
            (species_id,),
        )

    genome = [
        asset_record(item["path"], {"asset_type": item["asset_type"]})
        for item in genome_rows
    ]
    annotations = [
        asset_record(item["path"], {"format": item["format"]})
        for item in annotation_rows
    ]
    search_assets = [
        prefix_asset_record(
            item["index_prefix"],
            {
                "asset_type": item["asset_type"],
                "source_fasta": item["source_fasta"],
                "tool": item["tool"],
                "tool_version": item["tool_version"],
                "created_at": item["created_at"],
                "status": item["status"],
                "metadata": parse_json(item.pop("metadata_json", None)),
            },
        )
        for item in search_rows
    ]
    tracks = [
        asset_record(
            item["path"],
            {
                "sample_id": item["sample_id"],
                "track_type": item["track_type"],
                "source_bam": item["source_bam"],
                "bin_size": item["bin_size"],
                "normalization": item["normalization"],
                "created_at": item["created_at"],
                "status": item["status"],
                "metadata": parse_json(item.pop("metadata_json", None)),
            },
        )
        for item in track_rows
    ]
    groups = {
        "genome": genome,
        "annotation": annotations,
        "search": search_assets,
        "tracks": tracks,
        "expression": expression_summary,
    }
    return {
        "species_id": species_id,
        "species": dict(species),
        "ready": {"fasta", "fai"}.issubset({item["asset_type"] for item in genome})
        and "gff3_tabix" in {item["format"] for item in annotations}
        and bool({"gff3_tabix_tbi", "gff3_tabix_csi"} & {item["format"] for item in annotations})
        and {"blast_genome", "blast_protein", "diamond_protein"}.issubset({item["asset_type"] for item in search_assets})
        and all(item["exists"] for group in (genome, annotations, search_assets, tracks) for item in group)
        and all(item["status"] == "ready" for item in search_assets + tracks),
        "counts": {key: len(value) for key, value in groups.items() if isinstance(value, list)},
        "assets": groups,
    }
