from __future__ import annotations

from pathlib import Path

from fastapi import APIRouter, HTTPException

from ..config import ROOT
from ..runtime import load_runtime_json, load_species_manifest


router = APIRouter()


@router.get("/species/{species_id}/tracks")
def get_tracks(species_id: str) -> dict:
    try:
        return load_runtime_json(f"{species_id}/genome_browser_tracks.json")
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"species not found: {species_id}") from None


@router.get("/species/{species_id}/samples")
def get_samples(species_id: str) -> dict:
    try:
        manifest = load_species_manifest(species_id)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"species not found: {species_id}") from None
    return {"species_id": species_id, "samples": manifest["samples"]}


@router.get("/species/{species_id}/jbrowse-config")
def get_jbrowse_config(species_id: str) -> dict:
    path = ROOT / "data" / "jbrowse" / species_id / "config.json"
    if not path.exists():
        raise HTTPException(status_code=404, detail=f"JBrowse config not found: {species_id}")
    return {
        "species_id": species_id,
        "config_url": f"/files/data/jbrowse/{species_id}/config.json",
        "config": load_runtime_json(f"../jbrowse/{species_id}/config.json"),
    }


@router.get("/species/{species_id}/viewer-location/{gene_id}")
def get_viewer_location(species_id: str, gene_id: str) -> dict:
    try:
        manifest = load_species_manifest(species_id)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"species not found: {species_id}") from None

    # Gene coordinates are read from the runtime-independent API DB endpoint contract indirectly
    # by querying the SQLite-backed gene API would cause an HTTP round trip, so import locally.
    from ..db import connect, row

    with connect() as conn:
        gene = row(
            conn,
            "SELECT gene_id, seqid, start, end, strand FROM genes WHERE species_id = ? AND gene_id = ?",
            (species_id, gene_id),
        )
    if not gene:
        raise HTTPException(status_code=404, detail=f"gene not found: {species_id}/{gene_id}")

    padding = max(1000, int((gene["end"] - gene["start"] + 1) * 0.1))
    start = max(1, gene["start"] - padding)
    end = gene["end"] + padding
    return {
        "species_id": species_id,
        "gene_id": gene_id,
        "seqid": gene["seqid"],
        "start": start,
        "end": end,
        "strand": gene["strand"],
        "location": f"{gene['seqid']}:{start}-{end}",
        "config_url": f"/files/data/jbrowse/{species_id}/config.json",
        "assembly": manifest["species"]["display_name"],
    }
