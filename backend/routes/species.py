from __future__ import annotations

from fastapi import APIRouter, HTTPException

from ..runtime import load_species_index, load_species_manifest


router = APIRouter()


@router.get("/species")
def list_species() -> dict:
    return load_species_index()


@router.get("/species/{species_id}")
def get_species(species_id: str) -> dict:
    try:
        return load_species_manifest(species_id)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"species not found: {species_id}") from None


@router.get("/species/{species_id}/chromosomes")
def get_chromosomes(species_id: str) -> dict:
    try:
        manifest = load_species_manifest(species_id)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"species not found: {species_id}") from None
    return {"species_id": species_id, "chromosomes": manifest["chromosomes"]}

