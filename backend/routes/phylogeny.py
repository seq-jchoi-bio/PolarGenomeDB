from __future__ import annotations

from fastapi import APIRouter, HTTPException

from ..services.job_manager import job_manager
from ..services.phylogeny_jobs import PhylogenyRequest, create_phylogeny_job, execute_phylogeny_job, get_phylogeny_job


router = APIRouter()


@router.post("/phylogeny/tree")
def create_tree(request: PhylogenyRequest) -> dict:
    try:
        job = create_phylogeny_job(request)
        job_manager.submit("phylogeny", job["job_id"], execute_phylogeny_job, job["job_id"])
        return job
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from None


@router.get("/phylogeny/jobs/{job_id}")
def read_tree(job_id: str) -> dict:
    try:
        return get_phylogeny_job(job_id)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail="job not found") from None
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from None
