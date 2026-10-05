from __future__ import annotations

from fastapi import APIRouter, HTTPException

from ..services.primer_jobs import PrimerRequest, create_primer_job, execute_primer_job, get_primer_job
from ..services.job_manager import job_manager


router = APIRouter()


@router.post("/primer/design")
def design_primer(request: PrimerRequest) -> dict:
    try:
        job = create_primer_job(request)
        job_manager.submit("primer", job["job_id"], execute_primer_job, job["job_id"])
        return job
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"species not found: {request.species_id}") from None
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from None


@router.get("/primer/jobs/{job_id}")
def read_primer_job(job_id: str) -> dict:
    try:
        return get_primer_job(job_id)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"job not found: {job_id}") from None
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from None
