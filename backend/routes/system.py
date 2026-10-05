from __future__ import annotations

from fastapi import APIRouter

from ..api_limits import get_api_limits
from ..services.job_manager import job_manager


router = APIRouter()


@router.get("/system/limits")
def get_system_limits() -> dict:
    return get_api_limits().model_dump()


@router.get("/system/jobs")
def get_system_jobs() -> dict:
    return job_manager.status()
