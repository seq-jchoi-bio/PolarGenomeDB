from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Any

import yaml
from pydantic import BaseModel, Field

from .config import API_LIMITS_PATH


class SearchLimits(BaseModel):
    max_query_characters: int = Field(default=100_000, ge=1)
    max_blast_sequences_per_request: int = Field(default=10, ge=1)
    max_diamond_sequences_per_request: int = Field(default=100, ge=1)
    max_hits: int = Field(default=100, ge=1)
    timeout_seconds: int = Field(default=120, ge=1)
    max_concurrent_jobs: int = Field(default=2, ge=1)
    job_retention_days: int = Field(default=7, ge=1)


class JobLimits(BaseModel):
    max_concurrent_jobs: int = Field(default=2, ge=1)
    job_retention_days: int = Field(default=7, ge=1)


class GeneLimits(BaseModel):
    default_page_size: int = Field(default=50, ge=1)
    max_page_size: int = Field(default=500, ge=1)
    default_suggest_size: int = Field(default=20, ge=1)
    max_suggest_size: int = Field(default=100, ge=1)


class FeatureLimits(BaseModel):
    default_page_size: int = Field(default=500, ge=1)
    max_page_size: int = Field(default=2000, ge=1)


class ExpressionLimits(BaseModel):
    max_genes_per_matrix_request: int = Field(default=100, ge=1)


class GoLimits(BaseModel):
    max_similarity_rank: int = Field(default=20, ge=1)
    max_enrichment_genes: int = Field(default=1_000, ge=1)
    max_enrichment_background_genes: int = Field(default=100_000, ge=1)


class PrimerLimits(BaseModel):
    max_region_length: int = Field(default=100_000, ge=100)
    max_flank: int = Field(default=10_000, ge=0)
    max_offtarget_product_size: int = Field(default=5_000, ge=50)
    primer3_timeout_seconds: int = Field(default=60, ge=1)
    specificity_timeout_seconds: int = Field(default=180, ge=1)


class SimilarLimits(BaseModel):
    max_selected_isoforms: int = Field(default=40, ge=1)
    candidates_per_isoform: int = Field(default=5, ge=1, le=20)
    raw_hits_per_isoform: int = Field(default=25, ge=1, le=100)
    timeout_seconds: int = Field(default=120, ge=1)


class ApiLimits(BaseModel):
    jobs: JobLimits = Field(default_factory=JobLimits)
    search: SearchLimits = Field(default_factory=SearchLimits)
    genes: GeneLimits = Field(default_factory=GeneLimits)
    features: FeatureLimits = Field(default_factory=FeatureLimits)
    expression: ExpressionLimits = Field(default_factory=ExpressionLimits)
    go: GoLimits = Field(default_factory=GoLimits)
    primer: PrimerLimits = Field(default_factory=PrimerLimits)
    similar: SimilarLimits = Field(default_factory=SimilarLimits)


def read_limits_yaml(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise ValueError(f"API limits config must be a mapping: {path}")
    return data


@lru_cache(maxsize=1)
def get_api_limits() -> ApiLimits:
    return ApiLimits(**read_limits_yaml(API_LIMITS_PATH))
