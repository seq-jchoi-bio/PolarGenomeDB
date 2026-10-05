from __future__ import annotations

import asyncio
from contextlib import asynccontextmanager
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from .frontend import FrontendFiles

from .config import FRONTEND_DIST
from .services.job_manager import job_manager
from .routes import browser_files
from .routes import assets, domains, expression, features, genes, go, phylogeny, primer, search, similarity, species, system, tracks


@asynccontextmanager
async def lifespan(app: FastAPI):
    job_manager.start()
    try:
        yield
    finally:
        await asyncio.to_thread(job_manager.shutdown)


app = FastAPI(title="Polar Genome DB API", version="0.1.0", lifespan=lifespan)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=False,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/api/health")
def health() -> dict[str, str]:
    return {"status": "ok"}


app.include_router(species.router, prefix="/api", tags=["species"])
app.include_router(genes.router, prefix="/api", tags=["genes"])
app.include_router(features.router, prefix="/api", tags=["features"])
app.include_router(expression.router, prefix="/api", tags=["expression"])
app.include_router(domains.router, prefix="/api", tags=["domains"])
app.include_router(tracks.router, prefix="/api", tags=["tracks"])
app.include_router(assets.router, prefix="/api", tags=["assets"])
app.include_router(search.router, prefix="/api", tags=["search"])
app.include_router(similarity.router, prefix="/api", tags=["similarity"])
app.include_router(go.router, prefix="/api", tags=["go"])
app.include_router(primer.router, prefix="/api", tags=["primer"])
app.include_router(phylogeny.router, prefix="/api", tags=["phylogeny"])
app.include_router(system.router, prefix="/api", tags=["system"])
app.include_router(browser_files.router)

if FRONTEND_DIST.exists():
    app.mount("/", FrontendFiles(directory=str(FRONTEND_DIST), html=True), name="frontend")
