"""Expose only explicitly configured JBrowse assets, never the project directory."""
from __future__ import annotations

import json
import re
from pathlib import Path
from urllib.parse import unquote, urlsplit

from fastapi import APIRouter, HTTPException, Request
from fastapi.responses import FileResponse, Response, StreamingResponse

from ..config import ROOT
from ..runtime import load_species_index

router = APIRouter()
MAX_RANGE_BYTES = 8 * 1024 * 1024
INDEX_SUFFIXES = (".fai", ".gzi", ".tbi", ".csi")
DATA_SUFFIXES = (".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".gff.gz", ".gff3.gz", ".gtf.gz", ".bw", ".bigwig")


def safe_path(relative: str) -> Path | None:
    if not relative or relative.startswith("/") or "\\" in relative or ".." in relative.split("/"):
        return None
    path = (ROOT / relative).resolve()
    return path if path.is_relative_to(ROOT.resolve()) else None


def browser_allowlist() -> dict[str, tuple[Path, bool]]:
    allowed: dict[str, tuple[Path, bool]] = {}

    def add_locations(value):
        if isinstance(value, dict):
            uri = value.get("uri")
            if isinstance(uri, str):
                parsed = urlsplit(uri)
                if not parsed.scheme and not parsed.netloc and parsed.path.startswith("/files/"):
                    relative = unquote(parsed.path[len("/files/"):])
                    path = safe_path(relative)
                    if path and relative.lower().endswith(INDEX_SUFFIXES + DATA_SUFFIXES):
                        allowed[relative] = (path, not relative.lower().endswith(INDEX_SUFFIXES))
            for child in value.values():
                add_locations(child)
        elif isinstance(value, list):
            for child in value:
                add_locations(child)

    for species in load_species_index().get("species", []):
        species_id = species.get("species_id", "")
        if not re.fullmatch(r"[A-Za-z0-9_-]+", species_id):
            continue
        relative = f"data/jbrowse/{species_id}/config.json"
        path = safe_path(relative)
        if path is None or not path.is_file():
            continue
        config = json.loads(path.read_text())
        allowed[relative] = (path, False)
        # Only adapter locations are public; metadata paths are not authorization.
        for assembly in config.get("assemblies", []):
            add_locations(assembly.get("sequence", {}).get("adapter", {}))
        for track in config.get("tracks", []):
            add_locations(track.get("adapter", {}))
    return allowed


def range_bounds(header: str, size: int) -> tuple[int, int]:
    match = re.fullmatch(r"bytes=(\d+)-(\d+)", header)
    if not match:
        raise HTTPException(416, "A single bounded byte range is required")
    start, end = map(int, match.groups())
    if start >= size or end < start or end - start + 1 > MAX_RANGE_BYTES:
        raise HTTPException(416, "Requested byte range exceeds viewer limits")
    return start, min(end, size - 1)


@router.api_route("/files/{relative:path}", methods=["GET", "HEAD"])
def browser_file(relative: str, request: Request):
    record = browser_allowlist().get(relative)
    if record is None or not record[0].is_file():
        raise HTTPException(404, "File is not available")
    path, range_only = record
    if not range_only:
        return FileResponse(path)
    size = path.stat().st_size
    headers = {"Accept-Ranges": "bytes", "Cache-Control": "private, no-store"}
    if request.method == "HEAD":
        return Response(headers={**headers, "Content-Length": str(size)}, media_type="application/octet-stream")
    header = request.headers.get("range")
    if not header:
        raise HTTPException(403, "Bulk download is not available; viewer byte ranges only")
    start, end = range_bounds(header, size)

    def chunks():
        with path.open("rb") as handle:
            handle.seek(start)
            remaining = end - start + 1
            while remaining:
                block = handle.read(min(256 * 1024, remaining))
                if not block:
                    break
                remaining -= len(block)
                yield block

    # Do not honor If-Range by falling back to a full-file response.
    return StreamingResponse(chunks(), status_code=206, media_type="application/octet-stream", headers={
        **headers, "Content-Length": str(end - start + 1), "Content-Range": f"bytes {start}-{end}/{size}",
    })
