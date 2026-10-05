from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .config import RUNTIME_DIR


def load_runtime_json(relative_path: str, runtime_dir: Path = RUNTIME_DIR) -> dict[str, Any]:
    path = (runtime_dir / relative_path).resolve()
    if not path.exists():
        raise FileNotFoundError(path)
    return json.loads(path.read_text(encoding="utf-8"))


def load_species_index() -> dict[str, Any]:
    return load_runtime_json("species.json")


def load_species_manifest(species_id: str) -> dict[str, Any]:
    return load_runtime_json(f"{species_id}/manifest.json")

