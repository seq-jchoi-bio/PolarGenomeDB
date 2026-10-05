from __future__ import annotations

import os
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
DB_PATH = Path(os.environ.get("ANTARCTIC_DB", ROOT / "data" / "antarctic_plant.db")).resolve()
RUNTIME_DIR = Path(os.environ.get("ANTARCTIC_RUNTIME_DIR", ROOT / "data" / "runtime")).resolve()
JOB_DIR = Path(os.environ.get("ANTARCTIC_JOB_DIR", ROOT / "data" / "jobs")).resolve()
ORTHOLOG_BIN = Path(os.environ.get("ANTARCTIC_ORTHOLOG_BIN", Path(sys.executable).parent)).resolve()
FRONTEND_DIST = Path(os.environ.get("ANTARCTIC_FRONTEND_DIST", ROOT / "frontend" / "dist")).resolve()
API_LIMITS_PATH = Path(os.environ.get("ANTARCTIC_API_LIMITS", ROOT / "config" / "api_limits.yaml")).resolve()
PFAM_DB = Path(os.environ.get("ANTARCTIC_PFAM_DB", ROOT / "resources" / "pfam" / "Pfam-A.hmm")).resolve()
HMMSCAN = Path(os.environ.get("ANTARCTIC_HMMSCAN", ORTHOLOG_BIN / "hmmscan")).resolve()
GOFIGURE_PYTHON = Path(os.environ.get("ANTARCTIC_GOFIGURE_PYTHON", sys.executable)).expanduser().resolve()
GOFIGURE_SCRIPT = Path(os.environ.get("ANTARCTIC_GOFIGURE_SCRIPT", ROOT / "tools" / "gofigure.py")).expanduser().resolve()
