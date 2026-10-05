#!/usr/bin/env python3
"""Create BLAST/DIAMOND indexes for one species and register them in SQLite."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import shutil
import sqlite3
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = ROOT / "data" / "antarctic_plant.db"
DEFAULT_ASSETS_DIR = ROOT / "data" / "species_assets"
ORTHOLOG_BIN = Path(os.environ.get("ANTARCTIC_ORTHOLOG_BIN", Path(sys.executable).parent))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, help="Species manifest YAML file.")
    parser.add_argument("--db", default=str(DEFAULT_DB), help="SQLite DB path.")
    parser.add_argument("--assets-dir", default=str(DEFAULT_ASSETS_DIR), help="Runtime asset output directory.")
    parser.add_argument("--threads", type=int, default=max(1, min(8, os.cpu_count() or 1)), help="Threads for index tools.")
    parser.add_argument("--force", action="store_true", help="Rebuild indexes even if outputs already exist.")
    parser.add_argument(
        "--indexes",
        default="diamond_protein,blast_protein,blast_genome",
        help="Comma-separated index types: diamond_protein, blast_protein, blast_genome.",
    )
    parser.add_argument(
        "--require",
        choices=["none", "diamond", "makeblastdb", "all"],
        default="none",
        help="Fail if required external tool is unavailable.",
    )
    return parser.parse_args()


def parse_simple_yaml(path: Path) -> Dict[str, Any]:
    root: Dict[str, Any] = {}
    stack: List[Tuple[int, Dict[str, Any]]] = [(-1, root)]
    for line_no, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        indent = len(raw) - len(raw.lstrip(" "))
        if indent % 2:
            raise ValueError(f"{path}:{line_no}: indentation must use two-space levels")
        text = raw.strip()
        if ":" not in text:
            raise ValueError(f"{path}:{line_no}: expected key: value")
        key, value = text.split(":", 1)
        key = key.strip()
        value = value.strip()
        while stack and indent <= stack[-1][0]:
            stack.pop()
        parent = stack[-1][1]
        if value == "":
            child: Dict[str, Any] = {}
            parent[key] = child
            stack.append((indent, child))
        else:
            parent[key] = value.strip("\"'")
    return root


def require(mapping: Dict[str, Any], dotted: str) -> Any:
    current: Any = mapping
    for part in dotted.split("."):
        if not isinstance(current, dict) or part not in current:
            raise ValueError(f"manifest missing required key: {dotted}")
        current = current[part]
    return current


def resolve_path(manifest_path: Path, value: str) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = (manifest_path.parent / path).resolve()
    return path


def tool_path(env_name: str, *candidates: str) -> Optional[str]:
    configured = os.environ.get(env_name)
    if configured:
        return configured
    for candidate in candidates:
        if Path(candidate).exists():
            return candidate
        found = shutil.which(candidate)
        if found:
            return found
    return None


def command_version(command: str, args: List[str]) -> str:
    try:
        proc = subprocess.run([command, *args], check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    except OSError:
        return "unknown"
    return " ".join((proc.stdout or "").strip().splitlines()[:2]) or "unknown"


def init_db(conn: sqlite3.Connection) -> None:
    conn.executescript(
        """
        PRAGMA foreign_keys = ON;

        CREATE TABLE IF NOT EXISTS search_assets (
            species_id TEXT NOT NULL,
            asset_type TEXT NOT NULL,
            source_fasta TEXT NOT NULL,
            index_prefix TEXT NOT NULL,
            tool TEXT NOT NULL,
            tool_version TEXT,
            created_at TEXT NOT NULL,
            status TEXT NOT NULL,
            metadata_json TEXT,
            PRIMARY KEY (species_id, asset_type),
            FOREIGN KEY (species_id) REFERENCES species(species_id) ON DELETE CASCADE
        );
        """
    )


def run_command(cmd: List[str]) -> None:
    subprocess.run(cmd, check=True)


def diamond_ready(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 0


def blast_ready(prefix: Path, dbtype: str) -> bool:
    suffixes = [".phr", ".pin", ".psq"] if dbtype == "prot" else [".nhr", ".nin", ".nsq"]
    return all(Path(str(prefix) + suffix).exists() for suffix in suffixes)


def register_asset(
    conn: sqlite3.Connection,
    species_id: str,
    asset_type: str,
    source_fasta: Path,
    index_prefix: Path,
    tool: str,
    tool_version: str,
    status: str,
    metadata: Dict[str, Any],
) -> None:
    conn.execute(
        """
        INSERT OR REPLACE INTO search_assets
        (species_id, asset_type, source_fasta, index_prefix, tool, tool_version, created_at, status, metadata_json)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            species_id,
            asset_type,
            str(source_fasta),
            str(index_prefix),
            tool,
            tool_version,
            dt.datetime.now(dt.timezone.utc).isoformat(),
            status,
            json.dumps(metadata, ensure_ascii=False, sort_keys=True),
        ),
    )


def build_indexes(args: argparse.Namespace) -> Dict[str, Any]:
    manifest_path = Path(args.manifest).resolve()
    manifest = parse_simple_yaml(manifest_path)
    species_id = str(require(manifest, "species_id"))
    genome_fasta = resolve_path(manifest_path, str(require(manifest, "genome.fasta")))
    protein_fasta = resolve_path(manifest_path, str(require(manifest, "protein.fasta")))
    if not genome_fasta.exists():
        raise FileNotFoundError(f"genome FASTA not found: {genome_fasta}")
    if not protein_fasta.exists():
        raise FileNotFoundError(f"protein FASTA not found: {protein_fasta}")

    requested = {item.strip() for item in args.indexes.split(",") if item.strip()}
    unknown = requested - {"diamond_protein", "blast_protein", "blast_genome"}
    if unknown:
        raise ValueError(f"unknown index type(s): {', '.join(sorted(unknown))}")

    diamond = tool_path("DIAMOND", str(ORTHOLOG_BIN / "diamond"), "diamond")
    makeblastdb = tool_path("MAKEBLASTDB", str(ORTHOLOG_BIN / "makeblastdb"), "makeblastdb")
    if args.require in {"diamond", "all"} and not diamond:
        raise FileNotFoundError("DIAMOND is required but was not found")
    if args.require in {"makeblastdb", "all"} and not makeblastdb:
        raise FileNotFoundError("makeblastdb is required but was not found")

    db_path = Path(args.db).resolve()
    if not db_path.exists():
        raise FileNotFoundError(f"SQLite DB not found: {db_path}; run build_species.py first")
    search_dir = Path(args.assets_dir).resolve() / species_id / "search"
    search_dir.mkdir(parents=True, exist_ok=True)

    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA foreign_keys = ON")
    init_db(conn)
    if not conn.execute("SELECT 1 FROM species WHERE species_id = ?", (species_id,)).fetchone():
        raise ValueError(f"species {species_id!r} is not in the database; run build_species.py first")

    summary: Dict[str, Any] = {"species_id": species_id, "created": [], "reused": [], "skipped": [], "warnings": []}
    try:
        if "diamond_protein" in requested:
            out_prefix = search_dir / "diamond" / "protein"
            out_dmnd = Path(str(out_prefix) + ".dmnd")
            out_prefix.parent.mkdir(parents=True, exist_ok=True)
            if diamond:
                if args.force or not diamond_ready(out_dmnd):
                    run_command([diamond, "makedb", "--in", str(protein_fasta), "-d", str(out_prefix), "--threads", str(args.threads)])
                    summary["created"].append("diamond_protein")
                else:
                    summary["reused"].append("diamond_protein")
                register_asset(
                    conn,
                    species_id,
                    "diamond_protein",
                    protein_fasta,
                    out_dmnd,
                    "diamond",
                    command_version(diamond, ["version"]),
                    "ready",
                    {"threads": args.threads},
                )
            else:
                summary["skipped"].append("diamond_protein")
                summary["warnings"].append("diamond not found")

        if "blast_protein" in requested:
            out_prefix = search_dir / "blast" / "protein"
            out_prefix.parent.mkdir(parents=True, exist_ok=True)
            if makeblastdb:
                if args.force or not blast_ready(out_prefix, "prot"):
                    run_command([makeblastdb, "-in", str(protein_fasta), "-dbtype", "prot", "-out", str(out_prefix)])
                    summary["created"].append("blast_protein")
                else:
                    summary["reused"].append("blast_protein")
                register_asset(
                    conn,
                    species_id,
                    "blast_protein",
                    protein_fasta,
                    out_prefix,
                    "makeblastdb",
                    command_version(makeblastdb, ["-version"]),
                    "ready",
                    {"dbtype": "prot"},
                )
            else:
                summary["skipped"].append("blast_protein")
                summary["warnings"].append("makeblastdb not found")

        if "blast_genome" in requested:
            out_prefix = search_dir / "blast" / "genome"
            out_prefix.parent.mkdir(parents=True, exist_ok=True)
            if makeblastdb:
                if args.force or not blast_ready(out_prefix, "nucl"):
                    run_command([makeblastdb, "-in", str(genome_fasta), "-dbtype", "nucl", "-out", str(out_prefix)])
                    summary["created"].append("blast_genome")
                else:
                    summary["reused"].append("blast_genome")
                register_asset(
                    conn,
                    species_id,
                    "blast_genome",
                    genome_fasta,
                    out_prefix,
                    "makeblastdb",
                    command_version(makeblastdb, ["-version"]),
                    "ready",
                    {"dbtype": "nucl"},
                )
            else:
                summary["skipped"].append("blast_genome")
                summary["warnings"].append("makeblastdb not found")

        conn.commit()
    except Exception:
        conn.rollback()
        raise
    finally:
        conn.close()
    return summary


def main() -> int:
    args = parse_args()
    try:
        summary = build_indexes(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
