#!/usr/bin/env python3
"""Build RNA-seq expression assets and metadata for one species."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import os
import re
import shutil
import sqlite3
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = ROOT / "data" / "antarctic_plant.db"
DEFAULT_ASSETS_DIR = ROOT / "data" / "species_assets"
ORTHOLOG_BIN = Path(os.environ.get("ANTARCTIC_ORTHOLOG_BIN", Path(sys.executable).parent))
CONDA_BIN = Path(os.environ.get("ANTARCTIC_CONDA_BIN", Path(sys.executable).parent))
BAI_REFERENCE_LIMIT = 512 * 1024 * 1024


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, help="Species manifest YAML file.")
    parser.add_argument("--db", default=str(DEFAULT_DB), help="SQLite DB path.")
    parser.add_argument("--assets-dir", default=str(DEFAULT_ASSETS_DIR), help="Runtime asset output directory.")
    parser.add_argument("--threads", type=int, default=max(1, min(8, os.cpu_count() or 1)), help="Processors for bamCoverage.")
    parser.add_argument("--bin-size", type=int, default=25, help="bigWig bin size.")
    parser.add_argument("--normalization", default="CPM", help="bamCoverage normalization method.")
    parser.add_argument("--build-indexes", action="store_true", help="Create missing BAM indexes.")
    parser.add_argument("--build-bigwig", action="store_true", help="Create bigWig coverage tracks.")
    parser.add_argument("--load-counts", action="store_true", help="Load raw count matrix into SQLite.")
    parser.add_argument("--validate-only", action="store_true", help="Validate inputs without writing DB changes.")
    parser.add_argument("--limit-samples", type=int, help="Limit BAM/sample operations to the first N samples.")
    parser.add_argument("--force", action="store_true", help="Rebuild existing indexes/bigWigs and replace expression rows.")
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


def resolve_path(manifest_path: Path, value: Optional[str]) -> Optional[Path]:
    if not value:
        return None
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = (manifest_path.parent / path).resolve()
    return path


def as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


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


def run_command(cmd: List[str], log_path: Optional[Path] = None) -> None:
    if log_path:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with log_path.open("w", encoding="utf-8") as log:
            subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)
    else:
        subprocess.run(cmd, check=True)


def init_db(conn: sqlite3.Connection) -> None:
    conn.executescript(
        """
        PRAGMA foreign_keys = ON;

        CREATE TABLE IF NOT EXISTS samples (
            species_id TEXT NOT NULL,
            sample_id TEXT NOT NULL,
            display_name TEXT NOT NULL,
            tissue TEXT,
            condition TEXT,
            treatment TEXT,
            timepoint TEXT,
            replicate TEXT,
            source_accession TEXT,
            metadata_json TEXT,
            PRIMARY KEY (species_id, sample_id),
            FOREIGN KEY (species_id) REFERENCES species(species_id) ON DELETE CASCADE
        );

        CREATE TABLE IF NOT EXISTS raw_sample_metadata (
            species_id TEXT NOT NULL,
            sample_id TEXT NOT NULL,
            source_file TEXT NOT NULL,
            row_order INTEGER NOT NULL,
            metadata_json TEXT NOT NULL,
            PRIMARY KEY (species_id, sample_id, source_file),
            FOREIGN KEY (species_id) REFERENCES species(species_id) ON DELETE CASCADE
        );

        CREATE TABLE IF NOT EXISTS rnaseq_assets (
            species_id TEXT NOT NULL,
            sample_id TEXT NOT NULL,
            asset_type TEXT NOT NULL,
            path TEXT NOT NULL,
            status TEXT NOT NULL,
            metadata_json TEXT,
            PRIMARY KEY (species_id, sample_id, asset_type),
            FOREIGN KEY (species_id, sample_id) REFERENCES samples(species_id, sample_id) ON DELETE CASCADE
        );

        CREATE TABLE IF NOT EXISTS track_assets (
            species_id TEXT NOT NULL,
            sample_id TEXT NOT NULL,
            track_type TEXT NOT NULL,
            path TEXT NOT NULL,
            source_bam TEXT,
            bin_size INTEGER,
            normalization TEXT,
            created_at TEXT NOT NULL,
            status TEXT NOT NULL,
            metadata_json TEXT,
            PRIMARY KEY (species_id, sample_id, track_type, bin_size, normalization),
            FOREIGN KEY (species_id, sample_id) REFERENCES samples(species_id, sample_id) ON DELETE CASCADE
        );

        CREATE TABLE IF NOT EXISTS expression_counts (
            species_id TEXT NOT NULL,
            gene_id TEXT NOT NULL,
            sample_id TEXT NOT NULL,
            value_type TEXT NOT NULL,
            value REAL NOT NULL,
            source_file TEXT NOT NULL,
            PRIMARY KEY (species_id, gene_id, sample_id, value_type),
            FOREIGN KEY (species_id, gene_id) REFERENCES genes(species_id, gene_id) ON DELETE CASCADE,
            FOREIGN KEY (species_id, sample_id) REFERENCES samples(species_id, sample_id) ON DELETE CASCADE
        );

        CREATE INDEX IF NOT EXISTS idx_expression_gene ON expression_counts(species_id, gene_id);
        CREATE INDEX IF NOT EXISTS idx_expression_sample ON expression_counts(species_id, sample_id);
        """
    )


def sniff_delimiter(path: Path) -> str:
    sample = path.read_text(encoding="utf-8", errors="replace")[:4096]
    if "\t" in sample and sample.count("\t") >= sample.count(","):
        return "\t"
    return ","


def sample_id_from_bam_name(name: str) -> str:
    filename = Path(name).name
    return filename[:-4] if filename.endswith(".bam") else Path(filename).stem


def standard_sample_fields(row: Dict[str, str]) -> Dict[str, Optional[str]]:
    lower = {key.lower(): key for key in row}
    bam_value = row.get(lower.get("bam", ""), "")
    sample_id = (
        row.get(lower.get("sample", ""), "")
        or row.get(lower.get("sample_id", ""), "")
        or sample_id_from_bam_name(bam_value)
    )
    tissue = row.get(lower.get("tissue", ""), "") or None
    condition = row.get(lower.get("condition", ""), "") or row.get(lower.get("temperature", ""), "") or None
    replicate = row.get(lower.get("replicate", ""), "") or None
    source_accession = row.get(lower.get("accession", ""), "") or row.get(lower.get("sra", ""), "") or None
    return {
        "sample_id": sample_id,
        "bam_file": bam_value or None,
        "display_name": sample_id,
        "tissue": tissue,
        "condition": condition,
        "treatment": row.get(lower.get("treatment", ""), "") or None,
        "timepoint": row.get(lower.get("timepoint", ""), "") or None,
        "replicate": replicate,
        "source_accession": source_accession,
    }


def read_sample_metadata(path: Path) -> List[Dict[str, Any]]:
    delimiter = sniff_delimiter(path)
    rows: List[Dict[str, Any]] = []
    with path.open("rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"sample metadata has no header: {path}")
        for idx, row in enumerate(reader, start=1):
            clean_row = {str(key): (value or "").strip() for key, value in row.items()}
            fields = standard_sample_fields(clean_row)
            if not fields["sample_id"]:
                raise ValueError(f"cannot infer sample_id from metadata row {idx} in {path}")
            rows.append({"row_order": idx, "raw": clean_row, **fields})
    return rows


def bam_files_by_sample(bam_dir: Path) -> Dict[str, Path]:
    return {sample_id_from_bam_name(path.name): path for path in sorted(bam_dir.glob("*.bam"))}


def bam_files_by_name(bam_dir: Path) -> Dict[str, Path]:
    return {path.name: path for path in sorted(bam_dir.glob("*.bam"))}


def match_sample_to_bam(sample: Dict[str, Any], bam_by_sample: Dict[str, Path], bam_by_name: Dict[str, Path]) -> Optional[Path]:
    bam_file = sample.get("bam_file")
    if bam_file:
        bam_name = Path(str(bam_file)).name
        if bam_name in bam_by_name:
            return bam_by_name[bam_name]
        bam_stem = sample_id_from_bam_name(bam_name)
        if bam_stem in bam_by_sample:
            return bam_by_sample[bam_stem]
    return bam_by_sample.get(str(sample["sample_id"]))


def validate_bam_sample_links(samples: List[Dict[str, Any]], bam_dir: Path) -> Tuple[Dict[str, Path], List[str]]:
    bam_by_sample = bam_files_by_sample(bam_dir)
    bam_by_name = bam_files_by_name(bam_dir)
    matched: Dict[str, Path] = {}
    errors: List[str] = []
    for sample in samples:
        sample_id = str(sample["sample_id"])
        bam = match_sample_to_bam(sample, bam_by_sample, bam_by_name)
        if not bam:
            bam_hint = f" bam={sample['bam_file']}" if sample.get("bam_file") else ""
            errors.append(f"metadata sample {sample_id}{bam_hint} does not match any BAM in {bam_dir}")
            continue
        matched[sample_id] = bam
    return matched, errors


def max_reference_length(conn: sqlite3.Connection, species_id: str) -> int:
    row = conn.execute("SELECT MAX(length) FROM chromosomes WHERE species_id = ?", (species_id,)).fetchone()
    return int(row[0] or 0)


def desired_index_type(max_len: int) -> str:
    return "csi" if max_len >= BAI_REFERENCE_LIMIT else "bai"


def index_path_for_bam(bam: Path, index_type: str) -> Path:
    return Path(str(bam) + f".{index_type}")


def ensure_bam_index(samtools: str, bam: Path, index_type: str, force: bool) -> Tuple[Path, str]:
    index_path = index_path_for_bam(bam, index_type)
    if force and index_path.exists():
        index_path.unlink()
    if index_path.exists() and index_path.stat().st_mtime >= bam.stat().st_mtime:
        return index_path, "reused"
    cmd = [samtools, "index", "-c" if index_type == "csi" else "-b", str(bam)]
    run_command(cmd)
    return index_path, "created"


def build_bigwig(
    bam_coverage: str,
    bam: Path,
    output: Path,
    threads: int,
    bin_size: int,
    normalization: str,
    force: bool,
    log_path: Path,
) -> str:
    if force and output.exists():
        output.unlink()
    if output.exists() and output.stat().st_size > 0:
        return "reused"
    output.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        bam_coverage,
        "-b",
        str(bam),
        "-o",
        str(output),
        "--binSize",
        str(bin_size),
        "--normalizeUsing",
        normalization,
        "--numberOfProcessors",
        str(threads),
    ]
    env = os.environ.copy()
    env.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-cache")
    Path(env["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as log:
        subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT, env=env)
    return "created"


def read_count_matrix(path: Path, has_header: bool, sample_ids: List[str]) -> Tuple[List[str], List[Tuple[str, List[int]]]]:
    rows: List[Tuple[str, List[int]]] = []
    matrix_samples: Optional[List[str]] = None
    seen_genes = set()
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        for line_no, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if has_header and line_no == 1:
                if len(parts) < 2:
                    raise ValueError("count matrix header must include at least one sample column")
                matrix_samples = parts[1:]
                continue
            if len(parts) < 2:
                raise ValueError(f"{path}:{line_no}: expected gene_id plus one or more count columns")
            gene_id = parts[0].strip()
            if not gene_id or gene_id.startswith("__"):
                continue
            if gene_id in seen_genes:
                raise ValueError(f"{path}:{line_no}: duplicate gene_id {gene_id}")
            seen_genes.add(gene_id)
            values: List[int] = []
            for value in parts[1:]:
                value = value.strip()
                if value == "" or value.lower() in {"na", "nan", "null"}:
                    raise ValueError(f"{path}:{line_no}: NA/blank count is not allowed")
                if not re.fullmatch(r"[0-9]+", value):
                    raise ValueError(f"{path}:{line_no}: count must be a non-negative integer: {value!r}")
                values.append(int(value))
            rows.append((gene_id, values))
    if matrix_samples is None:
        if not rows:
            raise ValueError(f"count matrix has no rows: {path}")
        count_cols = len(rows[0][1])
        if count_cols > len(sample_ids):
            raise ValueError(f"count matrix has {count_cols} columns but metadata has only {len(sample_ids)} samples")
        matrix_samples = sample_ids[:count_cols]
    expected_cols = len(matrix_samples)
    for gene_id, values in rows:
        if len(values) != expected_cols:
            raise ValueError(f"gene {gene_id} has {len(values)} counts; expected {expected_cols}")
    return matrix_samples, rows


def load_counts(
    conn: sqlite3.Connection,
    species_id: str,
    count_path: Path,
    has_header: bool,
    value_type: str,
    sample_ids: List[str],
    force: bool,
) -> Dict[str, Any]:
    matrix_samples, rows = read_count_matrix(count_path, has_header, sample_ids)
    known_samples = {
        row[0]
        for row in conn.execute("SELECT sample_id FROM samples WHERE species_id = ?", (species_id,)).fetchall()
    }
    missing_samples = [sample_id for sample_id in matrix_samples if sample_id not in known_samples]
    if missing_samples:
        raise ValueError(f"count matrix sample(s) missing from samples table: {', '.join(missing_samples)}")
    known_genes = {
        row[0]
        for row in conn.execute("SELECT gene_id FROM genes WHERE species_id = ?", (species_id,)).fetchall()
    }
    missing_genes = [gene_id for gene_id, values in rows if gene_id not in known_genes]
    if missing_genes:
        preview = ", ".join(missing_genes[:10])
        raise ValueError(f"{len(missing_genes)} count matrix genes are not in annotation DB: {preview}")
    if force:
        conn.execute(
            "DELETE FROM expression_counts WHERE species_id = ? AND value_type = ? AND source_file = ?",
            (species_id, value_type, str(count_path)),
        )
    payload = []
    for gene_id, values in rows:
        for sample_id, value in zip(matrix_samples, values):
            payload.append((species_id, gene_id, sample_id, value_type, float(value), str(count_path)))
    conn.executemany(
        """
        INSERT OR REPLACE INTO expression_counts
        (species_id, gene_id, sample_id, value_type, value, source_file)
        VALUES (?, ?, ?, ?, ?, ?)
        """,
        payload,
    )
    library_sizes = {sample_id: 0 for sample_id in matrix_samples}
    for gene_id, values in rows:
        for sample_id, value in zip(matrix_samples, values):
            library_sizes[sample_id] += value
    cpm_payload = []
    for gene_id, values in rows:
        for sample_id, value in zip(matrix_samples, values):
            lib_size = library_sizes[sample_id]
            cpm = (value / lib_size * 1_000_000.0) if lib_size else 0.0
            cpm_payload.append((species_id, gene_id, sample_id, "cpm", cpm, str(count_path)))
    conn.executemany(
        """
        INSERT OR REPLACE INTO expression_counts
        (species_id, gene_id, sample_id, value_type, value, source_file)
        VALUES (?, ?, ?, ?, ?, ?)
        """,
        cpm_payload,
    )
    return {"count_rows": len(rows), "count_samples": len(matrix_samples), "library_sizes": library_sizes}


def build_expression(args: argparse.Namespace) -> Dict[str, Any]:
    manifest_path = Path(args.manifest).resolve()
    manifest = parse_simple_yaml(manifest_path)
    species_id = str(require(manifest, "species_id"))
    expression = manifest.get("expression")
    if expression is None:
        return {
            "species_id": species_id,
            "status": "skipped",
            "reason": "manifest has no expression section",
            "warnings": [],
        }
    if not isinstance(expression, dict):
        raise ValueError("manifest expression section must be a mapping")
    sample_metadata = resolve_path(manifest_path, expression.get("sample_metadata"))
    bam_dir = resolve_path(manifest_path, expression.get("bam_dir"))
    count_matrix = resolve_path(manifest_path, expression.get("count_matrix"))

    needs_samples = args.validate_only or args.load_counts or args.build_indexes or args.build_bigwig
    needs_bam = args.build_indexes or args.build_bigwig
    if not (args.validate_only or args.load_counts or args.build_indexes or args.build_bigwig):
        needs_samples = bool(sample_metadata)
        needs_bam = bool(bam_dir)

    if needs_samples and (not sample_metadata or not sample_metadata.exists()):
        raise FileNotFoundError(f"sample metadata not found: {sample_metadata}")
    if needs_bam and (not bam_dir or not bam_dir.exists()):
        raise FileNotFoundError(f"BAM directory not found: {bam_dir}")
    if args.load_counts and (not count_matrix or not count_matrix.exists()):
        raise FileNotFoundError(f"count matrix not found: {count_matrix}")
    if not sample_metadata or not sample_metadata.exists():
        sample_metadata = None
    if not bam_dir or not bam_dir.exists():
        bam_dir = None

    db_path = Path(args.db).resolve()
    if not db_path.exists():
        raise FileNotFoundError(f"SQLite DB not found: {db_path}")
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA foreign_keys = ON")
    init_db(conn)
    if not conn.execute("SELECT 1 FROM species WHERE species_id = ?", (species_id,)).fetchone():
        raise ValueError(f"species {species_id!r} is not in DB; run build_species.py first")

    samples = read_sample_metadata(sample_metadata) if sample_metadata else []
    if args.limit_samples:
        samples = samples[: args.limit_samples]
    bam_by_sample: Dict[str, Path] = {}
    if bam_dir:
        bam_by_sample, bam_errors = validate_bam_sample_links(samples, bam_dir)
        if bam_errors and needs_bam:
            raise ValueError("; ".join(bam_errors))
    max_len = max_reference_length(conn, species_id)
    index_type = desired_index_type(max_len)
    samtools = tool_path("SAMTOOLS", "/usr/local/bin/samtools", str(ORTHOLOG_BIN / "samtools"), "samtools")
    bam_coverage = tool_path("BAMCOVERAGE", str(CONDA_BIN / "bamCoverage"), "bamCoverage")
    summary: Dict[str, Any] = {
        "species_id": species_id,
        "samples": len(samples),
        "bam_dir": str(bam_dir) if bam_dir else None,
        "sample_metadata": str(sample_metadata) if sample_metadata else None,
        "bam_index_type": index_type,
        "max_reference_length": max_len,
        "indexes": {"created": 0, "reused": 0, "missing": 0},
        "bigwig": {"created": 0, "reused": 0, "skipped": 0},
        "counts": {},
        "warnings": [],
    }

    if args.validate_only:
        if bam_dir:
            for error in bam_errors:
                summary["warnings"].append(error)
        elif args.build_indexes or args.build_bigwig:
            summary["warnings"].append("BAM directory is not configured")
        conn.close()
        return summary

    try:
        conn.execute("BEGIN")
        for sample in samples:
            sample_id = sample["sample_id"]
            metadata_json = json.dumps(sample["raw"], ensure_ascii=False, sort_keys=True)
            conn.execute(
                """
                INSERT INTO samples
                (species_id, sample_id, display_name, tissue, condition, treatment, timepoint, replicate, source_accession, metadata_json)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                ON CONFLICT(species_id, sample_id) DO UPDATE SET
                    display_name = excluded.display_name,
                    tissue = excluded.tissue,
                    condition = excluded.condition,
                    treatment = excluded.treatment,
                    timepoint = excluded.timepoint,
                    replicate = excluded.replicate,
                    source_accession = excluded.source_accession,
                    metadata_json = excluded.metadata_json
                """,
                (
                    species_id,
                    sample_id,
                    sample["display_name"],
                    sample["tissue"],
                    sample["condition"],
                    sample["treatment"],
                    sample["timepoint"],
                    sample["replicate"],
                    sample["source_accession"],
                    metadata_json,
                ),
            )
            conn.execute(
                """
                INSERT OR REPLACE INTO raw_sample_metadata
                (species_id, sample_id, source_file, row_order, metadata_json)
                VALUES (?, ?, ?, ?, ?)
                """,
                (species_id, sample_id, str(sample_metadata), sample["row_order"], metadata_json),
            )

            if not bam_dir:
                summary["bigwig"]["skipped"] += 1
                continue
            bam = bam_by_sample.get(sample_id)
            if not bam:
                summary["warnings"].append(f"BAM missing for sample {sample_id}")
                continue
            index_path = index_path_for_bam(bam, index_type)
            index_status = "missing"
            if args.build_indexes:
                if not samtools:
                    raise FileNotFoundError("samtools not found")
                index_path, index_status = ensure_bam_index(samtools, bam, index_type, args.force)
            elif index_path.exists():
                index_status = "reused"
            else:
                summary["indexes"]["missing"] += 1
            if index_status in {"created", "reused"}:
                summary["indexes"][index_status] += 1
            conn.execute(
                """
                INSERT OR REPLACE INTO rnaseq_assets
                (species_id, sample_id, asset_type, path, status, metadata_json)
                VALUES (?, ?, ?, ?, ?, ?)
                """,
                (species_id, sample_id, "bam", str(bam), "ready", json.dumps({"index": str(index_path), "index_status": index_status}, sort_keys=True)),
            )
            if index_status in {"created", "reused"}:
                conn.execute(
                    """
                    INSERT OR REPLACE INTO rnaseq_assets
                    (species_id, sample_id, asset_type, path, status, metadata_json)
                    VALUES (?, ?, ?, ?, ?, ?)
                    """,
                    (species_id, sample_id, index_type, str(index_path), "ready", json.dumps({"source_bam": str(bam)}, sort_keys=True)),
                )

            if args.build_bigwig:
                if not bam_coverage:
                    raise FileNotFoundError("bamCoverage not found")
                if not index_path.exists():
                    raise FileNotFoundError(f"BAM index required for bamCoverage: {index_path}")
                output = Path(args.assets_dir).resolve() / species_id / "expression" / "bigwig" / f"{sample_id}.bin{args.bin_size}.{args.normalization.lower()}.bw"
                log_path = ROOT / "reports" / f"{species_id}_{sample_id}.bamCoverage_bin{args.bin_size}.log"
                status = build_bigwig(bam_coverage, bam, output, args.threads, args.bin_size, args.normalization, args.force, log_path)
                summary["bigwig"][status] += 1
                conn.execute(
                    """
                    INSERT OR REPLACE INTO track_assets
                    (species_id, sample_id, track_type, path, source_bam, bin_size, normalization, created_at, status, metadata_json)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        species_id,
                        sample_id,
                        "bigwig",
                        str(output),
                        str(bam),
                        args.bin_size,
                        args.normalization,
                        dt.datetime.now(dt.timezone.utc).isoformat(),
                        "ready",
                        json.dumps({"log": str(log_path)}, sort_keys=True),
                    ),
                )
            else:
                summary["bigwig"]["skipped"] += 1

        if args.load_counts:
            sample_ids = [sample["sample_id"] for sample in read_sample_metadata(sample_metadata)]
            has_header = as_bool(expression.get("count_has_header", False))
            value_type = str(expression.get("count_value_type", "raw_count"))
            summary["counts"] = load_counts(conn, species_id, count_matrix, has_header, value_type, sample_ids, args.force)

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
        summary = build_expression(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
