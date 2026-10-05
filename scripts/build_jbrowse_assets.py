#!/usr/bin/env python3
"""Prepare tabix-indexed annotation assets for JBrowse."""

from __future__ import annotations

import argparse
import json
import shutil
import sqlite3
import subprocess
from pathlib import Path
import os
import sys
from typing import Any, Dict, Iterable, List, Optional, Tuple


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = ROOT / "data" / "antarctic_plant.db"
DEFAULT_ASSETS = ROOT / "data" / "species_assets"
ORTHOLOG_BIN = Path(os.environ.get("ANTARCTIC_ORTHOLOG_BIN", Path(sys.executable).parent))
TBI_REFERENCE_LIMIT = 512 * 1024 * 1024


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", default=str(DEFAULT_DB), help="SQLite DB path.")
    parser.add_argument("--assets-dir", default=str(DEFAULT_ASSETS), help="Species asset directory.")
    parser.add_argument("--species", help="Optional single species_id.")
    parser.add_argument("--force", action="store_true", help="Rebuild existing compressed annotation and tabix index.")
    parser.add_argument("--keep-sorted", action="store_true", help="Keep the intermediate sorted annotation text file.")
    return parser.parse_args()


def rows(conn: sqlite3.Connection, query: str, params: Iterable[Any] = ()) -> List[Dict[str, Any]]:
    cursor = conn.execute(query, tuple(params))
    columns = [column[0] for column in cursor.description]
    return [dict(zip(columns, row)) for row in cursor.fetchall()]


def tool_path(name: str) -> str:
    candidates = [str(ORTHOLOG_BIN / name), name]
    for candidate in candidates:
        found = candidate if Path(candidate).exists() else shutil.which(candidate)
        if found:
            return str(found)
    raise FileNotFoundError(f"{name} not found")


def output_paths(species_id: str, annotation_format: str, source: Path, assets_dir: Path, index_type: str) -> Tuple[Path, Path, Path]:
    suffix = ".gff3"
    out_dir = assets_dir / species_id / "annotation"
    base = out_dir / f"{species_id}.annotation.sorted{suffix}"
    return base, Path(str(base) + ".gz"), Path(str(base) + f".gz.{index_type.lower()}")


def parse_gtf_attributes(value: str) -> Dict[str, str]:
    attributes: Dict[str, str] = {}
    for item in value.strip().rstrip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if " " in item:
            key, raw_value = item.split(" ", 1)
            attributes[key] = raw_value.strip().strip('"')
    return attributes


def gff3_escape(value: str) -> str:
    return value.replace("%", "%25").replace(";", "%3B").replace(",", "%2C").replace("=", "%3D").replace("&", "%26")


def gff3_attributes(items: Dict[str, str]) -> str:
    return ";".join(f"{key}={gff3_escape(value)}" for key, value in items.items() if value)


def gtf_to_gff3_line(parts: List[str]) -> Optional[str]:
    if len(parts) < 9:
        return None
    attrs = parse_gtf_attributes(parts[8])
    gene_id = attrs.get("gene_id")
    transcript_id = attrs.get("transcript_id")
    feature_type = parts[2]

    next_attrs: Dict[str, str] = {}
    if feature_type == "gene":
        if not gene_id:
            return None
        next_attrs["ID"] = gene_id
        next_attrs["Name"] = gene_id
    elif feature_type in {"transcript", "mRNA"}:
        if not transcript_id:
            return None
        parts[2] = "mRNA"
        next_attrs["ID"] = transcript_id
        if gene_id:
            next_attrs["Parent"] = gene_id
        next_attrs["Name"] = transcript_id
    else:
        parent = transcript_id or gene_id
        if parent:
            next_attrs["Parent"] = parent
        if gene_id:
            next_attrs["gene_id"] = gene_id
        if transcript_id:
            next_attrs["transcript_id"] = transcript_id

    if not next_attrs:
        return None
    parts[8] = gff3_attributes(next_attrs)
    return "\t".join(parts) + "\n"


def sortable_gtf_records(source: Path) -> Tuple[List[str], List[Tuple[str, int, int, str]]]:
    comments = ["##gff-version 3\n"]
    records: List[Tuple[str, int, int, str]] = []
    genes: Dict[str, Dict[str, Any]] = {}
    transcripts: Dict[str, Dict[str, Any]] = {}

    with source.open("rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            try:
                start = int(parts[3])
                end = int(parts[4])
            except ValueError:
                continue
            attrs = parse_gtf_attributes(parts[8])
            gene_id = attrs.get("gene_id")
            transcript_id = attrs.get("transcript_id")
            if not gene_id:
                continue

            gene = genes.setdefault(
                gene_id,
                {"seqid": parts[0], "source": parts[1], "start": start, "end": end, "strand": parts[6]},
            )
            gene["start"] = min(gene["start"], start)
            gene["end"] = max(gene["end"], end)

            if transcript_id:
                transcript = transcripts.setdefault(
                    transcript_id,
                    {
                        "gene_id": gene_id,
                        "seqid": parts[0],
                        "source": parts[1],
                        "start": start,
                        "end": end,
                        "strand": parts[6],
                    },
                )
                transcript["start"] = min(transcript["start"], start)
                transcript["end"] = max(transcript["end"], end)

            if parts[2] in {"gene", "transcript", "mRNA"}:
                continue
            converted = gtf_to_gff3_line(parts)
            if converted:
                records.append((parts[0], start, end, converted))

    for gene_id, gene in genes.items():
        attrs = gff3_attributes({"ID": gene_id, "Name": gene_id})
        line = f"{gene['seqid']}\t{gene['source']}\tgene\t{gene['start']}\t{gene['end']}\t.\t{gene['strand']}\t.\t{attrs}\n"
        records.append((gene["seqid"], int(gene["start"]), int(gene["end"]), line))
    for transcript_id, transcript in transcripts.items():
        attrs = gff3_attributes({"ID": transcript_id, "Parent": transcript["gene_id"], "Name": transcript_id})
        line = (
            f"{transcript['seqid']}\t{transcript['source']}\tmRNA\t{transcript['start']}\t{transcript['end']}\t.\t"
            f"{transcript['strand']}\t.\t{attrs}\n"
        )
        records.append((transcript["seqid"], int(transcript["start"]), int(transcript["end"]), line))

    records.sort(key=lambda item: (item[0], item[1], item[2]))
    return comments, records


def sortable_records(source: Path, annotation_format: str) -> Tuple[List[str], List[Tuple[str, int, int, str]]]:
    if annotation_format == "gtf":
        return sortable_gtf_records(source)
    comments: List[str] = []
    records: List[Tuple[str, int, int, str]] = []
    with source.open("rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip():
                continue
            if line.startswith("#"):
                if annotation_format != "gtf":
                    comments.append(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            try:
                start = int(parts[3])
                end = int(parts[4])
            except ValueError:
                continue
            records.append((parts[0], start, end, line if line.endswith("\n") else line + "\n"))
    records.sort(key=lambda item: (item[0], item[1], item[2]))
    return comments, records


def write_sorted_annotation(source: Path, sorted_path: Path, annotation_format: str) -> int:
    comments, records = sortable_records(source, annotation_format)
    sorted_path.parent.mkdir(parents=True, exist_ok=True)
    with sorted_path.open("wt", encoding="utf-8") as handle:
        for comment in comments:
            handle.write(comment if comment.endswith("\n") else comment + "\n")
        for record in records:
            handle.write(record[3])
    return len(records)


def bgzip_and_tabix(source: Path, sorted_path: Path, gz_path: Path, index_path: Path, index_type: str, annotation_format: str, force: bool, keep_sorted: bool) -> str:
    bgzip = tool_path("bgzip")
    tabix = tool_path("tabix")
    if force:
        for path in [sorted_path, gz_path, index_path]:
            if path.exists():
                path.unlink()
    if gz_path.exists() and index_path.exists():
        return "reused"
    write_sorted_annotation(source, sorted_path, annotation_format)
    with gz_path.open("wb") as output:
        subprocess.run([bgzip, "-c", str(sorted_path)], check=True, stdout=output)
    cmd = [tabix, "-f", "-p", "gff", str(gz_path)]
    if index_type == "CSI":
        cmd.insert(2, "-C")
    subprocess.run(cmd, check=True)
    if not index_path.exists():
        raise FileNotFoundError(f"tabix index was not created: {index_path}")
    if not keep_sorted and sorted_path.exists():
        sorted_path.unlink()
    return "created"


def register_annotation(conn: sqlite3.Connection, species_id: str, annotation_format: str, gz_path: Path, index_path: Path, index_type: str) -> None:
    runtime_format = "gff3" if annotation_format == "gtf" else annotation_format
    conn.execute(
        "INSERT OR REPLACE INTO annotation_assets(species_id, format, path) VALUES (?, ?, ?)",
        (species_id, f"{runtime_format}_tabix", str(gz_path)),
    )
    conn.execute(
        "INSERT OR REPLACE INTO annotation_assets(species_id, format, path) VALUES (?, ?, ?)",
        (species_id, f"{runtime_format}_tabix_{index_type.lower()}", str(index_path)),
    )


def annotation_index_type(conn: sqlite3.Connection, species_id: str) -> str:
    max_len = conn.execute("SELECT MAX(length) FROM chromosomes WHERE species_id = ?", (species_id,)).fetchone()[0] or 0
    return "CSI" if int(max_len) >= TBI_REFERENCE_LIMIT else "TBI"


def build_assets(args: argparse.Namespace) -> Dict[str, Any]:
    db_path = Path(args.db).resolve()
    assets_dir = Path(args.assets_dir).resolve()
    conn = sqlite3.connect(db_path)
    try:
        query = """
            SELECT species_id, format, path
            FROM annotation_assets
            WHERE format IN ('gff3', 'gtf')
        """
        params: Tuple[Any, ...] = ()
        if args.species:
            query += " AND species_id = ?"
            params = (args.species,)
        query += " ORDER BY species_id, format"
        annotations = rows(conn, query, params)
        summary = {"created": 0, "reused": 0, "species": []}
        conn.execute("BEGIN")
        for annotation in annotations:
            species_id = annotation["species_id"]
            annotation_format = annotation["format"]
            source = Path(annotation["path"])
            index_type = annotation_index_type(conn, species_id)
            sorted_path, gz_path, index_path = output_paths(species_id, annotation_format, source, assets_dir, index_type)
            status = bgzip_and_tabix(source, sorted_path, gz_path, index_path, index_type, annotation_format, args.force, args.keep_sorted)
            register_annotation(conn, species_id, annotation_format, gz_path, index_path, index_type)
            summary[status] += 1
            summary["species"].append(
                {
                    "species_id": species_id,
                    "format": annotation_format,
                    "status": status,
                    "index_type": index_type,
                    "annotation_gz": str(gz_path),
                    "tabix_index": str(index_path),
                }
            )
        conn.commit()
        return summary
    except Exception:
        conn.rollback()
        raise
    finally:
        conn.close()


def main() -> int:
    args = parse_args()
    try:
        summary = build_assets(args)
    except Exception as exc:
        print(f"ERROR: {exc}")
        return 1
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
