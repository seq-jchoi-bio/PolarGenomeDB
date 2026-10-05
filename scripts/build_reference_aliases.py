#!/usr/bin/env python3
"""Refresh queryable reference gene aliases without rebuilding the reference indexes."""

from __future__ import annotations

import argparse
import sqlite3
import sys
from pathlib import Path

from build_reference import (
    DEFAULT_DB,
    collect_gene_aliases,
    init_db,
    parse_gene_models,
    parse_simple_yaml,
    read_id_mapping,
    resolve_path,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, help="Reference manifest YAML file.")
    parser.add_argument("--db", default=str(DEFAULT_DB), help="SQLite DB path.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    manifest_path = Path(args.manifest).resolve()
    manifest = parse_simple_yaml(manifest_path)
    reference_id = str(manifest["reference_id"])
    annotation_path = resolve_path(manifest_path, manifest.get("annotation", {}).get("file"))
    optional = manifest.get("optional", {})
    mapping_path = resolve_path(manifest_path, optional.get("id_mapping") or optional.get("rap_msu_mapping"))
    if not annotation_path or not annotation_path.exists():
        raise FileNotFoundError(f"annotation not found: {annotation_path}")

    _genes, transcripts, _seqids, _transcript_to_gene = parse_gene_models(annotation_path)
    protein_to_gene, _protein_to_transcript = read_id_mapping(mapping_path)
    aliases = collect_gene_aliases(transcripts, protein_to_gene)

    db_path = Path(args.db).resolve()
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA foreign_keys = ON")
    try:
        init_db(conn)
        if not conn.execute("SELECT 1 FROM reference_species WHERE reference_id = ?", (reference_id,)).fetchone():
            raise ValueError(f"reference is not built in SQLite: {reference_id}")
        current_gene_ids = {
            row[0] for row in conn.execute(
                "SELECT gene_id FROM reference_genes WHERE reference_id = ?", (reference_id,)
            )
        }
        aliases = {item for item in aliases if item[0] in current_gene_ids}
        conn.execute("DELETE FROM reference_gene_aliases WHERE reference_id = ?", (reference_id,))
        conn.executemany(
            "INSERT INTO reference_gene_aliases VALUES (?, ?, ?, ?)",
            [(reference_id, gene_id, alias, alias_type) for gene_id, alias, alias_type in sorted(aliases)],
        )
        conn.commit()
    finally:
        conn.close()
    print(f"reference_id={reference_id} aliases={len(aliases)}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
