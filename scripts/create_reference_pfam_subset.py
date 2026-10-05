#!/usr/bin/env python3
"""Create a compact SQLite input containing one reference proteome for Pfam batching."""

from __future__ import annotations

import argparse
import sqlite3
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.build_reference import init_db


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-db", required=True)
    parser.add_argument("--output-db", required=True)
    parser.add_argument("--reference", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    output = Path(args.output_db).resolve()
    if output.exists():
        output.unlink()
    source = sqlite3.connect(Path(args.source_db).resolve())
    source.row_factory = sqlite3.Row
    destination = sqlite3.connect(output)
    destination.execute("PRAGMA foreign_keys = ON")
    try:
        init_db(destination)
        reference = source.execute("SELECT * FROM reference_species WHERE reference_id = ?", (args.reference,)).fetchone()
        if not reference:
            raise ValueError(f"reference not found: {args.reference}")
        columns = list(reference.keys())
        destination.execute(
            f"INSERT INTO reference_species ({','.join(columns)}) VALUES ({','.join('?' for _ in columns)})",
            [reference[column] for column in columns],
        )
        proteins = source.execute(
            "SELECT reference_id, protein_id, gene_id, transcript_id, description, length, sequence FROM reference_proteins WHERE reference_id = ?",
            (args.reference,),
        ).fetchall()
        destination.executemany(
            "INSERT INTO reference_proteins VALUES (?, ?, ?, ?, ?, ?, ?)",
            [tuple(record) for record in proteins],
        )
        destination.commit()
    finally:
        source.close()
        destination.close()
    print(f"reference_id={args.reference} proteins={len(proteins)} output_db={output}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
