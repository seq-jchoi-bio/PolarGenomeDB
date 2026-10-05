#!/usr/bin/env python3
"""Merge one reference Pfam table from a completed batch SQLite file."""

from __future__ import annotations

import argparse
import sqlite3
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", required=True, help="Local Polar Genome DB SQLite path.")
    parser.add_argument("--source-db", required=True, help="Completed reference Pfam batch SQLite path.")
    parser.add_argument("--reference", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    db_path = Path(args.db).resolve()
    source_path = Path(args.source_db).resolve()
    conn = sqlite3.connect(db_path)
    try:
        conn.execute("ATTACH DATABASE ? AS batch", (str(source_path),))
        count = conn.execute("SELECT COUNT(*) FROM batch.reference_pfam_domains WHERE reference_id = ?", (args.reference,)).fetchone()[0]
        conn.execute("DELETE FROM reference_pfam_domains WHERE reference_id = ?", (args.reference,))
        conn.execute(
            """
            INSERT INTO reference_pfam_domains
            (reference_id, protein_id, pfam_id, domain_name, start, end, evalue, score, source)
            SELECT reference_id, protein_id, pfam_id, domain_name, start, end, evalue, score, source
            FROM batch.reference_pfam_domains
            WHERE reference_id = ?
            """,
            (args.reference,),
        )
        conn.commit()
    finally:
        conn.close()
    print(f"reference_id={args.reference} imported_domains={count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
