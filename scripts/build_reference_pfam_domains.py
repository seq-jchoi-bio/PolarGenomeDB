#!/usr/bin/env python3
"""Build cached Pfam-A domains for one reference proteome."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from backend.config import DB_PATH
from backend.db import connect, rows
from backend.services.pfam import run_hmmscan
from build_reference import init_db


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, help="Reference ID in SQLite.")
    parser.add_argument("--db", default=str(DB_PATH), help="SQLite DB path.")
    parser.add_argument("--cpu", type=int, default=8, help="hmmscan CPUs per batch.")
    parser.add_argument("--batch-size", type=int, default=400, help="Protein records per hmmscan batch.")
    parser.add_argument("--replace", action="store_true", help="Replace existing domains for this reference.")
    return parser.parse_args()


def chunks(items: list[dict], size: int):
    for start in range(0, len(items), size):
        yield items[start:start + size]


def main() -> int:
    args = parse_args()
    if args.cpu < 1 or args.batch_size < 1:
        raise ValueError("--cpu and --batch-size must be positive")
    db_path = Path(args.db).resolve()
    with connect(db_path) as conn:
        init_db(conn)
        if not conn.execute("SELECT 1 FROM reference_species WHERE reference_id = ?", (args.reference,)).fetchone():
            raise ValueError(f"reference not found: {args.reference}")
        if args.replace:
            conn.execute("DELETE FROM reference_pfam_domains WHERE reference_id = ?", (args.reference,))
            conn.commit()
        proteins = rows(
            conn,
            "SELECT protein_id, sequence FROM reference_proteins WHERE reference_id = ? ORDER BY protein_id",
            (args.reference,),
        )
    inserted = 0
    for batch_index, batch in enumerate(chunks(proteins, args.batch_size), start=1):
        domains = run_hmmscan([(item["protein_id"], item["sequence"]) for item in batch], cpu=args.cpu)
        with connect(db_path) as conn:
            conn.executemany(
                """
                INSERT OR REPLACE INTO reference_pfam_domains
                (reference_id, protein_id, pfam_id, domain_name, start, end, evalue, score, source)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                [
                    (
                        args.reference,
                        domain["query"],
                        domain["pfam_id"],
                        domain.get("domain_name"),
                        domain["start"],
                        domain["end"],
                        domain.get("i_evalue"),
                        domain.get("score"),
                        "hmmscan_pfam_a",
                    )
                    for domain in domains
                ],
            )
            conn.commit()
        inserted += len(domains)
        print(f"[{args.reference}] batch={batch_index} proteins={min(batch_index * args.batch_size, len(proteins))}/{len(proteins)} domains={inserted}", flush=True)
    print(f"reference_id={args.reference} proteins={len(proteins)} inserted_domains={inserted}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
