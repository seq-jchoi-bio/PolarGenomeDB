#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from backend.db import connect, rows
from backend.services.pfam import ensure_pfam_schema, insert_pfam_domains, mark_pfam_scanned, run_hmmscan


def chunks(items: list[dict], size: int):
    for index in range(0, len(items), size):
        yield items[index : index + size]


def build_species(species_id: str, *, batch_size: int, cpu: int, refresh: bool, limit: int | None) -> dict:
    started = time.time()
    with connect() as conn:
        ensure_pfam_schema(conn)
        params: list[object] = [species_id]
        sql = """
            SELECT p.species_id, p.protein_id, p.gene_id, p.transcript_id, p.length, p.sequence
            FROM proteins p
            WHERE p.species_id = ?
        """
        if not refresh:
            sql += """
              AND NOT EXISTS (
                SELECT 1 FROM pfam_protein_status s
                WHERE s.species_id = p.species_id AND s.protein_id = p.protein_id
                  AND s.status = 'completed'
              )
            """
        sql += " ORDER BY p.protein_id"
        if limit:
            sql += " LIMIT ?"
            params.append(limit)
        proteins = rows(conn, sql, params)

    processed = 0
    inserted = 0
    failed_batches = 0
    for batch in chunks(proteins, batch_size):
        records = [(item["protein_id"], item["sequence"]) for item in batch]
        protein_rows = {item["protein_id"]: item for item in batch}
        try:
            domains = run_hmmscan(records, cpu=cpu)
        except Exception as exc:
            failed_batches += 1
            print(f"[WARN] {species_id} batch failed at protein {processed + 1}: {exc}", flush=True)
            continue
        with connect() as conn:
            inserted += insert_pfam_domains(conn, species_id, protein_rows, domains)
            mark_pfam_scanned(conn, species_id, protein_rows, domains)
            conn.commit()
        processed += len(batch)
        elapsed = time.time() - started
        print(
            f"[{species_id}] processed={processed}/{len(proteins)} inserted_domains={inserted} elapsed={elapsed/60:.1f}m",
            flush=True,
        )

    return {
        "species_id": species_id,
        "proteins_selected": len(proteins),
        "proteins_processed": processed,
        "domains_inserted": inserted,
        "failed_batches": failed_batches,
        "elapsed_sec": round(time.time() - started, 2),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Build query-side Pfam domains with hmmscan.")
    parser.add_argument("--species", action="append", help="Species ID. May be repeated. Defaults to all species in DB.")
    parser.add_argument("--batch-size", type=int, default=400)
    parser.add_argument("--cpu", type=int, default=4)
    parser.add_argument("--refresh", action="store_true", help="Re-run proteins even when cached domains exist.")
    parser.add_argument("--limit", type=int, default=None, help="Optional smoke-test protein limit per species.")
    args = parser.parse_args()

    with connect() as conn:
        species_ids = args.species or [item["species_id"] for item in rows(conn, "SELECT species_id FROM species ORDER BY species_id")]

    summaries = [
        build_species(species_id, batch_size=args.batch_size, cpu=args.cpu, refresh=args.refresh, limit=args.limit)
        for species_id in species_ids
    ]
    print("Pfam build summary")
    for summary in summaries:
        print(summary)
    if any(summary['failed_batches'] for summary in summaries):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
