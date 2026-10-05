#!/usr/bin/env python3
"""Precompute species-to-reference protein similarity hits."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import shutil
import sqlite3
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, Iterable, Optional


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = ROOT / "data" / "antarctic_plant.db"
ORTHOLOG_BIN = Path(os.environ.get("ANTARCTIC_ORTHOLOG_BIN", Path(sys.executable).parent))
DEFAULT_REFERENCE_BY_GROUP = {
    "dicot": "arabidopsis_tair10_ensembl62",
    "monocot": "rice_rapdb_irgsp1_2026",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species", required=True, help="Species ID to process.")
    parser.add_argument("--reference-id", help="Override reference ID. Defaults from species taxon_group.")
    parser.add_argument("--db", default=str(DEFAULT_DB), help="SQLite DB path.")
    parser.add_argument("--threads", type=int, default=max(1, min(8, os.cpu_count() or 1)), help="DIAMOND threads.")
    parser.add_argument("--max-target-seqs", type=int, default=10, help="DIAMOND max targets per protein.")
    parser.add_argument("--top-per-gene", type=int, default=5, help="Stored hits per gene.")
    parser.add_argument("--evalue", default="1e-5", help="DIAMOND e-value threshold.")
    parser.add_argument("--min-pident", type=float, default=20.0, help="Minimum percent identity to store.")
    parser.add_argument("--force", action="store_true", help="Replace existing hits for species/reference/method.")
    return parser.parse_args()


def rows(conn: sqlite3.Connection, query: str, params: Iterable[object] = ()) -> list[sqlite3.Row]:
    return conn.execute(query, tuple(params)).fetchall()


def diamond_path() -> str:
    for candidate in [str(ORTHOLOG_BIN / "diamond"), "diamond"]:
        if Path(candidate).exists():
            return candidate
        found = shutil.which(candidate)
        if found:
            return found
    raise FileNotFoundError("diamond not found")


def wrap_sequence(sequence: str, width: int = 60) -> str:
    return "\n".join(sequence[i : i + width] for i in range(0, len(sequence), width))


def write_species_proteins(conn: sqlite3.Connection, species_id: str, path: Path) -> Dict[str, str]:
    protein_to_gene: Dict[str, str] = {}
    with path.open("w", encoding="ascii") as handle:
        for protein_id, gene_id, sequence in rows(
            conn,
            """
            SELECT protein_id, gene_id, sequence
            FROM proteins
            WHERE species_id = ? AND sequence != ''
            ORDER BY gene_id, protein_id
            """,
            (species_id,),
        ):
            if not gene_id:
                continue
            protein_to_gene[str(protein_id)] = str(gene_id)
            handle.write(f">{protein_id}\n{wrap_sequence(str(sequence))}\n")
    if not protein_to_gene:
        raise ValueError(f"no proteins with gene links found for species: {species_id}")
    return protein_to_gene


def reference_for_species(conn: sqlite3.Connection, species_id: str, override: Optional[str]) -> str:
    if override:
        return override
    row = conn.execute("SELECT taxon_group FROM species WHERE species_id = ?", (species_id,)).fetchone()
    if not row:
        raise ValueError(f"species not found: {species_id}")
    reference_id = DEFAULT_REFERENCE_BY_GROUP.get(str(row[0]))
    if not reference_id:
        raise ValueError(f"no default reference configured for taxon_group: {row[0]}")
    return reference_id


def reference_diamond_db(conn: sqlite3.Connection, reference_id: str) -> Path:
    row = conn.execute(
        """
        SELECT index_prefix
        FROM reference_search_assets
        WHERE reference_id = ? AND asset_type = 'diamond_protein' AND status = 'ready'
        """,
        (reference_id,),
    ).fetchone()
    if not row:
        raise ValueError(f"reference DIAMOND asset not found: {reference_id}")
    path = Path(str(row[0]))
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def reference_gene_map(conn: sqlite3.Connection, reference_id: str) -> Dict[str, Optional[str]]:
    return {
        str(protein_id): (str(gene_id) if gene_id else None)
        for protein_id, gene_id in rows(
            conn,
            "SELECT protein_id, gene_id FROM reference_proteins WHERE reference_id = ?",
            (reference_id,),
        )
    }


def run_diamond(args: argparse.Namespace, query_fasta: Path, db_path: Path, out_path: Path) -> None:
    cmd = [
        diamond_path(),
        "blastp",
        "--query",
        str(query_fasta),
        "--db",
        str(db_path),
        "--out",
        str(out_path),
        "--outfmt",
        "6",
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "evalue",
        "bitscore",
        "--threads",
        str(args.threads),
        "--max-target-seqs",
        str(args.max_target_seqs),
        "--evalue",
        str(args.evalue),
        "--quiet",
    ]
    subprocess.run(cmd, check=True)


def parse_hits(path: Path, protein_to_gene: Dict[str, str], ref_to_gene: Dict[str, Optional[str]], min_pident: float) -> Dict[str, list[tuple]]:
    by_gene: Dict[str, list[tuple]] = {}
    with path.open("rt", encoding="utf-8") as handle:
        for raw in handle:
            qseqid, sseqid, pident, length, evalue, bitscore = raw.rstrip("\n").split("\t")
            pident_value = float(pident)
            if pident_value < min_pident:
                continue
            gene_id = protein_to_gene.get(qseqid)
            if not gene_id:
                continue
            by_gene.setdefault(gene_id, []).append(
                (
                    gene_id,
                    qseqid,
                    sseqid,
                    ref_to_gene.get(sseqid),
                    pident_value,
                    int(length),
                    float(evalue),
                    float(bitscore),
                )
            )
    for gene_id, hits in by_gene.items():
        hits.sort(key=lambda item: (-item[7], item[6], -item[4], item[2]))
    return by_gene


def build_similarity(args: argparse.Namespace) -> dict:
    db_path = Path(args.db).resolve()
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    try:
        reference_id = reference_for_species(conn, args.species, args.reference_id)
        db_path = reference_diamond_db(conn, reference_id)
        ref_to_gene = reference_gene_map(conn, reference_id)
        with tempfile.TemporaryDirectory(prefix=f"{args.species}_similarity_") as tempdir:
            query_fasta = Path(tempdir) / f"{args.species}.proteins.fa"
            out_path = Path(tempdir) / "diamond.tsv"
            protein_to_gene = write_species_proteins(conn, args.species, query_fasta)
            run_diamond(args, query_fasta, db_path, out_path)
            hits_by_gene = parse_hits(out_path, protein_to_gene, ref_to_gene, args.min_pident)

        payload = []
        now = dt.datetime.now(dt.timezone.utc).isoformat()
        for gene_id, hits in hits_by_gene.items():
            seen_refs = set()
            rank = 0
            for hit in hits:
                ref_key = hit[3] or hit[2]
                if ref_key in seen_refs:
                    continue
                seen_refs.add(ref_key)
                rank += 1
                if rank > args.top_per_gene:
                    break
                gene_id, protein_id, reference_protein_id, reference_gene_id, pident, length, evalue, bitscore = hit
                payload.append(
                    (
                        args.species,
                        gene_id,
                        protein_id,
                        reference_id,
                        reference_protein_id,
                        reference_gene_id,
                        pident,
                        length,
                        evalue,
                        bitscore,
                        rank,
                        "diamond",
                        now,
                    )
                )

        conn.execute("BEGIN")
        if args.force:
            conn.execute(
                "DELETE FROM similarity_hits WHERE species_id = ? AND reference_id = ? AND method = 'diamond'",
                (args.species, reference_id),
            )
        conn.executemany(
            """
            INSERT OR REPLACE INTO similarity_hits
            (species_id, gene_id, protein_id, reference_id, reference_protein_id, reference_gene_id,
             pident, alignment_length, evalue, bitscore, rank, method, created_at)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            payload,
        )
        conn.commit()
        return {
            "species_id": args.species,
            "reference_id": reference_id,
            "method": "diamond",
            "query_proteins": len(protein_to_gene),
            "genes_with_hits": len(hits_by_gene),
            "stored_hits": len(payload),
            "top_per_gene": args.top_per_gene,
            "min_pident": args.min_pident,
            "evalue": args.evalue,
        }
    except Exception:
        conn.rollback()
        raise
    finally:
        conn.close()


def main() -> int:
    args = parse_args()
    try:
        print(json.dumps(build_similarity(args), indent=2, sort_keys=True))
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
