from __future__ import annotations

import datetime as dt
import re
import subprocess
import tempfile
from pathlib import Path

from ..config import ORTHOLOG_BIN
from ..db import connect, row, rows


REFERENCE_BY_GROUP = {
    "dicot": "arabidopsis_tair10_ensembl62",
    "monocot": "rice_rapdb_irgsp1_2026",
}


GENE_SYMBOL_PATTERN = re.compile(r"(?:^|\s)gene_symbol:([^\s;]+)", re.IGNORECASE)


def reference_gene_symbol(description: object) -> str | None:
    """Return an explicitly annotated reference gene symbol, when available."""
    match = GENE_SYMBOL_PATTERN.search(str(description or ""))
    return match.group(1) if match else None


def annotation_symbols(conn, reference_id: str, gene_ids: list[str]) -> dict[str, str]:
    if not gene_ids:
        return {}
    placeholders = ",".join("?" for _ in gene_ids)
    output: dict[str, str] = {}
    for item in rows(conn, f"""
        SELECT gene_id, alias FROM reference_gene_aliases
        WHERE reference_id = ? AND alias_type = 'annotation_symbol'
          AND gene_id IN ({placeholders}) ORDER BY alias
    """, [reference_id, *gene_ids]):
        output.setdefault(str(item["gene_id"]), str(item["alias"]))
    return output


def ensure_schema(conn) -> None:
    conn.executescript("""
        CREATE TABLE IF NOT EXISTS protein_similarity_hits (
            species_id TEXT NOT NULL,
            gene_id TEXT NOT NULL,
            protein_id TEXT NOT NULL,
            reference_id TEXT NOT NULL,
            reference_protein_id TEXT NOT NULL,
            reference_gene_id TEXT,
            pident REAL, alignment_length INTEGER, evalue REAL, bitscore REAL,
            rank INTEGER NOT NULL, method TEXT NOT NULL, created_at TEXT NOT NULL,
            PRIMARY KEY (species_id, protein_id, reference_id, reference_protein_id, method)
        );
        CREATE INDEX IF NOT EXISTS idx_protein_similarity_lookup
          ON protein_similarity_hits(species_id, gene_id, protein_id, reference_id, rank);
    """)


def reference_id_for_protein(conn, species_id: str) -> str:
    species = row(conn, "SELECT taxon_group FROM species WHERE species_id = ?", (species_id,))
    reference_id = REFERENCE_BY_GROUP.get(str(species.get("taxon_group") if species else ""))
    if not reference_id:
        raise ValueError("no reference proteome is configured for this species")
    return reference_id


def result_rows(conn, species_id: str, gene_id: str, protein_id: str, reference_id: str, limit: int | None) -> list[dict]:
    hits = rows(conn, """
        SELECT h.reference_id, rs.display_name AS reference_display_name,
               h.reference_gene_id, rg.seqid AS reference_seqid, rg.start AS reference_start,
               rg.end AS reference_end, rg.strand AS reference_strand,
               h.reference_protein_id, rp.transcript_id AS reference_transcript_id,
               rp.description AS reference_protein_description,
               h.protein_id AS query_protein_id, h.pident, h.alignment_length,
               h.evalue, h.bitscore, h.rank, h.method, h.created_at
        FROM protein_similarity_hits h
        JOIN reference_species rs ON rs.reference_id = h.reference_id
        LEFT JOIN reference_proteins rp ON rp.reference_id = h.reference_id AND rp.protein_id = h.reference_protein_id
        LEFT JOIN reference_genes rg ON rg.reference_id = h.reference_id AND rg.gene_id = h.reference_gene_id
        WHERE h.species_id = ? AND h.gene_id = ? AND h.protein_id = ? AND h.reference_id = ?
        ORDER BY h.rank LIMIT ?
    """, (species_id, gene_id, protein_id, reference_id, -1 if limit is None else limit))
    symbols = annotation_symbols(conn, reference_id, sorted({str(hit["reference_gene_id"]) for hit in hits if hit.get("reference_gene_id")}))
    for hit in hits:
        hit["reference_gene_symbol"] = symbols.get(str(hit.get("reference_gene_id") or "")) or reference_gene_symbol(hit.get("reference_protein_description"))
    return hits


def protein_similarity(species_id: str, gene_id: str, protein_id: str, limit: int = 5) -> dict:
    with connect() as conn:
        ensure_schema(conn)
        protein = row(conn, "SELECT protein_id, sequence FROM proteins WHERE species_id = ? AND gene_id = ? AND protein_id = ?", (species_id, gene_id, protein_id))
        if not protein:
            raise ValueError("protein does not belong to gene")
        reference_id = reference_id_for_protein(conn, species_id)
        cached = result_rows(conn, species_id, gene_id, protein_id, reference_id, limit)
        if cached:
            return {"species_id": species_id, "gene_id": gene_id, "protein_id": protein_id, "hit_count": len(cached), "cached": True, "hits": cached}
        asset = row(conn, "SELECT index_prefix FROM reference_search_assets WHERE reference_id = ? AND asset_type = 'diamond_protein' AND status = 'ready'", (reference_id,))
        if not asset:
            raise ValueError("reference DIAMOND index is unavailable")
        sequence = str(protein["sequence"] or "").replace("*", "").strip()
        if not sequence:
            raise ValueError("protein sequence is empty")
        with tempfile.TemporaryDirectory(prefix="apdb_similarity_") as tempdir:
            query = Path(tempdir) / "query.fa"; output = Path(tempdir) / "hits.tsv"
            query.write_text(f">{protein_id}\n{sequence}\n", encoding="ascii")
            subprocess.run([str(ORTHOLOG_BIN / "diamond"), "blastp", "--query", str(query), "--db", str(asset["index_prefix"]), "--out", str(output), "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "--max-target-seqs", str(limit), "--evalue", "1e-5", "--threads", "2", "--quiet"], check=True, timeout=90)
            parsed = [line.rstrip("\n").split("\t") for line in output.read_text(encoding="utf-8").splitlines()]
        ref_genes = {item["protein_id"]: item["gene_id"] for item in rows(conn, "SELECT protein_id, gene_id FROM reference_proteins WHERE reference_id = ?", (reference_id,))}
        now = dt.datetime.now(dt.timezone.utc).isoformat()
        payload = [(species_id, gene_id, protein_id, reference_id, item[1], ref_genes.get(item[1]), float(item[2]), int(item[3]), float(item[4]), float(item[5]), index + 1, "diamond", now) for index, item in enumerate(parsed)]
        conn.execute("DELETE FROM protein_similarity_hits WHERE species_id = ? AND protein_id = ? AND reference_id = ? AND method = 'diamond'", (species_id, protein_id, reference_id))
        conn.executemany("INSERT INTO protein_similarity_hits VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", payload)
        conn.commit()
        hits = result_rows(conn, species_id, gene_id, protein_id, reference_id, limit)
    return {"species_id": species_id, "gene_id": gene_id, "protein_id": protein_id, "hit_count": len(hits), "cached": False, "hits": hits}
