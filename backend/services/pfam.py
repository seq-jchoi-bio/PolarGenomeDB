from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path
from typing import Any

from ..config import HMMSCAN, PFAM_DB


def ensure_pfam_schema(conn) -> None:
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS pfam_domains (
            species_id TEXT NOT NULL,
            protein_id TEXT NOT NULL,
            gene_id TEXT,
            transcript_id TEXT,
            pfam_id TEXT NOT NULL,
            domain_name TEXT,
            domain_length INTEGER,
            start INTEGER NOT NULL,
            end INTEGER NOT NULL,
            hmm_start INTEGER,
            hmm_end INTEGER,
            i_evalue REAL,
            c_evalue TEXT,
            full_evalue TEXT,
            score REAL,
            accuracy REAL,
            description TEXT,
            source TEXT DEFAULT 'hmmscan_pfam_a',
            PRIMARY KEY (species_id, protein_id, pfam_id, start, end)
        )
        """
    )
    conn.execute("CREATE INDEX IF NOT EXISTS idx_pfam_species_gene ON pfam_domains(species_id, gene_id)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_pfam_species_protein ON pfam_domains(species_id, protein_id)")
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS pfam_protein_status (
            species_id TEXT NOT NULL,
            protein_id TEXT NOT NULL,
            gene_id TEXT,
            transcript_id TEXT,
            status TEXT NOT NULL,
            domain_count INTEGER NOT NULL DEFAULT 0,
            source TEXT DEFAULT 'hmmscan_pfam_a',
            updated_at TEXT DEFAULT CURRENT_TIMESTAMP,
            PRIMARY KEY (species_id, protein_id)
        )
        """
    )
    conn.execute("CREATE INDEX IF NOT EXISTS idx_pfam_status_species_gene ON pfam_protein_status(species_id, gene_id)")


def fasta_record(identifier: str, sequence: str) -> str:
    clean = "".join(str(sequence or "").split())
    chunks = "\n".join(clean[index : index + 60] for index in range(0, len(clean), 60))
    return f">{identifier}\n{chunks}\n"


def parse_domtblout(path: Path, i_evalue_cutoff: float = 1e-3) -> list[dict[str, Any]]:
    domains: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split(maxsplit=22)
            if len(parts) < 22:
                continue
            try:
                i_evalue = float(parts[12])
            except ValueError:
                continue
            if i_evalue > i_evalue_cutoff:
                continue
            try:
                domains.append(
                    {
                        "domain_name": parts[0],
                        "pfam_id": parts[1],
                        "domain_length": int(parts[2]),
                        "query": parts[3],
                        "query_length": int(parts[5]),
                        "full_evalue": parts[6],
                        "full_score": float(parts[7]),
                        "domain_number": int(parts[9]),
                        "domain_count": int(parts[10]),
                        "c_evalue": parts[11],
                        "i_evalue": i_evalue,
                        "score": float(parts[13]),
                        "hmm_start": int(parts[15]),
                        "hmm_end": int(parts[16]),
                        "start": int(parts[17]),
                        "end": int(parts[18]),
                        "env_start": int(parts[19]),
                        "env_end": int(parts[20]),
                        "accuracy": float(parts[21]),
                        "description": parts[22] if len(parts) > 22 else "",
                    }
                )
            except ValueError:
                continue
    domains.sort(key=lambda item: (item["start"], item["i_evalue"], -item["score"]))
    return domains


def run_hmmscan(
    records: list[tuple[str, str]],
    *,
    cpu: int = 4,
    pfam_db: Path = PFAM_DB,
    hmmscan: Path = HMMSCAN,
    timeout: int | None = None,
) -> list[dict[str, Any]]:
    if not records:
        return []
    if not hmmscan.exists():
        raise FileNotFoundError(f"hmmscan not found: {hmmscan}")
    if not pfam_db.exists():
        raise FileNotFoundError(f"Pfam database not found: {pfam_db}")

    with tempfile.TemporaryDirectory(prefix="pgdb_pfam_") as tmpdir:
        tmp = Path(tmpdir)
        query_fasta = tmp / "query.fa"
        domtblout = tmp / "pfam.domtblout"
        with query_fasta.open("w", encoding="utf-8") as handle:
            for identifier, sequence in records:
                handle.write(fasta_record(identifier, sequence))
        subprocess.run(
            [
                str(hmmscan),
                "--cpu",
                str(cpu),
                "--noali",
                "--domtblout",
                str(domtblout),
                str(pfam_db),
                str(query_fasta),
            ],
            check=True,
            timeout=timeout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        return parse_domtblout(domtblout)


def insert_pfam_domains(conn, species_id: str, protein_rows: dict[str, dict[str, Any]], domains: list[dict[str, Any]]) -> int:
    ensure_pfam_schema(conn)
    count = 0
    for domain in domains:
      protein_id = str(domain.get("query") or "")
      protein = protein_rows.get(protein_id)
      if not protein:
          continue
      conn.execute(
          """
          INSERT OR REPLACE INTO pfam_domains
          (species_id, protein_id, gene_id, transcript_id, pfam_id, domain_name, domain_length,
           start, end, hmm_start, hmm_end, i_evalue, c_evalue, full_evalue, score, accuracy, description, source)
          VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
          """,
          (
              species_id,
              protein_id,
              protein.get("gene_id"),
              protein.get("transcript_id"),
              domain.get("pfam_id"),
              domain.get("domain_name"),
              domain.get("domain_length"),
              domain.get("start"),
              domain.get("end"),
              domain.get("hmm_start"),
              domain.get("hmm_end"),
              domain.get("i_evalue"),
              domain.get("c_evalue"),
              domain.get("full_evalue"),
              domain.get("score"),
              domain.get("accuracy"),
              domain.get("description"),
              "hmmscan_pfam_a",
          ),
      )
      count += 1
    return count


def mark_pfam_scanned(conn, species_id: str, protein_rows: dict[str, dict[str, Any]], domains: list[dict[str, Any]]) -> None:
    ensure_pfam_schema(conn)
    domain_counts: dict[str, int] = {}
    for domain in domains:
        protein_id = str(domain.get("query") or "")
        domain_counts[protein_id] = domain_counts.get(protein_id, 0) + 1
    for protein_id, protein in protein_rows.items():
        domain_count = domain_counts.get(protein_id, 0)
        conn.execute(
            """
            INSERT OR REPLACE INTO pfam_protein_status
            (species_id, protein_id, gene_id, transcript_id, status, domain_count, source, updated_at)
            VALUES (?, ?, ?, ?, ?, ?, ?, CURRENT_TIMESTAMP)
            """,
            (
                species_id,
                protein_id,
                protein.get("gene_id"),
                protein.get("transcript_id"),
                "completed",
                domain_count,
                "hmmscan_pfam_a",
            ),
        )
