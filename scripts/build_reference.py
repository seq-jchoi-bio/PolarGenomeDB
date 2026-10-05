#!/usr/bin/env python3
"""Build reference genome metadata and search indexes for 2nd-phase APIs."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import json
import os
import re
import shutil
import sqlite3
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Optional, Tuple
from urllib.parse import unquote


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = ROOT / "data" / "antarctic_plant.db"
DEFAULT_ASSETS_DIR = ROOT / "data" / "reference_assets"
ORTHOLOG_BIN = Path(os.environ.get("ANTARCTIC_ORTHOLOG_BIN", Path(sys.executable).parent))
ALLOWED_TAXON_GROUPS = {"monocot", "dicot", "other", "unknown"}
GENE_TYPES = {"gene", "pseudogene"}
TRANSCRIPT_TYPES = {"mRNA", "transcript", "lnc_RNA", "rRNA", "tRNA", "miRNA", "ncRNA"}
GENE_SYMBOL_ATTRIBUTE_KEYS = (
    "Oryzabase Gene Symbol Synonym(s)",
    "Oryzabase Gene Symbol",
    "gene_symbol",
    "gene_symbols",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, help="Reference manifest YAML file.")
    parser.add_argument("--db", default=str(DEFAULT_DB), help="SQLite DB path.")
    parser.add_argument("--assets-dir", default=str(DEFAULT_ASSETS_DIR), help="Reference asset output directory.")
    parser.add_argument("--threads", type=int, default=max(1, min(8, os.cpu_count() or 1)), help="Threads for DIAMOND.")
    parser.add_argument("--replace", action="store_true", help="Replace an existing reference entry.")
    parser.add_argument("--force-index", action="store_true", help="Rebuild search indexes even if outputs already exist.")
    parser.add_argument("--skip-indexes", action="store_true", help="Register metadata only; do not create DIAMOND/BLAST indexes.")
    parser.add_argument("--check-only", action="store_true", help="Validate inputs and report summary without writing SQLite.")
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
        while stack and indent <= stack[-1][0]:
            stack.pop()
        parent = stack[-1][1]
        key = key.strip()
        value = value.strip()
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


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def parse_attributes(raw: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    if "=" in raw:
        for item in raw.split(";"):
            if not item.strip() or "=" not in item:
                continue
            key, value = item.split("=", 1)
            attrs[key.strip()] = value.strip()
        return attrs
    for match in re.finditer(r'(\S+)\s+"([^"]+)"', raw):
        attrs[match.group(1)] = match.group(2)
    return attrs


def first_attr(attrs: Dict[str, str], keys: Iterable[str]) -> Optional[str]:
    for key in keys:
        value = attrs.get(key)
        if value:
            return value.split(",")[0]
    return None


def fasta_records(path: Path) -> Iterator[Tuple[str, str, str]]:
    record_id: Optional[str] = None
    description = ""
    chunks: List[str] = []
    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if record_id is not None:
                    yield record_id, description, "".join(chunks)
                description = line[1:]
                record_id = description.split()[0]
                chunks = []
            else:
                chunks.append(line)
    if record_id is not None:
        yield record_id, description, "".join(chunks)


def create_fai(fasta_path: Path, fai_path: Path) -> Path:
    if fai_path.exists() and fai_path.stat().st_mtime >= fasta_path.stat().st_mtime:
        return fai_path
    fai_path.parent.mkdir(parents=True, exist_ok=True)
    with fasta_path.open("rb") as src, fai_path.open("w", encoding="ascii") as out:
        name: Optional[str] = None
        length = 0
        offset = 0
        line_bases: Optional[int] = None
        line_width: Optional[int] = None
        while True:
            pos = src.tell()
            raw = src.readline()
            if not raw:
                break
            if raw.startswith(b">"):
                if name is not None:
                    out.write(f"{name}\t{length}\t{offset}\t{line_bases or 0}\t{line_width or 0}\n")
                name = raw[1:].decode("utf-8", errors="replace").strip().split()[0]
                length = 0
                offset = src.tell()
                line_bases = None
                line_width = None
                continue
            stripped = raw.rstrip(b"\r\n")
            if name is None or not stripped:
                continue
            if line_bases is None:
                line_bases = len(stripped)
                line_width = len(raw)
            length += len(stripped)
        if name is not None:
            out.write(f"{name}\t{length}\t{offset}\t{line_bases or 0}\t{line_width or 0}\n")
    return fai_path


def parse_fai(path: Path) -> List[Dict[str, Any]]:
    records = []
    with path.open("rt", encoding="ascii") as handle:
        for raw in handle:
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            records.append(
                {
                    "seqid": fields[0],
                    "length": int(fields[1]),
                    "offset": int(fields[2]),
                    "line_bases": int(fields[3]),
                    "line_width": int(fields[4]),
                }
            )
    return records


def parse_gene_models(path: Path) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], set[str], Dict[str, str]]:
    genes: List[Dict[str, Any]] = []
    transcripts: List[Dict[str, Any]] = []
    seqids: set[str] = set()
    transcript_to_gene: Dict[str, str] = {}
    synthetic_genes: Dict[str, Dict[str, Any]] = {}
    with open_text(path) as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            seqid, source, feature_type, start, end, _score, strand, _phase, attrs_raw = fields
            attrs = parse_attributes(attrs_raw)
            seqids.add(seqid)
            feature_id = first_attr(attrs, ["ID", "gene_id", "transcript_id", "Name"])
            parent_id = first_attr(attrs, ["Parent", "gene_id"])
            if feature_type in GENE_TYPES:
                gene_id = first_attr(attrs, ["ID", "gene_id", "Name"])
                if gene_id:
                    genes.append(
                        {
                            "gene_id": gene_id,
                            "seqid": seqid,
                            "start": int(start),
                            "end": int(end),
                            "strand": strand,
                            "source": source,
                            "attributes": attrs,
                        }
                    )
            elif feature_type in TRANSCRIPT_TYPES:
                transcript_id = first_attr(attrs, ["ID", "transcript_id", "Name"])
                gene_id = first_attr(attrs, ["Parent", "gene_id", "Locus_id", "locus_id"])
                if transcript_id:
                    if gene_id:
                        transcript_to_gene[transcript_id] = gene_id
                        if gene_id not in synthetic_genes:
                            synthetic_genes[gene_id] = {
                                "gene_id": gene_id,
                                "seqid": seqid,
                                "start": int(start),
                                "end": int(end),
                                "strand": strand,
                                "source": source,
                                "attributes": {"ID": gene_id, "derived_from": "transcript_span"},
                            }
                        else:
                            synthetic_genes[gene_id]["start"] = min(int(synthetic_genes[gene_id]["start"]), int(start))
                            synthetic_genes[gene_id]["end"] = max(int(synthetic_genes[gene_id]["end"]), int(end))
                    transcripts.append(
                        {
                            "transcript_id": transcript_id,
                            "gene_id": gene_id,
                            "seqid": seqid,
                            "start": int(start),
                            "end": int(end),
                            "strand": strand,
                            "source": source,
                            "attributes": attrs,
                        }
                    )
            elif feature_id and parent_id and parent_id in transcript_to_gene:
                continue
    if not genes and synthetic_genes:
        genes = sorted(synthetic_genes.values(), key=lambda item: (str(item["seqid"]), int(item["start"]), str(item["gene_id"])))
    return genes, transcripts, seqids, transcript_to_gene


def read_id_mapping(path: Optional[Path]) -> Tuple[Dict[str, str], Dict[str, str]]:
    protein_to_gene: Dict[str, str] = {}
    protein_to_transcript: Dict[str, str] = {}
    if not path:
        return protein_to_gene, protein_to_transcript
    if not path.exists():
        raise FileNotFoundError(f"id mapping not found: {path}")
    with open_text(path) as handle:
        first = handle.readline()
        handle.seek(0)
        if first.lower().startswith(("protein", "gene", "transcript", "original")):
            rows = read_optional_table(path)
            for item in rows:
                protein_id = item.get("protein_id") or item.get("Protein ID") or item.get("protein") or item.get("Original Protein ID")
                gene_id = item.get("gene_id") or item.get("Gene ID") or item.get("gene") or item.get("New Gene ID")
                transcript_id = item.get("transcript_id") or item.get("Transcript ID") or item.get("transcript") or item.get("New Transcript ID")
                if protein_id and gene_id:
                    protein_to_gene[protein_id] = gene_id
                if protein_id and transcript_id:
                    protein_to_transcript[protein_id] = transcript_id
            return protein_to_gene, protein_to_transcript

        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                continue
            source_id = fields[0].strip()
            for target_id in re.split(r"[,;\s]+", fields[1].strip()):
                if not target_id:
                    continue
                protein_to_gene[target_id] = source_id
                protein_to_transcript[target_id] = source_id
    return protein_to_gene, protein_to_transcript


def infer_protein_links(
    protein_id: str,
    description: str,
    gene_ids: set[str],
    transcript_ids: set[str],
    transcript_to_gene: Dict[str, str],
    protein_to_gene: Dict[str, str],
    protein_to_transcript: Dict[str, str],
) -> Tuple[Optional[str], Optional[str]]:
    attrs = parse_attributes(description)
    gene_id = protein_to_gene.get(protein_id) or first_attr(attrs, ["gene", "gene_id", "geneID", "locus"])
    transcript_id = protein_to_transcript.get(protein_id) or first_attr(attrs, ["transcript", "transcript_id", "transcriptID"])
    if not gene_id and protein_id in gene_ids:
        gene_id = protein_id
    if not transcript_id and protein_id in transcript_ids:
        transcript_id = protein_id
    if transcript_id and not gene_id:
        gene_id = transcript_to_gene.get(transcript_id)
    if not gene_id and "." in protein_id:
        prefix = protein_id.rsplit(".", 1)[0]
        if prefix in gene_ids:
            gene_id = prefix
    return gene_id, transcript_id


def annotation_aliases(value: Optional[str]) -> Iterator[str]:
    """Yield normalized annotation aliases while preserving the original spelling."""
    for item in unquote(value or "").split(","):
        alias = item.strip()
        if alias:
            yield alias


def collect_gene_aliases(
    transcripts: Iterable[Dict[str, Any]],
    protein_to_gene: Dict[str, str],
) -> set[Tuple[str, str, str]]:
    """Collect queryable gene aliases from reference annotations and RAP--MSU mappings."""
    aliases: set[Tuple[str, str, str]] = set()
    for transcript in transcripts:
        gene_id = transcript.get("gene_id")
        if not gene_id:
            continue
        attributes = transcript.get("attributes") or {}
        for key in GENE_SYMBOL_ATTRIBUTE_KEYS:
            for alias in annotation_aliases(attributes.get(key)):
                aliases.add((str(gene_id), alias, "annotation_symbol"))
    for mapped_id, gene_id in protein_to_gene.items():
        if not mapped_id or not gene_id or str(gene_id).casefold() == "none":
            continue
        aliases.add((str(gene_id), str(mapped_id), "msu_transcript"))
        if "." in mapped_id:
            aliases.add((str(gene_id), mapped_id.rsplit(".", 1)[0], "msu_gene"))
    return aliases


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


def command_version(command: str, args: List[str]) -> str:
    try:
        proc = subprocess.run([command, *args], check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    except OSError:
        return "unknown"
    return " ".join((proc.stdout or "").strip().splitlines()[:2]) or "unknown"


def run_command(cmd: List[str]) -> None:
    subprocess.run(cmd, check=True)


def diamond_ready(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 0


def blast_ready(prefix: Path, dbtype: str) -> bool:
    suffixes = [".phr", ".pin", ".psq"] if dbtype == "prot" else [".nhr", ".nin", ".nsq"]
    return all(Path(str(prefix) + suffix).exists() for suffix in suffixes)


def init_db(conn: sqlite3.Connection) -> None:
    conn.executescript(
        """
        PRAGMA foreign_keys = ON;

        CREATE TABLE IF NOT EXISTS reference_species (
            reference_id TEXT PRIMARY KEY,
            scientific_name TEXT NOT NULL,
            display_name TEXT NOT NULL,
            taxon_group TEXT NOT NULL,
            assembly_version TEXT,
            source_name TEXT,
            manifest_path TEXT NOT NULL,
            created_at TEXT NOT NULL,
            updated_at TEXT NOT NULL
        );

        CREATE TABLE IF NOT EXISTS reference_chromosomes (
            reference_id TEXT NOT NULL,
            seqid TEXT NOT NULL,
            length INTEGER NOT NULL,
            offset INTEGER,
            line_bases INTEGER,
            line_width INTEGER,
            PRIMARY KEY (reference_id, seqid),
            FOREIGN KEY (reference_id) REFERENCES reference_species(reference_id) ON DELETE CASCADE
        );

        CREATE TABLE IF NOT EXISTS reference_genes (
            reference_id TEXT NOT NULL,
            gene_id TEXT NOT NULL,
            seqid TEXT NOT NULL,
            start INTEGER NOT NULL,
            end INTEGER NOT NULL,
            strand TEXT,
            source TEXT,
            attributes_json TEXT,
            PRIMARY KEY (reference_id, gene_id),
            FOREIGN KEY (reference_id) REFERENCES reference_species(reference_id) ON DELETE CASCADE
        );

        CREATE TABLE IF NOT EXISTS reference_transcripts (
            reference_id TEXT NOT NULL,
            transcript_id TEXT NOT NULL,
            gene_id TEXT,
            seqid TEXT NOT NULL,
            start INTEGER NOT NULL,
            end INTEGER NOT NULL,
            strand TEXT,
            source TEXT,
            attributes_json TEXT,
            PRIMARY KEY (reference_id, transcript_id),
            FOREIGN KEY (reference_id) REFERENCES reference_species(reference_id) ON DELETE CASCADE
        );

        CREATE TABLE IF NOT EXISTS reference_proteins (
            reference_id TEXT NOT NULL,
            protein_id TEXT NOT NULL,
            gene_id TEXT,
            transcript_id TEXT,
            description TEXT,
            length INTEGER NOT NULL,
            sequence TEXT NOT NULL,
            PRIMARY KEY (reference_id, protein_id),
            FOREIGN KEY (reference_id) REFERENCES reference_species(reference_id) ON DELETE CASCADE
        );

        CREATE TABLE IF NOT EXISTS reference_gene_aliases (
            reference_id TEXT NOT NULL,
            gene_id TEXT NOT NULL,
            alias TEXT NOT NULL,
            alias_type TEXT NOT NULL,
            PRIMARY KEY (reference_id, gene_id, alias, alias_type),
            FOREIGN KEY (reference_id, gene_id) REFERENCES reference_genes(reference_id, gene_id) ON DELETE CASCADE
        );

        CREATE TABLE IF NOT EXISTS reference_assets (
            reference_id TEXT NOT NULL,
            asset_type TEXT NOT NULL,
            path TEXT NOT NULL,
            PRIMARY KEY (reference_id, asset_type),
            FOREIGN KEY (reference_id) REFERENCES reference_species(reference_id) ON DELETE CASCADE
        );

        CREATE TABLE IF NOT EXISTS reference_search_assets (
            reference_id TEXT NOT NULL,
            asset_type TEXT NOT NULL,
            source_fasta TEXT NOT NULL,
            index_prefix TEXT NOT NULL,
            tool TEXT NOT NULL,
            tool_version TEXT,
            created_at TEXT NOT NULL,
            status TEXT NOT NULL,
            metadata_json TEXT,
            PRIMARY KEY (reference_id, asset_type),
            FOREIGN KEY (reference_id) REFERENCES reference_species(reference_id) ON DELETE CASCADE
        );

        CREATE TABLE IF NOT EXISTS reference_go_annotations (
            reference_id TEXT NOT NULL,
            gene_id TEXT,
            protein_id TEXT,
            go_id TEXT NOT NULL,
            go_name TEXT,
            namespace TEXT,
            evidence TEXT,
            source TEXT,
            PRIMARY KEY (reference_id, gene_id, protein_id, go_id, evidence),
            FOREIGN KEY (reference_id) REFERENCES reference_species(reference_id) ON DELETE CASCADE
        );

        CREATE TABLE IF NOT EXISTS reference_pfam_domains (
            reference_id TEXT NOT NULL,
            protein_id TEXT NOT NULL,
            pfam_id TEXT NOT NULL,
            domain_name TEXT,
            start INTEGER,
            end INTEGER,
            evalue REAL,
            score REAL,
            source TEXT,
            PRIMARY KEY (reference_id, protein_id, pfam_id, start, end),
            FOREIGN KEY (reference_id) REFERENCES reference_species(reference_id) ON DELETE CASCADE
        );

        CREATE TABLE IF NOT EXISTS similarity_hits (
            species_id TEXT NOT NULL,
            gene_id TEXT NOT NULL,
            protein_id TEXT,
            reference_id TEXT NOT NULL,
            reference_protein_id TEXT NOT NULL,
            reference_gene_id TEXT,
            pident REAL,
            alignment_length INTEGER,
            evalue REAL,
            bitscore REAL,
            rank INTEGER NOT NULL,
            method TEXT NOT NULL,
            created_at TEXT NOT NULL,
            PRIMARY KEY (species_id, gene_id, reference_id, reference_protein_id, method),
            FOREIGN KEY (species_id) REFERENCES species(species_id) ON DELETE CASCADE,
            FOREIGN KEY (reference_id) REFERENCES reference_species(reference_id) ON DELETE CASCADE
        );

        CREATE INDEX IF NOT EXISTS idx_reference_proteins_gene ON reference_proteins(reference_id, gene_id);
        CREATE INDEX IF NOT EXISTS idx_reference_gene_aliases_lookup ON reference_gene_aliases(reference_id, alias COLLATE NOCASE);
        CREATE INDEX IF NOT EXISTS idx_reference_go_gene ON reference_go_annotations(reference_id, gene_id);
        CREATE INDEX IF NOT EXISTS idx_reference_pfam_protein ON reference_pfam_domains(reference_id, protein_id);
        CREATE INDEX IF NOT EXISTS idx_similarity_hits_species_gene ON similarity_hits(species_id, gene_id);
        """
    )


def read_optional_table(path: Optional[Path]) -> List[Dict[str, str]]:
    if not path:
        return []
    if not path.exists():
        raise FileNotFoundError(f"optional table not found: {path}")
    with open_text(path) as handle:
        sample = handle.read(4096)
        handle.seek(0)
        dialect = csv.excel_tab if sample.count("\t") >= sample.count(",") else csv.excel
        return list(csv.DictReader(handle, dialect=dialect))


def parse_go_cell(value: str) -> Iterator[Tuple[str, Optional[str], Optional[str]]]:
    for match in re.finditer(r"([^:,]+(?: [^:,]+)*):\s*([^()]+)\s*\((GO:\d+)\)", value):
        namespace = match.group(1).strip()
        go_name = match.group(2).strip()
        go_id = match.group(3).strip()
        yield go_id, go_name, namespace


def parse_go_cell_without_ids(value: str) -> Iterator[Tuple[str, Optional[str], Optional[str]]]:
    for go_id in re.findall(r"GO:\d+", value):
        yield go_id, None, None


def namespace_from_gaf_aspect(aspect: str) -> Optional[str]:
    return {
        "F": "molecular_function",
        "P": "biological_process",
        "C": "cellular_component",
    }.get(aspect)


def insert_gaf_go(conn: sqlite3.Connection, reference_id: str, path: Path) -> int:
    inserted = 0
    with open_text(path) as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("!"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 15:
                continue
            gene_id = fields[1].strip() or None
            go_id = fields[4].strip()
            if not gene_id or not go_id:
                continue
            conn.execute(
                """
                INSERT OR REPLACE INTO reference_go_annotations
                (reference_id, gene_id, protein_id, go_id, go_name, namespace, evidence, source)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    reference_id,
                    gene_id,
                    None,
                    go_id,
                    None,
                    namespace_from_gaf_aspect(fields[8].strip()),
                    fields[6].strip() or None,
                    fields[14].strip() or str(path),
                ),
            )
            inserted += 1
    return inserted


def insert_optional_go(conn: sqlite3.Connection, reference_id: str, path: Optional[Path]) -> int:
    if path and path.name.endswith((".gaf", ".gaf.gz")):
        return insert_gaf_go(conn, reference_id, path)
    rows = read_optional_table(path)
    inserted = 0
    for item in rows:
        gene_id = item.get("gene_id") or item.get("Gene_ID") or item.get("Locus_ID") or None
        protein_id = item.get("protein_id") or item.get("Protein_ID") or item.get("Transcript_ID") or None
        explicit_go = item.get("go_id") or item.get("GO_ID")
        go_terms = []
        if explicit_go:
            go_terms.append((explicit_go, item.get("go_name") or item.get("name") or None, item.get("namespace") or item.get("aspect") or None))
        elif item.get("GO"):
            go_terms = list(parse_go_cell(item["GO"]))
            if not go_terms:
                go_terms = list(parse_go_cell_without_ids(item["GO"]))
        for go_id, go_name, namespace in go_terms:
            conn.execute(
                """
                INSERT OR REPLACE INTO reference_go_annotations
                (reference_id, gene_id, protein_id, go_id, go_name, namespace, evidence, source)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (reference_id, gene_id, protein_id, go_id, go_name, namespace, item.get("evidence") or item.get("Evidence") or None, str(path) if path else None),
            )
            inserted += 1
    return inserted


def insert_optional_pfam(conn: sqlite3.Connection, reference_id: str, path: Optional[Path]) -> int:
    rows = read_optional_table(path)
    for item in rows:
        start = item.get("start") or item.get("ali_from")
        end = item.get("end") or item.get("ali_to")
        conn.execute(
            """
            INSERT OR REPLACE INTO reference_pfam_domains
            (reference_id, protein_id, pfam_id, domain_name, start, end, evalue, score, source)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                reference_id,
                item.get("protein_id") or item.get("target_name"),
                item.get("pfam_id") or item.get("accession") or item.get("hmm_acc"),
                item.get("domain_name") or item.get("hmm_name") or None,
                int(start) if start else None,
                int(end) if end else None,
                float(item["evalue"]) if item.get("evalue") else None,
                float(item["score"]) if item.get("score") else None,
                str(path) if path else None,
            ),
        )
    return len(rows)


def register_search_asset(conn: sqlite3.Connection, reference_id: str, asset_type: str, source_fasta: Path, index_prefix: Path, tool: str, version: str, metadata: Dict[str, Any]) -> None:
    conn.execute(
        """
        INSERT OR REPLACE INTO reference_search_assets
        (reference_id, asset_type, source_fasta, index_prefix, tool, tool_version, created_at, status, metadata_json)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (reference_id, asset_type, str(source_fasta), str(index_prefix), tool, version, utc_now(), "ready", json.dumps(metadata, sort_keys=True)),
    )


def build_reference(args: argparse.Namespace) -> Dict[str, Any]:
    manifest_path = Path(args.manifest).resolve()
    manifest = parse_simple_yaml(manifest_path)
    reference_id = str(require(manifest, "reference_id"))
    taxon_group = str(require(manifest, "taxon_group"))
    if taxon_group not in ALLOWED_TAXON_GROUPS:
        raise ValueError(f"taxon_group must be one of {', '.join(sorted(ALLOWED_TAXON_GROUPS))}")

    genome_fasta = resolve_path(manifest_path, str(require(manifest, "genome.fasta")))
    annotation_file = resolve_path(manifest_path, str(require(manifest, "annotation.file")))
    protein_fasta = resolve_path(manifest_path, str(require(manifest, "protein.fasta")))
    optional = manifest.get("optional", {})
    id_mapping_path = resolve_path(manifest_path, optional.get("id_mapping") or optional.get("rap_msu_mapping"))
    go_path = resolve_path(manifest_path, optional.get("go_annotation"))
    pfam_path = resolve_path(manifest_path, optional.get("pfam_annotation"))
    for label, path in [("genome FASTA", genome_fasta), ("annotation", annotation_file), ("protein FASTA", protein_fasta)]:
        if not path or not path.exists():
            raise FileNotFoundError(f"{label} not found: {path}")

    assets_dir = Path(args.assets_dir).resolve() / reference_id
    fai_path = create_fai(genome_fasta, assets_dir / "genome" / f"{reference_id}.genome.fa.fai")
    chromosomes = parse_fai(fai_path)
    chromosome_ids = {item["seqid"] for item in chromosomes}
    genes, transcripts, annotation_seqids, transcript_to_gene = parse_gene_models(annotation_file)
    gene_ids = {item["gene_id"] for item in genes}
    transcript_ids = {item["transcript_id"] for item in transcripts}
    protein_to_gene, protein_to_transcript = read_id_mapping(id_mapping_path)

    proteins = []
    linked_to_gene = 0
    linked_to_transcript = 0
    for protein_id, description, sequence in fasta_records(protein_fasta):
        gene_id, transcript_id = infer_protein_links(protein_id, description, gene_ids, transcript_ids, transcript_to_gene, protein_to_gene, protein_to_transcript)
        linked_to_gene += 1 if gene_id else 0
        linked_to_transcript += 1 if transcript_id else 0
        proteins.append(
            {
                "protein_id": protein_id,
                "gene_id": gene_id,
                "transcript_id": transcript_id,
                "description": description,
                "length": len(sequence),
                "sequence": sequence,
            }
        )
    gene_aliases = {
        item for item in collect_gene_aliases(transcripts, protein_to_gene)
        if item[0] in gene_ids
    }

    summary: Dict[str, Any] = {
        "reference_id": reference_id,
        "taxon_group": taxon_group,
        "chromosomes": len(chromosomes),
        "genes": len(genes),
        "transcripts": len(transcripts),
        "proteins": len(proteins),
        "proteins_linked_to_genes": linked_to_gene,
        "proteins_linked_to_transcripts": linked_to_transcript,
        "id_mapping_gene_links": len(protein_to_gene),
        "id_mapping_transcript_links": len(protein_to_transcript),
        "gene_aliases": len(gene_aliases),
        "annotation_seqids_missing_from_fasta": len(annotation_seqids - chromosome_ids),
        "created": [],
        "reused": [],
        "skipped": [],
        "warnings": [],
    }
    if annotation_seqids - chromosome_ids:
        summary["warnings"].append("some annotation seqids are missing from genome FASTA index")
    if proteins and linked_to_gene == 0:
        summary["warnings"].append("protein FASTA headers did not expose gene links; GO/similarity mapping may need an ID mapping table later")
    if args.check_only:
        return summary

    db_path = Path(args.db).resolve()
    db_path.parent.mkdir(parents=True, exist_ok=True)
    diamond = tool_path("DIAMOND", str(ORTHOLOG_BIN / "diamond"), "diamond")
    makeblastdb = tool_path("MAKEBLASTDB", str(ORTHOLOG_BIN / "makeblastdb"), "makeblastdb")

    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA foreign_keys = ON")
    init_db(conn)
    try:
        if args.replace:
            conn.execute("DELETE FROM reference_species WHERE reference_id = ?", (reference_id,))
        elif conn.execute("SELECT 1 FROM reference_species WHERE reference_id = ?", (reference_id,)).fetchone():
            raise ValueError(f"reference already exists: {reference_id}; use --replace")

        now = utc_now()
        conn.execute(
            """
            INSERT INTO reference_species
            (reference_id, scientific_name, display_name, taxon_group, assembly_version, source_name, manifest_path, created_at, updated_at)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                reference_id,
                str(require(manifest, "scientific_name")),
                str(require(manifest, "display_name")),
                taxon_group,
                manifest.get("assembly_version"),
                manifest.get("source_name"),
                str(manifest_path),
                now,
                now,
            ),
        )
        for asset_type, path in [("fasta", genome_fasta), ("fai", fai_path), ("annotation", annotation_file), ("protein_fasta", protein_fasta)]:
            conn.execute("INSERT INTO reference_assets VALUES (?, ?, ?)", (reference_id, asset_type, str(path)))
        if id_mapping_path:
            conn.execute("INSERT INTO reference_assets VALUES (?, ?, ?)", (reference_id, "id_mapping", str(id_mapping_path)))
        for item in chromosomes:
            conn.execute(
                "INSERT INTO reference_chromosomes VALUES (?, ?, ?, ?, ?, ?)",
                (reference_id, item["seqid"], item["length"], item["offset"], item["line_bases"], item["line_width"]),
            )
        for item in genes:
            conn.execute(
                "INSERT INTO reference_genes VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                (reference_id, item["gene_id"], item["seqid"], item["start"], item["end"], item["strand"], item["source"], json.dumps(item["attributes"], sort_keys=True)),
            )
        for item in transcripts:
            conn.execute(
                "INSERT INTO reference_transcripts VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (reference_id, item["transcript_id"], item["gene_id"], item["seqid"], item["start"], item["end"], item["strand"], item["source"], json.dumps(item["attributes"], sort_keys=True)),
            )
        for item in proteins:
            conn.execute(
                "INSERT INTO reference_proteins VALUES (?, ?, ?, ?, ?, ?, ?)",
                (reference_id, item["protein_id"], item["gene_id"], item["transcript_id"], item["description"], item["length"], item["sequence"]),
            )
        conn.executemany(
            "INSERT INTO reference_gene_aliases VALUES (?, ?, ?, ?)",
            [(reference_id, gene_id, alias, alias_type) for gene_id, alias, alias_type in sorted(gene_aliases)],
        )
        summary["go_annotations"] = insert_optional_go(conn, reference_id, go_path)
        summary["pfam_domains"] = insert_optional_pfam(conn, reference_id, pfam_path)

        if args.skip_indexes:
            summary["skipped"].append("search_indexes")
        else:
            search_dir = assets_dir / "search"
            if diamond:
                out_prefix = search_dir / "diamond" / "protein"
                out_dmnd = Path(str(out_prefix) + ".dmnd")
                out_prefix.parent.mkdir(parents=True, exist_ok=True)
                if args.force_index or not diamond_ready(out_dmnd):
                    run_command([diamond, "makedb", "--in", str(protein_fasta), "-d", str(out_prefix), "--threads", str(args.threads)])
                    summary["created"].append("diamond_protein")
                else:
                    summary["reused"].append("diamond_protein")
                register_search_asset(conn, reference_id, "diamond_protein", protein_fasta, out_dmnd, "diamond", command_version(diamond, ["version"]), {"threads": args.threads})
            else:
                summary["skipped"].append("diamond_protein")
                summary["warnings"].append("diamond not found")
            if makeblastdb:
                for asset_type, fasta_path, dbtype, suffix in [
                    ("blast_protein", protein_fasta, "prot", "protein"),
                    ("blast_genome", genome_fasta, "nucl", "genome"),
                ]:
                    out_prefix = search_dir / "blast" / suffix
                    out_prefix.parent.mkdir(parents=True, exist_ok=True)
                    if args.force_index or not blast_ready(out_prefix, dbtype):
                        run_command([makeblastdb, "-in", str(fasta_path), "-dbtype", dbtype, "-out", str(out_prefix)])
                        summary["created"].append(asset_type)
                    else:
                        summary["reused"].append(asset_type)
                    register_search_asset(conn, reference_id, asset_type, fasta_path, out_prefix, "makeblastdb", command_version(makeblastdb, ["-version"]), {"dbtype": dbtype})
            else:
                summary["skipped"].extend(["blast_protein", "blast_genome"])
                summary["warnings"].append("makeblastdb not found")
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
        summary = build_reference(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
