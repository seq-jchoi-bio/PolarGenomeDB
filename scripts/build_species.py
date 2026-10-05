#!/usr/bin/env python3
"""Build or append one species into the AntarcticPlantDB SQLite database."""

from __future__ import annotations

import argparse
import datetime as dt
import gzip
import json
import os
import re
import sqlite3
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Optional, Tuple
from zipfile import ZipFile


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = ROOT / "data" / "antarctic_plant.db"
DEFAULT_ASSETS_DIR = ROOT / "data" / "species_assets"
ALLOWED_TAXON_GROUPS = {"monocot", "dicot", "other", "unknown"}
FEATURE_TYPES = {"exon", "CDS", "cds", "five_prime_UTR", "three_prime_UTR", "UTR", "start_codon", "stop_codon"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, help="Species manifest YAML file.")
    parser.add_argument("--db", default=str(DEFAULT_DB), help="SQLite DB path.")
    parser.add_argument("--assets-dir", default=str(DEFAULT_ASSETS_DIR), help="Runtime asset output directory.")
    parser.add_argument("--reports-dir", default=str(ROOT / "reports"), help="Build report output directory.")
    parser.add_argument("--replace", action="store_true", help="Replace this species if it already exists.")
    return parser.parse_args()


def parse_simple_yaml(path: Path) -> Dict[str, Any]:
    """Parse the small manifest YAML subset used by this project."""
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
        key = key.strip()
        value = value.strip()
        while stack and indent <= stack[-1][0]:
            stack.pop()
        parent = stack[-1][1]
        if value == "":
            child: Dict[str, Any] = {}
            parent[key] = child
            stack.append((indent, child))
        else:
            parent[key] = value.strip("\"'")
    return root


def resolve_path(manifest_path: Path, value: Optional[str]) -> Optional[Path]:
    if not value:
        return None
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = (manifest_path.parent / path).resolve()
    return path


def require(mapping: Dict[str, Any], dotted: str) -> Any:
    current: Any = mapping
    for part in dotted.split("."):
        if not isinstance(current, dict) or part not in current:
            raise ValueError(f"manifest missing required key: {dotted}")
        current = current[part]
    return current


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def init_db(conn: sqlite3.Connection) -> None:
    conn.executescript(
        """
        PRAGMA foreign_keys = ON;

        CREATE TABLE IF NOT EXISTS species (
            species_id TEXT PRIMARY KEY,
            scientific_name TEXT NOT NULL,
            display_name TEXT NOT NULL,
            taxon_group TEXT NOT NULL,
            assembly_version TEXT,
            manifest_path TEXT NOT NULL,
            created_at TEXT NOT NULL,
            updated_at TEXT NOT NULL
        );

        CREATE TABLE IF NOT EXISTS chromosomes (
            species_id TEXT NOT NULL,
            seqid TEXT NOT NULL,
            length INTEGER NOT NULL,
            offset INTEGER,
            line_bases INTEGER,
            line_width INTEGER,
            PRIMARY KEY (species_id, seqid),
            FOREIGN KEY (species_id) REFERENCES species(species_id) ON DELETE CASCADE
        );

        CREATE TABLE IF NOT EXISTS genes (
            species_id TEXT NOT NULL,
            gene_id TEXT NOT NULL,
            seqid TEXT NOT NULL,
            start INTEGER NOT NULL,
            end INTEGER NOT NULL,
            strand TEXT,
            source TEXT,
            attributes_json TEXT,
            PRIMARY KEY (species_id, gene_id),
            FOREIGN KEY (species_id) REFERENCES species(species_id) ON DELETE CASCADE
        );

        CREATE TABLE IF NOT EXISTS transcripts (
            species_id TEXT NOT NULL,
            transcript_id TEXT NOT NULL,
            gene_id TEXT,
            seqid TEXT NOT NULL,
            start INTEGER NOT NULL,
            end INTEGER NOT NULL,
            strand TEXT,
            source TEXT,
            attributes_json TEXT,
            PRIMARY KEY (species_id, transcript_id),
            FOREIGN KEY (species_id) REFERENCES species(species_id) ON DELETE CASCADE
        );

        CREATE TABLE IF NOT EXISTS features (
            species_id TEXT NOT NULL,
            feature_id TEXT NOT NULL,
            feature_type TEXT NOT NULL,
            parent_id TEXT,
            gene_id TEXT,
            transcript_id TEXT,
            seqid TEXT NOT NULL,
            start INTEGER NOT NULL,
            end INTEGER NOT NULL,
            strand TEXT,
            phase TEXT,
            attributes_json TEXT,
            PRIMARY KEY (species_id, feature_id),
            FOREIGN KEY (species_id) REFERENCES species(species_id) ON DELETE CASCADE
        );

        CREATE TABLE IF NOT EXISTS proteins (
            species_id TEXT NOT NULL,
            protein_id TEXT NOT NULL,
            gene_id TEXT,
            transcript_id TEXT,
            description TEXT,
            length INTEGER NOT NULL,
            sequence TEXT NOT NULL,
            PRIMARY KEY (species_id, protein_id),
            FOREIGN KEY (species_id) REFERENCES species(species_id) ON DELETE CASCADE
        );

        CREATE TABLE IF NOT EXISTS genome_assets (
            species_id TEXT NOT NULL,
            asset_type TEXT NOT NULL,
            path TEXT NOT NULL,
            PRIMARY KEY (species_id, asset_type),
            FOREIGN KEY (species_id) REFERENCES species(species_id) ON DELETE CASCADE
        );

        CREATE TABLE IF NOT EXISTS annotation_assets (
            species_id TEXT NOT NULL,
            format TEXT NOT NULL,
            path TEXT NOT NULL,
            PRIMARY KEY (species_id, format, path),
            FOREIGN KEY (species_id) REFERENCES species(species_id) ON DELETE CASCADE
        );

        CREATE TABLE IF NOT EXISTS optional_assets (
            species_id TEXT NOT NULL,
            asset_key TEXT NOT NULL,
            path TEXT NOT NULL,
            exists_on_disk INTEGER NOT NULL,
            PRIMARY KEY (species_id, asset_key),
            FOREIGN KEY (species_id) REFERENCES species(species_id) ON DELETE CASCADE
        );

        CREATE TABLE IF NOT EXISTS id_mappings (
            species_id TEXT NOT NULL,
            source_id TEXT NOT NULL,
            source_type TEXT NOT NULL,
            target_id TEXT NOT NULL,
            target_type TEXT NOT NULL,
            source_file TEXT NOT NULL,
            PRIMARY KEY (species_id, source_id, source_type, target_id, target_type),
            FOREIGN KEY (species_id) REFERENCES species(species_id) ON DELETE CASCADE
        );

        CREATE TABLE IF NOT EXISTS build_runs (
            build_id INTEGER PRIMARY KEY AUTOINCREMENT,
            species_id TEXT NOT NULL,
            manifest_path TEXT NOT NULL,
            started_at TEXT NOT NULL,
            finished_at TEXT,
            status TEXT NOT NULL,
            summary_json TEXT,
            FOREIGN KEY (species_id) REFERENCES species(species_id) ON DELETE CASCADE
        );

        CREATE INDEX IF NOT EXISTS idx_genes_species_seqid ON genes(species_id, seqid, start, end);
        CREATE INDEX IF NOT EXISTS idx_transcripts_species_gene ON transcripts(species_id, gene_id);
        CREATE INDEX IF NOT EXISTS idx_features_species_parent ON features(species_id, parent_id);
        CREATE INDEX IF NOT EXISTS idx_features_species_gene ON features(species_id, gene_id, start, end);
        CREATE INDEX IF NOT EXISTS idx_features_species_transcript ON features(species_id, transcript_id);
        CREATE INDEX IF NOT EXISTS idx_features_species_seqid ON features(species_id, seqid, start, end);
        CREATE INDEX IF NOT EXISTS idx_proteins_species_gene ON proteins(species_id, gene_id);
        CREATE INDEX IF NOT EXISTS idx_proteins_species_transcript ON proteins(species_id, transcript_id);
        """
    )


def delete_species(conn: sqlite3.Connection, species_id: str) -> None:
    conn.execute("DELETE FROM species WHERE species_id = ?", (species_id,))


def create_fai(fasta_path: Path, fai_path: Path) -> Path:
    fai_path.parent.mkdir(parents=True, exist_ok=True)
    if fai_path.exists() and fai_path.stat().st_mtime >= fasta_path.stat().st_mtime:
        return fai_path

    with fasta_path.open("rb") as src, fai_path.open("w", encoding="ascii") as out:
        name: Optional[str] = None
        length = 0
        offset = 0
        line_bases = 0
        line_width = 0
        while True:
            pos = src.tell()
            line = src.readline()
            if not line:
                break
            if line.startswith(b">"):
                if name is not None:
                    out.write(f"{name}\t{length}\t{offset}\t{line_bases}\t{line_width}\n")
                name = line[1:].strip().split()[0].decode("utf-8", errors="replace")
                length = 0
                offset = src.tell()
                line_bases = 0
                line_width = 0
            else:
                stripped = line.rstrip(b"\r\n")
                if name is None or not stripped:
                    continue
                if line_bases == 0:
                    line_bases = len(stripped)
                    line_width = len(line)
                length += len(stripped)
        if name is not None:
            out.write(f"{name}\t{length}\t{offset}\t{line_bases}\t{line_width}\n")
    return fai_path


def read_fai(fai_path: Path) -> List[Tuple[str, int, int, int, int]]:
    rows = []
    with fai_path.open("rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 5:
                rows.append((parts[0], int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4])))
    return rows


def parse_gff3_attrs(text: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for item in text.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key.strip()] = value.strip()
    return attrs


def parse_gtf_attrs(text: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for key, value in re.findall(r'(\S+)\s+"([^"]*)"', text):
        attrs[key] = value
    return attrs


def gene_id_from_transcript_id(transcript_id: str) -> Optional[str]:
    match = re.match(r"^(.+)\.\d+$", transcript_id.strip())
    return match.group(1) if match else None


def infer_gene_id(attrs: Dict[str, str], raw_attr: str, feature_type: str, fmt: str) -> Optional[str]:
    if fmt == "gff3":
        if feature_type == "gene":
            return attrs.get("ID") or attrs.get("gene_id") or attrs.get("Name")
        return attrs.get("gene_id") or attrs.get("gene") or attrs.get("Parent")
    if attrs.get("gene_id"):
        return attrs["gene_id"]
    raw_id = raw_attr.strip()
    if feature_type == "gene" and raw_id:
        return raw_id
    if feature_type == "transcript" and raw_id:
        return gene_id_from_transcript_id(raw_id)
    return None


def infer_transcript_id(attrs: Dict[str, str], raw_attr: str, feature_type: str, fmt: str) -> Optional[str]:
    if fmt == "gff3":
        if feature_type in {"transcript", "mRNA"}:
            return attrs.get("ID") or attrs.get("transcript_id")
        return attrs.get("transcript_id") or attrs.get("Parent")
    return attrs.get("transcript_id") or (raw_attr.strip() if feature_type == "transcript" and raw_attr.strip() else None)


def parse_annotation(path: Path, fmt: str) -> Iterator[Dict[str, Any]]:
    with open_text(path) as handle:
        for line_no, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            seqid, source, feature_type, start, end, score, strand, phase, raw_attr = parts
            attrs = parse_gff3_attrs(raw_attr) if fmt == "gff3" else parse_gtf_attrs(raw_attr)
            yield {
                "line_no": line_no,
                "seqid": seqid,
                "source": source,
                "type": feature_type,
                "start": int(start),
                "end": int(end),
                "strand": strand,
                "phase": phase,
                "raw_attr": raw_attr,
                "attrs": attrs,
                "gene_id": infer_gene_id(attrs, raw_attr, feature_type, fmt),
                "transcript_id": infer_transcript_id(attrs, raw_attr, feature_type, fmt),
            }


def parse_proteins(path: Path) -> Iterator[Tuple[str, Optional[str], Optional[str], str, str]]:
    protein_id: Optional[str] = None
    description = ""
    chunks: List[str] = []

    def emit() -> Optional[Tuple[str, Optional[str], Optional[str], str, str]]:
        if protein_id is None:
            return None
        desc = description
        gene_match = re.search(r"(?:^|\s)gene=([^\s]+)", desc)
        tx_match = re.search(r"(?:^|\s)(?:transcript|transcript_id)=([^\s]+)", desc)
        transcript_id = tx_match.group(1) if tx_match else None
        gene_id = gene_match.group(1) if gene_match else None
        return protein_id, gene_id, transcript_id, desc, "".join(chunks)

    with open_text(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                row = emit()
                if row:
                    yield row
                header = line[1:].strip()
                protein_id = header.split()[0]
                description = header
                chunks = []
            else:
                chunks.append(line.strip())
        row = emit()
        if row:
            yield row


def xlsx_first_sheet_rows(path: Path) -> Iterator[List[str]]:
    ns = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with ZipFile(path) as zf:
        with zf.open("xl/worksheets/sheet1.xml") as handle:
            for event, elem in ET.iterparse(handle, events=("end",)):
                if not elem.tag.endswith("row"):
                    continue
                values = []
                for cell in elem.findall("x:c", ns):
                    value = cell.find("x:v", ns)
                    inline = cell.find("x:is", ns)
                    if value is not None:
                        values.append(value.text or "")
                    elif inline is not None:
                        texts = [text.text or "" for text in inline.findall(".//x:t", ns)]
                        values.append("".join(texts))
                    else:
                        values.append("")
                yield values
                elem.clear()


def delimited_rows(path: Path) -> Iterator[List[str]]:
    delimiter = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            yield [part.strip() for part in line.rstrip("\n").split(delimiter)]


def read_id_mappings(path: Optional[Path]) -> Tuple[Dict[str, str], Dict[str, str], List[Tuple[str, str, str, str]]]:
    protein_to_gene: Dict[str, str] = {}
    protein_to_transcript: Dict[str, str] = {}
    mapping_rows: List[Tuple[str, str, str, str]] = []
    if not path or not path.exists():
        return protein_to_gene, protein_to_transcript, mapping_rows

    rows = xlsx_first_sheet_rows(path) if path.suffix.lower() == ".xlsx" else delimited_rows(path)
    header: Optional[List[str]] = None
    for row in rows:
        normalized = [cell.strip() for cell in row]
        if not any(normalized):
            continue
        lower = [cell.lower() for cell in normalized]
        if any("gene" in cell for cell in lower) and any("transcript" in cell for cell in lower):
            header = normalized
            continue
        if header is None:
            continue
        record = {header[i]: normalized[i] for i in range(min(len(header), len(normalized)))}

        original_protein = record.get("Original Protein ID") or record.get("protein_id") or record.get("Protein ID")
        new_protein = record.get("New Protein ID")
        new_gene = record.get("New Gene ID") or record.get("new_gene_ID") or record.get("gene_id")
        new_transcript = record.get("New Transcript ID") or record.get("new_transcript_ID") or record.get("transcript_id")

        if original_protein and new_gene:
            protein_to_gene[original_protein] = new_gene
            mapping_rows.append((original_protein, "protein", new_gene, "gene"))
        if original_protein and new_transcript:
            protein_to_transcript[original_protein] = new_transcript
            mapping_rows.append((original_protein, "protein", new_transcript, "transcript"))
        if new_protein and new_gene:
            protein_to_gene[new_protein] = new_gene
            mapping_rows.append((new_protein, "protein", new_gene, "gene"))
        if new_protein and new_transcript:
            protein_to_transcript[new_protein] = new_transcript
            mapping_rows.append((new_protein, "protein", new_transcript, "transcript"))
        if new_transcript and new_gene:
            mapping_rows.append((new_transcript, "transcript", new_gene, "gene"))
    return protein_to_gene, protein_to_transcript, mapping_rows


def build_species(manifest_path: Path, db_path: Path, reports_dir: Path, assets_dir: Path, replace: bool) -> Dict[str, Any]:
    manifest = parse_simple_yaml(manifest_path)
    species_id = str(require(manifest, "species_id"))
    scientific_name = str(require(manifest, "scientific_name"))
    display_name = str(require(manifest, "display_name"))
    taxon_group = str(require(manifest, "taxon_group")).lower()
    if taxon_group not in ALLOWED_TAXON_GROUPS:
        raise ValueError(f"taxon_group must be one of {sorted(ALLOWED_TAXON_GROUPS)}")

    genome_fasta = resolve_path(manifest_path, require(manifest, "genome.fasta"))
    annotation_file = resolve_path(manifest_path, require(manifest, "annotation.file"))
    protein_fasta = resolve_path(manifest_path, require(manifest, "protein.fasta"))
    annotation_format = str(require(manifest, "annotation.format")).lower()
    if annotation_format == "gff":
        annotation_format = "gff3"
    if annotation_format not in {"gff3", "gtf"}:
        raise ValueError("annotation.format must be gff, gff3, or gtf")

    required_paths = {
        "genome.fasta": genome_fasta,
        "annotation.file": annotation_file,
        "protein.fasta": protein_fasta,
    }
    missing = [key for key, path in required_paths.items() if path is None or not path.exists()]
    if missing:
        raise FileNotFoundError(f"missing required input(s): {', '.join(missing)}")

    now = dt.datetime.now(dt.timezone.utc).isoformat()
    db_path.parent.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)
    species_assets_dir = assets_dir / species_id
    genome_assets_dir = species_assets_dir / "genome"

    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA foreign_keys = ON")
    init_db(conn)
    if conn.execute("SELECT 1 FROM species WHERE species_id = ?", (species_id,)).fetchone():
        if not replace:
            raise ValueError(f"species {species_id!r} already exists; pass --replace to rebuild it")
        delete_species(conn, species_id)
    conn.commit()

    build_id = None
    summary: Dict[str, Any] = {
        "species_id": species_id,
        "required_inputs": {key: str(path) for key, path in required_paths.items()},
        "annotation_format": annotation_format,
        "warnings": [],
    }
    try:
        conn.execute("BEGIN")
        conn.execute(
            """
            INSERT INTO species
            (species_id, scientific_name, display_name, taxon_group, assembly_version, manifest_path, created_at, updated_at)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                species_id,
                scientific_name,
                display_name,
                taxon_group,
                manifest.get("assembly_version"),
                str(manifest_path.resolve()),
                now,
                now,
            ),
        )
        conn.execute(
            "INSERT INTO build_runs (species_id, manifest_path, started_at, status) VALUES (?, ?, ?, ?)",
            (species_id, str(manifest_path.resolve()), now, "running"),
        )
        build_id = conn.execute("SELECT last_insert_rowid()").fetchone()[0]

        fai_path = create_fai(genome_fasta, genome_assets_dir / f"{species_id}.genome.fa.fai")
        chromosomes = read_fai(fai_path)
        conn.executemany(
            "INSERT INTO chromosomes (species_id, seqid, length, offset, line_bases, line_width) VALUES (?, ?, ?, ?, ?, ?)",
            [(species_id, *row) for row in chromosomes],
        )
        conn.execute("INSERT INTO genome_assets VALUES (?, ?, ?)", (species_id, "fasta", str(genome_fasta)))
        conn.execute("INSERT INTO genome_assets VALUES (?, ?, ?)", (species_id, "fai", str(fai_path)))
        masked = resolve_path(manifest_path, manifest.get("genome", {}).get("masked_fasta"))
        if masked:
            conn.execute("INSERT INTO genome_assets VALUES (?, ?, ?)", (species_id, "masked_fasta", str(masked)))
        conn.execute("INSERT INTO annotation_assets VALUES (?, ?, ?)", (species_id, annotation_format, str(annotation_file)))

        seqids = {row[0] for row in chromosomes}
        annotation_seqids = set()
        gene_count = transcript_count = feature_count = skipped_no_id = 0
        known_gene_ids = set()
        known_transcript_ids = set()
        seen_feature_counts: Dict[str, int] = {}
        for rec in parse_annotation(annotation_file, annotation_format):
            annotation_seqids.add(rec["seqid"])
            attrs_json = json.dumps(rec["attrs"], ensure_ascii=False, sort_keys=True)
            ftype = rec["type"]
            if ftype == "gene":
                gene_id = rec["gene_id"]
                if not gene_id:
                    skipped_no_id += 1
                    continue
                conn.execute(
                    """
                    INSERT OR REPLACE INTO genes
                    (species_id, gene_id, seqid, start, end, strand, source, attributes_json)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (species_id, gene_id, rec["seqid"], rec["start"], rec["end"], rec["strand"], rec["source"], attrs_json),
                )
                gene_count += 1
                known_gene_ids.add(gene_id)
            elif ftype in {"transcript", "mRNA"}:
                transcript_id = rec["transcript_id"]
                if not transcript_id:
                    skipped_no_id += 1
                    continue
                conn.execute(
                    """
                    INSERT OR REPLACE INTO transcripts
                    (species_id, transcript_id, gene_id, seqid, start, end, strand, source, attributes_json)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        species_id,
                        transcript_id,
                        rec["gene_id"],
                        rec["seqid"],
                        rec["start"],
                        rec["end"],
                        rec["strand"],
                        rec["source"],
                        attrs_json,
                    ),
                )
                transcript_count += 1
                known_transcript_ids.add(transcript_id)
            elif ftype in FEATURE_TYPES:
                base_id = rec["attrs"].get("ID") or f"{ftype}:{rec['seqid']}:{rec['start']}:{rec['end']}:{rec['line_no']}"
                idx = seen_feature_counts.get(base_id, 0)
                seen_feature_counts[base_id] = idx + 1
                feature_id = base_id if idx == 0 else f"{base_id}#{idx + 1}"
                parent_id = rec["attrs"].get("Parent") or rec["transcript_id"] or rec["gene_id"]
                conn.execute(
                    """
                    INSERT INTO features
                    (species_id, feature_id, feature_type, parent_id, gene_id, transcript_id, seqid, start, end, strand, phase, attributes_json)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        species_id,
                        feature_id,
                        ftype,
                        parent_id,
                        rec["gene_id"],
                        rec["transcript_id"],
                        rec["seqid"],
                        rec["start"],
                        rec["end"],
                        rec["strand"],
                        rec["phase"],
                        attrs_json,
                    ),
                )
                feature_count += 1

        optional = manifest.get("optional", {})
        mapping_path = None
        if isinstance(optional, dict):
            mapping_path = resolve_path(manifest_path, optional.get("id_mapping"))
        protein_gene_map, protein_transcript_map, id_mapping_rows = read_id_mappings(mapping_path)
        if id_mapping_rows and mapping_path:
            conn.executemany(
                """
                INSERT OR IGNORE INTO id_mappings
                (species_id, source_id, source_type, target_id, target_type, source_file)
                VALUES (?, ?, ?, ?, ?, ?)
                """,
                [(species_id, source_id, source_type, target_id, target_type, str(mapping_path)) for source_id, source_type, target_id, target_type in id_mapping_rows],
            )

        protein_count = protein_gene_links = protein_transcript_links = 0
        for protein_id, gene_id, transcript_id, description, sequence in parse_proteins(protein_fasta):
            gene_id = gene_id or protein_gene_map.get(protein_id)
            transcript_id = transcript_id or protein_transcript_map.get(protein_id)
            if not gene_id and protein_id in known_gene_ids:
                gene_id = protein_id
            if not transcript_id and protein_id in known_transcript_ids:
                transcript_id = protein_id
            if gene_id and gene_id in known_gene_ids:
                protein_gene_links += 1
            if transcript_id and transcript_id in known_transcript_ids:
                protein_transcript_links += 1
            conn.execute(
                """
                INSERT OR REPLACE INTO proteins
                (species_id, protein_id, gene_id, transcript_id, description, length, sequence)
                VALUES (?, ?, ?, ?, ?, ?, ?)
                """,
                (species_id, protein_id, gene_id, transcript_id, description, len(sequence), sequence),
            )
            protein_count += 1

        optional_assets = []
        if isinstance(optional, dict):
            for key, raw_path in optional.items():
                opt_path = resolve_path(manifest_path, raw_path)
                if opt_path:
                    optional_assets.append((species_id, key, str(opt_path), 1 if opt_path.exists() else 0))
        conn.executemany("INSERT INTO optional_assets VALUES (?, ?, ?, ?)", optional_assets)

        missing_seqids = sorted(annotation_seqids - seqids)
        if missing_seqids:
            preview = ", ".join(missing_seqids[:10])
            summary["warnings"].append(f"{len(missing_seqids)} annotation seqids not found in FASTA index: {preview}")
        if protein_gene_links == 0:
            summary["warnings"].append("No protein FASTA headers exposed gene= links; gene-protein mapping remains unresolved.")

        db_counts = {
            "chromosomes": conn.execute("SELECT COUNT(*) FROM chromosomes WHERE species_id = ?", (species_id,)).fetchone()[0],
            "genes": conn.execute("SELECT COUNT(*) FROM genes WHERE species_id = ?", (species_id,)).fetchone()[0],
            "transcripts": conn.execute("SELECT COUNT(*) FROM transcripts WHERE species_id = ?", (species_id,)).fetchone()[0],
            "features": conn.execute("SELECT COUNT(*) FROM features WHERE species_id = ?", (species_id,)).fetchone()[0],
            "proteins": conn.execute("SELECT COUNT(*) FROM proteins WHERE species_id = ?", (species_id,)).fetchone()[0],
        }
        summary.update(
            {
                "chromosomes": db_counts["chromosomes"],
                "genes": db_counts["genes"],
                "transcripts": db_counts["transcripts"],
                "features": db_counts["features"],
                "proteins": db_counts["proteins"],
                "annotation_gene_records": gene_count,
                "annotation_transcript_records": transcript_count,
                "annotation_feature_records": feature_count,
                "protein_records": protein_count,
                "protein_gene_links": protein_gene_links,
                "protein_transcript_links": protein_transcript_links,
                "id_mapping_rows": len(id_mapping_rows),
                "annotation_seqids": len(annotation_seqids),
                "annotation_seqids_missing_from_fasta": len(missing_seqids),
                "skipped_annotation_records_without_id": skipped_no_id,
                "optional_assets": len(optional_assets),
            }
        )
        finished = dt.datetime.now(dt.timezone.utc).isoformat()
        conn.execute(
            "UPDATE build_runs SET finished_at = ?, status = ?, summary_json = ? WHERE build_id = ?",
            (finished, "success", json.dumps(summary, ensure_ascii=False, sort_keys=True), build_id),
        )
        conn.commit()
    except Exception:
        conn.rollback()
        raise
    finally:
        conn.close()

    write_report(reports_dir / f"{species_id}_build_report.md", summary)
    return summary


def write_report(path: Path, summary: Dict[str, Any]) -> None:
    lines = [
        f"# {summary['species_id']} Build Report",
        "",
        "## Summary",
        "",
        f"- Annotation format: `{summary['annotation_format']}`",
        f"- Chromosomes/scaffolds: {summary['chromosomes']}",
        f"- Genes: {summary['genes']}",
        f"- Transcripts: {summary['transcripts']}",
        f"- Features: {summary['features']}",
        f"- Proteins: {summary['proteins']}",
        f"- Proteins linked to genes: {summary['protein_gene_links']}",
        f"- Proteins linked to transcripts: {summary['protein_transcript_links']}",
        f"- ID mapping rows loaded: {summary['id_mapping_rows']}",
        f"- Annotation seqids: {summary['annotation_seqids']}",
        f"- Annotation seqids missing from FASTA index: {summary['annotation_seqids_missing_from_fasta']}",
        f"- Optional assets recorded: {summary['optional_assets']}",
        "",
        "## Required Inputs",
        "",
    ]
    for key, value in summary["required_inputs"].items():
        lines.append(f"- `{key}`: {value}")
    lines.extend(["", "## Warnings", ""])
    if summary["warnings"]:
        lines.extend(f"- {warning}" for warning in summary["warnings"])
    else:
        lines.append("- None")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    try:
        summary = build_species(
            Path(args.manifest).resolve(),
            Path(args.db).resolve(),
            Path(args.reports_dir).resolve(),
            Path(args.assets_dir).resolve(),
            args.replace,
        )
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
