from __future__ import annotations

from pathlib import Path
from io import StringIO
from typing import Optional
from urllib.parse import quote

from fastapi import APIRouter, HTTPException, Query
from fastapi.responses import PlainTextResponse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

from ..api_limits import get_api_limits
from ..db import connect, parse_json, row, rows
from ..services.gene_service import build_gene_detail


router = APIRouter()
GENE_LIMITS = get_api_limits().genes


def escape_like(value: str) -> str:
    return value.replace("\\", "\\\\").replace("%", "\\%").replace("_", "\\_")


def reverse_complement(sequence: str) -> str:
    table = str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn", "TGCAYRMKSWVHDBNtgcayrmkswvhdbn")
    return sequence.translate(table)[::-1].upper()


def genome_paths(species_id: str) -> tuple[Path, Path]:
    with connect() as conn:
        data = rows(
            conn,
            "SELECT asset_type, path FROM genome_assets WHERE species_id = ? AND asset_type IN ('fasta', 'fai')",
            (species_id,),
        )
    values = {item["asset_type"]: Path(item["path"]) for item in data}
    if "fasta" not in values:
        raise HTTPException(status_code=404, detail=f"genome FASTA asset not found: {species_id}")
    fasta = values["fasta"]
    fai = values.get("fai", Path(str(fasta) + ".fai"))
    if not fasta.exists():
        raise HTTPException(status_code=404, detail=f"genome FASTA does not exist: {fasta}")
    if not fai.exists():
        raise HTTPException(status_code=404, detail=f"genome FASTA index does not exist: {fai}")
    return fasta, fai


def fetch_fasta_interval(species_id: str, seqid: str, start: int, end: int) -> tuple[str, int]:
    fasta, fai = genome_paths(species_id)
    fai_row = None
    with fai.open("rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[0] == seqid:
                fai_row = fields
                break
    if fai_row is None:
        raise HTTPException(status_code=404, detail=f"seqid not found in FASTA index: {seqid}")

    chrom_length, offset, line_bases, line_width = map(int, fai_row[1:5])
    bounded_start = max(1, start)
    bounded_end = min(chrom_length, end)
    if bounded_start > bounded_end:
        raise HTTPException(status_code=400, detail="requested interval is outside chromosome bounds")

    parts: list[str] = []
    with fasta.open("rb") as handle:
        pos = bounded_start
        while pos <= bounded_end:
            line_index = (pos - 1) // line_bases
            in_line = (pos - 1) % line_bases
            take = min(bounded_end - pos + 1, line_bases - in_line)
            byte_offset = offset + line_index * line_width + in_line
            handle.seek(byte_offset)
            parts.append(handle.read(take).decode("ascii"))
            pos += take
    return "".join(parts).upper(), chrom_length


def feature_to_display_interval(feature: dict, interval_start: int, interval_end: int, strand: str) -> dict:
    if strand == "-":
        start = interval_end - int(feature["end"]) + 1
        end = interval_end - int(feature["start"]) + 1
    else:
        start = int(feature["start"]) - interval_start + 1
        end = int(feature["end"]) - interval_start + 1
    return {
        **feature,
        "display_start": max(1, min(start, end)),
        "display_end": max(1, max(start, end)),
    }


def transcript_sequence_view(species_id: str, gene: dict, transcript: dict, features: list[dict]) -> dict:
    transcript_id = transcript["transcript_id"]
    exon_features = [
        item for item in features
        if item.get("transcript_id") == transcript_id and item.get("feature_type") == "exon"
    ]
    cds_features = [
        item for item in features
        if item.get("transcript_id") == transcript_id and item.get("feature_type") == "CDS"
    ]
    if not exon_features:
        return {
            "mode": "transcript",
            "transcript_id": transcript_id,
            "sequence": "",
            "features": [],
            "note": "No exon features are available for transcript sequence extraction.",
        }

    genomic_order = sorted(exon_features, key=lambda item: (int(item["start"]), int(item["end"])))
    parts = [
        fetch_fasta_interval(species_id, gene["seqid"], int(exon["start"]), int(exon["end"]))[0]
        for exon in genomic_order
    ]
    if gene.get("strand") == "-":
        ordered_exons = list(reversed(genomic_order))
        sequence = "".join(reverse_complement(part) for part in reversed(parts))
    else:
        ordered_exons = genomic_order
        sequence = "".join(parts)

    display_features: list[dict] = []
    cursor = 1
    for exon in ordered_exons:
        exon_length = int(exon["end"]) - int(exon["start"]) + 1
        exon_display = {
            **exon,
            "display_start": cursor,
            "display_end": cursor + exon_length - 1,
        }
        display_features.append(exon_display)
        for cds in cds_features:
            overlap_start = max(int(exon["start"]), int(cds["start"]))
            overlap_end = min(int(exon["end"]), int(cds["end"]))
            if overlap_start > overlap_end:
                continue
            if gene.get("strand") == "-":
                display_start = cursor + (int(exon["end"]) - overlap_end)
                display_end = cursor + (int(exon["end"]) - overlap_start)
            else:
                display_start = cursor + (overlap_start - int(exon["start"]))
                display_end = cursor + (overlap_end - int(exon["start"]))
            display_features.append({**cds, "display_start": display_start, "display_end": display_end})
        cursor += exon_length

    return {
        "mode": "transcript",
        "transcript_id": transcript_id,
        "seqid": gene["seqid"],
        "start": transcript["start"],
        "end": transcript["end"],
        "strand": gene.get("strand"),
        "sequence": sequence,
        "length": len(sequence),
        "feature_coordinate": "sequence",
        "features": display_features,
    }


@router.get("/species/{species_id}/genes")
def search_genes(
    species_id: str,
    q: Optional[str] = Query(default=None, min_length=1),
    seqid: Optional[str] = None,
    start: Optional[int] = Query(default=None, ge=1),
    end: Optional[int] = Query(default=None, ge=1),
    limit: int = Query(default=GENE_LIMITS.default_page_size, ge=1, le=GENE_LIMITS.max_page_size),
    offset: int = Query(default=0, ge=0),
) -> dict:
    clauses = ["species_id = ?"]
    params: list[object] = [species_id]
    if q:
        clauses.append("gene_id LIKE ?")
        params.append(f"%{q}%")
    if seqid:
        clauses.append("seqid = ?")
        params.append(seqid)
    if start is not None:
        clauses.append("end >= ?")
        params.append(start)
    if end is not None:
        clauses.append("start <= ?")
        params.append(end)
    params.extend([limit, offset])
    sql = f"""
        SELECT gene_id, seqid, start, end, strand, source, attributes_json
        FROM genes
        WHERE {' AND '.join(clauses)}
        ORDER BY seqid, start, end, gene_id
        LIMIT ? OFFSET ?
    """
    with connect() as conn:
        data = rows(conn, sql, params)
    for item in data:
        item["attributes"] = parse_json(item.pop("attributes_json", None))
    return {"species_id": species_id, "limit": limit, "offset": offset, "genes": data}


@router.get("/species/{species_id}/genes/suggest")
def suggest_genes(
    species_id: str,
    q: str = Query(min_length=1),
    limit: int = Query(default=GENE_LIMITS.default_suggest_size, ge=1, le=GENE_LIMITS.max_suggest_size),
) -> dict:
    value = escape_like(q.strip())
    prefix = f"{value}%"
    contains = f"%{value}%"
    sql = """
        SELECT gene_id, seqid, start, end, strand, source
        FROM genes
        WHERE species_id = ?
          AND (gene_id LIKE ? ESCAPE '\\' OR gene_id LIKE ? ESCAPE '\\')
        ORDER BY
          CASE WHEN gene_id LIKE ? ESCAPE '\\' THEN 0 ELSE 1 END,
          seqid,
          start,
          end,
          gene_id
        LIMIT ?
    """
    with connect() as conn:
        species_exists = conn.execute("SELECT 1 FROM species WHERE species_id = ?", (species_id,)).fetchone()
        if not species_exists:
            raise HTTPException(status_code=404, detail=f"species not found: {species_id}")
        data = rows(conn, sql, (species_id, prefix, contains, prefix, limit))

    suggestions = [
        {
            **item,
            "label": item["gene_id"],
            "location": f"{item['seqid']}:{item['start']}-{item['end']}",
            "match_type": "prefix" if item["gene_id"].startswith(q.strip()) else "contains",
        }
        for item in data
    ]
    return {"species_id": species_id, "q": q, "limit": limit, "suggestions": suggestions}


@router.get("/species/{species_id}/genes/{gene_id}")
def get_gene(species_id: str, gene_id: str) -> dict:
    with connect() as conn:
        detail = build_gene_detail(conn, species_id, gene_id)
        if not detail:
            raise HTTPException(status_code=404, detail=f"gene not found: {species_id}/{gene_id}")
    return detail


@router.get("/species/{species_id}/genes/{gene_id}/export")
def export_gene(species_id: str, gene_id: str, format: str = Query(pattern="^(fasta|gtf)$")) -> PlainTextResponse:
    detail = get_gene(species_id, gene_id)
    gene = detail["gene"]
    if format == "fasta":
        sequence, _ = fetch_fasta_interval(species_id, gene["seqid"], int(gene["start"]), int(gene["end"]))
        if gene.get("strand") == "-":
            sequence = reverse_complement(sequence)
        header = f">{gene_id} {gene['seqid']}:{gene['start']}-{gene['end']} strand={gene.get('strand', '.')} orientation=gene"
        content = header + "\n" + "\n".join(sequence[i:i + 80] for i in range(0, len(sequence), 80)) + "\n"
    else:
        def attribute(value: str) -> str:
            return str(value).replace("\\", "\\\\").replace('"', '\\"').replace("\n", " ").replace("\r", " ").replace("\t", " ")

        lines = []
        def add(item: dict, kind: str, transcript_id: str | None = None) -> None:
            attrs = f'gene_id "{attribute(gene_id)}";'
            if transcript_id:
                attrs += f' transcript_id "{attribute(transcript_id)}";'
            phase = item.get("phase") if kind == "CDS" else "."
            lines.append("\t".join(str(value) for value in [item["seqid"], item.get("source") or gene.get("source") or "PolarGenomeDB", kind, item["start"], item["end"], ".", item.get("strand") or ".", phase if phase is not None else ".", attrs]))

        add(gene, "gene")
        for transcript in detail["transcripts"]:
            add(transcript, "transcript", transcript["transcript_id"])
        for feature in detail["features"]:
            if feature["feature_type"] not in {"gene", "transcript", "mRNA"}:
                add(feature, feature["feature_type"], feature.get("transcript_id"))
        content = "\n".join(lines) + "\n"
    filename = quote(f"{species_id}_{gene_id}.{format}", safe="")
    return PlainTextResponse(content, headers={"Content-Disposition": f"attachment; filename*=UTF-8''{filename}"})


@router.get("/species/{species_id}/genes/{gene_id}/sequence.gb")
def export_sequence_genbank(
    species_id: str,
    gene_id: str,
    mode: str = Query(default="gene", pattern="^(gene|reference|transcript|protein)$"),
    upstream: int = Query(default=0, ge=0, le=10000),
    downstream: int = Query(default=0, ge=0, le=10000),
    transcript_id: Optional[str] = None,
) -> PlainTextResponse:
    payload = get_gene_sequence(species_id, gene_id, upstream, downstream, transcript_id)
    view = payload["views"][mode]
    if not view.get("sequence"):
        raise HTTPException(status_code=404, detail="No sequence available for this selection")
    identifier = view.get("protein_id") or view.get("transcript_id") or gene_id
    record = SeqRecord(Seq(view["sequence"]), id=identifier, name=identifier, description=f"{species_id} {gene_id}; {mode} sequence")
    record.annotations["molecule_type"] = "protein" if mode == "protein" else "DNA"
    record.annotations["topology"] = "linear"
    record.annotations["comment"] = f"Polar Genome DB export; mode={mode}; gene={gene_id}; transcript={view.get('transcript_id') or 'all'}."
    if view.get("seqid"):
        record.annotations["comment"] += f" Assembly locus {view['seqid']}:{view['start']}-{view['end']}; strand={view.get('strand', '.')} (1-based inclusive)."
    record.features.append(SeqFeature(FeatureLocation(0, len(record)), type="source", qualifiers={"note": [f"Species ID: {species_id}"]}))
    cds_groups: dict[str, list] = {}
    for feature in view.get("features", []):
        kind = feature.get("feature_type")
        if kind not in {"exon", "CDS"}:
            continue
        if mode == "reference":
            start = int(feature["start"]) - int(view["start"])
            end = int(feature["end"]) - int(view["start"]) + 1
            strand = -1 if feature.get("strand") == "-" else 1
        else:
            start = int(feature["display_start"]) - 1
            end = int(feature["display_end"])
            strand = 1
        if start < 0 or end > len(record) or start >= end:
            continue
        location = FeatureLocation(start, end, strand=strand)
        tid = feature.get("transcript_id") or ""
        qualifiers = {"gene": [gene_id]}
        if tid:
            qualifiers["transcript_id"] = [tid]
        if kind == "CDS":
            cds_groups.setdefault(tid, []).append((location, feature.get("phase"), qualifiers))
        else:
            record.features.append(SeqFeature(location, type="exon", qualifiers=qualifiers))
    for parts in cds_groups.values():
        parts.sort(key=lambda item: int(item[0].start), reverse=parts[0][0].strand == -1)
        locations = list(dict.fromkeys((int(p[0].start), int(p[0].end), p[0].strand) for p in parts))
        locations = [FeatureLocation(start, end, strand=strand) for start, end, strand in locations]
        qualifiers = dict(parts[0][2])
        if str(parts[0][1]) in {"0", "1", "2"}:
            qualifiers["codon_start"] = [str(int(parts[0][1]) + 1)]
        record.features.append(SeqFeature(locations[0] if len(locations) == 1 else CompoundLocation(locations), type="CDS", qualifiers=qualifiers))
    output = StringIO()
    SeqIO.write(record, output, "genbank")
    filename = quote(f"{species_id}_{identifier}_{mode}.gb", safe="")
    return PlainTextResponse(output.getvalue(), headers={"Content-Disposition": f"attachment; filename*=UTF-8''{filename}"})


@router.get("/species/{species_id}/genes/{gene_id}/sequence")
def get_gene_sequence(
    species_id: str,
    gene_id: str,
    upstream: int = Query(default=0, ge=0, le=10000),
    downstream: int = Query(default=0, ge=0, le=10000),
    transcript_id: Optional[str] = Query(default=None),
) -> dict:
    with connect() as conn:
        detail = build_gene_detail(conn, species_id, gene_id)
        if not detail:
            raise HTTPException(status_code=404, detail=f"gene not found: {species_id}/{gene_id}")
        proteins = rows(
            conn,
            """
            SELECT protein_id, transcript_id, length, sequence
            FROM proteins
            WHERE species_id = ? AND gene_id = ?
            ORDER BY transcript_id, protein_id
            """,
            (species_id, gene_id),
        )

    gene = detail["gene"]
    start = max(1, int(gene["start"]) - upstream)
    end = int(gene["end"]) + downstream
    reference_sequence, chrom_length = fetch_fasta_interval(species_id, gene["seqid"], start, end)
    end = min(end, chrom_length)
    gene_oriented = reverse_complement(reference_sequence) if gene.get("strand") == "-" else reference_sequence
    genomic_features = [
        feature_to_display_interval(item, start, end, str(gene.get("strand") or "+"))
        for item in detail["features"]
        if item.get("feature_type") in {"exon", "CDS"}
    ]

    transcripts = detail.get("transcripts") or []
    selected_transcript = None
    if transcript_id:
        selected_transcript = next((item for item in transcripts if item["transcript_id"] == transcript_id), None)
    if selected_transcript is None and transcripts:
        selected_transcript = transcripts[0]

    selected_transcript_id = selected_transcript["transcript_id"] if selected_transcript else None
    protein = next(
        (item for item in proteins if item.get("transcript_id") == selected_transcript_id),
        None,
    )

    transcript_view = (
        transcript_sequence_view(species_id, gene, selected_transcript, detail["features"])
        if selected_transcript
        else {"mode": "transcript", "sequence": "", "features": []}
    )
    protein_view = {
        "mode": "protein",
        "protein_id": protein["protein_id"] if protein else None,
        "transcript_id": selected_transcript_id,
        "sequence": protein["sequence"] if protein else "",
        "length": protein["length"] if protein else 0,
        "note": (
            None
            if protein
            else f"No translated protein sequence is available for {selected_transcript_id} in the current protein FASTA."
        ),
    }

    return {
        "species_id": species_id,
        "gene_id": gene_id,
        "upstream": upstream,
        "downstream": downstream,
        "views": {
            "gene": {
                "mode": "gene",
                "seqid": gene["seqid"],
                "start": start,
                "end": end,
                "strand": gene.get("strand"),
                "sequence": gene_oriented,
                "length": len(gene_oriented),
                "feature_coordinate": "sequence",
                "features": genomic_features,
            },
            "reference": {
                "mode": "reference",
                "seqid": gene["seqid"],
                "start": start,
                "end": end,
                "strand": "+",
                "sequence": reference_sequence,
                "length": len(reference_sequence),
                "feature_coordinate": "genome",
                "features": detail["features"],
            },
            "transcript": transcript_view,
            "protein": protein_view,
        },
    }
