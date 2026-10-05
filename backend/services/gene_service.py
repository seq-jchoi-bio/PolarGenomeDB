from __future__ import annotations

import sqlite3
from typing import Any, Dict, List, Optional

from ..db import parse_json, row, rows


def normalize_feature_json(items: List[Dict[str, Any]]) -> None:
    for item in items:
        item["attributes"] = parse_json(item.pop("attributes_json", None))


def gene_summary(conn: sqlite3.Connection, species_id: str, gene_id: str) -> Optional[Dict[str, Any]]:
    gene = row(
        conn,
        """
        SELECT gene_id, seqid, start, end, strand, source, attributes_json
        FROM genes
        WHERE species_id = ? AND gene_id = ?
        """,
        (species_id, gene_id),
    )
    if gene:
        gene["attributes"] = parse_json(gene.pop("attributes_json", None))
    return gene


def protein_summary(conn: sqlite3.Connection, species_id: str, protein_id: str) -> Optional[Dict[str, Any]]:
    return row(
        conn,
        """
        SELECT protein_id, gene_id, transcript_id, description, length
        FROM proteins
        WHERE species_id = ? AND protein_id = ?
        """,
        (species_id, protein_id),
    )


def overlapping_gene_summary(
    conn: sqlite3.Connection,
    species_id: str,
    seqid: str,
    start: int,
    end: int,
) -> Optional[Dict[str, Any]]:
    return row(
        conn,
        """
        SELECT gene_id, seqid, start, end, strand
        FROM genes
        WHERE species_id = ? AND seqid = ? AND start <= ? AND end >= ?
        ORDER BY (MIN(end, ?) - MAX(start, ?)) DESC, gene_id
        LIMIT 1
        """,
        (species_id, seqid, end, start, end, start),
    )


def transcript_feature_tree(conn: sqlite3.Connection, species_id: str, gene_id: str) -> Dict[str, Any]:
    transcripts = rows(
        conn,
        """
        SELECT transcript_id, gene_id, seqid, start, end, strand, source, attributes_json
        FROM transcripts
        WHERE species_id = ? AND gene_id = ?
        ORDER BY start, end, transcript_id
        """,
        (species_id, gene_id),
    )
    features = rows(
        conn,
        """
        SELECT feature_id, feature_type, parent_id, gene_id, transcript_id, seqid, start, end, strand, phase, attributes_json
        FROM features
        WHERE species_id = ? AND gene_id = ?
        ORDER BY start, end, feature_type, feature_id
        """,
        (species_id, gene_id),
    )
    proteins = rows(
        conn,
        """
        SELECT protein_id, gene_id, transcript_id, description, length
        FROM proteins
        WHERE species_id = ? AND gene_id = ?
        ORDER BY protein_id
        """,
        (species_id, gene_id),
    )
    normalize_feature_json(transcripts)
    normalize_feature_json(features)

    proteins_by_transcript: Dict[str, List[Dict[str, Any]]] = {}
    for protein in proteins:
        proteins_by_transcript.setdefault(str(protein.get("transcript_id") or ""), []).append(protein)

    features_by_transcript: Dict[str, List[Dict[str, Any]]] = {}
    gene_level_features: List[Dict[str, Any]] = []
    for feature in features:
        transcript_id = feature.get("transcript_id")
        parent_id = feature.get("parent_id")
        key = str(transcript_id or parent_id or "")
        if key:
            features_by_transcript.setdefault(key, []).append(feature)
        else:
            gene_level_features.append(feature)

    transcript_models = []
    for transcript in transcripts:
        transcript_id = str(transcript["transcript_id"])
        transcript_models.append(
            {
                **transcript,
                "features": features_by_transcript.get(transcript_id, []),
                "proteins": proteins_by_transcript.get(transcript_id, []),
            }
        )

    orphan_features = [
        feature
        for key, grouped in features_by_transcript.items()
        if key not in {str(transcript["transcript_id"]) for transcript in transcripts}
        for feature in grouped
    ]

    return {
        "transcripts": transcripts,
        "features": features,
        "proteins": proteins,
        "transcript_models": transcript_models,
        "gene_level_features": gene_level_features,
        "orphan_features": orphan_features,
    }


def build_gene_detail(conn: sqlite3.Connection, species_id: str, gene_id: str) -> Optional[Dict[str, Any]]:
    gene = gene_summary(conn, species_id, gene_id)
    if not gene:
        return None
    tree = transcript_feature_tree(conn, species_id, gene_id)
    return {
        "species_id": species_id,
        "gene": gene,
        **tree,
    }


def annotate_search_hit(conn: sqlite3.Connection, species_id: str, program: str, hit: Dict[str, Any]) -> Dict[str, Any]:
    target_id = str(hit.get("sseqid") or "")
    annotation: Dict[str, Any] = {"target_type": "unknown", "target_id": target_id}

    if program in {"blastp", "diamond"}:
        protein = protein_summary(conn, species_id, target_id)
        if protein:
            annotation.update(
                {
                    "target_type": "protein",
                    "protein_id": protein["protein_id"],
                    "gene_id": protein["gene_id"],
                    "transcript_id": protein["transcript_id"],
                    "detail_url": f"/species/{species_id}/genes/{protein['gene_id']}" if protein["gene_id"] else None,
                }
            )
            return {**hit, "annotation": annotation}

    gene = gene_summary(conn, species_id, target_id)
    if gene:
        annotation.update(
            {
                "target_type": "gene",
                "gene_id": gene["gene_id"],
                "seqid": gene["seqid"],
                "start": gene["start"],
                "end": gene["end"],
                "strand": gene["strand"],
                "detail_url": f"/species/{species_id}/genes/{gene['gene_id']}",
            }
        )
    else:
        hit_start = min(int(hit["sstart"]), int(hit["send"])) if program == "blastn" else None
        hit_end = max(int(hit["sstart"]), int(hit["send"])) if program == "blastn" else None
        region_gene = (
            overlapping_gene_summary(conn, species_id, target_id, hit_start, hit_end)
            if hit_start is not None and hit_end is not None
            else None
        )
        annotation.update(
            {
                "target_type": "gene" if region_gene else "genomic_region" if program == "blastn" else "unknown",
                "gene_id": region_gene["gene_id"] if region_gene else None,
                "seqid": target_id if program == "blastn" else None,
                "start": hit_start,
                "end": hit_end,
                "strand": region_gene["strand"] if region_gene else None,
                "detail_url": f"/species/{species_id}/genes/{region_gene['gene_id']}" if region_gene else None,
                "viewer_url": f"/viewer/{species_id}/{target_id}:{hit_start}-{hit_end}" if hit_start is not None and hit_end is not None else None,
            }
        )
    return {**hit, "annotation": annotation}
