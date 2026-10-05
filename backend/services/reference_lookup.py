from __future__ import annotations

from typing import Any

from ..db import connect, row, rows
from .protein_similarity import REFERENCE_BY_GROUP, reference_gene_symbol


MAX_REFERENCE_QUERIES = 30
MAX_CANDIDATES_PER_QUERY = 50


def normalize_queries(values: list[str]) -> list[str]:
    unique: list[str] = []
    seen: set[str] = set()
    for value in values:
        query = str(value).strip()
        key = query.casefold()
        if query and key not in seen:
            seen.add(key)
            unique.append(query)
    if not unique:
        raise ValueError("at least one reference gene, transcript, or symbol is required")
    if len(unique) > MAX_REFERENCE_QUERIES:
        raise ValueError(f"a maximum of {MAX_REFERENCE_QUERIES} reference queries is allowed")
    return unique


def reference_for_species(conn, species_id: str) -> dict[str, Any]:
    species = row(conn, "SELECT taxon_group FROM species WHERE species_id = ?", (species_id,))
    if not species:
        raise ValueError(f"species not found: {species_id}")
    reference_id = REFERENCE_BY_GROUP.get(str(species["taxon_group"]))
    if not reference_id:
        raise ValueError("no reference genome is configured for this species")
    reference = row(
        conn,
        "SELECT reference_id, display_name, scientific_name, assembly_version FROM reference_species WHERE reference_id = ?",
        (reference_id,),
    )
    if not reference:
        raise ValueError("configured reference genome is unavailable")
    return reference


def resolve_reference_query(conn, reference_id: str, query: str) -> list[dict[str, Any]]:
    """Resolve a reference gene, transcript/protein ID, or annotated symbol."""
    matched_gene_ids = {
        str(item["gene_id"])
        for item in rows(
            conn,
            "SELECT gene_id FROM reference_genes WHERE reference_id = ? AND lower(gene_id) = lower(?)",
            (reference_id, query),
        )
    }
    direct_protein_matches = rows(
        conn,
        """
        SELECT protein_id, gene_id, transcript_id, description, length
        FROM reference_proteins
        WHERE reference_id = ?
          AND (lower(protein_id) = lower(?) OR lower(transcript_id) = lower(?))
        """,
        (reference_id, query, query),
    )
    alias_matches = rows(
        conn,
        """
        SELECT DISTINCT gene_id
        FROM reference_gene_aliases
        WHERE reference_id = ? AND alias = ? COLLATE NOCASE
        """,
        (reference_id, query),
    )
    matched_gene_ids.update(str(item["gene_id"]) for item in alias_matches if item.get("gene_id"))
    symbol_candidates = rows(
        conn,
        """
        SELECT protein_id, gene_id, transcript_id, description, length
        FROM reference_proteins
        WHERE reference_id = ? AND lower(description) LIKE lower(?)
        """,
        (reference_id, f"%gene_symbol:{query}%"),
    )
    for item in symbol_candidates:
        symbol = reference_gene_symbol(item.get("description"))
        if symbol and symbol.casefold() == query.casefold() and item.get("gene_id"):
            matched_gene_ids.add(str(item["gene_id"]))
    if direct_protein_matches:
        proteins = direct_protein_matches
    elif not matched_gene_ids:
        return []
    else:
        placeholders = ",".join("?" for _ in matched_gene_ids)
        proteins = rows(
            conn,
            f"""
            SELECT protein_id, gene_id, transcript_id, description, length
            FROM reference_proteins
            WHERE reference_id = ? AND gene_id IN ({placeholders})
            ORDER BY gene_id, length DESC, protein_id
            """,
            [reference_id, *sorted(matched_gene_ids)],
        )
    protein_gene_ids = sorted({str(item["gene_id"]) for item in proteins if item.get("gene_id")})
    annotation_symbols: dict[str, str] = {}
    if protein_gene_ids:
        placeholders = ",".join("?" for _ in protein_gene_ids)
        for item in rows(
            conn,
            f"""
            SELECT gene_id, alias
            FROM reference_gene_aliases
            WHERE reference_id = ?
              AND alias_type = 'annotation_symbol'
              AND gene_id IN ({placeholders})
            ORDER BY alias
            """,
            [reference_id, *protein_gene_ids],
        ):
            annotation_symbols.setdefault(str(item["gene_id"]), str(item["alias"]))
    resolved = []
    for item in proteins:
        symbol = annotation_symbols.get(str(item.get("gene_id") or "")) or reference_gene_symbol(item.get("description"))
        resolved.append(
            {
                "reference_protein_id": item["protein_id"],
                "reference_gene_id": item["gene_id"],
                "reference_transcript_id": item.get("transcript_id") or item["protein_id"],
                "reference_gene_symbol": symbol,
                "reference_protein_description": item.get("description"),
                "length": item["length"],
            }
        )
    return [
        item for item in resolved
    ]


def resolve_reference_queries(species_id: str, query_values: list[str]) -> dict[str, Any]:
    queries = normalize_queries(query_values)
    with connect() as conn:
        reference = reference_for_species(conn, species_id)
        output = [
            {"query": query, "found": bool(resolved := resolve_reference_query(conn, reference["reference_id"], query)), "matches": resolved}
            for query in queries
        ]
    return {"species_id": species_id, "reference": reference, "query_count": len(queries), "queries": output}


def lookup_similar_genes(species_id: str, query_values: list[str]) -> dict[str, Any]:
    queries = normalize_queries(query_values)
    with connect() as conn:
        reference = reference_for_species(conn, species_id)
        output = []
        for query in queries:
            resolved = resolve_reference_query(conn, reference["reference_id"], query)
            gene_ids = sorted({str(item["reference_gene_id"]) for item in resolved if item.get("reference_gene_id")})
            transcript_ids = sorted({str(item["reference_protein_id"]) for item in resolved if item.get("reference_protein_id")})
            if not gene_ids and not transcript_ids:
                output.append({"query": query, "resolved": [], "hits": []})
                continue
            clauses = ["h.species_id = ?", "h.reference_id = ?"]
            params: list[Any] = [species_id, reference["reference_id"]]
            reference_clauses = []
            if gene_ids:
                reference_clauses.append(f"h.reference_gene_id IN ({','.join('?' for _ in gene_ids)})")
                params.extend(gene_ids)
            if transcript_ids:
                reference_clauses.append(f"h.reference_protein_id IN ({','.join('?' for _ in transcript_ids)})")
                params.extend(transcript_ids)
            params.append(MAX_CANDIDATES_PER_QUERY)
            hits = rows(
                conn,
                f"""
                SELECT h.gene_id, p.transcript_id, h.protein_id,
                       h.reference_gene_id, h.reference_protein_id,
                       h.pident, h.alignment_length, h.evalue, h.bitscore, h.rank
                FROM similarity_hits h
                LEFT JOIN proteins p
                  ON p.species_id = h.species_id AND p.protein_id = h.protein_id
                WHERE {' AND '.join(clauses)} AND ({' OR '.join(reference_clauses)})
                ORDER BY h.bitscore DESC, h.evalue ASC, h.pident DESC, h.gene_id, h.rank
                LIMIT ?
                """,
                params,
            )
            output.append({"query": query, "resolved": resolved, "hits": hits})
    return {
        "species_id": species_id,
        "reference": reference,
        "query_count": len(queries),
        "queries": output,
    }
