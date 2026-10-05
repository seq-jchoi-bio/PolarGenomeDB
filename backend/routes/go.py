from __future__ import annotations

from functools import lru_cache
import json
import math
import os
import re
import subprocess
import tempfile
from typing import Literal

from fastapi import APIRouter, HTTPException, Query
from pydantic import BaseModel, Field

from ..api_limits import get_api_limits
from ..config import GOFIGURE_PYTHON, GOFIGURE_SCRIPT, RUNTIME_DIR
from ..db import connect, parse_json, row, rows


router = APIRouter()
GO_LIMITS = get_api_limits().go

WANG_EDGE_WEIGHTS = {
    "is_a": 0.8,
    "part_of": 0.6,
    "regulates": 0.6,
    "positively_regulates": 0.6,
    "negatively_regulates": 0.6,
}

GO_NAMESPACES = {
    "biological_process": "BP",
    "molecular_function": "MF",
    "cellular_component": "CC",
}
GO_NAMESPACE_LABELS = {
    "biological_process": "Biological Process",
    "molecular_function": "Molecular Function",
    "cellular_component": "Cellular Component",
}
GO_CONFIDENCE_LEVELS = ("high", "medium", "low")
GO_INPUT_SPLIT = re.compile(r"[\s,;]+")


class GoEnrichmentRequest(BaseModel):
    genes: str = Field(min_length=1, max_length=250_000)
    background: str = Field(default="", max_length=2_500_000)
    min_count: int = Field(default=3, ge=1, le=1000)
    correction: Literal["benjamini-hochberg", "bonferroni", "holm"] = "benjamini-hochberg"
    confidence: list[Literal["high", "medium", "low"]] = Field(default_factory=lambda: ["high", "medium"])
    max_similarity_rank: int = Field(default=5, ge=1, le=20)


class GoSemanticTerm(BaseModel):
    go_id: str = Field(pattern=r"^GO:\d{7}$")
    pvalue: float = Field(ge=0)


class GoSemanticRequest(BaseModel):
    section: Literal["BP", "MF", "CC"]
    terms: list[GoSemanticTerm] = Field(min_length=2, max_length=30)
    palette: Literal["plasma", "viridis", "cividis", "magma", "Set2"] = "plasma"
    si: float = Field(default=0.5, ge=0.1, le=0.9)
    max_terms: int = Field(default=15, ge=5, le=30)
    label_mode: Literal["numbered", "description", "description-numbered"] = "numbered"
    bubble_size: Literal["members", "frequency"] = "members"
    font_size: Literal["small", "medium", "large"] = "small"
    legend_size: Literal["small", "medium", "large"] = "small"


def term_display(item: dict) -> dict:
    return {
        **item,
        "go_name": item.get("ontology_name") or item.get("go_name"),
        "namespace": item.get("ontology_namespace") or item.get("namespace"),
    }


@lru_cache(maxsize=1)
def load_ontology_graph() -> tuple[dict[str, dict], dict[str, list[tuple[str, str]]]]:
    with connect() as conn:
        terms = {
            item["go_id"]: item
            for item in rows(
                conn,
                """
                SELECT go_id, name, namespace, definition, is_obsolete, metadata_json
                FROM go_terms
                """,
            )
        }
        parents: dict[str, list[tuple[str, str]]] = {}
        for edge in rows(conn, "SELECT child_go_id, parent_go_id, relation_type FROM go_edges"):
            parents.setdefault(edge["child_go_id"], []).append((edge["parent_go_id"], edge["relation_type"]))
    return terms, parents


@lru_cache(maxsize=20000)
def wang_ancestor_weights(go_id: str) -> tuple[tuple[str, float], ...]:
    _terms, parents = load_ontology_graph()
    scores: dict[str, float] = {go_id: 1.0}
    queue = [go_id]
    while queue:
        child = queue.pop(0)
        child_score = scores[child]
        for parent, relation in parents.get(child, []):
            weight = WANG_EDGE_WEIGHTS.get(relation, 0.6)
            score = child_score * weight
            if score > scores.get(parent, 0.0):
                scores[parent] = score
                queue.append(parent)
    return tuple(sorted(scores.items()))


def wang_similarity(a: str, b: str) -> float:
    if a == b:
        return 1.0
    weights_a = dict(wang_ancestor_weights(a))
    weights_b = dict(wang_ancestor_weights(b))
    if not weights_a or not weights_b:
        return 0.0
    shared = set(weights_a) & set(weights_b)
    if not shared:
        return 0.0
    numerator = sum(weights_a[item] + weights_b[item] for item in shared)
    denominator = sum(weights_a.values()) + sum(weights_b.values())
    return numerator / denominator if denominator else 0.0


def ancestors_for(go_id: str) -> list[dict]:
    terms, _parents = load_ontology_graph()
    result = []
    for ancestor_id, weight in wang_ancestor_weights(go_id):
        if ancestor_id == go_id:
            continue
        term = terms.get(ancestor_id, {})
        result.append(
            {
                "go_id": ancestor_id,
                "name": term.get("name"),
                "namespace": term.get("namespace"),
                "semantic_weight": weight,
            }
        )
    result.sort(key=lambda item: (-item["semantic_weight"], item["go_id"]))
    return result


def split_gene_ids(value: str) -> list[str]:
    seen: set[str] = set()
    result: list[str] = []
    for token in GO_INPUT_SPLIT.split(value.strip()):
        if not token:
            continue
        key = token.casefold()
        if key in seen:
            continue
        seen.add(key)
        result.append(token)
    return result


def chunks(items: list[str], size: int = 900):
    for index in range(0, len(items), size):
        yield items[index:index + size]


def resolve_species_gene_ids(conn, species_id: str, submitted: list[str]) -> tuple[list[str], list[str]]:
    if not submitted:
        return [], []
    resolved: dict[str, str] = {}
    for part in chunks(submitted):
        placeholders = ",".join("?" for _ in part)
        for item in rows(
            conn,
            f"""
            SELECT gene_id
            FROM genes
            WHERE species_id = ? AND lower(gene_id) IN ({placeholders})
            """,
            [species_id, *[gene.casefold() for gene in part]],
        ):
            resolved[item["gene_id"].casefold()] = item["gene_id"]
    valid: list[str] = []
    invalid: list[str] = []
    for submitted_id in submitted:
        gene_id = resolved.get(submitted_id.casefold())
        if gene_id:
            valid.append(gene_id)
        else:
            invalid.append(submitted_id)
    return valid, invalid


def similarity_confidence(pident: float | None) -> str:
    value = float(pident or 0)
    if value >= 60:
        return "high"
    if value >= 40:
        return "medium"
    return "low"


def go_mapping_conditions(confidence: list[str], max_similarity_rank: int) -> tuple[str, list[object]]:
    selected = sorted(set(confidence) & set(GO_CONFIDENCE_LEVELS)) or ["high", "medium"]
    tier_clauses = []
    if "high" in selected:
        tier_clauses.append("h.pident >= 60")
    if "medium" in selected:
        tier_clauses.append("h.pident >= 40 AND h.pident < 60")
    if "low" in selected:
        tier_clauses.append("h.pident < 40")
    return f"h.rank <= ? AND ({' OR '.join(tier_clauses)})", [max_similarity_rank]


def go_mapping_sql(mode: Literal["direct", "all"], extra_clause: str = "") -> str:
    if mode == "all":
        term_id = "ga.ancestor_go_id"
        level = "ga.level"
        ontology_join = """
        JOIN go_ancestors ga ON ga.go_id = go.go_id
        LEFT JOIN go_terms t ON t.go_id = ga.ancestor_go_id
        """
    else:
        term_id = "go.go_id"
        level = "ga.level"
        ontology_join = """
        LEFT JOIN go_ancestors ga ON ga.go_id = go.go_id AND ga.ancestor_go_id = go.go_id
        LEFT JOIN go_terms t ON t.go_id = go.go_id
        """
    return f"""
        SELECT DISTINCT
          h.gene_id,
          {term_id} AS go_id,
          COALESCE(t.name, go.go_name, {term_id}) AS term,
          COALESCE(t.namespace, go.namespace) AS namespace,
          {level} AS level
        FROM similarity_hits h
        JOIN species s ON s.species_id = h.species_id
        JOIN reference_species rs ON rs.reference_id = h.reference_id
        JOIN reference_go_annotations go
          ON go.reference_id = h.reference_id
         AND go.gene_id = h.reference_gene_id
        {ontology_join}
        WHERE h.species_id = ?
          AND s.taxon_group = rs.taxon_group
          AND COALESCE(t.is_obsolete, 0) = 0
          AND COALESCE(t.namespace, go.namespace) IN ('biological_process', 'molecular_function', 'cellular_component')
          AND {{conditions}}
          {extra_clause}
    """


def population_statistics(
    conn,
    species_id: str,
    confidence: list[str],
    max_similarity_rank: int,
    background_ids: list[str] | None,
    mode: Literal["direct", "all"],
) -> tuple[dict[str, int], dict[tuple[str, str], int]]:
    conditions, condition_params = go_mapping_conditions(confidence, max_similarity_rank)
    background_join = ""
    if background_ids is not None:
        conn.execute("DROP TABLE IF EXISTS temp.go_enrichment_background")
        conn.execute("CREATE TEMP TABLE go_enrichment_background (gene_id TEXT PRIMARY KEY)")
        conn.executemany("INSERT OR IGNORE INTO go_enrichment_background (gene_id) VALUES (?)", [(gene_id,) for gene_id in background_ids])
        background_join = "JOIN go_enrichment_background bg ON bg.gene_id = h.gene_id"

    base = go_mapping_sql(mode).replace("{conditions}", conditions).replace("FROM similarity_hits h", f"FROM similarity_hits h {background_join}")
    population_rows = rows(
        conn,
        f"""
        SELECT namespace, COUNT(DISTINCT gene_id) AS gene_count
        FROM ({base}) mapped
        GROUP BY namespace
        """,
        [species_id, *condition_params],
    )
    term_rows = rows(
        conn,
        f"""
        SELECT
          namespace,
          go_id,
          COUNT(DISTINCT gene_id) AS gene_count
        FROM ({base}) mapped
        GROUP BY namespace, go_id
        """,
        [species_id, *condition_params],
    )
    if background_ids is not None:
        conn.execute("DROP TABLE IF EXISTS temp.go_enrichment_background")
    return (
        {item["namespace"]: int(item["gene_count"]) for item in population_rows},
        {(item["namespace"], item["go_id"]): int(item["gene_count"]) for item in term_rows},
    )


@lru_cache(maxsize=16)
def default_population_statistics(
    species_id: str,
    confidence: tuple[str, ...],
    max_similarity_rank: int,
    mode: Literal["direct", "all"],
) -> tuple[dict[str, int], dict[tuple[str, str], int]]:
    """Cache the static annotated population used when no custom background is supplied."""
    path = RUNTIME_DIR / species_id / "go_population.json"
    if confidence == ("high", "medium") and max_similarity_rank == 5 and path.exists():
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
            mode_payload = payload.get("modes", {}).get(mode, {})
            population_counts = {str(key): int(value) for key, value in mode_payload.get("population_counts", {}).items()}
            term_counts = {
                (str(item["namespace"]), str(item["go_id"])): int(item["gene_count"])
                for item in mode_payload.get("term_counts", [])
            }
            if population_counts and term_counts:
                return population_counts, term_counts
        except (OSError, ValueError, TypeError, KeyError):
            pass
    with connect() as conn:
        return population_statistics(conn, species_id, list(confidence), max_similarity_rank, None, mode)


def gene_go_terms(
    conn,
    species_id: str,
    gene_ids: list[str],
    confidence: list[str],
    max_similarity_rank: int,
    mode: Literal["direct", "all"],
) -> dict[str, dict[str, dict[str, dict]]]:
    if not gene_ids:
        return {}
    conditions, condition_params = go_mapping_conditions(confidence, max_similarity_rank)
    mapping: dict[str, dict[str, dict[str, dict]]] = {}
    for part in chunks(gene_ids):
        placeholders = ",".join("?" for _ in part)
        base = go_mapping_sql(mode, f"AND h.gene_id IN ({placeholders})").replace("{conditions}", conditions)
        for item in rows(
            conn,
            f"""
            SELECT gene_id, go_id, term, namespace, level
            FROM ({base}) mapped
            """,
            [species_id, *condition_params, *part],
        ):
            namespace = item["namespace"]
            mapping.setdefault(item["gene_id"], {}).setdefault(namespace, {})[item["go_id"]] = {
                "go_id": item["go_id"],
                "term": item["term"],
                "namespace": namespace,
                "level": item.get("level"),
            }
    return mapping


def hypergeom_tail(k: int, n: int, K: int, N: int) -> float:
    if k <= 0 or k > min(K, n) or n > N or K > N:
        return 1.0
    denominator = math.lgamma(N + 1) - math.lgamma(n + 1) - math.lgamma(N - n + 1)
    values = []
    for index in range(k, min(K, n) + 1):
        if n - index > N - K:
            continue
        values.append(
            math.lgamma(K + 1) - math.lgamma(index + 1) - math.lgamma(K - index + 1)
            + math.lgamma(N - K + 1) - math.lgamma(n - index + 1) - math.lgamma(N - K - (n - index) + 1)
            - denominator
        )
    if not values:
        return 1.0
    peak = max(values)
    return min(1.0, math.exp(peak) * sum(math.exp(value - peak) for value in values))


def add_multiple_test_corrections(items: list[dict], correction: str) -> None:
    count = len(items)
    if not count:
        return
    ordered = sorted(enumerate(items), key=lambda pair: pair[1]["pvalue"])
    bh = [1.0] * count
    previous = 1.0
    for rank in range(count, 0, -1):
        index, item = ordered[rank - 1]
        value = min(previous, item["pvalue"] * count / rank)
        bh[index] = min(1.0, value)
        previous = value
    holm = [1.0] * count
    running = 0.0
    for rank, (index, item) in enumerate(ordered, start=1):
        running = max(running, min(1.0, item["pvalue"] * (count - rank + 1)))
        holm[index] = running
    selected = {"benjamini-hochberg": "benjamini", "bonferroni": "bonferroni", "holm": "holm"}.get(correction, "benjamini")
    for index, item in enumerate(items):
        item["benjamini"] = bh[index]
        item["bonferroni"] = min(1.0, item["pvalue"] * count)
        item["holm"] = holm[index]
        item["adjusted_pvalue"] = item[selected]


def enrichment_categories(
    query_mapping: dict[str, dict[str, dict[str, dict]]],
    population_counts: dict[str, int],
    population_term_counts: dict[tuple[str, str], int],
    min_count: int,
    correction: str,
    mode: Literal["direct", "all"],
) -> tuple[list[dict], set[str]]:
    categories: list[dict] = []
    annotated: set[str] = set()
    for namespace, short_label in GO_NAMESPACES.items():
        term_genes: dict[str, set[str]] = {}
        term_info: dict[str, dict] = {}
        for gene_id, namespaces in query_mapping.items():
            for go_id, term in namespaces.get(namespace, {}).items():
                term_genes.setdefault(go_id, set()).add(gene_id)
                term_info[go_id] = term
        query_genes = set().union(*term_genes.values()) if term_genes else set()
        annotated.update(query_genes)
        population_total = population_counts.get(namespace, 0)
        rows: list[dict] = []
        for go_id, genes in term_genes.items():
            count = len(genes)
            population_hits = population_term_counts.get((namespace, go_id), 0)
            if count < min_count or not population_total or not population_hits or not query_genes:
                continue
            query_total = len(query_genes)
            pvalue = hypergeom_tail(count, query_total, population_hits, population_total)
            rows.append(
                {
                    "category": short_label,
                    "namespace": namespace,
                    "go_id": go_id,
                    "term": term_info[go_id]["term"],
                    "level": term_info[go_id].get("level"),
                    "count": count,
                    "percent": count / query_total * 100,
                    "pvalue": pvalue,
                    "list_total": query_total,
                    "pop_hits": population_hits,
                    "pop_total": population_total,
                    "fold_enrichment": (count / query_total) / (population_hits / population_total),
                    "genes": sorted(genes),
                }
            )
        add_multiple_test_corrections(rows, correction)
        rows.sort(key=lambda item: (item["pvalue"], -item["count"], item["go_id"]))
        categories.append(
            {
                "namespace": namespace,
                "label": GO_NAMESPACE_LABELS[namespace],
                "short_label": short_label,
                "query_annotated": len(query_genes),
                "background_annotated": population_total,
                "terms": rows[:1000],
                "term_count": len(rows),
                "mode": mode,
            }
        )
    return categories, annotated


@router.get("/references/{reference_id}/go/{gene_id}")
def get_reference_go(reference_id: str, gene_id: str) -> dict:
    with connect() as conn:
        reference = row(
            conn,
            "SELECT reference_id, display_name, taxon_group FROM reference_species WHERE reference_id = ?",
            (reference_id,),
        )
        if not reference:
            raise HTTPException(status_code=404, detail=f"reference not found: {reference_id}")
        terms = rows(
            conn,
            """
            SELECT
              go.go_id,
              COALESCE(t.name, go.go_name) AS go_name,
              COALESCE(t.namespace, go.namespace) AS go_namespace,
              t.definition,
              t.is_obsolete,
              go.evidence,
              go.source,
              COUNT(*) AS support_count
            FROM reference_go_annotations go
            LEFT JOIN go_terms t ON t.go_id = go.go_id
            WHERE go.reference_id = ? AND go.gene_id = ?
            GROUP BY go.go_id, go_name, go_namespace, t.definition, t.is_obsolete, go.evidence, go.source
            ORDER BY go_namespace, go.go_id, go.evidence
            """,
            (reference_id, gene_id),
        )
    for term in terms:
        term["namespace"] = term.pop("go_namespace", None)
    return {"reference_id": reference_id, "reference": reference, "gene_id": gene_id, "go_count": len(terms), "go_terms": terms}


@router.get("/species/{species_id}/go/status")
def get_species_go_status(species_id: str) -> dict:
    with connect() as conn:
        species = row(
            conn,
            "SELECT species_id, scientific_name, display_name, taxon_group FROM species WHERE species_id = ?",
            (species_id,),
        )
        if not species:
            raise HTTPException(status_code=404, detail=f"species not found: {species_id}")
        data = rows(
            conn,
            """
            SELECT
              h.reference_id,
              r.display_name AS reference_display_name,
              COUNT(DISTINCT h.gene_id) AS genes_with_similarity,
              COUNT(DISTINCT CASE WHEN go.go_id IS NOT NULL THEN h.gene_id END) AS genes_with_inferred_go,
              COUNT(DISTINCT go.go_id) AS distinct_go_terms,
              COUNT(go.go_id) AS inferred_go_links
            FROM similarity_hits h
            JOIN reference_species r ON r.reference_id = h.reference_id
            LEFT JOIN reference_go_annotations go
              ON go.reference_id = h.reference_id
             AND go.gene_id = h.reference_gene_id
            WHERE h.species_id = ?
            GROUP BY h.reference_id
            ORDER BY h.reference_id
            """,
            (species_id,),
        )
    return {"species_id": species_id, "species": species, "go": data}


@router.post("/species/{species_id}/go/enrichment")
def run_species_go_enrichment(species_id: str, request: GoEnrichmentRequest) -> dict:
    submitted = split_gene_ids(request.genes)
    background_submitted = split_gene_ids(request.background)
    if len(submitted) > GO_LIMITS.max_enrichment_genes:
        raise HTTPException(status_code=422, detail=f"Gene list accepts up to {GO_LIMITS.max_enrichment_genes:,} unique IDs.")
    if len(background_submitted) > GO_LIMITS.max_enrichment_background_genes:
        raise HTTPException(status_code=422, detail=f"Background list accepts up to {GO_LIMITS.max_enrichment_background_genes:,} unique IDs.")
    if not submitted:
        raise HTTPException(status_code=422, detail="Enter at least one gene ID.")
    if not request.confidence:
        raise HTTPException(status_code=422, detail="Select at least one homolog confidence level.")

    with connect() as conn:
        species = row(
            conn,
            "SELECT species_id, scientific_name, display_name, taxon_group FROM species WHERE species_id = ?",
            (species_id,),
        )
        if not species:
            raise HTTPException(status_code=404, detail=f"species not found: {species_id}")
        valid, invalid = resolve_species_gene_ids(conn, species_id, submitted)
        if not valid:
            raise HTTPException(status_code=422, detail="No valid gene IDs were found in the selected species.")
        background_valid: list[str] = []
        background_invalid: list[str] = []
        if background_submitted:
            background_valid, background_invalid = resolve_species_gene_ids(conn, species_id, background_submitted)
            if not background_valid:
                raise HTTPException(status_code=422, detail="No valid gene IDs were found in the background list.")

        if background_submitted:
            direct_population_counts, direct_population_term_counts = population_statistics(
                conn,
                species_id,
                request.confidence,
                request.max_similarity_rank,
                background_valid,
                "direct",
            )
            all_population_counts, all_population_term_counts = population_statistics(
                conn,
                species_id,
                request.confidence,
                request.max_similarity_rank,
                background_valid,
                "all",
            )
        else:
            direct_population_counts, direct_population_term_counts = default_population_statistics(
                species_id,
                tuple(sorted(request.confidence)),
                request.max_similarity_rank,
                "direct",
            )
            all_population_counts, all_population_term_counts = default_population_statistics(
                species_id,
                tuple(sorted(request.confidence)),
                request.max_similarity_rank,
                "all",
            )
        direct_query_mapping = gene_go_terms(conn, species_id, valid, request.confidence, request.max_similarity_rank, "direct")
        all_query_mapping = gene_go_terms(conn, species_id, valid, request.confidence, request.max_similarity_rank, "all")

    if background_valid:
        background_set = set(background_valid)
        direct_query_mapping = {gene_id: terms for gene_id, terms in direct_query_mapping.items() if gene_id in background_set}
        all_query_mapping = {gene_id: terms for gene_id, terms in all_query_mapping.items() if gene_id in background_set}
    direct_categories, direct_annotated = enrichment_categories(
        direct_query_mapping,
        direct_population_counts,
        direct_population_term_counts,
        request.min_count,
        request.correction,
        "direct",
    )
    all_categories, all_annotated = enrichment_categories(
        all_query_mapping,
        all_population_counts,
        all_population_term_counts,
        request.min_count,
        request.correction,
        "all",
    )
    valid_in_background = set(all_query_mapping)
    return {
        "found": True,
        "species_id": species_id,
        "species": species,
        "mode": "both",
        "method": "Reference homolog-inferred GO enrichment using one-sided hypergeometric testing.",
        "min_count": request.min_count,
        "correction": request.correction,
        "confidence": request.confidence,
        "max_similarity_rank": request.max_similarity_rank,
        "input_count": len(submitted),
        "valid_count": len(valid),
        "invalid": invalid[:100],
        "background_mode": "custom" if background_submitted else "species_annotated",
        "background_input_count": len(background_submitted),
        "background_valid_count": len(background_valid),
        "background_invalid": background_invalid[:100],
        "query_outside_background_count": len(set(valid) - valid_in_background) if background_submitted else 0,
        "unannotated": sorted(set(valid) - all_annotated)[:200],
        "unannotated_count": len(set(valid) - all_annotated),
        "categories": all_categories,
        "result_sets": {"direct": direct_categories, "all": all_categories},
    }


def semantic_go_figure(request: GoSemanticRequest) -> dict:
    """Run the same GO-Figure reduction path used by SmarLens for enriched GO terms."""
    if not GOFIGURE_PYTHON.is_file() or not GOFIGURE_SCRIPT.is_file():
        return {
            "found": False,
            "error": "GO-Figure is unavailable. Configure ANTARCTIC_GOFIGURE_PYTHON and ANTARCTIC_GOFIGURE_SCRIPT.",
        }

    ontology = {"BP": "bpo", "MF": "mfo", "CC": "cco"}[request.section]
    ontology_name = {
        "BP": "biological_process",
        "MF": "molecular_function",
        "CC": "cellular_component",
    }[request.section]
    cleaned: list[tuple[str, float]] = []
    seen: set[str] = set()
    for term in request.terms:
        go_id = term.go_id.strip()
        pvalue = float(term.pvalue)
        if go_id in seen or not math.isfinite(pvalue):
            continue
        seen.add(go_id)
        cleaned.append((go_id, max(pvalue, 1e-300)))
    cleaned.sort(key=lambda item: item[1])
    truncated = len(cleaned) > request.max_terms
    cleaned = cleaned[:request.max_terms]
    if len(cleaned) < 2:
        return {"found": False, "error": "Select at least two valid GO terms for a semantic plot."}

    gofigure_root = RUNTIME_DIR / "gofigure"
    cache_root = gofigure_root / "cache"
    cache_root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f"{request.section.lower()}_", dir=gofigure_root) as temp_dir:
        output_dir = os.path.abspath(temp_dir)
        input_path = os.path.join(output_dir, "gofigure_input.tsv")
        with open(input_path, "w", encoding="utf-8") as handle:
            handle.write("# GO\tPValue\n")
            for go_id, pvalue in cleaned:
                handle.write(f"{go_id}\t{pvalue:.12g}\n")

        appendix = f"polar_{request.section.lower()}"
        command = [
            str(GOFIGURE_PYTHON), str(GOFIGURE_SCRIPT),
            "-i", input_path,
            "-o", output_dir,
            "-j", "standard",
            "-n", ontology,
            "-q", "svg",
            "-w", appendix,
            "-a", str(request.max_terms),
            "-si", f"{request.si:.2f}",
            "-c", "log10-pval",
            "-p", request.palette,
            "-s", request.bubble_size,
            "-sr", "big" if request.font_size == "large" else request.font_size,
            "-u", request.label_mode,
            "-l", "description",
            "-m", str(request.max_terms),
            "-f", request.font_size,
            "-k", request.legend_size,
            "-t", f"{request.section} semantic GO",
        ]
        env = os.environ.copy()
        env.update({
            "MPLCONFIGDIR": str(cache_root / "matplotlib"),
            "XDG_CACHE_HOME": str(cache_root / "xdg"),
            "HOME": str(cache_root),
        })
        for directory in (env["MPLCONFIGDIR"], env["XDG_CACHE_HOME"]):
            os.makedirs(directory, exist_ok=True)
        try:
            completed = subprocess.run(
                command,
                cwd=output_dir,
                env=env,
                capture_output=True,
                text=True,
                timeout=90,
                check=False,
            )
        except subprocess.TimeoutExpired:
            return {"found": False, "error": "GO-Figure timed out. Try selecting fewer terms."}
        if completed.returncode:
            message = (completed.stderr or completed.stdout or "GO-Figure failed.").strip()
            return {"found": False, "error": message[-1000:]}

        svg_path = os.path.join(output_dir, f"{ontology_name}_{appendix}.svg")
        if not os.path.exists(svg_path):
            candidates = [
                os.path.join(output_dir, filename)
                for filename in os.listdir(output_dir)
                if filename.endswith(".svg")
            ]
            svg_path = candidates[0] if candidates else ""
        if not svg_path or not os.path.exists(svg_path):
            return {"found": False, "error": "GO-Figure finished without creating an SVG file."}

        with open(svg_path, encoding="utf-8") as handle:
            svg = handle.read()
        outputs: dict[str, str] = {"semantic_tsv": "", "semantic_full_tsv": ""}
        for key, filename in (
            ("semantic_tsv", f"{ontology_name}_{appendix}.tsv"),
            ("semantic_full_tsv", f"{ontology_name}_full_table_{appendix}.tsv"),
        ):
            path = os.path.join(output_dir, filename)
            if os.path.exists(path):
                with open(path, encoding="utf-8") as handle:
                    outputs[key] = handle.read()
        fallback_tsv = "GO\tPValue\n" + "\n".join(f"{go_id}\t{pvalue:.12g}" for go_id, pvalue in cleaned) + "\n"
        return {
            "found": True,
            "section": request.section,
            "terms_used": len(cleaned),
            "truncated": truncated,
            "svg": svg,
            "semantic_tsv": outputs["semantic_tsv"] or fallback_tsv,
            "semantic_full_tsv": outputs["semantic_full_tsv"] or fallback_tsv,
            "method": "GO-Figure semantic similarity reduction from enriched GO terms and raw p-values.",
        }


@router.post("/go/semantic")
def run_semantic_go(request: GoSemanticRequest) -> dict:
    result = semantic_go_figure(request)
    if not result.get("found"):
        raise HTTPException(status_code=400, detail=str(result.get("error") or "Semantic GO failed."))
    return result


@router.get("/species/{species_id}/go/{gene_id}")
def get_species_inferred_go(
    species_id: str,
    gene_id: str,
    reference_id: str | None = None,
    max_similarity_rank: int = Query(default=5, ge=1, le=GO_LIMITS.max_similarity_rank),
) -> dict:
    clauses = ["h.species_id = ?", "h.gene_id = ?", "h.rank <= ?"]
    params: list[object] = [species_id, gene_id, max_similarity_rank]
    if reference_id:
        clauses.append("h.reference_id = ?")
        params.append(reference_id)
    with connect() as conn:
        if not row(conn, "SELECT 1 FROM genes WHERE species_id = ? AND gene_id = ?", (species_id, gene_id)):
            raise HTTPException(status_code=404, detail=f"gene not found: {species_id}/{gene_id}")
        terms = rows(
            conn,
            f"""
            SELECT
              go.go_id,
              go.go_name,
              go.namespace,
              t.name AS ontology_name,
              t.namespace AS ontology_namespace,
              t.definition,
              go.evidence,
              go.source,
              h.reference_id,
              rs.display_name AS reference_display_name,
              h.reference_gene_id,
              h.reference_protein_id,
              h.rank AS similarity_rank,
              h.pident,
              h.evalue,
              h.bitscore
            FROM similarity_hits h
            JOIN reference_species rs ON rs.reference_id = h.reference_id
            JOIN reference_go_annotations go
              ON go.reference_id = h.reference_id
             AND go.gene_id = h.reference_gene_id
            LEFT JOIN go_terms t ON t.go_id = go.go_id
            WHERE {' AND '.join(clauses)}
            ORDER BY h.rank, go.namespace, go.go_id, go.evidence
            """,
            params,
        )

    grouped: dict[tuple, dict] = {}
    for term in terms:
        display = term_display(term)
        key = (display["go_id"], display["go_name"], display["namespace"])
        item = grouped.setdefault(
            key,
            {
                "go_id": display["go_id"],
                "go_name": display["go_name"],
                "namespace": display["namespace"],
                "definition": term.get("definition"),
                "support_count": 0,
                "support": [],
            },
        )
        item["support_count"] += 1
        item["support"].append(
            {
                "reference_id": term["reference_id"],
                "reference_display_name": term["reference_display_name"],
                "reference_gene_id": term["reference_gene_id"],
                "reference_protein_id": term["reference_protein_id"],
                "similarity_rank": term["similarity_rank"],
                "pident": term["pident"],
                "evalue": term["evalue"],
                "bitscore": term["bitscore"],
                "evidence": term["evidence"],
                "source": term["source"],
            }
        )
    go_terms = sorted(grouped.values(), key=lambda item: ((item["namespace"] or ""), item["go_id"]))
    return {
        "species_id": species_id,
        "gene_id": gene_id,
        "inference_method": "similarity_reference_go",
        "max_similarity_rank": max_similarity_rank,
        "go_count": len(go_terms),
        "go_terms": go_terms,
    }


@router.get("/go/terms/{go_id}")
def get_go_term(go_id: str, include_ancestors: bool = True) -> dict:
    with connect() as conn:
        term = row(
            conn,
            "SELECT go_id, name, namespace, definition, is_obsolete, metadata_json FROM go_terms WHERE go_id = ?",
            (go_id,),
        )
        if not term:
            alias = row(conn, "SELECT primary_go_id FROM go_alt_ids WHERE alt_go_id = ?", (go_id,))
            if alias:
                term = row(
                    conn,
                    "SELECT go_id, name, namespace, definition, is_obsolete, metadata_json FROM go_terms WHERE go_id = ?",
                    (alias["primary_go_id"],),
                )
        if not term:
            raise HTTPException(status_code=404, detail=f"GO term not found: {go_id}")
        parents = rows(
            conn,
            """
            SELECT e.parent_go_id AS go_id, t.name, t.namespace, e.relation_type
            FROM go_edges e
            JOIN go_terms t ON t.go_id = e.parent_go_id
            WHERE e.child_go_id = ?
            ORDER BY e.relation_type, e.parent_go_id
            """,
            (term["go_id"],),
        )
    term["metadata"] = parse_json(term.pop("metadata_json", None))
    payload = {"term": term, "parents": parents}
    if include_ancestors:
        payload["ancestors"] = ancestors_for(term["go_id"])
    return payload


@router.get("/species/{species_id}/go/{gene_id}/semantic")
def get_species_go_semantic(
    species_id: str,
    gene_id: str,
    reference_id: str | None = None,
    max_similarity_rank: int = Query(default=5, ge=1, le=GO_LIMITS.max_similarity_rank),
    similarity_threshold: float = Query(default=0.7, ge=0.0, le=1.0),
) -> dict:
    inferred = get_species_inferred_go(species_id, gene_id, reference_id, max_similarity_rank)
    terms = [item for item in inferred["go_terms"] if item.get("go_id")]
    terms.sort(key=lambda item: (-item["support_count"], item["go_id"]))
    clusters: list[dict] = []
    assigned: set[str] = set()
    for term in terms:
        if term["go_id"] in assigned:
            continue
        members = [term]
        assigned.add(term["go_id"])
        for candidate in terms:
            if candidate["go_id"] in assigned:
                continue
            if candidate.get("namespace") != term.get("namespace"):
                continue
            similarity = wang_similarity(term["go_id"], candidate["go_id"])
            if similarity >= similarity_threshold:
                member = {**candidate, "semantic_similarity_to_representative": similarity}
                members.append(member)
                assigned.add(candidate["go_id"])
        representative = max(members, key=lambda item: (item["support_count"], item["go_id"]))
        clusters.append(
            {
                "representative": representative,
                "member_count": len(members),
                "members": members,
            }
        )
    clusters.sort(key=lambda item: (-item["representative"]["support_count"], item["representative"]["go_id"]))
    return {
        "species_id": species_id,
        "gene_id": gene_id,
        "method": "wang_similarity",
        "similarity_threshold": similarity_threshold,
        "input_go_count": len(terms),
        "cluster_count": len(clusters),
        "clusters": clusters,
    }
