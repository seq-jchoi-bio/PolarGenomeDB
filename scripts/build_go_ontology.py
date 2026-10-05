#!/usr/bin/env python3
"""Load go-basic.obo into SQLite ontology tables."""

from __future__ import annotations

import argparse
import json
import sqlite3
import sys
from pathlib import Path
from collections import deque
from typing import Any, Dict, Iterator, List, Tuple


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = ROOT / "data" / "antarctic_plant.db"
DEFAULT_OBO = ROOT / "reference" / "go-basic.obo"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--obo", default=str(DEFAULT_OBO), help="go-basic.obo path.")
    parser.add_argument("--db", default=str(DEFAULT_DB), help="SQLite DB path.")
    parser.add_argument("--replace", action="store_true", help="Replace existing GO ontology rows.")
    return parser.parse_args()


def parse_obo(path: Path) -> Iterator[Dict[str, Any]]:
    current: Dict[str, Any] | None = None
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line == "[Term]":
                if current and current.get("id"):
                    yield current
                current = {"parents": [], "relations": [], "alt_ids": [], "subsets": []}
                continue
            if line.startswith("["):
                if current and current.get("id"):
                    yield current
                current = None
                continue
            if current is None or not line or line.startswith("!"):
                continue
            key, sep, value = line.partition(": ")
            if not sep:
                continue
            if key == "id":
                current["id"] = value
            elif key == "name":
                current["name"] = value
            elif key == "namespace":
                current["namespace"] = value
            elif key == "def":
                current["definition"] = value
            elif key == "is_obsolete":
                current["is_obsolete"] = value == "true"
            elif key == "alt_id":
                current["alt_ids"].append(value)
            elif key == "subset":
                current["subsets"].append(value)
            elif key == "is_a":
                current["parents"].append(("is_a", value.split(" ! ", 1)[0]))
            elif key == "relationship":
                parts = value.split()
                if len(parts) >= 2 and parts[0] in {"part_of", "regulates", "positively_regulates", "negatively_regulates"}:
                    current["parents"].append((parts[0], parts[1]))
    if current and current.get("id"):
        yield current


def init_db(conn: sqlite3.Connection) -> None:
    conn.executescript(
        """
        CREATE TABLE IF NOT EXISTS go_terms (
            go_id TEXT PRIMARY KEY,
            name TEXT NOT NULL,
            namespace TEXT NOT NULL,
            definition TEXT,
            is_obsolete INTEGER NOT NULL DEFAULT 0,
            metadata_json TEXT
        );

        CREATE TABLE IF NOT EXISTS go_edges (
            child_go_id TEXT NOT NULL,
            parent_go_id TEXT NOT NULL,
            relation_type TEXT NOT NULL,
            PRIMARY KEY (child_go_id, parent_go_id, relation_type),
            FOREIGN KEY (child_go_id) REFERENCES go_terms(go_id) ON DELETE CASCADE,
            FOREIGN KEY (parent_go_id) REFERENCES go_terms(go_id) ON DELETE CASCADE
        );

        CREATE TABLE IF NOT EXISTS go_alt_ids (
            alt_go_id TEXT PRIMARY KEY,
            primary_go_id TEXT NOT NULL,
            FOREIGN KEY (primary_go_id) REFERENCES go_terms(go_id) ON DELETE CASCADE
        );

        CREATE TABLE IF NOT EXISTS go_ancestors (
            go_id TEXT NOT NULL,
            ancestor_go_id TEXT NOT NULL,
            level INTEGER,
            distance INTEGER NOT NULL,
            PRIMARY KEY (go_id, ancestor_go_id),
            FOREIGN KEY (go_id) REFERENCES go_terms(go_id) ON DELETE CASCADE,
            FOREIGN KEY (ancestor_go_id) REFERENCES go_terms(go_id) ON DELETE CASCADE
        );

        CREATE INDEX IF NOT EXISTS idx_go_terms_namespace ON go_terms(namespace, go_id);
        CREATE INDEX IF NOT EXISTS idx_go_edges_child ON go_edges(child_go_id);
        CREATE INDEX IF NOT EXISTS idx_go_edges_parent ON go_edges(parent_go_id);
        CREATE INDEX IF NOT EXISTS idx_go_ancestors_ancestor ON go_ancestors(ancestor_go_id, go_id);
        CREATE INDEX IF NOT EXISTS idx_go_ancestors_level ON go_ancestors(level, ancestor_go_id);
        """
    )


def build_ancestor_rows(terms: List[Dict[str, Any]], edges: List[Tuple[str, str, str]]) -> List[Tuple[str, str, int | None, int]]:
    """Materialize transitive GO ancestry and shortest root depth for fast enrichment queries."""
    term_ids = {term["id"] for term in terms}
    parents: Dict[str, List[str]] = {term_id: [] for term_id in term_ids}
    children: Dict[str, List[str]] = {term_id: [] for term_id in term_ids}
    for child, parent, relation in edges:
        parents[child].append(parent)
        if relation in {"is_a", "part_of"}:
            children[parent].append(child)

    roots = {"GO:0008150", "GO:0003674", "GO:0005575"}
    levels: Dict[str, int] = {}
    queue: deque[Tuple[str, int]] = deque((root, 0) for root in roots if root in term_ids)
    while queue:
        term_id, level = queue.popleft()
        if term_id in levels and levels[term_id] <= level:
            continue
        levels[term_id] = level
        queue.extend((child, level + 1) for child in children.get(term_id, []))

    memo: Dict[str, Dict[str, int]] = {}
    visiting: set[str] = set()

    def ancestor_distances(term_id: str) -> Dict[str, int]:
        cached = memo.get(term_id)
        if cached is not None:
            return cached
        if term_id in visiting:
            return {term_id: 0}
        visiting.add(term_id)
        distances = {term_id: 0}
        for parent in parents.get(term_id, []):
            for ancestor_id, distance in ancestor_distances(parent).items():
                next_distance = distance + 1
                if next_distance < distances.get(ancestor_id, next_distance + 1):
                    distances[ancestor_id] = next_distance
        visiting.remove(term_id)
        memo[term_id] = distances
        return distances

    result: List[Tuple[str, str, int | None, int]] = []
    for term_id in sorted(term_ids):
        for ancestor_id, distance in ancestor_distances(term_id).items():
            result.append((term_id, ancestor_id, levels.get(ancestor_id), distance))
    return result


def build_ontology(args: argparse.Namespace) -> Dict[str, Any]:
    obo_path = Path(args.obo).resolve()
    if not obo_path.exists():
        raise FileNotFoundError(obo_path)
    terms = list(parse_obo(obo_path))
    db_path = Path(args.db).resolve()
    conn = sqlite3.connect(db_path)
    init_db(conn)
    try:
        conn.execute("BEGIN")
        if args.replace:
            conn.execute("DELETE FROM go_alt_ids")
            conn.execute("DELETE FROM go_edges")
            conn.execute("DELETE FROM go_terms")
        for term in terms:
            conn.execute(
                """
                INSERT OR REPLACE INTO go_terms
                (go_id, name, namespace, definition, is_obsolete, metadata_json)
                VALUES (?, ?, ?, ?, ?, ?)
                """,
                (
                    term["id"],
                    term.get("name", ""),
                    term.get("namespace", ""),
                    term.get("definition"),
                    1 if term.get("is_obsolete") else 0,
                    json.dumps({"subsets": term.get("subsets", [])}, sort_keys=True),
                ),
            )
        term_ids = {term["id"] for term in terms}
        edges: List[Tuple[str, str, str]] = []
        alt_ids: List[Tuple[str, str]] = []
        for term in terms:
            for relation, parent in term.get("parents", []):
                if parent in term_ids:
                    edges.append((term["id"], parent, relation))
            for alt_id in term.get("alt_ids", []):
                alt_ids.append((alt_id, term["id"]))
        conn.executemany("INSERT OR REPLACE INTO go_edges VALUES (?, ?, ?)", edges)
        conn.executemany("INSERT OR REPLACE INTO go_alt_ids VALUES (?, ?)", alt_ids)
        ancestor_rows = build_ancestor_rows(terms, edges)
        conn.execute("DELETE FROM go_ancestors")
        conn.executemany(
            "INSERT INTO go_ancestors (go_id, ancestor_go_id, level, distance) VALUES (?, ?, ?, ?)",
            ancestor_rows,
        )
        conn.commit()
        return {
            "obo": str(obo_path),
            "terms": len(terms),
            "edges": len(edges),
            "alt_ids": len(alt_ids),
            "ancestor_rows": len(ancestor_rows),
            "obsolete_terms": sum(1 for term in terms if term.get("is_obsolete")),
        }
    except Exception:
        conn.rollback()
        raise
    finally:
        conn.close()


def main() -> int:
    args = parse_args()
    try:
        print(json.dumps(build_ontology(args), indent=2, sort_keys=True))
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
