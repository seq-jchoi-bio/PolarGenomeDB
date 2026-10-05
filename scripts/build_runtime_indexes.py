#!/usr/bin/env python3
"""Export lightweight runtime JSON indexes from the AntarcticPlantDB SQLite DB."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import sqlite3
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = ROOT / "data" / "antarctic_plant.db"
DEFAULT_OUT = ROOT / "data" / "runtime"

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", default=str(DEFAULT_DB), help="SQLite DB path.")
    parser.add_argument("--out", default=str(DEFAULT_OUT), help="Runtime JSON output directory.")
    parser.add_argument("--species", help="Optional single species_id to export.")
    parser.add_argument("--pretty", action="store_true", help="Write indented JSON.")
    parser.add_argument("--skip-go-population", action="store_true", help="Skip the default GO enrichment population cache.")
    return parser.parse_args()


def rows(conn: sqlite3.Connection, query: str, params: Iterable[Any] = ()) -> List[Dict[str, Any]]:
    cursor = conn.execute(query, tuple(params))
    columns = [column[0] for column in cursor.description]
    return [dict(zip(columns, row)) for row in cursor.fetchall()]


def one(conn: sqlite3.Connection, query: str, params: Iterable[Any] = ()) -> Dict[str, Any]:
    result = rows(conn, query, params)
    return result[0] if result else {}


def parse_json(value: Optional[str]) -> Dict[str, Any]:
    if not value:
        return {}
    try:
        parsed = json.loads(value)
    except json.JSONDecodeError:
        return {"raw": value}
    return parsed if isinstance(parsed, dict) else {"value": parsed}


def relpath(path: Optional[str], base: Path) -> Optional[str]:
    if not path:
        return None
    resolved = Path(path)
    try:
        return str(resolved.resolve().relative_to(base.resolve()))
    except ValueError:
        return str(resolved)


def write_json(path: Path, payload: Dict[str, Any], pretty: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    kwargs = {"ensure_ascii": False, "sort_keys": True}
    if pretty:
        kwargs["indent"] = 2
    path.write_text(json.dumps(payload, **kwargs) + "\n", encoding="utf-8")


def count_value_types(conn: sqlite3.Connection, species_id: str) -> List[Dict[str, Any]]:
    return rows(
        conn,
        """
        SELECT value_type, COUNT(DISTINCT gene_id) AS genes, COUNT(DISTINCT sample_id) AS samples, COUNT(*) AS value_count
        FROM expression_counts
        WHERE species_id = ?
        GROUP BY value_type
        ORDER BY value_type
        """,
        (species_id,),
    )


def feature_counts(conn: sqlite3.Connection, species_id: str) -> Dict[str, int]:
    data = one(
        conn,
        """
        SELECT
          (SELECT COUNT(*) FROM genes WHERE species_id = ?) AS genes,
          (SELECT COUNT(*) FROM transcripts WHERE species_id = ?) AS transcripts,
          (SELECT COUNT(*) FROM proteins WHERE species_id = ?) AS proteins,
          (SELECT COUNT(*) FROM features WHERE species_id = ?) AS features,
          (SELECT COUNT(*) FROM samples WHERE species_id = ?) AS samples,
          (SELECT COUNT(*) FROM track_assets WHERE species_id = ? AND track_type = 'bigwig' AND status = 'ready') AS bigwig_tracks,
          (SELECT COUNT(*) FROM search_assets WHERE species_id = ? AND status = 'ready') AS search_assets
        """,
        (species_id, species_id, species_id, species_id, species_id, species_id, species_id),
    )
    return {key: int(value or 0) for key, value in data.items()}


def export_species(conn: sqlite3.Connection, species: Dict[str, Any], out_dir: Path, base: Path, pretty: bool) -> Dict[str, Any]:
    species_id = species["species_id"]
    species_dir = out_dir / species_id
    now = dt.datetime.now(dt.timezone.utc).isoformat()

    chromosomes = rows(
        conn,
        """
        SELECT seqid, length
        FROM chromosomes
        WHERE species_id = ?
        ORDER BY length DESC, seqid
        """,
        (species_id,),
    )
    genome_assets = rows(
        conn,
        "SELECT asset_type, path FROM genome_assets WHERE species_id = ? ORDER BY asset_type",
        (species_id,),
    )
    annotation_assets = rows(
        conn,
        "SELECT format, path FROM annotation_assets WHERE species_id = ? ORDER BY format, path",
        (species_id,),
    )
    search_assets = rows(
        conn,
        """
        SELECT asset_type, tool, tool_version, source_fasta, index_prefix, created_at, status, metadata_json
        FROM search_assets
        WHERE species_id = ?
        ORDER BY tool, asset_type
        """,
        (species_id,),
    )
    samples = rows(
        conn,
        """
        SELECT sample_id, display_name, tissue, condition, treatment, timepoint, replicate, source_accession, metadata_json
        FROM samples
        WHERE species_id = ?
        ORDER BY sample_id
        """,
        (species_id,),
    )
    tracks = rows(
        conn,
        """
        SELECT sample_id, track_type, path, source_bam, bin_size, normalization, created_at, status, metadata_json
        FROM track_assets
        WHERE species_id = ?
        ORDER BY sample_id, track_type, bin_size, normalization
        """,
        (species_id,),
    )

    for item in genome_assets:
        item["path"] = relpath(item["path"], base)
    for item in annotation_assets:
        item["path"] = relpath(item["path"], base)
    for item in search_assets:
        item["source_fasta"] = relpath(item["source_fasta"], base)
        item["index_prefix"] = relpath(item["index_prefix"], base)
        item["metadata"] = parse_json(item.pop("metadata_json", None))
    for item in samples:
        item["metadata"] = parse_json(item.pop("metadata_json", None))
    for item in tracks:
        item["path"] = relpath(item["path"], base)
        item["source_bam"] = relpath(item["source_bam"], base)
        item["metadata"] = parse_json(item.pop("metadata_json", None))

    overview = {
        "species_id": species_id,
        "scientific_name": species["scientific_name"],
        "display_name": species["display_name"],
        "taxon_group": species["taxon_group"],
        "assembly_version": species["assembly_version"],
        "counts": feature_counts(conn, species_id),
        "expression_value_types": count_value_types(conn, species_id),
        "runtime_manifest": f"{species_id}/manifest.json",
    }

    manifest = {
        "generated_at": now,
        "species": overview,
        "chromosomes": chromosomes,
        "assets": {
            "genome": genome_assets,
            "annotation": annotation_assets,
            "search": search_assets,
            "tracks": tracks,
        },
        "samples": samples,
    }
    browser_tracks = {
        "generated_at": now,
        "species_id": species_id,
        "chromosomes": chromosomes,
        "tracks": [
            {
                "track_id": f"{species_id}:{track['sample_id']}:{track['track_type']}:{track.get('bin_size') or 'native'}:{track.get('normalization') or 'none'}",
                "name": track["sample_id"],
                "type": track["track_type"],
                "path": track["path"],
                "sample_id": track["sample_id"],
                "bin_size": track["bin_size"],
                "normalization": track["normalization"],
                "status": track["status"],
            }
            for track in tracks
            if track["status"] == "ready"
        ],
    }
    expression = {
        "generated_at": now,
        "species_id": species_id,
        "samples": samples,
        "value_types": count_value_types(conn, species_id),
        "tracks": [track for track in tracks if track["track_type"] == "bigwig" and track["status"] == "ready"],
    }
    search = {
        "generated_at": now,
        "species_id": species_id,
        "assets": search_assets,
    }

    write_json(species_dir / "manifest.json", manifest, pretty)
    write_json(species_dir / "genome_browser_tracks.json", browser_tracks, pretty)
    write_json(species_dir / "expression_tracks.json", expression, pretty)
    write_json(species_dir / "search_assets.json", search, pretty)
    return overview


def export_go_population_cache(conn: sqlite3.Connection, species_id: str, out_dir: Path, pretty: bool) -> Dict[str, Any]:
    """Build the small, reusable population statistics needed by default GO enrichment."""
    from backend.routes.go import population_statistics

    confidence = ["high", "medium"]
    max_similarity_rank = 5
    modes: Dict[str, Any] = {}
    for mode in ("direct", "all"):
        population_counts, term_counts = population_statistics(
            conn,
            species_id,
            confidence,
            max_similarity_rank,
            None,
            mode,
        )
        modes[mode] = {
            "population_counts": population_counts,
            "term_counts": [
                {"namespace": namespace, "go_id": go_id, "gene_count": count}
                for (namespace, go_id), count in sorted(term_counts.items())
            ],
        }
    payload = {
        "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "species_id": species_id,
        "confidence": confidence,
        "max_similarity_rank": max_similarity_rank,
        "modes": modes,
    }
    write_json(out_dir / species_id / "go_population.json", payload, pretty)
    return {
        "path": f"{species_id}/go_population.json",
        "direct_terms": len(modes["direct"]["term_counts"]),
        "all_terms": len(modes["all"]["term_counts"]),
    }


def build_runtime_indexes(args: argparse.Namespace) -> Dict[str, Any]:
    db_path = Path(args.db).resolve()
    out_dir = Path(args.out).resolve()
    if not db_path.exists():
        raise FileNotFoundError(f"SQLite DB not found: {db_path}")

    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    try:
        species_query = """
            SELECT species_id, scientific_name, display_name, taxon_group, assembly_version
            FROM species
        """
        params: tuple[Any, ...] = ()
        if args.species:
            species_query += " WHERE species_id = ?"
            params = (args.species,)
        species_query += " ORDER BY species_id"
        species_rows = rows(conn, species_query, params)
        if args.species and not species_rows:
            raise ValueError(f"species not found: {args.species}")

        generated_at = dt.datetime.now(dt.timezone.utc).isoformat()
        overview = [export_species(conn, species, out_dir, ROOT, args.pretty) for species in species_rows]
        go_population = {}
        if not args.skip_go_population:
            go_population = {
                species["species_id"]: export_go_population_cache(conn, species["species_id"], out_dir, args.pretty)
                for species in species_rows
            }
        index = {
            "generated_at": generated_at,
            "db_path": relpath(str(db_path), ROOT),
            "species_count": len(overview),
            "species": overview,
        }
        write_json(out_dir / "species.json", index, args.pretty)
        return {
            "runtime_dir": str(out_dir),
            "species_count": len(overview),
            "species": [item["species_id"] for item in overview],
            "go_population": go_population,
        }
    finally:
        conn.close()


def main() -> int:
    args = parse_args()
    try:
        summary = build_runtime_indexes(args)
    except Exception as exc:
        print(f"ERROR: {exc}", flush=True)
        return 1
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
