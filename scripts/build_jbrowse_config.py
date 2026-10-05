#!/usr/bin/env python3
"""Build draft JBrowse 2 config files from runtime manifests."""

from __future__ import annotations

import argparse
import datetime as dt
import json
from pathlib import Path
from typing import Any, Dict, List, Optional


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RUNTIME = ROOT / "data" / "runtime"
DEFAULT_OUT = ROOT / "data" / "jbrowse"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runtime", default=str(DEFAULT_RUNTIME), help="Runtime JSON directory.")
    parser.add_argument("--out", default=str(DEFAULT_OUT), help="JBrowse config output directory.")
    parser.add_argument("--species", help="Optional single species_id.")
    parser.add_argument("--base-url", default="/files", help="URL prefix that serves project-root-relative files.")
    parser.add_argument("--pretty", action="store_true", help="Write indented JSON.")
    return parser.parse_args()


def read_json(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: Dict[str, Any], pretty: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    kwargs = {"ensure_ascii": False, "sort_keys": True}
    if pretty:
        kwargs["indent"] = 2
    path.write_text(json.dumps(payload, **kwargs) + "\n", encoding="utf-8")


def url_for(path: Optional[str], base_url: str) -> Optional[str]:
    if not path:
        return None
    if path.startswith(("http://", "https://", "/")):
        return path
    return f"{base_url.rstrip('/')}/{path}"


def first_asset(assets: List[Dict[str, Any]], asset_type: str) -> Optional[Dict[str, Any]]:
    for asset in assets:
        if asset.get("asset_type") == asset_type:
            return asset
    return None


def preferred_annotation_asset(annotations: List[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    for annotation_format in ("gff3_tabix", "gtf_tabix"):
        for annotation in annotations:
            if annotation.get("format") == annotation_format:
                return annotation
    return None


def annotation_track(species_id: str, annotation: Dict[str, Any], annotations: List[Dict[str, Any]], assembly_name: str, base_url: str) -> Optional[Dict[str, Any]]:
    fmt = annotation["format"]
    if fmt not in {"gff3_tabix", "gtf_tabix"}:
        return None
    index_base = fmt.rsplit("_", 1)[0]
    index = next((item for item in annotations if item["format"] in {f"{index_base}_tabix_tbi", f"{index_base}_tabix_csi"}), None)
    if not index:
        return None
    adapter_type = "Gff3TabixAdapter" if fmt == "gff3_tabix" else "GtfTabixAdapter"
    adapter: Dict[str, Any] = {
        "type": adapter_type,
        "gffGzLocation": {
            "uri": url_for(annotation["path"], base_url),
            "locationType": "UriLocation",
        },
        "index": {
            "indexType": "CSI" if index["format"].endswith("_csi") else "TBI",
            "location": {
                "uri": url_for(index["path"], base_url),
                "locationType": "UriLocation",
            },
        },
    }
    if fmt == "gtf_tabix":
        adapter["aggregateField"] = "gene_id"
    return {
        "type": "FeatureTrack",
        "trackId": f"{species_id}_gene_annotation",
        "name": "Gene annotation",
        "assemblyNames": [assembly_name],
        "adapter": adapter,
    }


def build_species_config(species_id: str, manifest: Dict[str, Any], base_url: str) -> Dict[str, Any]:
    species = manifest["species"]
    genome_assets = manifest["assets"]["genome"]
    annotation_assets = manifest["assets"]["annotation"]
    fasta = first_asset(genome_assets, "fasta")
    fai = first_asset(genome_assets, "fai")
    assembly_name = species["display_name"] or species_id
    tracks: List[Dict[str, Any]] = []

    for track in manifest["assets"]["tracks"]:
        if track["track_type"] != "bigwig" or track["status"] != "ready":
            continue
        tracks.append(
            {
                "type": "QuantitativeTrack",
                "trackId": f"{species_id}_{track['sample_id']}_bigwig",
                "name": f"{track['sample_id']} {track['normalization']} bin{track['bin_size']}",
                "assemblyNames": [assembly_name],
                "adapter": {
                    "type": "BigWigAdapter",
                    "bigWigLocation": {
                        "uri": url_for(track["path"], base_url),
                        "locationType": "UriLocation",
                    },
                },
            }
        )

    annotation = preferred_annotation_asset(annotation_assets)
    if annotation:
        track = annotation_track(species_id, annotation, annotation_assets, assembly_name, base_url)
        if track:
            tracks.insert(0, track)

    config = {
        "configuration": {},
        "assemblies": [
            {
                "name": assembly_name,
                "sequence": {
                    "type": "ReferenceSequenceTrack",
                    "trackId": f"{species_id}_reference_sequence",
                    "adapter": {
                        "type": "IndexedFastaAdapter",
                        "fastaLocation": {
                            "uri": url_for(fasta["path"], base_url) if fasta else None,
                            "locationType": "UriLocation",
                        },
                        "faiLocation": {
                            "uri": url_for(fai["path"], base_url) if fai else None,
                            "locationType": "UriLocation",
                        },
                    },
                },
                "aliases": [species_id, species["scientific_name"]],
            }
        ],
        "tracks": tracks,
        "metadata": {
            "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(),
            "species_id": species_id,
            "annotation_assets": annotation_assets,
            "notes": [
                "FASTA, tabix-indexed annotation, and bigWig tracks are represented directly when available.",
            ],
        },
    }
    return config


def main() -> int:
    args = parse_args()
    runtime_dir = Path(args.runtime).resolve()
    out_dir = Path(args.out).resolve()
    species_index = read_json(runtime_dir / "species.json")
    species_ids = [args.species] if args.species else [item["species_id"] for item in species_index["species"]]
    exported = []
    for species_id in species_ids:
        manifest_path = runtime_dir / species_id / "manifest.json"
        if not manifest_path.exists():
            raise FileNotFoundError(manifest_path)
        config = build_species_config(species_id, read_json(manifest_path), args.base_url)
        write_json(out_dir / species_id / "config.json", config, args.pretty)
        exported.append(species_id)
    write_json(
        out_dir / "index.json",
        {"generated_at": dt.datetime.now(dt.timezone.utc).isoformat(), "species": exported},
        args.pretty,
    )
    print(json.dumps({"jbrowse_dir": str(out_dir), "species": exported}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
