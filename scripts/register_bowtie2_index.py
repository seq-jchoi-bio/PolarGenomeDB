#!/usr/bin/env python3
"""Validate and register an externally built Bowtie2 index; never build it."""
import argparse
from contextlib import closing
import hashlib
import json
import os
from pathlib import Path
import shutil
import sqlite3
import subprocess
import sys

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.build_search_indexes import init_db, register_asset
from backend.services.bowtie_indexes import index_format


def inspect_lengths(output):
    lengths = {}
    for line in output.splitlines():
        fields = line.split('\t')
        if fields[0].startswith('Sequence-'):
            if len(fields) != 3 or fields[1] in lengths:
                raise ValueError('Invalid or duplicate sequence in Bowtie2 summary')
            lengths[fields[1]] = int(fields[2])
    if not lengths:
        raise ValueError('Bowtie2 index contains no reference sequences')
    return lengths


def register(manifest_path, db, validate_only=False):
    manifest_path = Path(manifest_path).resolve()
    manifest = yaml.safe_load(manifest_path.read_text())
    species = manifest['species_id']
    primer = manifest.get('primer') or {}
    configured = primer.get('bowtie2_index_prefix')
    if not Path(db).is_file():
        raise FileNotFoundError(db)
    with closing(sqlite3.connect(f'file:{Path(db).resolve()}?mode=ro', uri=True)) as conn:
        source = conn.execute("SELECT path FROM genome_assets WHERE species_id=? AND asset_type='fasta'", (species,)).fetchone()
        expected = dict(conn.execute('SELECT seqid,length FROM chromosomes WHERE species_id=?', (species,)))
    if not source or not expected:
        raise ValueError(f'Build species first: {species}')
    fasta = Path(source[0])
    if not fasta.is_absolute():
        fasta = ROOT / fasta
    manifest_fasta = Path(manifest['genome']['fasta']).expanduser()
    if not manifest_fasta.is_absolute():
        manifest_fasta = manifest_path.parent / manifest_fasta
    if manifest_fasta.resolve() != fasta.resolve():
        raise ValueError('Manifest FASTA differs from the registered genome')
    prefix = None
    metadata = {}
    version = ''
    if configured:
        prefix = Path(configured).expanduser()
        if not prefix.is_absolute():
            prefix = manifest_path.parent / prefix
        prefix = prefix.resolve()
        metadata['format'] = index_format(prefix)
        tool = os.environ.get('BOWTIE2_INSPECT') or shutil.which('bowtie2-inspect')
        if not tool:
            raise FileNotFoundError('bowtie2-inspect is required (PATH or BOWTIE2_INSPECT)')
        inspected = subprocess.run([tool, '-s', str(prefix)], check=True, capture_output=True, text=True, timeout=120)
        if inspect_lengths(inspected.stdout) != expected:
            raise ValueError('Bowtie2 chromosome names/lengths do not match the registered genome')
        version = subprocess.run([tool, '--version'], check=True, capture_output=True, text=True, timeout=30).stdout.splitlines()[0]
        metadata['validation'] = 'chromosome_names_and_lengths'
        checksum = primer.get('genome_sha256')
        if checksum:
            with fasta.open('rb') as handle:
                actual = hashlib.file_digest(handle, 'sha256').hexdigest()
            if actual.lower() != str(checksum).lower():
                raise ValueError('Genome FASTA SHA256 does not match the build-provenance checksum')
            metadata['genome_sha256'] = actual
        metadata['sequence_identity_verified'] = False
    status = 'ready' if prefix else 'disabled'
    if not validate_only:
        with closing(sqlite3.connect(db)) as conn, conn:
            init_db(conn)
            register_asset(conn, species, 'bowtie2_genome', fasta, prefix or '', 'bowtie2', version, status, metadata)
    return {'species_id': species, 'status': status, 'index_prefix': str(prefix) if prefix else None,
            'validation_only': validate_only, 'metadata': metadata}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', required=True)
    parser.add_argument('--db', default=str(ROOT / 'data/antarctic_plant.db'))
    parser.add_argument('--validate-only', action='store_true')
    args = parser.parse_args()
    print(json.dumps(register(args.manifest, args.db, args.validate_only), indent=2))


if __name__ == '__main__':
    main()
