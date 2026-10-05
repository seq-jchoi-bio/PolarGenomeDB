#!/usr/bin/env python3
"""Plan or run a species build in an isolated output directory."""
import argparse
from contextlib import closing
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import sqlite3
import subprocess
import sys

import yaml

ROOT = Path(__file__).resolve().parents[1]


def make_plan(manifest_path, output, threads, replace=False):
    manifest_path = Path(manifest_path).resolve()
    manifest = yaml.safe_load(manifest_path.read_text())
    species = manifest['species_id']
    if not re.fullmatch(r'[A-Za-z0-9_-]+', species):
        raise ValueError('species_id must use ASCII letters, digits, underscore or hyphen')
    if manifest.get('taxon_group') not in ('monocot', 'dicot'):
        raise ValueError('The full reference pipeline currently requires monocot or dicot')
    output = Path(output).resolve()
    db = str(output / 'antarctic_plant.db')
    assets = str(output / 'species_assets')
    runtime = str(output / 'runtime')
    common = ['--manifest', str(manifest_path), '--db', db]
    def step(name, script, args):
        return {'stage': name, 'command': [sys.executable, str(ROOT / 'scripts' / script), *args]}
    plan = [step('species', 'build_species.py', common + ['--assets-dir', assets, '--reports-dir', str(output / 'reports')] + (['--replace'] if replace else [])),
            step('search', 'build_search_indexes.py', common + ['--assets-dir', assets, '--threads', str(threads), '--require', 'all']),
            step('bowtie2', 'register_bowtie2_index.py', common),
            step('pfam', 'build_pfam_domains.py', ['--species', species, '--cpu', str(threads)]),
            step('similarity', 'build_similarity_hits.py', ['--species', species, '--db', db, '--threads', str(threads)])]
    expression = manifest.get('expression') or {}
    if expression:
        args = common + ['--assets-dir', assets, '--threads', str(threads)]
        if expression.get('count_matrix'):
            args += ['--load-counts']
        if expression.get('bam_dir'):
            args += ['--build-indexes', '--build-bigwig', '--bin-size', '25', '--normalization', 'CPM']
        plan.append(step('expression', 'build_expression.py', args))
    plan += [step('browser_assets', 'build_jbrowse_assets.py', ['--db', db, '--assets-dir', assets, '--species', species]),
             # Export every species so species.json never drops existing species.
             step('runtime', 'build_runtime_indexes.py', ['--db', db, '--out', runtime, '--pretty']),
             step('jbrowse', 'build_jbrowse_config.py', ['--runtime', runtime, '--out', str(output / 'jbrowse'), '--pretty'])]
    return manifest, plan


def check_inputs(manifest, path):
    pairs = [('genome', 'fasta'), ('annotation', 'file'), ('protein', 'fasta')]
    pairs += [('expression', key) for key in ('sample_metadata', 'count_matrix', 'bam_dir') if (manifest.get('expression') or {}).get(key)]
    for section, key in pairs:
        value = manifest[section][key]
        file = Path(value).expanduser()
        if not file.is_absolute():
            file = path.parent / file
        if not file.exists():
            raise FileNotFoundError(f'{section}.{key}: {file}')
    expression = manifest.get('expression') or {}
    if expression and not expression.get('sample_metadata'):
        raise ValueError('expression.sample_metadata is required when expression is configured')


def verify(db, species):
    with closing(sqlite3.connect(f'file:{db}?mode=ro', uri=True)) as conn:
        if conn.execute('PRAGMA quick_check').fetchone()[0] != 'ok':
            raise ValueError('SQLite quick_check failed')
        if conn.execute('PRAGMA foreign_key_check').fetchone():
            raise ValueError('SQLite foreign_key_check failed')
        counts = {table: conn.execute(f'SELECT count(*) FROM {table} WHERE species_id=?', (species,)).fetchone()[0]
                  for table in ('genes', 'transcripts', 'proteins')}
        if not all(counts.values()):
            raise ValueError(f'Empty core records: {counts}')
        missing = conn.execute("SELECT count(*) FROM proteins p WHERE species_id=? AND NOT EXISTS (SELECT 1 FROM pfam_protein_status s WHERE s.species_id=p.species_id AND s.protein_id=p.protein_id AND s.status='completed')", (species,)).fetchone()[0]
        if missing:
            raise ValueError(f'Pfam proteins not completed: {missing}')
        return counts


def check_shared_resources(source, manifest):
    reference = {'dicot': 'arabidopsis_tair10_ensembl62', 'monocot': 'rice_rapdb_irgsp1_2026'}[manifest['taxon_group']]
    with closing(sqlite3.connect(f'file:{source}?mode=ro', uri=True)) as conn:
        if not conn.execute('SELECT 1 FROM reference_species WHERE reference_id=?', (reference,)).fetchone():
            raise ValueError(f'Shared reference must be built first: {reference}')
        if not conn.execute('SELECT 1 FROM go_terms LIMIT 1').fetchone():
            raise ValueError('Shared GO ontology must be built first')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', required=True)
    parser.add_argument('--output', required=True, help='New staging directory, not the live data directory')
    parser.add_argument('--source-db', default=str(ROOT / 'data/antarctic_plant.db'), help='Existing DB containing shared reference/GO resources')
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--replace-species', action='store_true', help='Rebuild this species in the staging copy only')
    parser.add_argument('--execute', action='store_true', help='Without this flag, print the plan only')
    parser.add_argument('--start-at', help='Explicitly resume an existing staging build at this stage')
    args = parser.parse_args()
    if args.threads < 1:
        parser.error('--threads must be positive')
    path = Path(args.manifest).resolve()
    output = Path(args.output).resolve()
    source = Path(args.source_db).resolve()
    manifest, plan = make_plan(path, output, args.threads, args.replace_species)
    check_inputs(manifest, path)
    if not source.is_file():
        raise FileNotFoundError('Prepare a source DB with shared reference and GO data first')
    if output == source.parent or output in source.parents:
        raise ValueError('Output must be separate from the source database directory')
    stages = [s['stage'] for s in plan]
    if args.start_at and args.start_at not in stages:
        parser.error(f'--start-at must be one of {stages}')
    for stage in plan:
        print(f"[{stage['stage']}] {shlex.join(stage['command'])}", flush=True)
    if not args.execute:
        print('Plan only. No files or databases changed.')
        return
    check_shared_resources(source, manifest)
    config = {'manifest_sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'source_db': str(source), 'plan': plan}
    state_path = output / 'build-state.json'
    if args.start_at:
        state = json.loads(state_path.read_text())
        if state['config'] != config:
            raise ValueError('Resume configuration differs from the original build')
        first = stages.index(args.start_at)
        if any(state['stages'].get(name) != 'completed' for name in stages[:first]):
            raise ValueError('Cannot skip a prerequisite that has not completed')
    else:
        if output.exists():
            raise FileExistsError('Use a new output directory, or --start-at to resume an inspected build')
        output.mkdir(parents=True)
        with closing(sqlite3.connect(f'file:{source}?mode=ro', uri=True)) as src, closing(sqlite3.connect(output / 'antarctic_plant.db')) as dst:
            src.backup(dst)
        state = {'config': config, 'stages': {}}
        first = 0
    env = dict(os.environ, ANTARCTIC_DB=str(output / 'antarctic_plant.db'), ANTARCTIC_RUNTIME_DIR=str(output / 'runtime'))
    def save():
        state['updated_at'] = dt.datetime.now(dt.timezone.utc).isoformat()
        temporary = state_path.with_suffix('.tmp')
        temporary.write_text(json.dumps(state, indent=2) + '\n')
        temporary.replace(state_path)
    state['status'] = 'running'
    state.pop('error', None)
    state.pop('verification', None)
    for stage in plan[first:]:
        state['stages'][stage['stage']] = 'pending'
    save()
    try:
        for stage in plan[first:]:
            name = stage['stage']
            state['stages'][name] = 'running'
            save()
            with (output / f'{name}.log').open('a') as log:
                subprocess.run(stage['command'], cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
            state['stages'][name] = 'completed'
            save()
            print(f'Completed: {name}', flush=True)
        state['verification'] = verify(output / 'antarctic_plant.db', manifest['species_id'])
        state['status'] = 'completed'
    except BaseException as exc:
        state['status'] = 'failed'
        state['error'] = str(exc)
        for name, status in state['stages'].items():
            if status == 'running':
                state['stages'][name] = 'failed'
        raise
    finally:
        save()


if __name__ == '__main__':
    main()
