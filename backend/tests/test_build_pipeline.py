import json
from contextlib import closing
from pathlib import Path
import sqlite3
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import yaml

from backend.services.bowtie_indexes import index_format, PARTS
from scripts.build_all import make_plan, check_inputs
from scripts.register_bowtie2_index import register, inspect_lengths


class BuildPipelineTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.prefix = self.root / 'genome'
        self.manifest = self.root / 'species.yaml'
        self.fasta = self.root / 'genome.fa'
        self.fasta.write_text('>chr1\nACGT\n')
        self.data = {'species_id': 'Test', 'taxon_group': 'dicot',
                     'genome': {'fasta': 'genome.fa'}, 'annotation': {'file': 'genome.fa'},
                     'protein': {'fasta': 'genome.fa'}, 'primer': {'bowtie2_index_prefix': 'genome'}}
        self.write_manifest()
        self.db = self.root / 'db.sqlite'
        with closing(sqlite3.connect(self.db)) as conn, conn:
            conn.executescript('CREATE TABLE species(species_id TEXT PRIMARY KEY); INSERT INTO species VALUES ("Test"); CREATE TABLE chromosomes(species_id TEXT,seqid TEXT,length INTEGER); INSERT INTO chromosomes VALUES ("Test","chr1",4); CREATE TABLE genome_assets(species_id TEXT,asset_type TEXT,path TEXT);')
            conn.execute('INSERT INTO genome_assets VALUES (?,?,?)', ('Test', 'fasta', str(self.fasta)))

    def tearDown(self):
        self.tmp.cleanup()

    def write_manifest(self):
        self.manifest.write_text(yaml.safe_dump(self.data))

    def files(self, extension='bt2'):
        for part in PARTS:
            Path(f'{self.prefix}.{part}.{extension}').write_bytes(b'index')

    def test_formats_and_missing_empty_mixed(self):
        with self.assertRaises(ValueError):
            index_format(self.prefix)
        self.files('bt2l')
        self.assertEqual(index_format(self.prefix), 'bt2l')
        Path(f'{self.prefix}.1.bt2l').write_bytes(b'')
        with self.assertRaises(ValueError):
            index_format(self.prefix)
        self.files('bt2l')
        self.files('bt2')
        with self.assertRaises(ValueError):
            index_format(self.prefix)

    def test_register_readonly_mismatch_and_runtime_lookup(self):
        self.files()
        def tool(cmd, **kwargs):
            return subprocess.CompletedProcess(cmd, 0, stdout='version1\n' if '--version' in cmd else 'Sequence-1\tchr1\t4\n')
        with patch('scripts.register_bowtie2_index.shutil.which', return_value='inspect'), patch('scripts.register_bowtie2_index.subprocess.run', side_effect=tool):
            before = self.db.read_bytes()
            self.assertEqual(register(self.manifest, self.db, True)['status'], 'ready')
            self.assertEqual(before, self.db.read_bytes())
            register(self.manifest, self.db)
            from backend.services.primer_jobs import bowtie2_index_prefix
            with patch('backend.services.primer_jobs.connect', side_effect=lambda: closing(sqlite3.connect(self.db))):
                self.assertEqual(bowtie2_index_prefix('Test'), self.prefix.resolve())
            with closing(sqlite3.connect(self.db)) as conn, conn:
                conn.execute('UPDATE chromosomes SET length=5')
            with self.assertRaisesRegex(ValueError, 'names/lengths'):
                register(self.manifest, self.db)

    def test_no_prefix_disables_even_legacy_lookup(self):
        self.data.pop('primer')
        self.write_manifest()
        self.assertEqual(register(self.manifest, self.db)['status'], 'disabled')
        from backend.services.primer_jobs import bowtie2_index_prefix
        with patch('backend.services.primer_jobs.connect', side_effect=lambda: closing(sqlite3.connect(self.db))):
            self.assertIsNone(bowtie2_index_prefix('Test'))

    def test_duplicate_inspect_names_rejected(self):
        with self.assertRaises(ValueError):
            inspect_lengths('Sequence-1\tchr1\t4\nSequence-2\tchr1\t4')

    def test_plan_optional_expression_and_global_runtime(self):
        _, plan = make_plan(self.manifest, self.root / 'out', 2)
        self.assertNotIn('expression', [s['stage'] for s in plan])
        self.assertNotIn('--species', next(s['command'] for s in plan if s['stage'] == 'runtime'))
        self.data['expression'] = {'sample_metadata': 'meta.csv', 'count_matrix': 'count.txt'}
        self.write_manifest()
        _, plan = make_plan(self.manifest, self.root / 'out', 2)
        command = next(s['command'] for s in plan if s['stage'] == 'expression')
        self.assertIn('--load-counts', command)
        self.assertNotIn('--build-bigwig', command)
        self.data['expression']['bam_dir'] = 'bam'
        self.write_manifest()
        _, plan = make_plan(self.manifest, self.root / 'out', 2)
        self.assertIn('--build-bigwig', next(s['command'] for s in plan if s['stage'] == 'expression'))

    def test_missing_input_is_rejected(self):
        self.data['protein']['fasta'] = 'missing.fa'
        with self.assertRaises(FileNotFoundError):
            check_inputs(self.data, self.manifest)

    def test_orchestrator_stops_on_failure_and_can_resume(self):
        from scripts.build_all import main
        out = self.root / 'stage'
        plan = [{'stage': name, 'command': ['test-tool', name]} for name in ('species', 'search', 'runtime')]
        argv = ['build_all', '--manifest', str(self.manifest), '--output', str(out), '--source-db', str(self.db), '--execute']
        before = self.db.read_bytes()
        with patch('scripts.build_all.make_plan', return_value=(self.data, plan)), patch('scripts.build_all.check_shared_resources'), patch('sys.argv', argv):
            with patch('scripts.build_all.subprocess.run', side_effect=[None, subprocess.CalledProcessError(1, 'search')]) as run:
                with self.assertRaises(subprocess.CalledProcessError):
                    main()
                self.assertEqual(run.call_count, 2)
            state = json.loads((out / 'build-state.json').read_text())
            self.assertEqual(state['status'], 'failed')
            self.assertEqual(state['stages']['runtime'], 'pending')
            with patch('sys.argv', argv + ['--start-at', 'search']), patch('scripts.build_all.subprocess.run') as run, patch('scripts.build_all.verify', return_value={'genes': 1}):
                main()
                self.assertEqual(run.call_count, 2)
                self.assertEqual(run.call_args.kwargs['env']['ANTARCTIC_DB'], str(out.resolve() / 'antarctic_plant.db'))
            self.assertEqual(json.loads((out / 'build-state.json').read_text())['status'], 'completed')
        self.assertEqual(before, self.db.read_bytes())

    def test_checksum_mismatch_does_not_register(self):
        self.files()
        self.data['primer']['genome_sha256'] = '0' * 64
        self.write_manifest()
        with patch('scripts.register_bowtie2_index.shutil.which', return_value='inspect'), patch('scripts.register_bowtie2_index.subprocess.run', side_effect=[subprocess.CompletedProcess([], 0, 'Sequence-1\tchr1\t4\n'), subprocess.CompletedProcess([], 0, 'version\n')]):
            with self.assertRaisesRegex(ValueError, 'SHA256'):
                register(self.manifest, self.db)


if __name__ == '__main__':
    unittest.main()
