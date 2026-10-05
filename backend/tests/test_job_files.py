import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from backend.services.job_files import write_json


class JobFileTests(unittest.TestCase):
    def test_readers_see_previous_complete_record_until_replace(self):
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / 'status.json'
            write_json(target, {'status': 'queued'})
            original = Path.replace

            def inspect(temporary, destination):
                self.assertEqual(json.loads(target.read_text()), {'status': 'queued'})
                self.assertEqual(json.loads(temporary.read_text()), {'status': 'running'})
                self.assertEqual(temporary.parent, target.parent)
                return original(temporary, destination)

            with patch.object(Path, 'replace', inspect):
                write_json(target, {'status': 'running'})
            self.assertEqual(json.loads(target.read_text()), {'status': 'running'})
            self.assertEqual(list(target.parent.iterdir()), [target])

    def test_failed_replace_preserves_previous_record_and_cleans_temporary(self):
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / 'status.json'
            write_json(target, {'status': 'queued'})
            with patch.object(Path, 'replace', side_effect=OSError('failed')):
                with self.assertRaises(OSError):
                    write_json(target, {'status': 'running'})
            self.assertEqual(json.loads(target.read_text()), {'status': 'queued'})
            self.assertEqual(list(target.parent.iterdir()), [target])

    def test_invalid_payload_preserves_previous_record(self):
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / 'status.json'
            write_json(target, {'status': 'queued'})
            with self.assertRaises(TypeError):
                write_json(target, {'value': object()})
            self.assertEqual(json.loads(target.read_text()), {'status': 'queued'})
            self.assertEqual(list(target.parent.iterdir()), [target])
