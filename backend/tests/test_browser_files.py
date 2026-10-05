import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from fastapi import FastAPI
from fastapi.testclient import TestClient
from backend.routes import browser_files as files


class BrowserFileTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.root_patch = patch.object(files, "ROOT", self.root)
        self.root_patch.start()
        self.index_patch = patch.object(files, "load_species_index", return_value={"species": [{"species_id": "Cq"}]})
        self.index_patch.start()
        self.addCleanup(self.temp.cleanup)
        self.addCleanup(self.root_patch.stop)
        self.addCleanup(self.index_patch.stop)
        (self.root / "data/jbrowse/Cq").mkdir(parents=True)
        (self.root / "genome.fa").write_bytes(b"0123456789abcdef")
        (self.root / "genome.fa.fai").write_text("index")
        (self.root / "private.db").write_text("secret")
        (self.root / "private.bam").write_text("secret")
        self.config = self.root / "data/jbrowse/Cq/config.json"
        self.config.write_text(json.dumps({"assemblies": [{"sequence": {"adapter": {
            "fastaLocation": {"uri": "/files/genome.fa"},
            "faiLocation": {"uri": "/files/genome.fa.fai"},
        }}}], "metadata": {"uri": "/files/private.bam"}}))
        app = FastAPI()
        app.include_router(files.router)
        self.client = TestClient(app)

    def test_private_files_and_traversal_are_denied(self):
        for path in ["private.db", "private.bam", "missing.fa", "%2e%2e/private.db", "data/jbrowse/Cq/../config.json"]:
            for method in [self.client.get, self.client.head]:
                self.assertEqual(method(f"/files/{path}").status_code, 404)

    def test_config_indexes_and_bounded_ranges(self):
        self.assertEqual(self.client.get("/files/data/jbrowse/Cq/config.json").status_code, 200)
        self.assertEqual(self.client.get("/files/genome.fa.fai").text, "index")
        self.assertEqual(self.client.get("/files/genome.fa").status_code, 403)
        response = self.client.get("/files/genome.fa", headers={"Range": "bytes=2-5", "If-Range": '"invalid"'})
        self.assertEqual(response.status_code, 206)
        self.assertEqual(response.content, b"2345")
        self.assertEqual(response.headers["content-range"], "bytes 2-5/16")
        self.assertEqual(self.client.head("/files/genome.fa").headers["content-length"], "16")

    def test_range_limits(self):
        for value in ["bytes=0-", "bytes=-4", "bytes=0-2,4-8", "bytes=9-2", "bytes=30-40", f"bytes=0-{files.MAX_RANGE_BYTES}"]:
            self.assertEqual(self.client.get("/files/genome.fa", headers={"Range": value}).status_code, 416)

    def test_external_symlink_is_denied(self):
        with tempfile.TemporaryDirectory() as external:
            outside = Path(external) / "outside.fa"
            outside.write_text("secret")
            (self.root / "genome.fa").unlink()
            (self.root / "genome.fa").symlink_to(outside)
            self.assertEqual(self.client.get("/files/genome.fa", headers={"Range": "bytes=0-3"}).status_code, 404)
