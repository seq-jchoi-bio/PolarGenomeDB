import csv
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

from backend.services import similar_export as export


class SimilarExportTests(unittest.TestCase):
    def request(self):
        return {"species_id": "Da", "reference_id": "rice", "group": {
            "query_protein": {"protein_id": "query", "gene_id": "gene", "transcript_id": "tx", "length": 100},
            "candidates": [{"rank": 1, "protein_id": "p0", "gene_id": "g0", "transcript_id": "t0", "pident": 81.123456}],
        }, "raw_hits": [{"qseqid": "query", "sseqid": f"p{i}", "length": 90, "pident": 81.123456} for i in range(7)]}

    def test_all_raw_targets_are_compared_without_display_limit(self):
        compare = lambda *args: {"available": True, "shared_domains": [] if args[-1] == "p6" else [
            {"pfam_id": "PF001", "identity": 77.77, "alignment_error": None}]}
        with patch.object(export, "connect", return_value=MagicMock()), patch.object(export, "rows", return_value=[]), patch.object(export, "pairwise_domain_payload", side_effect=compare) as mocked:
            result = list(export.export_records(self.request()))
        self.assertEqual(mocked.call_count, 7)
        self.assertEqual(len(result), 15)
        self.assertEqual(result[0]["pident"], 81.123456)
        self.assertEqual(result[-1]["pairwise_status"], "no_shared_domains")

    def test_completed_tsv_has_header_and_preserves_tabs(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(export, "JOB_DIR", Path(directory)):
            job_id = "00000000-0000-0000-0000-000000000001"
            path = export.export_path(job_id)
            path.mkdir(parents=True)
            export.write_json(path / "status.json", {"species_id": "Da"})
            export.write_json(path / "request.json", {})
            with patch.object(export, "export_records", return_value=iter([{"record_type": "pairwise_domain", "domain_name": "A\tB"}])):
                export.execute_export(job_id)
            self.assertEqual(export.get_export("Da", job_id)["status"], "completed")
            status = export.read_json(path / "status.json")
            self.assertLessEqual(status["started_at"], status["finished_at"])
            self.assertEqual(status["updated_at"], status["finished_at"])
            with (path / "results.tsv").open() as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                self.assertEqual(reader.fieldnames, export.EXPORT_COLUMNS)
                self.assertEqual(next(reader)["domain_name"], "A\tB")
            with self.assertRaises(FileNotFoundError):
                export.get_export("Cq", job_id)

    def test_alignment_failure_is_not_silently_exported(self):
        with patch.object(export, "connect", return_value=MagicMock()), patch.object(export, "rows", return_value=[]), patch.object(export, "pairwise_domain_payload", return_value={"available": True, "shared_domains": [{"alignment_error": "timeout"}]}):
            with self.assertRaisesRegex(RuntimeError, "timeout"):
                list(export.export_records(self.request()))

    def test_failed_export_has_terminal_timestamp(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(export, "JOB_DIR", Path(directory)):
            job_id = "00000000-0000-0000-0000-000000000002"
            path = export.export_path(job_id)
            path.mkdir(parents=True)
            export.write_json(path / "status.json", {})
            export.write_json(path / "request.json", {})
            with patch.object(export, "export_records", side_effect=ValueError("test failure")):
                export.execute_export(job_id)
            status = export.read_json(path / "status.json")
            self.assertEqual(status["status"], "failed")
            self.assertEqual(status["failed_at"], status["finished_at"])
            self.assertEqual(status["updated_at"], status["finished_at"])


if __name__ == "__main__":
    unittest.main()
