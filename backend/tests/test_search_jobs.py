import unittest
import subprocess
from pathlib import Path
from unittest.mock import patch

from fastapi import FastAPI
from fastapi.testclient import TestClient

from backend.services import search_jobs as search


class SearchQueryLimitTests(unittest.TestCase):
    @staticmethod
    def fasta_records(count: int) -> str:
        return "\n".join(f">query_{index}\nMPEPTIDE" for index in range(1, count + 1))

    def test_plain_sequence_is_one_query(self):
        self.assertEqual(search.query_sequence_count("MPEPTIDE"), 1)

    def test_blast_accepts_ten_fasta_queries(self):
        query = self.fasta_records(10)
        self.assertEqual(search.query_sequence_count(query), 10)
        self.assertTrue(search.normalize_query(query, max_sequences=10, program="blastp").endswith("\n"))

    def test_blast_rejects_eleventh_fasta_query(self):
        with self.assertRaisesRegex(ValueError, "max 10"):
            search.normalize_query(self.fasta_records(11), max_sequences=10, program="blastp")

    def test_diamond_accepts_one_hundred_fasta_queries(self):
        query = self.fasta_records(100)
        self.assertEqual(search.query_sequence_count(query), 100)
        self.assertTrue(search.normalize_query(query, max_sequences=100, program="diamond").endswith("\n"))

    def test_diamond_rejects_one_hundred_and_first_fasta_query(self):
        with self.assertRaisesRegex(ValueError, "max 100"):
            search.normalize_query(self.fasta_records(101), max_sequences=100, program="diamond")

    def test_dna_iupac_and_whitespace(self):
        self.assertEqual(search.normalize_query("acgt rysw\n kmbdhvn", max_sequences=10, program="blastn"), ">query\nACGTRYSWKMBDHVN\n")

    def test_protein_ambiguity_and_stop(self):
        for program in ("blastp", "diamond"):
            self.assertEqual(search.normalize_query(">protein description\nmpeptideXbzjuo*", max_sequences=100, program=program), ">protein description\nMPEPTIDEXBZJUO*\n")

    def test_reject_invalid_characters_by_program(self):
        for program, query in (("blastn", "not_a_sequence_123!?"), ("blastn", "ACGTPE"), ("blastp", "MPEP-1"), ("diamond", "MPEP?")):
            with self.subTest(program=program, query=query), self.assertRaisesRegex(ValueError, "must contain only"):
                search.normalize_query(query, max_sequences=100, program=program)

    def test_reject_malformed_fasta(self):
        for query, message in ((">\nACGT", "needs a name"), (">empty\n>next\nACGT", "no sequence"), (">a first\nACGT\n>a second\nTGCA", "unique ID"), ("ACGT\n>a\nTGCA", "begin with")):
            with self.subTest(query=query), self.assertRaisesRegex(ValueError, message):
                search.normalize_query(query, max_sequences=10, program="blastn")

    def test_plain_and_fasta_have_equivalent_sequences(self):
        self.assertEqual(search.normalize_query("acgt\nacgt", max_sequences=10, program="blastn"), search.normalize_query(">query\nacgt\nacgt", max_sequences=10, program="blastn"))

    def test_invalid_request_never_creates_or_submits_job(self):
        from backend.routes.search import router, job_manager
        app = FastAPI()
        app.include_router(router)
        with TestClient(app) as client, patch.object(search, "asset_for_program") as asset, patch.object(search, "job_path") as job_path, patch.object(job_manager, "submit") as submit:
            for endpoint, program in (("blast", "blastn"), ("blast", "blastp"), ("diamond", "diamond")):
                response = client.post(f"/search/{endpoint}", json={"species_id": "Cq", "program": program, "query": "123!?"})
                self.assertEqual(response.status_code, 400)
                self.assertIn("must contain only", response.json()["detail"])
            asset.assert_not_called()
            job_path.assert_not_called()
            submit.assert_not_called()

    def test_failure_messages_do_not_expose_commands(self):
        error = subprocess.CalledProcessError(1, ["/private/tool", "/private/query.fa"])
        self.assertNotIn("/private", search.search_failure_message(error))
        self.assertIn("time limit", search.search_failure_message(subprocess.TimeoutExpired("blastn", 10)))
        self.assertIn("unavailable", search.search_failure_message(FileNotFoundError("/private/tool")))

    def test_result_limits_are_enforced_by_program(self):
        self.assertEqual(search.max_hits_for_program("blastn"), 10)
        self.assertEqual(search.max_hits_for_program("blastp"), 10)
        self.assertEqual(search.max_hits_for_program("diamond"), 1)

        asset = {"index_prefix": "data/test_index"}
        blast_command = search.command_for_request(
            search.SearchRequest(species_id="Cq", program="blastp", query="MPEPTIDE", max_hits=100),
            asset,
            Path("query.fa"),
            Path("results.tsv"),
        )
        diamond_command = search.command_for_request(
            search.SearchRequest(species_id="Cq", program="diamond", query="MPEPTIDE", max_hits=100),
            asset,
            Path("query.fa"),
            Path("results.tsv"),
        )
        self.assertEqual(blast_command[blast_command.index("-max_target_seqs") + 1], "10")
        self.assertEqual(diamond_command[diamond_command.index("--max-target-seqs") + 1], "1")

    def test_blast_raw_export_uses_blast_defaults_without_target_cap(self):
        command = search.command_for_request(
            search.SearchRequest(species_id="Cq", program="blastn", query="ATGC"),
            {"index_prefix": "data/test_index"},
            Path("query.fa"),
            Path("raw_results.tsv"),
            raw_export=True,
        )
        self.assertNotIn("-max_target_seqs", command)
        self.assertNotIn("-evalue", command)


if __name__ == "__main__":
    unittest.main()
