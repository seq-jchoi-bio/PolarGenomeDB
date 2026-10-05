import importlib.util
import random
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from backend.services import primer_jobs as p


def candidate(left=10, right=180, length=20):
    return {"rank": 1, "left_position": {"offset": left, "length": length},
            "right_position": {"offset": right, "length": length},
            "left_sequence": "ACGT" * 5, "right_sequence": "TGCA" * 5,
            "product_size": right - left + 1, "penalty": 1.0}


def target(strand="+"):
    segments = [
        {"number": 1, "start": 1, "end": 100, "genomic_start": 1001, "genomic_end": 1100},
        {"number": 2, "start": 101, "end": 200, "genomic_start": 2001, "genomic_end": 2100},
    ]
    if strand == "-":
        segments[0].update(genomic_start=2001, genomic_end=2100)
        segments[1].update(genomic_start=1001, genomic_end=1100)
    return {"target_id": "gene", "seqid": "chr1", "start": 1001, "end": 2100,
            "strand": strand, "genomic_region": "cdna", "junctions": [100], "exon_segments": segments}


class PrimerLogicTests(unittest.TestCase):
    def test_exon_only_excludes_cross_exon_product(self):
        c = candidate()
        self.assertTrue(p.candidate_product_spans_junction(c, [100]))
        self.assertEqual(p.filter_cdna_candidates([c], target(), "no_junction"), [])
        self.assertEqual(p.filter_cdna_candidates([c], target(), "allow_junction"), [c])

    def test_no_minimum_junction_overlap(self):
        for c in [candidate(left=93), candidate(left=99), candidate(right=100)]:
            self.assertEqual(p.filter_cdna_candidates([c], target(), "allow_junction"), [c])
            self.assertEqual(p.filter_cdna_candidates([c], target(), "no_junction"), [])

    def test_single_exon_is_excluded_from_multi_exon_mode(self):
        c = candidate(right=80)
        for t in [target(), {"junctions": []}]:
            self.assertEqual(p.filter_cdna_candidates([c], t, "no_junction"), [c])
            self.assertEqual(p.filter_cdna_candidates([c], t, "allow_junction"), [])

    def test_multi_exon_selection_has_no_single_exon_fallback(self):
        single, multi = candidate(right=80), candidate(left=60)
        multi["penalty"] = 10.0
        p.annotate_cdna_junctions([single, multi], target())
        selected = p.filter_cdna_candidates([single, multi], target(), "allow_junction")
        self.assertEqual(p.select_visible_candidates(selected), [multi])
        self.assertEqual(p.filter_cdna_candidates([single], target(), "allow_junction"), [])

    def test_expected_genomic_lengths_without_specificity(self):
        for strand in ["+", "-"]:
            multi, single, crossing = candidate(), candidate(right=80), candidate(left=99)
            p.annotate_cdna_junctions([multi, single, crossing], target(strand))
            self.assertEqual(multi["expected_genomic_product_size"], 1071)
            self.assertEqual(single["expected_genomic_product_size"], single["product_size"])
            self.assertIsNone(crossing["expected_genomic_product_size"])

    def test_plus_transcript_mapping(self):
        pair = p.expected_genome_pair(target(), candidate())
        self.assertEqual(pair["forward"], {"seqid": "chr1", "start": 1011, "end": 1030, "strand": "+"})
        self.assertEqual(pair["reverse"], {"seqid": "chr1", "start": 2062, "end": 2081, "strand": "-"})
        self.assertEqual(pair["span"], 1071)

    def test_minus_transcript_mapping(self):
        pair = p.expected_genome_pair(target("-"), candidate())
        self.assertEqual(pair["forward"], {"seqid": "chr1", "start": 2071, "end": 2090, "strand": "-"})
        self.assertEqual(pair["reverse"], {"seqid": "chr1", "start": 1020, "end": 1039, "strand": "+"})

    def test_gdna_reference_orientation_independent_of_gene_strand(self):
        t = {**target("-"), "genomic_region": "gene_body"}
        self.assertEqual(p.expected_genome_pair(t, candidate())["forward"]["strand"], "+")

    def test_junction_and_pasted_have_no_assigned_contiguous_pair(self):
        self.assertIsNone(p.expected_genome_pair(target(), candidate(left=93)))
        self.assertIsNone(p.expected_genome_pair({"seqid": None}, candidate()))

    def specificity(self, t, c, products, expected=None, limit=5000):
        response = {"engine": "bowtie2_exact_k2", "primer_counts": {c["left_sequence"]: 1, c["right_sequence"]: 1},
                    "products": [products], "expected_pairs": [expected], "screening_max_product_size": limit}
        request = p.PrimerRequest(species_id="Cq", gene_id="gene")
        with patch.object(p, "bowtie2_specificity", return_value=response):
            p.add_specificity(Path("/tmp"), request, t, [c])
        return c["specificity"]

    def test_single_hit_at_wrong_locus_is_not_on_target(self):
        c, t = candidate(), target()
        expected = p.expected_genome_pair(t, c)
        wrong = {"forward": {**expected["forward"], "seqid": "other"},
                 "reverse": {**expected["reverse"], "seqid": "other"}}
        s = self.specificity(t, c, [wrong], expected)
        self.assertEqual(s["on_target_product_count"], 0)
        self.assertEqual(s["status"], "off_target_detected")

    def test_expected_product_is_confirmed_by_both_coordinates(self):
        c, t = candidate(), target()
        expected = p.expected_genome_pair(t, c)
        self.assertEqual(self.specificity(t, c, [expected], expected)["status"], "target_confirmed")

    def test_zero_is_not_an_offtarget(self):
        c, t = candidate(), target()
        self.assertEqual(self.specificity(t, c, [], p.expected_genome_pair(t, c))["status"], "target_not_detected")

    def test_pasted_unique_does_not_imply_target_confirmation(self):
        s = self.specificity({"seqid": None}, candidate(), [{"forward": {}, "reverse": {}}])
        self.assertEqual(s["status"], "genome_matches_only")
        self.assertIsNone(s["on_target_product_count"])

    def test_out_of_range_target_is_not_reported_absent(self):
        c, t = candidate(), target()
        s = self.specificity(t, c, [], p.expected_genome_pair(t, c), limit=1000)
        self.assertEqual(s["status"], "target_out_of_range")

    def test_capped_counts_are_explicit(self):
        c, t = candidate(), target()
        expected = p.expected_genome_pair(t, c)
        other = {"forward": {**expected["forward"], "start": 999}, "reverse": expected["reverse"]}
        s = self.specificity(t, c, [expected, other], expected)
        self.assertEqual(s["pair_amplicon_count"], 2)
        self.assertTrue(s["pair_count_capped"])
        self.assertTrue(s["off_target_count_is_lower_bound"])

    def test_junction_genome_zero_is_separate_scope(self):
        c = candidate(left=93)
        s = self.specificity(target(), c, [])
        self.assertEqual(s["status"], "junction_genome_screen")
        self.assertIsNone(s["off_target_product_count"])

    def test_missing_bowtie_never_falls_back(self):
        with patch.object(p, "bowtie2_binary", return_value=None), patch.object(p, "bowtie2_index_prefix", return_value=Path("/tmp/test")):
            with self.assertRaisesRegex(ValueError, "requires Bowtie2"):
                p.require_bowtie2("Cq")
        with patch.object(p, "bowtie2_binary", return_value="bowtie2"), patch.object(p, "bowtie2_index_prefix", return_value=None):
            with self.assertRaisesRegex(ValueError, "complete Bowtie2 index"):
                p.require_bowtie2("Cq")

    def test_bowtie_error_propagates(self):
        with tempfile.TemporaryDirectory() as tmp, patch.object(p, "require_bowtie2", return_value=("bowtie2", Path("index"))), patch.object(p.subprocess, "run", side_effect=subprocess.CalledProcessError(1, "bowtie2")):
            with self.assertRaises(subprocess.CalledProcessError):
                p.bowtie2_specificity(Path(tmp), p.PrimerRequest(species_id="Cq", gene_id="gene"), target(), [candidate()])

    def test_screening_limit_expands_but_respects_policy(self):
        commands = []
        def run(command, **kwargs):
            commands.append(command)
            Path(command[command.index("-S") + 1]).write_text("@HD\tVN:1.0\n")
        with tempfile.TemporaryDirectory() as tmp, patch.object(p, "require_bowtie2", return_value=("bowtie2", Path("index"))), patch.object(p.subprocess, "run", side_effect=run):
            result = p.bowtie2_specificity(Path(tmp), p.PrimerRequest(species_id="Cq", gene_id="gene", product_size_range="1400-1600"), target(), [candidate()])
            self.assertEqual(result["screening_max_product_size"], 1600)
            self.assertEqual(commands[-1][commands[-1].index("-X") + 1], "1600")

    def test_sam_both_pair_orientations_and_exactness(self):
        for flags, positions in [((99, 147), (100, 180)), ((83, 163), (180, 100))]:
            with self.subTest(flags=flags), tempfile.TemporaryDirectory() as tmp:
                sam = Path(tmp) / "test.sam"
                sam.write_text("".join(f"pair0\t{flag}\tchr1\t{pos}\t255\t20M\t=\t{mate}\t0\t{'A'*20}\t*\tNM:i:0\n"
                                       for flag, pos, mate in zip(flags, positions, reversed(positions))))
                records = p.exact_sam_records(sam)
                products = p.sam_pair_products(records["pair0"])
                self.assertEqual(len(products), 1)
                self.assertEqual(products[0]["forward"]["start"], positions[0])
                sam.write_text(sam.read_text().replace("NM:i:0", "NM:i:1"))
                with self.assertRaisesRegex(ValueError, "non-exact"):
                    p.exact_sam_records(sam)


@unittest.skipUnless(importlib.util.find_spec("primer3"), "primer3-py required")
class CloningTests(unittest.TestCase):
    def test_old_out_of_range_example_is_not_silently_accepted(self):
        request = p.PrimerRequest(species_id="Cq", gene_id="Cq01g000020", template_type="cdna_cloning",
                                 primer_min_size=18, primer_max_size=24, primer_min_gc=55, primer_max_gc=65,
                                 primer_min_tm=57, primer_max_tm=63)
        # Preserve the real example's terminal cores without depending on the local DB.
        sequence = "ATGAGTACCTGTGATTACAAAAAT" + "ACGT" * 150 + "AGTGCATGGTGTGTTGATTGTTAA"
        self.assertEqual(p.full_length_cloning_pairs(sequence, request), [])

    def test_full_length_constraints_and_tail_direction(self):
        rng = random.Random(42)
        sequence = "".join(rng.choice("ACGT") for _ in range(400))
        request = p.PrimerRequest(species_id="Cq", gene_id="test", template_type="cdna_cloning",
                                 primer_min_size=18, primer_max_size=24,
                                 primer_min_tm=40, primer_max_tm=80, primer_min_gc=20, primer_max_gc=80,
                                 forward_tail="TAG", reverse_tail="GAT")
        pairs = p.full_length_cloning_pairs(sequence, request)
        self.assertTrue(pairs)
        p.apply_cloning_tails(pairs, request.forward_tail, request.reverse_tail)
        for pair in pairs:
            self.assertEqual(pair["left_position"]["offset"], 0)
            self.assertEqual(pair["right_position"]["offset"], len(sequence) - 1)
            self.assertEqual(pair["product_size"], len(sequence))
            self.assertEqual(pair["final_product_size"], len(sequence) + 6)
            self.assertTrue(pair["left_sequence"].startswith("TAG"))
            self.assertTrue(pair["right_sequence"].startswith("ATC"))
            self.assertLessEqual(max(pair["structure"].values()), 47)

    def test_cloning_defaults_and_explicit_overrides(self):
        base = dict(species_id="Cq", gene_id="test")
        clone = p.PrimerRequest(**base, template_type="cdna_cloning")
        self.assertEqual((clone.primer_min_size, clone.primer_max_size, clone.primer_max_tm,
                          clone.primer_min_gc, clone.primer_max_gc), (20, 40, 65, 30, 70))
        self.assertEqual(p.PrimerRequest(**base).primer_max_size, 24)
        self.assertEqual(p.PrimerRequest(**base, template_type="cdna_cloning", primer_min_gc=55).primer_min_gc, 55)
        with self.assertRaises(ValueError):
            p.PrimerRequest(**base, primer_max_size=40)

    def test_40_nt_and_tm_warning(self):
        from types import SimpleNamespace
        import primer3
        request = p.PrimerRequest(species_id="Cq", gene_id="test", template_type="cdna_cloning",
                                 primer_min_size=40, primer_max_size=40)
        sequence = "ACGT" * 20 + "AGCT" * 20
        with patch.object(primer3, "calc_tm", side_effect=[57.0, 64.0]), patch.object(primer3, "calc_hairpin", return_value=SimpleNamespace(tm=0)), patch.object(primer3, "calc_homodimer", return_value=SimpleNamespace(tm=0)), patch.object(primer3, "calc_heterodimer", return_value=SimpleNamespace(tm=0)):
            pairs = p.full_length_cloning_pairs(sequence, request)
        self.assertEqual(len(pairs), 1)
        self.assertEqual(len(pairs[0]["left_sequence"]), 40)
        self.assertEqual(pairs[0]["tm_difference"], 7)
        self.assertTrue(pairs[0]["warnings"])

    def test_empty_diagnostics_explain_both_ends(self):
        request = p.PrimerRequest(species_id="Cq", gene_id="test", template_type="cdna_cloning")
        diagnostics = {}
        self.assertEqual(p.full_length_cloning_pairs("AT" * 100, request, diagnostics), [])
        self.assertIn("Forward", diagnostics["empty_reason"])
        self.assertIn("Reverse", diagnostics["empty_reason"])
        self.assertIn("GC below minimum", diagnostics["empty_reason"])

    def test_structure_failure_explained(self):
        from types import SimpleNamespace
        import primer3
        request = p.PrimerRequest(species_id="Cq", gene_id="test", template_type="cdna_cloning")
        diagnostics = {}
        with patch.object(primer3, "calc_tm", return_value=60), patch.object(primer3, "calc_hairpin", return_value=SimpleNamespace(tm=0)), patch.object(primer3, "calc_homodimer", return_value=SimpleNamespace(tm=0)), patch.object(primer3, "calc_heterodimer", return_value=SimpleNamespace(tm=60)):
            self.assertEqual(p.full_length_cloning_pairs("ACGT" * 100, request, diagnostics), [])
        self.assertIn("Pair: heterodimer", diagnostics["empty_reason"])


if __name__ == "__main__":
    unittest.main()
