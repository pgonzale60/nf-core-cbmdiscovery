#!/usr/bin/env python3
"""Unit tests for fishtelomericanchors module."""

import pathlib
import sys
import unittest

# Ensure module lib is importable
CUR = pathlib.Path(__file__).resolve().parent
REPO_ROOT = CUR
while REPO_ROOT.parent != REPO_ROOT and not (REPO_ROOT / "modules").is_dir():
    REPO_ROOT = REPO_ROOT.parent

LIB_DIR = REPO_ROOT / "modules" / "local" / "fishtelomericanchors" / "resources" / "usr" / "lib"
if str(LIB_DIR) not in sys.path:
    sys.path.insert(0, str(LIB_DIR))

from teloanchor.cluster import (  # noqa: E402
    cluster_anchors,
    sequence_identity_anchored_3p,
)
from teloanchor.detector import (  # noqa: E402
    find_terminal_telomeric_tract,
    get_cyclic_permutations,
    revcomp,
)
from teloanchor.orient import (  # noqa: E402
    get_cyclic_prefixes,
    normalize_and_trim_flank,
)


class TestTeloAnchor(unittest.TestCase):
    def setUp(self):
        self.repeat = "TTAGGC"
        # 150 bp distinct somatic sequence
        self.soma_150 = (
            "GATTACACATTAGACCGGATACACAGATAGAGACAGATAGACAGATAGAC"
            "ACGTACGTACGTTTGCAATCGATCGATCGATCGATCGATCGATCGATCGA"
            "TGCTAGCTAGCTAGCTAAGGTCCAAATTGGCCAAGGTTAAGGCCTTGGAA"
        )
        self.assertTrue(len(self.soma_150) == 150)

    def test_revcomp(self):
        self.assertEqual(revcomp("TTAGGC"), "GCCTAA")
        self.assertEqual(revcomp("AATTCG"), "CGAATT")

    def test_cyclic_prefixes(self):
        prefixes = get_cyclic_prefixes("TTAGGC")
        self.assertEqual(prefixes, ["TTAGG", "TTAG", "TTA", "TT", "T"])

    def test_cyclic_permutations(self):
        perms = get_cyclic_permutations("TTAGGC")
        expected = sorted([
            "TTAGGC",
            "TAGGCT",
            "AGGCTT",
            "GGCTTA",
            "GCTTAG",
            "CTTAGG",
        ])
        self.assertEqual(perms, expected)

    def test_hifi_forward_read_orientation_and_trimming(self):
        # Flank + microhomology 'TTAG' + 12 telomeric repeats
        telo_array = "TTAGGC" * 12
        full_read = self.soma_150 + "TTAG" + telo_array + "A" * 20

        tract = find_terminal_telomeric_tract(
            full_read, repeat=self.repeat, platform="hifi"
        )
        self.assertIsNotNone(tract)
        self.assertEqual(tract["terminal_end"], "3p")
        self.assertEqual(tract["tract_orientation"], "fwd")

        trimmed = normalize_and_trim_flank(
            full_read, tract, repeat=self.repeat, max_flank_len=150
        )
        self.assertIsNotNone(trimmed)
        self.assertEqual(trimmed["trimmed_prefix"], "TTAG")
        self.assertEqual(trimmed["microhomology_trimmed_len"], 4)
        # Position 0 should be the exact 3' end of self.soma_150
        self.assertEqual(trimmed["anchor_seq"], self.soma_150)

    def test_hifi_reverse_read_orientation_and_trimming(self):
        # Bottom strand: 5' tip + 12 revcomp repeats (GCCTAA) + revcomp('TTAG')=CTAA + revcomp(self.soma_150)
        telo_rev_array = "GCCTAA" * 12
        full_read = (
            "T" * 20 + telo_rev_array + "CTAA" + revcomp(self.soma_150)
        )

        tract = find_terminal_telomeric_tract(
            full_read, repeat=self.repeat, platform="hifi"
        )
        self.assertIsNotNone(tract)
        self.assertEqual(tract["terminal_end"], "5p")
        self.assertEqual(tract["tract_orientation"], "rev")

        trimmed = normalize_and_trim_flank(
            full_read, tract, repeat=self.repeat, max_flank_len=150
        )
        self.assertIsNotNone(trimmed)
        self.assertEqual(trimmed["trimmed_prefix"], "TTAG")
        self.assertEqual(trimmed["anchor_seq"], self.soma_150)

    def test_illumina_read_detection(self):
        # Illumina: 75 bp somatic flank + 'TT' + 5 repeats
        flank_75 = self.soma_150[-75:]
        read_seq = flank_75 + "TT" + ("TTAGGC" * 5)

        tract = find_terminal_telomeric_tract(
            read_seq, repeat=self.repeat, platform="illumina"
        )
        self.assertIsNotNone(tract)
        self.assertEqual(tract["terminal_end"], "3p")

        trimmed = normalize_and_trim_flank(
            read_seq, tract, repeat=self.repeat, max_flank_len=150
        )
        self.assertIsNotNone(trimmed)
        self.assertEqual(trimmed["trimmed_prefix"], "TT")
        self.assertEqual(trimmed["anchor_seq"], flank_75)

    def test_clustering_and_consensus(self):
        # Create observations with small sequencing errors
        obs1_a = {
            "read_id": "read_1a",
            "anchor_seq": self.soma_150,
            "trimmed_prefix": "TTA",
        }
        obs1_b = {
            "read_id": "read_1b",
            "anchor_seq": self.soma_150,
            "trimmed_prefix": "TTA",
        }
        # 1 mismatch in 150 bp
        mutated_seq = (
            self.soma_150[:-10]
            + ("A" if self.soma_150[-10] != "A" else "C")
            + self.soma_150[-9:]
        )
        obs2 = {
            "read_id": "read_2",
            "anchor_seq": mutated_seq,
            "trimmed_prefix": "TTA",
        }

        ident = sequence_identity_anchored_3p(
            self.soma_150, mutated_seq
        )
        self.assertGreaterEqual(ident, 0.99)

        clusters = cluster_anchors([obs1_a, obs1_b, obs2], min_identity=0.97)
        self.assertEqual(len(clusters), 1)
        self.assertEqual(clusters[0]["read_count"], 3)
        self.assertEqual(clusters[0]["consensus_seq"], self.soma_150)

    def test_continuous_phase_microhomology_trimming(self):
        # Flank ending in CG + 'TT' prefix + 10 repeats
        soma = self.soma_150[:-4] + "CG"
        full_read = soma + "TT" + ("TTAGGC" * 10) + "A" * 10

        tract = find_terminal_telomeric_tract(
            full_read, repeat=self.repeat, platform="hifi"
        )
        self.assertIsNotNone(tract)
        trimmed = normalize_and_trim_flank(
            full_read, tract, repeat=self.repeat, max_flank_len=150
        )
        self.assertIsNotNone(trimmed)
        # Position 0 must be trimmed back to 'CG', excluding 'TT'
        self.assertEqual(trimmed["anchor_seq"], soma)
        self.assertEqual(trimmed["microhomology_trimmed_len"], 2)

    def test_inward_pointing_telomere_rejected(self):
        # 5p forward array: 'TTAGGC' * 12 + soma_150
        read_5p_fwd = ("TTAGGC" * 12) + self.soma_150
        # 3p reverse array: soma_150 + 'GCCTAA' * 12
        read_3p_rev = self.soma_150 + ("GCCTAA" * 12)

        tract_5p = find_terminal_telomeric_tract(
            read_5p_fwd, repeat=self.repeat, platform="hifi"
        )
        if tract_5p is not None:
            flank_5p = normalize_and_trim_flank(
                read_5p_fwd, tract_5p, repeat=self.repeat
            )
            self.assertIsNone(flank_5p)

        tract_3p = find_terminal_telomeric_tract(
            read_3p_rev, repeat=self.repeat, platform="hifi"
        )
        if tract_3p is not None:
            flank_3p = normalize_and_trim_flank(
                read_3p_rev, tract_3p, repeat=self.repeat
            )
            self.assertIsNone(flank_3p)

    def test_phase_trimming_across_arbitrary_somatic_endings(self):
        soma_base = self.soma_150[:100]
        # Case A: soma_base + "AA" + "TT" (2-bp prefix) + ("TTAGGC" * 10)
        read_a = soma_base + "AA" + "TT" + (self.repeat * 10)
        tract_a = find_terminal_telomeric_tract(
            read_a, repeat=self.repeat, platform="hifi"
        )
        self.assertIsNotNone(tract_a)
        trimmed_a = normalize_and_trim_flank(
            read_a, tract_a, repeat=self.repeat
        )
        self.assertIsNotNone(trimmed_a)
        self.assertTrue(trimmed_a["anchor_seq"].endswith("AA"))
        self.assertEqual(len(trimmed_a["anchor_seq"]), len(soma_base) + 2)
        self.assertEqual(trimmed_a["microhomology_trimmed_len"], 2)

        # Case B: soma_base + "CC" + "TTAG" (4-bp prefix) + ("TTAGGC" * 10)
        read_b = soma_base + "CC" + "TTAG" + (self.repeat * 10)
        tract_b = find_terminal_telomeric_tract(
            read_b, repeat=self.repeat, platform="hifi"
        )
        self.assertIsNotNone(tract_b)
        trimmed_b = normalize_and_trim_flank(
            read_b, tract_b, repeat=self.repeat
        )
        self.assertIsNotNone(trimmed_b)
        self.assertTrue(trimmed_b["anchor_seq"].endswith("CC"))
        self.assertEqual(trimmed_b["microhomology_trimmed_len"], 4)

        # Case C (5p reverse complement): ("GCCTAA" * 10) + "GCCT" + revcomp(soma_base + "CC" + "TTAG")
        rev_rep = revcomp(self.repeat)
        read_c = (rev_rep * 10) + "GCCT" + revcomp(soma_base + "CC" + "TTAG")
        tract_c = find_terminal_telomeric_tract(
            read_c, repeat=self.repeat, platform="hifi"
        )
        self.assertIsNotNone(tract_c)
        trimmed_c = normalize_and_trim_flank(
            read_c, tract_c, repeat=self.repeat
        )
        self.assertIsNotNone(trimmed_c)
        self.assertEqual(trimmed_c["anchor_seq"], soma_base + "CC")
        self.assertEqual(trimmed_c["microhomology_trimmed_len"], 4)

    def test_fast_prefilter_inside_detector(self):
        import time

        # 1000 non-telomeric reads + 5 telomeric reads
        non_telo = [self.soma_150] * 1000
        telo = [self.soma_150 + (self.repeat * 12)] * 5
        all_reads = non_telo + telo

        t0 = time.perf_counter()
        fished = 0
        for r in all_reads:
            if find_terminal_telomeric_tract(
                r, repeat=self.repeat, platform="hifi"
            ):
                fished += 1
        elapsed = time.perf_counter() - t0
        self.assertEqual(fished, 5)
        self.assertLess(elapsed, 0.05)


if __name__ == "__main__":
    unittest.main()
