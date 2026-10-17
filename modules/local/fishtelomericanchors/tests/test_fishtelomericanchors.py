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
        obs1 = {
            "read_id": "read_1",
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

        clusters = cluster_anchors([obs1, obs2], min_identity=0.97)
        self.assertEqual(len(clusters), 1)
        self.assertEqual(clusters[0]["read_count"], 2)
        self.assertEqual(clusters[0]["consensus_seq"], self.soma_150)


if __name__ == "__main__":
    unittest.main()
