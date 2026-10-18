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

    def test_continuous_phase_microhomology_trimming(self):
        # soma_150 ending in GCTT (shares 2 bp 'TT' microhomology with TTAGGC)
        soma = self.soma_150[:-4] + "GCTT"
        # Telomerase adds in-phase: AGGC + TTAGGC * 10
        telomerase_addition = "AGGC" + ("TTAGGC" * 10)
        full_read = soma + telomerase_addition + "A" * 10

        tract = find_terminal_telomeric_tract(
            full_read, repeat=self.repeat, platform="hifi"
        )
        self.assertIsNotNone(tract)
        trimmed = normalize_and_trim_flank(
            full_read, tract, repeat=self.repeat, max_flank_len=150
        )
        self.assertIsNotNone(trimmed)
        # Position 0 must be trimmed back to 'GC', excluding 'TT'
        self.assertEqual(trimmed["anchor_seq"], soma[:-2])
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

    def test_fast_prefilter_short_circuit(self):
        import time
        from teloanchor.detector import find_terminal_telomeric_tract

        # 1000 non-telomeric reads
        non_telo = [self.soma_150] * 1000
        # 5 telomeric reads
        telo = [self.soma_150 + ("TTAGGC" * 12)] * 5
        all_reads = non_telo + telo

        t0 = time.perf_counter()
        fished = 0
        rep3 = self.repeat * 3
        rc_rep3 = revcomp(self.repeat) * 3
        for r in all_reads:
            if rep3 not in r and rc_rep3 not in r:
                continue
            if find_terminal_telomeric_tract(r, repeat=self.repeat, platform="hifi"):
                fished += 1
        elapsed = time.perf_counter() - t0
        self.assertEqual(fished, 5)
        self.assertLess(elapsed, 0.05)


if __name__ == "__main__":
    unittest.main()
