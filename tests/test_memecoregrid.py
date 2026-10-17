#!/usr/bin/env python3
"""Unit tests for memecoregrid module."""

import pathlib
import sys
import unittest

# Ensure module lib is importable
CUR = pathlib.Path(__file__).resolve().parent
REPO_ROOT = CUR
while REPO_ROOT.parent != REPO_ROOT and not (REPO_ROOT / "modules").is_dir():
    REPO_ROOT = REPO_ROOT.parent

LIB_DIR = (
    REPO_ROOT
    / "modules"
    / "local"
    / "memecoregrid"
    / "resources"
    / "usr"
    / "lib"
)
if str(LIB_DIR) not in sys.path:
    sys.path.insert(0, str(LIB_DIR))

from memegrid.evaluator import (  # noqa: E402
    calculate_consensus_entropy,
    parse_meme_txt,
)
from memegrid.grid_runner import generate_grid_conditions  # noqa: E402
from memegrid.slicer import slice_window_at_radius  # noqa: E402


class TestMemeCoreGrid(unittest.TestCase):
    def test_window_slicing(self):
        # 300 bp test sequence: 150 'A's + 150 'C's
        seq = "A" * 150 + "C" * 150
        # Radius 30 should yield 30 'A's + 30 'C's = 60 bp
        sliced_30 = slice_window_at_radius(seq, radius=30, center_soma_len=150)
        self.assertEqual(len(sliced_30), 60)
        self.assertEqual(sliced_30, "A" * 30 + "C" * 30)

        # Radius 50 should yield 50 'A's + 50 'C's = 100 bp
        sliced_50 = slice_window_at_radius(seq, radius=50, center_soma_len=150)
        self.assertEqual(len(sliced_50), 100)
        self.assertEqual(sliced_50, "A" * 50 + "C" * 50)

    def test_consensus_entropy(self):
        # Homopolymer should have 0 entropy
        self.assertEqual(calculate_consensus_entropy("AAAAAAAAAA"), 0.0)
        # 4 bases balanced should have 2.0 bits
        self.assertAlmostEqual(calculate_consensus_entropy("ACGT"), 2.0, places=2)
        # Realistic motif TGTTCA
        ent = calculate_consensus_entropy("TGTTCA")
        self.assertGreater(ent, 1.0)

    def test_generate_grid_conditions(self):
        conds = generate_grid_conditions(radii=[30, 50])
        # 2 radii * 3 models * 2 strands * 3 width ranges = 36 conditions
        self.assertEqual(len(conds), 36)
        self.assertTrue(any(c["id"] == "R30_zoops_revcomp_w6_15" for c in conds))

    def test_parse_meme_txt(self):
        # Create a mock MEME output text
        mock_meme = """MEME version 4

MOTIF TGTTCA 1 width = 6 sites = 10 llr = 50 E-value = 1.2e-005

letter-probability matrix: alength= 4 w= 6 nsites= 10 E= 1.2e-005
  0.0 0.0 0.0 1.0
  0.0 0.0 1.0 0.0
  0.0 0.0 0.0 1.0
  0.0 0.0 0.0 1.0
  0.0 1.0 0.0 0.0
  1.0 0.0 0.0 0.0

Multilevel TGTTCA

sites sorted by position p-value
seq1  +  25  0.0001
seq2  +  28  0.0001
seq3  +  26  0.0001
seq4  +  27  0.0001
"""
        import tempfile

        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as tf:
            tf.write(mock_meme)
            tf_path = tf.name

        metrics = parse_meme_txt(tf_path, radius=30, total_sequences=4)
        self.assertIsNotNone(metrics)
        self.assertEqual(metrics["width"], 6)
        self.assertEqual(metrics["consensus"], "TGTTCA")
        # Pure discrete matrix should give maximum IC = 2.0 bits per position
        self.assertAlmostEqual(metrics["mean_ic"], 2.0, places=1)
        self.assertEqual(metrics["occupancy"], 1.0)
        self.assertGreater(metrics["composite_score"], 0.5)


if __name__ == "__main__":
    unittest.main()
