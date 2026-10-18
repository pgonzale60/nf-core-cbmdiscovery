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

    def test_core_30bp_beats_short_6bp_motif(self):
        import tempfile

        # Condition 1: 6 bp, mean_ic ~ 1.85, evalue 1e-6
        # matrix of 6 rows
        m6 = """MEME version 4
MOTIF SHORT 1 width = 6 sites = 20 llr = 50 E-value = 1e-006
letter-probability matrix: alength= 4 w= 6 nsites= 20 E= 1e-006
  0.9 0.05 0.05 0.0
  0.9 0.05 0.05 0.0
  0.9 0.05 0.05 0.0
  0.0 0.9 0.05 0.05
  0.0 0.05 0.9 0.05
  0.05 0.05 0.0 0.9
Multilevel AAACGT
sites sorted by position p-value
s1 + 25 1e-5
s2 + 25 1e-5
"""
        # Condition 2: 30 bp, mean_ic ~ 1.55, evalue 1e-40
        row_30 = "\n".join(["  0.8 0.1 0.1 0.0" if i % 2 == 0 else "  0.0 0.8 0.1 0.1" for i in range(30)])
        m30 = f"""MEME version 4
MOTIF LONGCORE 1 width = 30 sites = 20 llr = 200 E-value = 1e-040
letter-probability matrix: alength= 4 w= 30 nsites= 20 E= 1e-040
{row_30}
Multilevel {"AC" * 15}
sites sorted by position p-value
s1 + 25 1e-5
s2 + 25 1e-5
"""
        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f6:
            f6.write(m6)
            p6 = f6.name
        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f30:
            f30.write(m30)
            p30 = f30.name

        met6 = parse_meme_txt(p6, radius=30, total_sequences=2)
        met30 = parse_meme_txt(p30, radius=30, total_sequences=2)

        self.assertIsNotNone(met6)
        self.assertIsNotNone(met30)
        # Condition 2 (30 bp) must beat Condition 1 (6 bp)
        self.assertGreater(met30["composite_score"], met6["composite_score"])

    def test_insignificant_evalue_filtered_out(self):
        import tempfile

        # Motif with evalue = 2.5
        m_insig = """MEME version 4
MOTIF NOISE 1 width = 6 sites = 5 llr = 10 E-value = 2.5
letter-probability matrix: alength= 4 w= 6 nsites= 5 E= 2.5
  0.9 0.05 0.05 0.0
  0.9 0.05 0.05 0.0
  0.9 0.05 0.05 0.0
  0.0 0.9 0.05 0.05
  0.0 0.05 0.9 0.05
  0.05 0.05 0.0 0.9
Multilevel AAACGT
sites sorted by position p-value
s1 + 25 0.1
"""
        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as tf:
            tf.write(m_insig)
            p_tf = tf.name

        met = parse_meme_txt(p_tf, radius=30, total_sequences=1)
        self.assertIsNotNone(met)
        # Must have score 0 or be marked insignificant
        self.assertEqual(met["composite_score"], 0.0)

    def test_parallel_grid_runner_no_mpi_flag(self):
        import inspect
        from memegrid import grid_runner

        src = inspect.getsource(grid_runner)
        # Verify -p is NOT passed to meme command
        self.assertNotIn('"-p"', src)
        # Verify ProcessPoolExecutor is used
        self.assertIn("ProcessPoolExecutor", src)


if __name__ == "__main__":
    unittest.main()
