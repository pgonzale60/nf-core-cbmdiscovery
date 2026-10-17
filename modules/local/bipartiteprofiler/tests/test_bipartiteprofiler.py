#!/usr/bin/env python3
"""Unit tests for bipartiteprofiler module."""

import pathlib
import sys
import tempfile
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
    / "bipartiteprofiler"
    / "resources"
    / "usr"
    / "lib"
)
if str(LIB_DIR) not in sys.path:
    sys.path.insert(0, str(LIB_DIR))

from bipartite.masker import mask_sequences_with_hits  # noqa: E402
from bipartite.reporter import generate_cbm_discovery_report  # noqa: E402
from bipartite.spacing import (  # noqa: E402
    calculate_edge_to_edge_spacer,
    profile_bipartite_spacing,
)


class TestBipartiteProfiler(unittest.TestCase):
    def test_edge_to_edge_spacer(self):
        # Case 1: Embedded
        # Core: 100..120. Aux: 105..115
        sp_emb, typ_emb = calculate_edge_to_edge_spacer(100, 120, 105, 115)
        self.assertEqual(sp_emb, 0)
        self.assertEqual(typ_emb, "Embedded")

        # Case 2: Upstream
        # Core: 100..120. Aux: 70..80 -> spacer = 100 - 80 = 20 bp
        sp_up, typ_up = calculate_edge_to_edge_spacer(100, 120, 70, 80)
        self.assertEqual(sp_up, 20)
        self.assertEqual(typ_up, "Flank_Upstream")

        # Case 3: Downstream
        # Core: 100..120. Aux: 135..145 -> spacer = 135 - 120 = 15 bp
        sp_dn, typ_dn = calculate_edge_to_edge_spacer(100, 120, 135, 145)
        self.assertEqual(sp_dn, 15)
        self.assertEqual(typ_dn, "Flank_Downstream")

    def test_mask_sequences_with_hits(self):
        seqs = {"seq1": "AAAAACCCCCGGGGGTTTTT"}
        hits = [
            {"seq_name": "seq1", "start": 6, "stop": 10}  # CCCCC -> NNNNN
        ]
        masked = mask_sequences_with_hits(seqs, hits)
        self.assertEqual(masked["seq1"], "AAAAANNNNNGGGGGTTTTT")

    def test_bipartite_spacing_profiling(self):
        seq_names = ["border_1", "border_2"]
        # Core at 100..120 (+), pos0 at 150 -> offset = 150 - 120 = 30 bp
        core_hits = {
            "border_1": [{"start": 100, "stop": 120, "strand": "+"}],
            "border_2": [{"start": 100, "stop": 120, "strand": "+"}],
        }
        # Aux at 70..80 on same strand (+) -> spacer = 20 bp
        aux_hits = {
            "border_1": [{"start": 70, "stop": 80, "strand": "+"}],
            "border_2": [{"start": 105, "stop": 115, "strand": "+"}],  # Embedded
        }

        rows, summary = profile_bipartite_spacing(
            seq_names, core_hits, aux_hits, pos0_coord=150
        )
        self.assertEqual(len(rows), 2)
        self.assertEqual(summary["n_paired_occurrences"], 2)
        self.assertEqual(summary["pct_embedded_0bp"], 50.0)
        self.assertEqual(summary["pct_spacer_le30bp"], 100.0)
        self.assertEqual(summary["pct_same_strand"], 100.0)

    def test_report_generation(self):
        with tempfile.NamedTemporaryFile("w", suffix=".json") as f_census, \
             tempfile.NamedTemporaryFile("w", suffix=".json") as f_core, \
             tempfile.NamedTemporaryFile("w", suffix=".md") as f_rep:

            f_census.write('{"n_breaks": 12, "median_c_soma": 25.0, "median_c_germ": 5.0, "germline_fraction": 0.167}')
            f_census.flush()

            f_core.write('{"consensus": "TGTTCA", "width": 6, "mean_ic": 1.95, "evalue_str": "1e-5", "occupancy": 0.92}')
            f_core.flush()

            bipartite_summary = {
                "n_paired_occurrences": 10,
                "median_spacer_bp": 15.0,
                "iqr_spacer_bp": 5.0,
                "pct_embedded_0bp": 20.0,
                "pct_spacer_le30bp": 90.0,
                "pct_spacer_le50bp": 100.0,
                "pct_same_strand": 90.0,
            }

            best_aux = {
                "consensus": "WWTWATGAACA",
                "width": 11,
                "mean_ic": 1.45,
                "evalue_str": "1e-4",
                "proximal_sites_6_50bp": 8,
            }

            generate_cbm_discovery_report(
                sample_id="test_species",
                census_json_path=f_census.name,
                core_metrics_json_path=f_core.name,
                best_aux_motif=best_aux,
                bipartite_summary=bipartite_summary,
                output_report_md=f_rep.name,
            )

            report_text = pathlib.Path(f_rep.name).read_text(encoding="utf-8")
            self.assertIn("CBM Discovery & PDE Junction Report: `test_species`", report_text)
            self.assertIn("TGTTCA", report_text)
            self.assertIn("WWTWATGAACA", report_text)
            self.assertIn("Median Edge-to-Edge Spacer", report_text)


if __name__ == "__main__":
    unittest.main()
