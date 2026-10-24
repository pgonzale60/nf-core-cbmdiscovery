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

from bipartite.aux_grid import score_auxiliary_candidate  # noqa: E402
from bipartite.masker import (  # noqa: E402
    load_soma_len_map,
    mask_sequences_with_hits,
)
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
        self.assertEqual(summary["pct_retained_domain"], 100.0)
        self.assertEqual(summary["pct_eliminated_domain"], 0.0)
        self.assertEqual(rows[0]["partner_domain"], "Retained")

        # Test dynamic pos0 mapping (e.g. border_1 with pos0=93)
        pos0_map = {"border_1": 93, "border_2": 150}
        rows_dyn, _summary_dyn = profile_bipartite_spacing(
            seq_names, core_hits, aux_hits, pos0_coord=pos0_map
        )
        # In border_1 with pos0=93: Core at 100..120 (+) has offset 93 - 120 = -27 bp
        self.assertEqual(rows_dyn[0]["core_offset_to_pos0"], -27)
        self.assertEqual(rows[0]["core_offset_to_pos0"], 30)

        # Test load_soma_len_map from FASTA headers
        with tempfile.NamedTemporaryFile("w", suffix=".fasta") as tf_fa:
            tf_fa.write(">border_1 soma_len=93\nACGT\n>border_2\nACGT\n")
            tf_fa.flush()
            m_lens = load_soma_len_map(tf_fa.name, default=150)
            self.assertEqual(m_lens["border_1"], 93)
            self.assertEqual(m_lens["border_2"], 150)

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
                "pct_retained_domain": 20.0,
                "pct_eliminated_domain": 80.0,
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
            self.assertIn("Located in Retained Region", report_text)
            self.assertIn("Located in Eliminated Region", report_text)

    def test_score_auxiliary_candidate_favors_longer_palindrome_over_core_hexamer(self):
        # 11 bp imperfect palindrome (e.g. ATAGCGCTAT) across 15 sequences
        pal_11bp = {
            "width": 11,
            "nsites": 15,
            "evalue": 1.5e-4,
            "total_ic": 16.5,
            "mean_ic": 1.50,
            "consensus_entropy": 1.75,
            "unique_sequences_hit": 15,
            "proximal_sites_6_50bp": 15,
        }
        score_11bp = score_auxiliary_candidate(pal_11bp, total_sequences=15)

        # 6 bp central core (AGCAGC) with higher per-column IC and multiple ANR hits
        core_6bp = {
            "width": 6,
            "nsites": 22,
            "evalue": 9.2e-6,
            "total_ic": 10.8,
            "mean_ic": 1.80,
            "consensus_entropy": 1.45,
            "unique_sequences_hit": 15,
            "proximal_sites_6_50bp": 18,
        }
        score_6bp = score_auxiliary_candidate(core_6bp, total_sequences=15)

        self.assertGreater(score_11bp, score_6bp)
        self.assertGreater(score_11bp, 14.0)

        # Auxiliary scoring must NOT reward proximity to cleavage site
        pal_distal = dict(pal_11bp, proximal_sites_6_50bp=0)
        self.assertEqual(
            score_auxiliary_candidate(pal_distal, total_sequences=15),
            score_11bp,
        )

        # Telomeric repeat echo in auxiliary candidate must be hard-rejected (score 0.0)
        telo_aux = {
            "consensus": "TTAGGCTTAGGC",
            "width": 12,
            "nsites": 15,
            "evalue": 1e-12,
            "total_ic": 20.0,
            "mean_ic": 1.67,
            "consensus_entropy": 1.85,
            "unique_sequences_hit": 15,
            "proximal_sites_6_50bp": 15,
        }
        self.assertEqual(score_auxiliary_candidate(telo_aux, total_sequences=15), 0.0)
        self.assertEqual(telo_aux["rejection_flag"], "REJECTED_TELOMERIC_ECHO")

    def test_satellite_window_and_paralog_flank_deduplication(self):
        bait_lib = REPO_ROOT / "modules" / "local" / "baitgermlinejunctions" / "resources" / "usr" / "lib"
        if str(bait_lib) not in sys.path:
            sys.path.insert(0, str(bait_lib))
        from baitgerm.census import deconvolve_and_census, deduplicate_break_windows, is_tandem_satellite_window

        # 47 bp satellite monomer repeating across 300 bp window (C. angaria PS1015 style)
        monomer_47 = "GCCGAAGGCCCAGACTCAGATTAGGCTTAGGCTTAGGCTTAACGTAC"
        sat_win = (monomer_47 * 7)[:300]
        self.assertTrue(is_tandem_satellite_window(sat_win))

        # Two paralogous windows sharing >90% identity in window[:130] but diverging in germline
        flank_150 = (
            "GCTAGCTGATCGATCGTAGCTAGCTAGCATCGATCGTAGCTAGCATCGAT"
            "TGCATGCAAAGGTTCCAAGGTTCCAAGGTTCCTAGCTGATCGATAGTCGA"
            "CCTTAAGGCCTTAAGGACTGACTGACTGTACGTACGCATGCTAGCTAGTC"
        )
        germ_a = ("AACCTTGG" * 19)[:150]
        germ_b = ("GGTTAACC" * 19)[:150]
        loci = [
            {"anchor_id": "p1", "anchor_seq": flank_150, "break_window": flank_150 + germ_a, "c_telo": 20, "c_germ": 5, "coherence": 0.95},
            {"anchor_id": "p2", "anchor_seq": flank_150, "break_window": flank_150 + germ_b, "c_telo": 10, "c_germ": 3, "coherence": 0.95},
            {"anchor_id": "sat", "anchor_seq": sat_win[:150], "break_window": sat_win, "c_telo": 30, "c_germ": 8, "coherence": 0.99},
        ]
        deduped = deduplicate_break_windows(loci)
        # p1 and p2 collapse due to shared somatic flank [:130]
        self.assertEqual(len(deduped), 2)
        p_merged = next(l for l in deduped if l["anchor_id"] in ("p1", "p2"))
        self.assertEqual(p_merged["collapsed_copies"], 2)

        classified, census = deconvolve_and_census(deduped)
        class_map = {l["anchor_id"]: l["classification"] for l in classified}
        self.assertEqual(class_map["sat"], "TANDEM_SATELLITE_ARRAY")
        self.assertEqual(census["n_satellite_excluded"], 1)


if __name__ == "__main__":
    unittest.main()
