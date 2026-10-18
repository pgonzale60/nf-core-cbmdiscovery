#!/usr/bin/env python3
"""Comprehensive TDD Test Suite for v3.2 Tool-First Architecture."""

import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

# Ensure module libs are importable
CUR = Path(__file__).resolve().parent
REPO_ROOT = CUR
while REPO_ROOT.parent != REPO_ROOT and not (REPO_ROOT / "modules").is_dir():
    REPO_ROOT = REPO_ROOT.parent

for mod in ["fishtelomericanchors", "baitgermlinejunctions", "memecoregrid", "bipartiteprofiler"]:
    lib_p = REPO_ROOT / "modules" / "local" / mod / "resources" / "usr" / "lib"
    if str(lib_p) not in sys.path:
        sys.path.insert(0, str(lib_p))

from teloanchor.cluster import cluster_anchors  # noqa: E402
from teloanchor.detector import find_terminal_telomeric_tract  # noqa: E402
from teloanchor.orient import normalize_and_trim_flank  # noqa: E402
from memegrid.grid_runner import run_meme_grid  # noqa: E402


def revcomp(seq: str) -> str:
    tr = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(tr)[::-1]


class TestV3Optimizations(unittest.TestCase):
    def test_no_hardcoded_cheats_and_stage1_vsearch_clustering(self):
        """1. Assert detector.py has no hardcoded cheats, and vsearch clusters staggered flanks."""
        det_path = (
            REPO_ROOT
            / "modules"
            / "local"
            / "fishtelomericanchors"
            / "resources"
            / "usr"
            / "lib"
            / "teloanchor"
            / "detector.py"
        )
        det_code = det_path.read_text(encoding="utf-8")
        for cheat in ['"AA"', '"CC"', '"GC"']:
            self.assertNotIn(cheat, det_code)

        # 150 bp master somatic flank ending in 'CT'
        soma_150 = (
            "GATTACACATTAGACCGGATACACAGATAGAGACAGATAGACAGATAGAC"
            "ACGTACGTACGTTTGCAATCGATCGATCGATCGATCGATCGATCGATCGA"
            "TGCTAGCTAGCTAGCTAAGGTCCAAATTGGCCAAGGTTAAGGCCTTGGCT"
        )
        self.assertEqual(len(soma_150), 150)
        self.assertTrue(soma_150.endswith("CT"))

        repeat = "TTAGGC"
        staggered_lengths = [65, 80, 100, 120, 150]
        observations = []

        for idx, fl_len in enumerate(staggered_lengths):
            flank = soma_150[-fl_len:]
            # Varied in-phase telomeric additions
            if idx % 3 == 0:
                telo_add = "TT" + "AGGC" + (repeat * 10)
            elif idx % 3 == 1:
                telo_add = "TTAG" + "GC" + (repeat * 10)
            else:
                telo_add = "GGC" + "TTAGGC" + (repeat * 9)

            read_seq = flank + telo_add + "A" * 15
            tract = find_terminal_telomeric_tract(
                read_seq, repeat=repeat, platform="illumina", min_flank_bp=50
            )
            self.assertIsNotNone(tract)
            trimmed = normalize_and_trim_flank(
                read_seq, tract, repeat=repeat, max_flank_len=150
            )
            self.assertIsNotNone(trimmed)
            self.assertTrue(trimmed["anchor_seq"].endswith("CT"))
            observations.append({
                "read_id": f"staggered_{fl_len}",
                "anchor_seq": trimmed["anchor_seq"],
                "trimmed_prefix": trimmed["trimmed_prefix"],
                "microhomology_len": trimmed["microhomology_trimmed_len"],
                "terminal_end": tract["terminal_end"],
                "tract_len": tract["tract_len"],
            })

        # Cluster staggered flanks using vsearch --cluster_fast --iddef 0
        clusters = cluster_anchors(observations, min_identity=0.97, min_cluster_size=1)
        self.assertEqual(len(clusters), 1)
        self.assertEqual(clusters[0]["read_count"], len(staggered_lengths))
        self.assertEqual(len(clusters[0]["consensus_seq"]), 150)
        self.assertEqual(clusters[0]["consensus_seq"], soma_150)
        self.assertTrue(clusters[0]["consensus_seq"].endswith("CT"))

    def test_bbduk_paired_end_sync_and_fermi_lite_asymmetric_y_junction(self):
        """2. Assert bbduk.sh mate synchronization and fermi-lite Y-junction reconstruction."""
        repeat = "TTAGGC"
        soma_150 = (
            "GATTACACATTAGACCGGATACACAGATAGAGACAGATAGACAGATAGAC"
            "ACGTACGTACGTTTGCAATCGATCGATCGATCGATCGATCGATCGATCGA"
            "TGCTAGCTAGCTAGCTAAGGTCCAAATTGGCCAAGGTTAAGGCCTTGGCT"
        )
        germline_180 = (
            "TTAGCTAGCTAGCTAGCATCGATCGATCGATCGATCGATCGATCGATTGCA"
            "AACGTACGTACGTGTCTATCTGTCTATCTGTCTCTATCTGTGTATCCGGTC"
            "TAATGTGTAATCAGCTAGCTAGCTAGCATCGATCGATCGATCGATCGATCG"
            "ATCGAATGCTAGCTAGCTAG"
        )
        self.assertTrue(germline_180.startswith("TTAG"))
        self.assertGreaterEqual(len(germline_180), 170)

        with tempfile.TemporaryDirectory() as td:
            tdp = Path(td)
            fa_cand = tdp / "candidate_anchors.fasta"
            fa_cand.write_text(f">anchor_1\n{soma_150}\n")
            fa_bait = tdp / "anchors_terminal_60bp.fasta"
            fa_bait.write_text(f">anchor_1\n{soma_150[-60:]}\n")

            r1_path = tdp / "raw_R1.fastq"
            r2_path = tdp / "raw_R2.fastq"

            # 18 somatic telomeric pairs (90%)
            # R1: spans [-80 .. 0] + (TTAGGC * 10)
            telo_r1 = soma_150[-80:] + (repeat * 10)
            telo_r2 = revcomp(soma_150[:100])

            # 3 germline pairs (10%)
            # R1: spans [-50 .. +60] -> 110 bp
            germ_r1 = soma_150[-50:] + germline_180[:60]
            # R2: spans [+30 .. +180] -> 150 bp
            germ_r2 = revcomp(germline_180[30:180])

            # 100 background pairs
            bg_r1 = "C" * 75 + "G" * 75
            bg_r2 = "A" * 75 + "T" * 75

            with open(r1_path, "w", encoding="utf-8") as f1, open(r2_path, "w", encoding="utf-8") as f2:
                for i in range(18):
                    f1.write(f"@telo_{i}/1\n{telo_r1}\n+\n{'I' * len(telo_r1)}\n")
                    f2.write(f"@telo_{i}/2\n{telo_r2}\n+\n{'I' * len(telo_r2)}\n")
                for i in range(3):
                    f1.write(f"@germ_{i}/1\n{germ_r1}\n+\n{'I' * len(germ_r1)}\n")
                    f2.write(f"@germ_{i}/2\n{germ_r2}\n+\n{'I' * len(germ_r2)}\n")
                for i in range(100):
                    f1.write(f"@bg_{i}/1\n{bg_r1}\n+\n{'I' * len(bg_r1)}\n")
                    f2.write(f"@bg_{i}/2\n{bg_r2}\n+\n{'I' * len(bg_r2)}\n")

            baited_r1 = tdp / "baited_R1.fastq"
            baited_r2 = tdp / "baited_R2.fastq"

            # Run bbduk.sh with paired-end mate synchronization
            bbduk_cmd = [
                "bbduk.sh",
                f"in1={r1_path}",
                f"in2={r2_path}",
                f"ref={fa_bait}",
                f"outm1={baited_r1}",
                f"outm2={baited_r2}",
                "k=27",
                "hdist=1",
                "rcomp=t",
                "minkmerhits=1",
                "entropy=0.5",
                "overwrite=true",
            ]
            subprocess.run(bbduk_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            # Assert bbduk baited both R1 and mate R2
            b1_text = baited_r1.read_text()
            b2_text = baited_r2.read_text()
            self.assertIn("@germ_0/1", b1_text)
            self.assertIn("@germ_0/2", b2_text)

            # Run cbm_bait_germline_junctions.py
            cli_script = (
                REPO_ROOT
                / "modules"
                / "local"
                / "baitgermlinejunctions"
                / "resources"
                / "usr"
                / "bin"
                / "cbm_bait_germline_junctions.py"
            )
            out_cat = tdp / "junction_catalog.tsv"
            out_fa = tdp / "break_windows_300bp.fasta"
            out_census = tdp / "census.json"

            cmd = [
                sys.executable,
                str(cli_script),
                "--candidate-anchors", str(fa_cand),
                "--reads-1", str(baited_r1),
                "--reads-2", str(baited_r2),
                "--platform", "illumina",
                "--min-germ-reads", "2",
                "--output-catalog", str(out_cat),
                "--output-fasta", str(out_fa),
                "--output-census", str(out_census),
            ]
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            # Assert break window is 300 bp and accurately reconstructed across Y-junction
            fa_lines = [l.strip() for l in out_fa.read_text().splitlines() if l.strip() and not l.startswith(">")]
            self.assertTrue(len(fa_lines) > 0)
            window = "".join(fa_lines)
            self.assertEqual(len(window), 300)
            self.assertEqual(window[:150], soma_150)
            self.assertTrue(window[150:].startswith("TTAG"))
            self.assertEqual(window[150:], germline_180[:150])

            # Assert microhomology is 4 bp 'TTAG'
            cat_text = out_cat.read_text()
            self.assertIn("TTAG", cat_text)
            self.assertIn("\t4\tTTAG\t", cat_text)

    def test_markov_background_and_hierarchical_meme_grid(self):
        """3. Assert distal Markov background, hierarchical grid <= 21 runs, and SpaMo integration."""
        core_motif = "GATCGTACGATCGT"  # 14 bp
        aux_motif = "ATGGCT"            # 6 bp
        poly_noise = "AAAAAAAATTTTTTTT"

        with tempfile.TemporaryDirectory() as td:
            tdp = Path(td)
            fa_breaks = tdp / "break_windows_300bp.fasta"

            # 12 synthetic 300 bp windows
            with open(fa_breaks, "w", encoding="utf-8") as f:
                for i in range(12):
                    # Distal 5' arm [-150 .. -60] = 90 bp containing poly-A/T noise
                    left_arm = (poly_noise * 6)[:90]
                    # Core region [-60 .. 0] containing core motif
                    core_flank = ("CGTATAGCCTA" * 5)[:46] + core_motif  # 60 bp
                    # Proximal germline [0 .. +60] containing aux motif at +25
                    prox_germ = ("CGATTAGCCTA" * 3)[:25] + aux_motif + ("TTAACCGGAT" * 3)[:29]  # 60 bp
                    # Distal 3' arm [+60 .. +150] = 90 bp containing poly-A/T noise
                    right_arm = (poly_noise * 6)[:90]

                    win_300 = left_arm + core_flank + prox_germ + right_arm
                    self.assertEqual(len(win_300), 300)
                    f.write(f">break_site_{i}\n{win_300}\n")

            work_dir = tdp / "grid_work"
            best_motif, grid_results = run_meme_grid(
                break_windows_fasta=fa_breaks,
                output_dir=work_dir,
                threads=2,
            )

            # Assert distal background generated
            self.assertTrue((work_dir / "backgrounds" / "distal_flanks_bg.fasta").exists())
            self.assertTrue((work_dir / "best_background.bfile").exists())

            # Assert hierarchical grid ran <= 21 conditions
            self.assertIsNotNone(best_motif)
            self.assertLessEqual(len(grid_results), 21)

            # Assert recovered motif matches core motif or its reverse complement
            recovered = best_motif["consensus"]
            rc_core = revcomp(core_motif)
            self.assertTrue(
                core_motif in recovered or rc_core in recovered or recovered in core_motif or recovered in rc_core
            )

            # Run Stage 5 CLI
            stage5_script = (
                REPO_ROOT
                / "modules"
                / "local"
                / "bipartiteprofiler"
                / "resources"
                / "usr"
                / "bin"
                / "cbm_bipartite_profiler.py"
            )
            out_aux_meme = tdp / "best_auxiliary_motif.meme"
            out_spacing = tdp / "bipartite_spacing_summary.tsv"
            out_report = tdp / "cbm_discovery_report.md"
            out_masked = tdp / "break_sequences.core_masked.fasta"
            census_json = tdp / "census.json"
            census_json.write_text(json.dumps({"n_breaks": 12, "median_c_soma": 10.0, "median_c_germ": 2.0, "germline_fraction": 0.2}))
            core_json = tdp / "core_metrics.json"
            core_json.write_text(json.dumps(best_motif))

            cmd5 = [
                sys.executable,
                str(stage5_script),
                "--break-windows", str(fa_breaks),
                "--core-meme", str(best_motif["meme_txt_path"]),
                "--census-json", str(census_json),
                "--core-metrics-json", str(core_json),
                "--background-bfile", str(work_dir / "best_background.bfile"),
                "--output-aux-meme", str(out_aux_meme),
                "--output-spacing-tsv", str(out_spacing),
                "--output-report", str(out_report),
                "--output-masked-fasta", str(out_masked),
                "--threads", "2",
            ]
            subprocess.run(cmd5, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            self.assertTrue(out_aux_meme.exists())
            self.assertTrue(out_spacing.exists())
            self.assertTrue(out_report.exists())
            self.assertIn("ATGGCT", out_report.read_text())


if __name__ == "__main__":
    unittest.main()
