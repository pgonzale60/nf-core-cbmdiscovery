#!/usr/bin/env python3
"""Unit tests for baitgermlinejunctions module."""

import os
import pathlib
import subprocess
import sys
import tempfile
import unittest
import unittest.mock

# Ensure module lib is importable
CUR = pathlib.Path(__file__).resolve().parent
REPO_ROOT = CUR
while REPO_ROOT.parent != REPO_ROOT and not (REPO_ROOT / "modules").is_dir():
    REPO_ROOT = REPO_ROOT.parent

LIB_DIR = (
    REPO_ROOT
    / "modules"
    / "local"
    / "baitgermlinejunctions"
    / "resources"
    / "usr"
    / "lib"
)
if str(LIB_DIR) not in sys.path:
    sys.path.insert(0, str(LIB_DIR))

from baitgerm.baiter import (  # noqa: E402
    build_kmer_bait_index,
    revcomp,
)
from baitgerm.census import (  # noqa: E402
    deconvolve_and_census,
    deduplicate_break_windows,
)
from baitgerm.partitioner import partition_baited_read  # noqa: E402
from baitgerm.reconstructor import (  # noqa: E402
    annotate_microhomology,
    build_300bp_break_window,
    reconstruct_germline_consensus,
)


class TestBaitGermlineJunctions(unittest.TestCase):
    def setUp(self):
        self.repeat = "TTAGGC"
        # 150 bp somatic anchor ending at pos 0
        self.anchor_seq = (
            "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT"
            "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
            "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGCG"
        )
        self.assertTrue(len(self.anchor_seq) == 150)
        # 150 bp germline flank starting at +1 with 'TTAG' microhomology
        self.germ_flank = (
            "TTAG"
            + "CCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCC"
            + "GGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGG"
            + "CCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCC"
        )[:150]
        self.assertTrue(len(self.germ_flank) == 150)

    def test_bait_index_construction(self):
        anchors = {"anchor_1": self.anchor_seq}
        bait_index = build_kmer_bait_index(
            anchors, kmer_size=31, terminal_window=60
        )
        self.assertGreater(len(bait_index), 0)
        # Terminal 31-mer of anchor
        term_kmer = self.anchor_seq[-31:]
        self.assertIn(term_kmer, bait_index)
        self.assertIn("anchor_1", bait_index[term_kmer])
        # Revcomp should also be present
        self.assertIn(revcomp(term_kmer), bait_index)

    def test_partitioning_telo_vs_germline(self):
        # Telomeric read
        telo_read = self.anchor_seq[-60:] + (self.repeat * 8)
        p_telo = partition_baited_read(
            self.anchor_seq, telo_read, repeat=self.repeat
        )
        self.assertIsNotNone(p_telo)
        self.assertEqual(p_telo["branch"], "telo")

        # Germline read-through read
        germ_read = self.anchor_seq[-60:] + self.germ_flank[:80]
        p_germ = partition_baited_read(
            self.anchor_seq, germ_read, repeat=self.repeat
        )
        self.assertIsNotNone(p_germ)
        self.assertEqual(p_germ["branch"], "germline")
        self.assertEqual(p_germ["germline_ext"], self.germ_flank[:80])

        # Telomeric read with telomerase addition wobble (non-canonical overhang before repeats)
        wobble_read = self.anchor_seq[-60:] + "ACTTGGCTAC" + (self.repeat * 8)
        p_wobble = partition_baited_read(
            self.anchor_seq, wobble_read, repeat=self.repeat
        )
        self.assertIsNotNone(p_wobble)
        self.assertEqual(p_wobble["branch"], "telo")

        # Germline read with short 12 bp microhomology (2 repeats) followed by unique sequence
        mh_read = self.anchor_seq[-60:] + "TAGGCTTAGGCT" + ("GAGCCTATTCCGAAGAGCTCGGCTTC" * 3)
        p_mh = partition_baited_read(
            self.anchor_seq, mh_read, repeat=self.repeat
        )
        self.assertIsNotNone(p_mh)
        self.assertEqual(p_mh["branch"], "germline")

    def test_reconstruction_and_microhomology(self):
        # Multiple reads with small sequencing noise
        exts = [
            self.germ_flank,
            self.germ_flank[:-5] + "A" * 5,
            self.germ_flank,
        ]
        cons, coherence = reconstruct_germline_consensus(exts, max_len=150)
        self.assertEqual(cons, self.germ_flank)
        self.assertGreaterEqual(coherence, 0.95)

        # 1-bp homopolymer indel shift tolerance (HiFi noise)
        indel_ext = self.germ_flank[:25] + self.germ_flank[26:]
        _c_indel, coh_indel = reconstruct_germline_consensus(
            [self.germ_flank, indel_ext, self.germ_flank], max_len=150
        )
        self.assertGreaterEqual(coh_indel, 0.95)

        mh_len, mh_motif = annotate_microhomology(cons, repeat=self.repeat)
        self.assertEqual(mh_len, 4)
        self.assertEqual(mh_motif, "TTAG")

        win = build_300bp_break_window(self.anchor_seq, cons)
        self.assertEqual(len(win), 300)
        self.assertEqual(win[:150], self.anchor_seq)
        self.assertEqual(win[150:], self.germ_flank)

    def test_chimeric_noise_rejection(self):
        # Highly discordant reads
        exts = [
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
            "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
            "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
        ]
        _cons, coherence = reconstruct_germline_consensus(exts, max_len=40)
        self.assertLess(coherence, 0.85)

        loci = [{
            "anchor_id": "anchor_noise",
            "c_telo": 15,
            "c_germ": 3,
            "coherence": coherence,
            "break_window": "",
        }]
        classified, _census = deconvolve_and_census(
            loci, min_germ_reads=2, min_coherence=0.85
        )
        self.assertNotEqual(classified[0]["classification"], "PDE_BREAK_SITE")

    def test_deconvolve_and_census(self):
        loci = [
            {"anchor_id": "a_break", "c_telo": 24, "c_germ": 6, "coherence": 0.98, "break_window": "ACGT" * 75},
            {"anchor_id": "a_soft", "c_telo": 39, "c_germ": 18, "coherence": 0.78, "break_window": "TGCA" * 75},
            {"anchor_id": "a_native", "c_telo": 5, "c_germ": 0, "coherence": 0.0, "break_window": ""},
            {"anchor_id": "a_low_coh", "c_telo": 20, "c_germ": 4, "coherence": 0.45, "break_window": ""},
            {"anchor_id": "a_unbroken", "c_telo": 0, "c_germ": 67, "coherence": 0.99, "break_window": ""},
            {"anchor_id": "a_rdna_telo", "c_telo": 109, "c_germ": 5333, "coherence": 0.996, "break_window": ""},
        ]
        deduped = deduplicate_break_windows(loci)
        classified, census = deconvolve_and_census(deduped, min_germ_reads=2)

        class_map = {l["anchor_id"]: l["classification"] for l in classified}
        self.assertEqual(class_map["a_break"], "PDE_BREAK_SITE")
        self.assertEqual(class_map["a_soft"], "PDE_BREAK_SITE")
        self.assertEqual(class_map["a_native"], "GERMLINE_NATIVE_TELOMERE")
        self.assertEqual(class_map["a_low_coh"], "LOW_COHERENCE_JUNCTION")
        self.assertEqual(class_map["a_unbroken"], "UNBROKEN_GERMLINE_LOCUS")
        self.assertEqual(class_map["a_rdna_telo"], "GERMLINE_NATIVE_TELOMERE")

        self.assertEqual(census["n_breaks"], 2)
        self.assertEqual(census["n_germ_telo"], 2)
        self.assertEqual(census["median_c_soma"], 31.5)
        self.assertEqual(census["median_c_germ"], 12.0)
        self.assertAlmostEqual(census["germline_fraction"], 0.28, places=2)

    def test_length_shifted_break_window_deduplication_and_count_merging(self):
        soma_150 = self.anchor_seq
        germ_150 = self.germ_flank
        loci = [
            {
                "anchor_id": "anchor_1",
                "anchor_seq": soma_150,
                "germline_consensus": germ_150,
                "break_window": soma_150 + germ_150,
                "c_telo": 20,
                "c_germ": 4,
            },
            {
                "anchor_id": "anchor_2",
                "anchor_seq": soma_150[-120:],
                "germline_consensus": germ_150[:130],
                "break_window": soma_150[-120:] + germ_150[:130],
                "c_telo": 5,
                "c_germ": 2,
            },
        ]
        deduped = deduplicate_break_windows(loci, min_identity=0.95)
        self.assertEqual(len(deduped), 1)
        self.assertEqual(deduped[0]["c_telo"], 25)
        self.assertEqual(deduped[0]["c_germ"], 6)

    def test_illumina_short_prefix_and_mate_extension(self):
        # 150 bp somatic anchor ending at pos 0
        soma = self.anchor_seq
        # 210 bp germline sequence starting at +1
        germ_210 = (self.germ_flank + "ACGTACGT" * 10)[:210]

        # Read A: [-45 .. +105] -> 45 bp soma + 105 bp germ
        read_a = soma[-45:] + germ_210[:105]
        # Read B: [-20 .. +130] -> 20 bp soma (starts at -20) + 130 bp germ
        read_b = soma[-20:] + germ_210[:130]
        # Read C: [+60 .. +210] -> 150 bp germline mate
        read_c = germ_210[60:210]

        p_a = partition_baited_read(soma, read_a, repeat=self.repeat)
        p_b = partition_baited_read(soma, read_b, repeat=self.repeat)

        self.assertIsNotNone(p_a)
        self.assertIsNotNone(p_b)
        self.assertEqual(p_a["branch"], "germline")
        self.assertEqual(p_b["branch"], "germline")

        # Initial extensions from junction reads
        initial_exts = [p_a["germline_ext"], p_b["germline_ext"]]
        mock_raw = f"@utg0\n{germ_210}\n+\n{'I' * len(germ_210)}\n"
        with unittest.mock.patch("subprocess.run") as mock_sub:
            mock_sub.return_value = unittest.mock.Mock(stdout=mock_raw, returncode=0)
            cons, _coh = reconstruct_germline_consensus(
                initial_exts, max_len=150, mate_reads=[read_c]
            )
        self.assertEqual(len(cons), 150)
        self.assertEqual(cons, germ_210[:150])

    def test_end_to_end_cli_illumina_short_prefix_and_extension(self):
        soma = self.anchor_seq
        germ_210 = (self.germ_flank + "ACGTACGT" * 10)[:210]

        read_a = soma[-45:] + germ_210[:105]
        read_b = soma[-20:] + germ_210[:130]
        read_c = germ_210[60:210]

        with tempfile.TemporaryDirectory() as td:
            tdp = pathlib.Path(td)
            bin_dir = tdp / "bin"
            bin_dir.mkdir()
            mock_fml = bin_dir / "fml-asm"
            mock_fml.write_text(
                f"#!/bin/sh\nprintf '@utg0\\n{soma[-50:] + germ_210}\\n+\\n{'I' * (50 + len(germ_210))}\\n'\n"
            )
            mock_fml.chmod(0o755)

            fa_cand = tdp / "candidate_anchors.fasta"
            fa_cand.write_text(f">anchor_1\n{soma}\n>anchor_2\n{soma}\n")

            r1_fq = tdp / "r1.fastq"
            telo_reads = "\n".join([
                f"@telo_{i}/1\n{soma[-60:] + self.repeat * 8}\n+\n{'I' * (60 + len(self.repeat) * 8)}"
                for i in range(5)
            ])
            dummy_reads = "\n".join([
                f"@dummy_{i}\nAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTTAAAAACCCCC\n+\n" + ("I" * 50)
                for i in range(5)
            ])
            r1_content = (
                f"@read_a/1\n{read_a}\n+\n{'I' * len(read_a)}\n"
                f"@read_b/1\n{read_b}\n+\n{'I' * len(read_b)}\n"
                f"{telo_reads}\n"
                f"{dummy_reads}\n"
            )
            r1_fq.write_text(r1_content)

            r2_fq = tdp / "r2.fastq"
            r2_fq.write_text(f"@read_c/2\n{read_c}\n+\n{'I' * len(read_c)}\n")

            out_fa = tdp / "break_windows_300bp.fasta"
            out_tsv = tdp / "junction_catalog.tsv"
            out_json = tdp / "census.json"

            cli_script = (
                REPO_ROOT / "modules" / "local" / "baitgermlinejunctions" / "resources" / "usr" / "bin" / "cbm_bait_germline_junctions.py"
            ).resolve()

            cmd = [
                sys.executable, str(cli_script),
                "--candidate-anchors", str(fa_cand),
                "--reads-1", str(r1_fq), "--reads-2", str(r2_fq),
                "--platform", "illumina", "--min-germ-reads", "2",
                "--min-telo-reads", "5",
                "--threads", "2",
                "--output-catalog", str(out_tsv),
                "--output-fasta", str(out_fa),
                "--output-census", str(out_json),
            ]
            env = dict(os.environ, PATH=f"{bin_dir}:{os.environ.get('PATH', '')}")
            subprocess.run(cmd, capture_output=True, text=True, check=True, env=env)

            out_content = out_fa.read_text()
            lines = [l.strip() for l in out_content.splitlines() if l.strip() and not l.startswith(">")]
            seq = "".join(lines)
            self.assertEqual(len(seq), len(soma) + len(germ_210))
            self.assertEqual(seq[:150], soma)
            self.assertEqual(seq[150:], germ_210)

    def test_cli_illumina_single_end_rejected(self):
        cli_script = (
            pathlib.Path(__file__).parent.parent / "resources" / "usr" / "bin" / "cbm_bait_germline_junctions.py"
        ).resolve()

        with tempfile.TemporaryDirectory() as td:
            tdp = pathlib.Path(td)
            fa = tdp / "cands.fa"
            fa.write_text(">cand1\n" + self.anchor_seq + "\n")
            r1 = tdp / "r1.fq"
            r1.write_text("@r1\nACGT\n+\nIIII\n")

            cmd = [
                sys.executable, str(cli_script),
                "--candidate-anchors", str(fa),
                "--reads-1", str(r1),
                "--platform", "illumina",
            ]
            res = subprocess.run(cmd, capture_output=True, text=True)
            self.assertNotEqual(res.returncode, 0)
            self.assertIn("requires paired-end reads", res.stderr)


if __name__ == "__main__":
    unittest.main()

