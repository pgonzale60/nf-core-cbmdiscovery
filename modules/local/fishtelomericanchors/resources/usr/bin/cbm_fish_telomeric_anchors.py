#!/usr/bin/env python3
"""Stage 1: Fish telomeric reads and construct candidate somatic anchors."""

import argparse
import json
from pathlib import Path
import sys

# Ensure usr/lib is in python search path
_LIB_DIR = Path(__file__).resolve().parent.parent / "lib"
if str(_LIB_DIR) not in sys.path:
    sys.path.insert(0, str(_LIB_DIR))

from teloanchor.cluster import cluster_anchors  # noqa: E402
from teloanchor.detector import find_terminal_telomeric_tract  # noqa: E402
from teloanchor.orient import normalize_and_trim_flank  # noqa: E402
from teloanchor.reader import stream_fastq  # noqa: E402


def parse_args():
    parser = argparse.ArgumentParser(
        description="Fish telomeric reads and build candidate somatic anchors centered at pos 0"
    )
    parser.add_argument(
        "--reads-1", required=True, help="Path to FASTQ reads 1"
    )
    parser.add_argument(
        "--reads-2",
        default=None,
        help="Path to FASTQ reads 2 (optional, paired-end)",
    )
    parser.add_argument(
        "--platform",
        choices=["hifi", "illumina"],
        default="hifi",
        help="Sequencing platform (hifi or illumina)",
    )
    parser.add_argument(
        "--telomere-repeat",
        default="TTAGGC",
        help="Canonical telomeric repeat (default: TTAGGC)",
    )
    parser.add_argument(
        "--output-fasta",
        default="candidate_anchors.fasta",
        help="Output candidate anchors FASTA",
    )
    parser.add_argument(
        "--output-summary",
        default="pass1_summary.json",
        help="Output summary JSON",
    )
    parser.add_argument(
        "--min-cluster-size",
        type=int,
        default=1,
        help="Minimum reads per anchor cluster",
    )
    return parser.parse_args()


def process_read(qname, seq, platform, repeat, max_flank):
    tract = find_terminal_telomeric_tract(seq, repeat=repeat, platform=platform)
    if not tract:
        return None
    flank_res = normalize_and_trim_flank(
        seq, tract, repeat=repeat, max_flank_len=max_flank
    )
    if not flank_res:
        return None
    return {
        "read_id": qname,
        "anchor_seq": flank_res["anchor_seq"],
        "trimmed_prefix": flank_res["trimmed_prefix"],
        "microhomology_len": flank_res["microhomology_trimmed_len"],
        "terminal_end": tract["terminal_end"],
        "tract_len": tract["tract_len"],
    }


def main():
    args = parse_args()
    platform = args.platform.lower()
    repeat = args.telomere_repeat.upper()
    max_flank = 150 if platform == "hifi" else 150

    total_reads = 0
    fished_observations = []

    read_files = [args.reads_1]
    if args.reads_2:
        read_files.append(args.reads_2)

    for rf in read_files:
        for qname, seq, _qual in stream_fastq(rf):
            total_reads += 1
            obs = process_read(qname, seq, platform, repeat, max_flank)
            if obs:
                fished_observations.append(obs)

    clusters = cluster_anchors(
        fished_observations,
        min_identity=0.97,
        min_cluster_size=args.min_cluster_size,
    )

    # Write candidate anchors FASTA
    with open(args.output_fasta, "w", encoding="utf-8") as f_fa:
        for cl in clusters:
            f_fa.write(
                f">{cl['anchor_id']} read_count={cl['read_count']} len={len(cl['consensus_seq'])}\n"
            )
            f_fa.write(f"{cl['consensus_seq']}\n")

    summary = {
        "platform": platform,
        "telomere_repeat": repeat,
        "total_reads_scanned": total_reads,
        "fished_telomeric_reads": len(fished_observations),
        "candidate_anchors_count": len(clusters),
        "clusters": [
            {
                "anchor_id": cl["anchor_id"],
                "read_count": cl["read_count"],
                "anchor_len": len(cl["consensus_seq"]),
            }
            for cl in clusters
        ],
    }

    with open(args.output_summary, "w", encoding="utf-8") as f_json:
        json.dump(summary, f_json, indent=2)

    print(
        f"[FISHTELOMERICANCHORS] Processed {total_reads} reads -> Fished {len(fished_observations)} telomeric reads -> Clustered {len(clusters)} candidate anchors."
    )


if __name__ == "__main__":
    main()
