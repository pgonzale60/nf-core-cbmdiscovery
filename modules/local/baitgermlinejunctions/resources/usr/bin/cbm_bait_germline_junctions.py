#!/usr/bin/env python3
"""Stage 2 & 3: K-mer Baiting, Germline Reconstruction & Census."""

import argparse
import csv
import json
from pathlib import Path
import sys

# Ensure usr/lib is in python search path
_LIB_DIR = Path(__file__).resolve().parent.parent / "lib"
if str(_LIB_DIR) not in sys.path:
    sys.path.insert(0, str(_LIB_DIR))

from baitgerm.baiter import (  # noqa: E402
    build_kmer_bait_index,
    load_candidate_anchors,
    stream_paired_fastq_records,
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


def parse_args():
    parser = argparse.ArgumentParser(
        description="Bait germline reads, reconstruct break windows, and perform coverage/junction census"
    )
    parser.add_argument(
        "--candidate-anchors",
        required=True,
        help="Path to candidate_anchors.fasta",
    )
    parser.add_argument(
        "--reads-1", required=True, help="Path to FASTQ reads 1"
    )
    parser.add_argument(
        "--reads-2", default=None, help="Path to FASTQ reads 2 (optional)"
    )
    parser.add_argument(
        "--platform",
        choices=["hifi", "illumina"],
        default="hifi",
        help="Sequencing platform",
    )
    parser.add_argument(
        "--telomere-repeat", default="TTAGGC", help="Canonical telomere repeat"
    )
    parser.add_argument(
        "--kmer-size", type=int, default=31, help="Bait k-mer size"
    )
    parser.add_argument(
        "--min-germ-reads",
        type=int,
        default=2,
        help="Minimum germline read-through reads",
    )
    parser.add_argument(
        "--min-coherence",
        type=float,
        default=0.85,
        help="Minimum germline branch coherence",
    )
    parser.add_argument(
        "--output-catalog",
        default="junction_catalog.tsv",
        help="Output TSV junction catalog",
    )
    parser.add_argument(
        "--output-fasta",
        default="break_windows_300bp.fasta",
        help="Output FASTA of 300 bp break windows",
    )
    parser.add_argument(
        "--output-census",
        default="coverage_and_chromosome_census.json",
        help="Output census JSON",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads for reconstruction and deduplication",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    anchors = load_candidate_anchors(args.candidate_anchors)
    if not anchors:
        print(f"Warning: No anchors found in {args.candidate_anchors}")

    bait_k = min(args.kmer_size, 20)
    bait_index = build_kmer_bait_index(anchors, kmer_size=bait_k)

    def get_hits(read_seq: str, bait_idx: dict[str, list[str]], k_sz: int) -> set[str]:
        hits = set()
        if len(read_seq) < k_sz:
            return hits
        for i in range(len(read_seq) - k_sz + 1):
            km = read_seq[i : i + k_sz]
            if km in bait_idx:
                hits.update(bait_idx[km])
        return hits

    # Track extensions, local reads, and counts per anchor
    telo_counts = {a_id: 0 for a_id in anchors}
    germline_extensions = {a_id: [] for a_id in anchors}
    local_reads: dict[str, list[str]] = {a_id: [] for a_id in anchors}

    for r1_tuple, r2_tuple in stream_paired_fastq_records(args.reads_1, args.reads_2):
        s1 = r1_tuple[1] if r1_tuple else None
        s2 = r2_tuple[1] if r2_tuple else None
        hits_1 = get_hits(s1, bait_index, bait_k) if s1 else set()
        hits_2 = get_hits(s2, bait_index, bait_k) if s2 else set()
        all_hits = hits_1 | hits_2

        for a_id in all_hits:
            a_seq = anchors[a_id]
            p1 = partition_baited_read(a_seq, s1, repeat=args.telomere_repeat) if s1 else None
            p2 = partition_baited_read(a_seq, s2, repeat=args.telomere_repeat) if s2 else None

            if (p1 and p1["branch"] == "telo") or (p2 and p2["branch"] == "telo"):
                telo_counts[a_id] += 1

            if p1 and p1["branch"] == "germline":
                germline_extensions[a_id].append(p1["germline_ext"])
                local_reads[a_id].append(s1)
                if s2:
                    local_reads[a_id].append(s2)
            elif p2 and p2["branch"] == "germline":
                germline_extensions[a_id].append(p2["germline_ext"])
                local_reads[a_id].append(s2)
                local_reads[a_id].append(s1)

    loci = []
    for a_id, a_seq in anchors.items():
        c_telo = telo_counts[a_id]
        c_germ = len(germline_extensions[a_id])
        germ_cons = ""
        coherence = 0.0
        mh_len = 0
        mh_motif = ""
        break_window = ""

        if c_germ >= args.min_germ_reads:
            germ_cons, coherence = reconstruct_germline_consensus(
                germline_extensions[a_id],
                max_len=150,
                anchor_seq=a_seq,
                local_reads=local_reads.get(a_id, []),
                threads=args.threads,
            )
            mh_len, mh_motif = annotate_microhomology(
                germ_cons, repeat=args.telomere_repeat
            )
            break_window = build_300bp_break_window(a_seq, germ_cons)

        loci.append({
            "anchor_id": a_id,
            "anchor_seq": a_seq,
            "c_telo": c_telo,
            "c_germ": c_germ,
            "coherence": coherence,
            "microhomology_len": mh_len,
            "microhomology_motif": mh_motif,
            "germline_consensus": germ_cons,
            "break_window": break_window,
        })

    # Deduplicate allelic break windows via vsearch
    deduped_loci = deduplicate_break_windows(
        loci, min_identity=0.95, threads=args.threads
    )

    # Deconvolve into classes and compute census
    classified_loci, census = deconvolve_and_census(
        deduped_loci,
        min_germ_reads=args.min_germ_reads,
        min_coherence=args.min_coherence,
    )

    # Write junction catalog TSV
    fieldnames = [
        "anchor_id",
        "classification",
        "c_telo",
        "c_germ",
        "coherence",
        "microhomology_len",
        "microhomology_motif",
        "soma_flank_len",
        "germ_flank_len",
    ]
    with open(args.output_catalog, "w", newline="", encoding="utf-8") as f_tsv:
        writer = csv.DictWriter(f_tsv, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for loc in classified_loci:
            writer.writerow({
                "anchor_id": loc["anchor_id"],
                "classification": loc["classification"],
                "c_telo": loc["c_telo"],
                "c_germ": loc["c_germ"],
                "coherence": loc["coherence"],
                "microhomology_len": loc["microhomology_len"],
                "microhomology_motif": loc["microhomology_motif"],
                "soma_flank_len": min(150, len(loc["anchor_seq"])),
                "germ_flank_len": len(loc.get("germline_consensus", "")),
            })

    # Write break_windows_300bp.fasta for PDE_BREAK_SITE loci
    with open(args.output_fasta, "w", encoding="utf-8") as f_fa:
        for loc in classified_loci:
            if (
                loc["classification"] == "PDE_BREAK_SITE"
                and loc.get("break_window")
            ):
                hdr = (
                    f"{loc['anchor_id']} class=PDE_BREAK_SITE "
                    f"c_telo={loc['c_telo']} c_germ={loc['c_germ']} "
                    f"coherence={loc['coherence']} mh_len={loc['microhomology_len']}"
                )
                f_fa.write(f">{hdr}\n{loc['break_window']}\n")

    # Write census JSON
    with open(args.output_census, "w", encoding="utf-8") as f_json:
        json.dump(census, f_json, indent=2)

    print(
        f"[BAITGERMLINEJUNCTIONS] Census: N_breaks={census['n_breaks']}, "
        f"N_germ_telo={census['n_germ_telo']}, C_soma={census['median_c_soma']}, "
        f"C_germ={census['median_c_germ']}, Germline Fraction={census['germline_fraction']}"
    )


if __name__ == "__main__":
    main()
