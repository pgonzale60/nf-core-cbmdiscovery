#!/usr/bin/env python3
"""Stage 5: Core-Masking, Tier-2 Auxiliary Grid Search & Bipartite Profiler."""

import argparse
from collections import defaultdict
from pathlib import Path
import shutil
import sys

# Ensure usr/lib is in python search path
_LIB_DIR = Path(__file__).resolve().parent.parent / "lib"
if str(_LIB_DIR) not in sys.path:
    sys.path.insert(0, str(_LIB_DIR))

from bipartite.aux_grid import run_auxiliary_grid  # noqa: E402
from bipartite.masker import (  # noqa: E402
    load_fasta_dict,
    mask_core_motifs,
    run_fimo,
)
from bipartite.reporter import generate_cbm_discovery_report  # noqa: E402
from bipartite.spacing import (  # noqa: E402
    profile_bipartite_spacing,
    write_spacing_tsv,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Mask core motif, run auxiliary MEME grid search, profile bipartite spacing, and generate final report"
    )
    parser.add_argument(
        "--break-windows",
        required=True,
        help="Path to break_windows_300bp.fasta",
    )
    parser.add_argument(
        "--core-meme", required=True, help="Path to best_core_motif.meme"
    )
    parser.add_argument(
        "--census-json",
        required=True,
        help="Path to coverage_and_chromosome_census.json",
    )
    parser.add_argument(
        "--core-metrics-json",
        required=True,
        help="Path to core_motif_metrics.json",
    )
    parser.add_argument("--sample-id", default="sample", help="Sample ID")
    parser.add_argument(
        "--output-aux-meme",
        default="best_auxiliary_motif.meme",
        help="Output auxiliary MEME file",
    )
    parser.add_argument(
        "--output-spacing-tsv",
        default="bipartite_spacing_summary.tsv",
        help="Output TSV of bipartite spacing details",
    )
    parser.add_argument(
        "--output-report",
        default="cbm_discovery_report.md",
        help="Output markdown report",
    )
    parser.add_argument(
        "--output-masked-fasta",
        default="break_sequences.core_masked.fasta",
        help="Output core-masked FASTA",
    )
    parser.add_argument(
        "--threads", type=int, default=1, help="Threads for MEME/FIMO"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    raw_seqs = load_fasta_dict(args.break_windows)
    seq_names = list(raw_seqs.keys())

    # 1. Mask core motifs with N
    _masked_seqs, core_hits = mask_core_motifs(
        core_meme=args.core_meme,
        fasta_file=args.break_windows,
        output_masked_fasta=args.output_masked_fasta,
        fimo_out_dir="fimo_core_hits",
        threshold=1e-4,
    )

    core_hits_by_seq = defaultdict(list)
    for h in core_hits:
        core_hits_by_seq[h["seq_name"]].append(h)

    # 2. Auxiliary motif grid search
    best_aux, _grid_results = run_auxiliary_grid(
        masked_fasta=args.output_masked_fasta,
        out_dir="aux_grid_runs",
        threads=args.threads,
    )

    aux_hits_by_seq = defaultdict(list)
    if best_aux and best_aux.get("meme_txt_path"):
        shutil.copy(best_aux["meme_txt_path"], args.output_aux_meme)
        # Run FIMO with auxiliary motif
        aux_fimo_hits = run_fimo(
            meme_file=args.output_aux_meme,
            fasta_file=args.break_windows,
            out_dir="fimo_aux_hits",
            threshold=1e-3,
        )
        for h in aux_fimo_hits:
            aux_hits_by_seq[h["seq_name"]].append(h)
    else:
        # Create empty placeholder if no auxiliary motif detected
        with open(args.output_aux_meme, "w", encoding="utf-8") as f:
            f.write("MEME version 4\n\nALPHABET= ACGT\n")

    # 3. Profile bipartite spacing and cleavage offsets
    spacing_rows, spacing_summary = profile_bipartite_spacing(
        seq_names=seq_names,
        core_hits_by_seq=core_hits_by_seq,
        aux_hits_by_seq=aux_hits_by_seq if best_aux else None,
        pos0_coord=150,
    )
    write_spacing_tsv(spacing_rows, args.output_spacing_tsv)

    # 4. Generate CBM Discovery Markdown Report
    generate_cbm_discovery_report(
        sample_id=args.sample_id,
        census_json_path=args.census_json,
        core_metrics_json_path=args.core_metrics_json,
        best_aux_motif=best_aux,
        bipartite_summary=spacing_summary,
        output_report_md=args.output_report,
    )

    print(
        f"[BIPARTITEPROFILER] Completed for {args.sample_id}. "
        f"Auxiliary motif: {'Found (' + best_aux['consensus'] + ')' if best_aux else 'None'}. "
        f"Paired loci: {spacing_summary['n_paired_occurrences']}."
    )


if __name__ == "__main__":
    main()
