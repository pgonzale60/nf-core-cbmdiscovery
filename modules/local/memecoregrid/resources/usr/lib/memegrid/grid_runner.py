"""MEME parameter grid search runner and core motif selector."""

from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
from typing import Any

from .evaluator import parse_meme_txt
from .slicer import load_fasta, write_sliced_fasta

# Radii specified in architecture spec
RADII = [30, 50, 100, 150]

# Width intervals specified in architecture spec
WIDTH_RANGES = [
    {"label": "w6_15", "args": ["-minw", "6", "-maxw", "15"]},
    {"label": "w12_30", "args": ["-minw", "12", "-maxw", "30"]},
    {"label": "w15_50", "args": ["-minw", "15", "-maxw", "50"]},
]

MODELS = ["zoops", "anr", "oops"]


def generate_grid_conditions(radii: list[int] | None = None) -> list[dict[str, Any]]:
    """Generate parameter combinations across radius, model, strand, and width ranges."""
    selected_radii = radii or RADII
    conditions = []
    idx = 1
    for r in selected_radii:
        for mod in MODELS:
            for revcomp in [True, False]:
                for w in WIDTH_RANGES:
                    strand_lbl = "revcomp" if revcomp else "fwd"
                    cid = f"R{r}_{mod}_{strand_lbl}_{w['label']}"
                    conditions.append({
                        "id": cid,
                        "cond_idx": idx,
                        "radius": r,
                        "mod": mod,
                        "revcomp": revcomp,
                        "w_args": w["args"],
                        "w_label": w["label"],
                    })
                    idx += 1
    return conditions


def _run_single_condition(
    meme_bin: str | None,
    cond: dict[str, Any],
    r_fa_str: str,
    cond_dir_str: str,
    total_seqs: int,
) -> dict[str, Any] | None:
    """Run MEME for a single condition and parse metrics."""
    cond_dir = Path(cond_dir_str)
    meme_txt = cond_dir / "meme.txt"
    r = cond["radius"]
    cid = cond["id"]

    if not meme_txt.exists():
        if not meme_bin:
            raise RuntimeError("Required binary 'meme' not found in PATH")
        cmd = [
            meme_bin,
            r_fa_str,
            "-dna",
            "-nmotifs",
            "1",
            "-mod",
            cond["mod"],
            "-oc",
            cond_dir_str,
        ] + cond["w_args"]
        if cond["revcomp"]:
            cmd.append("-revcomp")

        subprocess.run(
            cmd,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True,
        )

    if meme_txt.exists():
        metrics = parse_meme_txt(
            str(meme_txt), radius=r, total_sequences=total_seqs
        )
        if metrics:
            metrics["condition_id"] = cid
            metrics["radius"] = r
            metrics["mod"] = cond["mod"]
            metrics["revcomp"] = cond["revcomp"]
            metrics["w_label"] = cond["w_label"]
            metrics["meme_txt_path"] = str(meme_txt)
            return metrics
    return None


def run_meme_grid(
    break_windows_fasta: str | Path,
    output_dir: str | Path,
    radii: list[int] | None = None,
    threads: int = 1,
) -> tuple[dict[str, Any] | None, list[dict[str, Any]]]:
    """Execute MEME parameter grid search across multi-radius sliced windows."""
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    records = load_fasta(break_windows_fasta)
    total_seqs = len(records)
    if total_seqs == 0:
        return None, []

    # 1. Slice FASTA for each radius
    sliced_dir = out_path / "sliced_windows"
    sliced_dir.mkdir(exist_ok=True)
    sliced_fastas = {}
    for r in radii or RADII:
        r_fa = sliced_dir / f"windows_radius_{r}bp.fasta"
        write_sliced_fasta(records, r, r_fa, center_soma_len=150)
        sliced_fastas[r] = r_fa

    # 2. Run MEME across grid conditions in parallel
    conditions = generate_grid_conditions(radii=radii)
    runs_dir = out_path / "runs"
    runs_dir.mkdir(exist_ok=True)

    meme_bin = shutil.which("meme")
    if not meme_bin and not any(runs_dir.glob("*/meme.txt")):
        raise RuntimeError("Required binary 'meme' not found in PATH")
    tasks = []
    for cond in conditions:
        cid = cond["id"]
        cond_dir = runs_dir / cid
        cond_dir.mkdir(exist_ok=True)
        r = cond["radius"]
        r_fa = sliced_fastas[r]
        tasks.append((cond, str(r_fa), str(cond_dir)))

    grid_results = []
    max_workers = max(1, threads)
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(
                _run_single_condition,
                meme_bin,
                cond,
                r_fa_str,
                c_dir_str,
                total_seqs,
            )
            for cond, r_fa_str, c_dir_str in tasks
        ]
        for f in as_completed(futures):
            res = f.result()
            if res:
                grid_results.append(res)

    # 3. Select best core motif maximizing composite score
    if not grid_results:
        return None, []

    grid_results.sort(
        key=lambda x: (x["composite_score"], -x["evalue"]), reverse=True
    )
    best_motif = grid_results[0]

    return best_motif, grid_results


def export_grid_outputs(
    best_motif: dict[str, Any] | None,
    grid_results: list[dict[str, Any]],
    output_meme: str | Path,
    summary_tsv: str | Path,
    metrics_json: str | Path,
):
    """Save best_core_motif.meme, core_grid_summary.tsv, and core_motif_metrics.json."""
    if best_motif and os.path.exists(best_motif["meme_txt_path"]):
        shutil.copy(best_motif["meme_txt_path"], output_meme)
    else:
        # Create empty placeholder if no motif found
        with open(output_meme, "w", encoding="utf-8") as f:
            f.write("MEME version 4\n\nALPHABET= ACGT\n")

    # Summary TSV
    fieldnames = [
        "condition_id",
        "radius",
        "mod",
        "revcomp",
        "w_label",
        "width",
        "nsites",
        "evalue",
        "consensus",
        "total_ic",
        "mean_ic",
        "consensus_entropy",
        "occupancy",
        "median_offset_to_pos0",
        "positional_iqr",
        "composite_score",
    ]
    with open(summary_tsv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(grid_results)

    # Metrics JSON
    with open(metrics_json, "w", encoding="utf-8") as f:
        json.dump(best_motif or {}, f, indent=2)
