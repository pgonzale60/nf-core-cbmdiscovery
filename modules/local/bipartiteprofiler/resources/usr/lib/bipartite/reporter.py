"""Markdown report generator for CBM Discovery pipeline."""

import json
import os
from pathlib import Path
from typing import Any


def classify_pde_regime(
    n_unique_breaks: int,
    core_metrics: dict[str, Any] | None = None,
    read_break_iqr: float = 0.0,
) -> tuple[str, str]:
    """Classify sample into universal DSB-PDE regime:

    - CLASS_I_PUNCTATE_MOTIF_PDE
    - CLASS_II_ZONAL_CBR_PDE
    - CLASS_III_HIGH_DENSITY_FRAGMENTATION_PDE
    - PDE_NEGATIVE
    Returns (regime_name, rationale_description).
    """
    if n_unique_breaks >= 100:
        return (
            "CLASS_III_HIGH_DENSITY_FRAGMENTATION_PDE",
            f"High-Density Micro-Chromosome Fragmentation PDE ({n_unique_breaks} breaks >= 100 with terminal microhomology).",
        )
    if n_unique_breaks < 4:
        return (
            "PDE_NEGATIVE",
            f"Negative for Programmed DNA Elimination ({n_unique_breaks} breaks < 4 after filtering).",
        )

    # For 4 <= n_breaks < 100, check Tier-1 Core Motif significance and spatial dispersion
    if core_metrics:
        evalue_val = core_metrics.get("evalue")
        if evalue_val is None:
            estr = core_metrics.get("evalue_str")
            if estr:
                try:
                    evalue_val = float(estr)
                except ValueError:
                    evalue_val = 1.0
            else:
                evalue_val = 1.0
        try:
            evalue = float(evalue_val)
        except (ValueError, TypeError):
            evalue = 1.0
        rejection_flag = core_metrics.get("rejection_flag", "PASS")
        try:
            pos_iqr = float(core_metrics.get("positional_iqr", 100.0))
        except (ValueError, TypeError):
            pos_iqr = 100.0
        is_pass = rejection_flag == "PASS"
        if is_pass and evalue < 0.01 and pos_iqr <= 60.0 and read_break_iqr < 50.0:
            return (
                "CLASS_I_PUNCTATE_MOTIF_PDE",
                f"Punctate CBM Programmed DNA Elimination with localized cleavage motif ({n_unique_breaks} breaks >= 4, E = {evalue:.2e} < 0.01, positional IQR = {pos_iqr:.1f} bp <= 60 bp).",
            )

    return (
        "CLASS_II_ZONAL_CBR_PDE",
        f"Zonal Chromosomal Breakage Region (CBR) Programmed DNA Elimination ({n_unique_breaks} breaks >= 4 without localized point motif).",
    )


def generate_cbm_discovery_report(
    sample_id: str,
    census_json_path: str | Path,
    core_metrics_json_path: str | Path,
    best_aux_motif: dict[str, Any] | None,
    bipartite_summary: dict[str, Any],
    output_report_md: str | Path,
):
    """Generate comprehensive cbm_discovery_report.md markdown summary."""
    census = {}
    if os.path.exists(census_json_path):
        with open(census_json_path, encoding="utf-8") as f:
            census = json.load(f)

    core_metrics = {}
    if os.path.exists(core_metrics_json_path):
        with open(core_metrics_json_path, encoding="utf-8") as f:
            core_metrics = json.load(f)

    n_breaks = int(census.get("n_breaks", 0))
    pde_regime, pde_rationale = classify_pde_regime(n_breaks, core_metrics)

    if core_metrics:
        core_metrics["pde_regime"] = pde_regime
        core_metrics["pde_regime_rationale"] = pde_rationale
        try:
            with open(core_metrics_json_path, "w", encoding="utf-8") as f_cm:
                json.dump(core_metrics, f_cm, indent=2)
        except OSError:
            pass

    md = [
        f"# CBM Discovery & PDE Junction Report: `{sample_id}`\n",
        "## 1. Executive Summary\n",
        f"- **Sample ID**: `{sample_id}`",
        f"- **Universal DSB-PDE Regime**: `{pde_regime}`",
        f"- **Regime Diagnosis**: {pde_rationale}",
        f"- **PDE Breakpoints Identified (N_breaks)**: `{census.get('n_breaks', 'N/A')}`",
        f"- **Germline Native Telomeres (N_germ_telo)**: `{census.get('n_germ_telo', 'N/A')}`",
        f"- **Estimated Germline Chromosome Count**: `{census.get('estimated_germline_chromosomes', 'N/A')}`",
        f"- **Median Somatic Coverage (C_soma)**: `{census.get('median_c_soma', 'N/A')}×`",
        f"- **Median Germline Coverage (C_germ)**: `{census.get('median_c_germ', 'N/A')}×`",
        f"- **Estimated Germline Fraction**: `{round(float(census.get('germline_fraction', 0)) * 100, 1)}%`",
        "\n---\n",
        "## 2. Tier-1 Core Motif Architecture\n",
    ]

    if core_metrics:
        md.extend([
            f"- **Core Consensus**: `{core_metrics.get('consensus', 'N/A')}`",
            f"- **Motif Width**: `{core_metrics.get('width', 'N/A')} bp`",
            f"- **Mean Information Content**: `{core_metrics.get('mean_ic', 'N/A')} bits/pos`",
            f"- **E-value**: `{core_metrics.get('evalue_str', core_metrics.get('evalue', 'N/A'))}`",
            f"- **Break Site Occupancy**: `{round(float(core_metrics.get('occupancy', 0)) * 100, 1)}%` ({core_metrics.get('unique_occupied_sites', 0)} sites)",
            f"- **Median Cleavage Offset (to Pos 0)**: `{core_metrics.get('median_offset_to_pos0', 'N/A')} bp` (IQR: `{core_metrics.get('positional_iqr', 'N/A')} bp`)",
            f"- **Grid Search Condition**: `{core_metrics.get('condition_id', 'N/A')}`",
        ])
    else:
        md.append("- No statistically significant core motif detected.\n")

    md.extend(["\n---\n", "## 3. Tier-2 Auxiliary & Bipartite Architecture\n"])

    if bipartite_summary.get("architecture") == "HOMOTYPIC_BIPARTITE" or (
        best_aux_motif and best_aux_motif.get("architecture") == "HOMOTYPIC_BIPARTITE"
    ):
        md.extend([
            "### Homotypic Bipartite Architecture Detected\n",
            "- **Architecture Type**: `HOMOTYPIC_BIPARTITE`\n",
            f"- **Primary Core Consensus**: `{core_metrics.get('consensus', 'N/A')}`",
            "- **Partner Type**: Dual core motif instances (Core–Core pairing)",
            "\n### Bipartite Spacing & Orientation Metrics\n",
            f"- **Paired Core-Core Loci**: `{bipartite_summary.get('n_paired_occurrences', 0)}`",
            f"- **Median Edge-to-Edge Spacer**: `{bipartite_summary.get('median_spacer_bp', 0)} bp` (IQR: `{bipartite_summary.get('iqr_spacer_bp', 0)} bp`)",
            f"- **Embedded Occurrences (0 bp)**: `{bipartite_summary.get('pct_embedded_0bp', 0)}%`",
            f"- **Tightly Spaced (≤ 30 bp)**: `{bipartite_summary.get('pct_spacer_le30bp', 0)}%`",
            f"- **Flanking (≤ 50 bp)**: `{bipartite_summary.get('pct_spacer_le50bp', 0)}%`",
            f"- **Strand Concordance (Same Strand)**: `{bipartite_summary.get('pct_same_strand', 0)}%`",
            f"- **Located in Retained Region (≤ Pos 0)**: `{bipartite_summary.get('pct_retained_domain', 0.0)}%`",
            f"- **Located in Eliminated Region (> Pos 0)**: `{bipartite_summary.get('pct_eliminated_domain', 0.0)}%`",
        ])
        if bipartite_summary.get("pct_spanning_domain", 0.0) > 0:
            md.append(f"- **Spanning Cleavage Junction**: `{bipartite_summary.get('pct_spanning_domain')}%`")
    elif best_aux_motif:
        md.extend([
            "### Tier-2 Auxiliary Motif Detected\n",
            f"- **Auxiliary Consensus**: `{best_aux_motif.get('consensus', 'N/A')}`",
            f"- **Width**: `{best_aux_motif.get('width', 'N/A')} bp`",
            f"- **E-value**: `{best_aux_motif.get('evalue_str', best_aux_motif.get('evalue', 'N/A'))}`",
            f"- **Mean Information Content**: `{best_aux_motif.get('mean_ic', 'N/A')} bits/pos`",
            f"- **Proximal Sites (6–50 bp from cleavage)**: `{best_aux_motif.get('proximal_sites_6_50bp', 0)}`",
            "\n### Bipartite Spacing & Orientation Metrics\n",
            f"- **Paired Core-Auxiliary Loci**: `{bipartite_summary.get('n_paired_occurrences', 0)}`",
            f"- **Median Edge-to-Edge Spacer**: `{bipartite_summary.get('median_spacer_bp', 0)} bp` (IQR: `{bipartite_summary.get('iqr_spacer_bp', 0)} bp`)",
            f"- **Embedded Occurrences (0 bp)**: `{bipartite_summary.get('pct_embedded_0bp', 0)}%`",
            f"- **Tightly Spaced (≤ 30 bp)**: `{bipartite_summary.get('pct_spacer_le30bp', 0)}%`",
            f"- **Flanking (≤ 50 bp)**: `{bipartite_summary.get('pct_spacer_le50bp', 0)}%`",
            f"- **Strand Concordance (Same Strand)**: `{bipartite_summary.get('pct_same_strand', 0)}%`",
            f"- **Located in Retained Region (≤ Pos 0)**: `{bipartite_summary.get('pct_retained_domain', 0.0)}%`",
            f"- **Located in Eliminated Region (> Pos 0)**: `{bipartite_summary.get('pct_eliminated_domain', 0.0)}%`",
        ])
        if bipartite_summary.get("pct_spanning_domain", 0.0) > 0:
            md.append(f"- **Spanning Cleavage Junction**: `{bipartite_summary.get('pct_spanning_domain')}%`")
    else:
        md.extend([
            "### Simple (Single Motif) Architecture\n",
            "- No statistically significant auxiliary motif was detected after core masking.",
            "- The cleavage system operates primarily via the single core recognition element.",
        ])

    md.append("\n---\n")

    with open(output_report_md, "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
