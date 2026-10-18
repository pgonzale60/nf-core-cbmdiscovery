"""Junction classification, deduplication, and coverage/chromosome census."""

import statistics
from typing import Any


def sequence_identity(s1: str, s2: str) -> float:
    """Calculate sequence identity between two equal-length or overlapping sequences."""
    l_min = min(len(s1), len(s2))
    if l_min == 0:
        return 0.0
    matches = sum(1 for a, b in zip(s1[:l_min], s2[:l_min], strict=False) if a == b)
    return matches / l_min


def are_loci_duplicates(
    loc1: dict[str, Any], loc2: dict[str, Any], min_identity: float = 0.95
) -> bool:
    """Check if two loci represent the same break site anchored at Position 0."""
    soma1 = loc1.get("anchor_seq", "")
    soma2 = loc2.get("anchor_seq", "")
    germ1 = loc1.get("germline_consensus", "")
    germ2 = loc2.get("germline_consensus", "")

    if soma1 and soma2:
        # Compare somatic anchors anchored at their 3' ends (Position 0)
        l_soma = min(len(soma1), len(soma2))
        m_soma = (
            sum(
                1
                for a, b in zip(soma1[-l_soma:], soma2[-l_soma:], strict=False)
                if a == b
            )
            if l_soma > 0
            else 0
        )
        id_soma = m_soma / l_soma if l_soma > 0 else 0.0

        if germ1 and germ2:
            # Compare germline consensus anchored at their 5' ends (Position +1)
            l_germ = min(len(germ1), len(germ2))
            m_germ = (
                sum(
                    1
                    for a, b in zip(germ1[:l_germ], germ2[:l_germ], strict=False)
                    if a == b
                )
                if l_germ > 0
                else 0
            )
            id_germ = m_germ / l_germ if l_germ > 0 else 0.0
            return (id_soma >= min_identity) and (id_germ >= min_identity)
        elif not germ1 and not germ2:
            return id_soma >= min_identity
        else:
            return False

    w1 = loc1.get("break_window", "")
    w2 = loc2.get("break_window", "")
    if w1 and w2:
        return sequence_identity(w1, w2) >= min_identity
    return False


def deduplicate_break_windows(
    loci: list[dict[str, Any]], min_identity: float = 0.95
) -> list[dict[str, Any]]:
    """Deduplicate allelic or identical break windows at >= min_identity anchored at Position 0."""
    # Sort loci by total read support descending
    sorted_loci = sorted(
        loci,
        key=lambda x: (x.get("c_telo", 0) + x.get("c_germ", 0)),
        reverse=True,
    )

    deduped: list[dict[str, Any]] = []
    for loc in sorted_loci:
        duplicate = False
        for kept in deduped:
            if are_loci_duplicates(loc, kept, min_identity=min_identity):
                duplicate = True
                # Merge read counts into the kept representative
                kept["c_telo"] = kept.get("c_telo", 0) + loc.get("c_telo", 0)
                kept["c_germ"] = kept.get("c_germ", 0) + loc.get("c_germ", 0)
                # Preserve the longer sequences
                if len(loc.get("anchor_seq", "")) > len(kept.get("anchor_seq", "")):
                    kept["anchor_seq"] = loc["anchor_seq"]
                if len(loc.get("germline_consensus", "")) > len(
                    kept.get("germline_consensus", "")
                ):
                    kept["germline_consensus"] = loc["germline_consensus"]
                if len(loc.get("break_window", "")) > len(kept.get("break_window", "")):
                    kept["break_window"] = loc["break_window"]
                break

        if not duplicate:
            deduped.append(dict(loc))

    return deduped


def deconvolve_and_census(
    loci: list[dict[str, Any]],
    min_germ_reads: int = 2,
    min_coherence: float = 0.85,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Classify anchors into PDE_BREAK_SITE, GERMLINE_NATIVE_TELOMERE, SHARED_NATIVE_TELOMERE,

    and compute coverage and chromosome census metrics.
    """
    # 1. First pass: identify initial PDE break sites passing germline support and coherence
    for loc in loci:
        c_germ = loc.get("c_germ", 0)
        coherence = loc.get("coherence", 0.0)
        c_telo = loc.get("c_telo", 0)

        if c_germ >= min_germ_reads and coherence >= min_coherence and c_telo > 0:
            loc["classification"] = "PDE_BREAK_SITE"
        else:
            loc["classification"] = "CANDIDATE_TELOMERE"

    pde_breaks = [l for l in loci if l["classification"] == "PDE_BREAK_SITE"]

    if pde_breaks:
        c_soma_vals = [l["c_telo"] for l in pde_breaks]
        c_germ_vals = [l["c_germ"] for l in pde_breaks]
        med_c_soma = float(statistics.median(c_soma_vals))
        med_c_germ = float(statistics.median(c_germ_vals))
    else:
        med_c_soma = 0.0
        med_c_germ = 0.0

    total_cov = med_c_soma + med_c_germ
    germ_fraction = (
        round(med_c_germ / total_cov, 4) if total_cov > 0 else 0.0
    )

    # 2. Second pass: classify telomeric loci without read-through (C_germ < min_germ_reads)
    cutoff = (
        med_c_germ * 1.5
        if med_c_germ > 0
        else (med_c_soma * 0.4 if med_c_soma > 0 else 5.0)
    )

    for loc in loci:
        if loc["classification"] == "PDE_BREAK_SITE":
            continue
        c_telo = loc.get("c_telo", 0)
        if c_telo <= cutoff:
            loc["classification"] = "GERMLINE_NATIVE_TELOMERE"
        else:
            loc["classification"] = "SHARED_NATIVE_TELOMERE"

    n_breaks = sum(1 for l in loci if l["classification"] == "PDE_BREAK_SITE")
    n_germ_telo = sum(
        1 for l in loci if l["classification"] == "GERMLINE_NATIVE_TELOMERE"
    )
    n_shared_telo = sum(
        1 for l in loci if l["classification"] == "SHARED_NATIVE_TELOMERE"
    )

    census = {
        "median_c_soma": round(med_c_soma, 2),
        "median_c_germ": round(med_c_germ, 2),
        "germline_fraction": germ_fraction,
        "n_breaks": n_breaks,
        "n_germ_telo": n_germ_telo,
        "estimated_germline_chromosomes": round(n_germ_telo / 2.0, 1),
        "n_shared_telo": n_shared_telo,
    }

    return loci, census
