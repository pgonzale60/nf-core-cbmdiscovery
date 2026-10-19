"""Junction classification, deduplication, and coverage/chromosome census."""

from pathlib import Path
import shutil
import statistics
import subprocess
import tempfile
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
    loci: list[dict[str, Any]], min_identity: float = 0.95, threads: int = 1
) -> list[dict[str, Any]]:
    """Deduplicate allelic or identical break windows at >= min_identity via vsearch --cluster_fast."""
    if not loci:
        return []
    if len(loci) == 1:
        return [dict(loci[0])]

    vsearch_bin = shutil.which("vsearch")
    if not vsearch_bin:
        raise RuntimeError("Required binary 'vsearch' not found in PATH")

    sorted_loci = sorted(
        loci,
        key=lambda x: (x.get("c_telo", 0) + x.get("c_germ", 0)),
        reverse=True,
    )

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        raw_fasta = tmp_path / "raw_windows.fasta"
        uc_file = tmp_path / "dedup.uc"

        loci_by_id: dict[str, dict[str, Any]] = {}
        written_count = 0
        with open(raw_fasta, "w", encoding="utf-8") as f_fa:
            for idx, loc in enumerate(sorted_loci):
                lid = f"loc_{idx}"
                loci_by_id[lid] = loc
                w = loc.get("break_window") or loc.get("anchor_seq", "")
                if w:
                    f_fa.write(f">{lid}\n{w}\n")
                    written_count += 1

        if written_count == 0:
            return [dict(loc) for loc in sorted_loci]

        cmd = [
            vsearch_bin,
            "--cluster_fast",
            str(raw_fasta),
            "--id",
            str(min_identity),
            "--iddef",
            "0",
            "--minseqlength",
            "1",
            "--qmask",
            "none",
            "--strand",
            "both",
            "--uc",
            str(uc_file),
            "--threads",
            str(max(1, threads)),
        ]

        subprocess.run(
            cmd,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True,
        )

        centroid_members: dict[str, list[dict[str, Any]]] = {}
        if uc_file.exists():
            with open(uc_file, encoding="utf-8") as f_uc:
                for line in f_uc:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split("\t")
                    rec_type = parts[0]
                    if rec_type == "S":
                        query_id = parts[8]
                        centroid_members.setdefault(query_id, []).append(loci_by_id[query_id])
                    elif rec_type == "H":
                        query_id = parts[8]
                        target_id = parts[9]
                        centroid_members.setdefault(target_id, []).append(loci_by_id[query_id])

        seen_loci_ids: set[int] = set()
        deduped: list[dict[str, Any]] = []
        for members in centroid_members.values():
            rep = dict(members[0])
            for m in members[1:]:
                rep["c_telo"] = rep.get("c_telo", 0) + m.get("c_telo", 0)
                rep["c_germ"] = rep.get("c_germ", 0) + m.get("c_germ", 0)
                if len(m.get("anchor_seq", "")) > len(rep.get("anchor_seq", "")):
                    rep["anchor_seq"] = m["anchor_seq"]
                if len(m.get("germline_consensus", "")) > len(rep.get("germline_consensus", "")):
                    rep["germline_consensus"] = m["germline_consensus"]
                if len(m.get("break_window", "")) > len(rep.get("break_window", "")):
                    rep["break_window"] = m["break_window"]
            for m in members:
                for idx, loc in enumerate(sorted_loci):
                    if loc is m:
                        seen_loci_ids.add(idx)
            deduped.append(rep)

        for idx, loc in enumerate(sorted_loci):
            if idx not in seen_loci_ids:
                deduped.append(dict(loc))

        return deduped


def deconvolve_and_census(
    loci: list[dict[str, Any]],
    min_germ_reads: int = 2,
    min_coherence: float = 0.85,
    min_telo_reads: int = 5,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Classify anchors into PDE_BREAK_SITE, GERMLINE_NATIVE_TELOMERE, SHARED_NATIVE_TELOMERE,

    UNBROKEN_GERMLINE_LOCUS, and compute coverage and chromosome census metrics.
    """
    # 1. First pass: find candidate break sites with initial dual-branch support
    cand_breaks = []
    for loc in loci:
        c_germ = loc.get("c_germ", 0)
        coherence = loc.get("coherence", 0.0)
        c_telo = loc.get("c_telo", 0)
        if (
            c_telo >= min_telo_reads
            and c_germ >= min_germ_reads
            and coherence >= 0.65
            and c_germ <= 5 * max(c_telo, min_telo_reads)
        ):
            cand_breaks.append(loc)

    if cand_breaks:
        c_soma_vals = [l["c_telo"] for l in cand_breaks]
        c_germ_vals = [l["c_germ"] for l in cand_breaks]
        med_c_soma = float(statistics.median(c_soma_vals))
        med_c_germ = float(statistics.median(c_germ_vals))
    else:
        med_c_soma = 0.0
        med_c_germ = 0.0

    # Adaptive telomeric bifurcation threshold based on median somatic coverage of candidate break sites
    telo_cutoff = (
        max(min_telo_reads, round(0.15 * med_c_soma))
        if med_c_soma > 0
        else min_telo_reads
    )

    for loc in loci:
        c_germ = loc.get("c_germ", 0)
        coherence = loc.get("coherence", 0.0)
        c_telo = loc.get("c_telo", 0)

        # High-copy repeat / rDNA capped native telomere filter:
        # Repetitive arrays terminating at germline telomeres exhibit massive c_germ (c_germ >> med_c_germ and c_germ >> c_telo)
        is_high_copy = (
            med_c_germ > 0
            and c_telo >= telo_cutoff
            and c_germ >= max(50, round(5.0 * med_c_germ))
            and c_germ >= 5 * c_telo
        )

        is_pde = False
        if not is_high_copy and c_telo >= telo_cutoff and c_germ >= min_germ_reads:
            if coherence >= min_coherence:
                is_pde = True
            elif c_germ >= min(5, min_germ_reads) and coherence >= 0.65:
                is_pde = True

        if is_pde:
            loc["classification"] = "PDE_BREAK_SITE"
        elif is_high_copy:
            loc["classification"] = "GERMLINE_NATIVE_TELOMERE"
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

    # 2. Second pass: classify non-PDE loci
    cutoff = (
        med_c_germ * 1.5
        if med_c_germ > 0
        else (med_c_soma * 0.4 if med_c_soma > 0 else 5.0)
    )

    for loc in loci:
        if loc["classification"] in ("PDE_BREAK_SITE", "GERMLINE_NATIVE_TELOMERE"):
            continue
        c_telo = loc.get("c_telo", 0)
        c_germ = loc.get("c_germ", 0)

        if c_telo < telo_cutoff and c_germ >= min_germ_reads:
            loc["classification"] = "UNBROKEN_GERMLINE_LOCUS"
        elif c_germ >= min_germ_reads:
            loc["classification"] = "LOW_COHERENCE_JUNCTION"
        elif c_telo <= cutoff:
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
