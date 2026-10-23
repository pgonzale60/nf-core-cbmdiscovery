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


def _is_short_microsatellite(sub: str) -> bool:
    """Return True if sequence is dominated by a short period 1..9 bp microsatellite."""
    n = len(sub)
    if n < 10:
        return True
    for p in range(1, 10):
        matches = sum(1 for i in range(n - p) if sub[i] == sub[i + p])
        if matches / (n - p) >= 0.80:
            return True
    return False


def is_tandem_satellite_window(
    seq: str, min_period: int = 10, max_period: int = 100, min_copies: int = 3
) -> bool:
    """Detect if a window contains a tandem satellite array of monomer length 10..100 bp repeating >= 3 times."""
    if not seq:
        return False
    s = seq.upper()
    n = len(s)
    if n < min_period * min_copies:
        return False

    upper_p = min(max_period, n // min_copies)
    for p in range(min_period, upper_p + 1):
        span = p * min_copies
        for start in range(0, n - span + 1, 4):
            m1 = s[start : start + p]
            m2 = s[start + p : start + 2 * p]
            m3 = s[start + 2 * p : start + 3 * p]
            id12 = sum(1 for a, b in zip(m1, m2, strict=False) if a == b) / p
            id23 = sum(1 for a, b in zip(m2, m3, strict=False) if a == b) / p
            if (
                id12 >= 0.82
                and id23 >= 0.82
                and not _is_short_microsatellite(m1)
                and not _is_short_microsatellite(m2)
            ):
                return True
    return False


def _vsearch_cluster_pass(
    loci: list[dict[str, Any]],
    seq_getter,
    min_identity: float,
    threads: int,
    vsearch_bin: str,
) -> list[dict[str, Any]]:
    """Run one pass of vsearch --cluster_fast using sequences from seq_getter."""
    if len(loci) <= 1:
        return [dict(l, collapsed_copies=l.get("collapsed_copies", 1)) for l in loci]

    sorted_loci = sorted(
        [dict(l, collapsed_copies=l.get("collapsed_copies", 1)) for l in loci],
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
                w = seq_getter(loc)
                if w and len(w) >= 20:
                    f_fa.write(f">{lid}\n{w}\n")
                    written_count += 1

        if written_count <= 1:
            return sorted_loci

        cmd = [
            vsearch_bin, "--cluster_fast", str(raw_fasta),
            "--id", str(min_identity), "--iddef", "0",
            "--minseqlength", "1", "--qmask", "none",
            "--strand", "both", "--uc", str(uc_file),
            "--threads", str(max(1, threads)),
        ]
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

        centroid_members: dict[str, list[dict[str, Any]]] = {}
        if uc_file.exists():
            with open(uc_file, encoding="utf-8") as f_uc:
                for line in f_uc:
                    parts = line.strip().split("\t")
                    if len(parts) < 10:
                        continue
                    if parts[0] == "S":
                        centroid_members.setdefault(parts[8], []).append(loci_by_id[parts[8]])
                    elif parts[0] == "H":
                        centroid_members.setdefault(parts[9], []).append(loci_by_id[parts[8]])

        seen_loci_ids: set[int] = set()
        deduped: list[dict[str, Any]] = []
        for members in centroid_members.values():
            rep = dict(members[0])
            rep["collapsed_copies"] = sum(m.get("collapsed_copies", 1) for m in members)
            for m in members[1:]:
                rep["c_telo"] = rep.get("c_telo", 0) + m.get("c_telo", 0)
                rep["c_germ"] = rep.get("c_germ", 0) + m.get("c_germ", 0)
                for k in ("anchor_seq", "germline_consensus", "break_window"):
                    if len(m.get(k, "")) > len(rep.get(k, "")):
                        rep[k] = m[k]
            for m in members:
                for idx, loc in enumerate(sorted_loci):
                    if loc is m:
                        seen_loci_ids.add(idx)
            deduped.append(rep)

        for idx, loc in enumerate(sorted_loci):
            if idx not in seen_loci_ids:
                deduped.append(dict(loc))

        return deduped


def deduplicate_break_windows(
    loci: list[dict[str, Any]], min_identity: float = 0.85, threads: int = 1
) -> list[dict[str, Any]]:
    """Collapse break windows sharing >= min_identity across full 300 bp window OR retained somatic flank [:130]."""
    if not loci:
        return []
    if len(loci) == 1:
        res = dict(loci[0])
        res.setdefault("collapsed_copies", 1)
        return [res]

    vsearch_bin = shutil.which("vsearch")
    if not vsearch_bin:
        raise RuntimeError("Required binary 'vsearch' not found in PATH")

    # Pass 1: Full 300 bp break window deduplication
    pass1 = _vsearch_cluster_pass(
        loci,
        seq_getter=lambda l: l.get("break_window") or l.get("anchor_seq", ""),
        min_identity=min_identity,
        threads=threads,
        vsearch_bin=vsearch_bin,
    )

    # Pass 2: Paralogous/satellite retained somatic flank window[:130] deduplication at >= 85% identity
    flank_id = min(min_identity, 0.85)
    return _vsearch_cluster_pass(
        pass1,
        seq_getter=lambda l: (l.get("anchor_seq") or l.get("break_window", "")[:150])[:130],
        min_identity=flank_id,
        threads=threads,
        vsearch_bin=vsearch_bin,
    )


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
        loc.setdefault("collapsed_copies", 1)
        w = loc.get("break_window") or loc.get("anchor_seq", "")
        loc["_is_satellite"] = is_tandem_satellite_window(w)
        c_germ = loc.get("c_germ", 0)
        coherence = loc.get("coherence", 0.0)
        c_telo = loc.get("c_telo", 0)
        if (
            not loc["_is_satellite"]
            and c_telo >= min_telo_reads
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
        is_sat = loc.pop("_is_satellite", False)

        # High-copy repeat / rDNA capped native telomere filter:
        # Repetitive arrays terminating at germline telomeres exhibit massive c_germ (c_germ >> med_c_germ and c_germ >> c_telo)
        is_high_copy = (
            med_c_germ > 0
            and c_telo >= telo_cutoff
            and c_germ >= max(50, round(5.0 * med_c_germ))
            and c_germ >= 5 * c_telo
        )

        is_pde = False
        if not is_sat and not is_high_copy and c_telo >= telo_cutoff and c_germ >= min_germ_reads:
            if coherence >= min_coherence:
                is_pde = True
            elif c_germ >= min(5, min_germ_reads) and coherence >= 0.65:
                is_pde = True

        if is_sat:
            loc["classification"] = "TANDEM_SATELLITE_ARRAY"
        elif is_pde:
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
        if loc["classification"] in (
            "PDE_BREAK_SITE",
            "GERMLINE_NATIVE_TELOMERE",
            "TANDEM_SATELLITE_ARRAY",
        ):
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
    n_satellite = sum(
        1 for l in loci if l["classification"] == "TANDEM_SATELLITE_ARRAY"
    )

    census = {
        "median_c_soma": round(med_c_soma, 2),
        "median_c_germ": round(med_c_germ, 2),
        "germline_fraction": germ_fraction,
        "n_breaks": n_breaks,
        "n_germ_telo": n_germ_telo,
        "estimated_germline_chromosomes": round(n_germ_telo / 2.0, 1),
        "n_shared_telo": n_shared_telo,
        "n_satellite_excluded": n_satellite,
    }

    return loci, census
