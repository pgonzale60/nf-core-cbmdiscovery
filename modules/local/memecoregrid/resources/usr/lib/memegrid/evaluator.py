"""MEME output parser, Information Content calculator, and motif scoring evaluator."""

from collections import Counter
import math
import os
from pathlib import Path
import re
from typing import Any

from cbm_defaults import SSOT


def calculate_consensus_entropy(consensus: str) -> float:
    """Calculate Shannon entropy (in bits) across the characters of a consensus sequence."""
    if not consensus:
        return 0.0
    cons_upper = consensus.upper()
    n = len(cons_upper)
    counts = Counter(cons_upper)
    entropy = 0.0
    for count in counts.values():
        p = count / n
        if p > 0:
            entropy -= p * math.log2(p)
    return round(entropy, 4)


def _revcomp(seq: str) -> str:
    tr = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(tr)[::-1]


def is_telomeric_repeat_echo(
    consensus: str, repeat: str = SSOT["telomere_repeat"]
) -> tuple[bool, float]:
    """Check if motif consensus is a telomeric repeat echo (>=50% exact 6-mers or >=70% in-phase 1-mismatch coverage)."""
    if not consensus or not repeat:
        return False, 0.0
    cons = consensus.upper().strip()
    rep = repeat.upper().strip()
    k = len(rep)
    n = len(cons)
    if n < k:
        return False, 0.0

    rc_rep = _revcomp(rep)
    fwd_perms = {(rep * 2)[i : i + k] for i in range(k)}
    rev_perms = {(rc_rep * 2)[i : i + k] for i in range(k)}

    n_kmers = n - k + 1
    fwd_exact = sum(1 for i in range(n_kmers) if cons[i : i + k] in fwd_perms)
    rev_exact = sum(1 for i in range(n_kmers) if cons[i : i + k] in rev_perms)
    exact_frac = max(fwd_exact, rev_exact) / n_kmers

    # Check in-phase 1-mismatch coverage across all cyclic phases of forward and reverse repeat
    best_cov_1mm = 0.0
    reps_needed = (n // k) + 3
    for base_rep in (rep, rc_rep):
        full_periodic = base_rep * reps_needed
        for phase in range(k):
            tmpl = full_periodic[phase : phase + n]
            covered = [False] * n
            mm1_count = 0
            for i in range(n_kmers):
                mm = sum(
                    1
                    for a, b in zip(cons[i : i + k], tmpl[i : i + k], strict=False)
                    if a != b
                )
                if mm <= 1:
                    mm1_count += 1
                    for j in range(i, i + k):
                        covered[j] = True
            if mm1_count >= 2:
                best_cov_1mm = max(best_cov_1mm, sum(covered) / n)

    is_echo = (exact_frac > SSOT["telomeric_echo_exact_kmer_cutoff"]) or (
        n >= 18 and best_cov_1mm >= SSOT["telomeric_echo_1mm_cov_cutoff"]
    )
    return is_echo, round(max(exact_frac, best_cov_1mm), 4)


def parse_meme_txt(
    txt_path: str | Path,
    radius: int,
    total_sequences: int = 1,
    telomere_repeat: str = SSOT["telomere_repeat"],
) -> dict[str, Any] | None:
    """Parse MEME text output and compute comprehensive motif metrics."""
    if not os.path.exists(txt_path):
        return None

    with open(txt_path, encoding="utf-8") as f:
        text = f.read()

    # Unwrap backslash-wrapped continuation lines if any
    unwrapped = re.sub(r"\\\s*\n\s*", "", text)

    # 1. Header regex for width, sites, E-value
    m = re.search(
        r"letter-probability matrix:\s+alength=\s*(\d+)\s+w=\s*(\d+)\s+nsites=\s*(\d+)\s+E=\s*(\S+)",
        unwrapped,
    )
    if not m:
        return None

    w = int(m.group(2))
    nsites = int(m.group(3))
    evalue_str = m.group(4)
    try:
        evalue = float(evalue_str)
    except ValueError:
        evalue = 1.0

    # 2. Consensus sequence
    cons_m = re.search(r"Multilevel\s+([A-Za-z]+)", unwrapped)
    consensus = cons_m.group(1).upper() if cons_m else ""

    # 3. Parse letter-probability matrix to compute Information Content
    idx = unwrapped.find(m.group(0)) + len(m.group(0))
    matrix_lines = unwrapped[idx:].strip().splitlines()
    ic_per_col = []
    for line in matrix_lines:
        line = line.strip()
        if not line or line.startswith("URL") or line.startswith("-"):
            break
        parts = [float(x) for x in line.split()]
        if len(parts) == 4:
            # IC = 2 - sum(-p * log2(p))
            h = sum(-p * math.log2(p) for p in parts if p > 0)
            ic_per_col.append(max(0.0, 2.0 - h))
        if len(ic_per_col) == w:
            break

    total_ic = sum(ic_per_col)
    mean_ic = total_ic / w if w > 0 else 0.0

    # 4. Parse individual site positions relative to Position 0
    # In sliced window of radius R, Position 0 is at 1-based index R
    in_sites = False
    site_records = []
    for line in unwrapped.splitlines():
        if "sites sorted by position p-value" in line:
            in_sites = True
            continue
        if in_sites:
            if line.startswith("---"):
                continue
            if "block diagrams" in line or line.startswith("MOTIF"):
                break
            parts = line.split()
            if len(parts) >= 4 and ("+" in parts or "-" in parts):
                seq_name = parts[0]
                strand = parts[1]
                try:
                    start_pos = int(parts[2])
                    # Distance relative to Position 0 (dynamic pos0 from sequence header, fallback to radius)
                    p0_match = re.search(r"_p0_(\d+)", seq_name)
                    pos0 = int(p0_match.group(1)) if p0_match else radius
                    # Distance from motif center to Position 0
                    motif_center = start_pos + (w - 1) / 2.0
                    dist_to_pos0 = motif_center - pos0
                    site_records.append({
                        "seq": seq_name,
                        "strand": strand,
                        "start": start_pos,
                        "dist_to_pos0": dist_to_pos0,
                    })
                except ValueError:
                    continue

    unique_seqs = len({s["seq"] for s in site_records})
    occupancy = unique_seqs / total_sequences if total_sequences > 0 else 0.0

    # 5. Positional IQR
    dists = sorted(s["dist_to_pos0"] for s in site_records)
    if len(dists) >= 4:
        q25 = dists[len(dists) // 4]
        q75 = dists[(3 * len(dists)) // 4]
        pos_iqr = float(q75 - q25)
        median_offset = float(dists[len(dists) // 2])
    elif dists:
        pos_iqr = float(max(dists) - min(dists))
        median_offset = float(dists[0])
    else:
        pos_iqr = 50.0
        median_offset = 0.0

    # 6. Homopolymer penalty & Composite Score
    cons_entropy = calculate_consensus_entropy(consensus)
    homopolymer_factor = (
        1.0
        if cons_entropy >= 1.2
        else (cons_entropy / 1.2 if cons_entropy >= 0.6 else 0.05)
    )

    # Score rewards mean IC, occupancy, and tight positional clustering normalized by search radius
    norm_iqr = pos_iqr / max(10, radius)
    iqr_penalty = 1.0 / (1.0 + math.log1p(norm_iqr * 5.0))
    is_telo_echo, telo_echo_frac = is_telomeric_repeat_echo(
        consensus, repeat=telomere_repeat
    )
    if is_telo_echo:
        composite_score = 0.0
        rejection_flag = "REJECTED_TELOMERIC_ECHO"
    elif evalue > SSOT["core_evalue_cutoff"]:
        composite_score = 0.0
        rejection_flag = "REJECTED_INSIGNIFICANT_EVALUE"
    else:
        composite_score = (
            (mean_ic * math.sqrt(w)) * occupancy * iqr_penalty * homopolymer_factor
        )
        rejection_flag = "PASS"

    return {
        "width": w,
        "nsites": nsites,
        "evalue": evalue,
        "evalue_str": evalue_str,
        "consensus": consensus,
        "total_ic": round(total_ic, 2),
        "mean_ic": round(mean_ic, 3),
        "consensus_entropy": cons_entropy,
        "unique_occupied_sites": unique_seqs,
        "occupancy": round(occupancy, 4),
        "median_offset_to_pos0": round(median_offset, 2),
        "positional_iqr": round(pos_iqr, 2),
        "telomeric_echo_fraction": telo_echo_frac,
        "rejection_flag": rejection_flag,
        "composite_score": round(composite_score, 4),
    }
