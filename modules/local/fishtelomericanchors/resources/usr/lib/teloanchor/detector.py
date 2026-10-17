"""Detector for terminal telomeric arrays in raw reads."""

from typing import Any


def revcomp(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    trans = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(trans)[::-1]


def get_cyclic_permutations(repeat: str) -> list[str]:
    """Return all unique cyclic permutations of a repeat."""
    n = len(repeat)
    doubled = repeat + repeat
    perms = [doubled[i : i + n] for i in range(n)]
    return sorted(set(perms))


def find_terminal_telomeric_tract(
    seq: str,
    repeat: str = "TTAGGC",
    platform: str = "hifi",
    min_repeats: int | None = None,
    tip_tolerance_bp: int | None = None,
    min_flank_bp: int | None = None,
    gap_tolerance_bp: int = 30,
) -> dict[str, Any] | None:
    """Detect terminal telomeric repeat array in a read sequence.

    Returns dict with keys:
        - terminal_end: '5p' or '3p'
        - tract_orientation: 'fwd' (TTAGGC-like) or 'rev' (GCCTAA-like)
        - tract_start: int
        - tract_stop: int
        - tract_len: int
        - non_telo_flank: str (raw orientation in read)
    or None if no terminal tract meets the criteria.
    """
    repeat = repeat.upper()
    rep_len = len(repeat)
    rev_repeat = revcomp(repeat)

    if platform.lower() == "hifi":
        eff_min_repeats = 10 if min_repeats is None else min_repeats
        eff_tip_tol = 100 if tip_tolerance_bp is None else tip_tolerance_bp
        eff_min_flank = 100 if min_flank_bp is None else min_flank_bp
    else:  # illumina
        eff_min_repeats = 4 if min_repeats is None else min_repeats
        eff_tip_tol = 30 if tip_tolerance_bp is None else tip_tolerance_bp
        eff_min_flank = 50 if min_flank_bp is None else min_flank_bp

    min_tract_len = eff_min_repeats * rep_len
    fwd_perms = get_cyclic_permutations(repeat)
    rev_perms = get_cyclic_permutations(rev_repeat)

    hits: list[tuple[int, int, str]] = []
    # Search for forward permutations
    for p in fwd_perms:
        idx = seq.find(p)
        while idx >= 0:
            hits.append((idx, idx + len(p), "fwd"))
            idx = seq.find(p, idx + 1)
    # Search for reverse complement permutations
    for p in rev_perms:
        idx = seq.find(p)
        while idx >= 0:
            hits.append((idx, idx + len(p), "rev"))
            idx = seq.find(p, idx + 1)

    if not hits:
        return None

    hits.sort(key=lambda x: (x[0], x[1]))

    # Merge into runs within gap_tolerance_bp having compatible orientation
    runs: list[dict[str, Any]] = []
    cur_start, cur_stop, cur_orients = hits[0][0], hits[0][1], [hits[0][2]]

    for start, stop, orient in hits[1:]:
        if start - cur_stop <= gap_tolerance_bp:
            cur_stop = max(cur_stop, stop)
            cur_orients.append(orient)
        else:
            fwd_count = cur_orients.count("fwd")
            rev_count = cur_orients.count("rev")
            dom_orient = "fwd" if fwd_count >= rev_count else "rev"
            runs.append({
                "start": cur_start,
                "stop": cur_stop,
                "len": cur_stop - cur_start,
                "orient": dom_orient,
            })
            cur_start, cur_stop, cur_orients = start, stop, [orient]

    fwd_count = cur_orients.count("fwd")
    rev_count = cur_orients.count("rev")
    dom_orient = "fwd" if fwd_count >= rev_count else "rev"
    runs.append({
        "start": cur_start,
        "stop": cur_stop,
        "len": cur_stop - cur_start,
        "orient": dom_orient,
    })

    read_len = len(seq)

    # 1. Check 3' terminal tract (most common for fwd somatic -> telo)
    for r in reversed(runs):
        dist_from_3p = read_len - r["stop"]
        if dist_from_3p <= eff_tip_tol and r["len"] >= min_tract_len:
            flank = seq[: r["start"]]
            if len(flank) >= eff_min_flank:
                return {
                    "terminal_end": "3p",
                    "tract_orientation": r["orient"],
                    "tract_start": r["start"],
                    "tract_stop": r["stop"],
                    "tract_len": r["len"],
                    "non_telo_flank": flank,
                }

    # 2. Check 5' terminal tract
    for r in runs:
        dist_from_5p = r["start"]
        if dist_from_5p <= eff_tip_tol and r["len"] >= min_tract_len:
            flank = seq[r["stop"] :]
            if len(flank) >= eff_min_flank:
                return {
                    "terminal_end": "5p",
                    "tract_orientation": r["orient"],
                    "tract_start": r["start"],
                    "tract_stop": r["stop"],
                    "tract_len": r["len"],
                    "non_telo_flank": flank,
                }

    return None
