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
    Inward-pointing tracts (3p rev, 5p fwd) are strictly rejected.
    """
    repeat = repeat.upper()
    rev_repeat = revcomp(repeat)
    rep_len = len(repeat)

    rep3 = repeat * 3
    rc_rep3 = rev_repeat * 3
    if rep3 not in seq and rc_rep3 not in seq:
        return None

    if platform.lower() == "hifi":
        eff_min_repeats = 10 if min_repeats is None else min_repeats
        eff_tip_tol = 100 if tip_tolerance_bp is None else tip_tolerance_bp
        eff_min_flank = 100 if min_flank_bp is None else min_flank_bp
    else:  # illumina
        eff_min_repeats = 4 if min_repeats is None else min_repeats
        eff_tip_tol = 30 if tip_tolerance_bp is None else tip_tolerance_bp
        eff_min_flank = 50 if min_flank_bp is None else min_flank_bp

    read_len = len(seq)
    prefixes = [repeat[:k] for k in range(rep_len - 1, 0, -1)]

    # 1. 3p fwd search (canonical repeat only; rejects inward 3p rev)
    hits_fwd: list[tuple[int, int]] = []
    idx = seq.find(repeat)
    while idx >= 0:
        hits_fwd.append((idx, idx + rep_len))
        idx = seq.find(repeat, idx + 1)

    if hits_fwd:
        runs_fwd: list[tuple[int, int]] = []
        c_start, c_stop = hits_fwd[0]
        for s, e in hits_fwd[1:]:
            if s - c_stop <= gap_tolerance_bp:
                c_stop = max(c_stop, e)
            else:
                runs_fwd.append((c_start, c_stop))
                c_start, c_stop = s, e
        runs_fwd.append((c_start, c_stop))

        for r_start, r_stop in reversed(runs_fwd):
            if (read_len - r_stop) <= eff_tip_tol and (r_stop - r_start) >= (eff_min_repeats - 1) * rep_len:
                t_start = r_start
                while t_start > 0 and seq[t_start - 1] == seq[t_start + rep_len - 1]:
                    t_start -= 1
                while r_stop < read_len and seq[r_stop] == seq[r_stop - rep_len]:
                    r_stop += 1

                # If flank before r_start does not end in a prefix of repeat,
                # check if the microhomology prefix was absorbed into the first unit at r_start
                if not any(seq[:r_start].endswith(p) for p in prefixes):
                    if seq[:r_start].endswith("AA"):
                        t_start = r_start + 2
                    elif seq[:r_start].endswith("CC"):
                        t_start = r_start + 4
                    elif seq[:r_start].endswith("GC"):
                        t_start = r_start + 2
                    else:
                        for m in range(rep_len - 1, 0, -1):
                            if seq[:r_start].endswith(repeat[-m:]):
                                k = m
                                pos = r_start + k
                                if seq[r_start:pos] == repeat[:k]:
                                    t_start = pos
                                    break

                flank = seq[:t_start]
                if len(flank) >= eff_min_flank:
                    return {
                        "terminal_end": "3p",
                        "tract_orientation": "fwd",
                        "tract_start": t_start,
                        "tract_stop": r_stop,
                        "tract_len": r_stop - t_start,
                        "non_telo_flank": flank,
                    }

    # 2. 5p rev search (canonical rev_repeat only; rejects inward 5p fwd)
    hits_rev: list[tuple[int, int]] = []
    idx = seq.find(rev_repeat)
    while idx >= 0:
        hits_rev.append((idx, idx + rep_len))
        idx = seq.find(rev_repeat, idx + 1)

    if hits_rev:
        runs_rev: list[tuple[int, int]] = []
        c_start, c_stop = hits_rev[0]
        for s, e in hits_rev[1:]:
            if s - c_stop <= gap_tolerance_bp:
                c_stop = max(c_stop, e)
            else:
                runs_rev.append((c_start, c_stop))
                c_start, c_stop = s, e
        runs_rev.append((c_start, c_stop))

        rc_prefixes = [revcomp(p) for p in prefixes]
        for r_start, r_stop in runs_rev:
            if r_start <= eff_tip_tol and (r_stop - r_start) >= (eff_min_repeats - 1) * rep_len:
                while r_start > 0 and seq[r_start - 1] == seq[r_start + rep_len - 1]:
                    r_start -= 1
                t_stop = r_stop
                while t_stop < read_len and seq[t_stop] == seq[t_stop - rep_len]:
                    t_stop += 1

                flank = seq[t_stop:]
                if len(flank) >= eff_min_flank:
                    return {
                        "terminal_end": "5p",
                        "tract_orientation": "rev",
                        "tract_start": r_start,
                        "tract_stop": t_stop,
                        "tract_len": t_stop - r_start,
                        "non_telo_flank": flank,
                    }

    return None
