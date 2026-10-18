"""Partitioner of baited reads at Position 0 into telomeric and germline branches."""

from typing import Any
from .baiter import revcomp


def get_cyclic_permutations(repeat: str) -> list[str]:
    """Return all unique cyclic permutations of repeat."""
    n = len(repeat)
    doubled = repeat + repeat
    return sorted({doubled[i : i + n] for i in range(n)})


def find_probe_alignment(
    probe: str, read_seq: str, max_mismatches: int = 2
) -> tuple[int, int, bool] | None:
    """Find probe or revcomp(probe) in read_seq allowing small mismatches.

    Returns (match_start, match_stop, is_revcomp) or None.
    """
    p_len = len(probe)
    if len(read_seq) < p_len:
        return None

    # 1. Exact match check (fast path)
    idx = read_seq.find(probe)
    if idx >= 0:
        return idx, idx + p_len, False

    rc_probe = revcomp(probe)
    idx_rc = read_seq.find(rc_probe)
    if idx_rc >= 0:
        return idx_rc, idx_rc + p_len, True

    # 2. Inexact match with 15-mer seeds
    seed_len = 15
    for is_rc, p in [(False, probe), (True, rc_probe)]:
        seed = p[:seed_len]
        s_idx = read_seq.find(seed)
        while s_idx >= 0:
            if s_idx + p_len <= len(read_seq):
                sub = read_seq[s_idx : s_idx + p_len]
                mismatches = sum(
                    1 for a, b in zip(sub, p, strict=False) if a != b
                )
                if mismatches <= max_mismatches:
                    return s_idx, s_idx + p_len, is_rc
            s_idx = read_seq.find(seed, s_idx + 1)

    return None


def is_telomeric_extension(
    flank: str, repeat: str = "TTAGGC", scan_window: int = 60
) -> bool:
    """Determine whether a sequence extension consists predominantly of telomeric repeats.

    Handles telomerase addition wobble and microhomology overhangs (5-15 bp) before
    the canonical repeat array by scanning for tandem repeats within the first 60 bp.
    """
    if not flank:
        return False

    perms = get_cyclic_permutations(repeat)
    rep_len = len(repeat)

    # 1. Tandem repeat scan: check for >= 2 consecutive repeats anywhere in first scan_window bp
    scan_sub = flank[:scan_window]
    for p in perms:
        if (p * 2) in scan_sub:
            return True

    # 2. Windowed telomeric density check across prefix (up to 48 bp)
    check_len = min(len(flank), 48)
    if check_len < rep_len:
        return False

    sub = flank[:check_len]
    telo_cov = 0
    i = 0
    while i <= len(sub) - rep_len:
        if sub[i : i + rep_len] in perms:
            telo_cov += rep_len
            i += rep_len
        else:
            i += 1

    return (telo_cov / len(sub)) >= 0.50


def partition_baited_read(
    anchor_seq: str,
    read_seq: str,
    repeat: str = "TTAGGC",
    probe_len: int = 40,
    min_germ_flank_bp: int = 15,
) -> dict[str, Any] | None:
    """Align anchor terminal probe against read_seq and partition into telo or germline branch.

    Returns dict with keys:
        - branch: 'telo' or 'germline'
        - germline_ext: str (if branch == 'germline')
        - telo_ext: str (if branch == 'telo')
    or None if probe does not align.
    """
    probe = (
        anchor_seq[-probe_len:]
        if len(anchor_seq) >= probe_len
        else anchor_seq
    )
    align_res = find_probe_alignment(probe, read_seq, max_mismatches=2)
    if not align_res:
        if len(anchor_seq) >= 20:
            short_probe = anchor_seq[-20:]
            rc_short = revcomp(short_probe)
            if read_seq.startswith(short_probe):
                match_stop = len(short_probe)
                oriented_read = read_seq
            elif read_seq.endswith(rc_short):
                oriented_read = revcomp(read_seq)
                match_stop = len(short_probe)
            else:
                return None
        else:
            return None
    else:
        match_start, match_stop, is_rc = align_res
        # Forward-orient the read so probe is on forward strand
        oriented_read = revcomp(read_seq) if is_rc else read_seq
        if is_rc:
            # Re-locate probe in oriented read
            idx = oriented_read.find(probe)
            if idx >= 0:
                match_stop = idx + len(probe)
            else:
                match_stop = len(read_seq) - match_start

    flank_after_pos0 = oriented_read[match_stop:]

    if is_telomeric_extension(flank_after_pos0, repeat=repeat):
        return {
            "branch": "telo",
            "telo_ext": flank_after_pos0,
        }

    if len(flank_after_pos0) >= min_germ_flank_bp:
        return {
            "branch": "germline",
            "germline_ext": flank_after_pos0,
        }

    return None
