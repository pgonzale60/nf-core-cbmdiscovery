"""Orientation normalizer and Position 0 trimmer for telomeric junction reads."""

from typing import Any
from .detector import revcomp


def get_cyclic_prefixes(repeat: str) -> list[str]:
    """Return cyclic prefixes of repeat ordered from longest to shortest.

    For 'TTAGGC', returns ['TTAGG', 'TTAG', 'TTA', 'TT', 'T'].
    """
    n = len(repeat)
    return [repeat[:k] for k in range(n - 1, 0, -1)]


def normalize_and_trim_flank(
    read_seq: str,
    tract_info: dict[str, Any],
    repeat: str = "TTAGGC",
    max_flank_len: int = 150,
) -> dict[str, Any] | None:
    """Standardize orientation to 5' -> 3' (flank on left, telomere on right),

    strip telomere array, trim cyclic prefixes of repeat to pos 0, and
    extract the terminal somatic anchor window [-max_flank_len .. 0].

    Returns dict with keys:
        - raw_flank: str
        - trimmed_flank: str
        - anchor_seq: str ([-max_flank_len .. 0])
        - microhomology_trimmed_len: int
        - trimmed_prefix: str
    or None if invalid.
    """
    repeat = repeat.upper()
    term_end = tract_info["terminal_end"]
    t_orient = tract_info["tract_orientation"]

    # 1. Extract raw flank in 5' -> 3' orientation (flank on left, telomere on right)
    if term_end == "3p" and t_orient == "fwd":
        raw_flank = read_seq[: tract_info["tract_start"]]
    elif term_end == "5p" and t_orient == "rev":
        raw_flank = revcomp(read_seq[tract_info["tract_stop"] :])
    else:
        # Inward-pointing tracts (3p rev, 5p fwd) are strictly rejected
        return None

    if not raw_flank:
        return None

    # Strip any trailing full telomere repeat units that might remain
    stripped = raw_flank
    while stripped.endswith(repeat):
        stripped = stripped[: -len(repeat)]

    # 2. Trim back any cyclic prefix of repeat (T, TT, TTA, TTAG, TTAGG)
    prefixes = get_cyclic_prefixes(repeat)
    trimmed_prefix = ""
    trimmed_flank = stripped

    for p in prefixes:
        if trimmed_flank.endswith(p):
            trimmed_prefix = p
            trimmed_flank = trimmed_flank[: -len(p)]
            break

    if not trimmed_flank:
        return None

    # Position 0 is defined as the 3' terminal base of trimmed_flank
    # Extract [-max_flank_len .. 0]
    anchor_seq = (
        trimmed_flank[-max_flank_len:]
        if len(trimmed_flank) >= max_flank_len
        else trimmed_flank
    )

    return {
        "raw_flank": raw_flank,
        "trimmed_flank": trimmed_flank,
        "anchor_seq": anchor_seq,
        "microhomology_trimmed_len": len(trimmed_prefix),
        "trimmed_prefix": trimmed_prefix,
    }
