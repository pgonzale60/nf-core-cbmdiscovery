"""Reconstruction of germline consensus, coherence scoring, and microhomology annotation."""

from collections import Counter


def get_cyclic_prefixes(repeat: str) -> list[str]:
    """Return prefixes of repeat from length L-1 down to 1."""
    n = len(repeat)
    return [repeat[:k] for k in range(n - 1, 0, -1)]


def reconstruct_germline_consensus(
    germline_extensions: list[str], max_len: int = 150
) -> tuple[str, float]:
    """Build majority-vote germline consensus from +1 to +max_len and calculate coherence.

    Returns (consensus_seq, coherence_fraction).
    """
    if not germline_extensions:
        return "", 0.0

    consensus_chars = []
    # Build consensus position by position
    for pos in range(max_len):
        bases = [
            ext[pos] for ext in germline_extensions if len(ext) > pos
        ]
        if not bases:
            break
        top_base = Counter(bases).most_common(1)[0][0]
        consensus_chars.append(top_base)

    consensus_seq = "".join(consensus_chars)
    if not consensus_seq:
        return "", 0.0

    # Calculate Branch Coherence: fraction of aligned bases agreeing with consensus
    total_aligned = 0
    total_matching = 0
    cons_len = len(consensus_seq)

    for ext in germline_extensions:
        overlap = min(len(ext), cons_len)
        for i in range(overlap):
            total_aligned += 1
            if ext[i] == consensus_seq[i]:
                total_matching += 1

    coherence = total_matching / total_aligned if total_aligned > 0 else 0.0
    return consensus_seq, round(coherence, 4)


def annotate_microhomology(
    germline_consensus: str, repeat: str = "TTAGGC"
) -> tuple[int, str]:
    """Compare bases +1 .. +5 of germline consensus against cyclic prefixes of repeat.

    Returns (microhomology_len, microhomology_motif).
    """
    if not germline_consensus:
        return 0, ""

    prefixes = get_cyclic_prefixes(repeat.upper())
    for p in prefixes:
        if germline_consensus.startswith(p):
            return len(p), p

    return 0, ""


def build_300bp_break_window(
    somatic_anchor: str,
    germline_consensus: str,
    soma_flank_target: int = 150,
    germ_flank_target: int = 150,
) -> str:
    """Construct 300 bp break window centered at Position 0:

    [-soma_flank_target .. 0] Retained Somatic + [+1 .. +germ_flank_target] Eliminated Germline.
    """
    soma_flank = (
        somatic_anchor[-soma_flank_target:]
        if len(somatic_anchor) >= soma_flank_target
        else somatic_anchor
    )
    germ_flank = germline_consensus[:germ_flank_target]
    return soma_flank + germ_flank
