import os
from pathlib import Path
import re


def load_fasta(fasta_path: str | Path) -> list[tuple[str, str]]:
    """Load FASTA into a list of (header, sequence) tuples."""
    records = []
    if not os.path.exists(fasta_path):
        return records
    with open(fasta_path, encoding="utf-8") as f:
        cur_header = None
        cur_seq = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_header:
                    records.append((cur_header, "".join(cur_seq)))
                cur_header = line[1:].strip()
                cur_seq = []
            else:
                cur_seq.append(line.upper())
        if cur_header:
            records.append((cur_header, "".join(cur_seq)))
    return records


def extract_soma_len_from_header(header: str, default: int = 150) -> int:
    """Extract soma_len from FASTA header if present, else fallback to default."""
    m = re.search(r"\bsoma_len=(\d+)", header)
    if m:
        return int(m.group(1))
    return default


def _revcomp(seq: str) -> str:
    tr = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(tr)[::-1]


def mask_residual_telomeric_tracts(
    seq: str, center_soma_len: int = 150, repeat: str = "TTAGGC"
) -> str:
    """Mask residual tandem (TTAGGC)2+ tracts (>=10 bp exact or 1-mismatch dimer) except [-6..+6] around pos0."""
    if not seq or not repeat or len(seq) < 10:
        return seq
    rep = repeat.upper()
    rc_rep = _revcomp(rep)
    k = len(rep)
    exact_10mers = {(rep * 3)[i : i + 10] for i in range(k)} | {
        (rc_rep * 3)[i : i + 10] for i in range(k)
    }
    dimer_12mers = {(rep * 3)[i : i + 12] for i in range(k)} | {
        (rc_rep * 3)[i : i + 12] for i in range(k)
    }

    n = len(seq)
    s_up = seq.upper()
    to_mask = [False] * n

    for i in range(n - 9):
        if s_up[i : i + 10] in exact_10mers:
            for j in range(i, i + 10):
                to_mask[j] = True

    for i in range(n - 11):
        sub12 = s_up[i : i + 12]
        if any(
            sum(1 for a, b in zip(sub12, t, strict=False) if a != b) <= 1
            for t in dimer_12mers
        ):
            for j in range(i, i + 12):
                to_mask[j] = True

    # Preserve [-6 .. +6] around Position 0 for cleavage microhomology
    prot_start = max(0, center_soma_len - 6)
    prot_end = min(n, center_soma_len + 6)
    for j in range(prot_start, prot_end):
        to_mask[j] = False

    return "".join("N" if to_mask[j] else seq[j] for j in range(n))


def slice_window_at_radius(
    seq: str, radius: int, center_soma_len: int = 150, repeat: str = "TTAGGC"
) -> str:
    """Slice sequence centered at boundary between somatic and germline flank.

    Masks residual tandem telomeric tracts outside [-6..+6] and extracts [-radius .. +radius].
    """
    masked = mask_residual_telomeric_tracts(
        seq, center_soma_len=center_soma_len, repeat=repeat
    )
    soma_start = max(0, center_soma_len - radius)
    soma_slice = masked[soma_start:center_soma_len]

    germ_end = min(len(masked), center_soma_len + radius)
    germ_slice = masked[center_soma_len:germ_end]

    return soma_slice + germ_slice


def write_sliced_fasta(
    records: list[tuple[str, str]],
    radius: int,
    output_path: str | Path,
    center_soma_len: int = 150,
) -> int:
    """Write sliced records to output_path. Returns count of written records."""
    count = 0
    with open(output_path, "w", encoding="utf-8") as f:
        for header, seq in records:
            s_len = extract_soma_len_from_header(header, default=center_soma_len)
            sliced = slice_window_at_radius(seq, radius, center_soma_len=s_len)
            pos0 = min(s_len, radius)
            if len(sliced) >= 20:  # Minimum valid window size
                # Clean header for MEME compatibility and embed pos0 coordinate
                clean_header = header.split()[0].replace(":", "_")
                f.write(f">{clean_header}_r{radius}_p0_{pos0}\n{sliced}\n")
                count += 1
    return count

