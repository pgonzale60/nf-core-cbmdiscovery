"""Multi-radius window slicing for 300 bp break windows centered at Position 0."""

import os
from pathlib import Path


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


def slice_window_at_radius(
    seq: str, radius: int, center_soma_len: int = 150
) -> str:
    """Slice sequence centered at boundary between somatic and germline flank.

    Extracts [-radius .. 0] somatic bases and [+1 .. +radius] germline bases.
    """
    # pos 0 is at index center_soma_len - 1
    soma_start = max(0, center_soma_len - radius)
    soma_slice = seq[soma_start:center_soma_len]

    germ_end = min(len(seq), center_soma_len + radius)
    germ_slice = seq[center_soma_len:germ_end]

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
            sliced = slice_window_at_radius(
                seq, radius, center_soma_len=center_soma_len
            )
            if len(sliced) >= 20:  # Minimum valid window size
                # Clean header for MEME compatibility
                clean_header = header.split()[0].replace(":", "_")
                f.write(f">{clean_header}_r{radius}\n{sliced}\n")
                count += 1
    return count
