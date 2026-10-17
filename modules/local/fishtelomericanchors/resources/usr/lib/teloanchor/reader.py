"""FastQ reader supporting plain and gzipped single-end and paired-end files."""

from collections.abc import Iterator
import gzip
import os
from pathlib import Path


def open_file_maybe_gzip(path: str | Path, mode: str = "rt"):
    """Open a file with gzip if ends with .gz, else plain open."""
    path_str = str(path)
    if path_str.endswith(".gz"):
        return gzip.open(path_str, mode, encoding="utf-8", errors="replace")
    return open(path_str, mode, encoding="utf-8", errors="replace")


def stream_fastq(fastq_path: str | Path) -> Iterator[tuple[str, str, str]]:
    """Yield (header, sequence, quality) from a FASTQ file."""
    if not os.path.exists(fastq_path):
        raise FileNotFoundError(f"FASTQ file not found: {fastq_path}")

    with open_file_maybe_gzip(fastq_path) as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline()
            _plus = fh.readline()
            qual = fh.readline()
            if not qual:
                break
            qname = header.strip()[1:].split()[0]
            yield qname, seq.strip().upper(), qual.strip()


def stream_paired_fastq(
    r1_path: str | Path, r2_path: str | Path
) -> Iterator[tuple[tuple[str, str, str], tuple[str, str, str]]]:
    """Yield ((qname1, seq1, qual1), (qname2, seq2, qual2)) from paired FASTQ files."""
    iter1 = stream_fastq(r1_path)
    iter2 = stream_fastq(r2_path)
    for rec1, rec2 in zip(iter1, iter2, strict=False):
        yield rec1, rec2
