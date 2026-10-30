"""K-mer bait index builder and FASTQ streaming baiter."""

from collections import defaultdict
from collections.abc import Iterator
import gzip
import itertools
import os
from pathlib import Path
from typing import Any


def revcomp(seq: str) -> str:
    """Return reverse complement of a DNA sequence."""
    trans = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(trans)[::-1]


def load_candidate_anchors(fasta_path: str | Path) -> dict[str, str]:
    """Load candidate anchors from FASTA into a dictionary: anchor_id -> sequence."""
    anchors = {}
    if not os.path.exists(fasta_path):
        return anchors
    with open(fasta_path, encoding="utf-8") as f:
        cur_id = None
        cur_seq = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id:
                    anchors[cur_id] = "".join(cur_seq)
                cur_id = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line.upper())
        if cur_id:
            anchors[cur_id] = "".join(cur_seq)
    return anchors


def build_kmer_bait_index(
    anchors: dict[str, str], kmer_size: int = 20, terminal_window: int = 60
) -> dict[str, list[str]]:
    """Build bait dictionary: kmer -> list of anchor_ids from the terminal window of each anchor."""
    bait_index = defaultdict(list)
    for anchor_id, seq in anchors.items():
        subseq = seq[-terminal_window:] if len(seq) >= terminal_window else seq
        if len(subseq) < kmer_size:
            continue
        for i in range(len(subseq) - kmer_size + 1):
            kmer = subseq[i : i + kmer_size]
            if "N" in kmer:
                continue
            bait_index[kmer].append(anchor_id)
            rc_kmer = revcomp(kmer)
            if rc_kmer != kmer:
                bait_index[rc_kmer].append(anchor_id)
    # Deduplicate anchor lists
    return {k: sorted(set(v)) for k, v in bait_index.items()}


def stream_fastq_records(fastq_path: str | Path) -> Iterator[tuple[str, str]]:
    """Stream FASTQ yielding (header, seq)."""
    path_str = str(fastq_path)
    opener = (
        gzip.open(path_str, "rt", encoding="utf-8", errors="replace")
        if path_str.endswith(".gz")
        else open(path_str, "r", encoding="utf-8", errors="replace")
    )
    with opener as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            s = fh.readline()
            _p = fh.readline()
            _q = fh.readline()
            yield h.strip()[1:].split()[0], s.strip().upper()


def stream_paired_fastq_records(
    fastq_1: str | Path,
    fastq_2: str | Path | None = None,
) -> Iterator[tuple[tuple[str, str] | None, tuple[str, str] | None]]:
    """Stream paired FASTQ records in lockstep."""
    r1_stream = stream_fastq_records(fastq_1)
    if fastq_2 and os.path.exists(fastq_2):
        r2_stream = stream_fastq_records(fastq_2)
        yield from itertools.zip_longest(r1_stream, r2_stream, fillvalue=None)
    else:
        for r1 in r1_stream:
            yield r1, None


def bait_reads(
    fastq_1: str | Path,
    fastq_2: str | Path | None = None,
    bait_index: dict[str, list[str]] | None = None,
    kmer_size: int = 20,
) -> Iterator[dict[str, Any]]:
    """Stream raw reads and yield baited records: {'anchor_id': ..., 'read_id': ..., 'seq': ...}."""
    if not bait_index:
        return

    def get_hits(read_seq: str) -> set[str]:
        hits = set()
        if len(read_seq) < kmer_size:
            return hits
        for i in range(len(read_seq) - kmer_size + 1):
            km = read_seq[i : i + kmer_size]
            if km in bait_index:
                hits.update(bait_index[km])
        return hits

    files = [fastq_1]
    if fastq_2 and os.path.exists(fastq_2):
        files.append(fastq_2)

    for fq in files:
        if not fq:
            continue
        for qname, seq in stream_fastq_records(fq):
            hits = get_hits(seq)
            for a_id in hits:
                yield {"anchor_id": a_id, "read_id": qname, "seq": seq}
