"""Reconstruction of germline consensus, coherence scoring, and microhomology annotation."""

from collections import Counter
from pathlib import Path
import shutil
import subprocess
import tempfile


def get_cyclic_prefixes(repeat: str) -> list[str]:
    """Return prefixes of repeat from length L-1 down to 1."""
    n = len(repeat)
    return [repeat[:k] for k in range(n - 1, 0, -1)]


def revcomp(seq: str) -> str:
    """Reverse-complement a DNA sequence."""
    tr = str.maketrans(
        "ACGTURYKMSWBDHVNacgturykmswbdhvn", "TGCAAYRMKWSVHDBNtgcaayrmkwsvhdbn"
    )
    return seq.translate(tr)[::-1]


def banded_alignment_identity(
    query: str, target: str, max_band: int = 3
) -> float:
    """Compute alignment identity allowing small indels (banded semi-global alignment).

    Tolerates 1-2 bp homopolymer indels common in long reads (PacBio HiFi).
    """
    m = len(query)
    n = len(target)
    if m == 0 or n == 0:
        return 0.0

    min_l = min(m, n)
    raw_matches = sum(1 for i in range(min_l) if query[i] == target[i])
    if raw_matches / min_l >= 0.95:
        return raw_matches / min_l

    dp = [0] * (n + 1)
    for i in range(1, m + 1):
        prev = dp[:]
        j_start = max(1, i - max_band)
        j_end = min(n, i + max_band)
        for j in range(j_start, j_end + 1):
            match = 1 if query[i - 1] == target[j - 1] else 0
            best = prev[j - 1] + match
            if j > j_start and dp[j - 1] > best:
                best = dp[j - 1]
            if prev[j] > best:
                best = prev[j]
            dp[j] = best

    end_j_start = max(1, m - max_band)
    end_j_end = min(n, m + max_band)
    best_matches = (
        max(dp[j] for j in range(end_j_start, end_j_end + 1))
        if end_j_end >= end_j_start
        else dp[-1]
    )
    align_len = max(m, n)
    return best_matches / align_len if align_len > 0 else 0.0


def reconstruct_germline_consensus(
    germline_extensions: list[str],
    max_len: int = 150,
    anchor_seq: str = "",
    local_reads: list[str] | None = None,
    mate_reads: list[str] | None = None,
    threads: int = 1,
) -> tuple[str, float]:
    """Reconstruct [+1 .. +max_len] germline consensus via vsearch (Case A) or fermi-lite (Case B).

    Returns (consensus_seq, branch_coherence).
    """
    if not germline_extensions:
        return "", 0.0

    vsearch_bin = shutil.which("vsearch")
    fml_bin = shutil.which("fml-asm")
    if not vsearch_bin:
        raise RuntimeError("Required binary 'vsearch' not found in PATH")
    if not fml_bin:
        raise RuntimeError("Required binary 'fml-asm' not found in PATH")

    all_mates = list(local_reads or []) + list(mate_reads or [])
    consensus_seq = ""

    # Case A: Direct Full-Length Overhangs (>= 2 extensions with length >= max_len)
    long_exts = [e for e in germline_extensions if len(e) >= max_len]
    if len(long_exts) >= 2:
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            fa_path = tmp_path / "long_exts.fa"
            cons_path = tmp_path / "cons.fa"
            with open(fa_path, "w", encoding="utf-8") as f_fa:
                for idx, ext in enumerate(long_exts):
                    f_fa.write(f">ext_{idx}\n{ext[:max_len]}\n")

            cmd = [
                vsearch_bin,
                "--cluster_fast",
                str(fa_path),
                "--id",
                "0.95",
                "--iddef",
                "0",
                "--minseqlength",
                "1",
                "--qmask",
                "none",
                "--strand",
                "plus",
                "--consout",
                str(cons_path),
                "--sizeout",
                "--threads",
                str(max(1, threads)),
            ]
            subprocess.run(
                cmd,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=True,
            )
            if cons_path.exists():
                lines = [
                    line.strip()
                    for line in cons_path.read_text(encoding="utf-8").splitlines()
                    if line.strip() and not line.startswith(">")
                ]
                if lines:
                    c = "".join(lines).replace("-", "").replace(".", "").upper()
                    if len(c) >= max_len:
                        consensus_seq = c[:max_len]

    # Case B: Assembly via fermi-lite (fml-asm)
    if len(consensus_seq) < max_len and (anchor_seq or all_mates):
        reads_to_assemble = list(all_mates)
        for ext in germline_extensions[:200]:
            r_full = (anchor_seq[-50:] if anchor_seq else "") + ext
            if r_full and r_full not in reads_to_assemble:
                reads_to_assemble.append(r_full)
        if len(reads_to_assemble) > 500:
            reads_to_assemble = reads_to_assemble[:500]

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            fq_path = tmp_path / "local_pool.fastq"
            with open(fq_path, "w", encoding="utf-8") as f_fq:
                if anchor_seq:
                    a_seed = (
                        anchor_seq[-150:]
                        if len(anchor_seq) >= 150
                        else anchor_seq
                    )
                    for s_idx in range(3):
                        f_fq.write(
                            f"@anchor_seed_{s_idx}\n{a_seed}\n+\n{'I' * len(a_seed)}\n"
                        )
                for r_idx, r_seq in enumerate(reads_to_assemble):
                    f_fq.write(
                        f"@read_{r_idx}\n{r_seq}\n+\n{'I' * len(r_seq)}\n"
                    )

            unitigs: list[str] = []
            for fml_args in [
                ["-e", "-1", "-l", "35", "-c", "1,1", "-O"],
                ["-e", "21", "-l", "35", "-c", "2", "-t", "1"],
            ]:
                try:
                    res = subprocess.run(
                        [fml_bin] + fml_args + [str(fq_path)],
                        capture_output=True,
                        text=True,
                        timeout=30,
                    )
                    raw_out = res.stdout.strip()
                    if raw_out:
                        u_lines = raw_out.splitlines()
                        for i in range(1, len(u_lines), 4):
                            unitigs.append(u_lines[i].strip().upper())
                        if unitigs:
                            break
                except subprocess.TimeoutExpired:
                    continue

            if anchor_seq:
                for p_len in [50, 40, 30, 20]:
                    if len(anchor_seq) < p_len:
                        continue
                    probe = anchor_seq[-p_len:].upper()
                    for utg in unitigs:
                        idx = utg.find(probe)
                        if idx >= 0:
                            after = utg[idx + p_len :]
                            if len(after) > len(consensus_seq):
                                consensus_seq = after[:max_len]
                        rc_utg = revcomp(utg)
                        idx_rc = rc_utg.find(probe)
                        if idx_rc >= 0:
                            after = rc_utg[idx_rc + p_len :]
                            if len(after) > len(consensus_seq):
                                consensus_seq = after[:max_len]
                    if len(consensus_seq) >= max_len:
                        break
            elif unitigs:
                # Without anchor_seq, pick unitig matching germline extension prefix
                ext0 = germline_extensions[0][: min(20, len(germline_extensions[0]))]
                for utg in unitigs:
                    if ext0 in utg:
                        idx = utg.find(ext0)
                        sub = utg[idx : idx + max_len]
                        if len(sub) > len(consensus_seq):
                            consensus_seq = sub
                    elif ext0 in revcomp(utg):
                        rc_utg = revcomp(utg)
                        idx = rc_utg.find(ext0)
                        sub = rc_utg[idx : idx + max_len]
                        if len(sub) > len(consensus_seq):
                            consensus_seq = sub

    # Fallback to column-wise consensus if needed
    if not consensus_seq:
        chars = []
        for pos in range(max_len):
            bases = [ext[pos] for ext in germline_extensions if len(ext) > pos]
            if not bases:
                break
            top_base = Counter(bases).most_common(1)[0][0]
            chars.append(top_base)
        consensus_seq = "".join(chars)

    if not consensus_seq:
        return "", 0.0

    # Calculate Branch Coherence against consensus_seq with indel tolerance
    cons_len = len(consensus_seq)
    identities = []

    for ext in germline_extensions:
        overlap = min(len(ext), cons_len)
        if overlap == 0:
            continue
        identities.append(
            banded_alignment_identity(
                ext[:overlap], consensus_seq[:overlap], max_band=3
            )
        )

    coherence = sum(identities) / len(identities) if identities else 0.0
    return consensus_seq[:max_len], round(coherence, 4)


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
