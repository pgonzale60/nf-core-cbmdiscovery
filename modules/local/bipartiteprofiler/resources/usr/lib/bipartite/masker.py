"""Core motif instance locator (via FIMO) and sequence masker."""

from collections import defaultdict
import csv
import os
from pathlib import Path
import re
import shutil
import subprocess
from typing import Any

from cbm_defaults import SSOT


def load_soma_len_map(
    fasta_path: str | Path, default: int | None = None
) -> dict[str, int]:
    """Map clean sequence name to soma_len from FASTA header (fallback to seq midpoint or SSOT)."""
    soma_map = {}
    if not os.path.exists(fasta_path):
        return soma_map
    with open(fasta_path, encoding="utf-8") as f:
        cur_id = None
        cur_hdr = ""
        cur_len = 0
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None:
                    m = re.search(r"\bsoma_len=(\d+)", cur_hdr)
                    if not m:
                        m = re.search(r"_p0_(\d+)", cur_hdr)
                    soma_map[cur_id] = (
                        int(m.group(1))
                        if m
                        else (
                            default
                            if default is not None
                            else (cur_len // 2 or SSOT["max_soma_flank_len"])
                        )
                    )
                cur_id = line[1:].split()[0].replace(":", "_")
                cur_hdr = line
                cur_len = 0
            else:
                cur_len += len(line)
        if cur_id is not None:
            m = re.search(r"\bsoma_len=(\d+)", cur_hdr)
            if not m:
                m = re.search(r"_p0_(\d+)", cur_hdr)
            soma_map[cur_id] = (
                int(m.group(1))
                if m
                else (
                    default
                    if default is not None
                    else (cur_len // 2 or SSOT["max_soma_flank_len"])
                )
            )
    return soma_map


def slice_window_at_radius(
    seq: str,
    radius: int,
    pos0: int,
) -> tuple[str, int]:
    """Slice sequence to [-radius .. +radius] around pos0.

    Returns:
        (sliced_seq, new_pos0)
    """
    s = seq.upper()
    n = len(s)
    if n <= 2 * radius + 10:
        return s, pos0

    soma_start = max(0, pos0 - radius)
    germ_end = min(n, pos0 + radius)
    sliced = s[soma_start:germ_end]
    new_pos0 = pos0 - soma_start
    return sliced, new_pos0


def write_sliced_fasta(
    raw_seqs: dict[str, str],
    soma_lens: dict[str, int],
    radius: int,
    output_fasta: str | Path,
) -> tuple[dict[str, str], dict[str, int]]:
    """Slice sequences to radius and write FASTA with embedded soma_len."""
    out_path = Path(output_fasta)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    sliced_seqs = {}
    new_soma_lens = {}

    with open(out_path, "w", encoding="utf-8") as f:
        for name, seq in raw_seqs.items():
            pos0 = soma_lens.get(name, len(seq) // 2)
            sliced, new_pos0 = slice_window_at_radius(seq, radius, pos0)
            sliced_seqs[name] = sliced
            new_soma_lens[name] = new_pos0
            f.write(f">{name} soma_len={new_pos0}\n{sliced}\n")

    return sliced_seqs, new_soma_lens


def load_fasta_dict(fasta_path: str | Path) -> dict[str, str]:
    """Load FASTA into dict: header -> sequence."""
    seqs = {}
    if not os.path.exists(fasta_path):
        return seqs
    with open(fasta_path, encoding="utf-8") as f:
        cur_h = None
        cur_s = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_h:
                    seqs[cur_h] = "".join(cur_s)
                # Keep ID clean (e.g. up to space)
                cur_h = line[1:].split()[0].replace(":", "_")
                cur_s = []
            else:
                cur_s.append(line.upper())
        if cur_h:
            seqs[cur_h] = "".join(cur_s)
    return seqs


def has_motifs(meme_file: str | Path) -> bool:
    """Check if MEME file contains at least one MOTIF definition."""
    if not os.path.exists(meme_file):
        return False
    with open(meme_file, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("MOTIF "):
                return True
    return False


def run_fimo(
    meme_file: str | Path,
    fasta_file: str | Path,
    out_dir: str | Path,
    threshold: float = SSOT["fimo_core_threshold"],
) -> list[dict[str, Any]]:
    """Execute FIMO to find motif occurrences passing threshold."""
    if not has_motifs(meme_file):
        return []

    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    fimo_tsv = out_path / "fimo.tsv"

    if not fimo_tsv.exists():
        fimo_bin = shutil.which("fimo")
        if not fimo_bin:
            raise RuntimeError("Required binary 'fimo' not found in PATH")
        cmd = [
            fimo_bin,
            "--thresh",
            str(threshold),
            "--oc",
            str(out_path),
            str(meme_file),
            str(fasta_file),
        ]
        subprocess.run(
            cmd,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True,
        )

    hits = []
    if fimo_tsv.exists():
        with open(fimo_tsv, encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for r in reader:
                s = r.get("sequence_name")
                if not s or s.startswith("#"):
                    continue
                try:
                    hits.append({
                        "seq_name": s.replace(":", "_"),
                        "motif_id": r.get("motif_id", ""),
                        "start": int(r["start"]),  # 1-based inclusive
                        "stop": int(r["stop"]),  # 1-based inclusive
                        "strand": r.get("strand", "+"),
                        "score": float(r.get("score", 0.0)),
                        "p_value": float(r.get("p-value", 1.0)),
                        "matched_seq": r.get("matched_sequence", ""),
                    })
                except (ValueError, KeyError):
                    continue

    return hits


def mask_sequences_with_hits(
    seq_dict: dict[str, str], hits: list[dict[str, Any]]
) -> dict[str, str]:
    """Mask hits with 'N' in each sequence."""
    hits_by_seq = defaultdict(list)
    for h in hits:
        hits_by_seq[h["seq_name"]].append(h)

    masked_seqs = {}
    for name, seq in seq_dict.items():
        char_list = list(seq)
        for h in hits_by_seq.get(name, []):
            start_0 = max(0, h["start"] - 1)
            stop_1 = min(len(char_list), h["stop"])
            for i in range(start_0, stop_1):
                char_list[i] = "N"
        masked_seqs[name] = "".join(char_list)

    return masked_seqs


def _revcomp(seq: str) -> str:
    tr = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(tr)[::-1]


def mask_residual_telomeric_tracts(
    seq: str,
    center_soma_len: int = SSOT["max_soma_flank_len"],
    repeat: str = SSOT["telomere_repeat"],
) -> str:
    """Mask residual tandem (TTAGGC)2+ tracts (>=10 bp exact or 1-mismatch dimer) except [-6..+6] around pos0."""
    if not seq or not repeat or len(seq) < 10:
        return seq
    rep = repeat.upper()
    rc_rep = _revcomp(rep)
    k = len(rep)
    min_exact = min(10, 2 * k)
    min_dimer = 2 * k
    exact_kmers = {(rep * 3)[i : i + min_exact] for i in range(k)} | {
        (rc_rep * 3)[i : i + min_exact] for i in range(k)
    }
    dimer_kmers = {(rep * 3)[i : i + min_dimer] for i in range(k)} | {
        (rc_rep * 3)[i : i + min_dimer] for i in range(k)
    }

    n = len(seq)
    s_up = seq.upper()
    to_mask = [False] * n

    for i in range(n - min_exact + 1):
        if s_up[i : i + min_exact] in exact_kmers:
            for j in range(i, i + min_exact):
                to_mask[j] = True

    for i in range(n - min_dimer + 1):
        sub_d = s_up[i : i + min_dimer]
        if any(
            sum(1 for a, b in zip(sub_d, t, strict=False) if a != b) <= 1
            for t in dimer_kmers
        ):
            for j in range(i, i + min_dimer):
                to_mask[j] = True

    prot_r = SSOT["cleavage_mh_protection_radius"]
    prot_start = max(0, center_soma_len - prot_r)
    prot_end = min(n, center_soma_len + prot_r)
    for j in range(prot_start, prot_end):
        to_mask[j] = False

    return "".join("N" if to_mask[j] else seq[j] for j in range(n))


def mask_core_motifs(
    core_meme: str | Path,
    fasta_file: str | Path,
    output_masked_fasta: str | Path,
    fimo_out_dir: str | Path,
    threshold: float = SSOT["fimo_core_threshold"],
    repeat: str = SSOT["telomere_repeat"],
) -> tuple[dict[str, str], list[dict[str, Any]]]:
    """Run FIMO with core motif on FASTA, mask instances and residual telomeric tracts with N, and save masked FASTA."""
    seq_dict = load_fasta_dict(fasta_file)
    soma_lens = load_soma_len_map(fasta_file)
    hits = run_fimo(
        core_meme, fasta_file, out_dir=fimo_out_dir, threshold=threshold
    )
    masked_dict = mask_sequences_with_hits(seq_dict, hits)
    for name, mseq in list(masked_dict.items()):
        pos0 = soma_lens.get(name, len(mseq) // 2)
        masked_dict[name] = mask_residual_telomeric_tracts(
            mseq, center_soma_len=pos0, repeat=repeat
        )

    with open(output_masked_fasta, "w", encoding="utf-8") as f:
        for name, mseq in masked_dict.items():
            f.write(f">{name}\n{mseq}\n")

    return masked_dict, hits
