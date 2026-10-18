"""Clustering of somatic anchor sequences at >= 97% identity via vsearch."""

from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from typing import Any


def sequence_identity_anchored_3p(s1: str, s2: str) -> float:
    """Calculate sequence identity between two sequences anchored at their 3' ends (pos 0)."""
    l_min = min(len(s1), len(s2))
    if l_min == 0:
        return 0.0
    u1 = s1[-l_min:]
    u2 = s2[-l_min:]
    matches = sum(1 for a, b in zip(u1, u2, strict=False) if a == b)
    return matches / l_min


def cluster_anchors(
    observations: list[dict[str, Any]],
    min_identity: float = 0.98,
    min_cluster_size: int = 1,
    threads: int = 1,
) -> list[dict[str, Any]]:
    """Cluster somatic anchor sequences using vsearch --cluster_fast.

    Returns list of clusters, each containing:
        - anchor_id: str ("anchor_1", "anchor_2", ...)
        - consensus_seq: str
        - read_count: int
        - read_ids: list[str]
        - observations: list[dict]
    """
    if not observations:
        return []

    vsearch_bin = shutil.which("vsearch")
    if not vsearch_bin:
        raise RuntimeError("Required binary 'vsearch' not found in PATH")

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        raw_fasta = tmp_path / "raw_flanks.fasta"
        cons_fasta = tmp_path / "vsearch_cons.fasta"
        uc_file = tmp_path / "vsearch_clusters.uc"

        obs_by_id = {}
        with open(raw_fasta, "w", encoding="utf-8") as f_fa:
            for idx, obs in enumerate(observations):
                rec_id = f"obs_{idx}"
                obs_by_id[rec_id] = obs
                seq = obs["anchor_seq"].upper().strip()
                f_fa.write(f">{rec_id}\n{seq}\n")

        cmd = [
            vsearch_bin,
            "--cluster_fast",
            str(raw_fasta),
            "--id",
            str(min_identity),
            "--iddef",
            "0",
            "--minseqlength",
            "1",
            "--qmask",
            "none",
            "--strand",
            "plus",
            "--consout",
            str(cons_fasta),
            "--uc",
            str(uc_file),
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

        cluster_members: dict[int, list[dict[str, Any]]] = {}
        centroid_to_cluster: dict[str, int] = {}

        if uc_file.exists():
            with open(uc_file, encoding="utf-8") as f_uc:
                for line in f_uc:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split("\t")
                    rec_type = parts[0]
                    if rec_type == "S":
                        cl_idx = int(parts[1])
                        query_id = parts[8]
                        cluster_members.setdefault(cl_idx, []).append(obs_by_id[query_id])
                        centroid_to_cluster[query_id] = cl_idx
                    elif rec_type == "H":
                        cl_idx = int(parts[1])
                        query_id = parts[8]
                        cluster_members.setdefault(cl_idx, []).append(obs_by_id[query_id])

        clusters: list[tuple[str, str]] = []
        if cons_fasta.exists():
            cur_header = None
            cur_seq: list[str] = []
            with open(cons_fasta, encoding="utf-8") as f_c:
                for line in f_c:
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        if cur_header:
                            clusters.append((cur_header, "".join(cur_seq)))
                        cur_header = line[1:]
                        cur_seq = []
                    else:
                        cur_seq.append(line)
                if cur_header:
                    clusters.append((cur_header, "".join(cur_seq)))

        final_clusters = []
        for header, c_seq in clusters:
            m_size = re.search(r"size=(\d+)", header)
            size = int(m_size.group(1)) if m_size else 1
            if size < min_cluster_size:
                continue

            m_cent = re.search(r"centroid=([^;]+)", header)
            if m_cent:
                centroid_id = m_cent.group(1)
            else:
                centroid_id = header.split(";")[0]

            cl_idx = centroid_to_cluster.get(centroid_id)
            if cl_idx is not None and cl_idx in cluster_members:
                obs_list = cluster_members[cl_idx]
            elif centroid_id in obs_by_id:
                obs_list = [obs_by_id[centroid_id]]
            else:
                obs_list = []

            cleaned_seq = c_seq.replace("-", "").replace(".", "").upper()

            final_clusters.append({
                "consensus_seq": cleaned_seq,
                "read_count": len(obs_list) if obs_list else size,
                "read_ids": [o["read_id"] for o in obs_list],
                "observations": obs_list,
            })

        final_clusters.sort(key=lambda c: c["read_count"], reverse=True)
        for idx, cl in enumerate(final_clusters, start=1):
            cl["anchor_id"] = f"anchor_{idx}"

        return final_clusters
