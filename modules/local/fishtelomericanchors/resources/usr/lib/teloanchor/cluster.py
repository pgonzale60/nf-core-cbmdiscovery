"""Clustering of somatic anchor sequences at >= 97% identity."""

from collections import Counter
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
    min_identity: float = 0.97,
    min_cluster_size: int = 1,
) -> list[dict[str, Any]]:
    """Greedy centroid clustering of anchor sequences anchored at 3' end (pos 0).

    Returns list of clusters, each containing:
        - anchor_id: str
        - consensus_seq: str
        - read_count: int
        - read_ids: list[str]
        - observations: list[dict]
    """
    if not observations:
        return []

    # Sort observations by sequence length descending so longest anchor serves as seed
    sorted_obs = sorted(
        observations, key=lambda x: len(x["anchor_seq"]), reverse=True
    )

    clusters: list[dict[str, Any]] = []

    for obs in sorted_obs:
        seq = obs["anchor_seq"]
        matched_cluster = None
        for cl in clusters:
            ident = sequence_identity_anchored_3p(cl["rep_seq"], seq)
            if ident >= min_identity:
                matched_cluster = cl
                break

        if matched_cluster is not None:
            matched_cluster["observations"].append(obs)
        else:
            clusters.append({
                "rep_seq": seq,
                "observations": [obs],
            })

    # Filter by min_cluster_size and build consensus
    final_clusters = []
    for cl in clusters:
        obs_list = cl["observations"]
        if len(obs_list) < min_cluster_size:
            continue

        # Build 3'-anchored consensus sequence
        max_len = max(len(o["anchor_seq"]) for o in obs_list)
        consensus_chars = []
        for offset_from_3p in range(max_len, 0, -1):
            bases = []
            for o in obs_list:
                s = o["anchor_seq"]
                if len(s) >= offset_from_3p:
                    bases.append(s[-offset_from_3p])
            if bases:
                top_base = Counter(bases).most_common(1)[0][0]
                consensus_chars.append(top_base)

        consensus_seq = "".join(consensus_chars)
        final_clusters.append({
            "consensus_seq": consensus_seq,
            "read_count": len(obs_list),
            "read_ids": [o["read_id"] for o in obs_list],
            "observations": obs_list,
        })

    # Sort clusters by read_count descending
    final_clusters.sort(key=lambda c: c["read_count"], reverse=True)

    for idx, cl in enumerate(final_clusters, start=1):
        cl["anchor_id"] = f"anchor_{idx}"

    return final_clusters
