#!/usr/bin/env python3
# =============================================================
# Script 08: Comparative phylogenetic metrics
# =============================================================
# Computes all 9 metrics for Table 3:
#   (A) Mantel test (rho, p-value; 9,999 permutations)
#   (B) Robinson-Foulds distance and normalised RF
#   (C) Clustering Information Distance (CID)
#   (D) Bootstrap support: median, %>=70, %>=90
#   (E) Parsimony-informative sites (PIS, count and %)
#   (F) Resolution index
#
# Usage:
#   python 08_statistics.py \
#       --tree_6gene tree_6gene_620.treefile \
#       --tree_7gene tree_7gene_620.treefile \
#       --tree_wgsnp wgsnp_620.treefile \
#       --align_6gene concat_6gene_620.fasta \
#       --align_7gene concat_7gene_620.fasta \
#       --output results/
#
# Requirements: Biopython, NumPy, SciPy, dendropy, treeswift
# =============================================================

import argparse
import csv
import warnings
from pathlib import Path

import numpy as np
from scipy.stats import pearsonr
from Bio import Phylo, SeqIO, AlignIO
from io import StringIO
from collections import defaultdict

warnings.filterwarnings("ignore")


# ═══════════════════════════════════════════════════════════
# A. Mantel Test
# ═══════════════════════════════════════════════════════════
def patristic_distances(tree):
    """Compute pairwise patristic distance matrix from tree."""
    tips = [c.name for c in tree.get_terminals()]
    n    = len(tips)
    idx  = {t: i for i, t in enumerate(tips)}
    dist = np.zeros((n, n))

    for i, t1 in enumerate(tips):
        for j, t2 in enumerate(tips):
            if i < j:
                d = tree.distance(t1, t2)
                dist[i, j] = dist[j, i] = d
    return tips, dist


def mantel_test(m1, m2, n_perm=9999, seed=42):
    """
    Mantel test: correlation between upper triangles of
    two symmetric distance matrices.
    Returns: rho (Pearson r), p-value (one-tailed, > obs)
    """
    rng = np.random.default_rng(seed)
    idx = np.triu_indices(len(m1), k=1)
    v1, v2 = m1[idx], m2[idx]
    obs_r, _ = pearsonr(v1, v2)

    count = 0
    for _ in range(n_perm):
        perm = rng.permutation(len(m1))
        pm   = m1[np.ix_(perm, perm)]
        r, _ = pearsonr(pm[idx], v2)
        if r >= obs_r:
            count += 1

    p_val = (count + 1) / (n_perm + 1)
    return obs_r, p_val


# ═══════════════════════════════════════════════════════════
# B. Robinson-Foulds distance
# ═══════════════════════════════════════════════════════════
def get_bipartitions(tree):
    """Get set of bipartitions (frozensets of leaf names)."""
    tips = frozenset(c.name for c in tree.get_terminals())
    biparts = set()
    for clade in tree.find_clades(order="level"):
        if clade == tree.root:
            continue
        leaves = frozenset(c.name for c in clade.get_terminals())
        if len(leaves) > 1 and leaves != tips:
            biparts.add(leaves)
    return biparts


def rf_distance(tree1, tree2):
    """Robinson-Foulds distance and normalised RF."""
    b1 = get_bipartitions(tree1)
    b2 = get_bipartitions(tree2)
    rf = len(b1.symmetric_difference(b2))
    norm_rf = rf / (len(b1) + len(b2)) if (len(b1) + len(b2)) > 0 else 0
    return rf, norm_rf


# ═══════════════════════════════════════════════════════════
# C. Clustering Information Distance (CID)
# Simplified implementation; use TreeDist R package for
# the full implementation used in the paper.
# ═══════════════════════════════════════════════════════════
def clustering_info_distance(tree1, tree2):
    """
    Approximate CID based on Shannon entropy of bipartitions.
    For full implementation, see Smith (2020) Bioinformatics
    and the R package TreeDist.
    """
    import math

    def entropy_biparts(tree):
        tips = [c.name for c in tree.get_terminals()]
        n = len(tips)
        if n == 0:
            return 0, {}
        H = 0
        bp_probs = {}
        for clade in tree.find_clades():
            leaves = [c.name for c in clade.get_terminals()]
            k = len(leaves)
            if 1 < k < n:
                p = k / n
                h = -p * math.log2(p) - (1-p) * math.log2(1-p)
                H += h
                bp_probs[frozenset(leaves)] = h
        return H, bp_probs

    H1, bp1 = entropy_biparts(tree1)
    H2, bp2 = entropy_biparts(tree2)

    # Mutual information approximation
    shared = set(bp1.keys()) & set(bp2.keys())
    MI = sum(min(bp1[k], bp2[k]) for k in shared)

    # CID
    joint_H = H1 + H2 - MI
    cid = 1 - (MI / joint_H) if joint_H > 0 else 0
    return cid


# ═══════════════════════════════════════════════════════════
# D. Bootstrap support
# ═══════════════════════════════════════════════════════════
def bootstrap_stats(tree):
    """Median BS, %>=70, %>=90."""
    bs_vals = []
    for clade in tree.find_clades():
        if clade.confidence is not None:
            bs_vals.append(clade.confidence)
    if not bs_vals:
        return None, None, None
    bs = np.array(bs_vals)
    return (np.median(bs),
            np.mean(bs >= 70) * 100,
            np.mean(bs >= 90) * 100)


# ═══════════════════════════════════════════════════════════
# E. Parsimony-informative sites (PIS)
# ═══════════════════════════════════════════════════════════
def count_pis(fasta_file):
    """Count parsimony-informative sites in alignment."""
    aln = AlignIO.read(fasta_file, "fasta")
    n_seq  = len(aln)
    n_cols = aln.get_alignment_length()
    pis = 0

    for col in range(n_cols):
        column = [str(aln[i].seq[col]).upper()
                  for i in range(n_seq)]
        # Count character states (excluding gaps and N)
        counts = defaultdict(int)
        for c in column:
            if c not in ("-", "N", "?"):
                counts[c] += 1
        # PIS: >= 2 states, each in >= 2 sequences
        informative = [v for v in counts.values() if v >= 2]
        if len(informative) >= 2:
            pis += 1

    return pis, n_cols, pis / n_cols * 100 if n_cols > 0 else 0


# ═══════════════════════════════════════════════════════════
# F. Resolution index
# ═══════════════════════════════════════════════════════════
def resolution_index(tree):
    """RI = n_internal / (n_leaves - 1)."""
    n_leaves   = len(list(tree.get_terminals()))
    n_internal = len([c for c in tree.find_clades()
                      if not c.is_terminal()])
    max_internal = n_leaves - 1
    return n_internal, max_internal, \
           n_internal / max_internal if max_internal > 0 else 0


# ═══════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════
def main():
    parser = argparse.ArgumentParser(
        description="Compute comparative phylogenetic metrics")
    parser.add_argument("--tree_6gene", required=True)
    parser.add_argument("--tree_7gene", required=True)
    parser.add_argument("--tree_wgsnp", required=True)
    parser.add_argument("--align_6gene", required=True)
    parser.add_argument("--align_7gene", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--n_perm", type=int, default=9999)
    args = parser.parse_args()

    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=== Comparative Phylogenetic Metrics ===")
    print(f"Mantel permutations: {args.n_perm:,}")

    # Load trees
    print("\nLoading trees...")
    t6 = Phylo.read(args.tree_6gene, "newick")
    t7 = Phylo.read(args.tree_7gene, "newick")
    tw = Phylo.read(args.tree_wgsnp, "newick")

    # Get common taxa
    tips6 = set(c.name for c in t6.get_terminals())
    tips7 = set(c.name for c in t7.get_terminals())
    tipsw = set(c.name for c in tw.get_terminals())
    common = tips6 & tips7 & tipsw
    print(f"Common taxa: {len(common)}")

    # Prune to common taxa
    def prune(tree, keep):
        for tip in list(tree.get_terminals()):
            if tip.name not in keep:
                tree.prune(tip.name)
        return tree

    t6 = prune(t6, common)
    t7 = prune(t7, common)
    tw = prune(tw, common)

    # Compute patristic distance matrices
    print("Computing patristic distance matrices...")
    tips6, d6 = patristic_distances(t6)
    tips7, d7 = patristic_distances(t7)
    tipsw, dw = patristic_distances(tw)

    # Align matrices (same tip order)
    tip_order = sorted(common)
    idx6 = [tips6.index(t) for t in tip_order]
    idx7 = [tips7.index(t) for t in tip_order]
    idxw = [tipsw.index(t) for t in tip_order]
    d6 = d6[np.ix_(idx6, idx6)]
    d7 = d7[np.ix_(idx7, idx7)]
    dw = dw[np.ix_(idxw, idxw)]

    results = {}

    # ── Mantel test ──────────────────────────────────────────
    print("\nMantel tests (may take a few minutes)...")
    r6w, p6w = mantel_test(d6, dw, args.n_perm)
    r7w, p7w = mantel_test(d7, dw, args.n_perm)
    print(f"  6-gene vs wgSNP: rho={r6w:.4f}, p={p6w:.4f}")
    print(f"  7-gene vs wgSNP: rho={r7w:.4f}, p={p7w:.4f}")
    results["mantel_rho_6_wgsnp"] = r6w
    results["mantel_p_6_wgsnp"]   = p6w
    results["mantel_rho_7_wgsnp"] = r7w
    results["mantel_p_7_wgsnp"]   = p7w

    # ── RF distance ──────────────────────────────────────────
    print("\nRF distances...")
    rf6w,  nrf6w  = rf_distance(t6, tw)
    rf7w,  nrf7w  = rf_distance(t7, tw)
    rf67,  nrf67  = rf_distance(t6, t7)
    print(f"  6-gene vs wgSNP: RF={rf6w}, nRF={nrf6w:.4f}")
    print(f"  7-gene vs wgSNP: RF={rf7w}, nRF={nrf7w:.4f}")
    print(f"  6-gene vs 7-gene: RF={rf67}, nRF={nrf67:.4f}")
    results["RF_6_wgsnp"]  = rf6w;  results["nRF_6_wgsnp"]  = nrf6w
    results["RF_7_wgsnp"]  = rf7w;  results["nRF_7_wgsnp"]  = nrf7w
    results["RF_6_7"]      = rf67;  results["nRF_6_7"]      = nrf67

    # ── CID ──────────────────────────────────────────────────
    print("\nClustering information distance...")
    cid6w = clustering_info_distance(t6, tw)
    cid7w = clustering_info_distance(t7, tw)
    print(f"  6-gene vs wgSNP: CID={cid6w:.4f}")
    print(f"  7-gene vs wgSNP: CID={cid7w:.4f}")
    results["CID_6_wgsnp"] = cid6w
    results["CID_7_wgsnp"] = cid7w

    # ── Bootstrap ────────────────────────────────────────────
    print("\nBootstrap support statistics...")
    for name, tree in [("6gene", t6), ("7gene", t7), ("wgsnp", tw)]:
        med, pct70, pct90 = bootstrap_stats(tree)
        print(f"  {name}: median={med:.1f}%, "
              f">=70%: {pct70:.1f}%, >=90%: {pct90:.1f}%")
        results[f"bs_median_{name}"]  = med
        results[f"bs_pct70_{name}"]   = pct70
        results[f"bs_pct90_{name}"]   = pct90

    # ── PIS ──────────────────────────────────────────────────
    print("\nParsimony-informative sites...")
    pis6, len6, pct6 = count_pis(args.align_6gene)
    pis7, len7, pct7 = count_pis(args.align_7gene)
    print(f"  6-gene: {pis6} PIS / {len6} bp ({pct6:.2f}%)")
    print(f"  7-gene: {pis7} PIS / {len7} bp ({pct7:.2f}%)")
    results["PIS_6"] = pis6; results["PIS_pct_6"] = pct6
    results["PIS_7"] = pis7; results["PIS_pct_7"] = pct7

    # ── Resolution index ─────────────────────────────────────
    print("\nResolution index...")
    for name, tree in [("6gene", t6), ("7gene", t7), ("wgsnp", tw)]:
        ni, max_ni, ri = resolution_index(tree)
        print(f"  {name}: {ni}/{max_ni} ({ri*100:.1f}%)")
        results[f"RI_{name}"] = ri

    # ── Save results ─────────────────────────────────────────
    out_csv = out_dir / "comparison_metrics.csv"
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Metric", "Value"])
        for k, v in results.items():
            w.writerow([k, f"{v:.4f}" if isinstance(v, float) else v])

    print(f"\n✓ Results saved: {out_csv}")


if __name__ == "__main__":
    main()
