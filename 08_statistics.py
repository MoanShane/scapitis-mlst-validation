#!/usr/bin/env python3
# =============================================================
# Script 08: Three-way statistical comparison
#   (6-gene MLST vs 7-gene MLST vs wgSNP, validated against cgMLST)
# =============================================================
# Usage:
#   python 08_statistics.py \
#       --tree_6gene PATH --tree_7gene PATH --tree_wgsnp PATH \
#       --st_6gene PATH --st_7gene PATH --cgmlst PATH \
#       --output DIR
#
# Computes the 9 head-to-head metrics reported in the manuscript:
#   1. Mantel correlation (patristic distance, vs wgSNP)
#   2. Robinson-Foulds distance (vs wgSNP)
#   3. Normalised RF distance (vs wgSNP)
#   4. Clustering information distance, CID (vs wgSNP)
#   5. Internal node count
#   6. Resolution index
#   7. Median bootstrap support
#   8. Proportion of nodes with bootstrap >= 70%
#   9. Adjusted Rand Index (ARI) vs cgMLST cluster assignment
#
# IMPORTANT — two pitfalls fixed in this version:
#
#   (a) Missing-data exclusion: ST assignment CSVs produced by
#       04_extract_mlst.py mark strains with incomplete locus
#       coverage as ST == "NA" (has_missing_data == "yes"). These
#       strains MUST be excluded before computing ARI — including
#       them (even inadvertently) fragments true alleles and
#       collapses ARI from ~0.93 down to ~0.75.
#
#   (b) pandas read_csv() pitfall: by default, pandas treats the
#       literal string "NA" as a missing-value marker and silently
#       converts it to float NaN. A filter such as
#       `df[df["ST"] != "NA"]` will therefore match EVERYTHING
#       (since NaN != "NA" is always True), silently failing to
#       exclude anything. This script reads ST assignment CSVs with
#       keep_default_na=False, na_values=[] specifically to avoid
#       this, then filters on the literal string "NA" explicitly.
#
# Requirements: Biopython 1.79, NumPy 1.21, SciPy 1.7,
#               scikit-learn 1.0, pandas 1.3
# =============================================================

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import adjusted_rand_score
from Bio import Phylo

sys.setrecursionlimit(10000)


def load_st_assignments(csv_path, st_col="ST"):
    """
    Load an ST assignment CSV produced by 04_extract_mlst.py,
    correctly handling the literal string "NA" used to flag
    strains with missing locus data.

    Returns a DataFrame with rows for has_missing_data == "yes"
    already removed.
    """
    df = pd.read_csv(csv_path, keep_default_na=False, na_values=[])
    df["strain"] = df["strain"].astype(str).str.strip()

    if "has_missing_data" in df.columns:
        n_total = len(df)
        df = df[df["has_missing_data"] != "yes"].copy()
        n_excluded = n_total - len(df)
        if n_excluded > 0:
            print(f"  {csv_path.name}: excluded {n_excluded} strain(s) "
                  f"with missing locus data; {len(df)} remain")
    else:
        # Fallback for CSVs without the has_missing_data column —
        # still guard against the literal "NA" string in the ST column.
        df = df[df[st_col] != "NA"].copy()

    return df


def patristic_distances(tree_path, target_set):
    """Compute pairwise patristic distance matrix for a pruned tree."""
    tree = Phylo.read(tree_path, "newick")
    all_t = [t.name for t in tree.get_terminals()]
    for name in all_t:
        if name not in target_set:
            tree.prune(name)

    terminals = tree.get_terminals()
    names = [t.name for t in terminals]
    depths = {}
    ancestors = {}

    def build(clade, path, depth):
        bl = clade.branch_length or 0.0
        d = depth + bl
        depths[id(clade)] = d
        new_path = path + [id(clade)]
        if clade.is_terminal():
            ancestors[clade.name] = new_path
        else:
            for c in clade.clades:
                build(c, new_path, d)

    build(tree.root, [], 0.0)

    n = len(names)
    dist = np.zeros((n, n))
    for i, na in enumerate(names):
        pa = set(ancestors[na])
        da = depths[ancestors[na][-1]]
        for j, nb in enumerate(names):
            if j <= i:
                continue
            pb = ancestors[nb]
            lca = next((x for x in reversed(pb) if x in pa), None)
            db = depths[ancestors[nb][-1]]
            lca_d = depths[lca] if lca else 0.0
            d = da + db - 2 * lca_d
            dist[i, j] = dist[j, i] = d

    idx = {nm: i for i, nm in enumerate(names)}
    return names, dist, idx


def reorder_matrix(names, dist, order):
    idx = {nm: i for i, nm in enumerate(names)}
    n = len(order)
    m = np.zeros((n, n))
    for i, a in enumerate(order):
        for j, b in enumerate(order):
            m[i, j] = dist[idx[a]][idx[b]]
    return m


def mantel_test(m1, m2, perms=9999, seed=42):
    """Mantel test with permutation-based p-value."""
    tri_idx = np.triu_indices(len(m1), k=1)
    v1, v2 = m1[tri_idx], m2[tri_idx]
    if v1.std() == 0 or v2.std() == 0:
        return float("nan"), float("nan")
    obs_r, _ = stats.pearsonr(v1, v2)
    rng = np.random.default_rng(seed)
    count = 0
    for _ in range(perms):
        p = rng.permutation(len(m1))
        mp = m1[np.ix_(p, p)]
        r, _ = stats.pearsonr(mp[tri_idx], v2)
        if r >= obs_r:
            count += 1
    return obs_r, (count + 1) / (perms + 1)


def get_bipartitions(tree_path, common_set):
    tree = Phylo.read(tree_path, "newick")
    all_t = [t.name for t in tree.get_terminals()]
    for nm in all_t:
        if nm not in common_set:
            tree.prune(nm)
    leaves = frozenset(t.name for t in tree.get_terminals())
    n = len(leaves)
    bips = []
    for clade in tree.find_clades():
        if clade.is_terminal():
            continue
        sub = frozenset(t.name for t in clade.get_terminals())
        comp = leaves - sub
        if len(sub) > 1 and len(comp) > 1:
            p = len(sub) / n
            bips.append((frozenset([sub, comp]), p))
    return bips, n


def rf_distance(tree_path_1, tree_path_2, common_set):
    """Robinson-Foulds distance and normalised RF."""
    bips1, _ = get_bipartitions(tree_path_1, common_set)
    bips2, _ = get_bipartitions(tree_path_2, common_set)
    b1 = {b for b, _ in bips1}
    b2 = {b for b, _ in bips2}
    rf = len(b1.symmetric_difference(b2))
    norm = rf / (len(b1) + len(b2)) if (len(b1) + len(b2)) else 0
    return rf, norm


def cid_distance(tree_path_1, tree_path_2, common_set):
    """Clustering Information Distance (Smith 2020)."""
    bips1, n = get_bipartitions(tree_path_1, common_set)
    bips2, _ = get_bipartitions(tree_path_2, common_set)
    set1 = {b for b, _ in bips1}
    set2 = {b for b, _ in bips2}
    shared = set1 & set2

    def entropy(bips):
        h = 0
        for _, p in bips:
            if 0 < p < 1:
                h += -(p * np.log2(p) + (1 - p) * np.log2(1 - p))
        return h

    h1 = entropy(bips1)
    h2 = entropy(bips2)
    bips1_dict = {b: p for b, p in bips1}
    bips2_dict = {b: p for b, p in bips2}

    h_shared = 0
    for b in shared:
        p1 = bips1_dict[b]
        p2 = bips2_dict[b]
        p_avg = (p1 + p2) / 2
        if 0 < p_avg < 1:
            h_shared += -(p_avg * np.log2(p_avg) +
                          (1 - p_avg) * np.log2(1 - p_avg))

    if h1 + h2 == 0:
        return 0.0
    return max(0.0, min(1.0, 1 - 2 * h_shared / (h1 + h2)))


def bootstrap_stats(contree_path, common_set):
    tree = Phylo.read(contree_path, "newick")
    all_t = [t.name for t in tree.get_terminals()]
    for name in all_t:
        if name not in common_set:
            tree.prune(name)
    bs = np.array([
        float(c.confidence)
        for c in tree.find_clades()
        if not c.is_terminal() and c.confidence is not None
    ])
    if len(bs) == 0:
        return {"n": 0, "median": 0, "gt70": 0, "gt90": 0}
    return {
        "n": len(bs),
        "median": np.median(bs),
        "gt70": (bs >= 70).sum() / len(bs) * 100,
        "gt90": (bs >= 90).sum() / len(bs) * 100,
    }


def resolution_index(treefile, common_set):
    tree = Phylo.read(treefile, "newick")
    all_t = [t.name for t in tree.get_terminals()]
    for nm in all_t:
        if nm not in common_set:
            tree.prune(nm)
    n_leaves = len(tree.get_terminals())
    n_int = len(tree.get_nonterminals())
    max_int = n_leaves - 1
    ri = n_int / max_int * 100 if max_int > 0 else 0
    return n_int, max_int, ri


def main():
    parser = argparse.ArgumentParser(
        description="Three-way MLST/wgSNP/cgMLST statistical comparison")
    parser.add_argument("--tree_6gene", required=True)
    parser.add_argument("--tree_7gene", required=True)
    parser.add_argument("--tree_wgsnp", required=True)
    parser.add_argument("--contree_6gene", default=None,
                         help="6-gene .contree (bootstrap); defaults to "
                              "tree_6gene with .treefile->.contree")
    parser.add_argument("--contree_7gene", default=None)
    parser.add_argument("--contree_wgsnp", default=None)
    parser.add_argument("--st_6gene", required=True,
                         help="6-gene st_assignments_final.csv from script 04")
    parser.add_argument("--st_7gene", required=True,
                         help="7-gene st_assignments_final.csv from script 04")
    parser.add_argument("--cgmlst", required=True,
                         help="Wang et al. (2022) cgMLST cluster/ST assignment CSV "
                              "with columns: Accession, Wang_Cluster, Wang_ST")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    def default_contree(path):
        return path.replace(".treefile", ".contree")

    c6 = args.contree_6gene or default_contree(args.tree_6gene)
    c7 = args.contree_7gene or default_contree(args.tree_7gene)
    cw = args.contree_wgsnp or default_contree(args.tree_wgsnp)

    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("THREE-WAY STATISTICAL COMPARISON")
    print("=" * 60)

    # ── Step 1: common taxa ──────────────────────────────────────
    print("\n[1] Checking taxa...")
    s6 = set(t.name for t in Phylo.read(args.tree_6gene, "newick").get_terminals())
    s7 = set(t.name for t in Phylo.read(args.tree_7gene, "newick").get_terminals())
    sw = set(t.name for t in Phylo.read(args.tree_wgsnp, "newick").get_terminals())
    sw -= {"Reference"}
    common = sorted(s6 & s7 & sw)
    common_set = set(common)
    print(f"  6-gene taxa : {len(s6)}")
    print(f"  7-gene taxa : {len(s7)}")
    print(f"  wgSNP taxa  : {len(sw)}")
    print(f"  Common      : {len(common)}")

    # ── Step 2: patristic distances ──────────────────────────────
    print("\n[2] Computing patristic distances...")
    t0 = time.time()
    pn6, pd6, _ = patristic_distances(args.tree_6gene, common_set)
    print(f"  6-gene: {time.time()-t0:.0f}s")
    t0 = time.time()
    pn7, pd7, _ = patristic_distances(args.tree_7gene, common_set)
    print(f"  7-gene: {time.time()-t0:.0f}s")
    t0 = time.time()
    pnw, pdw, _ = patristic_distances(args.tree_wgsnp, common_set)
    print(f"  wgSNP:  {time.time()-t0:.0f}s")

    pd6_r = reorder_matrix(pn6, pd6, common)
    pd7_r = reorder_matrix(pn7, pd7, common)
    pdw_r = reorder_matrix(pnw, pdw, common)

    # ── Step 3: Mantel test ──────────────────────────────────────
    print("\n[3] Mantel test (9,999 permutations)...")
    r_6w, p_6w = mantel_test(pd6_r, pdw_r)
    r_7w, p_7w = mantel_test(pd7_r, pdw_r)
    print(f"  6-gene vs wgSNP: rho={r_6w:.4f}, p={p_6w:.4f}")
    print(f"  7-gene vs wgSNP: rho={r_7w:.4f}, p={p_7w:.4f}")

    # ── Step 4: RF distance ──────────────────────────────────────
    print("\n[4] Robinson-Foulds distance...")
    rf_6w, nrf_6w = rf_distance(args.tree_6gene, args.tree_wgsnp, common_set)
    rf_7w, nrf_7w = rf_distance(args.tree_7gene, args.tree_wgsnp, common_set)
    print(f"  6-gene vs wgSNP: RF={rf_6w}, norm={nrf_6w:.4f}")
    print(f"  7-gene vs wgSNP: RF={rf_7w}, norm={nrf_7w:.4f}")

    # ── Step 5: CID ───────────────────────────────────────────────
    print("\n[5] Clustering Information Distance...")
    cid_6w = cid_distance(args.tree_6gene, args.tree_wgsnp, common_set)
    cid_7w = cid_distance(args.tree_7gene, args.tree_wgsnp, common_set)
    print(f"  6-gene vs wgSNP: CID={cid_6w:.4f}")
    print(f"  7-gene vs wgSNP: CID={cid_7w:.4f}")

    # ── Step 6: Bootstrap support ────────────────────────────────
    print("\n[6] Bootstrap support...")
    bs6 = bootstrap_stats(c6, common_set)
    bs7 = bootstrap_stats(c7, common_set)
    bsw = bootstrap_stats(cw, common_set)
    print(f"  {'Metric':<20} {'6-gene':>10} {'7-gene':>10} {'wgSNP':>10}")
    for k, label in [("n", "Internal nodes"), ("median", "Median BS (%)"),
                      ("gt70", "Nodes BS>=70 (%)"), ("gt90", "Nodes BS>=90 (%)")]:
        fmt = ".0f" if k == "n" else ".1f"
        print(f"  {label:<20} {bs6[k]:>10{fmt}} {bs7[k]:>10{fmt}} {bsw[k]:>10{fmt}}")

    # ── Step 7: Resolution Index ─────────────────────────────────
    print("\n[7] Resolution Index...")
    ni6, mx6, ri6 = resolution_index(args.tree_6gene, common_set)
    ni7, mx7, ri7 = resolution_index(args.tree_7gene, common_set)
    niw, mxw, riw = resolution_index(args.tree_wgsnp, common_set)
    print(f"  6-gene: RI={ri6:.1f}%   7-gene: RI={ri7:.1f}%   wgSNP: RI={riw:.1f}%")

    # ── Step 8: ARI vs cgMLST cluster ───────────────────────────────
    print("\n[8] ARI vs cgMLST cluster...")
    wang = pd.read_csv(args.cgmlst)
    wang["Accession"] = wang["Accession"].astype(str).str.strip()

    your6 = load_st_assignments(Path(args.st_6gene))
    m6 = pd.merge(your6, wang, left_on="strain", right_on="Accession")
    ari_6_cg = adjusted_rand_score(
        m6["Wang_Cluster"].astype(str), m6["ST"].astype(str))
    print(f"  6-gene: {len(your6)} strains with complete locus data; "
          f"{len(m6)} overlap with cgMLST")
    print(f"  6-gene vs cgMLST Cluster ARI: {ari_6_cg:.4f}")

    your7 = load_st_assignments(Path(args.st_7gene))
    m7 = pd.merge(your7, wang, left_on="strain", right_on="Accession")
    ari_7_cg = adjusted_rand_score(
        m7["Wang_Cluster"].astype(str), m7["ST"].astype(str))
    print(f"  7-gene: {len(your7)} strains with complete locus data; "
          f"{len(m7)} overlap with cgMLST")
    print(f"  7-gene vs cgMLST Cluster ARI: {ari_7_cg:.4f}")

    # ── Final table ───────────────────────────────────────────────
    print(f"\n{'='*70}")
    print("FINAL NINE-METRIC TABLE")
    print(f"{'='*70}")
    rows = [
        ("Common taxa for comparison", len(common), len(common), len(common)),
        ("Mantel rho (vs wgSNP)", r_6w, r_7w, "ref"),
        ("RF distance (vs wgSNP)", rf_6w, rf_7w, "ref"),
        ("Norm. RF (vs wgSNP)", nrf_6w, nrf_7w, "ref"),
        ("CID (vs wgSNP)", cid_6w, cid_7w, "ref"),
        ("Internal nodes", bs6["n"], bs7["n"], bsw["n"]),
        ("Resolution Index (%)", ri6, ri7, riw),
        ("Median bootstrap (%)", bs6["median"], bs7["median"], bsw["median"]),
        ("Nodes BS>=70 (%)", bs6["gt70"], bs7["gt70"], bsw["gt70"]),
        ("ARI vs cgMLST cluster", ari_6_cg, ari_7_cg, "n/a"),
    ]
    print(f"  {'Metric':<30} {'6-gene':>12} {'7-gene':>12} {'wgSNP':>12}")
    for label, v6, v7, vw in rows:
        fmt_f = lambda v: f"{v:.4f}" if isinstance(v, float) else str(v)
        print(f"  {label:<30} {fmt_f(v6):>12} {fmt_f(v7):>12} {fmt_f(vw):>12}")

    df_out = pd.DataFrame(
        [(r[0], str(r[1]), str(r[2]), str(r[3])) for r in rows],
        columns=["Metric", "6-gene", "7-gene", "wgSNP"])
    out_csv = out_dir / "comparison_metrics.csv"
    df_out.to_csv(out_csv, index=False)
    print(f"\n✓ Saved: {out_csv}")
    print("  Compare against expected_output/comparison_metrics_expected.csv "
          "to verify your re-run reproduces the published numbers.")


if __name__ == "__main__":
    main()
