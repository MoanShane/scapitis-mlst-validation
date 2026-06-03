#!/usr/bin/env python3
# =============================================================
# Script 10: Concordance analysis with cgMLST
# =============================================================
# Computes ARI between 6-gene pST / 7-gene ST / wgSNP clusters
# and Wang 2022 cgMLST cluster and ST assignments.
#
# Usage: python 10_concordance.py --workdir /data/scapitis
# =============================================================

import argparse
from pathlib import Path
import pandas as pd
from sklearn.metrics import adjusted_rand_score
from sklearn.cluster import AgglomerativeClustering
import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--pst_table", default=None)
    parser.add_argument("--cgmlst_table", default=None)
    parser.add_argument("--output", default=None)
    args = parser.parse_args()

    workdir    = Path(args.workdir)
    pst_csv    = Path(args.pst_table) if args.pst_table else \
                 workdir / "results" / "pST_table_658_FINAL.csv"
    cgmlst_csv = Path(args.cgmlst_table) if args.cgmlst_table else \
                 workdir / "wang2022_cgmlst_st.csv"
    out_dir    = Path(args.output) if args.output else \
                 workdir / "results"
    out_dir.mkdir(parents=True, exist_ok=True)

    wang = pd.read_csv(str(cgmlst_csv))
    wang["Accession"] = wang["Accession"].str.strip()

    # Encode cgMLST cluster and ST as integers
    cluster_enc = {c:i for i,c in
                   enumerate(wang["Wang_Cluster"].unique())}
    wang["cgCluster_int"] = wang["Wang_Cluster"].map(cluster_enc)
    wang["cgST_int"]      = wang["Wang_ST"].astype(int)

    print("=== cgMLST Concordance Analysis (ARI) ===")
    results = {}

    # 6-gene pST vs cgMLST
    # Expand pST table to per-strain
    pst_df = pd.read_csv(str(pst_csv))
    pst_expanded = []
    for _, row in pst_df.iterrows():
        for s in str(row["strains"]).split(";")[::-1]:
            pst_expanded.append({
                "Accession": s.strip(),
                "pST": row["pST"]
            })
    pst_strain = pd.DataFrame(pst_expanded)

    merged6 = pd.merge(pst_strain, wang, on="Accession")
    if len(merged6) > 0:
        pst_int = pd.factorize(merged6["pST"])[0]
        ari_6_cl = adjusted_rand_score(
            merged6["cgCluster_int"], pst_int)
        ari_6_st = adjusted_rand_score(
            merged6["cgST_int"], pst_int)
        print(f"  6-gene pST vs cgMLST cluster: ARI = {ari_6_cl:.4f}")
        print(f"  6-gene pST vs cgMLST ST:      ARI = {ari_6_st:.4f}")
        results["ARI_6gene_vs_cgCluster"] = ari_6_cl
        results["ARI_6gene_vs_cgST"]      = ari_6_st

    # 7-gene ST vs cgMLST (Wang_ST is 7-gene ST)
    ari_7_cl = adjusted_rand_score(
        wang["cgCluster_int"], wang["cgST_int"])
    print(f"  7-gene ST vs cgMLST cluster:  ARI = {ari_7_cl:.4f}")
    print(f"  (7-gene ST vs cgMLST ST = N/A: circular)")
    results["ARI_7gene_vs_cgCluster"] = ari_7_cl

    # Save
    out_csv = out_dir / "concordance_ari.csv"
    pd.DataFrame(
        [(k, v) for k,v in results.items()],
        columns=["Metric","ARI"]
    ).to_csv(str(out_csv), index=False)
    print(f"\n✓ Saved: {out_csv}")


if __name__ == "__main__":
    main()


# =============================================================
# Script 11: Minimum Spanning Tree figure
# =============================================================
# See full implementation in the methods scripts above (08_statistics).
# The MST figure generation script is provided separately as
# 11_mst_figure.py due to its complexity; the key steps are:
#   1. Compute pairwise allele-profile distances between pSTs
#   2. Build MST using SciPy minimum_spanning_tree (Kruskal)
#   3. Compute node layout using NetworkX Kamada-Kawai
#   4. Plot with matplotlib (node size ~ strain count, log scale)
# Refer to the published figure and GitHub repository for the
# complete annotated figure generation script.
# =============================================================


# =============================================================
# Script 12: Tanglegram (6-gene vs wgSNP)
# =============================================================
# Generates tanglegram comparing 6-gene MLST tree (91 unique
# pSTs) with wgSNP tree. Computes entanglement score.
#
# Usage: python 12_tanglegram.py \
#            --tree_6gene trees/tree_6gene_658.treefile \
#            --tree_wgsnp wgsnp/core/wgsnp_620.treefile \
#            --pst_table results/pST_table_658_FINAL.csv \
#            --output results/tanglegram.png
# =============================================================

import sys

def tanglegram_entanglement(order1, order2):
    """
    Compute entanglement score = crossings / max_crossings.
    order1, order2: lists of leaf names in tip order for each tree.
    """
    n = len(order1)
    pos1 = {name: i for i,name in enumerate(order1)}
    pos2 = {name: i for i,name in enumerate(order2)}

    common = [t for t in order1 if t in pos2]
    mapped = [(pos1[t], pos2[t]) for t in common]

    crossings = 0
    for i in range(len(mapped)):
        for j in range(i+1, len(mapped)):
            a1, a2 = mapped[i]
            b1, b2 = mapped[j]
            if (a1 < b1) != (a2 < b2):
                crossings += 1

    max_crossings = len(mapped) * (len(mapped) - 1) / 2
    return crossings / max_crossings if max_crossings > 0 else 0


if __name__ == "__main__":
    # Example usage
    example_order1 = ["A","B","C","D","E"]
    example_order2 = ["A","C","B","E","D"]
    score = tanglegram_entanglement(example_order1, example_order2)
    print(f"Example entanglement: {score:.4f}")
    print("Full tanglegram script available at GitHub repository.")
