#!/usr/bin/env python3
# =============================================================
# Script 09: Provisional sequence type (pST) definition
# =============================================================
# Defines 91 pSTs from unique 6-gene allele combinations
# across the 658-strain dataset.
#
# Usage: python 09_provisional_st.py --workdir /data/scapitis
# =============================================================

import argparse
from pathlib import Path
from collections import defaultdict, Counter
import csv
import pandas as pd
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--concat_658",
                        default=None,
                        help="Path to 658-strain concat FASTA")
    parser.add_argument("--cgmlst_table",
                        default=None,
                        help="Wang 2022 cgMLST assignment table (CSV)")
    parser.add_argument("--output", default=None)
    args = parser.parse_args()

    workdir    = Path(args.workdir)
    concat_658 = Path(args.concat_658) if args.concat_658 else \
                 workdir / "concat_6gene" / "concat_6gene_658.fasta"
    cgmlst_csv = Path(args.cgmlst_table) if args.cgmlst_table else \
                 workdir / "wang2022_cgmlst_st.csv"
    out_dir    = Path(args.output) if args.output else \
                 workdir / "results"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Gene positions in 658-strain concat (from partition file)
    GENE_POS = {
        "femA": (0,   492),
        "ftsZ": (492, 1010),
        "gap":  (1010,1531),
        "pyrH": (1531,2005),
        "rpoB": (2005,2539),
        "tuf":  (2539,3032),
    }
    GENES = ["femA","ftsZ","gap","pyrH","rpoB","tuf"]

    # Read sequences
    print("Reading 658-strain concat...")
    strain_seqs = {}
    for rec in SeqIO.parse(str(concat_658), "fasta"):
        seq = str(rec.seq).upper()
        strain_seqs[rec.id] = {
            g: seq[s:e] for g,(s,e) in GENE_POS.items()}
    print(f"  Strains: {len(strain_seqs)}")

    # Build allele catalogue (frequency-based numbering)
    print("Building allele catalogue...")
    allele_cat = {}
    for gene in GENES:
        cnt = Counter(strain_seqs[s][gene]
                      for s in strain_seqs)
        allele_cat[gene] = {
            seq: rank
            for rank,(seq,_) in enumerate(cnt.most_common(), 1)}
        print(f"  {gene}: {len(cnt)} unique alleles")

    # Assign allele profiles
    strain_profile = {
        s: tuple(allele_cat[g][seqs[g]] for g in GENES)
        for s,seqs in strain_seqs.items()
    }

    profile_strains = defaultdict(list)
    for s,p in strain_profile.items():
        profile_strains[p].append(s)

    # Load cgMLST data
    wang = pd.read_csv(str(cgmlst_csv))
    wang["Accession"] = wang["Accession"].str.strip()
    nrcs_a  = set(wang[wang["Wang_Cluster"]=="A"]["Accession"])
    # NOTE: Wang et al. (2022) assigned no cluster label to the L-clone,
    # so this set is empty with the published reference file. The L-clone
    # is defined instead in 11_cgmlst_concordance.py, by cgMLST allelic
    # distance to LNZR-1 (GCA_000712995.1).
    l_clone = set(wang[wang["Wang_ST"]==6]["Accession"])

    # Sort profiles by cgMLST priority
    def sort_key(item):
        prof, strains = item
        n_A   = sum(1 for s in strains if s in nrcs_a)
        n_L   = sum(1 for s in strains if s in l_clone)
        lnzr  = "GCA_000712995.1" in strains
        return (-n_A, -(n_L + int(lnzr)), -len(strains))

    sorted_profiles = sorted(profile_strains.items(),
                             key=sort_key)

    # Build pST table
    print("\nBuilding pST catalogue...")
    rows = []
    per_strain_rows = []
    for pst_num, (profile, strains) in \
            enumerate(sorted_profiles, 1):
        wang_sub = wang[wang["Accession"].isin(strains)]
        n_over   = len(wang_sub)

        cg_top = "—"; st7_top = "—"
        if n_over > 0:
            cg_cnt  = wang_sub["Wang_Cluster"].value_counts()
            st7_cnt = wang_sub["Wang_ST"].value_counts()
            cg_top  = cg_cnt.index[0]
            st7_top = f"ST{st7_cnt.index[0]}"

        n_A  = sum(1 for s in strains if s in nrcs_a)
        lnzr = "GCA_000712995.1" in strains
        lin  = ("NRCS-A" if n_A == n_over and n_over > 0
                else "L-clone" if lnzr else "")

        rows.append({
            "pST":           f"pST{pst_num}",
            **{f"allele_{g}": str(profile[i])
               for i,g in enumerate(GENES)},
            "N_total":        len(strains),
            "N_Wang_overlap": n_over,
            "cgMLST_cluster": cg_top,
            "7gene_ST":       st7_top,
            "Lineage":        lin,
            # Representative strains for display in Supplementary Table S1.
            # TRUNCATED ON PURPOSE - never use this column for analysis.
            "Rep_strains":    ",".join(sorted(strains)[:5]),
        })
        # Full per-strain record, written to a separate file below.
        for s in sorted(strains):
            per_strain_rows.append({
                "accession": s,
                "pST":       f"pST{pst_num}",
                **{g: profile[i] for i, g in enumerate(GENES)},
            })

    df = pd.DataFrame(rows)
    out_csv = out_dir / "pST_table_658_FINAL.csv"
    df.to_csv(str(out_csv), index=False)

    # ---------------------------------------------------------------
    # Complete per-strain assignment table.
    #
    # This file - NOT the Rep_strains column above - is the input for
    # every downstream concordance calculation. An earlier version of
    # this pipeline expanded the truncated representative-strain list
    # instead, which silently reduced the analysis set and made the
    # six-gene and seven-gene schemes non-comparable.
    # ---------------------------------------------------------------
    per_strain = pd.DataFrame(per_strain_rows)
    per_csv = out_dir / "pST_per_strain_658.csv"
    per_strain.to_csv(str(per_csv), index=False)
    assert len(per_strain) == len(strain_seqs), (
        f"per-strain table has {len(per_strain)} rows but "
        f"{len(strain_seqs)} strains were read - refusing to continue")
    print(f"  Per-strain assignments: {len(per_strain)} strains")

    print(f"  Total pSTs: {len(df)}")
    print(f"  pST1 (NRCS-A): "
          f"{df[df['pST']=='pST1']['N_total'].values[0]} strains")
    print(f"\n✓ Saved: {out_csv}")
    print(f"✓ Saved: {per_csv}")


if __name__ == "__main__":
    main()
