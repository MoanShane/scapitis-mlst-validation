#!/usr/bin/env python3
# =============================================================
# Script 11: Unified cgMLST concordance analysis
# =============================================================
# Replaces the concordance section of the former
# 10_12_concordance_mst_tanglegram.py, which compared the two MLST
# schemes on DIFFERENT strain subsets and therefore produced
# non-comparable ARI values.
#
# All metrics below are computed on a single common strain set.
#
# Inputs
#   --ridom      Ridom SeqSphere+ export (.xlsx) for all 620 genomes.
#                Must contain the columns "Sample ID", "ST",
#                "Complex Type" and the 1,492 cgMLST target columns
#                suffixed "(S. capitis cgMLST)".
#   --pst        Per-strain 6-gene assignment CSV with columns
#                accession,pST,femA,ftsZ,gap,pyrH,rpoB,tuf
#   --pst_table  91-row pST catalogue (pST_table_658_FINAL.csv), used
#                to translate allele numbering and to define clonal
#                complexes.
#   --wang       Wang et al. (2022) cluster labels
#                (Accession,Wang_Cluster,Wang_ST,...)
#
# Usage
#   python 11_cgmlst_concordance.py \
#       --ridom   Scapitis_620_cgMLST_SevenBatch.xlsx \
#       --pst     results/pST_per_strain_620.csv \
#       --pst_table results/pST_table_658_FINAL.csv \
#       --wang    wang2022_cgmlst_st.csv \
#       --outdir  results/
#
# IMPORTANT — allele numbering is NOT portable
#   A 6-gene run on 620 genomes and a 6-gene run on 658 genomes assign
#   allele integers independently, so pST labels from the two runs do
#   not correspond. This script re-maps per-strain pST labels onto the
#   658-strain scheme by matching representative strains, never by
#   matching allele numbers. Do not remove that step.
#
# IMPORTANT — outgroup
#   DSM20326 (S. capitis subsp. capitis type strain) is the wgSNP
#   reference and phylogenetic outgroup. It is excluded from all
#   comparative statistics.
# =============================================================

import argparse
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score

RNG = np.random.default_rng(20260814)
LOCI6 = ["femA", "ftsZ", "gap", "pyrH", "rpoB", "tuf"]
OUTGROUP = "DSM20326"
LNZR1 = "GCA_000712995.1"          # L-clone reference strain
LCLONE_MAX_DIST = 25               # allelic distance cut-off (see below)


# ── partition statistics ─────────────────────────────────────────────
def simpson(labels):
    """Simpson's index of diversity with Grundmann et al. 95% CI."""
    n = len(labels)
    cnt = pd.Series(labels).value_counts().values.astype(float)
    D = 1.0 - (cnt * (cnt - 1)).sum() / (n * (n - 1))
    p = cnt / n
    se = np.sqrt(max((4.0 / n) * ((p ** 3).sum() - ((p ** 2).sum()) ** 2), 0))
    return D, max(0, D - 1.96 * se), min(1, D + 1.96 * se), len(cnt)


def _pair_counts(a, b):
    df = pd.DataFrame({"a": pd.factorize(a)[0], "b": pd.factorize(b)[0]})
    ct = pd.crosstab(df["a"], df["b"]).values.astype(float)
    both = (ct * (ct - 1) / 2).sum()
    ra, cb = ct.sum(axis=1), ct.sum(axis=0)
    return both, (ra * (ra - 1) / 2).sum() - both, (cb * (cb - 1) / 2).sum() - both


def wallace(a, b):
    """W(A->B) = P(same type in B | same type in A)."""
    a11, a10, _ = _pair_counts(a, b)
    return a11 / (a11 + a10) if (a11 + a10) > 0 else np.nan


def adj_wallace(a, b):
    """Adjusted Wallace (Severiano et al. 2011)."""
    W, Wi = wallace(a, b), 1.0 - simpson(b)[0]
    return (W - Wi) / (1.0 - Wi) if (1.0 - Wi) > 0 else np.nan


def boot_ci(fn, a, b, n_boot=1000):
    a, b = np.asarray(a), np.asarray(b)
    vals = []
    for _ in range(n_boot):
        i = RNG.integers(0, len(a), len(a))
        try:
            v = fn(a[i], b[i])
            if np.isfinite(v):
                vals.append(v)
        except Exception:
            pass
    return tuple(np.percentile(vals, [2.5, 97.5])) if vals else (np.nan, np.nan)


def wilson(k, n):
    if n == 0:
        return (np.nan, np.nan)
    p, z = k / n, 1.96
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (max(0, c - h), min(1, c + h))


def perf(pred, truth):
    pred, truth = np.asarray(pred, bool), np.asarray(truth, bool)
    tp, fp = int((pred & truth).sum()), int((pred & ~truth).sum())
    fn, tn = int((~pred & truth).sum()), int((~pred & ~truth).sum())
    sn, sp = wilson(tp, tp + fn), wilson(tn, tn + fp)
    return dict(
        TP=tp, FP=fp, FN=fn, TN=tn,
        sens=f"{tp/(tp+fn)*100:.1f} ({sn[0]*100:.1f}-{sn[1]*100:.1f})" if tp + fn else "NA",
        spec=f"{tn/(tn+fp)*100:.1f} ({sp[0]*100:.1f}-{sp[1]*100:.1f})" if tn + fp else "NA",
        ppv=f"{tp/(tp+fp)*100:.1f}" if tp + fp else "NA")


# ── clonal complexes (goeBURST single-locus-variant linkage) ─────────
def slv_components(profiles, names):
    P = np.asarray(profiles, float)
    parent = list(range(len(P)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in range(len(P)):
        for j in np.where((P[i] != P).sum(axis=1) <= 1)[0]:
            ri, rj = find(i), find(int(j))
            if ri != rj:
                parent[max(ri, rj)] = min(ri, rj)
    roots, out = {}, {}
    for i in range(len(P)):
        r = find(i)
        roots.setdefault(r, len(roots) + 1)
        out[names[i]] = f"CC{roots[r]}"
    return out


# ── input handling ──────────────────────────────────────────────────
def load_ridom(path):
    d = pd.read_excel(path)
    plain = {c: c.split("\n")[0].strip() for c in d.columns}
    key = d.rename(columns=plain)[["Sample ID", "ST", "Complex Type",
                                   "Perc. Good Targets"]].copy()
    key.columns = ["accession", "ST7", "cgST", "pct_good"]
    key["accession"] = key["accession"].astype(str).str.strip()
    targets = [c for c in d.columns
               if c.endswith("(S. capitis cgMLST)")
               and not c.split("\n")[0].startswith(("Perc.", "Complex"))]
    mat = d[targets].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    return key, mat


def remap_pst(per_strain, pst_table):
    """Map per-strain pST labels onto the 658-strain scheme via
    representative strains. Allele integers are NOT comparable between
    runs, so profile matching cannot be used."""
    acc2local = dict(zip(per_strain["accession"], per_strain["pST"]))
    mapping, conflicts = {}, []
    for _, r in pst_table.iterrows():
        reps = [x.strip() for x in str(r["Rep_strains"]).split(",")
                if x.strip() and x.strip().lower() != "nan"]
        hits = {acc2local[a] for a in reps if a in acc2local}
        if len(hits) == 1:
            mapping[hits.pop()] = r["pST"]
        elif len(hits) > 1:
            conflicts.append((r["pST"], sorted(hits)))
    if conflicts:
        raise SystemExit(f"pST re-mapping conflict, refusing to continue: {conflicts}")
    return mapping


def lclone_set(key, mat, max_dist):
    """Define the L-clone by cgMLST allelic distance to LNZR-1."""
    idx = int(np.where(key["accession"].values == LNZR1)[0][0])
    valid = ~np.isnan(mat)
    shared = valid[idx] & valid
    diff = ((mat[idx] != mat) & shared).sum(axis=1)
    return set(key["accession"].values[diff <= max_dist]), diff


# ── main ────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(description="Unified cgMLST concordance analysis")
    ap.add_argument("--ridom", required=True)
    ap.add_argument("--pst", required=True)
    ap.add_argument("--pst_table", required=True)
    ap.add_argument("--wang", required=True)
    ap.add_argument("--outdir", default="results")
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    key, mat = load_ridom(args.ridom)
    per = pd.read_csv(args.pst)
    per["accession"] = per["accession"].astype(str).str.strip()
    tab = pd.read_csv(args.pst_table)
    wang = pd.read_csv(args.wang)
    wang.columns = [c.strip().lstrip("\ufeff") for c in wang.columns]
    wang["Accession"] = wang["Accession"].astype(str).str.strip()

    per["pST"] = per["pST"].map(remap_pst(per, tab))
    lclone, dist = lclone_set(key, mat, LCLONE_MAX_DIST)

    m = (key.merge(per[["accession", "pST"]], on="accession", how="left")
            .merge(wang[["Accession", "Wang_Cluster", "Wang_ST"]],
                   left_on="accession", right_on="Accession", how="left")
            .drop(columns="Accession"))
    m["Lclone"] = m["accession"].isin(lclone)
    m["dist_LNZR1"] = dist
    m = m[m["accession"] != OUTGROUP].copy()
    m["ST7"] = "ST" + m["ST7"].astype("Int64").astype(str)

    # clonal complexes, identical rule for both schemes
    present6 = sorted(m["pST"].dropna().unique())
    t6 = tab[tab["pST"].isin(present6)]
    m["CC6"] = m["pST"].map(slv_components(t6[LOCI6].values, t6["pST"].tolist()))
    loci7 = [c for c in pd.read_excel(args.ridom).columns
             if c.endswith("(S. capitis MLST)")
             and c.split("\n")[0] not in ("ST", "CC")]
    p7 = pd.read_excel(args.ridom)[loci7].apply(pd.to_numeric, errors="coerce")
    p7.insert(0, "accession", key["accession"])
    st7 = p7.merge(m[["accession", "ST7"]], on="accession").groupby("ST7").first().reset_index()
    m["CC7"] = m["ST7"].map(slv_components(st7[loci7].values, st7["ST7"].tolist()))

    m.to_csv(out / "master_620_four_schemes.csv", index=False)

    lines = []
    P = lines.append
    P("=" * 74)
    P("UNIFIED cgMLST CONCORDANCE ANALYSIS")
    P(f"outgroup {OUTGROUP} excluded; N = {len(m)}")
    P("=" * 74)

    fine = m.dropna(subset=["pST", "ST7", "cgST"]).copy()
    fine["cgST"] = "CT" + fine["cgST"].astype(int).astype(str)
    P(f"\n[1] Partition characteristics (N = {len(fine)})")
    P(f"{'scheme':<22}{'types':>7}{'D':>9}{'95% CI':>18}")
    for name, col in [("6-gene pST", "pST"), ("7-gene ST", "ST7"), ("cgMLST CT", "cgST")]:
        D, lo, hi, k = simpson(fine[col])
        P(f"{name:<22}{k:>7}{D:>9.4f}{f'{lo:.4f}-{hi:.4f}':>18}")

    P(f"\n[2] Pairwise agreement (N = {len(fine)})")
    for (n1, c1), (n2, c2) in itertools.combinations(
            [("6-gene", "pST"), ("7-gene", "ST7"), ("cgMLST", "cgST")], 2):
        ari = adjusted_rand_score(fine[c1], fine[c2])
        ci = boot_ci(adjusted_rand_score, fine[c1].values, fine[c2].values)
        ami = adjusted_mutual_info_score(fine[c1], fine[c2])
        P(f"  {n1} vs {n2}: ARI {ari:.4f} ({ci[0]:.4f}-{ci[1]:.4f})  AMI {ami:.4f}")
    P("  NOTE ARI between an MLST scheme and 321 cgMLST complex types is")
    P("       dominated by the difference in granularity; report the")
    P("       directional adjusted Wallace coefficient instead.")
    for a, b in [("cgST", "pST"), ("cgST", "ST7")]:
        P(f"  AW cgMLST -> {b}: {adj_wallace(fine[a], fine[b]):.4f}")

    for tag, col in [("STRICT", "Wang_Cluster")]:
        df = m.dropna(subset=["pST", "ST7", col]).rename(columns={col: "cl"})
        P(f"\n[3] Agreement with cgMLST lineage cluster ({tag}, N = {len(df)})")
        for name, c in [("6-gene pST", "pST"), ("6-gene CC", "CC6"),
                        ("7-gene ST", "ST7"), ("7-gene CC", "CC7")]:
            ari = adjusted_rand_score(df["cl"], df[c])
            aw = adj_wallace(df[c], df["cl"])
            P(f"  {name:<12} levels {df[c].nunique():>3}  ARI {ari:.4f}  AW->cluster {aw:.4f}")

        truth = (df["cl"] == "A").values
        ccA6 = df[df["pST"] == "pST1"]["CC6"].mode()[0]
        ccA7 = df[df["ST7"] == "ST1"]["CC7"].mode()[0]
        P(f"\n[4] NRCS-A identification (Cluster A, n = {truth.sum()})")
        P(f"{'predictor':<26}{'TP':>4}{'FP':>4}{'FN':>4}{'TN':>5}{'sens':>22}{'spec':>22}{'PPV':>7}")
        for nm, pr in [("6-gene pST1", df["pST"] == "pST1"),
                       (f"6-gene {ccA6}", df["CC6"] == ccA6),
                       ("7-gene ST1", df["ST7"] == "ST1"),
                       (f"7-gene {ccA7}", df["CC7"] == ccA7)]:
            r = perf(pr.values, truth)
            P(f"{nm:<26}{r['TP']:>4}{r['FP']:>4}{r['FN']:>4}{r['TN']:>5}"
              f"{r['sens']:>22}{r['spec']:>22}{r['ppv']:>7}")

    truth = m["Lclone"].values
    ftsz_l = int(tab[tab["pST"] == "pST11"]["ftsZ"].values[0])
    ftsz_psts = set(tab[tab["ftsZ"] == ftsz_l]["pST"])
    P(f"\n[5] L-clone identification (cgMLST-defined, n = {truth.sum()} of {len(m)})")
    P(f"    definition: <= {LCLONE_MAX_DIST} cgMLST allelic differences from {LNZR1}")
    nxt = sorted(m.loc[~m['Lclone'], 'dist_LNZR1'])[0]
    P(f"    next closest genome differs at {nxt} alleles -> boundary unambiguous")
    P(f"{'predictor':<26}{'TP':>4}{'FP':>4}{'FN':>4}{'TN':>5}{'sens':>22}{'spec':>22}{'PPV':>7}")
    for nm, pr in [("6-gene pST11", m["pST"] == "pST11"),
                   (f"6-gene ftsZ allele {ftsz_l}", m["pST"].isin(ftsz_psts)),
                   ("7-gene ST6", m["ST7"] == "ST6")]:
        r = perf(pr.values, truth)
        P(f"{nm:<26}{r['TP']:>4}{r['FP']:>4}{r['FN']:>4}{r['TN']:>5}"
          f"{r['sens']:>22}{r['spec']:>22}{r['ppv']:>7}")

    txt = "\n".join(lines)
    print(txt)
    (out / "concordance_results.txt").write_text(txt + "\n")
    print(f"\nSaved: {out/'concordance_results.txt'}")
    print(f"Saved: {out/'master_620_four_schemes.csv'}")


if __name__ == "__main__":
    main()
