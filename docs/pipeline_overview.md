# Pipeline Overview

This document describes each step of the *S. capitis* 6-gene vs 7-gene
MLST validation pipeline in more detail than the top-level README, and
documents two non-obvious pitfalls that were encountered and fixed
during development. If you adapt any part of this pipeline, please read
the "Missing-data handling" section carefully — it affects step 04 and
step 08 directly, and any downstream concordance/ARI calculation
indirectly.

## Step-by-step description

| Step | Script | What it does |
|---|---|---|
| 01 | `01_download_sra.sh` | Downloads raw reads for all public strains (BioProject PRJNA493527 and additional accessions) via SRA Toolkit |
| 02 | `02_assemble_genomes.sh` | De novo assembly with Shovill (SPAdes backend) |
| 03 | `03_qc_check.sh` | Assembly quality control (contig count, N50, total length) |
| 04 | `04_extract_mlst.py` | BLAST-based extraction of the 6 (or 7) MLST loci from each assembly, followed by allele numbering and ST assignment |
| 05 | `05_align_concat.sh` | MAFFT alignment per locus, trimAl trimming, concatenation into a single multi-gene alignment |
| 06 | `06_wgsnp_snippy.sh` | Whole-genome SNP calling against the DSM20326 reference using Snippy, followed by core-genome alignment |
| 07 | `07_phylogeny_iqtree.sh` | Maximum-likelihood tree inference with IQ-TREE2 (ModelFinder + 1,000 ultrafast bootstrap) for all three datasets (6-gene, 7-gene, wgSNP) |
| 08 | `08_statistics.py` | Nine-metric three-way comparison (Mantel, RF, normalised RF, CID, bootstrap support, resolution index, ARI vs cgMLST) |
| 09 | `09_provisional_st.py` | Defines provisional sequence types (pST) from unique 6-gene allele combinations across all 658 strains |
| 10–12 | `10_12_concordance_mst_tanglegram.py` | cgMLST concordance heatmap, minimum spanning tree figure, and 6-gene-vs-wgSNP tanglegram |

## Missing-data handling

### The problem

Not every public genome has complete sequence coverage at every one of
the 6 (or 7) MLST loci. In this dataset, 41 strains (the `ERR3378xxx`
series) are missing 6 of the 7 seven-gene loci, and the reference strain
DSM20326 is missing `rluB`. This is genuine biological/technical missing
data — incomplete assembly coverage at that specific locus — not a real
"short allele."

### What went wrong (and how it was fixed)

An earlier version of the seven-gene ST-assignment logic treated *any*
sequence record, including one consisting entirely of gap/N characters
(i.e., 0 bp of real sequence after stripping), as a valid — if
unusually short — allele. Feeding a 0-bp "allele" into a multi-locus
profile is catastrophic: it fragments what should be a single shared
allele into many spurious distinct alleles, one for every strain with
that pattern of missing loci, and inflates the apparent number of
sequence types far beyond what the data actually support.

The concrete symptom: the seven-gene ARI against the Wang et al. (2022)
cgMLST cluster assignment collapsed from an expected ~0.93 down to
~0.75–0.76 — an 18-percentage-point drop that, on inspection, traced
directly back to those 41+1 strains being assigned spurious "novel"
alleles instead of being excluded.

### The fix

`04_extract_mlst.py` in this repository implements the correct logic:

1. A locus is only recorded as present for a strain if BLAST found a
   valid hit AND the extraction succeeded. Absence is absence — it is
   never backfilled with an empty-string placeholder that could later
   be mistaken for a real allele.
2. The allele catalogue for each gene (the mapping from unique sequence
   to allele number) is built using **only strains with complete
   7/7-locus data**. A missing locus must never define or match an
   allele.
3. Any strain missing ≥1 locus is recorded in the output CSV with
   `ST = "NA"` and `has_missing_data = "yes"`, with the specific
   missing loci listed in `missing_loci`. It is excluded from ST
   assignment entirely, not force-assigned a number.

`08_statistics.py` then loads these CSVs via `load_st_assignments()`,
which explicitly drops `has_missing_data == "yes"` rows before any ARI
or concordance calculation.

### A second, independent pitfall: pandas and the string "NA"

Even with the correct upstream logic, a downstream filter can silently
fail. `pandas.read_csv()` has built-in logic that treats certain literal
strings — including `"NA"` — as missing-value markers, and converts them
to `NaN` on read, *by default*, regardless of what your CSV actually
contains in that cell.

This means a filter written naively as:

```python
df = pd.read_csv("st_assignments_final.csv")
df = df[df["ST"] != "NA"]   # BUG: this does nothing useful
```

will not work as intended. Because pandas already converted the string
`"NA"` to `NaN` during the read, the column no longer contains the
string `"NA"` anywhere — it contains `NaN`. The comparison
`NaN != "NA"` evaluates to `True` for every row, so the filter matches
*everything*, including the rows you meant to exclude. No error is
raised; the bug is completely silent.

The fix used throughout this pipeline:

```python
df = pd.read_csv(
    "st_assignments_final.csv",
    keep_default_na=False,   # do not auto-convert "NA"/"NaN"/etc. to NaN
    na_values=[],             # do not treat any additional strings as NaN
)
df = df[df["ST"] != "NA"]    # now this comparison is against the literal string
```

If you write any new code that reads these CSVs, use this pattern.

## Reproducibility checklist

Before reporting results from a re-run of this pipeline, confirm:

- [ ] `st_assignments_final.csv` for both schemes contains a
      `has_missing_data` column, and the count of `yes` rows matches
      what you expect given known incomplete strains in your dataset
      (41 + 1 = 42 for the published 620-strain dataset).
- [ ] Any code reading these CSVs uses
      `keep_default_na=False, na_values=[]`.
- [ ] `comparison_metrics.csv` from step 08 matches
      `expected_output/comparison_metrics_expected.csv` within rounding
      tolerance. In particular, seven-gene ARI vs cgMLST cluster should
      be in the 0.92–0.93 range, not 0.75–0.80 — a value in the latter
      range indicates the missing-data exclusion did not take effect.
