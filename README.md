# scapitis-mlst-validation

Six-gene versus seven-gene MLST comparison pipeline for *Staphylococcus capitis*

## Overview

Source code for:

> Tsai M-S, Ling TZ. *Six-Gene versus Seven-Gene MLST for Staphylococcus capitis*:
> Performance Comparison Against Whole-Genome SNP and Core-Genome MLST, with
> Detection of Putative L-Clone-Affiliated Isolates in Taiwan. *Journal of
> Clinical Microbiology* (submitted).

This pipeline compares the Song et al. (2019) six-gene MLST scheme for
*Staphylococcus capitis* with the Wang et al. (2025) seven-gene scheme,
using whole-genome SNP (wgSNP) analysis and core-genome MLST (cgMLST) as
reference standards, across 658 genomes (620 public sequences + 29 YGH
clinical isolates + 9 Song 2019 reference strains).

The two schemes are **complementary rather than hierarchically ranked**.
The six-gene scheme corresponds more closely to the wgSNP reference
phylogeny; the seven-gene scheme agrees more closely with cgMLST lineage
clusters. Which scheme appears "better" depends on the reference standard
chosen, so both are reported.

---

## ?? Read before running: three ways this analysis can silently go wrong

### 1. Comparing the two schemes on different strain subsets

Concordance statistics are only comparable when computed on the **same
set of strains**. An earlier version of this pipeline computed six-gene
concordance on one subset and seven-gene concordance on another, which
made the two ARI values non-comparable and produced misleading results.

`11_cgmlst_concordance.py` enforces a single common strain set for every
metric. Do not remove that behaviour.

### 2. Allele numbering is not portable between runs

A six-gene run on 620 genomes and a six-gene run on 658 genomes assign
allele integers **independently**, because numbering is frequency-ranked
within whichever dataset was processed. The same allele profile can
therefore receive different pST labels in the two runs.

`11_cgmlst_concordance.py` re-maps pST labels by matching
**representative strains**, never by matching allele numbers, and aborts
if it detects a conflict. Never merge outputs from two runs by allele
number.

### 3. Truncated strain lists

`pST_table_658_FINAL.csv` contains a `Rep_strains` column listing at most
five representative strains per pST. That column exists for display in
Supplementary Table S1 only. Expanding it into a per-strain table
silently discards most of the dataset.

Use `pST_per_strain_658.csv`, written by `09_provisional_st.py`, for all
downstream analysis.

### 4. Missing data at MLST loci (BLAST-based pipeline only)

Some genomes have genuinely incomplete locus coverage at one or more of
the seven MLST genes ??most notably 41 strains in the `ERR3378xxx`
series, which are missing 6 of the 7 seven-gene loci, and DSM20326,
which is missing `rluB`. **A strain with missing data at any locus must
be excluded from ST assignment for that scheme, not assigned a spurious
allele derived from an empty sequence.** Earlier iterations conflated
"0 bp after stripping gaps/N" with "a real zero-length allele", which
fragmented true alleles into many spurious ones. `04_extract_mlst.py`
and `08_statistics.py` implement the correct exclusion logic (see the
`has_missing_data` / `missing_loci` columns).

**Note:** the seven-gene sequence types reported in the manuscript were
not produced by this BLAST-based route. They come from allele calling
performed from genome assemblies within MBioSEQ Ridom Typer, which
returned a complete seven-locus profile for all 620 genomes and
therefore recovered the 42 strains that this pipeline excludes. The
exclusion logic above remains correct for the BLAST route.

### 5. pandas and the string "NA"

`pandas.read_csv()` treats the literal string `"NA"` as a missing-value
marker by default and silently converts it to `NaN`. If you filter ST
assignment tables on the string `"NA"`, read with
`keep_default_na=False, na_values=[]` first, or the filter will silently
match nothing.

---

## Repository structure

```
scapitis-mlst-validation/
??? README.md
??? LICENSE
??? environment.yml                  # conda environment (pinned versions)
??? requirements_python.txt          # pip packages (pinned)
??? 01_download_sra.sh               # Download reads from SRA
??? 02_assemble_genomes.sh           # de novo assembly (Shovill/SPAdes)
??? 03_qc_check.sh                   # Assembly QC
??? 04_extract_mlst.py               # BLAST-based MLST allele extraction
??                                      (6-gene and 7-gene; handles missing data)
??? 05_align_concat.sh               # MAFFT alignment + trimAl + concat
??? 06_wgsnp_snippy.sh               # Whole-genome SNP analysis (Snippy)
??? 07_phylogeny_iqtree.sh           # IQ-TREE2 ML trees
??? 08_statistics.py                 # Mantel / RF / CID / PIS / bootstrap
??? 09_provisional_st.py             # pST definition; writes the pST catalogue
??                                      AND the per-strain assignment table
??? 11_cgmlst_concordance.py         # Unified cgMLST concordance analysis:
??                                      Simpson's D, ARI, adjusted Wallace,
??                                      clonal complexes, NRCS-A and L-clone
??                                      identification performance
??? docs/
??  ??? pipeline_overview.md
??? expected_output/
    ??? comparison_metrics_expected.csv
```

The cgMLST step itself is performed in MBioSEQ Ridom Typer, a commercial
GUI application, and is therefore not scripted. The manuscript records
the software version and scheme; `11_cgmlst_concordance.py` takes the
Ridom Typer Excel export as input.

---

## Dependencies

### Software

| Software                 | Version         | Use                                 |
| ------------------------ | --------------- | ----------------------------------- |
| SRA Toolkit              | 3.0.5           | Download SRA reads                  |
| Shovill                  | 1.1.0           | Genome assembly                     |
| SPAdes                   | 3.15            | Assembler (via Shovill)             |
| BLAST+                   | 2.13.0          | MLST allele extraction              |
| MAFFT                    | 7.490           | Multiple sequence alignment         |
| trimAl                   | 1.4.1           | Alignment trimming                  |
| Snippy                   | 4.6.0           | Whole-genome SNP analysis           |
| IQ-TREE2                 | 2.2.0           | Phylogenetic inference              |
| MBioSEQ Ridom Typer      | 12.0.5 (2026/05)| cgMLST typing (commercial, GUI)     |
| Python                   | 3.9             | Statistical analysis                |

MBioSEQ Ridom Typer (Ridom GmbH, a Bruker company, M羹nster, Germany) was
formerly marketed as Ridom SeqSphere+. The *S. capitis* cgMLST scheme
used here comprises 1,492 targets.

### Python packages (pinned)

```
biopython==1.79
numpy==1.21.6
scipy==1.7.3
scikit-learn==1.0.2
pandas==1.3.5
matplotlib==3.5.3
networkx==2.8.8
openpyxl==3.0.10
```

**Why pinned, not `>=`:** an upgrade-sensitive default in
`pandas.read_csv()` (silent string-to-`NaN` conversion of `"NA"`,
described above) once caused a result to regress with no error or
warning. `scikit-learn`'s `adjusted_rand_score` has also had minor
behavioural changes across versions. Pinning exact versions is the only
way to guarantee that a re-run reproduces the published numbers.

---

## Installation

```bash
git clone https://github.com/MoanShane/scapitis-mlst-validation.git
cd scapitis-mlst-validation

conda env create -f environment.yml
conda activate scapitis_mlst
```

---

## Usage

```bash
WORKDIR=/data/scapitis
mkdir -p ${WORKDIR}/{reads,assemblies,mlst_6gene,mlst_7gene,wgsnp,trees,results}

bash 01_download_sra.sh    ${WORKDIR}
bash 02_assemble_genomes.sh ${WORKDIR}
bash 03_qc_check.sh        ${WORKDIR}

python 04_extract_mlst.py --assembly_dir ${WORKDIR}/assemblies \
    --scheme 6gene --ref_dir reference_alleles \
    --output_dir ${WORKDIR}/mlst_6gene
python 04_extract_mlst.py --assembly_dir ${WORKDIR}/assemblies \
    --scheme 7gene --ref_dir reference_alleles \
    --output_dir ${WORKDIR}/mlst_7gene

bash 05_align_concat.sh ${WORKDIR} 6gene
bash 05_align_concat.sh ${WORKDIR} 7gene
bash 06_wgsnp_snippy.sh ${WORKDIR}
bash 07_phylogeny_iqtree.sh ${WORKDIR}

python 08_statistics.py \
    --tree_6gene ${WORKDIR}/trees/tree_6gene_620.treefile \
    --tree_7gene ${WORKDIR}/trees/tree_7gene_620.treefile \
    --tree_wgsnp ${WORKDIR}/wgsnp/wgsnp_620.treefile \
    --st_6gene   ${WORKDIR}/mlst_6gene/st_assignments_final.csv \
    --st_7gene   ${WORKDIR}/mlst_7gene/st_assignments_final.csv \
    --output     ${WORKDIR}/results/

python 09_provisional_st.py --workdir ${WORKDIR}

# cgMLST typing is performed separately in MBioSEQ Ridom Typer;
# export the project to .xlsx, then:
python 11_cgmlst_concordance.py \
    --ridom     Scapitis_620_cgMLST_SevenBatch.xlsx \
    --pst       ${WORKDIR}/results/pST_per_strain_658.csv \
    --pst_table ${WORKDIR}/results/pST_table_658_FINAL.csv \
    --wang      wang2022_cgmlst_st.csv \
    --outdir    ${WORKDIR}/results/
```

### Verifying your re-run

Compare your output against
`expected_output/comparison_metrics_expected.csv`, which lists the value
of every metric **together with the N on which it was computed**. If your
N differs, the metric is not comparable with the published value, whatever
the number itself looks like.

Key values to check:

| Metric                                | Six-gene | Seven-gene | N   |
| ------------------------------------- | -------- | ---------- | --- |
| Mantel ? vs wgSNP                     | 0.9449   | 0.9334     | 619 |
| Simpson's D                           | 0.598    | 0.608      | 613 |
| ARI vs cgMLST cluster                 | 0.886    | 0.934      | 466 |
| Adjusted Wallace ??cluster            | 0.927    | 0.997      | 466 |
| NRCS-A sensitivity (clonal complex)   | 100.0%   | 99.4%      | 466 |
| L-clone sensitivity                   | 100.0%   | 100.0%     | 619 |
| L-clone specificity                   | 98.0%    | 100.0%     | 619 |

---

## Input data

- Public genomes: BioProject PRJNA493527 and additional GenBank accessions
  (Supplementary Table S2)
- YGH clinical isolates: GenBank accessions PZ469898?Z470071
- Reference genome (wgSNP): *S. capitis* subsp. *capitis* DSM20326
  (GCF_040739495.1)
- L-clone anchor strain: LNZR-1 (GCA_000712995.1)

---

## Output files

| File                                    | Description                                             |
| --------------------------------------- | ------------------------------------------------------- |
| `concat_6gene_658.fasta`                | 6-gene alignment, 658 strains, 3,032 bp untrimmed       |
| `concat_7gene_620.fasta`                | 7-gene alignment, 2,780 bp                              |
| `st_assignments_final.csv`              | Per-strain ST with `has_missing_data` / `missing_loci`  |
| `tree_6gene_620.treefile`               | ML tree, 6-gene scheme                                  |
| `tree_7gene_620.treefile`               | ML tree, 7-gene scheme                                  |
| `wgsnp_620.treefile`                    | ML tree, wgSNP (619 taxa after outgroup exclusion)      |
| `results/comparison_metrics.csv`        | Topological and bootstrap metrics                       |
| `results/pST_table_658_FINAL.csv`       | 91 provisional sequence types (display table)           |
| `results/pST_per_strain_658.csv`        | **Per-strain pST ??use this for analysis**              |
| `results/concordance_results.txt`       | Full concordance report                                 |
| `results/master_620_four_schemes.csv`   | Per-strain pST, ST, cgST, cluster and L-clone status    |

---

## Citation

> Tsai M-S, Ling TZ. *Six-Gene versus Seven-Gene MLST for Staphylococcus
> capitis*: Performance Comparison Against Whole-Genome SNP and
> Core-Genome MLST, with Detection of Putative L-Clone-Affiliated
> Isolates in Taiwan. *Journal of Clinical Microbiology* (submitted).
> DOI: [to be added on acceptance]

Please cite the specific code version via its GitHub Release tag
(`v2.0-jcm-revision`) rather than the `main` branch.

---

## License

MIT License. See [LICENSE](LICENSE).

---

## Contact

Moan-Shane Tsai, MD
Division of Infectious Diseases
Yuan's General Hospital
Kaohsiung, Taiwan
