# scapitis-mlst-validation

Six-gene MLST validation pipeline for *Staphylococcus capitis*

# *S. capitis* 6-gene vs 7-gene MLST Validation Pipeline

## Overview

Source code for:

> Tsai M-S, Ling TZ. *Six-Gene versus Seven-Gene MLST for Staphylococcus capitis*:
> Performance Validation Against Whole-Genome SNP and Core-Genome MLST, with
> Detection of Putative L-Clone-Affiliated Isolates in Taiwan. *Journal of
> Clinical Microbiology* (submitted).

This pipeline validates the Song et al. (2019) six-gene MLST scheme for
*Staphylococcus capitis* against the Wang et al. (2025) seven-gene scheme,
whole-genome SNP (wgSNP) analysis, and core-genome MLST, using 658 genomes
(620 public sequences + 29 YGH clinical isolates + 9 Song 2019 reference
strains).

---

## ?? Important: missing-data handling (read before running)

This dataset contains genomes with genuinely incomplete locus coverage at
one or more of the 7 MLST genes (most notably 41 strains in the
`ERR3378xxx` series, which are missing 6 of the 7 seven-gene loci, and
DSM20326, which is missing `rluB`). **A strain with missing data at any
locus must be excluded from ST assignment at that scheme, not assigned a
spurious allele derived from an empty sequence.** Earlier iterations of
this pipeline conflated "0 bp after stripping gaps/N" with "a real
zero-length allele," which fragmented true alleles into many spurious
ones and collapsed the seven-gene ARI vs cgMLST cluster from the expected
~0.93 down to ~0.75. `04_extract_mlst.py` and `08_statistics.py` in this
repository already implement the correct exclusion logic (see
`has_missing_data` / `missing_loci` columns in the ST assignment output).
If you adapt this code, preserve that behaviour.

A second, unrelated pitfall: `pandas.read_csv()` treats the literal
string `"NA"` as a missing-value marker by default and silently converts
it to `NaN`. If you filter ST assignment tables on the string `"NA"`,
read with `keep_default_na=False, na_values=[]` first, or the filter will
silently match nothing and missing-data strains will leak back into your
analysis.

---

## Repository Structure

```
scapitis-mlst-validation/
??? README.md
??? LICENSE
??? environment.yml                       # conda environment (pinned versions)
??? requirements_python.txt               # pip-installable Python packages (pinned)
??? 01_download_sra.sh                    # Download reads from SRA
??? 02_assemble_genomes.sh                # de novo assembly (Shovill/SPAdes)
??? 03_qc_check.sh                        # Assembly QC
??? 04_extract_mlst.py                    # BLAST-based MLST allele extraction
??                                           (6-gene and 7-gene; handles missing data)
??? 05_align_concat.sh                    # MAFFT alignment + trimAl + concat
??? 06_wgsnp_snippy.sh                    # Whole-genome SNP analysis (Snippy)
??? 07_phylogeny_iqtree.sh                # IQ-TREE2 ML trees
??? 08_statistics.py                      # Mantel/RF/CID/ARI/PIS/Bootstrap
??                                           (9-metric three-way comparison)
??? 09_provisional_st.py                  # pST definition and catalogue
??? 10_12_concordance_mst_tanglegram.py   # cgMLST concordance + MST + tanglegram
??? docs/
??  ??? pipeline_overview.md              # Detailed pipeline description
??? expected_output/
    ??? comparison_metrics_expected.csv   # Reference values reported in the paper,
                                             for verifying a successful re-run
```

**Note:** steps 10??2 (cgMLST concordance, minimum spanning tree figure,
tanglegram) are implemented in a single combined script,
`10_12_concordance_mst_tanglegram.py`, run with a `--step` flag to select
which output to generate (see Usage below).

---

## Dependencies

### Software (versions used in this study)

| Software    | Version | Use                          |
| ----------- | ------- | ----------------------------- |
| SRA Toolkit | 3.0.5   | Download SRA reads            |
| Shovill     | 1.1.0   | Genome assembly                |
| SPAdes      | 3.15    | Assembler (via Shovill)        |
| BLAST+      | 2.13.0  | MLST allele extraction         |
| MAFFT       | 7.490   | Multiple sequence alignment    |
| trimAl      | 1.4.1   | Alignment trimming             |
| Snippy      | 4.6.0   | Whole-genome SNP analysis      |
| IQ-TREE2    | 2.2.0   | Phylogenetic inference         |
| Python      | 3.9     | Statistical analysis           |

### Python packages (pinned ??see note on version sensitivity below)

```
biopython==1.79
numpy==1.21.6
scipy==1.7.3
scikit-learn==1.0.2
pandas==1.3.5
matplotlib==3.5.3
networkx==2.8.8
```

**Why pinned, not `>=`:** during development, an upgrade-sensitive default
in `pandas.read_csv()` (silent string-to-`NaN` conversion of `"NA"`,
described above) caused a result to silently regress without any error
or warning. `scikit-learn`'s `adjusted_rand_score` and
`AgglomerativeClustering` have also had minor behavioural changes across
minor versions. Pinning exact versions is the only way to guarantee a
re-run reproduces the published numbers; `environment.yml` and
`requirements_python.txt` both pin exact versions for this reason.

---

## Installation

```bash
# Clone repository
git clone https://github.com/MoanShane/scapitis-mlst-validation.git
cd scapitis-mlst-validation

# Create conda environment (recommended ??pins both Python and external tools)
conda env create -f environment.yml
conda activate scapitis_mlst

# Or install Python packages only (you must install BLAST+, MAFFT, trimAl,
# Snippy, IQ-TREE2, Shovill, and SRA Toolkit separately if using this route)
pip install -r requirements_python.txt
```

---

## Usage

### Quick start (run all steps)

```bash
# Set working directory
WORKDIR=/data/scapitis
mkdir -p ${WORKDIR}/{reads,assemblies,mlst_6gene,mlst_7gene,wgsnp,trees,results}

# Step 1: Download SRA reads
bash 01_download_sra.sh ${WORKDIR}

# Step 2: Assemble genomes
bash 02_assemble_genomes.sh ${WORKDIR}

# Step 3: QC check
bash 03_qc_check.sh ${WORKDIR}

# Step 4: Extract MLST alleles (6-gene and 7-gene)
# NOTE: strains with incomplete locus coverage are flagged
# has_missing_data=yes in the output and excluded from ST assignment ??# do not force-assign them an ST downstream.
python 04_extract_mlst.py \
    --assembly_dir ${WORKDIR}/assemblies \
    --scheme 6gene \
    --ref_dir reference_alleles \
    --output_dir ${WORKDIR}/mlst_6gene

python 04_extract_mlst.py \
    --assembly_dir ${WORKDIR}/assemblies \
    --scheme 7gene \
    --ref_dir reference_alleles \
    --output_dir ${WORKDIR}/mlst_7gene

# Step 5: Align and concatenate
bash 05_align_concat.sh ${WORKDIR} 6gene
bash 05_align_concat.sh ${WORKDIR} 7gene

# Step 6: wgSNP analysis
bash 06_wgsnp_snippy.sh ${WORKDIR}

# Step 7: Phylogenetic inference
bash 07_phylogeny_iqtree.sh ${WORKDIR}

# Step 8: Statistical comparison (9-metric table; reads ST assignment CSVs
# with keep_default_na=False to avoid the pandas "NA" pitfall)
python 08_statistics.py \
    --tree_6gene ${WORKDIR}/trees/tree_6gene_620.treefile \
    --tree_7gene ${WORKDIR}/trees/tree_7gene_620.treefile \
    --tree_wgsnp ${WORKDIR}/wgsnp/wgsnp_620.treefile \
    --st_6gene   ${WORKDIR}/mlst_6gene/st_assignments_final.csv \
    --st_7gene   ${WORKDIR}/mlst_7gene/st_assignments_final.csv \
    --cgmlst     wang2022_cgmlst_st.csv \
    --output     ${WORKDIR}/results/

# Step 9: Provisional ST catalogue
python 09_provisional_st.py --workdir ${WORKDIR}

# Steps 10-12: cgMLST concordance, MST figure, tanglegram
python 10_12_concordance_mst_tanglegram.py --workdir ${WORKDIR} --step concordance
python 10_12_concordance_mst_tanglegram.py --workdir ${WORKDIR} --step mst
python 10_12_concordance_mst_tanglegram.py --workdir ${WORKDIR} --step tanglegram
```

### Verifying your re-run

After Step 8, compare `${WORKDIR}/results/comparison_metrics.csv` against
`expected_output/comparison_metrics_expected.csv`. The two should match
within rounding tolerance (Mantel/RF/CID/ARI to 3 decimal places). If your
seven-gene ARI vs cgMLST cluster comes out far below ~0.92 (e.g. in the
0.7??.8 range), you have very likely hit the missing-data pitfall
described above ??check that strains in `has_missing_data=yes` were
correctly excluded before computing ARI.

---

## Input Data

- Public genomes: BioProject PRJNA493527 and additional GenBank accessions
  (see Supplementary Table S2 in the manuscript)
- YGH clinical isolates: GenBank accessions PZ469898?Z470071
  (six separate BankIt submissions, one per gene: femA, ftsZ, gap, rpoB,
  pyrH, tuf)
- Reference genome (wgSNP): *S. capitis* subsp. *capitis* DSM20326
  (GCF_040739495.1)

---

## Output Files

| File                                       | Description                                            |
| ------------------------------------------- | -------------------------------------------------------- |
| `concat_6gene_658.fasta`                    | Concatenated 6-gene alignment (3,032 bp, 658 strains)     |
| `concat_7gene_620.fasta`                     | Concatenated 7-gene alignment (2,780 bp, 620 strains)     |
| `st_assignments_final.csv`                   | Per-strain ST with `has_missing_data` / `missing_loci`    |
| `tree_6gene_620.treefile`                    | ML tree, 6-gene scheme (620 strains)                       |
| `tree_7gene_620.treefile`                    | ML tree, 7-gene scheme (620 strains)                       |
| `wgsnp_620.treefile`                         | ML tree, wgSNP (619 strains)                               |
| `results/comparison_metrics.csv`             | All 9 comparative metrics                                  |
| `results/pST_table_658_FINAL.csv`            | 91 provisional sequence types                               |
| `results/MST_6gene_pST.png`                  | Minimum spanning tree figure                                |

---

## Citation

If you use this code, please cite:

> Tsai M-S, Ling TZ. *Six-Gene versus Seven-Gene MLST for Staphylococcus
> capitis*: Performance Validation Against Whole-Genome SNP and
> Core-Genome MLST, with Detection of Putative L-Clone-Affiliated Isolates
> in Taiwan. *Journal of Clinical Microbiology* (submitted). DOI: [to be
> added on acceptance]

Please also cite the specific code version used via the GitHub Release
tag (e.g. `v1.0-jcm-submission`) rather than the `main` branch, since
`main` may be updated after publication.

---

## License

MIT License. See [LICENSE](LICENSE) for details.

---

## Contact

Moan-Shane Tsai, MD
Division of Infectious Diseases
Yuan's General Hospital
Kaohsiung, Taiwan
