# scapitis-mlst-validation
Six-gene MLST validation pipeline for Staphylococcus capitis
# S. capitis 6-gene vs 7-gene MLST Validation Pipeline

## Overview

Source code for:

> Tsai M-S, Ling TZ. *Six-Gene versus Seven-Gene MLST for Staphylococcus capitis*: Performance Validation Against Whole-Genome SNP and Core-Genome MLST, with Detection of Putative L-Clone-Affiliated Isolates in Taiwan. *Journal of Clinical Microbiology* (submitted).

This pipeline validates the [Song et al. (2019)] six-gene MLST scheme for *Staphylococcus capitis* against the [Wang et al. (2025)](https://doi.org/10.1186/s12866-025-04339-z) seven-gene scheme, whole-genome SNP (wgSNP) analysis, and core-genome MLST, using 658 genomes from 620 public sequences + 29 YGH clinical isolates + 9 Song 2019 reference strains.

---

## Repository Structure

```
scapitis_mlst_code/
??? README.md
??? envs/
??  ??? environment.yml          # conda environment
??  ??? requirements_python.txt # Python packages
??? scripts/
??  ??? 01_download_sra.sh       # Download reads from SRA
??  ??? 02_assemble_genomes.sh   # de novo assembly (Shovill/SPAdes)
??  ??? 03_qc_check.sh           # Assembly QC
??  ??? 04_extract_mlst.py       # BLAST-based MLST allele extraction
??  ??? 05_align_concat.sh       # MAFFT alignment + trimAl + concat
??  ??? 06_wgsnp_snippy.sh       # Whole-genome SNP (Snippy)
??  ??? 07_phylogeny_iqtree.sh   # IQ-TREE2 ML trees
??  ??? 08_statistics.py         # Mantel/RF/CID/ARI/PIS/Bootstrap
??  ??? 09_provisional_st.py     # pST definition and catalogue
??  ??? 10_concordance.py        # cgMLST concordance analysis
??  ??? 11_mst_figure.py         # Minimum spanning tree figure
??  ??? 12_tanglegram.py         # Tanglegram (6-gene vs wgSNP)
??? docs/
    ??? pipeline_overview.md     # Detailed pipeline description
```

---

## Dependencies

### Software (versions used in this study)

| Software | Version | Use |
|----------|---------|-----|
| SRA Toolkit | 3.0.5 | Download SRA reads |
| Shovill | 1.1.0 | Genome assembly |
| SPAdes | 3.15 | Assembler (via Shovill) |
| BLAST+ | 2.13.0 | MLST allele extraction |
| MAFFT | 7.490 | Multiple sequence alignment |
| trimAl | 1.4.1 | Alignment trimming |
| Snippy | 4.6.0 | Whole-genome SNP analysis |
| IQ-TREE2 | 2.2.0 | Phylogenetic inference |
| Python | 3.9 | Statistical analysis |

### Python packages

```
biopython>=1.79
numpy>=1.21
scipy>=1.7
scikit-learn>=1.0
pandas>=1.3
matplotlib>=3.5
networkx>=2.8
```

---

## Installation

```bash
# Clone repository
git clone https://github.com/[your-username]/scapitis-mlst-validation.git
cd scapitis-mlst-validation

# Create conda environment
conda env create -f envs/environment.yml
conda activate scapitis_mlst

# Or install Python packages only
pip install -r envs/requirements_python.txt
```

---

## Usage

### Quick start (run all steps)

```bash
# Set working directory
WORKDIR=/data/scapitis
mkdir -p ${WORKDIR}/{reads,assemblies,mlst_6gene,mlst_7gene,wgsnp,trees,results}

# Step 1: Download SRA reads
bash scripts/01_download_sra.sh ${WORKDIR}

# Step 2: Assemble genomes
bash scripts/02_assemble_genomes.sh ${WORKDIR}

# Step 3: QC check
bash scripts/03_qc_check.sh ${WORKDIR}

# Step 4: Extract MLST alleles (6-gene and 7-gene)
python scripts/04_extract_mlst.py \
    --assembly_dir ${WORKDIR}/assemblies \
    --scheme 6gene \
    --output_dir ${WORKDIR}/mlst_6gene
python scripts/04_extract_mlst.py \
    --assembly_dir ${WORKDIR}/assemblies \
    --scheme 7gene \
    --output_dir ${WORKDIR}/mlst_7gene

# Step 5: Align and concatenate
bash scripts/05_align_concat.sh ${WORKDIR} 6gene
bash scripts/05_align_concat.sh ${WORKDIR} 7gene

# Step 6: wgSNP analysis
bash scripts/06_wgsnp_snippy.sh ${WORKDIR}

# Step 7: Phylogenetic inference
bash scripts/07_phylogeny_iqtree.sh ${WORKDIR}

# Step 8: Statistical comparison
python scripts/08_statistics.py \
    --tree_6gene ${WORKDIR}/trees/tree_6gene_620.treefile \
    --tree_7gene ${WORKDIR}/trees/tree_7gene_620.treefile \
    --tree_wgsnp ${WORKDIR}/wgsnp/wgsnp_620.treefile \
    --output ${WORKDIR}/results/

# Step 9-12: Additional analyses
python scripts/09_provisional_st.py --workdir ${WORKDIR}
python scripts/10_concordance.py    --workdir ${WORKDIR}
python scripts/11_mst_figure.py     --workdir ${WORKDIR}
python scripts/12_tanglegram.py     --workdir ${WORKDIR}
```

---

## Input Data

- Public genomes: BioProject PRJNA493527 and additional GenBank accessions (see Supplementary Table S2)
- YGH clinical isolates: GenBank accession numbers [XXXXXXXX?XXXXXXX]
- Reference genome (wgSNP): *S. capitis* subsp. *capitis* DSM20326 (GCF_040739495.1)

---

## Output Files

| File | Description |
|------|-------------|
| `concat_6gene_658.fasta` | Concatenated 6-gene alignment (3,032 bp, 658 strains) |
| `concat_7gene_620.fasta` | Concatenated 7-gene alignment (2,780 bp, 620 strains) |
| `tree_6gene_620.treefile` | ML tree, 6-gene scheme (620 strains) |
| `tree_7gene_620.treefile` | ML tree, 7-gene scheme (620 strains) |
| `wgsnp_620.treefile` | ML tree, wgSNP (619 strains) |
| `results/comparison_metrics.csv` | All 9 comparative metrics |
| `results/pST_table_658_FINAL.csv` | 91 provisional sequence types |
| `results/MST_6gene_pST.png` | Minimum spanning tree figure |

---

## Citation

If you use this code, please cite:

> Tsai M-S, Ling TZ. *Six-Gene versus Seven-Gene MLST for Staphylococcus capitis*...
> *Journal of Clinical Microbiology* (submitted). DOI: [to be added]

---

## License

MIT License. See [LICENSE](LICENSE) for details.

---

## Contact

Moan-Shane Tsai, MD  
Division of Infectious Diseases  
Yuan's General Hospital  
Kaohsiung, Taiwan  
