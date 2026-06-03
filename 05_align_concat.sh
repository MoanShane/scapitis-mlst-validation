#!/usr/bin/env bash
# =============================================================
# Script 05: Multiple sequence alignment, trimming, concatenation
# =============================================================
# Usage: bash 05_align_concat.sh <WORKDIR> <SCHEME: 6gene|7gene>
#
# Steps:
#   1. MAFFT alignment (per gene)
#   2. trimAl trimming (-automated1)
#   3. Concatenation in fixed gene order
#   4. Partition file generation for IQ-TREE2
#
# Requirements: MAFFT v7.490, trimAl v1.4.1, Python 3.9
# =============================================================

set -euo pipefail

WORKDIR="${1:-/data/scapitis}"
SCHEME="${2:-6gene}"
THREADS=8

case "${SCHEME}" in
    6gene)
        GENES="femA ftsZ gap pyrH rpoB tuf"
        ALN_DIR="${WORKDIR}/aligned_6gene"
        TRM_DIR="${WORKDIR}/trimmed_6gene"
        CAT_DIR="${WORKDIR}/concat_6gene"
        ;;
    7gene)
        GENES="atpB_2 carB clpP hisS mntC phoA rluB"
        ALN_DIR="${WORKDIR}/aligned_7gene"
        TRM_DIR="${WORKDIR}/trimmed_7gene"
        CAT_DIR="${WORKDIR}/concat_7gene"
        ;;
    *)
        echo "Error: scheme must be 6gene or 7gene" >&2
        exit 1
        ;;
esac

RAW_DIR="${WORKDIR}/extracted_${SCHEME}/alleles"
mkdir -p "${ALN_DIR}" "${TRM_DIR}" "${CAT_DIR}"

echo "=== Alignment and Concatenation (${SCHEME}) ==="
echo "Genes: ${GENES}"

# ── Step 1: MAFFT alignment ──────────────────────────────────
echo ""
echo "Step 1: MAFFT alignment"
for gene in ${GENES}; do
    INPUT="${RAW_DIR}/${gene}_all.fasta"
    OUTPUT="${ALN_DIR}/${gene}_aligned.fasta"

    if [ ! -f "${INPUT}" ]; then
        echo "  [WARN] Missing: ${INPUT}"
        continue
    fi

    echo -n "  ${gene}... "
    mafft \
        --auto \
        --thread ${THREADS} \
        --quiet \
        "${INPUT}" > "${OUTPUT}"
    echo "$(grep -c '^>' ${OUTPUT}) sequences aligned"
done

# ── Step 2: trimAl trimming ──────────────────────────────────
echo ""
echo "Step 2: trimAl trimming"
for gene in ${GENES}; do
    INPUT="${ALN_DIR}/${gene}_aligned.fasta"
    OUTPUT="${TRM_DIR}/${gene}_trimmed.fasta"

    [ -f "${INPUT}" ] || continue

    echo -n "  ${gene}... "
    trimal \
        -in "${INPUT}" \
        -out "${OUTPUT}" \
        -automated1 \
        -fasta
    LEN=$(awk 'NR==2{print length($0); exit}' "${OUTPUT}")
    echo "${LEN} bp after trimming"
done

# ── Step 3 & 4: Concatenation + partition file ───────────────
echo ""
echo "Step 3-4: Concatenation + partition file"

python3 - "${TRM_DIR}" "${CAT_DIR}" "${SCHEME}" \
         $(echo ${GENES}) << 'PYEOF'
import sys
from pathlib import Path
from collections import defaultdict

trim_dir = Path(sys.argv[1])
cat_dir  = Path(sys.argv[2])
scheme   = sys.argv[3]
genes    = sys.argv[4:]

# Read all trimmed sequences
strain_seqs = defaultdict(dict)
gene_len    = {}

for gene in genes:
    fasta = trim_dir / f"{gene}_trimmed.fasta"
    if not fasta.exists():
        print(f"  [WARN] Missing trimmed: {fasta}")
        continue
    seq = ""
    current = None
    with open(fasta) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current and seq:
                    strain_seqs[current][gene] = seq
                current = line[1:].strip()
                seq = ""
            else:
                seq += line.upper()
        if current and seq:
            strain_seqs[current][gene] = seq

    # All sequences should be same length after trimAl
    lens = [len(strain_seqs[s][gene])
            for s in strain_seqs if gene in strain_seqs[s]]
    if lens:
        gene_len[gene] = lens[0]
        print(f"  {gene}: {gene_len[gene]} bp, "
              f"{len(lens)} sequences")

# Get strains with ALL genes
complete = [s for s in strain_seqs
            if all(g in strain_seqs[s] for g in genes
                   if g in gene_len)]
print(f"\n  Strains with all {len(genes)} genes: {len(complete)}")

# Write concatenated FASTA
n_strains = len(complete)
cat_file  = cat_dir / f"concat_{scheme}_{n_strains}.fasta"
with open(cat_file, "w") as f:
    for strain in sorted(complete):
        concat = "".join(
            strain_seqs[strain].get(g, "N" * gene_len.get(g, 0))
            for g in genes)
        f.write(f">{strain}\n{concat}\n")
print(f"\n  Concatenated FASTA: {cat_file}")
print(f"  Total length: {sum(gene_len.values())} bp")

# Write IQ-TREE2 partition file
part_file = cat_dir / f"partition_{scheme}_{n_strains}.txt"
pos = 1
with open(part_file, "w") as f:
    for gene in genes:
        if gene not in gene_len:
            continue
        end = pos + gene_len[gene] - 1
        f.write(f"DNA, {gene} = {pos}-{end}\n")
        pos = end + 1
print(f"  Partition file: {part_file}")
PYEOF

echo ""
echo "=== Alignment and concatenation complete ==="
