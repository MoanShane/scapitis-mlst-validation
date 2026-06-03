#!/usr/bin/env bash
# =============================================================
# Script 07: Phylogenetic inference using IQ-TREE2
# =============================================================
# Usage: bash 07_phylogeny_iqtree.sh <WORKDIR>
#
# Runs IQ-TREE2 for:
#   (a) 6-gene partitioned MLST tree (620 and 658 strains)
#   (b) 7-gene partitioned MLST tree (620 strains)
#   (c) wgSNP tree (619 strains)
#
# All trees use:
#   - ModelFinder (-m MFP) for model selection
#   - 1,000 ultrafast bootstrap replicates (-B 1000)
#   - DSM20326 as outgroup
#
# Requirements: IQ-TREE2 v2.2.0
# =============================================================

set -euo pipefail

WORKDIR="${1:-/data/scapitis}"
THREADS="${2:-AUTO}"
LOG="${WORKDIR}/logs/07_phylogeny.log"

mkdir -p "${WORKDIR}/trees" "${WORKDIR}/logs"

echo "=== Phylogenetic inference (IQ-TREE2 v2.2.0) ===" | tee "${LOG}"
echo "Start: $(date)" | tee -a "${LOG}"

# ── Function to run IQ-TREE2 ─────────────────────────────────
run_iqtree() {
    local NAME="$1"
    local FASTA="$2"
    local PARTITION="$3"   # empty = no partition
    local OUTGROUP="$4"
    local PREFIX="${WORKDIR}/trees/${NAME}"

    echo "" | tee -a "${LOG}"
    echo "  Running: ${NAME}" | tee -a "${LOG}"

    if [ -f "${PREFIX}.treefile" ]; then
        echo "  [SKIP] ${NAME} — treefile exists" | tee -a "${LOG}"
        return
    fi

    local CMD=(
        iqtree2
        -s "${FASTA}"
        -m MFP
        -B 1000
        -T "${THREADS}"
        --prefix "${PREFIX}"
        -o "${OUTGROUP}"
        --safe
    )

    if [ -n "${PARTITION}" ]; then
        CMD+=(-p "${PARTITION}")
    else
        CMD+=(-m GTR+G)
    fi

    "${CMD[@]}" >> "${LOG}" 2>&1

    echo "  [OK] ${NAME}: $(cat ${PREFIX}.log | \
         grep 'Log-likelihood' | tail -1)" | tee -a "${LOG}"
}

# ── (a) 6-gene MLST trees ────────────────────────────────────
# 620-strain public dataset
run_iqtree \
    "tree_6gene_620" \
    "${WORKDIR}/concat_6gene/concat_6gene_620.fasta" \
    "${WORKDIR}/concat_6gene/partition_6gene_620.txt" \
    "DSM20326"

# 658-strain full dataset (620 public + 29 YGH + 9 Song refs)
run_iqtree \
    "tree_6gene_658" \
    "${WORKDIR}/concat_6gene/concat_6gene_658.fasta" \
    "${WORKDIR}/concat_6gene/partition_6gene_658.txt" \
    "DSM20326"

# ── (b) 7-gene MLST tree (620 strains) ───────────────────────
run_iqtree \
    "tree_7gene_620" \
    "${WORKDIR}/concat_7gene/concat_7gene_620.fasta" \
    "${WORKDIR}/concat_7gene/partition_7gene_620.txt" \
    "DSM20326"

# ── (c) wgSNP tree (619 strains, excl. DSM20326 reference) ──
run_iqtree \
    "wgsnp_620" \
    "${WORKDIR}/wgsnp/core/core.aln" \
    "" \
    "DSM20326"

echo "" | tee -a "${LOG}"
echo "=== All trees complete ===" | tee -a "${LOG}"
echo "Done: $(date)" | tee -a "${LOG}"
echo "" | tee -a "${LOG}"
echo "Output files:" | tee -a "${LOG}"
ls "${WORKDIR}/trees/"*.treefile | tee -a "${LOG}"
