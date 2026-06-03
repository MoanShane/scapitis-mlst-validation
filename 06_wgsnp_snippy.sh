#!/usr/bin/env bash
# =============================================================
# Script 06: Whole-genome SNP analysis using Snippy
# =============================================================
# Usage: bash 06_wgsnp_snippy.sh <WORKDIR>
#
# Reference genome: S. capitis subsp. capitis DSM20326
# (GCF_040739495.1) — used as both mapping reference and
# phylogenetic outgroup.
#
# Note on reference choice: DSM20326 was chosen over CR01
# (GCA_000499705.1) because (i) CR01 was suppressed from NCBI,
# and (ii) DSM20326 is the officially designated type strain
# providing a phylogenetically appropriate outgroup that avoids
# NRCS-A-specific mapping bias.
#
# Requirements: Snippy v4.6.0
# =============================================================

set -euo pipefail

WORKDIR="${1:-/data/scapitis}"
ASSEMBLIES_DIR="${WORKDIR}/assemblies"
WGSNP_DIR="${WORKDIR}/wgsnp"
REF="${WORKDIR}/reference/DSM20326.fa"  # GCF_040739495.1
THREADS=8
LOG="${WORKDIR}/logs/06_wgsnp.log"

mkdir -p "${WGSNP_DIR}" "${WORKDIR}/logs"

echo "=== Whole-genome SNP analysis (Snippy v4.6.0) ===" | tee "${LOG}"
echo "Reference: DSM20326 (GCF_040739495.1)" | tee -a "${LOG}"
echo "Start: $(date)" | tee -a "${LOG}"

# ── Per-strain SNP calling (assembly-based) ──────────────────
TOTAL=$(ls "${ASSEMBLIES_DIR}"/*.fasta | wc -l)
COUNT=0

for assembly in "${ASSEMBLIES_DIR}"/*.fasta; do
    strain=$(basename "${assembly}" .fasta)
    outdir="${WGSNP_DIR}/${strain}"
    COUNT=$((COUNT + 1))

    if [ -f "${outdir}/snps.vcf" ]; then
        echo "  [SKIP ${COUNT}/${TOTAL}] ${strain}" | tee -a "${LOG}"
        continue
    fi

    echo "  [${COUNT}/${TOTAL}] ${strain}" | tee -a "${LOG}"
    snippy \
        --ctgs "${assembly}" \
        --ref  "${REF}" \
        --outdir "${outdir}" \
        --cpus "${THREADS}" \
        --force \
        --quiet \
        2>>"${LOG}"
done

# ── Core SNP alignment ───────────────────────────────────────
echo "" | tee -a "${LOG}"
echo "Generating core SNP alignment..." | tee -a "${LOG}"

OUTDIR_LIST=$(ls -d "${WGSNP_DIR}"/*/  | \
              grep -v "^${WGSNP_DIR}/core" | \
              tr '\n' ' ')

snippy-core \
    --ref   "${REF}" \
    --prefix "${WGSNP_DIR}/core/core" \
    ${OUTDIR_LIST} \
    2>>"${LOG}"

# Report
NSNPS=$(grep -v '^#' "${WGSNP_DIR}/core/core.tab" | wc -l)
NSTRAINS=$(grep -c '^>' "${WGSNP_DIR}/core/core.aln" || true)
echo "" | tee -a "${LOG}"
echo "Core SNP alignment:" | tee -a "${LOG}"
echo "  Strains: ${NSTRAINS}" | tee -a "${LOG}"
echo "  SNP sites: ${NSNPS}" | tee -a "${LOG}"
echo "  File: ${WGSNP_DIR}/core/core.aln" | tee -a "${LOG}"
echo "Done: $(date)" | tee -a "${LOG}"
