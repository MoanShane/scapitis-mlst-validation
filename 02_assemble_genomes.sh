#!/usr/bin/env bash
# =============================================================
# Script 02: de novo genome assembly using Shovill/SPAdes
# =============================================================
# Usage: bash 02_assemble_genomes.sh <WORKDIR>
#
# Requirements: Shovill v1.1.0, SPAdes v3.15
# =============================================================

set -euo pipefail

WORKDIR="${1:-/data/scapitis}"
READS_DIR="${WORKDIR}/reads"
ASSEMBLIES_DIR="${WORKDIR}/assemblies"
THREADS=8
MIN_CONTIG=200
LOG="${WORKDIR}/logs/02_assembly.log"

mkdir -p "${ASSEMBLIES_DIR}" "${WORKDIR}/logs"

echo "=== Genome Assembly (Shovill v1.1.0) ===" | tee "${LOG}"
echo "Start: $(date)" | tee -a "${LOG}"

for acc_dir in "${READS_DIR}"/*/; do
    acc=$(basename "${acc_dir}")
    outdir="${ASSEMBLIES_DIR}/${acc}"

    # Skip if already assembled
    if [ -f "${outdir}/contigs.fa" ]; then
        echo "  [SKIP] ${acc}" | tee -a "${LOG}"
        continue
    fi

    # Find R1/R2 reads
    r1=$(ls "${acc_dir}"*_1.fastq.gz 2>/dev/null | head -1)
    r2=$(ls "${acc_dir}"*_2.fastq.gz 2>/dev/null | head -1)

    if [ -z "${r1}" ] || [ -z "${r2}" ]; then
        echo "  [WARN] No paired reads for ${acc}" | tee -a "${LOG}"
        continue
    fi

    echo "  [ASM] ${acc}" | tee -a "${LOG}"
    shovill \
        --R1 "${r1}" \
        --R2 "${r2}" \
        --outdir "${outdir}" \
        --assembler spades \
        --cpus "${THREADS}" \
        --minlen "${MIN_CONTIG}" \
        --force \
        2>>"${LOG}"

    # Copy contigs with strain name
    cp "${outdir}/contigs.fa" \
       "${ASSEMBLIES_DIR}/${acc}.fasta"

    echo "  [OK] ${acc}: $(grep -c '^>' ${ASSEMBLIES_DIR}/${acc}.fasta) contigs" \
        | tee -a "${LOG}"
done

echo "Done: $(date)" | tee -a "${LOG}"
echo "Total assemblies: $(ls ${ASSEMBLIES_DIR}/*.fasta | wc -l)" \
    | tee -a "${LOG}"
