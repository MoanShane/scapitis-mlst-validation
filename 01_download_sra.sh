#!/usr/bin/env bash
# =============================================================
# Script 01: Download S. capitis genomes from SRA/NCBI
# =============================================================
# Usage: bash 01_download_sra.sh <WORKDIR>
#
# Downloads all 620 public S. capitis genomes used in:
# Tsai & Ling (submitted) - JCM
#
# Requirements: SRA Toolkit v3.0.5, parallel (optional)
# =============================================================

set -euo pipefail

WORKDIR="${1:-/data/scapitis}"
READS_DIR="${WORKDIR}/reads"
ASSEMBLIES_DIR="${WORKDIR}/assemblies"
ACCESSION_LIST="${WORKDIR}/accession_list.txt"
THREADS=8
LOG="${WORKDIR}/logs/01_download.log"

mkdir -p "${READS_DIR}" "${ASSEMBLIES_DIR}" \
         "${WORKDIR}/logs"

echo "=== Download S. capitis genomes ===" | tee "${LOG}"
echo "Start: $(date)" | tee -a "${LOG}"

# ── SRR/ERR/DRR accessions: download as raw reads ──────────
while IFS= read -r acc; do
    [[ "$acc" =~ ^(SRR|ERR|DRR) ]] || continue
    outdir="${READS_DIR}/${acc}"
    if [ -d "${outdir}" ] && ls "${outdir}"/*.fastq.gz &>/dev/null; then
        echo "  [SKIP] ${acc} already downloaded"
        continue
    fi
    echo "  [DL] ${acc}" | tee -a "${LOG}"
    prefetch "${acc}" --output-directory "${READS_DIR}" \
             --max-size 10GB 2>>"${LOG}"
    fasterq-dump "${READS_DIR}/${acc}" \
        --outdir "${outdir}" \
        --threads "${THREADS}" \
        --split-files 2>>"${LOG}"
    # Compress
    gzip -f "${outdir}"/*.fastq 2>>/dev/null || true
done < "${ACCESSION_LIST}"

# ── GCA accessions: download as assembled genomes ──────────
while IFS= read -r acc; do
    [[ "$acc" =~ ^GCA ]] || continue
    outdir="${ASSEMBLIES_DIR}/${acc}"
    if [ -f "${outdir}/${acc}.fasta" ]; then
        echo "  [SKIP] ${acc} already downloaded"
        continue
    fi
    mkdir -p "${outdir}"
    echo "  [DL] ${acc}" | tee -a "${LOG}"
    datasets download genome accession "${acc}" \
        --filename "${outdir}/${acc}.zip" \
        --no-progressbar 2>>"${LOG}"
    unzip -q "${outdir}/${acc}.zip" -d "${outdir}/tmp"
    find "${outdir}/tmp" -name "*.fna" -exec \
        cp {} "${outdir}/${acc}.fasta" \;
    rm -rf "${outdir}/tmp" "${outdir}/${acc}.zip"
done < "${ACCESSION_LIST}"

echo "Done: $(date)" | tee -a "${LOG}"
echo "=== Download complete ===" | tee -a "${LOG}"
