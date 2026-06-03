#!/usr/bin/env bash
# =============================================================
# Script 03: Assembly quality control
# =============================================================
# Usage: bash 03_qc_check.sh <WORKDIR>
#
# Computes N50, total length, contig count for each assembly.
# Filters out assemblies outside expected S. capitis range.
# Expected: total 2.3–2.7 Mb, N50 > 10,000 bp
# =============================================================

set -euo pipefail

WORKDIR="${1:-/data/scapitis}"
ASSEMBLIES_DIR="${WORKDIR}/assemblies"
QC_OUT="${WORKDIR}/results/assembly_qc.tsv"
MIN_SIZE=2300000   # 2.3 Mb
MAX_SIZE=2700000   # 2.7 Mb
MIN_N50=10000

mkdir -p "${WORKDIR}/results"

echo -e "Strain\tTotal_bp\tContigs\tN50\tGC_pct\tStatus" \
    > "${QC_OUT}"

python3 - << 'PYEOF'
import os, sys
from pathlib import Path

WORKDIR = sys.argv[1] if len(sys.argv) > 1 else "/data/scapitis"
ASM_DIR = Path(WORKDIR) / "assemblies"
QC_OUT  = Path(WORKDIR) / "results" / "assembly_qc.tsv"
MIN_SIZE, MAX_SIZE, MIN_N50 = 2300000, 2700000, 10000

rows = []
for fasta in sorted(ASM_DIR.glob("*.fasta")):
    strain = fasta.stem
    seqs, lengths = [], []
    total_gc = 0
    with open(fasta) as f:
        seq = ""
        for line in f:
            if line.startswith(">"):
                if seq:
                    lengths.append(len(seq))
                    total_gc += seq.count("G") + seq.count("C")
                seq = ""
            else:
                seq += line.strip().upper()
        if seq:
            lengths.append(len(seq))
            total_gc += seq.count("G") + seq.count("C")

    total_bp = sum(lengths)
    n_contigs = len(lengths)
    gc_pct = total_gc / total_bp * 100 if total_bp else 0

    # N50
    lengths_sorted = sorted(lengths, reverse=True)
    cumsum = 0
    n50 = 0
    for l in lengths_sorted:
        cumsum += l
        if cumsum >= total_bp / 2:
            n50 = l
            break

    # QC status
    ok = (MIN_SIZE <= total_bp <= MAX_SIZE and n50 >= MIN_N50)
    status = "PASS" if ok else "FAIL"

    rows.append((strain, total_bp, n_contigs, n50,
                 f"{gc_pct:.1f}", status))
    print(f"  {strain}: {total_bp:,} bp, N50={n50:,}, "
          f"GC={gc_pct:.1f}%, {status}")

with open(QC_OUT, "a") as f:
    for row in rows:
        f.write("\t".join(str(x) for x in row) + "\n")

pass_n = sum(1 for r in rows if r[-1] == "PASS")
fail_n = sum(1 for r in rows if r[-1] == "FAIL")
print(f"\nQC summary: {pass_n} PASS, {fail_n} FAIL")
print(f"QC report: {QC_OUT}")
PYEOF "${WORKDIR}"

echo "=== Assembly QC complete ==="
