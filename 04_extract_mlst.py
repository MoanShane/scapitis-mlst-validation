#!/usr/bin/env python3
# =============================================================
# Script 04: BLAST-based MLST allele extraction
# =============================================================
# Usage:
#   python 04_extract_mlst.py --assembly_dir DIR
#                             --scheme 6gene|7gene
#                             --ref_dir DIR
#                             --output_dir DIR
#
# Extracts MLST allele sequences from genome assemblies using
# BLASTn. Handles both 6-gene (Song 2019) and 7-gene (Wang
# 2025) MLST schemes.
#
# Requirements: BLAST+ v2.13.0, Biopython v1.79
# =============================================================

import argparse
import subprocess
import os
from pathlib import Path
from Bio import SeqIO
from collections import defaultdict

# ── MLST scheme definitions ──────────────────────────────────
SCHEMES = {
    "6gene": {
        "genes": ["femA", "ftsZ", "gap", "pyrH", "rpoB", "tuf"],
        "desc": "Song et al. (2019) 6-gene MLST",
    },
    "7gene": {
        "genes": ["atpB_2", "carB", "clpP", "hisS",
                  "mntC", "phoA", "rluB"],
        "desc": "Wang et al. (2025) 7-gene MLST",
    }
}

# BLAST parameters
PERC_IDENTITY  = 70
QCOV_HSP_PERC  = 70
MAX_TARGETS    = 5
OUTFMT = "6 sseqid pident length qlen sstart send sstrand bitscore"


def build_blast_db(assembly_path, db_path):
    """Build BLAST database from assembly."""
    cmd = [
        "makeblastdb",
        "-in",    str(assembly_path),
        "-dbtype","nucl",
        "-out",   str(db_path),
        "-parse_seqids"
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"makeblastdb failed: {result.stderr}")


def blast_gene(query_fasta, db_path):
    """Run BLASTn for one gene against one assembly."""
    cmd = [
        "blastn",
        "-query",          str(query_fasta),
        "-db",             str(db_path),
        "-outfmt",         OUTFMT,
        "-perc_identity",  str(PERC_IDENTITY),
        "-qcov_hsp_perc",  str(QCOV_HSP_PERC),
        "-max_target_seqs",str(MAX_TARGETS),
        "-num_threads",    "4",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.stdout.strip()


def extract_hit(db_path, sseqid, start, end, strand):
    """Extract the matched region from the assembly."""
    cmd = [
        "blastdbcmd",
        "-db",     str(db_path),
        "-entry",  sseqid,
        "-range",  f"{start}-{end}",
        "-outfmt", "%s",
    ]
    if strand == "minus":
        cmd += ["-strand", "minus"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.stdout.strip()


def process_strain(strain, assembly_path, ref_dir, scheme,
                   output_dir, genes):
    """Extract all MLST alleles from one strain assembly."""
    db_path = output_dir / "blast_dbs" / strain

    # Build BLAST DB
    build_blast_db(assembly_path, db_path)

    results = {}
    for gene in genes:
        query = ref_dir / f"{gene}_ref.fasta"
        if not query.exists():
            print(f"    [WARN] Reference not found: {query}")
            continue

        blast_out = blast_gene(query, db_path)
        if not blast_out:
            print(f"    [MISS] {strain}: {gene} — no BLAST hit")
            continue

        # Take best hit (highest bitscore)
        best = None
        best_score = 0
        for line in blast_out.split("\n"):
            parts = line.split("\t")
            if len(parts) < 8:
                continue
            sseqid   = parts[0]
            pident   = float(parts[1])
            length   = int(parts[2])
            qlen     = int(parts[3])
            sstart   = int(parts[4])
            send     = int(parts[5])
            strand   = parts[6]
            bitscore = float(parts[7])

            # Validate hit length (±50% of query)
            if not (qlen * 0.5 <= length <= qlen * 1.5):
                continue

            if bitscore > best_score:
                best_score = bitscore
                best = (sseqid, min(sstart, send),
                        max(sstart, send), strand)

        if best is None:
            print(f"    [MISS] {strain}: {gene} — no valid hit")
            continue

        sseqid, start, end, strand = best
        seq = extract_hit(db_path, sseqid, start, end, strand)
        if seq:
            results[gene] = seq
        else:
            print(f"    [ERR] {strain}: {gene} — extraction failed")

    return results


def main():
    parser = argparse.ArgumentParser(
        description="BLAST-based MLST allele extraction")
    parser.add_argument("--assembly_dir", required=True,
                        help="Directory containing *.fasta assemblies")
    parser.add_argument("--scheme", required=True,
                        choices=["6gene","7gene"],
                        help="MLST scheme")
    parser.add_argument("--ref_dir",
                        default="reference_alleles",
                        help="Directory with per-gene reference FASTAs")
    parser.add_argument("--output_dir", required=True,
                        help="Output directory")
    parser.add_argument("--strains", default=None,
                        help="Optional: text file with strain IDs (one per line)")
    args = parser.parse_args()

    asm_dir    = Path(args.assembly_dir)
    ref_dir    = Path(args.ref_dir) / args.scheme
    output_dir = Path(args.output_dir)
    genes      = SCHEMES[args.scheme]["genes"]

    (output_dir / "blast_dbs").mkdir(parents=True, exist_ok=True)
    (output_dir / "alleles").mkdir(parents=True, exist_ok=True)

    print(f"=== MLST Allele Extraction ({args.scheme}) ===")
    print(f"Scheme: {SCHEMES[args.scheme]['desc']}")
    print(f"Genes:  {', '.join(genes)}")

    # Get strain list
    if args.strains:
        with open(args.strains) as f:
            strains = [l.strip() for l in f if l.strip()]
        assemblies = {s: asm_dir / f"{s}.fasta" for s in strains}
    else:
        assemblies = {f.stem: f
                      for f in sorted(asm_dir.glob("*.fasta"))}

    print(f"Strains: {len(assemblies)}")

    # Per-gene output files
    gene_files = {
        gene: open(output_dir / "alleles" / f"{gene}_all.fasta", "w")
        for gene in genes
    }

    # Stats
    success = defaultdict(int)
    total   = len(assemblies)

    for i, (strain, asm_path) in enumerate(assemblies.items(), 1):
        if not asm_path.exists():
            print(f"  [{i}/{total}] {strain}: assembly not found")
            continue
        print(f"  [{i}/{total}] {strain}", end="", flush=True)

        results = process_strain(
            strain, asm_path, ref_dir,
            args.scheme, output_dir, genes)

        for gene in genes:
            if gene in results:
                gene_files[gene].write(
                    f">{strain}\n{results[gene]}\n")
                success[gene] += 1
        print(f" — {len(results)}/{len(genes)} genes")

    for f in gene_files.values():
        f.close()

    # Summary
    print(f"\n=== Extraction Summary ===")
    print(f"{'Gene':<10} {'Success':>8} {'Failed':>8} {'%':>6}")
    print("─" * 36)
    for gene in genes:
        n = success[gene]
        print(f"{gene:<10} {n:>8} {total-n:>8} "
              f"{n/total*100:>5.1f}%")


if __name__ == "__main__":
    main()
