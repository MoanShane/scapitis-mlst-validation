#!/usr/bin/env python3
# =============================================================
# Script 04: BLAST-based MLST allele extraction and ST assignment
# =============================================================
# Usage:
#   python 04_extract_mlst.py --assembly_dir DIR \
#       --scheme 6gene|7gene \
#       --ref_dir DIR \
#       --output_dir DIR
#
# Extracts MLST allele sequences from genome assemblies using
# BLASTn, then assigns sequence types (ST) per strain. Handles
# both 6-gene (Song 2019) and 7-gene (Wang 2025) MLST schemes.
#
# IMPORTANT — missing-data handling:
#   A strain for which BLAST fails to find a valid hit at one or
#   more loci has genuinely missing data at that locus. Such a
#   strain MUST be excluded from ST assignment (has_missing_data
#   = "yes"), not assigned a spurious allele derived from an
#   absent/empty sequence. Conflating "no hit" with "a real
#   zero-length allele" fragments true alleles into many spurious
#   ones and corrupts downstream concordance metrics (this caused
#   the seven-gene ARI vs cgMLST cluster to collapse from ~0.93 to
#   ~0.75 during development of this pipeline — see README).
#
# Requirements: BLAST+ v2.13.0, Biopython v1.79
# =============================================================

import argparse
import subprocess
import os
from pathlib import Path
from Bio import SeqIO
from collections import defaultdict
import csv

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
PERC_IDENTITY = 70
QCOV_HSP_PERC = 70
MAX_TARGETS = 5
OUTFMT = "6 sseqid pident length qlen sstart send sstrand bitscore"


def build_blast_db(assembly_path, db_path):
    """Build BLAST database from assembly."""
    cmd = [
        "makeblastdb",
        "-in", str(assembly_path),
        "-dbtype", "nucl",
        "-out", str(db_path),
        "-parse_seqids",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"makeblastdb failed: {result.stderr}")


def blast_gene(query_fasta, db_path):
    """Run BLASTn for one gene against one assembly."""
    cmd = [
        "blastn",
        "-query", str(query_fasta),
        "-db", str(db_path),
        "-outfmt", OUTFMT,
        "-perc_identity", str(PERC_IDENTITY),
        "-qcov_hsp_perc", str(QCOV_HSP_PERC),
        "-max_target_seqs", str(MAX_TARGETS),
        "-num_threads", "4",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.stdout.strip()


def extract_hit(db_path, sseqid, start, end, strand):
    """Extract the matched region from the assembly."""
    cmd = [
        "blastdbcmd",
        "-db", str(db_path),
        "-entry", sseqid,
        "-range", f"{start}-{end}",
        "-outfmt", "%s",
    ]
    if strand == "minus":
        cmd += ["-strand", "minus"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.stdout.strip()


def process_strain(strain, assembly_path, ref_dir, scheme,
                    output_dir, genes):
    """
    Extract all MLST alleles from one strain assembly.

    Returns a dict {gene: sequence} containing ONLY genes with a
    valid BLAST hit and successful extraction. Genes absent from
    this dict for a given strain represent genuinely missing data
    (no hit, or extraction failure) — they must never be
    backfilled with a placeholder sequence.
    """
    db_path = output_dir / "blast_dbs" / strain
    build_blast_db(assembly_path, db_path)

    results = {}
    for gene in genes:
        query = ref_dir / f"{gene}_ref.fasta"
        if not query.exists():
            print(f"  [WARN] Reference not found: {query}")
            continue

        blast_out = blast_gene(query, db_path)
        if not blast_out:
            print(f"  [MISS] {strain}: {gene} — no BLAST hit")
            continue

        # Take best hit (highest bitscore)
        best = None
        best_score = 0
        for line in blast_out.split("\n"):
            parts = line.split("\t")
            if len(parts) < 8:
                continue
            sseqid = parts[0]
            pident = float(parts[1])
            length = int(parts[2])
            qlen = int(parts[3])
            sstart = int(parts[4])
            send = int(parts[5])
            strand = parts[6]
            bitscore = float(parts[7])

            # Validate hit length (±50% of query) — a hit far
            # shorter or longer than expected is treated as no
            # valid hit, not as a short/long "allele"
            if not (qlen * 0.5 <= length <= qlen * 1.5):
                continue

            if bitscore > best_score:
                best_score = bitscore
                best = (sseqid, min(sstart, send),
                        max(sstart, send), strand)

        if best is None:
            print(f"  [MISS] {strain}: {gene} — no valid hit")
            continue

        sseqid, start, end, strand = best
        seq = extract_hit(db_path, sseqid, start, end, strand)
        if seq:
            results[gene] = seq
        else:
            print(f"  [ERR] {strain}: {gene} — extraction failed")

    return results


def assign_sequence_types(allele_seqs_per_gene, genes):
    """
    Assign ST numbers from per-gene allele sequences.

    allele_seqs_per_gene: {gene: {strain: sequence}}
        Only strains with a successfully extracted sequence for
        that gene appear in the inner dict — absence means
        missing data at that locus for that strain.

    Returns: list of dicts with keys
        strain, ST, is_novel, has_missing_data, missing_loci
    """
    all_strains = set()
    for gene in genes:
        all_strains |= set(allele_seqs_per_gene[gene].keys())

    # Identify strains missing at least one locus
    strain_missing_gene = defaultdict(list)
    for gene in genes:
        seqs_this_gene = allele_seqs_per_gene[gene]
        for strain in all_strains:
            if strain not in seqs_this_gene or not seqs_this_gene[strain]:
                strain_missing_gene[strain].append(gene)

    incomplete_strains = set(strain_missing_gene.keys())
    complete_strains = all_strains - incomplete_strains

    print(f"\n  Total strains          : {len(all_strains)}")
    print(f"  Complete (all loci)     : {len(complete_strains)}")
    print(f"  Incomplete (>=1 missing): {len(incomplete_strains)}")

    # Build allele catalogue from complete strains ONLY — a
    # missing/empty sequence must never define or match an allele
    gene_allele_map = {}
    for gene in genes:
        valid_seqs = [allele_seqs_per_gene[gene][s]
                      for s in complete_strains
                      if allele_seqs_per_gene[gene].get(s)]
        unique_seqs = list(dict.fromkeys(valid_seqs))
        gene_allele_map[gene] = {seq: i + 1
                                  for i, seq in enumerate(unique_seqs)}
        print(f"  {gene}: {len(unique_seqs)} unique alleles "
              f"(from {len(valid_seqs)} complete-strain sequences)")

    # Build per-strain profile for complete strains
    strain_profile = {}
    for strain in sorted(complete_strains):
        profile = []
        for gene in genes:
            seq = allele_seqs_per_gene[gene][strain]
            allele_num = gene_allele_map[gene].get(seq, -1)
            profile.append(allele_num)
        strain_profile[strain] = tuple(profile)

    # Assign ST per unique profile
    profile_to_st = {}
    next_st = 1
    profile_count = defaultdict(int)
    for strain in strain_profile:
        profile_count[strain_profile[strain]] += 1

    records = []
    for strain in sorted(all_strains):
        if strain in incomplete_strains:
            records.append({
                "strain": strain,
                "ST": "NA",
                "is_novel": "NA",
                "has_missing_data": "yes",
                "missing_loci": ";".join(strain_missing_gene[strain]),
            })
            continue

        profile = strain_profile[strain]
        if profile not in profile_to_st:
            profile_to_st[profile] = next_st
            next_st += 1
        st = profile_to_st[profile]
        is_novel = "yes" if profile_count[profile] == 1 else "no"
        records.append({
            "strain": strain,
            "ST": st,
            "is_novel": is_novel,
            "has_missing_data": "no",
            "missing_loci": "",
        })

    return records


def main():
    parser = argparse.ArgumentParser(
        description="BLAST-based MLST allele extraction and ST assignment")
    parser.add_argument("--assembly_dir", required=True,
                         help="Directory containing *.fasta assemblies")
    parser.add_argument("--scheme", required=True,
                         choices=["6gene", "7gene"],
                         help="MLST scheme")
    parser.add_argument("--ref_dir",
                         default="reference_alleles",
                         help="Directory with per-gene reference FASTAs")
    parser.add_argument("--output_dir", required=True,
                         help="Output directory")
    parser.add_argument("--strains", default=None,
                         help="Optional: text file with strain IDs (one per line)")
    args = parser.parse_args()

    asm_dir = Path(args.assembly_dir)
    ref_dir = Path(args.ref_dir) / args.scheme
    output_dir = Path(args.output_dir)
    genes = SCHEMES[args.scheme]["genes"]

    (output_dir / "blast_dbs").mkdir(parents=True, exist_ok=True)
    (output_dir / "alleles").mkdir(parents=True, exist_ok=True)

    print(f"=== MLST Allele Extraction ({args.scheme}) ===")
    print(f"Scheme: {SCHEMES[args.scheme]['desc']}")
    print(f"Genes: {', '.join(genes)}")

    # Get strain list
    if args.strains:
        with open(args.strains) as f:
            strains = [l.strip() for l in f if l.strip()]
        assemblies = {s: asm_dir / f"{s}.fasta" for s in strains}
    else:
        assemblies = {f.stem: f
                       for f in sorted(asm_dir.glob("*.fasta"))}

    print(f"Strains: {len(assemblies)}")

    # Per-gene output FASTA files
    gene_files = {
        gene: open(output_dir / "alleles" / f"{gene}_all.fasta", "w")
        for gene in genes
    }

    # In-memory store for ST assignment: {gene: {strain: sequence}}
    allele_seqs_per_gene = {gene: {} for gene in genes}

    success = defaultdict(int)
    total = len(assemblies)

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
                gene_files[gene].write(f">{strain}\n{results[gene]}\n")
                allele_seqs_per_gene[gene][strain] = results[gene]
                success[gene] += 1
            # NOTE: genes absent from `results` are intentionally
            # left absent from allele_seqs_per_gene[gene] — this
            # is how missing data is represented downstream.

        print(f" — {len(results)}/{len(genes)} genes")

    for f in gene_files.values():
        f.close()

    # ── Extraction summary ──────────────────────────────────────
    print(f"\n=== Extraction Summary ===")
    print(f"{'Gene':<10} {'Success':>8} {'Failed':>8} {'%':>6}")
    print("-" * 36)
    for gene in genes:
        n = success[gene]
        print(f"{gene:<10} {n:>8} {total-n:>8} "
              f"{n/total*100:>5.1f}%")

    # ── ST assignment (with correct missing-data exclusion) ────
    print(f"\n=== Assigning Sequence Types ===")
    records = assign_sequence_types(allele_seqs_per_gene, genes)

    out_csv = output_dir / "st_assignments_final.csv"
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["strain", "ST", "is_novel",
                           "has_missing_data", "missing_loci"])
        writer.writeheader()
        writer.writerows(records)

    n_assigned = sum(1 for r in records if r["has_missing_data"] == "no")
    n_excluded = sum(1 for r in records if r["has_missing_data"] == "yes")
    print(f"\n✓ Written: {out_csv}")
    print(f"  ST assigned        : {n_assigned}")
    print(f"  Excluded (missing) : {n_excluded}")
    if n_excluded > 0:
        print(f"  IMPORTANT: {n_excluded} strain(s) with missing locus data "
              f"are recorded with ST='NA' and must be excluded — not "
              f"force-assigned — in any downstream ARI/concordance "
              f"calculation. See README for the pandas 'NA' pitfall when "
              f"reading this file back in.")


if __name__ == "__main__":
    main()
