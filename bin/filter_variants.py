#!/usr/bin/env python

# Copyright 2024 - GHGA
# Author: Kuebra Narci - @kubranarci
#!/usr/bin/env python3
#!/usr/bin/env python3
"""
filter_variants.py

Filters a CSV of variants by:
1. Optional overlap with one or two BED files (1-based inclusive)
2. Presence in a VCF file (matching CHROM, POS, REF, ALT)

Usage:
    python filter_variants.py --csv variants_to_extract.csv \
                              --vcf input_vcf.vcf.gz \
                              --bed regions.bed \
                              --targets targets.bed \
                              --out filtered.csv
"""

import argparse
import pandas as pd
import pysam

def load_bed(bed_file):
    """Load BED regions into a dictionary of intervals per chromosome."""
    bed_dict = {}
    if not bed_file:
        return bed_dict
    with open(bed_file) as f:
        for line in f:
            if line.strip() == '':
                continue
            chrom, start, end = line.strip().split()[:3]
            start = int(start)
            end = int(end)
            if chrom not in bed_dict:
                bed_dict[chrom] = []
            bed_dict[chrom].append((start, end))
    # Sort intervals for each chromosome
    for chrom in bed_dict:
        bed_dict[chrom].sort()
    return bed_dict

def merge_beds(*bed_dicts):
    """Merge multiple BED dictionaries into one."""
    merged = {}
    for bd in bed_dicts:
        for chrom, intervals in bd.items():
            if chrom not in merged:
                merged[chrom] = []
            merged[chrom].extend(intervals)
    # Sort merged intervals
    for chrom in merged:
        merged[chrom].sort()
    return merged

def pos_in_bed(chrom, pos, bed_dict):
    """Check if a 1-based position is inside any interval of BED (1-based inclusive)."""
    if not bed_dict:
        return True  # No BED provided → accept all positions
    if chrom not in bed_dict:
        return False
    for start, end in bed_dict[chrom]:
        if start <= pos <= end:
            return True
    return False

def load_vcf_variants(vcf_file):
    """Load VCF variants into a set of (CHROM, POS, REF, ALT)."""
    vcf_set = set()
    vcf = pysam.VariantFile(vcf_file)
    for rec in vcf.fetch():
        for alt in rec.alts:
            vcf_set.add((rec.chrom, rec.pos, rec.ref, alt))
    return vcf_set

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True, help="Input CSV file with variants")
    parser.add_argument("--vcf", required=True, help="Input bgzipped + indexed VCF")
    parser.add_argument("--bed", default=None, help="Optional BED file (1-based inclusive)")
    parser.add_argument("--targets", default=None, help="Optional second BED file (1-based inclusive)")
    parser.add_argument("--out", required=True, help="Output filtered CSV file")
    args = parser.parse_args()

    # Load CSV
    df = pd.read_csv(args.csv)
    df['POS'] = df['POS'].astype(int)

    # Load BED(s)
    bed1 = load_bed(args.bed)
    bed2 = load_bed(args.targets)
    bed_dict = merge_beds(bed1, bed2)

    # Load VCF variants
    vcf_variants = load_vcf_variants(args.vcf)

    # Filter by BED (if any)
    df['in_bed'] = df.apply(lambda row: pos_in_bed(row['CHROM'], row['POS'], bed_dict), axis=1)
    df_bed_filtered = df[df['in_bed']]

    # Filter by VCF
    df_final = df_bed_filtered[df_bed_filtered.apply(
        lambda row: (row['CHROM'], row['POS'], row['REF'], row['ALT']) in vcf_variants, axis=1
    )]

    # Drop helper column and save
    df_final.drop(columns=['in_bed'], inplace=True)
    df_final.to_csv(args.out, index=False)
    print(f"Filtered CSV saved to {args.out}, {len(df_final)} variants retained.")

if __name__ == "__main__":
    main()
