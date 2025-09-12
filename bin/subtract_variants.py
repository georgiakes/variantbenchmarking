#!/usr/bin/env python

# Copyright 2024 - GHGA
# Author: Kuebra Narci - @kubranarci

import argparse
import sys
import gzip

def parse_vcf(file_path):
    """
    Parses a VCF file and returns a set of unique variants.
    Variants are identified by a tuple of (chromosome, position, reference_allele, alternate_allele).
    It returns a set of unique variant tuples.
    """
    variants = set()
    try:
        if file_path.endswith('.gz'):
            vcf_file = gzip.open(file_path, 'rt')
        else:
            vcf_file = open(file_path, 'r')

        for line in vcf_file:
            # Skip header and comment lines
            if line.startswith('#'):
                continue

            # Split the line into columns
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                # Extract the key fields for variant identification
                # CHROM, POS, REF, ALT
                chrom = parts[0]
                pos = parts[1]
                ref = parts[3]
                alt = parts[4]

                # Create a unique identifier for the variant
                variant_id = (chrom, pos, ref, alt)
                variants.add(variant_id)

        vcf_file.close()

    except FileNotFoundError:
        print(f"Error: File not found at {file_path}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred while parsing {file_path}: {e}", file=sys.stderr)
        sys.exit(1)

    return variants

def parse_bed_file(file_path):
    """
    Parses a BED file and returns a dictionary of regions.
    The dictionary keys are chromosome names, and values are lists of (start, end) tuples.
    """
    regions = {}
    try:
        with open(file_path, 'r') as bed_file:
            for line in bed_file:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    if chrom not in regions:
                        regions[chrom] = []
                    regions[chrom].append((start, end))
    except FileNotFoundError:
        print(f"Error: BED file not found at {file_path}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred while parsing {file_path}: {e}", file=sys.stderr)
        sys.exit(1)

    return regions

def is_in_region(chrom, pos, bed_regions):
    """
    Checks if a given position is within any of the regions for a chromosome.
    The BED format is 0-indexed, but VCF is 1-indexed. We'll handle this by
    checking if the position is within [start+1, end].

    Args:
        chrom (str): The chromosome of the variant.
        pos (int): The 1-based position of the variant.
        bed_regions (dict): The dictionary of BED regions.

    Returns:
        bool: True if the position is in a region, False otherwise.
    """
    if chrom in bed_regions:
        for start, end in bed_regions[chrom]:
            # VCF position is 1-based, BED start is 0-based, BED end is exclusive.
            # A common interpretation is to check if pos >= start+1 and pos <= end.
            if pos >= start + 1 and pos <= end:
                return True
    return False

def get_sample_name_from_header(header_line):
    """
    Extracts the sample name from a VCF header line.

    Args:
        header_line (str): A VCF header line starting with '#CHROM'

    Returns:
        str: The sample name or None if not found
    """
    if not header_line.startswith('#CHROM'):
        return None

    parts = header_line.strip().split('\t')
    if len(parts) >= 10:  # VCF format has at least 9 fixed columns + sample columns
        return parts[9]  # First sample column
    return None

def update_sample_name_in_header(header_line, new_sample_name):
    """
    Updates the sample name in a VCF header line.

    Args:
        header_line (str): A VCF header line starting with '#CHROM'
        new_sample_name (str): The new sample name to use

    Returns:
        str: The updated header line
    """
    if not header_line.startswith('#CHROM'):
        return header_line

    parts = header_line.strip().split('\t')
    if len(parts) >= 10:
        parts[9] = new_sample_name  # Update the first sample column
        return '\t'.join(parts) + '\n'

    return header_line

def subtract_vcf_files(primary_vcf, to_subtract_vcf, output_vcf, bed_file=None, zip_output=False, sample_name=None):
    """
    Subtracts variants from the primary VCF that exist in the second VCF.
    It writes the filtered variants to a new output VCF file.

    Args:
        primary_vcf (str): Path to the main VCF file.
        to_subtract_vcf (str): Path to the VCF file with variants to subtract.
        output_vcf (str): Path to the output VCF file.
        bed_file (str): Optional path to a BED file for region-based filtering.
        zip_output (bool): If True, the output VCF will be compressed with gzip.
        sample_name (str): Optional new sample name for the output VCF.
    """
    print(f"Parsing variants from {to_subtract_vcf}...")
    variants_to_subtract = parse_vcf(to_subtract_vcf)
    print(f"Found {len(variants_to_subtract)} variants to subtract.")

    bed_regions = None
    if bed_file:
        print(f"Parsing regions from BED file {bed_file}...")
        bed_regions = parse_bed_file(bed_file)
        print(f"Found regions on {len(bed_regions)} chromosomes in BED file.")

    print(f"Filtering variants from {primary_vcf}...")

    try:
        if primary_vcf.endswith('.gz'):
            primary_file = gzip.open(primary_vcf, 'rt')
        else:
            primary_file = open(primary_vcf, 'r')

        # Use gzip.open for writing if the zip_output flag is set
        if zip_output:
            out_file = gzip.open(output_vcf, 'wt')
        else:
            out_file = open(output_vcf, 'w')

        with out_file:
            for line in primary_file:
                # Handle header lines
                if line.startswith('#'):
                    # If this is the sample header line and we want to rename the sample
                    if sample_name and line.startswith('#CHROM'):
                        updated_line = update_sample_name_in_header(line, sample_name)
                        out_file.write(updated_line)
                    else:
                        out_file.write(line)
                    continue

                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    chrom = parts[0]
                    pos = int(parts[1]) # Convert position to integer for comparison
                    ref = parts[3]
                    alt = parts[4]
                    variant_id = (chrom, str(pos), ref, alt)

                    # Check if the variant is in the to_subtract list
                    if variant_id not in variants_to_subtract:
                        # If a BED file was provided, check if the variant is in a region
                        if bed_regions:
                            if is_in_region(chrom, pos, bed_regions):
                                out_file.write(line)
                        else:
                            # If no BED file, write the line directly
                            out_file.write(line)

        primary_file.close()
        print(f"Successfully created filtered VCF at {output_vcf}.")

    except FileNotFoundError:
        print(f"Error: File not found at {primary_vcf}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred while filtering {primary_vcf}: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    """
    Main function to handle command-line arguments and run the subtraction logic.
    """
    parser = argparse.ArgumentParser(description="Subtracts variants from one VCF file if they exist in another, with an optional BED file for region filtering.")
    parser.add_argument("primary_vcf", help="The VCF file to filter.")
    parser.add_argument("to_subtract_vcf", help="The VCF file containing variants to remove .")
    parser.add_argument("output_vcf", help="The name of the output VCF file.")
    parser.add_argument("--bed-file", help="Optional BED file to restrict the output to specific regions.", default=None)
    parser.add_argument("--zip-output", action="store_true", help="Compress the output file with gzip.")
    parser.add_argument("--sample-name", help="Optional new sample name for the output VCF.", default=None)

    args = parser.parse_args()

    subtract_vcf_files(args.primary_vcf, args.to_subtract_vcf, args.output_vcf, args.bed_file, args.zip_output, args.sample_name)

if __name__ == "__main__":
    main()
