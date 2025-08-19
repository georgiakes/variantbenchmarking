#!/usr/bin/env python
# A script to parse, analyze, and plot structural variant data from VCF files.

# Copyright 2025 - GHGA
# Author: Kuebra Narci - @kubranarci

import matplotlib.pyplot as plt
import os
import argparse
import numpy as np
import gzip  # Import the gzip module for handling compressed files
import math

def get_sample_name(vcf_file):
    """
    Extracts the sample name from the VCF file name.
    """
    # Get the base name of the file (e.g., 'test1.HG002.manta.tp-base.vcf').
    base_name = os.path.basename(vcf_file)
    # Split the name by the dot '.' and take the second part, which is the sample name in this case.
    # For example, 'test1.HG002.manta.tp-base.vcf' becomes 'manta'.
    return base_name.split('.')[2]

def format_bp_label(value, pos):
    """
    Formats a numeric value into a human-readable string with 'bp', 'kbp', or 'Mbp' suffixes.
    This is used for the x-axis tick labels.
    """
    value = abs(value)
    if value >= 1e6:
        return f'{value / 1e6:.1f} Mbp'
    elif value >= 1e3:
        return f'{value / 1e3:.1f} kbp'
    elif value > 0:
        return f'{int(value)} bp'
    else:
        return '0'

def parse_svlen_data(vcf_files):
    """
    Parses multiple VCF files, extracts SVLEN values, and determines the maximum SV length.

    Returns:
        tuple: A tuple containing a dictionary of SV data and the maximum SV length found.
    """
    svlen_data = {}
    max_svlen = 0

    # Iterate through each VCF file provided.
    for vcf_file in vcf_files:
        print(f"Processing file: {vcf_file}")

        file_name = get_sample_name(vcf_file)
        svlen_data[file_name] = {'positive': [], 'negative': []}

        try:
            # Check if the file is gzipped based on the file extension.
            if vcf_file.endswith('.gz'):
                file_handle = gzip.open(vcf_file, 'rt')
            else:
                file_handle = open(vcf_file, 'r')

            with file_handle as f:
                for line in f:
                    # Skip header lines.
                    if line.startswith('#'):
                        continue

                    parts = line.strip().split('\t')
                    info_column = parts[7]
                    svlen_val = None
                    sv_type = None

                    # --- 1. First priority: Check for SVLEN field ---
                    if 'SVLEN=' in info_column:
                        info_fields = info_column.split(';')
                        for field in info_fields:
                            if field.startswith('SVLEN='):
                                svlen_str = field.split('=')[1]
                                try:
                                    if ',' in svlen_str:
                                        svlen_str = svlen_str.split(',')[0]
                                    svlen_val = int(svlen_str)
                                except ValueError:
                                    print(f"Could not parse SVLEN value: {svlen_str} in {vcf_file}")

                    # --- 2. Second priority: If SVLEN is missing, check for SVTYPE ---
                    if svlen_val is None and 'SVTYPE=' in info_column:
                        info_fields = info_column.split(';')
                        for field in info_fields:
                            if field.startswith('SVTYPE='):
                                sv_type = field.split('=')[1]

                        if sv_type == 'DEL':
                            if 'END=' in info_column:
                                for field in info_fields:
                                    if field.startswith('END='):
                                        end_str = field.split('=')[1]
                                        try:
                                            end_pos = int(end_str)
                                            svlen_val = end_pos - int(parts[1])
                                        except ValueError:
                                            print(f"Could not parse END value: {end_str} in {vcf_file}")
                                            svlen_val = None

                            if svlen_val is None:
                                svlen_val = len(parts[4]) - len(parts[3])

                        elif sv_type == 'INS':
                            svlen_val = len(parts[4]) - len(parts[3])

                    # --- 3. Third priority (ultimate fallback): Infer from allele lengths ---
                    if svlen_val is None:
                        ref_len = len(parts[3])
                        alt_len = len(parts[4])
                        if ref_len > alt_len:
                            svlen_val = -(ref_len - alt_len)
                        elif alt_len > ref_len:
                            svlen_val = alt_len - ref_len

                    if svlen_val is not None and svlen_val != 0:
                        abs_svlen = abs(svlen_val)
                        if abs_svlen > max_svlen:
                            max_svlen = abs_svlen
                        if svlen_val > 0:
                            svlen_data[file_name]['positive'].append(svlen_val)
                        else:
                            svlen_data[file_name]['negative'].append(svlen_val)

        except FileNotFoundError:
            print(f"Error: File not found at {vcf_file}")
            continue

    return svlen_data, max_svlen

def plot_svlen_distributions(svlen_data, max_svlen, output_file, plot_title):
    """
    Plots the SVLEN distributions based on the parsed data using a mirrored linear scale and a log y-axis.
    """
    plt.figure(figsize=(12, 8))

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    line_styles = ['-', '--', ':', '-.']

    # We will use the same number of bins for both positive and negative sides.
    num_bins = 50
    bins = np.linspace(-max_svlen, max_svlen, num_bins * 2)

    # Iterate through the collected data and plot each distribution.
    for i, (file_name, lengths_dict) in enumerate(svlen_data.items()):
        current_color = colors[i % len(colors)]
        current_linestyle = line_styles[i % len(line_styles)]

        # Combine positive and negative lengths into a single array for plotting.
        all_lengths = lengths_dict['positive'] + lengths_dict['negative']

        # Plot the histogram as a line plot.
        if all_lengths:
            counts, bin_edges = np.histogram(all_lengths, bins=bins)
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

            plt.plot(bin_centers, counts, label=file_name,
                     color=current_color, linestyle=current_linestyle)

    plt.yscale('log')
    plt.xlabel("Deletions   |   Insertions")
    plt.ylabel("Count (Log10)")
    plt.title(plot_title)
    plt.legend(title="Method")
    plt.grid(True, which="both", linestyle='--', alpha=0.7)

    # Set the x-limits to be symmetric around zero
    plt.xlim(-max_svlen, max_svlen)

    # Use the custom formatter for the x-axis tick labels
    plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(format_bp_label))

    # Save the plot to a file.
    plt.savefig(output_file)
    print(f"Plot saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot SVLEN distributions from one or more VCF files.")
    parser.add_argument('vcf_files', nargs='+', help="One or more VCF files to process.")
    parser.add_argument('--output', '-o', dest='output_file', default='svlen_distributions.png',
                        help="The name of the output image file (default: svlen_distributions.png).")
    parser.add_argument('--title', '-t', dest='plot_title', default='Structural Variant Length Distributions by Type',
                        help="The title for the plot.")

    args = parser.parse_args()

    sv_data, max_len = parse_svlen_data(args.vcf_files)
    plot_svlen_distributions(sv_data, max_len, args.output_file, args.plot_title)
