#!/usr/bin/env python

# Copyright 2025 - GHGA
# Author: Kuebra Narci - @kubranarci

import matplotlib.pyplot as plt
import os
import argparse
import numpy as np

def plot_svlen_distributions(vcf_files, output_file="svlen_distributions.png"):
    """
    Parses multiple VCF files, extracts SVLEN values, and plots their cumulative distributions
    on a single augmented plot.

    Args:
        vcf_files (list): A list of file paths to the VCF files.
        output_file (str): The name of the output image file for the plot.
    """
    # Create a dictionary to store SVLEN values for each file.
    svlen_data = {}

    # Iterate through each VCF file provided.
    for vcf_file in vcf_files:
        print(f"Processing file: {vcf_file}")

        # Get the base name of the file (e.g., 'manta' or 'lumpy') for the plot legend.
        file_name = os.path.basename(vcf_file).split('.')[2]
        svlen_data[file_name] = []

        # Open and read the VCF file line by line.
        try:
            with open(vcf_file, 'r') as f:
                for line in f:
                    # Skip header lines.
                    if line.startswith('#'):
                        continue

                    parts = line.strip().split('\t')

                    # The INFO column is the 8th column (index 7).
                    info_column = parts[7]

                    # Look for the 'SVLEN=' tag.
                    if 'SVLEN=' in info_column:
                        info_fields = info_column.split(';')
                        for field in info_fields:
                            if field.startswith('SVLEN='):
                                svlen_str = field.split('=')[1]
                                try:
                                    svlen_val = int(svlen_str)
                                    svlen_data[file_name].append(svlen_val)
                                except ValueError:
                                    # Handle cases where SVLEN might be a comma-separated list.
                                    if ',' in svlen_str:
                                        svlen_val = int(svlen_str.split(',')[0])
                                        svlen_data[file_name].append(svlen_val)
                                    else:
                                        print(f"Could not parse SVLEN value: {svlen_str} in {vcf_file}")
        except FileNotFoundError:
            print(f"Error: File not found at {vcf_file}")
            continue

    # Set up the single plot.
    plt.figure(figsize=(12, 8))

    # Define a custom color palette for the plots.
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']

    # Iterate through the collected data and plot each distribution on the same axes.
    for i, (file_name, lengths) in enumerate(svlen_data.items()):
        if not lengths:
            print(f"No SVLEN data found for file: {file_name}")
            continue

        # Sort the SVLEN values for the cumulative distribution plot.
        sorted_lengths = np.sort(lengths)
        y_values = np.arange(len(sorted_lengths)) / float(len(sorted_lengths))

        # Plot the cumulative distribution with a distinct color and label.
        plt.plot(sorted_lengths, y_values, label=file_name, color=colors[i % len(colors)])

    # Set plot properties for the single figure.
    plt.xscale('log')
    plt.xlabel("Structural Variant Length (log scale)")
    plt.ylabel("Cumulative Fraction")
    plt.title("Augmented Structural Variant Length Distributions")
    plt.legend(title="VCF File")
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.xlim([10, 1000000]) # Set a consistent x-axis for comparison.

    # Save the plot to a file.
    plt.savefig(output_file)
    print(f"Plot saved to {output_file}")

# Example usage:
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot SVLEN distributions from one or more VCF files.")
    parser.add_argument('vcf_files', nargs='+', help="One or more VCF files to process.")
    parser.add_argument('--output', '-o', dest='output_file', default='svlen_distributions.png',
                        help="The name of the output image file (default: svlen_distributions.png).")

    args = parser.parse_args()

    # Call the function with the command-line arguments.
    plot_svlen_distributions(args.vcf_files, args.output_file)
