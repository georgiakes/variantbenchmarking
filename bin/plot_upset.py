#!/usr/bin/env python

# Copyright 2025 - GHGA
# Author: Kuebra Narci - @kubranarci

"""
Variant Upset Plot Generator
Command-line tool to create upset plots from FP/FN/TP results
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_contents
import sys
import os

def detect_tools_from_header(header):
    """Detect tool names from VCF header columns using simple pattern logic"""
    tools = []

    for col in header:
        if any(x in col for x in ['_GT', '.GT']):
            # Split by both dots and underscores to handle both patterns
            parts = []
            for separator in ['.', '_']:
                if separator in col:
                    parts = col.split(separator)
                    break

            if not parts:  # If no separators found, skip
                continue

            # Remove GT and any empty parts
            clean_parts = [part for part in parts if part and part != 'GT']

            # Apply the pattern logic
            if len(clean_parts) >= 3:
                # Pattern: "test1.HG002.delly.falsenegatives.vcf_GT" -> use 3rd element
                tool_name = clean_parts[2]
            elif len(clean_parts) >= 2:
                # Pattern: "HG002.delly_GT" -> use 2nd element
                tool_name = clean_parts[1]
            else:
                # Fallback: use first element if only one exists
                tool_name = clean_parts[0] if clean_parts else None

            if tool_name and tool_name not in tools:
                tools.append(tool_name)

    return sorted(tools)

def parse_vcf_file(file_path, category_name):
    """Parse VCF file and extract tool detection information"""
    try:
        if not os.path.exists(file_path):
            print(f"Error: File not found: {file_path}")
            return pd.DataFrame()

        # Read the file
        with open(file_path, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]

        if not lines:
            return pd.DataFrame()

        # Find header line
        header_line = None
        for line in lines:
            if line.startswith('CHROM,POS') or ('CHROM' in line and 'POS' in line and 'REF' in line and 'ALT' in line):
                header_line = line
                break

        if not header_line:
            print(f"Error: No header found in {file_path}")
            return pd.DataFrame()

        header = header_line.split(',')
        tools = detect_tools_from_header(header)

        if not tools:
            print(f"Warning: No tools detected in {file_path}")
            return pd.DataFrame()

        print(f"Detected tools in {category_name}: {', '.join(tools)}")

        # Parse data
        data = []
        for line in lines:
            if line.startswith('CHROM,POS') or line.startswith('#'):
                continue
            fields = line.split(',')
            if len(fields) >= len(header):
                data.append(fields)

        df = pd.DataFrame(data, columns=header)

        # Create variant ID
        df['variant_id'] = df['CHROM'] + ':' + df['POS'] + ':' + df['REF'] + ':' + df['ALT']

        # Check tool detection
        tool_detection = []
        for _, row in df.iterrows():
            variant_id = row['variant_id']
            detected_tools = []
            for tool in tools:
                # Find GT column for this tool
                gt_columns = [col for col in header if tool in col and ('_GT' in col or '.GT' in col)]
                if gt_columns:
                    gt_value = row[gt_columns[0]]
                    # Check if variant was detected (not 0/0, ./. or similar)
                    if (gt_value not in ['0/0', './.', '0|0', '.|.', ''] and
                        '0' not in str(gt_value) and gt_value != '0'):
                        detected_tools.append(tool)

            if detected_tools:
                tool_detection.append({'variant_id': variant_id, 'tools': detected_tools})

        if not tool_detection:
            return pd.DataFrame()

        # Create indicator matrix
        tool_df = pd.DataFrame(tool_detection)
        exploded = tool_df.explode('tools')
        indicator_matrix = pd.crosstab(exploded['variant_id'], exploded['tools'])
        indicator_matrix = (indicator_matrix > 0).astype(int)
        indicator_matrix['category'] = category_name

        return indicator_matrix

    except Exception as e:
        print(f"Error parsing {file_path}: {str(e)}")
        return pd.DataFrame()

def create_grouped_upset_plots(dataframes, output_file=None, title="Variant Detection Upset Plots"):
    """Create two upset plots: TP_base+FN and TP_comp+FP in the same PNG"""
    if not dataframes:
        print("Error: No data to plot")
        return

    # Combine all data
    combined_df = pd.concat(dataframes)

    # Get all unique tools
    all_tools = [col for col in combined_df.columns if col != 'category']

    # Create data for first group: TP_base + FN
    group1_data = {}
    group1_categories = ['TP_Base', 'FN']

    for category in group1_categories:
        category_df = combined_df[combined_df['category'] == category]
        if not category_df.empty:
            for tool in all_tools:
                if tool in category_df.columns:
                    variants = set(category_df[category_df[tool] > 0].index)
                    if variants:
                        group1_data[f"{category}_{tool}"] = variants

    # Create data for second group: TP_comp + FP
    group2_data = {}
    group2_categories = ['TP_Comp', 'FP']

    for category in group2_categories:
        category_df = combined_df[combined_df['category'] == category]
        if not category_df.empty:
            for tool in all_tools:
                if tool in category_df.columns:
                    variants = set(category_df[category_df[tool] > 0].index)
                    if variants:
                        group2_data[f"{category}_{tool}"] = variants

    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))

    # Plot first group (TP_base + FN)
    if group1_data:
        # Create a separate figure for the first upset plot
        fig1, ax_temp = plt.subplots(figsize=(10, 6))
        upset1 = UpSet(from_contents(group1_data),
                      subset_size='count',
                      show_counts=True,
                      sort_by='cardinality')
        upset1.plot(fig=fig1)
        fig1.canvas.draw()

        # Convert the plot to an image and display it on the first subplot
        ax1.imshow(fig1.canvas.renderer.buffer_rgba())
        ax1.axis('off')
        ax1.set_title('TP_Base + FN', fontsize=14, fontweight='bold')
        plt.close(fig1)
    else:
        ax1.text(0.5, 0.5, 'No data for TP_Base + FN',
                horizontalalignment='center', verticalalignment='center',
                transform=ax1.transAxes, fontsize=12)
        ax1.set_title('TP_Base + FN (No Data)', fontsize=14, fontweight='bold')
        ax1.set_xticks([])
        ax1.set_yticks([])

    # Plot second group (TP_comp + FP)
    if group2_data:
        # Create a separate figure for the second upset plot
        fig2, ax_temp = plt.subplots(figsize=(10, 6))
        upset2 = UpSet(from_contents(group2_data),
                      subset_size='count',
                      show_counts=True,
                      sort_by='cardinality')
        upset2.plot(fig=fig2)
        fig2.canvas.draw()

        # Convert the plot to an image and display it on the second subplot
        ax2.imshow(fig2.canvas.renderer.buffer_rgba())
        ax2.axis('off')
        ax2.set_title('TP_Comp + FP', fontsize=14, fontweight='bold')
        plt.close(fig2)
    else:
        ax2.text(0.5, 0.5, 'No data for TP_Comp + FP',
                horizontalalignment='center', verticalalignment='center',
                transform=ax2.transAxes, fontsize=12)
        ax2.set_title('TP_Comp + FP (No Data)', fontsize=14, fontweight='bold')
        ax2.set_xticks([])
        ax2.set_yticks([])

    plt.suptitle(title, fontsize=16, y=0.98)
    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {output_file}")
    else:
        plt.show()

def main():
    parser = argparse.ArgumentParser(description='Create upset plots from variant benchmarking results')
    parser.add_argument('--fp', help='False positives VCF file')
    parser.add_argument('--fn', help='False negatives VCF file')
    parser.add_argument('--tp-base', help='True positives base VCF file')
    parser.add_argument('--tp-comp', help='True positives comparison VCF file')
    parser.add_argument('--output', '-o', help='Output plot file (optional)')
    parser.add_argument('--title', help='Plot title', default='Variant Detection Upset Plots')

    args = parser.parse_args()

    if not any([args.fp, args.fn, args.tp_base, args.tp_comp]):
        parser.print_help()
        sys.exit(1)

    # Parse all provided files
    dataframes = []

    if args.fp:
        print(f"Processing false positives: {args.fp}")
        fp_df = parse_vcf_file(args.fp, 'FP')
        if not fp_df.empty:
            dataframes.append(fp_df)

    if args.fn:
        print(f"Processing false negatives: {args.fn}")
        fn_df = parse_vcf_file(args.fn, 'FN')
        if not fn_df.empty:
            dataframes.append(fn_df)

    if args.tp_base:
        print(f"Processing true positives base: {args.tp_base}")
        tp_base_df = parse_vcf_file(args.tp_base, 'TP_Base')
        if not tp_base_df.empty:
            dataframes.append(tp_base_df)

    if args.tp_comp:
        print(f"Processing true positives comparison: {args.tp_comp}")
        tp_comp_df = parse_vcf_file(args.tp_comp, 'TP_Comp')
        if not tp_comp_df.empty:
            dataframes.append(tp_comp_df)

    if not dataframes:
        print("Error: No valid data found in any input files")
        sys.exit(1)

    # Create the grouped upset plots
    create_grouped_upset_plots(dataframes, args.output, args.title)

if __name__ == "__main__":
    main()
