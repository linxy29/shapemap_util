#!/usr/bin/env python3
"""
Analysis functions for bowtie2 alignment pipeline.
These functions are called by run_bowtie_parameter_sweep.py to analyze alignment results.
"""

import pysam
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from typing import List, Tuple


def parse_raw_lengths(stats_file: str) -> Tuple[List[int], List[int]]:
    """Parse raw read length distribution from samtools stats (RL section).
    
    Args:
        stats_file: Path to the samtools stats output file
        
    Returns:
        Tuple of (lengths, counts) lists
    """
    raw_lengths = []
    raw_counts = []
    with open(stats_file) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 3 and parts[0] == "RL":
                raw_lengths.append(int(parts[1]))
                raw_counts.append(int(parts[2]))
    return raw_lengths, raw_counts


def calculate_mapped_lengths(bam_file: str) -> List[int]:
    """Calculate mapped lengths from BAM file using CIGAR strings.
    
    Args:
        bam_file: Path to the sorted BAM file
        
    Returns:
        List of mapped read lengths
    """
    mapped_lengths = []
    
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            if read.is_unmapped:
                continue
            
            # Calculate mapped length from CIGAR
            mapped_length = 0
            for operation, length in read.cigartuples or []:
                if operation in [0, 1, 7, 8]:  # M, I, =, X (consume query)
                    mapped_length += length
            
            mapped_lengths.append(mapped_length)
    
    return mapped_lengths


def plot_distributions(raw_lengths: List[int], raw_counts: List[int], 
                       mapped_lengths: List[int], output_file: str = None) -> None:
    """Plot raw read lengths and mapped lengths distributions using histograms.
    
    Args:
        raw_lengths: List of raw read lengths
        raw_counts: List of counts for each raw read length
        mapped_lengths: List of mapped read lengths
        output_file: Optional path to save the plot
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Convert raw length distribution to expanded list for histogram
    raw_expanded = []
    if raw_lengths and raw_counts:
        for length, count in zip(raw_lengths, raw_counts):
            raw_expanded.extend([length] * count)
    
    # Determine common bin range for consistent x-axis
    all_lengths = []
    if raw_expanded:
        all_lengths.extend(raw_expanded)
    if mapped_lengths:
        all_lengths.extend(mapped_lengths)
    
    if all_lengths:
        min_len, max_len = min(all_lengths), max(all_lengths)
        bins = range(min_len, max_len + 2)
    
    # Plot raw read lengths
    if raw_expanded:
        if not all_lengths:
            bins = range(min(raw_expanded), max(raw_expanded) + 2)
        ax1.hist(raw_expanded, bins=bins, alpha=0.7, color='blue', 
                edgecolor='black', linewidth=0.5)
        ax1.set_xlabel("Length (bp)")
        ax1.set_ylabel("Count")
        ax1.set_title("Raw Read Lengths")
        ax1.grid(True, alpha=0.3)
    else:
        ax1.text(0.5, 0.5, 'No raw length data', ha='center', va='center', 
                transform=ax1.transAxes, fontsize=12)
        ax1.set_title("Raw Read Lengths")
    
    # Plot mapped lengths
    if mapped_lengths:
        if not all_lengths:
            bins = range(min(mapped_lengths), max(mapped_lengths) + 2)
        ax2.hist(mapped_lengths, bins=bins, alpha=0.7, color='orange', 
                edgecolor='black', linewidth=0.5)
        ax2.set_xlabel("Length (bp)")
        ax2.set_ylabel("Count")
        ax2.set_title("Mapped Lengths (CIGAR-based)")
        ax2.grid(True, alpha=0.3)
    else:
        ax2.text(0.5, 0.5, 'No mapped length data', ha='center', va='center', 
                transform=ax2.transAxes, fontsize=12)
        ax2.set_title("Mapped Lengths (CIGAR-based)")
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {output_file}")
    else:
        plt.show()


def save_statistics_to_csv(raw_lengths: List[int], raw_counts: List[int], 
                          mapped_lengths: List[int], output_prefix: str) -> None:
    """Save length distribution data and detailed summary statistics to CSV files.
    
    Args:
        raw_lengths: List of raw read lengths
        raw_counts: List of counts for each raw read length
        mapped_lengths: List of mapped read lengths
        output_prefix: Prefix for output CSV files
    """
    # Save distribution data
    all_data = []
    
    # Add raw read length distribution
    if raw_lengths and raw_counts:
        for length, count in zip(raw_lengths, raw_counts):
            all_data.append({
                'type': 'raw',
                'length': length,
                'count': count
            })
    
    # Add mapped read lengths distribution
    if mapped_lengths:
        mapped_counts = pd.Series(mapped_lengths).value_counts().sort_index()
        for length, count in mapped_counts.items():
            all_data.append({
                'type': 'mapped',
                'length': length,
                'count': count
            })
    
    # Save distribution data to CSV
    if all_data:
        df_distributions = pd.DataFrame(all_data)
        dist_csv_file = f"{output_prefix}_length_distributions.csv"
        df_distributions.to_csv(dist_csv_file, index=False)
        print(f"Length distributions saved to {dist_csv_file}")
    
    # Prepare detailed summary statistics
    summary_data = []
    
    if raw_lengths and raw_counts:
        expanded_raw = []
        for length, count in zip(raw_lengths, raw_counts):
            expanded_raw.extend([length] * count)
        
        if expanded_raw:
            raw_array = np.array(expanded_raw)
            summary_data.append({
                'data_type': 'raw_reads',
                'total_reads': len(expanded_raw),
                'mean_length': np.mean(raw_array),
                'median_length': np.median(raw_array),
                'std_length': np.std(raw_array),
                'min_length': np.min(raw_array),
                'max_length': np.max(raw_array),
                'q5': np.percentile(raw_array, 5),
                'q25': np.percentile(raw_array, 25),
                'q75': np.percentile(raw_array, 75),
                'q95': np.percentile(raw_array, 95),
                'mode_length': pd.Series(expanded_raw).mode().iloc[0] if len(expanded_raw) > 0 else None,
                'length_range': np.max(raw_array) - np.min(raw_array),
                'coefficient_of_variation': np.std(raw_array) / np.mean(raw_array) * 100
            })
    
    if mapped_lengths:
        mapped_array = np.array(mapped_lengths)
        summary_data.append({
            'data_type': 'mapped_reads',
            'total_reads': len(mapped_lengths),
            'mean_length': np.mean(mapped_array),
            'median_length': np.median(mapped_array),
            'std_length': np.std(mapped_array),
            'min_length': np.min(mapped_array),
            'max_length': np.max(mapped_array),
            'q5': np.percentile(mapped_array, 5),
            'q25': np.percentile(mapped_array, 25),
            'q75': np.percentile(mapped_array, 75),
            'q95': np.percentile(mapped_array, 95),
            'mode_length': pd.Series(mapped_lengths).mode().iloc[0] if len(mapped_lengths) > 0 else None,
            'length_range': np.max(mapped_array) - np.min(mapped_array),
            'coefficient_of_variation': np.std(mapped_array) / np.mean(mapped_array) * 100
        })
    
    # Save detailed summary statistics to separate CSV
    if summary_data:
        df_summary = pd.DataFrame(summary_data)
        summary_csv_file = f"{output_prefix}_summary_statistics.csv"
        df_summary.to_csv(summary_csv_file, index=False)
        print(f"Detailed summary statistics saved to {summary_csv_file}")
        
        # Also print summary to console
        print("\nSummary Statistics:")
        print("=" * 50)
        for stats in summary_data:
            print(f"{stats['data_type'].replace('_', ' ').title()}:")
            print(f"  Total reads: {stats['total_reads']:,}")
            print(f"  Mean length: {stats['mean_length']:.2f} bp")
            print(f"  Median length: {stats['median_length']:.2f} bp")
            print(f"  Length range: {stats['min_length']:.0f}-{stats['max_length']:.0f} bp")
            print(f"  Standard deviation: {stats['std_length']:.2f} bp")
            print(f"  Coefficient of variation: {stats['coefficient_of_variation']:.2f}%")
            print()