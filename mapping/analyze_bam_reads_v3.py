#!/usr/bin/env python3
"""
BAM Read Alignment Analyzer - Memory Efficient Version
Xinyi Lin, 202509

A memory-efficient tool for analyzing BAM file read alignments with focus on:
- Multi-mapping read analysis with minimal memory footprint
- Transcript assignment statistics using online computation
- Read filtering by alignment quality
- Extraction of uniquely mapped reads
- Extract the primary alignments

Version: v3.0 (Memory Efficient - Stream Processing)

Key improvements over v2:
- Stream processing: Process one read at a time, write immediately, discard
- Online statistics: No array storage, constant memory usage
- Auto-check BAM sort order (must be sorted by read name)
- ~99.9% memory reduction for large files
- Same speed as v2, optimized for throughput

Usage:
    # Check if BAM is sorted by name first:
    samtools view -H input.bam | grep '@HD'
    # If not sorted by name: samtools sort -n input.bam -o input.sorted.bam

    # Statistics only (minimal memory)
    python analyze_bam_reads_v3_memory_efficient.py input.sorted.bam -o output_prefix

    # Full analysis with per-read TSV
    python analyze_bam_reads_v3_memory_efficient.py input.sorted.bam -o output_prefix --full-analysis

    # Extract unique reads
    python analyze_bam_reads_v3_memory_efficient.py input.sorted.bam -o output_prefix --extract-unique
"""

# ============================================================================
# IMPORTS
# ============================================================================

import argparse
import csv
import logging
import sys
import gc
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set, TextIO

import matplotlib.pyplot as plt
import numpy as np
import pysam

# ============================================================================
# CONSTANTS
# ============================================================================

# Default parameters
DEFAULT_OUTPUT_PREFIX = 'read_alignment_analysis'
DEFAULT_MIN_LENGTH = 0
PROGRESS_INTERVAL = 1000000

# Plot settings
PLOT_DPI = 150
PLOT_FIGSIZE = (12, 8)

# ============================================================================
# LOGGING CONFIGURATION
# ============================================================================

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# ============================================================================
# DATA STRUCTURES
# ============================================================================

@dataclass
class ReadAlignment:
    """
    Stores alignment information for a single read.
    Memory optimized: Only stores transcript counts, no full list.
    """
    read_id: str
    alignment_count: int = 0
    transcript_counts: Counter = field(default_factory=Counter)

    def add_alignment(self, transcript: str) -> None:
        """Add an alignment record for this read."""
        self.alignment_count += 1
        self.transcript_counts[transcript] += 1

    def get_max_transcript(self) -> Tuple[str, int]:
        """Get the transcript with the most alignments."""
        if self.transcript_counts:
            return self.transcript_counts.most_common(1)[0]
        return '', 0

    def get_unique_transcripts(self) -> List[str]:
        """Get list of unique transcripts this read maps to."""
        return list(self.transcript_counts.keys())

    def get_max_ratio(self) -> float:
        """Calculate the ratio of max transcript count to total alignments."""
        if self.alignment_count > 0:
            _, max_count = self.get_max_transcript()
            return max_count / self.alignment_count
        return 0.0

    def is_uniquely_mapped(self) -> bool:
        """Check if this read maps to only one unique transcript."""
        return len(self.transcript_counts) == 1


@dataclass
class OnlineStatistics:
    """
    Online statistics computation without storing all values.
    Computes running statistics with constant memory usage.
    """
    # Basic counters
    total_reads: int = 0
    single_aligned: int = 0
    multi_aligned: int = 0
    uniquely_mapped: int = 0

    # Ratio statistics (online computation)
    ratio_sum: float = 0.0
    ratio_min: float = float('inf')
    ratio_max: float = float('-inf')

    # Category counts
    ratio_eq_1: int = 0  # Exactly 1.0
    ratio_08_1: int = 0  # [0.8, 1.0)
    ratio_05_08: int = 0  # [0.5, 0.8)
    ratio_0_05: int = 0  # (0, 0.5)

    # Alignment distribution
    alignment_distribution: Counter = field(default_factory=Counter)

    def update(self, read_data: ReadAlignment) -> None:
        """Update statistics with a new read."""
        self.total_reads += 1

        # Alignment type counts
        if read_data.alignment_count == 1:
            self.single_aligned += 1
        else:
            self.multi_aligned += 1

        # Unique mapping
        if read_data.is_uniquely_mapped():
            self.uniquely_mapped += 1

        # Ratio statistics
        ratio = read_data.get_max_ratio()
        self.ratio_sum += ratio
        self.ratio_min = min(self.ratio_min, ratio)
        self.ratio_max = max(self.ratio_max, ratio)

        # Ratio categories
        if ratio == 1.0:
            self.ratio_eq_1 += 1
        elif 0.8 <= ratio < 1.0:
            self.ratio_08_1 += 1
        elif 0.5 <= ratio < 0.8:
            self.ratio_05_08 += 1
        elif 0 < ratio < 0.5:
            self.ratio_0_05 += 1

        # Alignment distribution
        self.alignment_distribution[read_data.alignment_count] += 1

    @property
    def mean_ratio(self) -> float:
        """Calculate mean ratio."""
        return self.ratio_sum / self.total_reads if self.total_reads > 0 else 0.0

    def get_summary_dict(self) -> Dict:
        """Get summary statistics as a dictionary."""
        total = self.total_reads
        return {
            'total_reads': total,
            'single_aligned': self.single_aligned,
            'multi_aligned': self.multi_aligned,
            'uniquely_mapped': self.uniquely_mapped,
            'single_aligned_pct': self._safe_percentage(self.single_aligned, total),
            'multi_aligned_pct': self._safe_percentage(self.multi_aligned, total),
            'uniquely_mapped_pct': self._safe_percentage(self.uniquely_mapped, total),
            'mean_ratio': self.mean_ratio,
            'min_ratio': self.ratio_min if total > 0 else 0,
            'max_ratio': self.ratio_max if total > 0 else 0,
        }

    @staticmethod
    def _safe_percentage(count: int, total: int) -> float:
        """Calculate percentage safely."""
        return (count / total * 100) if total > 0 else 0.0

    def print_summary(self) -> None:
        """Print formatted summary statistics."""
        summary = self.get_summary_dict()

        print("\n" + "=" * 60)
        print("ALIGNMENT SUMMARY STATISTICS")
        print("=" * 60)

        print(f"\nTotal unique reads: {summary['total_reads']:,}")

        if summary['total_reads'] > 0:
            print(f"\nAlignment types:")
            print(f"  Single-aligned reads:     {summary['single_aligned']:,} ({summary['single_aligned_pct']:.2f}%)")
            print(f"  Multi-aligned reads:      {summary['multi_aligned']:,} ({summary['multi_aligned_pct']:.2f}%)")
            print(f"  Uniquely mapped reads:    {summary['uniquely_mapped']:,} ({summary['uniquely_mapped_pct']:.2f}%)")

            print(f"\nRatio statistics:")
            print(f"  Mean ratio:      {summary['mean_ratio']:.4f}")
            print(f"  Min ratio:       {summary['min_ratio']:.4f}")
            print(f"  Max ratio:       {summary['max_ratio']:.4f}")

            print(f"\nDominance categories:")
            total = summary['total_reads']
            print(f"  Uniquely mapped (ratio = 1.0):  {self.ratio_eq_1:,} ({self._safe_percentage(self.ratio_eq_1, total):.2f}%)")
            print(f"  High dominance [0.8-1.0):       {self.ratio_08_1:,} ({self._safe_percentage(self.ratio_08_1, total):.2f}%)")
            print(f"  Moderate dominance [0.5-0.8):   {self.ratio_05_08:,} ({self._safe_percentage(self.ratio_05_08, total):.2f}%)")
            print(f"  Low dominance (0-0.5):          {self.ratio_0_05:,} ({self._safe_percentage(self.ratio_0_05, total):.2f}%)")

            if self.alignment_distribution:
                print("\nTop 10 alignment count distribution:")
                for count, freq in self.alignment_distribution.most_common(10):
                    print(f"  {count:3d} alignments: {freq:,} reads")

        print("=" * 60)


# ============================================================================
# BAM SORT ORDER CHECKER
# ============================================================================

class BAMSortChecker:
    """Check if BAM file is sorted by query name."""

    @staticmethod
    def check_sort_order(bam_file: Path) -> Tuple[bool, str]:
        """
        Check if BAM file is sorted by query name.

        Returns:
            Tuple of (is_sorted_by_name, sort_order_description)
        """
        try:
            with pysam.AlignmentFile(str(bam_file), "rb") as bamfile:
                header = bamfile.header

                # Check header for sort order
                if 'HD' in header:
                    hd = header['HD']
                    if 'SO' in hd:
                        sort_order = hd['SO']
                        if sort_order == 'queryname':
                            return True, "queryname (correctly sorted)"
                        else:
                            return False, f"{sort_order} (NOT sorted by query name)"

                # If no sort order in header, check first few reads
                logger.warning("No sort order in BAM header, checking read order...")
                return BAMSortChecker._verify_read_order(bamfile)

        except Exception as e:
            logger.error(f"Error checking BAM sort order: {e}")
            return False, "unknown (error reading header)"

    @staticmethod
    def _verify_read_order(bamfile, check_limit: int = 10000) -> Tuple[bool, str]:
        """
        Verify read order by checking if consecutive reads with same name are grouped.

        Args:
            bamfile: Open BAM file handle
            check_limit: Number of reads to check

        Returns:
            Tuple of (likely_sorted_by_name, description)
        """
        prev_read_name = None
        read_name_changes = 0
        total_checked = 0

        for alignment in bamfile:
            current_read_name = alignment.query_name

            if prev_read_name is not None and current_read_name != prev_read_name:
                read_name_changes += 1

            prev_read_name = current_read_name
            total_checked += 1

            if total_checked >= check_limit:
                break

        if total_checked == 0:
            return False, "empty BAM file"

        # If reads change names very frequently, likely sorted by name
        # (each read's alignments are grouped together)
        if total_checked > 100:
            # Heuristic: if we see many name changes, likely sorted by name
            if read_name_changes > total_checked * 0.5:
                return True, "likely queryname (verified by sampling)"
            else:
                return False, "likely coordinate sorted (not queryname)"

        return False, "insufficient data to verify sort order"


# ============================================================================
# STREAMING BAM PROCESSOR
# ============================================================================

class StreamingBAMProcessor:
    """
    Memory-efficient BAM processor using streaming approach.
    Processes one read at a time, writes immediately, discards from memory.
    """

    def __init__(self, bam_file: str, min_mapped_length: int = 0):
        """
        Initialize streaming BAM processor.

        Args:
            bam_file: Path to input BAM file (must be sorted by query name)
            min_mapped_length: Minimum mapped length threshold
        """
        self.bam_file = Path(bam_file)
        self.min_mapped_length = min_mapped_length
        self.stats = OnlineStatistics()
        self.processing_stats = {
            'total_alignments': 0,
            'filtered_alignments': 0,
            'unmapped_alignments': 0
        }

        if not self.bam_file.exists():
            raise FileNotFoundError(f"BAM file not found: {bam_file}")

        # Check sort order
        self._check_sort_order()

        logger.info(f"Initialized streaming BAM processor for: {self.bam_file}")

    def _check_sort_order(self) -> None:
        """Check if BAM is sorted by query name and raise error if not."""
        is_sorted, sort_desc = BAMSortChecker.check_sort_order(self.bam_file)

        logger.info(f"BAM sort order: {sort_desc}")

        if not is_sorted:
            error_msg = f"""
BAM file is NOT sorted by query name (current: {sort_desc})

This tool requires name-sorted BAM for memory-efficient streaming.

To sort your BAM file by read name, run:
    samtools sort -n {self.bam_file} -o {self.bam_file.stem}.sorted.bam

Then use the sorted file with this tool.
"""
            raise ValueError(error_msg)

    def process_streaming(self, tsv_writer: Optional[csv.DictWriter] = None,
                         unique_bam_writer: Optional[pysam.AlignmentFile] = None,
                         progress_interval: int = PROGRESS_INTERVAL) -> None:
        """
        Process BAM file in streaming mode.

        Args:
            tsv_writer: Optional CSV writer for per-read output
            unique_bam_writer: Optional BAM writer for unique reads
            progress_interval: Progress update interval
        """
        logger.info(f"Processing BAM file in streaming mode")
        logger.info(f"Min length filter: {self.min_mapped_length} bp")

        current_read = None
        current_read_id = None
        reads_processed = 0

        try:
            with pysam.AlignmentFile(str(self.bam_file), "rb") as bamfile:
                for i, alignment in enumerate(bamfile):
                    self.processing_stats['total_alignments'] += 1

                    # Progress reporting
                    if i % progress_interval == 0 and i > 0:
                        logger.info(f"Processed {i:,} alignments, {reads_processed:,} reads...")

                    # Skip unmapped reads
                    if not alignment.reference_name:
                        self.processing_stats['unmapped_alignments'] += 1
                        continue

                    # Apply length filter
                    mapped_length = alignment.reference_length if alignment.reference_length is not None else 0
                    if mapped_length < self.min_mapped_length:
                        self.processing_stats['filtered_alignments'] += 1
                        continue

                    # Check if we've moved to a new read
                    read_id = alignment.query_name
                    if read_id != current_read_id:
                        # Process and write previous read
                        if current_read is not None:
                            self._process_completed_read(current_read, tsv_writer, unique_bam_writer)
                            reads_processed += 1

                        # Start new read
                        current_read = ReadAlignment(read_id)
                        current_read_id = read_id

                    # Add alignment to current read
                    current_read.add_alignment(alignment.reference_name)

                    # Store alignment for potential BAM writing
                    if unique_bam_writer is not None:
                        if not hasattr(current_read, '_alignments'):
                            current_read._alignments = []
                        current_read._alignments.append(alignment)

                # Process last read
                if current_read is not None:
                    self._process_completed_read(current_read, tsv_writer, unique_bam_writer)
                    reads_processed += 1

        except Exception as e:
            logger.error(f"Error processing BAM file: {e}")
            raise

        logger.info(f"Streaming processing complete: {reads_processed:,} reads processed")
        self._log_processing_summary()

    def _process_completed_read(self, read_data: ReadAlignment,
                               tsv_writer: Optional[csv.DictWriter],
                               unique_bam_writer: Optional[pysam.AlignmentFile]) -> None:
        """
        Process a completed read: update stats, write to TSV/BAM, then discard.

        Args:
            read_data: Completed read alignment data
            tsv_writer: Optional TSV writer
            unique_bam_writer: Optional BAM writer for unique reads
        """
        # Update online statistics
        self.stats.update(read_data)

        # Write to TSV if requested
        if tsv_writer is not None:
            self._write_read_to_tsv(read_data, tsv_writer)

        # Write to unique BAM if requested and read is uniquely mapped
        if unique_bam_writer is not None and read_data.is_uniquely_mapped():
            self._write_unique_read_to_bam(read_data, unique_bam_writer)

        # Read is now processed and can be discarded (Python GC will handle it)

    def _write_read_to_tsv(self, read_data: ReadAlignment, writer: csv.DictWriter) -> None:
        """Write read data to TSV."""
        max_transcript, max_count = read_data.get_max_transcript()
        unique_transcripts = read_data.get_unique_transcripts()

        writer.writerow({
            'read_id': read_data.read_id,
            'alignment_count': read_data.alignment_count,
            'unique_transcript_count': len(unique_transcripts),
            'is_uniquely_mapped': 'Yes' if read_data.is_uniquely_mapped() else 'No',
            'max_transcript': max_transcript,
            'max_count': max_count,
            'max_ratio': f"{read_data.get_max_ratio():.4f}",
            'unique_transcripts': ','.join(unique_transcripts)
        })

    def _write_unique_read_to_bam(self, read_data: ReadAlignment,
                                  bam_writer: pysam.AlignmentFile) -> None:
        """Write uniquely mapped read's primary alignment to BAM."""
        if hasattr(read_data, '_alignments'):
            for alignment in read_data._alignments:
                # Only write primary alignments
                if not alignment.is_secondary and not alignment.is_supplementary:
                    bam_writer.write(alignment)

    def _log_processing_summary(self) -> None:
        """Log processing summary statistics."""
        stats = self.processing_stats
        total = stats['total_alignments']

        logger.info(f"Processing summary:")
        logger.info(f"  Total alignments: {total:,}")

        if total > 0:
            filtered_pct = stats['filtered_alignments'] / total * 100
            unmapped_pct = stats['unmapped_alignments'] / total * 100
            passed = total - stats['filtered_alignments'] - stats['unmapped_alignments']
            passed_pct = passed / total * 100

            logger.info(f"  Unmapped: {stats['unmapped_alignments']:,} ({unmapped_pct:.2f}%)")
            logger.info(f"  Filtered (< {self.min_mapped_length} bp): {stats['filtered_alignments']:,} ({filtered_pct:.2f}%)")
            logger.info(f"  Passed: {passed:,} ({passed_pct:.2f}%)")


# ============================================================================
# FILE I/O MODULE
# ============================================================================

class FileHandler:
    """Handles file I/O operations."""

    @staticmethod
    def get_tsv_writer(output_file: Path) -> Tuple[TextIO, csv.DictWriter]:
        """
        Create TSV writer for streaming output.

        Returns:
            Tuple of (file_handle, csv_writer)
        """
        fieldnames = [
            'read_id',
            'alignment_count',
            'unique_transcript_count',
            'is_uniquely_mapped',
            'max_transcript',
            'max_count',
            'max_ratio',
            'unique_transcripts'
        ]

        file_handle = open(output_file, 'w', newline='')
        writer = csv.DictWriter(file_handle, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        return file_handle, writer

    @staticmethod
    def write_summary_stats(stats: OnlineStatistics, output_file: Path) -> None:
        """Write summary statistics to file."""
        logger.info(f"Writing summary statistics to: {output_file}")

        try:
            with open(output_file, 'w') as f:
                summary = stats.get_summary_dict()

                f.write("ALIGNMENT SUMMARY STATISTICS\n")
                f.write("=" * 60 + "\n\n")

                f.write(f"Total unique reads: {summary['total_reads']:,}\n\n")

                if summary['total_reads'] > 0:
                    f.write("Alignment types:\n")
                    f.write(f"  Single-aligned reads:     {summary['single_aligned']:,} ({summary['single_aligned_pct']:.2f}%)\n")
                    f.write(f"  Multi-aligned reads:      {summary['multi_aligned']:,} ({summary['multi_aligned_pct']:.2f}%)\n")
                    f.write(f"  Uniquely mapped reads:    {summary['uniquely_mapped']:,} ({summary['uniquely_mapped_pct']:.2f}%)\n\n")

                    f.write("Ratio statistics:\n")
                    f.write(f"  Mean ratio:      {summary['mean_ratio']:.4f}\n")
                    f.write(f"  Min ratio:       {summary['min_ratio']:.4f}\n")
                    f.write(f"  Max ratio:       {summary['max_ratio']:.4f}\n\n")

                    f.write("Dominance categories:\n")
                    total = summary['total_reads']
                    f.write(f"  Uniquely mapped (ratio = 1.0):  {stats.ratio_eq_1:,} ({stats._safe_percentage(stats.ratio_eq_1, total):.2f}%)\n")
                    f.write(f"  High dominance [0.8-1.0):       {stats.ratio_08_1:,} ({stats._safe_percentage(stats.ratio_08_1, total):.2f}%)\n")
                    f.write(f"  Moderate dominance [0.5-0.8):   {stats.ratio_05_08:,} ({stats._safe_percentage(stats.ratio_05_08, total):.2f}%)\n")
                    f.write(f"  Low dominance (0-0.5):          {stats.ratio_0_05:,} ({stats._safe_percentage(stats.ratio_0_05, total):.2f}%)\n\n")

                    if stats.alignment_distribution:
                        f.write("Top 10 alignment count distribution:\n")
                        for count, freq in stats.alignment_distribution.most_common(10):
                            f.write(f"  {count:3d} alignments: {freq:,} reads\n")

            logger.info("Successfully wrote summary statistics")

        except Exception as e:
            logger.error(f"Error writing summary file: {e}")
            raise


# ============================================================================
# VISUALIZATION MODULE
# ============================================================================

class Visualizer:
    """Creates visualization plots with online statistics (no histograms)."""

    @staticmethod
    def plot_summary(stats: OnlineStatistics, output_prefix: str) -> None:
        """
        Create summary plots from online statistics.

        Args:
            stats: OnlineStatistics object
            output_prefix: Prefix for output files
        """
        if stats.total_reads == 0:
            logger.warning("No data available for plotting")
            return

        logger.info("Creating summary visualization...")

        try:
            fig, axes = plt.subplots(2, 2, figsize=PLOT_FIGSIZE)

            # Plot 1: Alignment type distribution
            Visualizer._plot_alignment_types(axes[0, 0], stats)

            # Plot 2: Dominance categories
            Visualizer._plot_dominance_categories(axes[0, 1], stats)

            # Plot 3: Top alignment counts
            Visualizer._plot_alignment_distribution(axes[1, 0], stats)

            # Plot 4: Summary statistics table
            Visualizer._plot_summary_table(axes[1, 1], stats)

            plt.tight_layout()
            output_file = f'{output_prefix}_summary.png'
            plt.savefig(output_file, dpi=PLOT_DPI, bbox_inches='tight')
            plt.close()

            logger.info(f"Saved summary plot to: {output_file}")

        except Exception as e:
            logger.error(f"Error creating plots: {e}")

    @staticmethod
    def _plot_alignment_types(ax, stats: OnlineStatistics) -> None:
        """Plot alignment type distribution."""
        categories = ['Single-aligned', 'Multi-aligned']
        counts = [stats.single_aligned, stats.multi_aligned]
        colors = ['steelblue', 'darkorange']

        bars = ax.bar(categories, counts, color=colors, alpha=0.7, edgecolor='black')
        ax.set_ylabel('Number of Reads')
        ax.set_title('Alignment Type Distribution')
        ax.grid(True, alpha=0.3, axis='y')

        # Add value labels on bars
        for bar, count in zip(bars, counts):
            height = bar.get_height()
            pct = stats._safe_percentage(count, stats.total_reads)
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{count:,}\n({pct:.1f}%)',
                   ha='center', va='bottom', fontsize=9)

    @staticmethod
    def _plot_dominance_categories(ax, stats: OnlineStatistics) -> None:
        """Plot dominance category distribution."""
        categories = ['Unique\n(1.0)', 'High\n[0.8-1.0)', 'Moderate\n[0.5-0.8)', 'Low\n(0-0.5)']
        counts = [stats.ratio_eq_1, stats.ratio_08_1, stats.ratio_05_08, stats.ratio_0_05]
        colors = ['lightgreen', 'lightyellow', 'lightcoral', 'lightgray']

        bars = ax.bar(categories, counts, color=colors, alpha=0.7, edgecolor='black')
        ax.set_ylabel('Number of Reads')
        ax.set_title('Transcript Dominance Categories')
        ax.grid(True, alpha=0.3, axis='y')

        # Add value labels on bars
        for bar, count in zip(bars, counts):
            height = bar.get_height()
            pct = stats._safe_percentage(count, stats.total_reads)
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{count:,}\n({pct:.1f}%)',
                   ha='center', va='bottom', fontsize=8)

    @staticmethod
    def _plot_alignment_distribution(ax, stats: OnlineStatistics) -> None:
        """Plot top alignment count distribution."""
        top_10 = stats.alignment_distribution.most_common(10)

        if top_10:
            counts, freqs = zip(*top_10)
            bars = ax.bar(range(len(counts)), freqs, alpha=0.7, color='purple', edgecolor='black')
            ax.set_xticks(range(len(counts)))
            ax.set_xticklabels([str(c) for c in counts])
            ax.set_xlabel('Number of Alignments per Read')
            ax.set_ylabel('Number of Reads')
            ax.set_title('Top 10 Alignment Count Distribution')
            ax.grid(True, alpha=0.3, axis='y')

            # Add value labels on bars (only for bars tall enough)
            for bar in bars:
                height = bar.get_height()
                if height > max(freqs) * 0.05:  # Only label if >5% of max
                    ax.text(bar.get_x() + bar.get_width()/2., height,
                           f'{int(height):,}',
                           ha='center', va='bottom', fontsize=7)
        else:
            ax.text(0.5, 0.5, 'No data available',
                   ha='center', va='center', transform=ax.transAxes)

    @staticmethod
    def _plot_summary_table(ax, stats: OnlineStatistics) -> None:
        """Plot summary statistics as a table."""
        ax.axis('off')

        summary = stats.get_summary_dict()

        table_data = [
            ['Metric', 'Value'],
            ['Total Reads', f"{summary['total_reads']:,}"],
            ['Uniquely Mapped', f"{summary['uniquely_mapped']:,}\n({summary['uniquely_mapped_pct']:.2f}%)"],
            ['', ''],
            ['Mean Ratio', f"{summary['mean_ratio']:.4f}"],
            ['Min Ratio', f"{summary['min_ratio']:.4f}"],
            ['Max Ratio', f"{summary['max_ratio']:.4f}"],
        ]

        table = ax.table(cellText=table_data, cellLoc='left', loc='center',
                        colWidths=[0.6, 0.4])
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1, 2)

        # Style header row
        for i in range(2):
            table[(0, i)].set_facecolor('#4472C4')
            table[(0, i)].set_text_props(weight='bold', color='white')

        ax.set_title('Summary Statistics', fontsize=12, weight='bold', pad=20)


# ============================================================================
# MAIN APPLICATION
# ============================================================================

class BAMAnalyzerApp:
    """Main application controller."""

    def __init__(self, args: argparse.Namespace):
        """Initialize the application."""
        self.args = args
        self.processor = None

    def run(self) -> None:
        """Run the complete analysis pipeline."""
        try:
            # Initialize processor (will check sort order)
            self.processor = StreamingBAMProcessor(
                self.args.bam_file,
                self.args.min_length
            )

            # Run appropriate mode
            if self.args.full_analysis:
                self._run_full_analysis()
            elif self.args.extract_unique:
                self._run_unique_extraction()
            else:
                self._run_statistics_only()

            logger.info("Analysis complete!")

        except Exception as e:
            logger.error(f"Analysis failed: {e}")
            raise

    def _run_full_analysis(self) -> None:
        """Run full analysis mode with TSV output."""
        logger.info("Running FULL ANALYSIS mode (with per-read TSV)")

        # Open TSV file for streaming write
        output_file = Path(f'{self.args.output}_read_alignment_analysis.tsv')
        logger.info(f"Writing per-read analysis to: {output_file}")

        tsv_handle, tsv_writer = FileHandler.get_tsv_writer(output_file)

        try:
            # Process with TSV writing
            self.processor.process_streaming(tsv_writer=tsv_writer)
        finally:
            tsv_handle.close()

        # Generate reports
        self._generate_reports()

    def _run_unique_extraction(self) -> None:
        """Run unique read extraction mode."""
        logger.info("Running UNIQUE EXTRACTION mode")

        # Open BAM file for writing unique reads
        output_bam = Path(f'{self.args.output}.bam')
        logger.info(f"Writing unique reads to: {output_bam}")

        with pysam.AlignmentFile(str(self.processor.bam_file), "rb") as template:
            with pysam.AlignmentFile(str(output_bam), "wb", template=template) as bam_writer:
                # Process with BAM writing
                self.processor.process_streaming(unique_bam_writer=bam_writer)

        # Generate reports
        self._generate_reports()

    def _run_statistics_only(self) -> None:
        """Run statistics only mode."""
        logger.info("Running STATISTICS mode")

        # Process without writing
        self.processor.process_streaming()

        # Generate reports
        self._generate_reports()

    def _generate_reports(self) -> None:
        """Generate statistics and visualizations."""
        # Print summary to console
        self.processor.stats.print_summary()

        # Write summary to file
        summary_file = Path(f'{self.args.output}_summary_statistics.txt')
        FileHandler.write_summary_stats(self.processor.stats, summary_file)

        # Create plots if requested
        if not self.args.no_plot:
            logger.info("Generating visualization plots...")
            Visualizer.plot_summary(self.processor.stats, self.args.output)


# ============================================================================
# COMMAND-LINE INTERFACE
# ============================================================================

def parse_arguments() -> argparse.Namespace:
    """Parse and validate command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Memory-efficient BAM file read alignment analyzer (requires name-sorted BAM)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    # Required arguments
    parser.add_argument(
        'bam_file',
        help='Input BAM file path (must be sorted by query name)'
    )

    # Optional arguments
    parser.add_argument(
        '-o', '--output',
        default=DEFAULT_OUTPUT_PREFIX,
        help=f'Output prefix for all generated files (default: {DEFAULT_OUTPUT_PREFIX})'
    )

    parser.add_argument(
        '-m', '--min-length',
        type=int,
        default=DEFAULT_MIN_LENGTH,
        help=f'Minimum mapped length threshold in bp (default: {DEFAULT_MIN_LENGTH})'
    )

    # Feature flags
    parser.add_argument(
        '--extract-unique',
        action='store_true',
        help='Extract reads mapped to unique transcripts and save to separate BAM'
    )

    parser.add_argument(
        '--full-analysis',
        action='store_true',
        help='Generate detailed per-read TSV output (uses streaming write)'
    )

    parser.add_argument(
        '--no-plot',
        action='store_true',
        help='Skip generating visualization plots'
    )

    # Logging options
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )

    parser.add_argument(
        '-q', '--quiet',
        action='store_true',
        help='Suppress all but error messages'
    )

    args = parser.parse_args()

    # Validate arguments
    if not Path(args.bam_file).exists():
        parser.error(f"BAM file not found: {args.bam_file}")

    if args.min_length < 0:
        parser.error("Minimum length must be non-negative")

    # Adjust logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    elif args.quiet:
        logging.getLogger().setLevel(logging.ERROR)

    return args


# ============================================================================
# ENTRY POINT
# ============================================================================

def main():
    """Main entry point."""
    try:
        args = parse_arguments()
        app = BAMAnalyzerApp(args)
        app.run()

    except KeyboardInterrupt:
        logger.info("Analysis interrupted by user")
        sys.exit(1)
    except FileNotFoundError as e:
        logger.error(f"File not found: {e}")
        sys.exit(1)
    except ValueError as e:
        logger.error(str(e))
        sys.exit(1)
    except PermissionError as e:
        logger.error(f"Permission denied: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
