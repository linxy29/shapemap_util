#!/usr/bin/env python3
"""
BAM Read Alignment Analyzer
Xinyi Lin, 202509

A comprehensive tool for analyzing BAM file read alignments with focus on:
- Multi-mapping read analysis
- Transcript assignment statistics
- Read filtering by alignment quality
- Extraction of uniquely mapped reads (v2)
- Extract the primary alignments (updated in v2.1)

Version: v2.1

Usage:
    python analyze_bam_reads_v2.py input.bam -o output_prefix
    python analyze_bam_reads_v2.py input.bam -o output_prefix -m 50 --extract-unique
"""

# ============================================================================
# IMPORTS
# ============================================================================

import argparse
import csv
import logging
import sys
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set

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
PLOT_FIGSIZE = (14, 10)
HISTOGRAM_BINS = 50

# Ratio categories for analysis
RATIO_THRESHOLDS = {
    'single': (1.0, 1.0),
    'highly_dominant': (0.8, 1.0),
    'moderately_dominant': (0.5, 0.8),
    'low_dominant': (0.0, 0.5)
}

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

    Attributes:
        read_id: Unique identifier for the read
        alignment_count: Total number of alignments for this read
        transcripts: List of all transcript alignments (may contain duplicates)
        transcript_counts: Counter of transcript occurrences
    """
    read_id: str
    alignment_count: int = 0
    transcripts: List[str] = field(default_factory=list)
    transcript_counts: Counter = field(default_factory=Counter)

    def add_alignment(self, transcript: str) -> None:
        """Add an alignment record for this read."""
        self.alignment_count += 1
        self.transcripts.append(transcript)
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

# ============================================================================
# STATISTICS MODULE
# ============================================================================

class AlignmentStatistics:
    """
    Calculates and stores comprehensive alignment statistics.
    """

    def __init__(self, read_data: Dict[str, ReadAlignment]):
        self.read_data = read_data
        self.total_reads = len(read_data)
        self._statistics_calculated = False
        self._calculate_statistics()

    def _calculate_statistics(self) -> None:
        """Calculate all summary statistics."""
        if self.total_reads == 0:
            self._set_empty_statistics()
            return

        # Basic counts
        self.single_aligned = sum(
            1 for data in self.read_data.values()
            if data.alignment_count == 1
        )
        self.multi_aligned = sum(
            1 for data in self.read_data.values()
            if data.alignment_count > 1
        )

        # Unique mapping counts
        self.uniquely_mapped = sum(
            1 for data in self.read_data.values()
            if data.is_uniquely_mapped()
        )

        # Alignment distribution
        self.alignment_distribution = Counter(
            data.alignment_count for data in self.read_data.values()
        )

        # Max count ratios
        self.max_count_ratios = [
            data.get_max_ratio() for data in self.read_data.values()
        ]

        self._statistics_calculated = True

    def _set_empty_statistics(self) -> None:
        """Set default values for empty dataset."""
        self.single_aligned = 0
        self.multi_aligned = 0
        self.uniquely_mapped = 0
        self.alignment_distribution = Counter()
        self.max_count_ratios = []
        self._statistics_calculated = True

    def get_summary_dict(self) -> Dict:
        """Get summary statistics as a dictionary."""
        if not self._statistics_calculated:
            self._calculate_statistics()

        return {
            'total_reads': self.total_reads,
            'single_aligned': self.single_aligned,
            'multi_aligned': self.multi_aligned,
            'uniquely_mapped': self.uniquely_mapped,
            'single_aligned_pct': self._safe_percentage(self.single_aligned, self.total_reads),
            'multi_aligned_pct': self._safe_percentage(self.multi_aligned, self.total_reads),
            'uniquely_mapped_pct': self._safe_percentage(self.uniquely_mapped, self.total_reads)
        }

    def get_ratio_statistics(self) -> Dict:
        """Calculate detailed ratio statistics."""
        ratios = np.array(self.max_count_ratios)

        if len(ratios) == 0:
            return self._empty_ratio_statistics()

        return {
            'total': len(ratios),
            'mean': np.mean(ratios),
            'median': np.median(ratios),
            'std': np.std(ratios),
            'min': np.min(ratios),
            'max': np.max(ratios),
            'single_aligned': np.sum(ratios == 1.0),
            'highly_dominant': np.sum((ratios >= 0.8) & (ratios < 1.0)),
            'moderately_dominant': np.sum((ratios >= 0.5) & (ratios < 0.8)),
            'low_dominant': np.sum((ratios > 0) & (ratios < 0.5)),
            'percentiles': {p: np.percentile(ratios, p) for p in [10, 25, 50, 75, 90, 95, 99]}
        }

    def _empty_ratio_statistics(self) -> Dict:
        """Return empty ratio statistics structure."""
        return {
            'total': 0,
            'mean': 0,
            'median': 0,
            'std': 0,
            'min': 0,
            'max': 0,
            'single_aligned': 0,
            'highly_dominant': 0,
            'moderately_dominant': 0,
            'low_dominant': 0,
            'percentiles': {p: 0 for p in [10, 25, 50, 75, 90, 95, 99]}
        }

    @staticmethod
    def _safe_percentage(count: int, total: int) -> float:
        """Calculate percentage safely, handling division by zero."""
        return (count / total * 100) if total > 0 else 0.0

    def print_summary(self) -> None:
        """Print formatted summary statistics to console."""
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

            if self.alignment_distribution:
                print("\nTop 10 alignment count distribution:")
                for count, freq in self.alignment_distribution.most_common(10):
                    print(f"  {count:3d} alignments: {freq:,} reads")

        print("=" * 60)

# ============================================================================
# BAM PROCESSING MODULE
# ============================================================================

class BAMProcessor:
    """
    Handles all BAM file processing operations.
    """

    def __init__(self, bam_file: str, min_mapped_length: int = 0):
        """
        Initialize BAM processor.

        Args:
            bam_file: Path to input BAM file
            min_mapped_length: Minimum mapped length threshold

        Raises:
            FileNotFoundError: If BAM file doesn't exist
        """
        self.bam_file = Path(bam_file)
        self.min_mapped_length = min_mapped_length
        self.read_data: Dict[str, ReadAlignment] = {}
        self.processing_stats = {
            'total_alignments': 0,
            'filtered_alignments': 0,
            'unmapped_alignments': 0
        }

        if not self.bam_file.exists():
            raise FileNotFoundError(f"BAM file not found: {bam_file}")

        logger.info(f"Initialized BAM processor for: {self.bam_file}")

    def process_alignments(self, progress_interval: int = PROGRESS_INTERVAL) -> None:
        """
        Process all alignments in the BAM file.

        Args:
            progress_interval: Number of alignments between progress updates
        """
        logger.info(f"Processing BAM file with min length filter: {self.min_mapped_length} bp")

        try:
            with pysam.AlignmentFile(str(self.bam_file), "rb") as bamfile:
                self._process_bam_records(bamfile, progress_interval)
        except Exception as e:
            logger.error(f"Error processing BAM file: {e}")
            raise

        self._log_processing_summary()

    def _process_bam_records(self, bamfile, progress_interval: int) -> None:
        """Process individual BAM records."""
        for i, alignment in enumerate(bamfile):
            self.processing_stats['total_alignments'] += 1

            # Progress reporting
            if i % progress_interval == 0 and i > 0:
                logger.info(f"Processed {i:,} alignments...")

            # Skip unmapped reads
            if not alignment.reference_name:
                self.processing_stats['unmapped_alignments'] += 1
                continue

            # Apply length filter
            mapped_length = self._get_mapped_length(alignment)
            if mapped_length < self.min_mapped_length:
                self.processing_stats['filtered_alignments'] += 1
                continue

            # Store alignment data
            self._store_alignment(alignment)

    def _get_mapped_length(self, alignment) -> int:
        """Get the mapped length of an alignment safely."""
        return alignment.reference_length if alignment.reference_length is not None else 0

    def _store_alignment(self, alignment) -> None:
        """Store alignment information."""
        read_id = alignment.query_name
        if read_id not in self.read_data:
            self.read_data[read_id] = ReadAlignment(read_id)

        self.read_data[read_id].add_alignment(alignment.reference_name)

    def _log_processing_summary(self) -> None:
        """Log processing summary statistics."""
        stats = self.processing_stats
        total = stats['total_alignments']

        logger.info(f"Processing complete:")
        logger.info(f"  Total alignments: {total:,}")

        if total > 0:
            filtered_pct = stats['filtered_alignments'] / total * 100
            unmapped_pct = stats['unmapped_alignments'] / total * 100
            passed = total - stats['filtered_alignments'] - stats['unmapped_alignments']
            passed_pct = passed / total * 100

            logger.info(f"  Unmapped: {stats['unmapped_alignments']:,} ({unmapped_pct:.2f}%)")
            logger.info(f"  Filtered (< {self.min_mapped_length} bp): {stats['filtered_alignments']:,} ({filtered_pct:.2f}%)")
            logger.info(f"  Passed: {passed:,} ({passed_pct:.2f}%)")

        logger.info(f"  Unique reads processed: {len(self.read_data):,}")

    def extract_unique_transcript_reads(self) -> Dict[str, ReadAlignment]:
        """
        Extract reads that map to only one unique transcript.

        Returns:
            Dictionary of read_id to ReadAlignment for uniquely mapped reads
        """
        unique_reads = {
            read_id: alignment
            for read_id, alignment in self.read_data.items()
            if alignment.is_uniquely_mapped()
        }

        logger.info(f"Extracted {len(unique_reads):,} uniquely mapped reads")
        return unique_reads

    def get_statistics(self) -> AlignmentStatistics:
        """Get alignment statistics object."""
        return AlignmentStatistics(self.read_data)

# ============================================================================
# FILE I/O MODULE
# ============================================================================

class FileHandler:
    """
    Handles all file I/O operations.
    """

    @staticmethod
    def write_analysis_tsv(read_data: Dict[str, ReadAlignment], output_file: Path) -> None:
        """
        Write detailed analysis results to TSV file.

        Args:
            read_data: Dictionary of read alignment data
            output_file: Output file path
        """
        logger.info(f"Writing analysis results to: {output_file}")

        try:
            with open(output_file, 'w', newline='') as csvfile:
                fieldnames = [
                    'read_id',
                    'alignment_count',
                    'unique_transcript_count',
                    'is_uniquely_mapped',
                    'max_transcript',
                    'max_count',
                    'max_ratio',
                    'all_transcripts',
                    'unique_transcripts'
                ]
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
                writer.writeheader()

                for read_id, data in read_data.items():
                    max_transcript, max_count = data.get_max_transcript()
                    unique_transcripts = data.get_unique_transcripts()

                    writer.writerow({
                        'read_id': read_id,
                        'alignment_count': data.alignment_count,
                        'unique_transcript_count': len(unique_transcripts),
                        'is_uniquely_mapped': 'Yes' if data.is_uniquely_mapped() else 'No',
                        'max_transcript': max_transcript,
                        'max_count': max_count,
                        'max_ratio': f"{data.get_max_ratio():.4f}",
                        'all_transcripts': ','.join(data.transcripts),
                        'unique_transcripts': ','.join(unique_transcripts)
                    })

            logger.info(f"Successfully wrote {len(read_data)} records")

        except Exception as e:
            logger.error(f"Error writing TSV file: {e}")
            raise

    @staticmethod
    def write_unique_reads_bam(
        unique_read_ids: Set[str],
        input_bam: Path,
        output_bam: Path,
        min_mapped_length: int = 0
    ) -> None:
        """
        Write reads with unique transcript mappings to BAM file.
        Only primary alignments are extracted (v2.1: skips secondary and supplementary alignments).

        Args:
            unique_read_ids: Set of read IDs that map to single unique transcript
            input_bam: Input BAM file path
            output_bam: Output BAM file path
            min_mapped_length: Minimum mapped length threshold
        """
        logger.info(f"Writing unique reads (primary alignments only) to BAM: {output_bam}")

        try:
            with pysam.AlignmentFile(str(input_bam), "rb") as infile:
                with pysam.AlignmentFile(str(output_bam), "wb", template=infile) as outfile:
                    written_count = 0
                    total_checked = 0
                    skipped_secondary = 0
                    skipped_supplementary = 0

                    for alignment in infile:
                        total_checked += 1

                        # Apply same filtering as in processing
                        if not alignment.reference_name:
                            continue

                        mapped_length = alignment.reference_length if alignment.reference_length is not None else 0
                        if mapped_length < min_mapped_length:
                            continue

                        # Write if in unique set AND is primary alignment
                        if alignment.query_name in unique_read_ids:
                            # Skip secondary alignments (flag 0x100)
                            if alignment.is_secondary:
                                skipped_secondary += 1
                                continue

                            # Skip supplementary alignments (flag 0x800)
                            if alignment.is_supplementary:
                                skipped_supplementary += 1
                                continue

                            # Write primary alignment
                            outfile.write(alignment)
                            written_count += 1

                        # Progress reporting
                        if total_checked % PROGRESS_INTERVAL == 0:
                            logger.info(f"Checked {total_checked:,} alignments, written {written_count:,}")

            logger.info(f"Wrote {written_count:,} primary alignment records from {len(unique_read_ids):,} unique reads")
            logger.info(f"Skipped {skipped_secondary:,} secondary and {skipped_supplementary:,} supplementary alignments")

        except Exception as e:
            logger.error(f"Error writing BAM file: {e}")
            raise

# ============================================================================
# VISUALIZATION MODULE
# ============================================================================

class Visualizer:
    """
    Creates all visualization plots.
    """

    @staticmethod
    def plot_ratio_distribution(statistics: AlignmentStatistics, output_prefix: str) -> None:
        """
        Create comprehensive ratio distribution plots.

        Args:
            statistics: AlignmentStatistics object
            output_prefix: Prefix for output files
        """
        ratios = np.array(statistics.max_count_ratios)

        if len(ratios) == 0:
            logger.warning("No data available for plotting")
            return

        plot_file = f'{output_prefix}_ratio_distribution.png'

        try:
            Visualizer._create_ratio_plots(ratios, plot_file)
            logger.info(f"Saved ratio distribution plot to: {plot_file}")

            # Print detailed statistics
            Visualizer._print_ratio_statistics(statistics)

        except Exception as e:
            logger.error(f"Error creating plots: {e}")

    @staticmethod
    def _create_ratio_plots(ratios: np.ndarray, output_file: str) -> None:
        """Create the actual ratio distribution plots."""
        _, axes = plt.subplots(2, 2, figsize=PLOT_FIGSIZE)

        # Main histogram
        Visualizer._plot_main_histogram(axes[0, 0], ratios)

        # Multi-aligned reads histogram
        Visualizer._plot_multi_aligned_histogram(axes[0, 1], ratios)

        # Cumulative distribution
        Visualizer._plot_cumulative_distribution(axes[1, 0], ratios)

        # Category boxplot
        Visualizer._plot_category_boxplot(axes[1, 1], ratios)

        plt.tight_layout()
        plt.savefig(output_file, dpi=PLOT_DPI, bbox_inches='tight')
        plt.close()

    @staticmethod
    def _plot_main_histogram(ax, ratios: np.ndarray) -> None:
        """Plot the main ratio histogram."""
        ax.hist(ratios, bins=HISTOGRAM_BINS, edgecolor='black', alpha=0.7, color='steelblue')
        ax.set_xlabel('Max Count / Alignment Count Ratio')
        ax.set_ylabel('Number of Reads')
        ax.set_title('Distribution of Transcript Dominance Ratio')
        ax.grid(True, alpha=0.3)
        ax.axvline(x=1.0, color='red', linestyle='--', alpha=0.5, label='Uniquely mapped')
        ax.legend()

    @staticmethod
    def _plot_multi_aligned_histogram(ax, ratios: np.ndarray) -> None:
        """Plot histogram for multi-aligned reads only."""
        non_unique = ratios[ratios < 1.0]

        if len(non_unique) > 0:
            ax.hist(non_unique, bins=HISTOGRAM_BINS, edgecolor='black',
                   alpha=0.7, color='darkorange')
            ax.set_xlabel('Ratio (Multi-aligned Reads Only)')
            ax.set_ylabel('Number of Reads')
            ax.set_title(f'Multi-aligned Reads Distribution (n={len(non_unique):,})')
            ax.grid(True, alpha=0.3)
        else:
            ax.text(0.5, 0.5, 'No multi-aligned reads',
                   ha='center', va='center', transform=ax.transAxes, fontsize=12)
            ax.set_title('Multi-aligned Reads Distribution')

    @staticmethod
    def _plot_cumulative_distribution(ax, ratios: np.ndarray) -> None:
        """Plot cumulative distribution."""
        sorted_ratios = np.sort(ratios)
        cumulative_pct = np.arange(1, len(sorted_ratios) + 1) / len(sorted_ratios) * 100

        ax.plot(sorted_ratios, cumulative_pct, 'b-', linewidth=2)
        ax.set_xlabel('Max Count / Alignment Count Ratio')
        ax.set_ylabel('Cumulative Percentage (%)')
        ax.set_title('Cumulative Distribution')
        ax.grid(True, alpha=0.3)

        # Add reference lines
        ax.axhline(y=50, color='red', linestyle='--', alpha=0.5, label='50th percentile')
        ax.axhline(y=90, color='green', linestyle='--', alpha=0.5, label='90th percentile')
        ax.axvline(x=1.0, color='purple', linestyle=':', alpha=0.5, label='Uniquely mapped')
        ax.legend()

    @staticmethod
    def _plot_category_boxplot(ax, ratios: np.ndarray) -> None:
        """Create boxplot for ratio categories."""
        categories = []
        labels = []

        # Define categories
        category_data = [
            (ratios[ratios == 1.0], 'Unique\nMapped'),
            (ratios[(ratios >= 0.8) & (ratios < 1.0)], 'High\nDominance\n[0.8-1.0)'),
            (ratios[(ratios >= 0.5) & (ratios < 0.8)], 'Moderate\nDominance\n[0.5-0.8)'),
            (ratios[(ratios > 0) & (ratios < 0.5)], 'Low\nDominance\n(0-0.5)')
        ]

        for data, label in category_data:
            if len(data) > 0:
                categories.append(data)
                labels.append(f'{label}\n(n={len(data):,})')

        if categories:
            bp = ax.boxplot(categories, labels=labels, patch_artist=True)

            # Color the boxes
            colors = ['lightgreen', 'lightyellow', 'lightcoral', 'lightgray']
            for patch, color in zip(bp['boxes'], colors[:len(bp['boxes'])]):
                patch.set_facecolor(color)

            ax.set_ylabel('Ratio Value')
            ax.set_title('Distribution by Dominance Category')
            ax.grid(True, alpha=0.3, axis='y')
        else:
            ax.text(0.5, 0.5, 'No data to display',
                   ha='center', va='center', transform=ax.transAxes)

    @staticmethod
    def _print_ratio_statistics(statistics: AlignmentStatistics) -> None:
        """Print detailed ratio statistics."""
        stats = statistics.get_ratio_statistics()

        print("\n" + "=" * 60)
        print("TRANSCRIPT DOMINANCE RATIO STATISTICS")
        print("=" * 60)

        print(f"\nBasic Statistics:")
        print(f"  Total reads:     {stats['total']:,}")
        print(f"  Mean ratio:      {stats['mean']:.4f}")
        print(f"  Median ratio:    {stats['median']:.4f}")
        print(f"  Std deviation:   {stats['std']:.4f}")
        print(f"  Min ratio:       {stats['min']:.4f}")
        print(f"  Max ratio:       {stats['max']:.4f}")

        print(f"\nDominance Categories:")
        total = stats['total']
        if total > 0:
            categories = [
                ("Uniquely mapped (ratio = 1.0)", stats['single_aligned']),
                ("High dominance [0.8-1.0)", stats['highly_dominant']),
                ("Moderate dominance [0.5-0.8)", stats['moderately_dominant']),
                ("Low dominance (0-0.5)", stats['low_dominant'])
            ]

            for category_name, count in categories:
                pct = count / total * 100
                print(f"  {category_name:30s}: {count:,} ({pct:.2f}%)")

        print(f"\nKey Percentiles:")
        for p in [10, 25, 50, 75, 90, 95, 99]:
            value = stats['percentiles'][p]
            print(f"  {p:2d}th percentile: {value:.4f}")

        print("=" * 60)

# ============================================================================
# MAIN APPLICATION
# ============================================================================

class BAMAnalyzerApp:
    """
    Main application controller.
    """

    def __init__(self, args: argparse.Namespace):
        """Initialize the application with command-line arguments."""
        self.args = args
        self.processor = None
        self.statistics = None

    def run(self) -> None:
        """Run the complete analysis pipeline."""
        try:
            # Initialize processor
            self.processor = BAMProcessor(
                self.args.bam_file,
                self.args.min_length
            )

            # Process alignments
            logger.info("Starting alignment processing...")
            self.processor.process_alignments()

            # Write main analysis results
            self._write_main_results()

            # Extract unique reads if requested
            if self.args.extract_unique:
                self._extract_unique_reads()

            # Generate statistics and visualizations
            self._generate_reports()

            logger.info("Analysis complete!")

        except Exception as e:
            logger.error(f"Analysis failed: {e}")
            raise

    def _write_main_results(self) -> None:
        """Write the main analysis results."""
        output_file = Path(f'{self.args.output}_read_alignment_analysis.tsv')
        FileHandler.write_analysis_tsv(
            self.processor.read_data,
            output_file
        )

    def _extract_unique_reads(self) -> None:
        """Extract and save uniquely mapped reads."""
        logger.info("Extracting uniquely mapped reads...")

        unique_reads = self.processor.extract_unique_transcript_reads()
        unique_read_ids = set(unique_reads.keys())

        # Write to BAM
        output_bam = Path(f'{self.args.output}.bam')
        FileHandler.write_unique_reads_bam(
            unique_read_ids,
            self.processor.bam_file,
            output_bam,
            self.processor.min_mapped_length
        )

        # Write summary TSV
        summary_file = Path(f'{self.args.output}_unique_reads_summary.tsv')
        FileHandler.write_analysis_tsv(unique_reads, summary_file)

    def _generate_reports(self) -> None:
        """Generate statistics and visualizations."""
        self.statistics = self.processor.get_statistics()

        # Print summary
        self.statistics.print_summary()

        # Create plots if requested
        if not self.args.no_plot:
            logger.info("Generating visualization plots...")
            Visualizer.plot_ratio_distribution(
                self.statistics,
                self.args.output
            )

# ============================================================================
# COMMAND-LINE INTERFACE
# ============================================================================

def parse_arguments() -> argparse.Namespace:
    """Parse and validate command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Comprehensive BAM file read alignment analyzer',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    # Required arguments
    parser.add_argument(
        'bam_file',
        help='Input BAM file path'
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
    """Main entry point for the application."""
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
    except PermissionError as e:
        logger.error(f"Permission denied: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()