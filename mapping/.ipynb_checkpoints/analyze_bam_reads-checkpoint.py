#!/usr/bin/env python3
"""
BAM Read Alignment Analyzer
Xinyi Lin, 202509

Analyze BAM file to extract and visualize read alignment information:
1. Count alignment records per read ID
2. Track transcripts each read aligns to
3. Calculate maximum occurrences of the same transcript per read
4. Generate visualization of alignment patterns

Example usage:
    python analyze_bam_reads.py input.bam -o output_prefix
    python analyze_bam_reads.py input.bam -o output_prefix -m 50
"""

import argparse
import csv
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pysam


@dataclass
class ReadAlignment:
    """Data structure for storing read alignment information"""
    read_id: str
    alignment_count: int = 0
    transcripts: List[str] = field(default_factory=list)
    transcript_counts: Counter = field(default_factory=Counter)
    
    def add_alignment(self, transcript: str) -> None:
        """Add an alignment to this read"""
        self.alignment_count += 1
        self.transcripts.append(transcript)
        self.transcript_counts[transcript] += 1
    
    def get_max_transcript(self) -> Tuple[str, int]:
        """Get the transcript with maximum occurrences"""
        if self.transcript_counts:
            return self.transcript_counts.most_common(1)[0]
        return '', 0
    
    def get_unique_transcripts(self) -> List[str]:
        """Get list of unique transcripts"""
        return list(set(self.transcripts))
    
    def get_max_ratio(self) -> float:
        """Calculate max_count/alignment_count ratio"""
        if self.alignment_count > 0:
            _, max_count = self.get_max_transcript()
            return max_count / self.alignment_count
        return 0.0


class AlignmentStatistics:
    """Class for calculating and storing alignment statistics"""
    
    def __init__(self, read_data: Dict[str, ReadAlignment]):
        self.read_data = read_data
        self.total_reads = len(read_data)
        self._calculate_statistics()
    
    def _calculate_statistics(self) -> None:
        """Calculate summary statistics"""
        self.single_aligned = sum(
            1 for data in self.read_data.values() 
            if data.alignment_count == 1
        )
        self.multi_aligned = sum(
            1 for data in self.read_data.values() 
            if data.alignment_count > 1
        )
        
        self.alignment_distribution = Counter(
            data.alignment_count for data in self.read_data.values()
        )
        
        self.max_count_ratios = [
            data.get_max_ratio() for data in self.read_data.values()
        ]
    
    def print_summary(self) -> None:
        """Print summary statistics to console"""
        print("\n=== Summary Statistics ===")
        print(f"Total unique reads: {self.total_reads}")
        
        if self.total_reads > 0:
            single_pct = self.single_aligned / self.total_reads * 100
            multi_pct = self.multi_aligned / self.total_reads * 100
            print(f"Single-aligned reads: {self.single_aligned} ({single_pct:.2f}%)")
            print(f"Multi-aligned reads: {self.multi_aligned} ({multi_pct:.2f}%)")
            
            print("\nAlignment count distribution (top 10):")
            for count, freq in self.alignment_distribution.most_common(10):
                print(f"  {count} alignments: {freq} reads")
    
    def get_ratio_statistics(self) -> Dict:
        """Calculate detailed ratio statistics"""
        ratios = np.array(self.max_count_ratios)
        
        stats = {
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
        
        return stats


class BAMAnalyzer:
    """Main class for BAM file analysis"""
    
    def __init__(self, bam_file: str, min_mapped_length: int = 0):
        """
        Initialize BAM analyzer
        
        Args:
            bam_file: Path to input BAM file
            min_mapped_length: Minimum mapped length threshold for filtering
        """
        self.bam_file = Path(bam_file)
        self.min_mapped_length = min_mapped_length
        self.read_data: Dict[str, ReadAlignment] = {}
        self.total_alignments = 0
        self.filtered_alignments = 0
        
        if not self.bam_file.exists():
            raise FileNotFoundError(f"BAM file not found: {bam_file}")
    
    def process_alignments(self, progress_interval: int = 100000) -> None:
        """
        Process all alignments in the BAM file
        
        Args:
            progress_interval: Number of alignments between progress updates
        """
        print(f"Opening BAM file: {self.bam_file}")
        print(f"Minimum mapped length filter: {self.min_mapped_length} bp")
        
        with pysam.AlignmentFile(str(self.bam_file), "rb") as bamfile:
            for i, alignment in enumerate(bamfile):
                self.total_alignments += 1
                
                if i % progress_interval == 0 and i > 0:
                    print(f"Processed {i:,} alignments...")
                
                if not alignment.reference_name:
                    continue
                
                mapped_length = alignment.reference_length if alignment.reference_length else 0
                # print(f"Alignment {i}: Read {alignment.query_name}, Mapped Length: {mapped_length}")
                if mapped_length < self.min_mapped_length:
                    self.filtered_alignments += 1
                    continue
                
                read_id = alignment.query_name
                if read_id not in self.read_data:
                    self.read_data[read_id] = ReadAlignment(read_id)
                
                self.read_data[read_id].add_alignment(alignment.reference_name)
        
        self._print_processing_summary()
    
    def _print_processing_summary(self) -> None:
        """Print processing summary"""
        print(f"Total alignments in BAM: {self.total_alignments:,}")
        
        if self.total_alignments > 0:
            filtered_pct = self.filtered_alignments / self.total_alignments * 100
            passed = self.total_alignments - self.filtered_alignments
            passed_pct = passed / self.total_alignments * 100
            
            print(f"Alignments filtered (< {self.min_mapped_length} bp): "
                  f"{self.filtered_alignments:,} ({filtered_pct:.2f}%)")
            print(f"Alignments passing filter: {passed:,} ({passed_pct:.2f}%)")
        
        print(f"Total unique reads processed: {len(self.read_data):,}")
    
    def get_statistics(self) -> AlignmentStatistics:
        """Get alignment statistics object"""
        return AlignmentStatistics(self.read_data)


class ResultWriter:
    """Class for writing analysis results"""
    
    @staticmethod
    def write_tsv(read_data: Dict[str, ReadAlignment], output_file: Path) -> None:
        """
        Write results to TSV file
        
        Args:
            read_data: Dictionary of read alignment data
            output_file: Output file path
        """
        print(f"Writing results to: {output_file}")
        
        with open(output_file, 'w', newline='') as csvfile:
            fieldnames = [
                'read_id', 'alignment_count', 'all_transcripts',
                'unique_transcripts', 'max_transcript', 'max_count'
            ]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            
            for read_id, data in read_data.items():
                max_transcript, max_count = data.get_max_transcript()
                unique_transcripts = data.get_unique_transcripts()
                
                writer.writerow({
                    'read_id': read_id,
                    'alignment_count': data.alignment_count,
                    'all_transcripts': ','.join(data.transcripts),
                    'unique_transcripts': ','.join(unique_transcripts),
                    'max_transcript': max_transcript,
                    'max_count': max_count
                })


class Visualizer:
    """Class for creating visualizations"""
    
    @staticmethod
    def plot_ratio_distribution(statistics: AlignmentStatistics, output_prefix: str) -> None:
        """
        Create distribution plots for max_count/alignment_count ratios
        
        Args:
            statistics: AlignmentStatistics object
            output_prefix: Prefix for output files
        """
        ratios = np.array(statistics.max_count_ratios)
        if len(ratios) == 0:
            print("Warning: No data to plot")
            return
        
        plot_file = f'{output_prefix}_max_count_ratio_distribution.png'
        
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        # Main histogram
        axes[0, 0].hist(ratios, bins=50, edgecolor='black', alpha=0.7)
        axes[0, 0].set_xlabel('Max Count / Alignment Count Ratio')
        axes[0, 0].set_ylabel('Number of Reads')
        axes[0, 0].set_title('Distribution of Max Count to Alignment Count Ratio')
        axes[0, 0].grid(True, alpha=0.3)
        
        # Histogram for multi-aligned reads
        non_one_ratios = ratios[ratios < 1.0]
        if len(non_one_ratios) > 0:
            axes[0, 1].hist(non_one_ratios, bins=50, edgecolor='black', 
                           alpha=0.7, color='orange')
            axes[0, 1].set_xlabel('Max Count / Alignment Count Ratio (excluding 1.0)')
            axes[0, 1].set_ylabel('Number of Reads')
            axes[0, 1].set_title('Distribution for Multi-aligned Reads (Ratio < 1.0)')
            axes[0, 1].grid(True, alpha=0.3)
        else:
            axes[0, 1].text(0.5, 0.5, 'No multi-aligned reads', 
                           ha='center', va='center', transform=axes[0, 1].transAxes)
        
        # Cumulative distribution
        sorted_ratios = np.sort(ratios)
        cumulative_pct = np.arange(1, len(sorted_ratios) + 1) / len(sorted_ratios) * 100
        axes[1, 0].plot(sorted_ratios, cumulative_pct, 'b-', linewidth=2)
        axes[1, 0].set_xlabel('Max Count / Alignment Count Ratio')
        axes[1, 0].set_ylabel('Cumulative Percentage of Reads (%)')
        axes[1, 0].set_title('Cumulative Distribution of Ratio')
        axes[1, 0].grid(True, alpha=0.3)
        axes[1, 0].axhline(y=50, color='r', linestyle='--', alpha=0.5, label='50th percentile')
        axes[1, 0].axhline(y=90, color='g', linestyle='--', alpha=0.5, label='90th percentile')
        axes[1, 0].axvline(x=1.0, color='purple', linestyle=':', alpha=0.5, label='Ratio = 1.0')
        axes[1, 0].legend()
        
        # Box plot for ratio categories
        Visualizer._create_category_boxplot(axes[1, 1], ratios)
        
        plt.tight_layout()
        plt.savefig(plot_file, dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"\nMax count ratio distribution plot saved to: {plot_file}")
        
        # Print ratio statistics
        Visualizer._print_ratio_statistics(statistics)
    
    @staticmethod
    def _create_category_boxplot(ax, ratios: np.ndarray) -> None:
        """Create boxplot for ratio categories"""
        categories = []
        labels = []
        
        single = ratios[ratios == 1.0]
        high = ratios[(ratios >= 0.8) & (ratios < 1.0)]
        moderate = ratios[(ratios >= 0.5) & (ratios < 0.8)]
        low = ratios[(ratios > 0) & (ratios < 0.5)]
        
        for data, label in [
            (single, f'Single\n(n={len(single)})'),
            (high, f'High\n[0.8-1.0)\n(n={len(high)})'),
            (moderate, f'Moderate\n[0.5-0.8)\n(n={len(moderate)})'),
            (low, f'Low\n[0-0.5)\n(n={len(low)})')
        ]:
            if len(data) > 0:
                categories.append(data)
                labels.append(label)
        
        if categories:
            ax.boxplot(categories, labels=labels)
            ax.set_ylabel('Max Count / Alignment Count Ratio')
            ax.set_title('Ratio Distribution by Dominance Category')
            ax.grid(True, alpha=0.3, axis='y')
        else:
            ax.text(0.5, 0.5, 'No data to display', 
                   ha='center', va='center', transform=ax.transAxes)
    
    @staticmethod
    def _print_ratio_statistics(statistics: AlignmentStatistics) -> None:
        """Print detailed ratio statistics"""
        stats = statistics.get_ratio_statistics()
        
        print("\n=== Max Count / Alignment Count Ratio Statistics ===")
        print(f"Total reads analyzed: {stats['total']}")
        print(f"Mean ratio: {stats['mean']:.4f}")
        print(f"Median ratio: {stats['median']:.4f}")
        print(f"Standard deviation: {stats['std']:.4f}")
        print(f"Min ratio: {stats['min']:.4f}")
        print(f"Max ratio: {stats['max']:.4f}")
        
        print("\nRatio Categories:")
        total = stats['total']
        if total > 0:
            for category, count in [
                ("Reads with ratio = 1.0 (single-aligned or all to same transcript)", 
                 stats['single_aligned']),
                ("Reads with ratio [0.8-1.0) (highly dominant transcript)", 
                 stats['highly_dominant']),
                ("Reads with ratio [0.5-0.8) (moderately dominant transcript)", 
                 stats['moderately_dominant']),
                ("Reads with ratio (0-0.5) (low dominant transcript)", 
                 stats['low_dominant'])
            ]:
                pct = count / total * 100
                print(f"  {category}: {count} ({pct:.2f}%)")
        
        print("\nPercentiles:")
        for p, value in stats['percentiles'].items():
            print(f"  {p}th percentile: {value:.4f}")


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Analyze BAM file read alignments',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument(
        'bam_file',
        help='Input BAM file path'
    )
    
    parser.add_argument(
        '-o', '--output',
        default='read_alignment_analysis',
        help='Output prefix for files (default: read_alignment_analysis)'
    )
    
    parser.add_argument(
        '-m', '--min-length',
        type=int,
        default=0,
        help='Minimum mapped length threshold in bp (default: 0, no filtering)'
    )
    
    parser.add_argument(
        '--no-plot',
        action='store_true',
        help='Skip generating visualization plots'
    )
    
    return parser.parse_args()


def main():
    """Main execution function"""
    try:
        args = parse_arguments()
        
        # Initialize analyzer
        analyzer = BAMAnalyzer(args.bam_file, args.min_length)
        
        # Process alignments
        analyzer.process_alignments()
        
        # Write results
        output_file = Path(f'{args.output}_read_alignment_analysis.tsv')
        ResultWriter.write_tsv(analyzer.read_data, output_file)
        
        # Get and print statistics
        statistics = analyzer.get_statistics()
        statistics.print_summary()
        
        # Create visualizations
        if not args.no_plot:
            Visualizer.plot_ratio_distribution(statistics, args.output)
        
        print("Analysis complete!")
        
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except PermissionError as e:
        print(f"Permission denied: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()