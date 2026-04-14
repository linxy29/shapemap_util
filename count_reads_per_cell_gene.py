#!/usr/bin/env python3
"""
Count reads per cell and gene from BAM file with cell barcodes and gene tags.
Output format: genes as rows, cells as columns (TSV format)
Usage: count_reads_per_cell_gene.py -i input.bam -o output.tsv [options]
"""

import pysam
import pandas as pd
from collections import defaultdict
from tqdm import tqdm
import argparse
import sys

def count_reads_per_cell_gene(bam_file, output_tsv, cb_tag="CB", gene_tag="XT", 
                               min_mapq=0, verbose=True):
    """
    Count reads per cell and gene from BAM file.
    
    Args:
        bam_file: Path to BAM file
        output_tsv: Output TSV file name
        cb_tag: Cell barcode tag (default: CB)
        gene_tag: Gene tag (default: XT)
        min_mapq: Minimum mapping quality (default: 0)
        verbose: Print progress information
    """
    counts = defaultdict(lambda: defaultdict(int))
    
    if verbose:
        print(f"Reading BAM file: {bam_file}")
        print(f"Cell barcode tag: {cb_tag}")
        print(f"Gene tag: {gene_tag}")
        print(f"Minimum MAPQ: {min_mapq}")
        print("-" * 50)
    
    try:
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            total_reads = 0
            counted_reads = 0
            skipped_unmapped = 0
            skipped_no_tags = 0
            skipped_low_mapq = 0
            
            # Get total for progress bar (optional, can be slow for large files)
            try:
                total = bam.count() if verbose else None
                bam.reset()  # Reset file pointer after count()
            except:
                total = None
            
            iterator = tqdm(bam, total=total, desc="Processing reads") if verbose else bam
            
            for read in iterator:
                total_reads += 1
                
                # Skip unmapped reads
                if read.is_unmapped:
                    skipped_unmapped += 1
                    continue
                
                # Skip low quality mappings
                if read.mapping_quality < min_mapq:
                    skipped_low_mapq += 1
                    continue
                
                try:
                    cell_bc = read.get_tag(cb_tag)
                    gene = read.get_tag(gene_tag)
                    counts[gene][cell_bc] += 1
                    counted_reads += 1
                except KeyError:
                    # Skip reads without CB or XT tags
                    skipped_no_tags += 1
                    continue
    
    except FileNotFoundError:
        print(f"Error: BAM file '{bam_file}' not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading BAM file: {e}", file=sys.stderr)
        sys.exit(1)
    
    if verbose:
        print(f"\nProcessing summary:")
        print(f"  Total reads: {total_reads:,}")
        print(f"  Counted reads: {counted_reads:,}")
        print(f"  Skipped (unmapped): {skipped_unmapped:,}")
        print(f"  Skipped (low MAPQ): {skipped_low_mapq:,}")
        print(f"  Skipped (no tags): {skipped_no_tags:,}")
    
    if counted_reads == 0:
        print("Warning: No reads were counted. Check your BAM file and tag names.", 
              file=sys.stderr)
        sys.exit(1)
    
    # Convert to DataFrame (genes as rows, cells as columns)
    if verbose:
        print("\nConverting to matrix...")
    
    df = pd.DataFrame(counts).fillna(0).astype(int)
    
    # Sort genes alphabetically and cells alphabetically
    df = df.sort_index()  # Sort gene names (rows)
    df = df[sorted(df.columns)]  # Sort cell barcodes (columns)
    
    if verbose:
        print(f"Matrix shape: {df.shape[0]} genes × {df.shape[1]} cells")
        print(f"Total counts in matrix: {df.sum().sum():,}")
    
    # Save to TSV with 'gene' as the index name
    if verbose:
        print(f"\nSaving to {output_tsv}...")
    
    df.index.name = 'gene'
    df.to_csv(output_tsv, sep='\t')
    
    if verbose:
        print(f"Done! Output saved to {output_tsv}")
        print("\n" + "="*50)
        print("Matrix Statistics:")
        print("="*50)
        print(f"Number of genes: {df.shape[0]}")
        print(f"Number of cells: {df.shape[1]}")
        print(f"\nTop 10 genes by total counts:")
        print(df.sum(axis=1).sort_values(ascending=False).head(10))
        print(f"\nTop 10 cells by total counts:")
        print(df.sum(axis=0).sort_values(ascending=False).head(10))
    
    return df

def main():
    parser = argparse.ArgumentParser(
        description='Count reads per cell and gene from BAM file with barcodes',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  %(prog)s -i input.bam -o output.tsv
  
  # With custom tags
  %(prog)s -i input.bam -o output.tsv --cb-tag CR --gene-tag GN
  
  # With mapping quality filter
  %(prog)s -i input.bam -o output.tsv --min-mapq 10
  
  # Quiet mode (no progress output)
  %(prog)s -i input.bam -o output.tsv --quiet
        """
    )
    
    parser.add_argument('-i', '--input', required=True,
                        help='Input BAM file with cell barcodes and gene tags')
    parser.add_argument('-o', '--output', required=True,
                        help='Output TSV file for count matrix')
    parser.add_argument('--cb-tag', default='CB',
                        help='Cell barcode tag in BAM file (default: CB)')
    parser.add_argument('--gene-tag', default='XT',
                        help='Gene tag in BAM file (default: XT)')
    parser.add_argument('--min-mapq', type=int, default=0,
                        help='Minimum mapping quality (default: 0)')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='Suppress progress output')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    
    args = parser.parse_args()
    
    # Run the counting
    count_reads_per_cell_gene(
        bam_file=args.input,
        output_tsv=args.output,
        cb_tag=args.cb_tag,
        gene_tag=args.gene_tag,
        min_mapq=args.min_mapq,
        verbose=not args.quiet
    )

if __name__ == "__main__":
    main()