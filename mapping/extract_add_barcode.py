#!/usr/bin/env python3

# Example Usage:
# python extract_add_barcode.py add_tags input.bam output_tagged.bam
# python extract_add_barcode.py add_tags input.bam output.bam --barcode-tag BC --umi-tag UM
# python extract_add_barcode.py extract_barcodes merge2_DMSO_tagged.bam barcodes.txt
# python extract_add_barcode.py extract_barcodes merge2_DMSO_tagged.bam unique_barcodes.txt --unique

import pysam
import argparse
import re
from collections import defaultdict

def parse_read_name(read_name):
    """
    Parse read name to extract barcode and UMI.
    Expected format: read_name_BARCODE_UMI
    
    Example: LH00504:260:232GTNLT3:1:2284:14488:25291_TCAGAGAAATGAAAGAGTACGAGCGTGG_GATAGC
    Returns: (TCAGAGAAATGAAAGAGTACGAGCGTGG, GATAGC)
    """
    # Split by underscore and take the last two parts as UMI and barcode
    parts = read_name.split('_')
    
    if len(parts) >= 3:
        # Last part is UMI, second to last is barcode
        umi = parts[-1]
        barcode = parts[-2]
        return barcode, umi
    else:
        # If format doesn't match, return None values
        return None, None

def add_barcode_umi_tags(input_bam, output_bam, barcode_tag='CB', umi_tag='UM'):
    """
    Add barcode and UMI tags to BAM file based on read names.
    
    Args:
        input_bam: Input BAM file path
        output_bam: Output BAM file path
        barcode_tag: BAM tag for barcode (default: CB for Cell Barcode)
        umi_tag: BAM tag for UMI (default: UM for UMI Barcode)
    """
    
    stats = {
        'total_reads': 0,
        'reads_with_barcode': 0,
        'reads_without_barcode': 0,
        'unique_barcodes': set(),
        'unique_umis': set()
    }
    
    print(f"Processing BAM file: {input_bam}")
    print(f"Adding tags: {barcode_tag} (barcode), {umi_tag} (UMI)")
    
    with pysam.AlignmentFile(input_bam, "rb") as inbam:
        with pysam.AlignmentFile(output_bam, "wb", template=inbam) as outbam:
            
            for read in inbam:
                stats['total_reads'] += 1
                
                # Parse barcode and UMI from read name
                barcode, umi = parse_read_name(read.query_name)
                
                if barcode and umi:
                    # Add tags to the read
                    read.set_tag(barcode_tag, barcode)
                    read.set_tag(umi_tag, umi)
                    
                    stats['reads_with_barcode'] += 1
                    stats['unique_barcodes'].add(barcode)
                    stats['unique_umis'].add(umi)
                else:
                    stats['reads_without_barcode'] += 1
                    if stats['reads_without_barcode'] <= 5:  # Show first few examples
                        print(f"Warning: Could not parse barcode/UMI from: {read.query_name}")
                
                # Write the read to output
                outbam.write(read)
                
                # Progress update
                if stats['total_reads'] % 100000 == 0:
                    print(f"Processed {stats['total_reads']:,} reads...")
    
    # Print statistics
    print("\n" + "="*50)
    print("PROCESSING SUMMARY")
    print("="*50)
    print(f"Total reads processed: {stats['total_reads']:,}")
    print(f"Reads with barcode/UMI: {stats['reads_with_barcode']:,}")
    print(f"Reads without barcode/UMI: {stats['reads_without_barcode']:,}")
    print(f"Unique barcodes found: {len(stats['unique_barcodes']):,}")
    print(f"Unique UMIs found: {len(stats['unique_umis']):,}")
    
    # Index the output BAM file
    print(f"\nIndexing output BAM file...")
    pysam.index(output_bam)
    print(f"Done! Output written to {output_bam}")

def extract_barcode_umi_list(input_bam, output_file):
    """
    Extract a list of all barcodes and UMIs from the BAM file.
    """
    barcode_umi_counts = defaultdict(lambda: defaultdict(int))
    total_reads = 0
    
    print(f"Extracting barcode/UMI combinations from {input_bam}...")
    
    with pysam.AlignmentFile(input_bam, "rb") as inbam:
        for read in inbam:
            total_reads += 1
            barcode, umi = parse_read_name(read.query_name)
            
            if barcode and umi:
                barcode_umi_counts[barcode][umi] += 1
            
            if total_reads % 100000 == 0:
                print(f"Processed {total_reads:,} reads...")
    
    # Write results to file
    with open(output_file, 'w') as f:
        f.write("barcode\tumi\tcount\n")
        for barcode in sorted(barcode_umi_counts.keys()):
            for umi in sorted(barcode_umi_counts[barcode].keys()):
                count = barcode_umi_counts[barcode][umi]
                f.write(f"{barcode}\t{umi}\t{count}\n")
    
    print(f"Barcode/UMI list written to {output_file}")
    print(f"Found {len(barcode_umi_counts)} unique barcodes")
    total_combinations = sum(len(umis) for umis in barcode_umi_counts.values())
    print(f"Found {total_combinations} unique barcode/UMI combinations")

def extract_barcodes_only(bam_file, output_file, barcode_tag='CB', unique_only=False):
    """
    Extract only barcodes from BAM file tags, one per line.
    
    Args:
        bam_file: Input BAM file with barcode tags
        output_file: Output file (one barcode per line)
        barcode_tag: BAM tag for barcode (default: CB)
        unique_only: If True, only output unique barcodes
    """
    barcodes_seen = set()
    total_reads = 0
    barcodes_extracted = 0
    reads_without_tags = 0
    
    print(f"Extracting barcodes from {bam_file}...")
    
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        with open(output_file, 'w') as outfile:
            
            for read in bam:
                total_reads += 1
                
                try:
                    barcode = read.get_tag(barcode_tag)
                    
                    # If unique_only is True, check if we've seen this barcode before
                    if unique_only:
                        if barcode not in barcodes_seen:
                            outfile.write(barcode + '\n')
                            barcodes_seen.add(barcode)
                            barcodes_extracted += 1
                    else:
                        # Write all barcodes (including duplicates)
                        outfile.write(barcode + '\n')
                        barcodes_extracted += 1
                        
                except KeyError:
                    # Read doesn't have barcode tag
                    reads_without_tags += 1
                    continue
                
                # Progress update
                if total_reads % 100000 == 0:
                    print(f"Processed {total_reads:,} reads...")
    
    print(f"\nCompleted!")
    print(f"Total reads processed: {total_reads:,}")
    print(f"Reads with barcode tags: {total_reads - reads_without_tags:,}")
    print(f"Reads without barcode tags: {reads_without_tags:,}")
    print(f"Barcodes extracted: {barcodes_extracted:,}")
    if unique_only:
        print(f"Unique barcodes: {len(barcodes_seen):,}")
    print(f"Output written to: {output_file}")

def validate_tags(bam_file, barcode_tag='CB', umi_tag='UM', num_examples=10):
    """
    Validate that tags were added correctly by showing examples.
    """
    print(f"\nValidating tags in {bam_file}...")
    print(f"Showing first {num_examples} examples:")
    print("-" * 80)
    print(f"{'Read Name':<50} {'Barcode':<30} {'UMI':<10}")
    print("-" * 80)
    
    count = 0
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            if count >= num_examples:
                break
                
            try:
                barcode = read.get_tag(barcode_tag)
                umi = read.get_tag(umi_tag)
                print(f"{read.query_name:<50} {barcode:<30} {umi:<10}")
                count += 1
            except KeyError:
                # Skip reads without tags
                continue

def main():
    parser = argparse.ArgumentParser(description='Add barcode and UMI tags to BAM file from read names')
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Add tags command
    add_parser = subparsers.add_parser('add_tags', help='Add barcode and UMI tags to BAM file')
    add_parser.add_argument('input_bam', help='Input BAM file')
    add_parser.add_argument('output_bam', help='Output BAM file with tags')
    add_parser.add_argument('--barcode-tag', default='CB', help='BAM tag for barcode (default: CB)')
    add_parser.add_argument('--umi-tag', default='UM', help='BAM tag for UMI (default: UM)')
    
    # Extract list command
    extract_parser = subparsers.add_parser('extract_list', help='Extract list of barcodes and UMIs')
    extract_parser.add_argument('input_bam', help='Input BAM file')
    extract_parser.add_argument('output_file', help='Output TSV file with barcode/UMI list')
    
    # Extract barcodes command
    extract_barcodes_parser = subparsers.add_parser('extract_barcodes', help='Extract only barcodes from tagged BAM file')
    extract_barcodes_parser.add_argument('bam_file', help='Input BAM file with barcode tags')
    extract_barcodes_parser.add_argument('output_file', help='Output file (one barcode per line)')
    extract_barcodes_parser.add_argument('--barcode-tag', default='CB', help='BAM tag for barcode (default: CB)')
    extract_barcodes_parser.add_argument('--unique', action='store_true', help='Only output unique barcodes')
    
    # Validate command
    validate_parser = subparsers.add_parser('validate', help='Validate tags in BAM file')
    validate_parser.add_argument('bam_file', help='BAM file to validate')
    validate_parser.add_argument('--barcode-tag', default='CB', help='BAM tag for barcode (default: CB)')
    validate_parser.add_argument('--umi-tag', default='UM', help='BAM tag for UMI (default: UM)')
    validate_parser.add_argument('--examples', type=int, default=10, help='Number of examples to show')
    
    args = parser.parse_args()
    
    if args.command == 'add_tags':
        add_barcode_umi_tags(args.input_bam, args.output_bam, args.barcode_tag, args.umi_tag)
    elif args.command == 'extract_list':
        extract_barcode_umi_list(args.input_bam, args.output_file)
    elif args.command == 'extract_barcodes':
        extract_barcodes_only(args.bam_file, args.output_file, args.barcode_tag, args.unique)
    elif args.command == 'validate':
        validate_tags(args.bam_file, args.barcode_tag, args.umi_tag, args.examples)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()