#!/home/users/astar/gis/gaog/.conda/envs/main/bin/python

## This script is getting from Gao Ge.
## Version: v5.0
## Parallel processing using tabix random access.
## Each worker fetches and processes genes independently.
## Key optimizations: parallel I/O, vectorized cell calculations, no repeated scans.

import argparse
import pysam
import gzip
import numpy as np
from multiprocessing import Pool, cpu_count
from functools import partial

__doc__="""
Calculate sliding window mutation rates from cellSNP VCF files (Parallel Version).

Description:
    This script reads a cellSNP VCF file using pysam with tabix index,
    processes genes in parallel using random access, and outputs a combined
    VCF file with all cells as sample columns.

    Key features:
    - Parallel gene processing using tabix random access (no repeated scans)
    - Vectorized calculations across all cells simultaneously
    - Memory efficient: only one gene per worker in memory

    To index vcf file: tabix -p vcf cellSNP.cells.vcf.gz

Input:
    cellSNP VCF file (BGZF compressed with .tbi index) containing per-cell mutation data.
    Expected FORMAT fields: GT:AD:DP:OTH:PL:ALL

Output:
    Gzipped VCF file with one row per window, all cells as sample columns.
    FORMAT: DP:AD:MR (coverage, mutant, mutrate)

Usage:
    python cellSNP2window_v3.py -i cellSNP.cells.vcf.gz -o output.vcf.gz -w 10 -s 10 --coverage-threshold 200 --mutrate-threshold 0.2 -p 8

Arguments:
    -i, --input              Input cellSNP VCF file (BGZF compressed with .tbi index)
    -o, --output             Output VCF file (gzipped)
    -w, --window_length      Window length in bp (default: 10)
    -s, --step               Step size in bp (default: 10)
    --coverage-threshold     Minimum coverage sum per window (default: 200)
    --mutrate-threshold      Maximum mutrate per position to include (default: 0.2)
    -p, --processes          Number of parallel processes (default: number of CPUs)
"""
__version__="v5.0"
__author__="linxy"
__last_modify__="06-Jan-2026"


def get_vcf_info(vcf_path):
    """
    Get VCF metadata (samples, contigs) without loading all data.

    Parameters:
        vcf_path: str
            Path to the BGZF compressed VCF file

    Returns:
        samples: list
            List of sample names (cell barcodes)
        contigs: list
            List of contig names from VCF header
    """
    try:
        vcf = pysam.VariantFile(vcf_path)
    except NotImplementedError as e:
        if "seek not implemented" in str(e):
            raise SystemExit(
                f"\nError: Input file '{vcf_path}' is compressed with regular gzip, not BGZF.\n"
                f"pysam requires BGZF compression for random access.\n\n"
                f"Please convert your file using one of these commands:\n"
                f"  zcat {vcf_path} | bgzip -c > output.vcf.bgz.gz\n"
                f"  gunzip -c {vcf_path} | bgzip > output.vcf.bgz.gz\n"
            )
        raise

    samples = list(vcf.header.samples)
    contigs = list(vcf.header.contigs)
    vcf.close()
    return samples, contigs


def process_gene(gene, vcf_path, samples, n_samples, win_len, step,
                 coverage_threshold, mutrate_threshold):
    """
    Process a single gene: fetch data via tabix, calculate windows.

    Parameters:
        gene: str
            Gene/contig name to process
        vcf_path: str
            Path to the VCF file
        samples: list
            List of sample names
        n_samples: int
            Number of samples
        win_len: int
            Window length in bp
        step: int
            Step size in bp
        coverage_threshold: int
            Minimum coverage per window
        mutrate_threshold: float
            Maximum mutrate per position to include

    Returns:
        gene: str
        windows: list of tuples
            [(win_start, win_end, dp_array, ad_array, mr_array, valid_mask), ...]
    """
    vcf = pysam.VariantFile(vcf_path)

    # Collect data for this gene using tabix random access
    positions = []
    dp_list = []
    ad_list = []

    try:
        for record in vcf.fetch(contig=gene):
            positions.append(record.pos)
            dp_row = np.zeros(n_samples, dtype=np.int32)
            ad_row = np.zeros(n_samples, dtype=np.int32)

            for s_idx, sample in enumerate(samples):
                sample_data = record.samples[sample]
                dp = sample_data.get('DP', None)
                ad = sample_data.get('AD', None)

                if dp is not None and ad is not None:
                    dp_row[s_idx] = dp
                    ad_row[s_idx] = ad

            dp_list.append(dp_row)
            ad_list.append(ad_row)
    except ValueError:
        # Gene not in index
        vcf.close()
        return gene, []

    vcf.close()

    if not positions:
        return gene, []

    # Convert to numpy arrays
    positions = np.array(positions, dtype=np.int32)
    dp_matrix = np.vstack(dp_list)  # shape: (n_positions, n_samples)
    ad_matrix = np.vstack(ad_list)

    # Sort by position (should already be sorted, but ensure)
    sort_idx = np.argsort(positions)
    positions = positions[sort_idx]
    dp_matrix = dp_matrix[sort_idx]
    ad_matrix = ad_matrix[sort_idx]

    # Calculate windows
    windows = calculate_windows(
        positions, dp_matrix, ad_matrix, n_samples,
        win_len, step, coverage_threshold, mutrate_threshold
    )

    return gene, windows


def calculate_windows(positions, dp_matrix, ad_matrix, n_samples,
                      win_len, step, coverage_threshold, mutrate_threshold):
    """
    Calculate sliding windows using vectorized operations.

    Parameters:
        positions: np.array
            Sorted position array
        dp_matrix: np.array
            Coverage matrix (n_positions, n_samples)
        ad_matrix: np.array
            Alternate counts matrix (n_positions, n_samples)
        n_samples: int
            Number of samples
        win_len: int
            Window length in bp
        step: int
            Step size in bp
        coverage_threshold: int
            Minimum coverage per window
        mutrate_threshold: float
            Maximum mutrate per position to include

    Returns:
        windows: list of tuples
    """
    if len(positions) == 0:
        return []

    min_pos = int(positions[0])
    max_pos = int(positions[-1])

    windows = []
    win_start = min_pos

    while win_start <= max_pos:
        win_end = win_start + win_len - 1

        # Binary search for window boundaries
        left_idx = np.searchsorted(positions, win_start, side='left')
        right_idx = np.searchsorted(positions, win_end, side='right')

        if left_idx < right_idx:  # Window has positions
            # Slice arrays for this window (views, no copy)
            win_dp = dp_matrix[left_idx:right_idx]
            win_ad = ad_matrix[left_idx:right_idx]

            # Calculate per-position mutrate and create mask
            with np.errstate(divide='ignore', invalid='ignore'):
                win_mutrates = np.where(win_dp > 0, win_ad / win_dp, 0)
            win_mask = win_mutrates <= mutrate_threshold

            # Vectorized masked sum across all cells
            dp_sum = (win_dp * win_mask).sum(axis=0)
            ad_sum = (win_ad * win_mask).sum(axis=0)

            # Create validity mask based on coverage threshold
            valid_mask = dp_sum >= coverage_threshold

            # Only add window if at least one sample passes
            if valid_mask.any():
                # Calculate mutation rates
                with np.errstate(divide='ignore', invalid='ignore'):
                    mr = np.where(dp_sum > 0, np.round(ad_sum / dp_sum, 5), 0)

                windows.append((
                    win_start,
                    win_end,
                    dp_sum.astype(np.int32),
                    ad_sum.astype(np.int32),
                    mr.astype(np.float32),
                    valid_mask
                ))

        win_start += step

    return windows


def process_gene_wrapper(args):
    """Wrapper for multiprocessing."""
    return process_gene(*args)


def write_vcf_header(output_path, samples, contigs):
    """Write VCF header to output file."""
    with gzip.open(output_path, 'wt') as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##source=cellSNP2window_v5.0_parallel\n")
        f.write('##INFO=<ID=WIN_END,Number=1,Type=Integer,Description="Window end position">\n')
        f.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total coverage in window">\n')
        f.write('##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Total alternate counts in window">\n')
        f.write('##FORMAT=<ID=MR,Number=1,Type=Float,Description="Mutation rate (AD/DP)">\n')

        for contig in contigs:
            f.write(f'##contig=<ID={contig}>\n')

        header_cols = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
        header_cols.extend(samples)
        f.write('\t'.join(header_cols) + '\n')


def write_vcf_results(output_path, results, n_samples):
    """
    Append results to VCF file.

    Parameters:
        output_path: str
            Path to output file
        results: list
            List of (gene, windows) tuples
        n_samples: int
            Number of samples
    """
    lines = []
    for gene, windows in results:
        for win_start, win_end, dp_arr, ad_arr, mr_arr, valid_mask in windows:
            sample_cols = []
            for i in range(n_samples):
                if valid_mask[i]:
                    sample_cols.append(f"{dp_arr[i]}:{ad_arr[i]}:{mr_arr[i]}")
                else:
                    sample_cols.append('.:.:.')

            row = [
                gene,
                str(win_start),
                '.',
                'N',
                '.',
                '.',
                'PASS',
                f'WIN_END={win_end}',
                'DP:AD:MR'
            ] + sample_cols

            lines.append('\t'.join(row) + '\n')

    with gzip.open(output_path, 'at') as f:
        f.writelines(lines)


def main(args):
    """
    Main function: parallel processing of genes using tabix random access.
    """
    print(f"Reading VCF metadata: {args.input_file}")
    samples, contigs = get_vcf_info(args.input_file)
    n_samples = len(samples)
    genes = contigs  # Process all contigs
    print(f"Found {n_samples} cells and {len(genes)} genes/chromosomes")

    # Determine number of processes
    n_processes = args.processes if args.processes else cpu_count()
    n_processes = min(n_processes, len(genes))
    print(f"Using {n_processes} parallel processes")

    # Write VCF header
    print(f"Writing output to: {args.output_file}")
    write_vcf_header(args.output_file, samples, contigs)

    # Prepare arguments for parallel processing
    process_args = [
        (gene, args.input_file, samples, n_samples, args.win_len, args.step,
         args.coverage_threshold, args.mutrate_threshold)
        for gene in genes
    ]

    # Process genes in parallel
    total_windows = 0
    genes_with_windows = 0

    print("Processing genes in parallel...")
    with Pool(processes=n_processes) as pool:
        # Use imap for ordered results and progress tracking
        results_buffer = []
        buffer_size = 100  # Write every N genes

        for i, (gene, windows) in enumerate(pool.imap(process_gene_wrapper, process_args)):
            if windows:
                results_buffer.append((gene, windows))
                total_windows += len(windows)
                genes_with_windows += 1

            # Periodic write to reduce memory
            if len(results_buffer) >= buffer_size:
                write_vcf_results(args.output_file, results_buffer, n_samples)
                results_buffer = []

            # Progress update
            if (i + 1) % 500 == 0:
                print(f"  Processed {i + 1}/{len(genes)} genes...")

        # Write remaining results
        if results_buffer:
            write_vcf_results(args.output_file, results_buffer, n_samples)

    print(f"Done! Total: {total_windows} windows across {genes_with_windows} genes")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", dest='input_file', required=True,
                        help="Input cellSNP VCF file (BGZF compressed with .tbi index)")
    parser.add_argument("-o", "--output", dest='output_file', required=True,
                        help="Output VCF file (gzipped)")
    parser.add_argument("-w", "--window_length", dest='win_len', type=int, default=10,
                        help="Length of window in bp (default: 10)")
    parser.add_argument("-s", "--step", dest='step', default=10, type=int,
                        help="Step size in bp (default: 10)")
    parser.add_argument("--coverage-threshold", dest='coverage_threshold', default=200, type=int,
                        help="Minimum coverage sum per window (default: 200)")
    parser.add_argument("--mutrate-threshold", dest='mutrate_threshold', default=0.2, type=float,
                        help="Maximum mutrate per position to include (default: 0.2)")
    parser.add_argument("-p", "--processes", dest='processes', type=int, default=None,
                        help="Number of parallel processes (default: number of CPUs)")

    args = parser.parse_args()
    main(args)
