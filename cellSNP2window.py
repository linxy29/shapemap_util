#!/home/users/astar/gis/gaog/.conda/envs/main/bin/python

## This script is getting from Gao Ge.
## Version: v3.0
## Change input to VCF file format using pysam for efficient reading.
## Process each cell and output combined VCF with sliding window results.
## Order dataframe by gene and pos before processing.

## Important: Might need to add the step for filtering based on mutation rates.


import argparse
import pysam
import gzip
from collections import defaultdict

__doc__="""
Calculate sliding window mutation rates from cellSNP VCF files.

Description:
    This script reads a cellSNP VCF file using pysam, calculates mutation rates
    within sliding windows for each cell, and outputs a combined VCF file with
    all cells as sample columns.

Input:
    cellSNP VCF file (gzipped, with .tbi index) containing per-cell mutation data.
    Expected FORMAT fields: GT:AD:DP:OTH:PL:ALL

Output:
    Gzipped VCF file with one row per window, all cells as sample columns.
    FORMAT: DP:AD:MR (coverage, mutant, mutrate)

Usage:
    python cellSNP2window.py -i cellSNP.cells.vcf.gz -o output.vcf.gz -w 10 -s 10 --coverage-threshold 200 --mutrate-threshold 0.2

Arguments:
    -i, --input              Input cellSNP VCF file (gzipped)
    -o, --output             Output VCF file (gzipped)
    -w, --window_length      Window length in bp (default: 10)
    -s, --step               Step size in bp (default: 10)
    --coverage-threshold     Minimum coverage sum per window (default: 200)
    --mutrate-threshold      Maximum mutrate per position to include (default: 0.2)
"""
__version__="v3.0"
__author__="linxy"
__last_modify__="26-Nov-2025"


def read_vcf_pysam(vcf_path):
    """
    Read cellSNP VCF file using pysam for efficient parsing.

    Parameters:
        vcf_path: str
            Path to the gzipped VCF file

    Returns:
        data: dict
            {gene: {pos: {sample_name: {'DP': int, 'AD': int}}}}
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

    # Structure: {gene: {pos: {sample: {'DP': val, 'AD': val}}}}
    data = defaultdict(lambda: defaultdict(dict))

    for record in vcf:
        chrom = record.chrom
        pos = record.pos

        for sample in samples:
            sample_data = record.samples[sample]
            # Handle missing data (./.)
            dp = sample_data.get('DP', None)
            ad = sample_data.get('AD', None)

            if dp is not None and ad is not None:
                data[chrom][pos][sample] = {
                    'DP': dp,
                    'AD': ad
                }

    vcf.close()
    return dict(data), samples, contigs


def calculate_windows_per_cell(gene_data, samples, win_len, step, coverage_threshold, mutrate_threshold):
    """
    Calculate sliding windows for all cells in a gene.

    Parameters:
        gene_data: dict
            {pos: {sample: {'DP': int, 'AD': int}}}
        samples: list
            List of sample names
        win_len: int
            Window length in bp
        step: int
            Step size in bp
        coverage_threshold: int
            Minimum coverage per window
        mutrate_threshold: float
            Maximum mutrate per position to include

    Returns:
        windows: list of dict
            [{win_start, win_end, sample_data: {sample: {DP, AD, MR}}}]
    """
    if not gene_data:
        return []

    positions = sorted(gene_data.keys())
    min_pos = positions[0]
    max_pos = positions[-1]

    windows = []
    win_start = min_pos

    while win_start <= max_pos:
        win_end = win_start + win_len - 1

        # Get positions within this window
        window_positions = [p for p in positions if win_start <= p <= win_end]

        if window_positions:
            sample_data = {}

            for sample in samples:
                dp_sum = 0
                ad_sum = 0

                for pos in window_positions:
                    if sample in gene_data[pos]:
                        cell_data = gene_data[pos][sample]
                        dp = cell_data['DP']
                        ad = cell_data['AD']

                        # Calculate position mutrate for filtering
                        pos_mutrate = ad / dp if dp > 0 else 0

                        if pos_mutrate <= mutrate_threshold:
                            dp_sum += dp
                            ad_sum += ad

                # Apply coverage threshold per cell
                if dp_sum >= coverage_threshold:
                    mutrate = round(ad_sum / dp_sum, 5) if dp_sum > 0 else 0
                    sample_data[sample] = {
                        'DP': dp_sum,
                        'AD': ad_sum,
                        'MR': mutrate
                    }
                else:
                    # Mark as missing data
                    sample_data[sample] = None

            # Only add window if at least one cell passes threshold
            if any(v is not None for v in sample_data.values()):
                windows.append({
                    'win_start': win_start,
                    'win_end': win_end,
                    'sample_data': sample_data
                })

        win_start += step

    return windows


def write_vcf_output(output_path, windows_by_gene, samples, contigs):
    """
    Write sliding window results as a VCF file.

    Parameters:
        output_path: str
            Path to output gzipped VCF file
        windows_by_gene: dict
            {gene: [window_dicts]}
        samples: list
            List of sample names
        contigs: list
            List of contig names
    """
    with gzip.open(output_path, 'wt') as f:
        # Write VCF header
        f.write("##fileformat=VCFv4.2\n")
        f.write("##source=window.py_v3.0\n")
        f.write('##INFO=<ID=WIN_END,Number=1,Type=Integer,Description="Window end position">\n')
        f.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total coverage in window">\n')
        f.write('##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Total alternate counts in window">\n')
        f.write('##FORMAT=<ID=MR,Number=1,Type=Float,Description="Mutation rate (AD/DP)">\n')

        # Write contig lines
        for contig in contigs:
            f.write(f'##contig=<ID={contig}>\n')

        # Write header line
        header_cols = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
        header_cols.extend(samples)
        f.write('\t'.join(header_cols) + '\n')

        # Write data rows
        for gene in windows_by_gene:
            for window in windows_by_gene[gene]:
                win_start = window['win_start']
                win_end = window['win_end']
                sample_data = window['sample_data']

                # Build row
                row = [
                    gene,                    # CHROM
                    str(win_start),          # POS
                    '.',                     # ID
                    'N',                     # REF (placeholder)
                    '.',                     # ALT
                    '.',                     # QUAL
                    'PASS',                  # FILTER
                    f'WIN_END={win_end}',    # INFO
                    'DP:AD:MR'               # FORMAT
                ]

                # Add sample columns
                for sample in samples:
                    if sample_data.get(sample) is not None:
                        sd = sample_data[sample]
                        row.append(f"{sd['DP']}:{sd['AD']}:{sd['MR']}")
                    else:
                        row.append('.:.:.')

                f.write('\t'.join(row) + '\n')


def main(args):
    '''
    Read VCF file, calculate sliding windows for each cell, and output combined VCF.
    '''
    print(f"Reading VCF file: {args.input_file}")
    data, samples, contigs = read_vcf_pysam(args.input_file)
    print(f"Found {len(samples)} cells and {len(data)} genes/chromosomes")

    # Calculate windows for each gene
    windows_by_gene = {}
    total_windows = 0

    for gene in data:
        windows = calculate_windows_per_cell(
            data[gene],
            samples,
            args.win_len,
            args.step,
            args.coverage_threshold,
            args.mutrate_threshold
        )
        if windows:
            windows_by_gene[gene] = windows
            total_windows += len(windows)

    print(f"Calculated {total_windows} windows across {len(windows_by_gene)} genes")

    # Write output VCF
    write_vcf_output(args.output_file, windows_by_gene, samples, contigs)
    print(f"Output written to: {args.output_file}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", dest='input_file', required=True,
                        help="Input cellSNP VCF file (gzipped)")
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

    args = parser.parse_args()
    main(args)