#!/usr/bin/env python3
"""
Fast FASTQ modifier for UMI-tools extracted files.
Moves barcode and UMI from read name to the 5' end of the sequence.

v2: Automatically detects numeric CID barcodes (e.g., 10236:14346) and
    converts them to nucleotide sequences before prepending.
    Read names are always kept unchanged.

v3: Uses pigz subprocesses for fast parallel (de)compression instead of
    Python's built-in gzip module (~10x faster on large files).
    Removed multiprocessing Pool (the computation is lightweight; pigz
    already parallelises the I/O bottleneck).

Encoding for numeric CID (each character -> dinucleotide):
  0->AA, 1->AC, 2->AG, 3->AT, 4->CA, 5->CC, 6->CG, 7->CT, 8->GA, 9->GC, :->GT

Usage: python read_name_to_sequence_v2.py <input.fastq[.gz]> <output.fastq[.gz]> [quality_char] [threads]
"""

import sys
import subprocess

CHAR_TO_SEQ = {
    '0': 'AA', '1': 'AC', '2': 'AG', '3': 'AT',
    '4': 'CA', '5': 'CC', '6': 'CG', '7': 'CT',
    '8': 'GA', '9': 'GC', ':': 'GT'
}

NUMERIC_CHARS = set(CHAR_TO_SEQ.keys())


def is_numeric_cid(barcode):
    """Check if barcode is a numeric CID (digits and colons only)."""
    return all(c in NUMERIC_CHARS for c in barcode)


def encode_cid(cid_str):
    """Convert numeric CID string (e.g., '10236:14346') to nucleotide sequence."""
    return ''.join(CHAR_TO_SEQ[c] for c in cid_str)


def extract_barcode_umi(read_name):
    """
    Extract barcode and UMI from read name.
    Expected format: @READ_ID_BARCODE_UMI or similar

    Returns: (barcode, umi) or (None, None) if not found
    """
    parts = read_name.split(" ")[0].split('_')
    if len(parts) >= 2:
        umi = parts[-1]
        barcode = parts[-2]
        return barcode, umi
    return None, None


def modify_fastq(input_file, output_file, high_qual='I', threads=4):
    """
    Modify FASTQ file using pigz subprocesses for fast parallel I/O.
    """
    # Open input: pigz for .gz, plain open otherwise
    if input_file.endswith('.gz'):
        decomp = subprocess.Popen(
            ['pigz', '-dc', '-p', str(threads), input_file],
            stdout=subprocess.PIPE, bufsize=1024*1024
        )
        fin = decomp.stdout
    else:
        fin = open(input_file, 'rb')
        decomp = None

    # Open output: pigz for .gz, plain open otherwise
    out_fh = None
    if output_file.endswith('.gz'):
        out_fh = open(output_file, 'wb')
        comp = subprocess.Popen(
            ['pigz', '-p', str(threads), '-c'],
            stdin=subprocess.PIPE, stdout=out_fh,
            bufsize=1024*1024
        )
        fout = comp.stdin
    else:
        fout = open(output_file, 'wb')
        comp = None

    processed = 0
    modified = 0
    numeric_converted = 0
    high_qual_bytes = high_qual.encode()

    try:
        while True:
            header = fin.readline()
            if not header:
                break

            sequence = fin.readline().rstrip(b'\n')
            plus = fin.readline()
            quality = fin.readline().rstrip(b'\n')

            processed += 1

            barcode, umi = extract_barcode_umi(header.decode().strip())

            if barcode and umi:
                if is_numeric_cid(barcode):
                    barcode = encode_cid(barcode)
                    numeric_converted += 1

                bc_umi = barcode + umi
                bc_umi_bytes = bc_umi.encode()
                qual_add = high_qual_bytes * len(bc_umi)

                fout.write(header)
                fout.write(bc_umi_bytes + sequence + b'\n')
                fout.write(plus)
                fout.write(qual_add + quality + b'\n')
                modified += 1
            else:
                fout.write(header)
                fout.write(sequence + b'\n')
                fout.write(plus)
                fout.write(quality + b'\n')

            if processed % 1000000 == 0:
                print(f"Processed {processed:,} reads...", file=sys.stderr)

    finally:
        fin.close()
        if fout is not None:
            fout.close()
        if decomp is not None:
            decomp.wait()
        if comp is not None:
            comp.wait()
        if out_fh is not None:
            out_fh.close()

    print(f"\nTotal reads processed: {processed:,}", file=sys.stderr)
    print(f"Reads modified: {modified:,}", file=sys.stderr)
    if numeric_converted > 0:
        print(f"Numeric CID converted: {numeric_converted:,}", file=sys.stderr)


def main():
    if len(sys.argv) < 3:
        print("Usage: python script.py <input.fastq[.gz]> <output.fastq[.gz]> [quality_char] [threads]")
        print("\nQuality character options:")
        print("  I = Q40 (default, 99.99% accuracy)")
        print("\nThreads: default = 4 (used by pigz for parallel compression)")
        print("\nNote: Output will be automatically gzipped")
        print("\nv3: Uses pigz for fast parallel (de)compression.")
        print("    Numeric CID barcodes (e.g., 10236:14346) are automatically")
        print("    converted to nucleotide sequences. Read names are unchanged.")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    quality_char = sys.argv[3] if len(sys.argv) > 3 else 'I'
    threads = int(sys.argv[4]) if len(sys.argv) > 4 else 4

    # Ensure output is gzipped
    if not output_file.endswith('.gz'):
        output_file += '.gz'
        print(f"Note: Output will be gzipped as {output_file}", file=sys.stderr)

    print(f"Input: {input_file}", file=sys.stderr)
    print(f"Output: {output_file}", file=sys.stderr)
    print(f"Quality character: {quality_char} (Q{ord(quality_char) - 33})", file=sys.stderr)
    print(f"Threads: {threads}", file=sys.stderr)
    print("", file=sys.stderr)

    modify_fastq(input_file, output_file, quality_char, threads)

    print("\nDone!", file=sys.stderr)


if __name__ == "__main__":
    main()
