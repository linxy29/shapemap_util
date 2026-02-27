#!/usr/bin/env python3
"""
Fast multi-threaded FASTQ modifier for UMI-tools extracted files.
Moves barcode and UMI from read name to the 5' end of the sequence.
Usage: python read_name_to_sequence.py <input.fastq[.gz]> <output.fastq[.gz]> [quality_char] [threads]
"""

import sys
import gzip
from multiprocessing import Pool, cpu_count
from functools import partial
import io

def open_fastq(filepath):
    """Open FASTQ file, handling both gzipped and uncompressed files."""
    if filepath.endswith('.gz'):
        return gzip.open(filepath, 'rt')
    return open(filepath, 'r')

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

def process_record(record, high_qual='I'):
    """
    Process a single FASTQ record (4 lines).

    Args:
        record: Tuple of (header, sequence, plus, quality)
        high_qual: Phred quality character

    Returns:
        Modified record as a string
    """
    header, sequence, plus, quality = record

    # Extract barcode and UMI from header
    barcode, umi = extract_barcode_umi(header.strip())

    if barcode and umi:
        # Prepend barcode and UMI to sequence
        bc_umi = barcode + umi
        new_sequence = bc_umi + sequence

        # Add high quality scores for barcode and UMI
        qual_addition = high_qual * len(bc_umi)
        new_quality = qual_addition + quality

        return f"{header}{new_sequence}\n{plus}{new_quality}\n"
    else:
        # Return original record if barcode/UMI not found
        return f"{header}{sequence}\n{plus}{quality}\n"

def process_chunk(records, high_qual='I'):
    """Process a chunk of FASTQ records."""
    results = []
    for record in records:
        results.append(process_record(record, high_qual))
    return ''.join(results)

def read_fastq_chunks(filepath, chunk_size=10000):
    """
    Read FASTQ file in chunks for parallel processing.

    Args:
        filepath: Path to FASTQ file
        chunk_size: Number of records per chunk

    Yields:
        Lists of FASTQ records (tuples of 4 lines)
    """
    chunk = []
    with open_fastq(filepath) as f:
        while True:
            header = f.readline()
            if not header:
                if chunk:
                    yield chunk
                break

            sequence = f.readline()
            plus = f.readline()
            quality = f.readline()

            chunk.append((header, sequence.strip(), plus, quality.strip()))

            if len(chunk) >= chunk_size:
                yield chunk
                chunk = []

def modify_fastq_parallel(input_file, output_file, high_qual='I', threads=None, chunk_size=10000):
    """
    Modify FASTQ file using parallel processing.

    Args:
        input_file: Input FASTQ file path
        output_file: Output FASTQ file path
        high_qual: Phred quality character (default 'I' = Q40)
        threads: Number of threads (default: CPU count)
        chunk_size: Records per chunk for parallel processing
    """

    if threads is None:
        threads = cpu_count()

    print(f"Using {threads} threads", file=sys.stderr)

    processed = 0

    with gzip.open(output_file, 'wt', compresslevel=6) as f_out:
        if threads == 1:
            # Single-threaded mode
            for chunk in read_fastq_chunks(input_file, chunk_size):
                result = process_chunk(chunk, high_qual)
                f_out.write(result)
                processed += len(chunk)
                if processed % 100000 == 0:
                    print(f"Processed {processed} reads...", file=sys.stderr)
        else:
            # Multi-threaded mode
            with Pool(threads) as pool:
                process_func = partial(process_chunk, high_qual=high_qual)

                for result in pool.imap(process_func, read_fastq_chunks(input_file, chunk_size), chunksize=1):
                    f_out.write(result)
                    processed += chunk_size  # Approximate
                    if processed % 100000 == 0:
                        print(f"Processed ~{processed} reads...", file=sys.stderr)

    print(f"\nTotal reads processed: ~{processed}", file=sys.stderr)

def modify_fastq_streaming(input_file, output_file, high_qual='I'):
    """
    Single-threaded streaming version with optimized I/O.
    Often faster for smaller files or when I/O is the bottleneck.
    """

    processed = 0
    modified = 0

    # Use larger buffer for better I/O performance
    with open_fastq(input_file) as f_in, \
         gzip.open(output_file, 'wt', compresslevel=6) as f_out:

        # Pre-allocate string buffer
        write_buffer = []
        buffer_size = 1000

        while True:
            header = f_in.readline()
            if not header:
                # Flush remaining buffer
                if write_buffer:
                    f_out.write(''.join(write_buffer))
                break

            sequence = f_in.readline().strip()
            plus = f_in.readline()
            quality = f_in.readline().strip()

            processed += 1

            # Extract barcode and UMI from header
            barcode, umi = extract_barcode_umi(header.strip())

            if barcode and umi:
                bc_umi = barcode + umi
                new_sequence = bc_umi + sequence
                qual_addition = high_qual * len(bc_umi)
                new_quality = qual_addition + quality

                write_buffer.append(f"{header}{new_sequence}\n{plus}{new_quality}\n")
                modified += 1
            else:
                write_buffer.append(f"{header}{sequence}\n{plus}{quality}\n")

            # Flush buffer periodically
            if len(write_buffer) >= buffer_size:
                f_out.write(''.join(write_buffer))
                write_buffer = []

            if processed % 100000 == 0:
                print(f"Processed {processed} reads...", file=sys.stderr)

    print(f"\nTotal reads processed: {processed}", file=sys.stderr)
    print(f"Reads modified: {modified}", file=sys.stderr)

def main():
    if len(sys.argv) < 3:
        print("Usage: python script.py <input.fastq[.gz]> <output.fastq[.gz]> [quality_char] [threads]")
        print("\nQuality character options:")
        print("  I = Q40 (default, 99.99% accuracy)")
        print(f"\nThreads: default = {cpu_count()} (use 1 for single-threaded)")
        print("\nNote: Output will be automatically gzipped")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    quality_char = sys.argv[3] if len(sys.argv) > 3 else 'I'
    threads = int(sys.argv[4]) if len(sys.argv) > 4 else None

    # Ensure output is gzipped
    if not output_file.endswith('.gz'):
        output_file += '.gz'
        print(f"Note: Output will be gzipped as {output_file}", file=sys.stderr)

    print(f"Input: {input_file}", file=sys.stderr)
    print(f"Output: {output_file}", file=sys.stderr)
    print(f"Quality character: {quality_char} (Q{ord(quality_char) - 33})", file=sys.stderr)
    print("", file=sys.stderr)

    # Use streaming for single thread, parallel for multiple threads
    if threads == 1:
        print("Using optimized single-threaded mode", file=sys.stderr)
        modify_fastq_streaming(input_file, output_file, quality_char)
    else:
        modify_fastq_parallel(input_file, output_file, quality_char, threads)

    print("\nDone!", file=sys.stderr)

if __name__ == "__main__":
    main()
