#!/usr/bin/env python3
"""
Trim a transcriptome-mapped BAM header to only include references (genes)
that have at least one read mapped to them.

This is useful when the BAM was aligned to a full transcriptome but only a
subset of genes received reads, resulting in an unnecessarily large header.

Usage:
    python trim_bam_header.py -i input.bam -o output.bam
    python trim_bam_header.py -i input.bam -o output.bam --gene-list genes.txt
"""

import argparse
import pysam


def get_mapped_genes(bam_path: str) -> set[str]:
    """Return the set of reference names that have at least one mapped read."""
    genes = set()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if not read.is_unmapped:
                genes.add(read.reference_name)
    return genes


def trim_header(bam_path: str, output_path: str, gene_list_path: str | None = None):
    with pysam.AlignmentFile(bam_path, "rb") as bam_in:
        old_header = bam_in.header.to_dict()

        # Step 1: get gene list
        mapped_genes = get_mapped_genes(bam_path)
        print(f"Genes with mapped reads: {len(mapped_genes)}")

        # Optionally save gene list
        if gene_list_path:
            with open(gene_list_path, "w") as f:
                for gene in sorted(mapped_genes):
                    f.write(gene + "\n")
            print(f"Gene list saved to: {gene_list_path}")

        # Step 2: filter @SQ lines
        old_sq = old_header.get("SQ", [])
        new_sq = [sq for sq in old_sq if sq["SN"] in mapped_genes]
        print(f"Header @SQ lines: {len(old_sq)} -> {len(new_sq)}")

        new_header = dict(old_header)
        new_header["SQ"] = new_sq

        # Step 3: write new BAM with trimmed header
        # Build old-to-new reference ID mapping
        old_name_to_new_id = {sq["SN"]: i for i, sq in enumerate(new_sq)}

        with pysam.AlignmentFile(output_path, "wb", header=new_header) as bam_out:
            for read in bam_in.fetch(until_eof=True):
                if read.is_unmapped:
                    bam_out.write(read)
                    continue
                ref_name = read.reference_name
                if ref_name in old_name_to_new_id:
                    read.reference_id = old_name_to_new_id[ref_name]
                    bam_out.write(read)

    # Index the output BAM
    pysam.sort("-o", output_path + ".tmp.bam", output_path)
    import shutil
    shutil.move(output_path + ".tmp.bam", output_path)
    pysam.index(output_path)
    print(f"Output BAM: {output_path}")
    print(f"Output index: {output_path}.bai")


def main():
    parser = argparse.ArgumentParser(
        description="Trim BAM header to only include genes with mapped reads."
    )
    parser.add_argument("-i", "--input", required=True, help="Input BAM file")
    parser.add_argument("-o", "--output", required=True, help="Output BAM file")
    parser.add_argument(
        "--gene-list",
        default=None,
        help="Optional: save the list of mapped genes to this file",
    )
    args = parser.parse_args()

    trim_header(args.input, args.output, args.gene_list)


if __name__ == "__main__":
    main()
