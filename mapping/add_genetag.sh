#!/bin/bash

# Script to add gene tags to BAM files
# Usage: ./add_gene_tag.sh <input.bam> <output.bam>

# Check if correct number of arguments provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 <input.bam> <output.bam>"
    echo "Example: $0 input.bam output_genetag.bam"
    exit 1
fi

# Assign arguments to variables
INPUT_BAM=$1
OUTPUT_BAM=$2

# Check if input file exists
if [ ! -f "$INPUT_BAM" ]; then
    echo "Error: Input file '$INPUT_BAM' not found!"
    exit 1
fi

# Run the pipeline
echo "Processing $INPUT_BAM..."
samtools view -h "$INPUT_BAM" | \
awk 'BEGIN{OFS="\t"}
    /^@/ {print; next}                   # print header lines unchanged
    {   gene=$3;                         # get gene name from column 3
        print $0 "\tXT:Z:" gene          # append gene tag
    }' | \
samtools view -bS - > "$OUTPUT_BAM"

# Check if output was created successfully
if [ -f "$OUTPUT_BAM" ]; then
    echo "Success! Output written to $OUTPUT_BAM"
else
    echo "Error: Failed to create output file"
    exit 1
fi