#!/bin/bash

# Script to extract all pass.fastq.gz files and concatenate them into one file
# Usage: ./extract_pass_fastq.sh [output_filename]

# Set default output filename if not provided
OUTPUT_FILE=${1:-"combined_pass_reads.fastq.gz"}

# Get the current directory (should be the nanopore directory)
NANOPORE_DIR=$(pwd)

echo "Starting extraction of pass.fastq.gz files from: $NANOPORE_DIR"
echo "Output file: $OUTPUT_FILE"

# Check if output file already exists
if [ -f "$OUTPUT_FILE" ]; then
    echo "Warning: Output file $OUTPUT_FILE already exists."
    read -p "Do you want to overwrite it? (y/n): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Operation cancelled."
        exit 1
    fi
    rm "$OUTPUT_FILE"
fi

# Find all pass.fastq.gz files and count them
PASS_FILES=($(find . -name "*_pass.fastq.gz" -type f))
TOTAL_FILES=${#PASS_FILES[@]}

echo "Found $TOTAL_FILES pass.fastq.gz files"

if [ $TOTAL_FILES -eq 0 ]; then
    echo "No pass.fastq.gz files found!"
    exit 1
fi

# List all files that will be processed
echo "Files to be processed:"
for file in "${PASS_FILES[@]}"; do
    echo "  $file"
done

echo
read -p "Proceed with concatenation? (y/n): " -n 1 -r
echo

if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Operation cancelled."
    exit 1
fi

# Process each file
echo "Starting concatenation..."
PROCESSED=0

for file in "${PASS_FILES[@]}"; do
    PROCESSED=$((PROCESSED + 1))
    echo "Processing file $PROCESSED/$TOTAL_FILES: $file"
    
    # Get file size for progress indication
    SIZE=$(du -h "$file" | cut -f1)
    echo "  File size: $SIZE"
    
    # Concatenate the file (since they're gzipped, we can just cat them)
    cat "$file" >> "$OUTPUT_FILE"
    
    if [ $? -eq 0 ]; then
        echo "  ✓ Successfully added"
    else
        echo "  ✗ Error processing file: $file"
        exit 1
    fi
done

echo
echo "Concatenation completed successfully!"
echo "Output file: $OUTPUT_FILE"
echo "Final file size: $(du -h "$OUTPUT_FILE" | cut -f1)"

# Verify the output file
if [ -f "$OUTPUT_FILE" ] && [ -s "$OUTPUT_FILE" ]; then
    echo "✓ Output file created and is not empty"
    
    # Count reads in the final file (optional verification)
    echo "Counting reads in output file (this may take a while)..."
    READ_COUNT=$(zcat "$OUTPUT_FILE" | wc -l | awk '{print $1/4}')
    echo "Total reads in combined file: $READ_COUNT"
else
    echo "✗ Error: Output file is missing or empty"
    exit 1
fi

echo "Done!"