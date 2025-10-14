#!/bin/bash

# Script to extract all pass.fastq.gz files and concatenate them into separate files
# based on whether the folder contains "DMSO" or "2A3"
# Usage: ./extract_pass_fastq.sh [--fast]
# Options:
#   --fast: Use file sizes for statistics instead of counting reads (much faster)

# Check for --fast option
FAST_MODE=false
if [[ "$1" == "--fast" ]]; then
    FAST_MODE=true
    echo "Running in fast mode (using file sizes for statistics)"
fi

# Set output filenames
DMSO_OUTPUT="combined_DMSO_pass_reads.fastq.gz"
F2A3_OUTPUT="combined_2A3_pass_reads.fastq.gz"
STATS_FILE="sequencing_statistics.txt"

# Get the current directory (should be the nanopore directory)
NANOPORE_DIR=$(pwd)

echo "Starting extraction of pass.fastq.gz files from: $NANOPORE_DIR"
echo "DMSO output file: $DMSO_OUTPUT"
echo "2A3 output file: $F2A3_OUTPUT"
echo "Statistics file: $STATS_FILE"

# Check if output files already exist
for OUTPUT_FILE in "$DMSO_OUTPUT" "$F2A3_OUTPUT"; do
    if [ -f "$OUTPUT_FILE" ]; then
        echo "Warning: Output file $OUTPUT_FILE already exists."
        read -p "Do you want to overwrite existing files? (y/n): " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            echo "Operation cancelled."
            exit 1
        fi
        rm -f "$DMSO_OUTPUT" "$F2A3_OUTPUT"
        break
    fi
done

# Find all pass.fastq.gz and fail.fastq.gz files and separate them by DMSO/2A3
DMSO_PASS_FILES=()
F2A3_PASS_FILES=()
DMSO_FAIL_FILES=()
F2A3_FAIL_FILES=()

# Find pass files
while IFS= read -r -d '' file; do
    if [[ "$file" == *"DMSO"* ]]; then
        DMSO_PASS_FILES+=("$file")
    elif [[ "$file" == *"2A3"* ]]; then
        F2A3_PASS_FILES+=("$file")
    else
        echo "Warning: Pass file doesn't contain DMSO or 2A3 in path: $file"
    fi
done < <(find . -name "*_pass.fastq.gz" -type f -print0)

# Find fail files
while IFS= read -r -d '' file; do
    if [[ "$file" == *"DMSO"* ]]; then
        DMSO_FAIL_FILES+=("$file")
    elif [[ "$file" == *"2A3"* ]]; then
        F2A3_FAIL_FILES+=("$file")
    else
        echo "Warning: Fail file doesn't contain DMSO or 2A3 in path: $file"
    fi
done < <(find . -name "*_fail.fastq.gz" -type f -print0)

TOTAL_DMSO_PASS=${#DMSO_PASS_FILES[@]}
TOTAL_2A3_PASS=${#F2A3_PASS_FILES[@]}
TOTAL_DMSO_FAIL=${#DMSO_FAIL_FILES[@]}
TOTAL_2A3_FAIL=${#F2A3_FAIL_FILES[@]}
TOTAL_PASS=$((TOTAL_DMSO_PASS + TOTAL_2A3_PASS))
TOTAL_FAIL=$((TOTAL_DMSO_FAIL + TOTAL_2A3_FAIL))

echo "Found files:"
echo "  - DMSO pass files: $TOTAL_DMSO_PASS"
echo "  - DMSO fail files: $TOTAL_DMSO_FAIL"
echo "  - 2A3 pass files: $TOTAL_2A3_PASS"
echo "  - 2A3 fail files: $TOTAL_2A3_FAIL"
echo "  - Total pass files: $TOTAL_PASS"
echo "  - Total fail files: $TOTAL_FAIL"

if [ $TOTAL_PASS -eq 0 ]; then
    echo "No pass.fastq.gz files found!"
    exit 1
fi

# List all files that will be processed
echo
echo "DMSO pass files to be processed:"
for file in "${DMSO_PASS_FILES[@]}"; do
    echo "  $file"
done

echo
echo "2A3 pass files to be processed:"
for file in "${F2A3_PASS_FILES[@]}"; do
    echo "  $file"
done

echo
read -p "Proceed with concatenation? (y/n): " -n 1 -r
echo

if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Operation cancelled."
    exit 1
fi

# Initialize statistics variables
DMSO_PASS_SIZE=0
DMSO_FAIL_SIZE=0
F2A3_PASS_SIZE=0
F2A3_FAIL_SIZE=0
DMSO_PASS_READS=0
DMSO_FAIL_READS=0
F2A3_PASS_READS=0
F2A3_FAIL_READS=0

# Calculate total file sizes for fail files
echo
echo "Calculating fail file sizes..."
for file in "${DMSO_FAIL_FILES[@]}"; do
    SIZE=$(stat -c%s "$file" 2>/dev/null || stat -f%z "$file" 2>/dev/null)
    DMSO_FAIL_SIZE=$((DMSO_FAIL_SIZE + SIZE))
done

for file in "${F2A3_FAIL_FILES[@]}"; do
    SIZE=$(stat -c%s "$file" 2>/dev/null || stat -f%z "$file" 2>/dev/null)
    F2A3_FAIL_SIZE=$((F2A3_FAIL_SIZE + SIZE))
done

# Process DMSO pass files
if [ $TOTAL_DMSO_PASS -gt 0 ]; then
    echo
    echo "Processing DMSO pass files..."
    PROCESSED=0

    for file in "${DMSO_PASS_FILES[@]}"; do
        PROCESSED=$((PROCESSED + 1))
        echo "Processing DMSO pass file $PROCESSED/$TOTAL_DMSO_PASS: $file"

        # Get file size for progress indication and statistics
        SIZE=$(stat -c%s "$file" 2>/dev/null || stat -f%z "$file" 2>/dev/null)
        DMSO_PASS_SIZE=$((DMSO_PASS_SIZE + SIZE))
        echo "  File size: $(du -h "$file" | cut -f1)"

        # Concatenate the file (since they're gzipped, we can just cat them)
        cat "$file" >> "$DMSO_OUTPUT"

        if [ $? -eq 0 ]; then
            echo "  ✓ Successfully added"
        else
            echo "  ✗ Error processing file: $file"
            exit 1
        fi
    done

    echo "DMSO concatenation completed!"
    echo "DMSO output file size: $(du -h "$DMSO_OUTPUT" | cut -f1)"
else
    echo "No DMSO pass files to process."
fi

# Process 2A3 pass files
if [ $TOTAL_2A3_PASS -gt 0 ]; then
    echo
    echo "Processing 2A3 pass files..."
    PROCESSED=0

    for file in "${F2A3_PASS_FILES[@]}"; do
        PROCESSED=$((PROCESSED + 1))
        echo "Processing 2A3 pass file $PROCESSED/$TOTAL_2A3_PASS: $file"

        # Get file size for progress indication and statistics
        SIZE=$(stat -c%s "$file" 2>/dev/null || stat -f%z "$file" 2>/dev/null)
        F2A3_PASS_SIZE=$((F2A3_PASS_SIZE + SIZE))
        echo "  File size: $(du -h "$file" | cut -f1)"

        # Concatenate the file (since they're gzipped, we can just cat them)
        cat "$file" >> "$F2A3_OUTPUT"

        if [ $? -eq 0 ]; then
            echo "  ✓ Successfully added"
        else
            echo "  ✗ Error processing file: $file"
            exit 1
        fi
    done

    echo "2A3 concatenation completed!"
    echo "2A3 output file size: $(du -h "$F2A3_OUTPUT" | cut -f1)"
else
    echo "No 2A3 pass files to process."
fi

echo
echo "All concatenations completed successfully!"

# Count reads or use file sizes for statistics
echo
echo "Generating statistics..."

if [ "$FAST_MODE" = true ]; then
    echo "Using file sizes for approximation..."
else
    # Count actual reads in output files
    echo "Counting reads (this may take a while)..."

    if [ -f "$DMSO_OUTPUT" ] && [ -s "$DMSO_OUTPUT" ]; then
        echo "  Counting DMSO pass reads..."
        DMSO_PASS_READS=$(zcat "$DMSO_OUTPUT" | wc -l | awk '{print $1/4}')
    fi

    if [ -f "$F2A3_OUTPUT" ] && [ -s "$F2A3_OUTPUT" ]; then
        echo "  Counting 2A3 pass reads..."
        F2A3_PASS_READS=$(zcat "$F2A3_OUTPUT" | wc -l | awk '{print $1/4}')
    fi

    # Count fail reads if not in fast mode
    if [ $TOTAL_DMSO_FAIL -gt 0 ]; then
        echo "  Counting DMSO fail reads..."
        for file in "${DMSO_FAIL_FILES[@]}"; do
            COUNT=$(zcat "$file" | wc -l | awk '{print $1/4}')
            DMSO_FAIL_READS=$((DMSO_FAIL_READS + COUNT))
        done
    fi

    if [ $TOTAL_2A3_FAIL -gt 0 ]; then
        echo "  Counting 2A3 fail reads..."
        for file in "${F2A3_FAIL_FILES[@]}"; do
            COUNT=$(zcat "$file" | wc -l | awk '{print $1/4}')
            F2A3_FAIL_READS=$((F2A3_FAIL_READS + COUNT))
        done
    fi
fi

# Generate statistics report
echo
echo "Creating statistics report..."

# Function to format large numbers with commas
format_number() {
    echo $1 | sed ':a;s/\B[0-9]\{3\}\>/,&/;ta'
}

# Function to convert bytes to human readable
bytes_to_human() {
    local bytes=$1
    if [ $bytes -gt 1073741824 ]; then
        echo "$(echo "scale=2; $bytes/1073741824" | bc) GB"
    elif [ $bytes -gt 1048576 ]; then
        echo "$(echo "scale=2; $bytes/1048576" | bc) MB"
    else
        echo "$(echo "scale=2; $bytes/1024" | bc) KB"
    fi
}

# Write statistics to file
cat > "$STATS_FILE" << EOF
========================================
Sequencing Statistics Report
========================================
Generated: $(date)
Directory: $NANOPORE_DIR
Mode: $(if [ "$FAST_MODE" = true ]; then echo "Fast (file size approximation)"; else echo "Accurate (read counting)"; fi)

========================================
FILE COUNTS
========================================
DMSO:
  Pass files: $TOTAL_DMSO_PASS
  Fail files: $TOTAL_DMSO_FAIL

2A3:
  Pass files: $TOTAL_2A3_PASS
  Fail files: $TOTAL_2A3_FAIL

Total:
  Pass files: $TOTAL_PASS
  Fail files: $TOTAL_FAIL

========================================
FILE SIZES
========================================
DMSO:
  Pass size: $(bytes_to_human $DMSO_PASS_SIZE) ($(format_number $DMSO_PASS_SIZE) bytes)
  Fail size: $(bytes_to_human $DMSO_FAIL_SIZE) ($(format_number $DMSO_FAIL_SIZE) bytes)
  Total size: $(bytes_to_human $((DMSO_PASS_SIZE + DMSO_FAIL_SIZE)))
  Pass proportion: $(echo "scale=2; $DMSO_PASS_SIZE * 100 / ($DMSO_PASS_SIZE + $DMSO_FAIL_SIZE)" | bc)%

2A3:
  Pass size: $(bytes_to_human $F2A3_PASS_SIZE) ($(format_number $F2A3_PASS_SIZE) bytes)
  Fail size: $(bytes_to_human $F2A3_FAIL_SIZE) ($(format_number $F2A3_FAIL_SIZE) bytes)
  Total size: $(bytes_to_human $((F2A3_PASS_SIZE + F2A3_FAIL_SIZE)))
  Pass proportion: $(echo "scale=2; $F2A3_PASS_SIZE * 100 / ($F2A3_PASS_SIZE + $F2A3_FAIL_SIZE)" | bc)%

Overall:
  Total pass size: $(bytes_to_human $((DMSO_PASS_SIZE + F2A3_PASS_SIZE)))
  Total fail size: $(bytes_to_human $((DMSO_FAIL_SIZE + F2A3_FAIL_SIZE)))
  Total size: $(bytes_to_human $((DMSO_PASS_SIZE + F2A3_PASS_SIZE + DMSO_FAIL_SIZE + F2A3_FAIL_SIZE)))
  Pass proportion: $(echo "scale=2; ($DMSO_PASS_SIZE + $F2A3_PASS_SIZE) * 100 / ($DMSO_PASS_SIZE + $F2A3_PASS_SIZE + $DMSO_FAIL_SIZE + $F2A3_FAIL_SIZE)" | bc)%

EOF

# Add read counts if available
if [ "$FAST_MODE" = false ]; then
    cat >> "$STATS_FILE" << EOF
========================================
READ COUNTS
========================================
DMSO:
  Pass reads: $(format_number $DMSO_PASS_READS)
  Fail reads: $(format_number $DMSO_FAIL_READS)
  Total reads: $(format_number $((DMSO_PASS_READS + DMSO_FAIL_READS)))
  Pass rate: $(if [ $((DMSO_PASS_READS + DMSO_FAIL_READS)) -gt 0 ]; then echo "scale=2; $DMSO_PASS_READS * 100 / ($DMSO_PASS_READS + $DMSO_FAIL_READS)" | bc; else echo "N/A"; fi)%

2A3:
  Pass reads: $(format_number $F2A3_PASS_READS)
  Fail reads: $(format_number $F2A3_FAIL_READS)
  Total reads: $(format_number $((F2A3_PASS_READS + F2A3_FAIL_READS)))
  Pass rate: $(if [ $((F2A3_PASS_READS + F2A3_FAIL_READS)) -gt 0 ]; then echo "scale=2; $F2A3_PASS_READS * 100 / ($F2A3_PASS_READS + $F2A3_FAIL_READS)" | bc; else echo "N/A"; fi)%

Overall:
  Total pass reads: $(format_number $((DMSO_PASS_READS + F2A3_PASS_READS)))
  Total fail reads: $(format_number $((DMSO_FAIL_READS + F2A3_FAIL_READS)))
  Total reads: $(format_number $((DMSO_PASS_READS + F2A3_PASS_READS + DMSO_FAIL_READS + F2A3_FAIL_READS)))
  Pass rate: $(if [ $((DMSO_PASS_READS + F2A3_PASS_READS + DMSO_FAIL_READS + F2A3_FAIL_READS)) -gt 0 ]; then echo "scale=2; ($DMSO_PASS_READS + $F2A3_PASS_READS) * 100 / ($DMSO_PASS_READS + $F2A3_PASS_READS + $DMSO_FAIL_READS + $F2A3_FAIL_READS)" | bc; else echo "N/A"; fi)%

EOF
else
    cat >> "$STATS_FILE" << EOF
========================================
ESTIMATED READ COUNTS (based on file size)
========================================
Note: These are rough estimates assuming average compression ratio.
For accurate counts, run without --fast option.

Estimated reads (assuming ~1000 bytes per read compressed):
  DMSO pass: ~$(format_number $((DMSO_PASS_SIZE / 1000)))
  DMSO fail: ~$(format_number $((DMSO_FAIL_SIZE / 1000)))
  2A3 pass: ~$(format_number $((F2A3_PASS_SIZE / 1000)))
  2A3 fail: ~$(format_number $((F2A3_FAIL_SIZE / 1000)))

EOF
fi

cat >> "$STATS_FILE" << EOF
========================================
OUTPUT FILES
========================================
DMSO combined: $DMSO_OUTPUT
2A3 combined: $F2A3_OUTPUT
Statistics report: $STATS_FILE
========================================
EOF

# Display statistics summary
echo
echo "✓ Statistics report saved to: $STATS_FILE"
echo
echo "Summary:"
cat "$STATS_FILE" | grep -E "Pass rate:|Pass proportion:" | head -6

echo
echo "Done!"