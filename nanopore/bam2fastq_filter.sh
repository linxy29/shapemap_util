#!/bin/bash
BAM_DIR="./bam"
FASTQ_DIR="./fastq"

# Create fastq directory if it doesn't exist
mkdir -p "$FASTQ_DIR"

# Part 1: Convert all BAM to FASTQ
echo "Converting BAM files to FASTQ..."
for bamfile in "$BAM_DIR"/*.bam; do
    basename=$(basename "$bamfile" .bam)
    outfile="$FASTQ_DIR/${basename}.fastq.gz"
    
    echo "Converting $bamfile -> $outfile"
    samtools fastq -@ 8 "$bamfile" | gzip > "$outfile"
done
echo "✅ All BAM files converted to FASTQ.gz files"

# Part 2: Filter all FASTQ files
echo "Starting FASTQ filtering..."
for fastqfile in "$FASTQ_DIR"/*.fastq.gz; do
    # Skip files that are already filtered (end with _pass.fastq.gz or _fail.fastq.gz)
    if [[ $fastqfile == *"_pass.fastq.gz" ]] || [[ $fastqfile == *"_fail.fastq.gz" ]]; then
        continue
    fi
    
    basename=$(basename "$fastqfile" .fastq.gz)
    pass_file="$FASTQ_DIR/${basename}_pass.fastq.gz"
    fail_file="$FASTQ_DIR/${basename}_fail.fastq.gz"
    
    echo "Filtering reads in $fastqfile"
    python /home/users/astar/gis/linxy/code/shapemap_util/nanopore/filter_reads_V2.py "$fastqfile" "$pass_file" "$fail_file"
done

echo "✅ All FASTQ files have been filtered in: $FASTQ_DIR"
