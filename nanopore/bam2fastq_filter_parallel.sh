#!/bin/bash
#Usage: bash /home/users/astar/gis/linxy/code/shapemap_util/nanopore/bam2fastq_filter_parallel.sh

BAM_DIR="./bam"
FASTQ_DIR="./fastq"
# Set number of parallel jobs (adjust based on your CPU cores)
MAX_JOBS=2

# Create fastq directory if it doesn't exist
mkdir -p "$FASTQ_DIR"

# Function to convert BAM to FASTQ
convert_bam() {
    local bamfile=$1
    local basename=$(basename "$bamfile" .bam)
    local outfile="$FASTQ_DIR/${basename}.fastq.gz"
    
    echo "Converting $bamfile -> $outfile"
    samtools fastq -@ 2 "$bamfile" | gzip > "$outfile"
    echo "✅ Completed: $outfile"
}

# Function to filter FASTQ
filter_fastq() {
    local fastqfile=$1
    
    # Skip files that are already filtered
    if [[ $fastqfile == *"_pass.fastq.gz" ]] || [[ $fastqfile == *"_fail.fastq.gz" ]]; then
        return
    fi
    
    local basename=$(basename "$fastqfile" .fastq.gz)
    local pass_file="$FASTQ_DIR/${basename}_pass.fastq.gz"
    local fail_file="$FASTQ_DIR/${basename}_fail.fastq.gz"
    
    echo "Filtering reads in $fastqfile"
    python /home/users/astar/gis/linxy/code/shapemap_util/nanopore/filter_reads_V2.py "$fastqfile" "$pass_file" "$fail_file"
    echo "✅ Completed filtering: $fastqfile"
}

# Export functions so they can be used by parallel processes
export -f convert_bam
export -f filter_fastq
export FASTQ_DIR

# Part 1: Convert all BAM to FASTQ in parallel
echo "Converting BAM files to FASTQ in parallel (${MAX_JOBS} jobs)..."

# Method 1: Using GNU Parallel (recommended if available)
if command -v parallel &> /dev/null; then
    find "$BAM_DIR" -name "*.bam" | parallel -j "$MAX_JOBS" convert_bam {}
else
    # Method 2: Using background jobs with job control
    job_count=0
    for bamfile in "$BAM_DIR"/*.bam; do
        convert_bam "$bamfile" &
        
        ((job_count++))
        # Wait when we reach MAX_JOBS
        if (( job_count >= MAX_JOBS )); then
            wait -n  # Wait for any one job to finish
            ((job_count--))
        fi
    done
    wait  # Wait for all remaining jobs
fi

echo "✅ All BAM files converted to FASTQ.gz files"

# Part 2: Filter all FASTQ files in parallel
echo "Starting FASTQ filtering in parallel (${MAX_JOBS} jobs)..."

if command -v parallel &> /dev/null; then
    find "$FASTQ_DIR" -name "*.fastq.gz" ! -name "*_pass.fastq.gz" ! -name "*_fail.fastq.gz" | \
        parallel -j "$MAX_JOBS" filter_fastq {}
else
    job_count=0
    for fastqfile in "$FASTQ_DIR"/*.fastq.gz; do
        # Skip files that are already filtered
        if [[ $fastqfile == *"_pass.fastq.gz" ]] || [[ $fastqfile == *"_fail.fastq.gz" ]]; then
            continue
        fi
        
        filter_fastq "$fastqfile" &
        
        ((job_count++))
        if (( job_count >= MAX_JOBS )); then
            wait -n
            ((job_count--))
        fi
    done
    wait
fi

echo "✅ All FASTQ files have been filtered in: $FASTQ_DIR"
