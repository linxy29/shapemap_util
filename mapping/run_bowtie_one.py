#!/usr/bin/env python3
"""
Script to run bowtie pipeline once with specified parameters.
Results are saved to specific folders for each parameter combination.
"""

# Not sure whether it can be used in sbatch, need to test.

import subprocess
import itertools
from pathlib import Path
import sys
import argparse
from run_bowtie_parameter_sweep import run_full_pipeline_with_params

def main():
    """Main function to run bowtie pipeline with specified parameters."""
    parser = argparse.ArgumentParser(description='Run bowtie pipeline with specified parameters.')
    parser.add_argument('--reference_index', type=str, required=False, help='Path to the bowtie reference index')
    parser.add_argument('--input_fastq', type=str, required=False, help='Path to the input FASTQ file')
    parser.add_argument('--base_output_dir', type=str, required=False, help='Base directory for output files')
    parser.add_argument('--bowtie_params', type=str, required=False, help='Additional bowtie parameters as a single string')
    parser.add_argument('--threads', type=int, default=8, help='Number of threads to use')
    args = parser.parse_args()

    run_full_pipeline_with_params(
        args.reference_index,
        args.input_fastq,
        args.base_output_dir,
        args.bowtie_params,
        args.threads
    )

def test():
    """Test function with hardcoded parameters."""
    sys.argv = ['run_bowtie_one.py']  # Clear sys.argv to avoid interference
    
    reference_index = "/home/users/astar/gis/linxy/scratch/data/references/transcriptome/hg38_ncbi_whole/hg38_ncbi_selected_transcriptome.rmdup"
    input_fastq = "/home/users/astar/gis/linxy/scratch/data/singlecell_shapemap/RHH_only/DMSO/trim/subsample_1mreads_RHH9430_trimmed.fq.gz"
    #input_fastq = "/home/users/astar/gis/linxy/scratch/data/singlecell_shapemap/PIPseq_merge/DMSO/PIP2_PIP_DMSO_trimmed.fq.gz"
    base_output_dir = "/home/users/astar/gis/linxy/scratch/data/singlecell_shapemap/RHH_only/DMSO/mapping/parameter_sweep"
    #base_output_dir = "/home/users/astar/gis/linxy/scratch/data/singlecell_shapemap/PIPseq_merge/DMSO/mapping/"
    threads = 8
    bowtie_params = "--sensitive-local --mp 6,2 --rdg 5,3 --rfg 5,3 -k 4"
    run_full_pipeline_with_params(reference_index, input_fastq, base_output_dir, bowtie_params, threads)

if __name__ == "__main__":
    main()
    #test()