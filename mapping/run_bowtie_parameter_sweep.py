#!/usr/bin/env python3
"""
Script to run bowtie pipeline with different parameter combinations.
Results are saved to specific folders for each parameter combination.
"""

import subprocess
import itertools
from pathlib import Path
import sys
import os

def modify_bowtie_pipeline(reference_index, input_fastq, output_prefix, threads=4, 
                          bowtie_params="--very-sensitive-local --mp 3,1 --rdg 5,1 --rfg 5,1"):
    """
    Simple step-by-step bowtie2 alignment pipeline.
    """
    
    # Output files
    bam_file = f"{output_prefix}.bam"
    temp_bam = f"{output_prefix}.temp.bam"
    unmapped_file = f"{output_prefix}.unmapped.fq.gz"
    log_file = f"{output_prefix}.log"
    stats_file = f"{output_prefix}.stats"
    flagstat_file = f"{output_prefix}.flagstat"
    
    # Create output directory
    os.makedirs(os.path.dirname(output_prefix), exist_ok=True)
    
    # Step 1: Bowtie2 alignment -> BAM
    print("Step 1: Running alignment...")
    cmd1 = f"bowtie2 {bowtie_params} -p {threads} -x {reference_index} -U {input_fastq} --un-gz {unmapped_file} 2> {log_file} | samtools view -bS > {temp_bam}"
    subprocess.run(cmd1, shell=True, check=True)
    
    # Step 2: Sort BAM
    print("Step 2: Sorting BAM...")
    subprocess.run(['samtools', 'sort', '-@', str(threads), '-o', bam_file, temp_bam], check=True)
    
    # Step 3: Generate stats
    print("Step 3: Generating statistics...")
    subprocess.run(f'samtools flagstat {bam_file} > {flagstat_file}', shell=True, check=True)
    subprocess.run(f'samtools stats {bam_file} > {stats_file}', shell=True, check=True)
    
    # Step 4: Index BAM
    print("Step 4: Indexing BAM...")
    subprocess.run(['samtools', 'index', bam_file], check=True)
    
    # Clean up
    os.remove(temp_bam)
    print("Done!")
    
    return {
        'bam': bam_file,
        'stats': stats_file,
        'flagstat': flagstat_file,
        'log': log_file,
        'unmapped': unmapped_file
    }

def run_full_pipeline_with_params(reference_index, input_fastq, base_output_dir, 
                                 bowtie_params, threads=4):
    """Run the complete pipeline with specified parameters.
    
    Args:
        bowtie_params: String containing all bowtie2 parameters (e.g., "--very-sensitive-local --mp 3,1 --rdg 5,1 --rfg 5,1")
    """
    
    # Create parameter-specific output directory name from the parameter string
    param_name = bowtie_params.replace('--', '').replace(' ', '_').replace(',', '-')
    output_dir = Path(f'{base_output_dir}/bowtie_{param_name}')
    #output_dir = Path(base_output_dir) / param_name
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Set output prefix
    input_basename = Path(input_fastq).stem.replace('.fq', '').replace('.fastq', '')
    output_prefix = output_dir / input_basename
    
    print(f"\n{'='*60}")
    print(f"Running pipeline with parameters:")
    print(f"  Parameters: {bowtie_params}")
    print(f"  Output directory: {output_dir}")
    print(f"{'='*60}")
    
    try:
        # Run bowtie2 alignment
        output_files = modify_bowtie_pipeline(
            reference_index, input_fastq, str(output_prefix), threads,
            bowtie_params
        )
        
        # Import and run analysis functions from original pipeline
        sys.path.append(str(Path(__file__).parent))
        from bowtie_pipeline import parse_raw_lengths, calculate_mapped_lengths, plot_distributions, save_statistics_to_csv
        
        # Calculate lengths
        print("Calculating mapped read lengths...")
        raw_lengths, raw_counts = parse_raw_lengths(output_files['stats'])
        mapped_lengths = calculate_mapped_lengths(output_files['bam'])
        
        # Generate plot
        print("Generating plot...")
        plot_file = f"{output_prefix}_length_distributions.png"
        plot_distributions(raw_lengths, raw_counts, mapped_lengths, plot_file)
        
        # Save statistics
        print("Saving statistics to CSV...")
        save_statistics_to_csv(raw_lengths, raw_counts, mapped_lengths, str(output_prefix))
        
        print(f"Pipeline completed successfully for {param_name}")
        return True
        
    except Exception as e:
        print(f"ERROR: Pipeline failed for {param_name}: {str(e)}")
        return False

def main():
    """Main function to run parameter sweep."""
    
    # Configuration - modify these paths as needed
    reference_index = "/home/users/astar/gis/linxy/scratch/data/references/transcriptome/hg38_ncbi_whole/hg38_ncbi_selected_transcriptome.rmdup"
    input_fastq = "/home/users/astar/gis/linxy/scratch/data/singlecell_shapemap/RHH_only/DMSO/trim/subsample_1mreads_RHH9430_trimmed.fq.gz"
    base_output_dir = "/home/users/astar/gis/linxy/scratch/data/singlecell_shapemap/RHH_only/DMSO/mapping/parameter_sweep"
    threads = 8
    
    # Parameter combinations
    sensitivity_options = [
        "--very-sensitive-local",
        "--very-sensitive", 
        "--sensitive-local",
        "--sensitive"
    ]
    
    mp_options = ["3,1", "6,2"]
    rdg_options = ["5,1", "5,3"]
    rfg_options = ["5,1", "5,3"]
    
    # Create base output directory
    Path(base_output_dir).mkdir(parents=True, exist_ok=True)
    
    # Generate all parameter combinations as strings
    param_combinations = []
    for sens, mp, rdg, rfg in itertools.product(sensitivity_options, mp_options, rdg_options, rfg_options):
        param_string = f"{sens} --mp {mp} --rdg {rdg} --rfg {rfg}"
        param_combinations.append(param_string)
    
    print(f"Total parameter combinations to test: {len(param_combinations)}")
    print("\nParameter combinations:")
    for i, params in enumerate(param_combinations, 1):
        print(f"  {i:2d}. {params}")
    
    # Run all combinations
    successful_runs = 0
    failed_runs = 0
    
    for i, bowtie_params in enumerate(param_combinations, 1):
        print(f"\n\nRunning combination {i}/{len(param_combinations)}")
        
        # Check if output folder already exists
        param_name = bowtie_params.replace('--', '').replace(' ', '_').replace(',', '-')
        output_dir = Path(base_output_dir) / param_name
        
        if output_dir.exists():
            print(f"Output directory {output_dir} already exists. Skipping...")
            successful_runs += 1
            continue
        
        success = run_full_pipeline_with_params(
            reference_index, input_fastq, base_output_dir,
            bowtie_params, threads
        )
        
        if success:
            successful_runs += 1
        else:
            failed_runs += 1
    
    # Summary
    print(f"\n\n{'='*60}")
    print("PARAMETER SWEEP SUMMARY")
    print(f"{'='*60}")
    print(f"Total combinations: {len(param_combinations)}")
    print(f"Successful runs: {successful_runs}")
    print(f"Failed runs: {failed_runs}")
    print(f"Results saved in: {base_output_dir}")
    
    if failed_runs > 0:
        print(f"\nWarning: {failed_runs} runs failed. Check the output for error messages.")

if __name__ == "__main__":
    main()