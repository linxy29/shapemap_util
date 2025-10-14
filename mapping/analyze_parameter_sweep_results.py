#!/usr/bin/env python3
"""
Script to collect and visualize results from bowtie parameter sweep.
Collects data from .log files and summary_statistics.csv files.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
import re
import glob

def parse_log_file(log_file):
    """Parse bowtie2 log file to extract alignment statistics."""
    stats = {}
    
    try:
        with open(log_file, 'r') as f:
            content = f.read()
        
        # Extract total reads
        total_match = re.search(r'(\d+) reads; of these:', content)
        if total_match:
            stats['total_reads'] = int(total_match.group(1))
        
        # Extract alignment rates
        overall_match = re.search(r'([\d.]+)% overall alignment rate', content)
        if overall_match:
            stats['overall_alignment_rate'] = float(overall_match.group(1))
        
        # Extract aligned exactly 1 time
        exact_1_match = re.search(r'(\d+) \([\d.]+%\) aligned exactly 1 time', content)
        if exact_1_match:
            stats['aligned_exactly_1'] = int(exact_1_match.group(1))
        
        # Extract aligned >1 times
        multi_match = re.search(r'(\d+) \([\d.]+%\) aligned >1 times', content)
        if multi_match:
            stats['aligned_multiple'] = int(multi_match.group(1))
        
        # Extract unaligned reads
        unaligned_match = re.search(r'(\d+) \([\d.]+%\) were unaligned', content)
        if unaligned_match:
            stats['unaligned'] = int(unaligned_match.group(1))
        
        # Calculate additional metrics
        if 'total_reads' in stats:
            if 'aligned_exactly_1' in stats and 'aligned_multiple' in stats:
                stats['total_aligned'] = stats['aligned_exactly_1'] + stats['aligned_multiple']
                stats['alignment_rate_calculated'] = (stats['total_aligned'] / stats['total_reads']) * 100
                stats['aligned_exactly_1_rate'] = (stats['aligned_exactly_1'] / stats['total_reads']) * 100
                stats['aligned_multiple_rate'] = (stats['aligned_multiple'] / stats['total_reads']) * 100
        
    except Exception as e:
        print(f"Error parsing log file {log_file}: {e}")
    
    return stats

def parse_summary_statistics(csv_file):
    """Parse summary_statistics.csv file."""
    try:
        df = pd.read_csv(csv_file)
        return df
    except Exception as e:
        print(f"Error parsing CSV file {csv_file}: {e}")
        return pd.DataFrame()

def extract_parameters_from_path(path):
    """Extract bowtie2 parameters from the directory path."""
    folder_name = Path(path).name
    
    # Parse folder name to extract parameters
    # Expected format: sensitivity_mp_X-Y_rdg_A-B_rfg_C-D
    parts = folder_name.split('_')
    
    params = {}
    
    # Find parameter indices by looking for the parameter labels
    mp_idx = -1
    rdg_idx = -1
    rfg_idx = -1
    
    for i, part in enumerate(parts):
        if part == 'mp':
            mp_idx = i
        elif part == 'rdg':
            rdg_idx = i
        elif part == 'rfg':
            rfg_idx = i
    
    # Extract sensitivity (everything before 'mp')
    if mp_idx > 0:
        sensitivity_parts = parts[:mp_idx]
        params['sensitivity'] = '-'.join(sensitivity_parts)
    
    # Extract MP parameter (value after 'mp' label)
    if mp_idx >= 0 and mp_idx + 1 < len(parts):
        params['mp'] = parts[mp_idx + 1].replace('-', ',')
    
    # Extract RDG parameter (value after 'rdg' label)
    if rdg_idx >= 0 and rdg_idx + 1 < len(parts):
        params['rdg'] = parts[rdg_idx + 1].replace('-', ',')
    
    # Extract RFG parameter (value after 'rfg' label)
    if rfg_idx >= 0 and rfg_idx + 1 < len(parts):
        params['rfg'] = parts[rfg_idx + 1].replace('-', ',')
    
    return params

def collect_all_results(base_dir):
    """Collect all results from parameter sweep directories."""
    base_path = Path(base_dir)
    all_data = []
    
    # Find all subdirectories
    for param_dir in base_path.glob('*'):
        if param_dir.is_dir():
            print(f"Processing directory: {param_dir.name}")
            
            # Extract parameters from directory name
            params = extract_parameters_from_path(param_dir)
            
            # Find log files
            log_files = list(param_dir.glob('*.log'))
            summary_files = list(param_dir.glob('*summary_statistics.csv'))
            
            for log_file in log_files:
                # Parse log file
                log_stats = parse_log_file(log_file)
                
                # Find corresponding summary file
                base_name = log_file.stem
                summary_file = param_dir / f"{base_name}_summary_statistics.csv"
                
                summary_stats = {}
                if summary_file.exists():
                    summary_df = parse_summary_statistics(summary_file)
                    if not summary_df.empty:
                        # Extract key statistics for raw and mapped reads
                        for _, row in summary_df.iterrows():
                            data_type = row['data_type']
                            prefix = 'raw' if 'raw' in data_type else 'mapped'
                            summary_stats[f'{prefix}_total_reads'] = row['total_reads']
                            summary_stats[f'{prefix}_mean_length'] = row['mean_length']
                            summary_stats[f'{prefix}_median_length'] = row['median_length']
                            summary_stats[f'{prefix}_std_length'] = row['std_length']
                            summary_stats[f'{prefix}_min_length'] = row['min_length']
                            summary_stats[f'{prefix}_max_length'] = row['max_length']
                
                # Combine all data
                combined_data = {
                    'parameter_combination': param_dir.name,
                    'log_file': str(log_file),
                    **params,
                    **log_stats,
                    **summary_stats
                }
                
                all_data.append(combined_data)
    
    return pd.DataFrame(all_data)

def create_visualizations(df, output_dir):
    """Create comprehensive visualizations of the results."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Set style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # 1. Overall alignment rates by sensitivity
    if 'overall_alignment_rate' in df.columns and 'sensitivity' in df.columns:
        plt.figure(figsize=(12, 8))
        sns.boxplot(data=df, x='sensitivity', y='overall_alignment_rate')
        plt.title('Overall Alignment Rate by Sensitivity Setting')
        plt.xlabel('Sensitivity Setting')
        plt.ylabel('Overall Alignment Rate (%)')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(output_path / 'alignment_rate_by_sensitivity.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    # 2. Alignment rates by MP parameter
    if 'overall_alignment_rate' in df.columns and 'mp' in df.columns:
        plt.figure(figsize=(10, 6))
        sns.boxplot(data=df, x='mp', y='overall_alignment_rate')
        plt.title('Overall Alignment Rate by MP Parameter')
        plt.xlabel('MP Parameter')
        plt.ylabel('Overall Alignment Rate (%)')
        plt.tight_layout()
        plt.savefig(output_path / 'alignment_rate_by_mp.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    # 3. Alignment rates by RDG parameter
    if 'overall_alignment_rate' in df.columns and 'rdg' in df.columns:
        plt.figure(figsize=(10, 6))
        sns.boxplot(data=df, x='rdg', y='overall_alignment_rate')
        plt.title('Overall Alignment Rate by RDG Parameter')
        plt.xlabel('RDG Parameter')
        plt.ylabel('Overall Alignment Rate (%)')
        plt.tight_layout()
        plt.savefig(output_path / 'alignment_rate_by_rdg.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    # 4. Alignment rates by RFG parameter
    if 'overall_alignment_rate' in df.columns and 'rfg' in df.columns:
        plt.figure(figsize=(10, 6))
        sns.boxplot(data=df, x='rfg', y='overall_alignment_rate')
        plt.title('Overall Alignment Rate by RFG Parameter')
        plt.xlabel('RFG Parameter')
        plt.ylabel('Overall Alignment Rate (%)')
        plt.tight_layout()
        plt.savefig(output_path / 'alignment_rate_by_rfg.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    # 5. Heatmap of alignment rates (Sensitivity vs MP)
    if all(col in df.columns for col in ['sensitivity', 'mp', 'overall_alignment_rate']):
        pivot_df = df.pivot_table(values='overall_alignment_rate', 
                                index='sensitivity', 
                                columns='mp', 
                                aggfunc='mean')
        
        plt.figure(figsize=(10, 8))
        sns.heatmap(pivot_df, annot=True, fmt='.2f', cmap='YlOrRd', 
                   cbar_kws={'label': 'Overall Alignment Rate (%)'})
        plt.title('Alignment Rate Heatmap: Sensitivity vs MP Parameter')
        plt.xlabel('MP Parameter')
        plt.ylabel('Sensitivity Setting')
        plt.tight_layout()
        plt.savefig(output_path / 'alignment_rate_heatmap_sensitivity_mp.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    # 6. Heatmap of alignment rates (RDG vs RFG)
    if all(col in df.columns for col in ['rdg', 'rfg', 'overall_alignment_rate']):
        pivot_df = df.pivot_table(values='overall_alignment_rate', 
                                index='rdg', 
                                columns='rfg', 
                                aggfunc='mean')
        
        plt.figure(figsize=(10, 8))
        sns.heatmap(pivot_df, annot=True, fmt='.2f', cmap='YlOrRd', 
                   cbar_kws={'label': 'Overall Alignment Rate (%)'})
        plt.title('Alignment Rate Heatmap: RDG vs RFG Parameter')
        plt.xlabel('RFG Parameter')
        plt.ylabel('RDG Parameter')
        plt.tight_layout()
        plt.savefig(output_path / 'alignment_rate_heatmap_rdg_rfg.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    # 7. Read length distributions comparison
    length_cols = [col for col in df.columns if 'mean_length' in col]
    if length_cols:
        fig, axes = plt.subplots(1, len(length_cols), figsize=(6*len(length_cols), 6))
        if len(length_cols) == 1:
            axes = [axes]
        
        for i, col in enumerate(length_cols):
            data_type = col.split('_')[0]
            sns.boxplot(data=df, x='sensitivity', y=col, ax=axes[i])
            axes[i].set_title(f'{data_type.title()} Read Mean Length by Sensitivity')
            axes[i].set_xlabel('Sensitivity Setting')
            axes[i].set_ylabel('Mean Length (bp)')
            axes[i].tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        plt.savefig(output_path / 'read_length_comparison.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    # 8. Parameter combination performance summary
    if 'overall_alignment_rate' in df.columns:
        # Sort by alignment rate
        df_sorted = df.sort_values('overall_alignment_rate', ascending=False)
        
        plt.figure(figsize=(15, 8))
        x_pos = range(len(df_sorted))
        plt.bar(x_pos, df_sorted['overall_alignment_rate'])
        plt.title('Overall Alignment Rate by Parameter Combination')
        plt.xlabel('Parameter Combination')
        plt.ylabel('Overall Alignment Rate (%)')
        plt.xticks(x_pos, df_sorted['parameter_combination'], rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(output_path / 'parameter_combination_performance.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    # 9. Alignment rates by parameters (concatenated plots)
    if all(col in df.columns for col in ['aligned_exactly_1_rate', 'aligned_multiple_rate']):
        # By sensitivity
        if 'sensitivity' in df.columns:
            fig, axes = plt.subplots(1, 2, figsize=(15, 6))
            
            sns.boxplot(data=df, x='sensitivity', y='aligned_exactly_1_rate', ax=axes[0])
            axes[0].set_title('Aligned Exactly 1 Rate by Sensitivity')
            axes[0].set_xlabel('Sensitivity Setting')
            axes[0].set_ylabel('Aligned Exactly 1 Rate (%)')
            axes[0].tick_params(axis='x', rotation=45)
            
            sns.boxplot(data=df, x='sensitivity', y='aligned_multiple_rate', ax=axes[1])
            axes[1].set_title('Aligned Multiple Rate by Sensitivity')
            axes[1].set_xlabel('Sensitivity Setting')
            axes[1].set_ylabel('Aligned Multiple Rate (%)')
            axes[1].tick_params(axis='x', rotation=45)
            
            plt.tight_layout()
            plt.savefig(output_path / 'alignment_rates_by_sensitivity.png', dpi=300, bbox_inches='tight')
            plt.close()
        
        # By RDG parameter
        if 'rdg' in df.columns:
            fig, axes = plt.subplots(1, 2, figsize=(12, 6))
            
            sns.boxplot(data=df, x='rdg', y='aligned_exactly_1_rate', ax=axes[0])
            axes[0].set_title('Aligned Exactly 1 Rate by RDG Parameter')
            axes[0].set_xlabel('RDG Parameter')
            axes[0].set_ylabel('Aligned Exactly 1 Rate (%)')
            
            sns.boxplot(data=df, x='rdg', y='aligned_multiple_rate', ax=axes[1])
            axes[1].set_title('Aligned Multiple Rate by RDG Parameter')
            axes[1].set_xlabel('RDG Parameter')
            axes[1].set_ylabel('Aligned Multiple Rate (%)')
            
            plt.tight_layout()
            plt.savefig(output_path / 'alignment_rates_by_rdg.png', dpi=300, bbox_inches='tight')
            plt.close()
        
        # By RFG parameter
        if 'rfg' in df.columns:
            fig, axes = plt.subplots(1, 2, figsize=(12, 6))
            
            sns.boxplot(data=df, x='rfg', y='aligned_exactly_1_rate', ax=axes[0])
            axes[0].set_title('Aligned Exactly 1 Rate by RFG Parameter')
            axes[0].set_xlabel('RFG Parameter')
            axes[0].set_ylabel('Aligned Exactly 1 Rate (%)')
            
            sns.boxplot(data=df, x='rfg', y='aligned_multiple_rate', ax=axes[1])
            axes[1].set_title('Aligned Multiple Rate by RFG Parameter')
            axes[1].set_xlabel('RFG Parameter')
            axes[1].set_ylabel('Aligned Multiple Rate (%)')
            
            plt.tight_layout()
            plt.savefig(output_path / 'alignment_rates_by_rfg.png', dpi=300, bbox_inches='tight')
            plt.close()
    # 10. Correlation matrix of numeric variables
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    if len(numeric_cols) > 1:
        plt.figure(figsize=(12, 10))
        correlation_matrix = df[numeric_cols].corr()
        sns.heatmap(correlation_matrix, annot=True, fmt='.2f', cmap='coolwarm', center=0,
                   square=True, linewidths=0.5, cbar_kws={"shrink": .8})
        plt.title('Correlation Matrix of Numeric Variables')
        plt.tight_layout()
        plt.savefig(output_path / 'correlation_matrix.png', dpi=300, bbox_inches='tight')
        plt.close()

def generate_summary_report(df, output_dir):
    """Generate a summary report of the results."""
    output_path = Path(output_dir)
    
    with open(output_path / 'summary_report.txt', 'w') as f:
        f.write("BOWTIE2 PARAMETER SWEEP ANALYSIS SUMMARY\n")
        f.write("=" * 50 + "\n\n")
        
        f.write(f"Total parameter combinations tested: {len(df)}\n\n")
        
        if 'overall_alignment_rate' in df.columns:
            f.write("ALIGNMENT RATE STATISTICS:\n")
            f.write(f"Mean alignment rate: {df['overall_alignment_rate'].mean():.2f}%\n")
            f.write(f"Median alignment rate: {df['overall_alignment_rate'].median():.2f}%\n")
            f.write(f"Best alignment rate: {df['overall_alignment_rate'].max():.2f}%\n")
            f.write(f"Worst alignment rate: {df['overall_alignment_rate'].min():.2f}%\n\n")
        
        if 'aligned_exactly_1_rate' in df.columns:
            f.write("ALIGNED EXACTLY 1 RATE STATISTICS:\n")
            f.write(f"Mean aligned exactly 1 rate: {df['aligned_exactly_1_rate'].mean():.2f}%\n")
            f.write(f"Median aligned exactly 1 rate: {df['aligned_exactly_1_rate'].median():.2f}%\n")
            f.write(f"Maximum aligned exactly 1 rate: {df['aligned_exactly_1_rate'].max():.2f}%\n")
            f.write(f"Minimum aligned exactly 1 rate: {df['aligned_exactly_1_rate'].min():.2f}%\n\n")
        
        if 'aligned_multiple_rate' in df.columns:
            f.write("ALIGNED MULTIPLE RATE STATISTICS:\n")
            f.write(f"Mean aligned multiple rate: {df['aligned_multiple_rate'].mean():.2f}%\n")
            f.write(f"Median aligned multiple rate: {df['aligned_multiple_rate'].median():.2f}%\n")
            f.write(f"Maximum aligned multiple rate: {df['aligned_multiple_rate'].max():.2f}%\n")
            f.write(f"Minimum aligned multiple rate: {df['aligned_multiple_rate'].min():.2f}%\n\n")
            
            # Best performing combination
            best_idx = df['overall_alignment_rate'].idxmax()
            best_combo = df.loc[best_idx]
            f.write("BEST PERFORMING COMBINATION:\n")
            f.write(f"Parameter combination: {best_combo['parameter_combination']}\n")
            f.write(f"Alignment rate: {best_combo['overall_alignment_rate']:.2f}%\n")
            if 'sensitivity' in best_combo:
                f.write(f"Sensitivity: {best_combo['sensitivity']}\n")
            if 'mp' in best_combo:
                f.write(f"MP: {best_combo['mp']}\n")
            if 'rdg' in best_combo:
                f.write(f"RDG: {best_combo['rdg']}\n")
            if 'rfg' in best_combo:
                f.write(f"RFG: {best_combo['rfg']}\n")
            f.write("\n")
        
        # Parameter-specific analysis
        if 'sensitivity' in df.columns and 'overall_alignment_rate' in df.columns:
            f.write("ALIGNMENT RATE BY SENSITIVITY:\n")
            sensitivity_stats = df.groupby('sensitivity')['overall_alignment_rate'].agg(['mean', 'std', 'count'])
            for sens, stats in sensitivity_stats.iterrows():
                f.write(f"{sens}: {stats['mean']:.2f}% ± {stats['std']:.2f}% (n={stats['count']})\n")
            f.write("\n")
        
        if 'mp' in df.columns and 'overall_alignment_rate' in df.columns:
            f.write("ALIGNMENT RATE BY MP PARAMETER:\n")
            mp_stats = df.groupby('mp')['overall_alignment_rate'].agg(['mean', 'std', 'count'])
            for mp, stats in mp_stats.iterrows():
                f.write(f"{mp}: {stats['mean']:.2f}% ± {stats['std']:.2f}% (n={stats['count']})\n")
            f.write("\n")
        
        if 'rdg' in df.columns and 'overall_alignment_rate' in df.columns:
            f.write("ALIGNMENT RATE BY RDG PARAMETER:\n")
            rdg_stats = df.groupby('rdg')['overall_alignment_rate'].agg(['mean', 'std', 'count'])
            for rdg, stats in rdg_stats.iterrows():
                f.write(f"{rdg}: {stats['mean']:.2f}% ± {stats['std']:.2f}% (n={stats['count']})\n")
            f.write("\n")
        
        if 'rfg' in df.columns and 'overall_alignment_rate' in df.columns:
            f.write("ALIGNMENT RATE BY RFG PARAMETER:\n")
            rfg_stats = df.groupby('rfg')['overall_alignment_rate'].agg(['mean', 'std', 'count'])
            for rfg, stats in rfg_stats.iterrows():
                f.write(f"{rfg}: {stats['mean']:.2f}% ± {stats['std']:.2f}% (n={stats['count']})\n")
            f.write("\n")
        
        # Parameter-specific analysis for alignment rates
        if 'aligned_exactly_1_rate' in df.columns:
            if 'sensitivity' in df.columns:
                f.write("ALIGNED EXACTLY 1 RATE BY SENSITIVITY:\n")
                sens_stats = df.groupby('sensitivity')['aligned_exactly_1_rate'].agg(['mean', 'std', 'count'])
                for sens, stats in sens_stats.iterrows():
                    f.write(f"{sens}: {stats['mean']:.2f}% ± {stats['std']:.2f}% (n={stats['count']})\n")
                f.write("\n")
            
            if 'rdg' in df.columns:
                f.write("ALIGNED EXACTLY 1 RATE BY RDG PARAMETER:\n")
                rdg_stats = df.groupby('rdg')['aligned_exactly_1_rate'].agg(['mean', 'std', 'count'])
                for rdg, stats in rdg_stats.iterrows():
                    f.write(f"{rdg}: {stats['mean']:.2f}% ± {stats['std']:.2f}% (n={stats['count']})\n")
                f.write("\n")
            
            if 'rfg' in df.columns:
                f.write("ALIGNED EXACTLY 1 RATE BY RFG PARAMETER:\n")
                rfg_stats = df.groupby('rfg')['aligned_exactly_1_rate'].agg(['mean', 'std', 'count'])
                for rfg, stats in rfg_stats.iterrows():
                    f.write(f"{rfg}: {stats['mean']:.2f}% ± {stats['std']:.2f}% (n={stats['count']})\n")
                f.write("\n")
        
        if 'aligned_multiple_rate' in df.columns:
            if 'sensitivity' in df.columns:
                f.write("ALIGNED MULTIPLE RATE BY SENSITIVITY:\n")
                sens_stats = df.groupby('sensitivity')['aligned_multiple_rate'].agg(['mean', 'std', 'count'])
                for sens, stats in sens_stats.iterrows():
                    f.write(f"{sens}: {stats['mean']:.2f}% ± {stats['std']:.2f}% (n={stats['count']})\n")
                f.write("\n")
            
            if 'rdg' in df.columns:
                f.write("ALIGNED MULTIPLE RATE BY RDG PARAMETER:\n")
                rdg_stats = df.groupby('rdg')['aligned_multiple_rate'].agg(['mean', 'std', 'count'])
                for rdg, stats in rdg_stats.iterrows():
                    f.write(f"{rdg}: {stats['mean']:.2f}% ± {stats['std']:.2f}% (n={stats['count']})\n")
                f.write("\n")
            
            if 'rfg' in df.columns:
                f.write("ALIGNED MULTIPLE RATE BY RFG PARAMETER:\n")
                rfg_stats = df.groupby('rfg')['aligned_multiple_rate'].agg(['mean', 'std', 'count'])
                for rfg, stats in rfg_stats.iterrows():
                    f.write(f"{rfg}: {stats['mean']:.2f}% ± {stats['std']:.2f}% (n={stats['count']})\n")
                f.write("\n")

def main():
    """Main function to run the analysis."""
    # Configuration
    base_dir = "/home/users/astar/gis/linxy/scratch/data/singlecell_shapemap/RHH_only/DMSO/mapping/parameter_sweep"
    output_dir = "/home/users/astar/gis/linxy/scratch/data/singlecell_shapemap/RHH_only/DMSO/mapping/parameter_sweep_analysis"
    
    print("Collecting results from parameter sweep...")
    df = collect_all_results(base_dir)
    
    if df.empty:
        print("No data found. Please check the base directory path.")
        return
    
    print(f"Found {len(df)} parameter combinations")
    print(f"Columns available: {list(df.columns)}")
    
    # Save combined data
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path / 'combined_results.csv', index=False)
    print(f"Combined data saved to {output_path / 'combined_results.csv'}")
    
    # Create visualizations
    print("Creating visualizations...")
    create_visualizations(df, output_dir)
    
    # Generate summary report
    print("Generating summary report...")
    generate_summary_report(df, output_dir)
    
    print(f"Analysis complete! Results saved in: {output_dir}")
    print(f"Check the following files:")
    print(f"  - combined_results.csv: Raw data")
    print(f"  - summary_report.txt: Summary statistics")
    print(f"  - *.png: Visualization plots")

if __name__ == "__main__":
    main()