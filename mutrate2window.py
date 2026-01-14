#!/home/users/astar/gis/gaog/.conda/envs/main/bin/python

## This script is getting from Gao Ge. 
## Version: v2.0 
## Change the input file from gz to txt, and it contains header.
## Order dataframe by gene and pos before processing.
## Add test() function to test the script.

## Important: Might need to add the step for filtering based on mutation rates. In this case, need to consider how to make sure the position is the same between control and signal files.


import argparse
import pandas as pd
import numpy as np
import os
import sys

__doc__="""
    The script is used for getting fixed windows from *.mutrate.txt.gz.
    Usage: python mutrate2window.py -i mutrate.txt.gz -o win.gz -w 10 -s 10 --coverage-threshold 100
    Criteria:
    - Input: mutrate.txt.gz
    - Output: 10nt.csv.gz
    """
__version__="v2.0"
__author__="linxy"
__last_modify__="Dec-2025"

def gene_rolling(df, win_len, step, coverage_threshold, mutrate_threshold=0.2):
    '''
    intro: Create sliding windows based on actual genomic positions (not consecutive indices)
    
    Parameters: df: pd.DataFrame
                    DataFrame of single gene
                win_len: int
                    Length of window in base pairs
                step: int, default None
                    Step size for sliding window in base pairs
                threshold: cut off wins whose coverage less than the threshold
    
    return: DataFrame of single gene with win per row
                win_start: start position of the window
                win_end: end position of the window
                coverage: sum of bases coverage in the window
                mutant: sum of bases mutation in the window
    '''
    
    # step1: get sub df and sort by position
    df = df[['gene','pos','coverage','mutant','mutrate']].copy()
    df = df.sort_values(by=['pos'])
    
    # Get the range of positions
    min_pos = df['pos'].min()
    max_pos = df['pos'].max()
    
    # Create windows based on actual genomic positions
    windows = []
    win_start = min_pos
    
    while win_start <= max_pos:
        win_end = win_start + win_len - 1
        
        # Get all positions within this window
        window_data = df[(df['pos'] >= win_start) & (df['pos'] <= win_end) & (df['mutrate'] <= mutrate_threshold)]
        
        if len(window_data) > 0:  # Only process windows with data
            coverage_sum = window_data['coverage'].sum()
            mutant_sum = window_data['mutant'].sum()
            
            # Apply coverage threshold
            if coverage_sum >= coverage_threshold:
                mutrate = round(mutant_sum / coverage_sum, 5) if coverage_sum > 0 else 0
                
                windows.append({
                    'win': win_start,
                    'win_end': win_end,
                    'coverage': coverage_sum,
                    'mutant': mutant_sum,
                    'mutrate': mutrate
                })
        
        win_start += step
    
    # Convert to DataFrame
    result_df = pd.DataFrame(windows)
    
    # Add window number for compatibility
    if len(result_df) > 0:
        result_df = result_df[['win', 'coverage', 'mutant', 'mutrate', 'win_end']]
    else:
        # Return empty DataFrame with correct columns
        result_df = pd.DataFrame(columns=['win', 'coverage', 'mutant', 'mutrate', 'win_end'])
    
    return result_df

def gene_rolling_withpad(df, win_len, step, threshold):
    '''
    intro: optimize df.rolling and specify it, and cut off sparse windows using given coverage threshold
    
    Parameters: df: pd.DataFrame
                    DataFrame of single gene
                win_len: int
                    Length of window
                step: int, default None
                    Evaluate the window at every step result
                threshold: cut off wins whose coverage less than the threshold
    
    return: DataFrame of single gene with win per row
                coverage: sum of bases coverage in the window
                mutant: sum of bases mutation in the window
    '''
    
    # step1: get sub df
    df = df[['gene','pos','coverage','mutant']]
    df = df.sort_values(by=['gene','pos'])
    df.set_index('pos', inplace=True, drop=True)
    # gene = df.iloc[0,0]  # gene id
    
    # step2: padding at both ends(optimize df.rolling()):
    # add a same row as the first row before the 1st row
    new_row = df.iloc[0].to_frame().T
    df = pd.concat([new_row, df]).reset_index(drop=True)  # the first row is fake row, pos still is 1-base!
    # add zero rows at the tail so that the fade-pos end is integer multiple of win-size
    n = win_len - (len(df) - 1 ) % win_len  # rows num added
    # global df_na
    df_na_concat = pd.concat([df_na] * n)
    df = pd.concat([df, df_na_concat]).reset_index(drop=True)
    
    # step3: rolling
    # gene = df.iloc[0,0]
    df = df[['coverage', 'mutant']].rolling(window=win_len, step=step, min_periods=0).sum().iloc[1:]  # without gene name
    # df['gene'] = gene  # get varible gene in step1
    df.reset_index(inplace=True)
    df.rename(columns={'index':'win'}, inplace=True)
    
    # step4: cut off
    df = df[df['coverage']>=threshold]
    
    # step5: get mutrate
    df['mutrate'] = round(df['mutant'] / df['coverage'], 5)
    # return df[['gene', 'win', 'coverage', 'mutant', 'mutrate']]
    return df[['win', 'coverage', 'mutant', 'mutrate']]

# def main(win_len, step, threshold, input_file, output_file):
def main(args):
    '''
    mutrate: win_mutant / win_cov
    '''

    # Input validation
    if not os.path.exists(args.input_file):
        print(f"[ERROR] Input file not found: {args.input_file}")
        sys.exit(1)

    # Print parameters
    print(f"[INFO] Input file: {args.input_file}")
    print(f"[INFO] Output file: {args.output_file}")
    print(f"[INFO] Window length: {args.win_len} bp")
    print(f"[INFO] Step size: {args.step} bp")
    print(f"[INFO] Coverage threshold: {args.coverage_threshold}")
    print(f"[INFO] Mutrate threshold: {args.mutrate_threshold}")
    print(f"[INFO] Reading input file...")

    col_names = ["gene", "pos", "end", "gene.position", "mutrate", "strand", "coverage", "coverage_withIndel", "mutant", "normalized_cov", "g_readcount", "refnt", "detail"]
    #mtr_df = pd.read_csv(args.input_file, sep="\t", names=col_names, index_col = "gene.position", compression='gzip') # Original Gao Ge's version
    try:
        mtr_df = pd.read_csv(args.input_file, sep="\t", names=col_names, header=None, index_col="gene.position")
    except Exception as e:
        print(f"[ERROR] Failed to read input file: {e}")
        sys.exit(1)

    print(f"[INFO] Loaded {len(mtr_df)} positions from input file")

    grouped = mtr_df.groupby('gene')
    print(f"[INFO] Found {len(grouped)} genes to process")
    dfs = []  # gene : wins
    
    #global df_na
    #df_na = pd.DataFrame({'gene': None, 'coverage': 0, 'mutant': 0}, index=[0])  # will be used in gene_rolling_withpad()
    
    for gene, group in grouped:
        group.reset_index(inplace=True)
        win = gene_rolling(group, args.win_len, args.step, args.coverage_threshold, args.mutrate_threshold)
        if len(win) > 0:  # Only add if there are windows
            win.insert(0, 'gene', gene) 
            dfs.append(win)
    
    if dfs:  # Only concatenate if there are DataFrames
        mtr_win_df = pd.concat(dfs, axis=0)
        mtr_win_df.reset_index(drop=True, inplace=True)
    else:
        # Create empty DataFrame with correct columns
        mtr_win_df = pd.DataFrame(columns=['gene', 'win', 'coverage', 'mutant', 'mutrate', 'win_end'])
    
    # write the results to the output file as a tab-separated file
    mtr_win_df.to_csv(args.output_file, sep='\t', index=False, compression='gzip')
    
    print(f"Processed {len(mtr_win_df)} windows across {len(dfs)} genes")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i","--input", dest='input_file', help="inputfile")
    parser.add_argument("-o","--output", dest='output_file', default="-", help="outputfile")
    parser.add_argument("-w","--window_length", dest='win_len', type=int, default=10,help="the length of window, default 10")
    parser.add_argument("-s","--step",dest='step',default=10,type=int,help="evaluate the window at every step result, default 10")
    parser.add_argument("--coverage-threshold",dest='coverage_threshold',default=200,type=int,help="threshold of coverage sum of one window, default 200")
    parser.add_argument("--mutrate-threshold",dest='mutrate_threshold',default=0.2,type=float,help="threshold of mutrate for one position, default 0.2")
    # parser.add_argument("-","--", dest='',default="",help="")

    args = parser.parse_args()
    main(args)
    print("Windows' mutation rates calculation completed successfully.")
    # Uncomment below for testing
    #main(argparse.Namespace(input_file="/home/users/astar/gis/linxy/scratch/data/spatial_shapemap/4th_Round/2A3/basecount_per_cluster/1.0_counts_mutrate.txt", output_file="/home/users/astar/gis/linxy/scratch/data/spatial_shapemap/4th_Round/2A3/basecount_per_cluster/1.0_counts_mutrate.10nt.csv.gz", win_len=10, step=10, threshold=100))