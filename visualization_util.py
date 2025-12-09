## This file contains utility functions for visualizing mutation rates and reactivity data.

def normalize_reactivity(reactivity, percentile=95):
    """
    Normalize reactivity values using a specified percentile.
    
    Parameters:
    -----------
    reactivity : array-like
        Array of reactivity values to normalize
    percentile : float, default 90
        Percentile value to use for normalization
    
    Returns:
    --------
    numpy.ndarray
        Normalized reactivity values
    """
    import numpy as np
    
    # Calculate the specified percentile value
    P05 = np.percentile(reactivity, 5)
    P95 = np.percentile(reactivity, 95)
    
    # Avoid division by zero
    if P95 - P05 == 0:
        print("Warning: P95 equals P05, returning original reactivity values.")
        return reactivity
    
    # Normalize reactivity values
    normalized = np.zeros_like(reactivity)
    for i, val in enumerate(reactivity):
        if val < P05:
            normalized[i] = 0
        elif val > P95:
            normalized[i] = 1
        else:
            normalized[i] = (val - P05) / (P95 - P05)
    
    return normalized

def read_mutrate_window_files(directory_path='.', pattern = "_counts_mutrate.10nt.csv.gz", prefix ='', output_file=None, coverage_threshold = 200):
    """
    Process all windowed mutrate files in the specified directory.
    
    Parameters:
    -----------
    directory_path : str, default '.'
        Path to the directory containing the mutrate files
    pattern : str, default '10nt.csv.gz'
        Pattern to match windowed mutrate files
    output_file : str, optional
        If provided, save the concatenated results to this file
    prefix : str, default ''
        Prefix to add to sample names (e.g., 'C' for clarity)
    
    Returns:
    --------
    pandas.DataFrame
        Concatenated dataframe with columns: sample, gene.position, mutrate, coverage, mutant
    """
    import os
    import glob
    import pandas as pd
    
    # Find all files ending with _mutrate.txt
    path_pattern = os.path.join(directory_path, f'*{pattern}')
    mutrate_files = glob.glob(path_pattern)
    
    if not mutrate_files:
        print(f"No files ending with '_mutrate.txt' found in {directory_path}")
        return pd.DataFrame()
    
    print(f"Found {len(mutrate_files)} mutrate files:")
    for file in mutrate_files:
        print(f"  - {os.path.basename(file)}")
    
    all_data = []
    
    for file_path in mutrate_files:
        try:
            # Extract sample name from filename (remove pattern part)
            filename = os.path.basename(file_path)
            sample_name = filename.replace(f"{pattern}", '')

            # Extra step: Add C in sample name for clarity (optional)
            sample_name = f"{prefix}{sample_name}"
            
            print(f"Processing {filename} (sample: {sample_name})...")
            
            # Read the file
            df = pd.read_csv(file_path, sep='\t', compression='gzip')
            
            df['sample'] = sample_name
            
            all_data.append(df)
            print(f"  - Added {len(df)} rows from {sample_name}")
            
        except Exception as e:
            print(f"Error processing {file_path}: {str(e)}")
            continue
    
    if not all_data:
        print("No data was successfully processed.")
        return pd.DataFrame()
    
    # Concatenate all dataframes
    concat_df = pd.concat(all_data, ignore_index=True)

    # Filter dataframes based on coverage
    final_df = concat_df[concat_df['coverage'] >= coverage_threshold]
    
    print(f"\nFinal concatenated dataframe:")
    print(f"  - Total rows: {len(final_df)}")
    print(f"  - Samples: {sorted(final_df['sample'].unique())}")
    print(f"  - Columns: {list(final_df.columns)}")
    
    # Save to file if requested
    if output_file:
        final_df.to_csv(output_file, sep='\t', index=False)
        print(f"  - Saved to: {output_file}")
    
    return final_df

def calculate_reactivity(signal_df, control_df, index_cols = ['gene', 'pos']):
    """
    Calculate reactivity as the difference between signal and control mutational rates.
    Require mutrate and coverage columns in both dataframes as well as the index columns.
    
    Parameters:
    -----------
    signal_df : pandas.DataFrame
        DataFrame containing mutational rates for the signal condition (e.g., treated)
    control_df : pandas.DataFrame
        DataFrame containing mutational rates for the control condition (e.g., untreated)
    index_cols : list of str, default ['gene', 'pos']
        Columns to use as index for merging the two DataFrames
    
    Returns:
    --------
    pandas.DataFrame
        DataFrame with columns: gene, pos, signal_mutrate, control_mutrate, reactivity
    """
    import pandas as pd

    # Merge the two DataFrames on the specified index columns
    merged_df = pd.merge(signal_df, control_df, on=index_cols, suffixes=('_signal', '_control'))

    # Give warning if there is NA values in mutrate columns after merging
    signal_mutrate_na = merged_df['mutrate_signal'].isna().sum()
    control_mutrate_na = merged_df['mutrate_control'].isna().sum()
    if signal_mutrate_na > 0 or control_mutrate_na > 0:
        print(f"Warning: There are {signal_mutrate_na} NA values in signal mutrate and {control_mutrate_na} NA values in control mutrate after merging.")

    # Exclude rows with NA values in mutrate columns
    merged_df = merged_df.dropna(subset=['mutrate_signal', 'mutrate_control'])
    
    # Calculate reactivity
    merged_df['reactivity'] = merged_df['mutrate_signal'] - merged_df['mutrate_control']
    
    # Select relevant columns to return
    result_df = merged_df[index_cols + ['mutrate_signal', 'mutrate_control', 'reactivity', 'coverage_signal', 'coverage_control']]
    
    return result_df

def gene_normalization_1nt(crude_reac, bottom_winsorize=0.0, top_winsorize=0.0, 
                          transcriptome_winsorize=False, gene_winsorize=False):
    ## transcriptome-wide winsorization
    from scipy.stats import stats
    from scipy.stats.mstats import winsorize
    import pandas as pd
    import numpy as np

    crude_reac_array = crude_reac.to_numpy()
    if transcriptome_winsorize:
        crude_reac_array = winsorize(crude_reac_array, limits=[bottom_winsorize, top_winsorize], 
                                   axis=0, nan_policy="omit")
    winsorized_reac = pd.DataFrame(crude_reac_array, index=crude_reac.index, columns=crude_reac.columns)
    
    concat_genes = []
    for gene, gene_crude_reac in winsorized_reac.groupby(level="gene"):
        ## gene winsorization
        crude_reac_array = gene_crude_reac.to_numpy()
        if gene_winsorize:
            crude_reac_array = winsorize(crude_reac_array, limits=[bottom_winsorize, top_winsorize], 
                                       axis=0, nan_policy="omit")
        gene_winsorized_reac = pd.DataFrame(crude_reac_array, index=gene_crude_reac.index, 
                                          columns=gene_crude_reac.columns)
        
        concat_stages = []
        for stage, stage_crude_reac in gene_winsorized_reac.items():
            
            iqr = stats.iqr(stage_crude_reac)
            if stage_crude_reac.index.size >= 500:
                pct95_val = np.nanpercentile(stage_crude_reac, 99)
            else:
                pct95_val = np.nanpercentile(stage_crude_reac, 95)
            filter_outlies_series = stage_crude_reac.loc[stage_crude_reac < max(pct95_val, 1.5*iqr)]
            pct90_val_filtered = np.nanpercentile(filter_outlies_series, 90)
            avg_pct90_filtered_vals = filter_outlies_series.loc[filter_outlies_series>pct90_val_filtered].mean()
            normalized_reac = stage_crude_reac/avg_pct90_filtered_vals
            concat_stages.append(normalized_reac)
        gene_norm_reac = pd.concat(concat_stages, axis=1)
        concat_genes.append(gene_norm_reac)
    gene_norm_reac = pd.concat(concat_genes, axis=0)
    gene_norm_reac.columns.names = crude_reac.columns.names
    return gene_norm_reac

def normalize_reactivity_zy(df, reactivity_columns=None, bottom_winsorize=0.0, top_winsorize=0.0,
                              transcriptome_winsorize=False, gene_winsorize=False, index_col='win'):
    """
    Wrapper function for gene_normalization_1nt that prepares the input dataframe.

    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe with 'gene' as a column and reactivity measurements
    reactivity_columns : list, optional
        List of column names containing reactivity data. If None, will auto-detect
        numeric columns excluding index columns
    bottom_winsorize : float, default 0.0
        Bottom percentile for winsorization
    top_winsorize : float, default 0.0
        Top percentile for winsorization
    transcriptome_winsorize : bool, default False
        Whether to apply transcriptome-wide winsorization
    gene_winsorize : bool, default False
        Whether to apply gene-level winsorization
    index_col : str or list, default 'win'
        Column(s) to use as index along with 'gene'. Can be a single column name
        or a list of column names for multi-level indexing.

    Returns:
    --------
    pd.DataFrame
        DataFrame with original columns plus normalized reactivity columns
    """
    import pandas as pd
    import numpy as np

    # Make a copy to avoid modifying the original dataframe
    df_work = df.copy()

    # Handle index_col as list or string
    if isinstance(index_col, str):
        index_cols = [index_col]
    else:
        index_cols = list(index_col)

    # Auto-detect reactivity columns if not specified
    if reactivity_columns is None:
        # Get numeric columns, excluding index columns
        numeric_cols = df_work.select_dtypes(include=[np.number]).columns.tolist()
        # Remove all index columns from numeric columns
        for col in index_cols:
            if col in numeric_cols:
                numeric_cols.remove(col)
        reactivity_columns = numeric_cols
        print(f"Auto-detected reactivity columns: {reactivity_columns}")

    # Ensure reactivity_columns is always a list
    if isinstance(reactivity_columns, str):
        reactivity_columns = [reactivity_columns]

    # Set multi-level index with gene and all index columns
    index_to_set = ['gene'] + index_cols
    crude_reac = df_work.set_index(index_to_set)[reactivity_columns]
    
    # Ensure we have a DataFrame, not a Series
    if isinstance(crude_reac, pd.Series):
        crude_reac = crude_reac.to_frame()
    
    # Call the original normalization function
    normalized_result = gene_normalization_1nt(
        crude_reac=crude_reac,
        bottom_winsorize=bottom_winsorize,
        top_winsorize=top_winsorize,
        transcriptome_winsorize=transcriptome_winsorize,
        gene_winsorize=gene_winsorize
    )
    
    # Reset index to match original dataframe structure
    normalized_result = normalized_result.reset_index()

    # Start with the complete original dataframe
    result = df.copy()

    # Create mapping for normalized column names
    if len(reactivity_columns) == 1:
        # If single column, rename to 'normalized_reactivity'
        normalized_columns_map = {reactivity_columns[0]: 'normalized_reactivity'}
    else:
        # If multiple columns, add '_normalized' suffix
        normalized_columns_map = {col: f"{col}_normalized" for col in reactivity_columns}

    # Rename normalized columns
    normalized_result = normalized_result.rename(columns=normalized_columns_map)

    # Get the new normalized column names
    new_normalized_cols = list(normalized_columns_map.values())

    # Merge the normalized columns back using all index columns
    merge_cols = ['gene'] + index_cols
    normalized_data = normalized_result[merge_cols + new_normalized_cols]
    result = pd.merge(result, normalized_data, on=merge_cols, how='left')

    return result

def get_df_with_normalized_reactivity(signal_df, control_df, index_cols=['gene', 'pos'],
                                     sample_col="group", percentile=95, normalization_method="normalize_reactivity_zy"):
    """
    Add a column of normalized reactivity values to the input DataFrame.
    
    Parameters:
    -----------
    signal_df : pandas.DataFrame
        DataFrame containing mutational rates for the signal condition (e.g., treated)
    control_df : pandas.DataFrame
        DataFrame containing mutational rates for the control condition (e.g., untreated)
    index_cols : list of str, default ['gene', 'pos', 'sample']
        Columns to use as index for merging the two DataFrames
    sample_col : str, default 'sample'
        Column name representing different samples or conditions
    percentile : float, default 95
        Percentile value to use for normalization
    normalization_method : str, default "normalize_reactivity_zy"
        Choose normalization method: "normalize_reactivity_zy" or "normalize_reactivity"
    
    Returns:
    --------
    pandas.DataFrame
        DataFrame with an additional column 'normalized_reactivity'
    """
    import pandas as pd

    # Validate normalization method
    valid_methods = ["normalize_reactivity_zy", "normalize_reactivity"]
    if normalization_method not in valid_methods:
        raise ValueError(f"normalization_method must be one of {valid_methods}")
    
    # Calculate reactivity
    print("Calculating reactivity...")
    reactivity_df = calculate_reactivity(signal_df, control_df, index_cols=index_cols + [sample_col])
    print("The number of rows in reactivity_df:", len(reactivity_df))
    
    # Normalize reactivity values based on chosen method
    if normalization_method == "normalize_reactivity":
        print("Normalizing reactivity using normalize_reactivity function...")
        normalized_reactivity = reactivity_df.groupby(sample_col)['reactivity'].transform(
            lambda x: normalize_reactivity(x, percentile=percentile)
        )
        reactivity_df['normalized_reactivity'] = normalized_reactivity
        
    elif normalization_method == "normalize_reactivity_zy":
        print("Normalizing reactivity using normalize_reactivity_zy function group by group...")
        normalized_groups = []
        for group_name, group_data in reactivity_df.groupby(sample_col):
            print(f"Processing group: {group_name}")
            # Apply normalize_reactivity_zy to this group
            normalized_group = normalize_reactivity_zy(
                df=group_data, 
                reactivity_columns="reactivity",
                index_col = [col for col in index_cols if col != "gene"]
            )
            normalized_groups.append(normalized_group)
        
        # Concatenate all normalized groups back together
        reactivity_df = pd.concat(normalized_groups, ignore_index=True)
    
    return reactivity_df

def plot_correlation(dataframe, index_col, group_col, value_col, output_file, title_prefix=''):
    """
    Generate a correlation heatmap with hierarchical clustering.
    
    Parameters:
    -----------
    dataframe : pandas.DataFrame
        Input dataframe with data to correlate
    index_col : str
        Column to use as index (e.g., 'gene_position')
    group_col : str
        Column to use for grouping/samples (e.g., 'sample')
    value_col : str
        Column containing values to correlate (e.g., 'mutrate')
    output_file : str
        Path to save the correlation plot
    """
    import seaborn as sns
    from matplotlib import pyplot as plt
    import pandas as pd

    df_wide = dataframe.pivot_table(index=index_col, columns=group_col, values=value_col)
    df_wide = df_wide.dropna()
    
    if df_wide.empty:
        print("Error: No data remaining after pivot and dropna")
        return
    
    print(f"Pivot table shape: {df_wide.shape}")
    print(f"Samples in correlation: {list(df_wide.columns)}")
    
    corr_df = df_wide.corr(method='pearson')
    
    # Create clustermap
    g = sns.clustermap(corr_df, cmap="coolwarm", annot=True, method='ward', 
                       figsize=(5, 5), fmt='.2f', linewidths=0.5,
                       cbar_kws={'label': f'Pearson Correlation'})
    
    # Adjust title position - include number of data points
    n_datapoints = df_wide.shape[0]
    g.fig.suptitle(f'Pearson Correlation of {value_col} (n={n_datapoints}) {title_prefix}', y=1.02)
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Correlation plot saved to: {output_file}")
    
    plt.show()

def read_mutrate_files(directory_path='.', prefix ='', output_file=None, header=None, coverage_threshold=100):
    """
    Process all mutrate.txt and mutrate.txt.gz files in the specified directory.

    Parameters:
    -----------
    directory_path : str, default '.'
        Path to the directory containing the mutrate files
    output_file : str, optional
        If provided, save the concatenated results to this file
    prefix : str, default ''
        Prefix to add to sample names (e.g., 'C' for clarity)
    header : list, optional
        Column names to use if the file has no header. Default is:
        ["gene", "pos", "pos+1", "gene.position", "mutrate", "strand", "coverage",
         "coverage_withIndel", "mutant", "normalized_cov", "g_readcount", "refnt", "detail"]
    coverage_threshold : int, default 100
        Minimum coverage required. Rows with coverage below this threshold will be filtered out.

    Returns:
    --------
    pandas.DataFrame
        Concatenated dataframe with columns: sample, gene.position, mutrate, coverage, mutant
    """
    import os
    import glob
    import pandas as pd
    import gzip

    # Set default header if not provided
    if header is None:
        header = ["gene", "pos", "pos+1", "gene.position", "mutrate", "strand", "coverage",
                  "coverage_withIndel", "mutant", "normalized_cov", "g_readcount", "refnt", "detail"]

    # Find all files ending with mutrate.txt or mutrate.txt.gz
    pattern_txt = os.path.join(directory_path, '*mutrate.txt')
    pattern_gz = os.path.join(directory_path, '*mutrate.txt.gz')
    mutrate_files = glob.glob(pattern_txt) + glob.glob(pattern_gz)
    
    if not mutrate_files:
        print(f"No files ending with 'mutrate.txt' or 'mutrate.txt.gz' found in {directory_path}")
        return pd.DataFrame()

    print(f"Found {len(mutrate_files)} mutrate files:")
    for file in mutrate_files:
        print(f"  - {os.path.basename(file)}")

    all_data = []

    for file_path in mutrate_files:
        try:
            # Extract sample name from filename
            filename = os.path.basename(file_path)
            if filename.endswith('.gz'):
                sample_name = filename[:-len('mutrate.txt.gz')]
            else:
                sample_name = filename[:-len('mutrate.txt')]

            # Remove trailing _ or . if present
            while sample_name and sample_name[-1] in ('_', '.'):
                sample_name = sample_name[:-1]

            # Extra step: Add prefix to sample name for clarity (optional)
            sample_name = f"{prefix}{sample_name}"

            print(f"Processing {filename} (sample: {sample_name})...")

            # Read the file (handle both .txt and .txt.gz)
            # First, check if the file has a header by examining the first line
            if file_path.endswith('.gz'):
                with gzip.open(file_path, 'rt') as f:
                    first_line = f.readline().strip()
                    f.seek(0)  # Reset to beginning

                    # Check if first line looks like a header (contains 'gene' or similar)
                    if 'gene' in first_line.lower() or 'pos' in first_line.lower():
                        df = pd.read_csv(f, sep='\t')
                    else:
                        # No header, use provided header
                        f.seek(0)  # Reset again
                        df = pd.read_csv(f, sep='\t', header=None, names=header)
            else:
                with open(file_path, 'r') as f:
                    first_line = f.readline().strip()

                # Check if first line looks like a header
                if 'gene' in first_line.lower() or 'pos' in first_line.lower():
                    df = pd.read_csv(file_path, sep='\t')
                else:
                    # No header, use provided header
                    df = pd.read_csv(file_path, sep='\t', header=None, names=header)

            # Check if required columns exist
            required_cols = ['gene','pos', 'mutrate', 'coverage', 'mutant']
            missing_cols = [col for col in required_cols if col not in df.columns]
            
            if missing_cols:
                print(f"Warning: Missing columns in {filename}: {missing_cols}")
                continue
            
            # Extract required columns and add sample column
            subset_df = df[required_cols].copy()
            subset_df['sample'] = sample_name

            # Filter by coverage threshold
            rows_before = len(subset_df)
            subset_df = subset_df[subset_df['coverage'] >= coverage_threshold]
            rows_after = len(subset_df)

            # Reorder columns to put sample first
            subset_df = subset_df[['sample', 'gene', 'pos', 'mutrate', 'coverage', 'mutant']]

            all_data.append(subset_df)
            print(f"  - Added {rows_after} rows from {sample_name} (filtered {rows_before - rows_after} rows with coverage < {coverage_threshold})")
            
        except Exception as e:
            print(f"Error processing {file_path}: {str(e)}")
            continue
    
    if not all_data:
        print("No data was successfully processed.")
        return pd.DataFrame()
    
    # Concatenate all dataframes
    final_df = pd.concat(all_data, ignore_index=True)
    
    print(f"\nFinal concatenated dataframe:")
    print(f"  - Total rows: {len(final_df)}")
    print(f"  - Samples: {sorted(final_df['sample'].unique())}")
    print(f"  - Columns: {list(final_df.columns)}")
    
    # Save to file if requested
    if output_file:
        final_df.to_csv(output_file, sep='\t', index=False)
        print(f"  - Saved to: {output_file}")
    
    return final_df