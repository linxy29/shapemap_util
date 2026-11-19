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

def calculate_reactivity(signal_df, control_df, index_cols = ['gene', 'pos'], keep_cols=[],
                         trim_5=0, trim_3=0, pos_col='pos', gene_col='gene'):
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
    keep_cols : list of str, default []
        Extra columns to carry through into the returned DataFrame.
    trim_5 : int, default 0
        Drop positions within this many units of the 5' end of each transcript
        (i.e. pos < min_pos(gene) + trim_5). Units are whatever pos_col uses
        (nucleotides for pos_col='pos', windows for pos_col='win'). The 5' / 3'
        ends of a transcript often carry artifactual reactivity from RT priming
        bias and library-prep edge effects; trimming masks them before any
        downstream normalization.
    trim_3 : int, default 0
        Drop positions within this many units of the 3' end of each transcript
        (i.e. pos > max_pos(gene) - trim_3). Same units as trim_5.
    pos_col : str, default 'pos'
        Name of the position column used for edge trimming.
    gene_col : str, default 'gene'
        Name of the gene/transcript column used to scope edge trimming.

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

    # Trim 5' / 3' edge positions per transcript using observed min/max in pos_col
    if trim_5 > 0 or trim_3 > 0:
        if gene_col not in merged_df.columns or pos_col not in merged_df.columns:
            raise ValueError(
                f"trim_5/trim_3 require '{gene_col}' and '{pos_col}' columns in the merged DataFrame; "
                f"got columns {list(merged_df.columns)}"
            )
        bounds = merged_df.groupby(gene_col)[pos_col].agg(_trim_min='min', _trim_max='max')
        merged_df = merged_df.merge(bounds, left_on=gene_col, right_index=True)
        keep_mask = (merged_df[pos_col] >= merged_df['_trim_min'] + trim_5) & \
                    (merged_df[pos_col] <= merged_df['_trim_max'] - trim_3)
        n_dropped = (~keep_mask).sum()
        if n_dropped:
            print(f"Edge trim (5'={trim_5}, 3'={trim_3}) dropped {n_dropped} positions across {bounds.shape[0]} {gene_col}(s).")
        merged_df = merged_df[keep_mask].drop(columns=['_trim_min', '_trim_max'])

    # Resolve keep_cols that may have been suffixed by the merge (present in both inputs).
    # Prefer the unsuffixed column if it survived the merge, otherwise fall back to the
    # _signal copy (assumed equal to _control for per-sample attributes like 'group').
    resolved_keep_cols = []
    for col in keep_cols:
        if col in merged_df.columns:
            resolved_keep_cols.append(col)
        elif f"{col}_signal" in merged_df.columns:
            merged_df = merged_df.rename(columns={f"{col}_signal": col})
            resolved_keep_cols.append(col)
        else:
            raise KeyError(
                f"keep_cols entry '{col}' not found in merged DataFrame "
                f"(also looked for '{col}_signal'). Available columns: {list(merged_df.columns)}"
            )

    # Select relevant columns to return
    result_df = merged_df[index_cols + ['mutrate_signal', 'mutrate_control', 'reactivity', 'coverage_signal', 'coverage_control'] + resolved_keep_cols]

    return result_df

def gene_normalization_1nt(crude_reac, bottom_winsorize=0.0, top_winsorize=0.0,
                          transcriptome_winsorize=False, gene_winsorize=False):
    """
    Normalize reactivity values at 1nt resolution using gene-level normalization.

    Parameters:
    -----------
    crude_reac : pd.DataFrame
        DataFrame with multi-index (gene, position) and reactivity columns
    bottom_winsorize : float, default 0.0
        Bottom percentile for winsorization (0.0 to 1.0)
    top_winsorize : float, default 0.0
        Top percentile for winsorization (0.0 to 1.0)
    transcriptome_winsorize : bool, default False
        Whether to apply transcriptome-wide winsorization
    gene_winsorize : bool, default False
        Whether to apply gene-level winsorization

    Returns:
    --------
    pd.DataFrame
        DataFrame with normalized reactivity values
    """
    import numpy as np
    import pandas as pd
    from scipy.stats import stats
    from scipy.stats.mstats import winsorize

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

def get_df_with_normalized_reactivity(signal_df, control_df, index_cols=['gene', 'pos'], keep_cols=[],
                                     sample_col="group", percentile=95, normalization_method="normalize_reactivity_zy",
                                     trim_5=0, trim_3=0, pos_col='pos', gene_col='gene'):
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
    trim_5, trim_3 : int, default 0
        5' / 3' edge positions to drop per transcript before normalization. See
        calculate_reactivity for semantics. Units follow pos_col.
    pos_col : str, default 'pos'
        Position column used for edge trimming (e.g. 'pos' for nt, 'win' for windows).
    gene_col : str, default 'gene'
        Gene/transcript column used to scope edge trimming.

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

    # Calculate reactivity (with optional 5'/3' edge trimming)
    print("Calculating reactivity...")
    reactivity_df = calculate_reactivity(signal_df, control_df, index_cols=index_cols + [sample_col],
                                         keep_cols=keep_cols,
                                         trim_5=trim_5, trim_3=trim_3,
                                         pos_col=pos_col, gene_col=gene_col)
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

def read_mutrate_files(directory_path='.', prefix ='', output_file=None, header=None, coverage_threshold=100, mutrate_threshold=0.2):
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
    mutrate_threshold : float, default 0.2
        Maximum mutrate allowed. Rows with mutrate above this threshold will be filtered out.

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
            rows_after_coverage = len(subset_df)

            # Filter by mutrate threshold
            subset_df = subset_df[subset_df['mutrate'] <= mutrate_threshold]
            rows_after = len(subset_df)

            # Reorder columns to put sample first
            subset_df = subset_df[['sample', 'gene', 'pos', 'mutrate', 'coverage', 'mutant']]

            all_data.append(subset_df)
            print(f"  - Added {rows_after} rows from {sample_name} (filtered {rows_before - rows_after_coverage} rows with coverage < {coverage_threshold}, {rows_after_coverage - rows_after} rows with mutrate > {mutrate_threshold})")
            
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

def lookup_window_region(regions_df, gene, window_start, window_end, symbol_col='symbol'):
    """
    Assign one window to the 5'UTR, CDS or 3'UTR it overlaps most, and return that region's extent.

    Coordinates are 1-based inclusive TRANSCRIPT coordinates on both sides, which is what
    gencode.vM38.selected.transcript.regions.csv stores and what a BAM aligned straight to the
    selected transcript FASTA already reports, so no genomic projection happens here.

    Parameters:
    -----------
    regions_df : pandas.DataFrame
        Transcript region table, one row per transcript, with symbol_col, 'biotype', 'length' and the
        'utr5_start'/'utr5_end', 'cds_start'/'cds_end', 'utr3_start'/'utr3_end' pairs. A region the
        transcript does not have is left blank.
    gene : str
        Transcript to look up, matched against symbol_col.
    window_start, window_end : int
        Window extent, 1-based inclusive.
    symbol_col : str, default 'symbol'
        Column of regions_df holding the transcript name.

    Returns:
    --------
    dict with:
        region : "5'UTR", 'CDS', "3'UTR", 'noncoding', 'unassigned' or 'no_annotation'
        region_start, region_end : extent of that region, or None when there is no extent to give
        overlap : nucleotides the window shares with the winning region
        n_regions_touched : how many of the three the window overlaps at all
        is_boundary_window : True when the window straddles more than one region
        biotype, length : copied from regions_df, or None when the transcript is absent
    """
    import pandas as pd

    gene_rows = regions_df[regions_df[symbol_col] == gene]
    if gene_rows.empty:
        return {'region': 'no_annotation', 'region_start': None, 'region_end': None, 'overlap': 0,
                'n_regions_touched': 0, 'is_boundary_window': False, 'biotype': None, 'length': None}
    row = gene_rows.iloc[0]

    ## Overlap in nucleotides between the window and each region, on the same 1-based inclusive
    ## convention; a region the transcript does not have (blank start/end) contributes 0.
    overlaps = []
    for region_name, start_col, end_col in [("5'UTR", 'utr5_start', 'utr5_end'),
                                            ('CDS', 'cds_start', 'cds_end'),
                                            ("3'UTR", 'utr3_start', 'utr3_end')]:
        region_start, region_end = row[start_col], row[end_col]
        if pd.isna(region_start) or pd.isna(region_end):
            overlaps.append((region_name, 0, None, None))
            continue
        region_start, region_end = int(region_start), int(region_end)
        overlap = min(window_end, region_end) - max(window_start, region_start) + 1
        overlaps.append((region_name, max(overlap, 0), region_start, region_end))

    ## MAX overlap wins, the same rule the bin annotation in the notebooks uses, so a title can never
    ## disagree with the table. Ties break toward the earlier region in 5' -> 3' order, which is what
    ## max() over a list in that order already does.
    region_name, overlap, region_start, region_end = max(overlaps, key=lambda item: item[1])
    n_regions_touched = sum(1 for _, region_overlap, _, _ in overlaps if region_overlap > 0)

    if overlap == 0:
        ## No overlap with any of the three: a noncoding transcript has no CDS to sit in, while a
        ## coding one means the window runs past the annotated end of the transcript.
        region_name = 'noncoding' if row['biotype'] == 'noncoding' else 'unassigned'
        region_start = region_end = None

    return {'region': region_name, 'region_start': region_start, 'region_end': region_end,
            'overlap': overlap, 'n_regions_touched': n_regions_touched,
            'is_boundary_window': n_regions_touched > 1,
            'biotype': row['biotype'], 'length': row['length']}


## Colours of the 5'UTR / CDS / 3'UTR track drawn along the top of a window panel. Deliberately
## clear of the gold window band and of the red/grey reactivity lines, so the track reads as
## structure and not as another series.
REGION_TRACK_COLORS = {"5'UTR": '#4C72B0', 'CDS': '#55A868', "3'UTR": '#8172B2'}


def transcript_region_extents(regions_df, gene, symbol_col='symbol'):
    """
    Return the 5'UTR / CDS / 3'UTR extents of one transcript, 5' -> 3', skipping the ones it lacks.

    Same table and same 1-based inclusive transcript coordinates as lookup_window_region, which
    assigns a single window to one region; this one hands back all of them so a panel can draw the
    whole structure rather than only the region its window sits in.

    Returns:
    --------
    list of (region_name, start, end), ordered by start. Empty when the transcript is absent from
    regions_df or carries no annotated region (a noncoding transcript with blank CDS/UTR columns).
    """
    import pandas as pd

    gene_rows = regions_df[regions_df[symbol_col] == gene]
    if gene_rows.empty:
        return []
    row = gene_rows.iloc[0]

    extents = []
    for region_name, start_col, end_col in [("5'UTR", 'utr5_start', 'utr5_end'),
                                            ('CDS', 'cds_start', 'cds_end'),
                                            ("3'UTR", 'utr3_start', 'utr3_end')]:
        region_start, region_end = row[start_col], row[end_col]
        if pd.isna(region_start) or pd.isna(region_end):
            continue
        extents.append((region_name, int(region_start), int(region_end)))
    return sorted(extents, key=lambda item: item[1])


def plot_window_reactivity(reactivity_df, low_corr_df, row_idx, flank=50,
                           group_col='cell_type', control_label='bulk', method_label='',
                           y_col='smoothed_reactivity', unit_col='biorep_celltype',
                           pos_col='pos', gene_col='gene',
                           coverage_col='coverage_signal',
                           zoom_range=None, views=('whole', 'zoom'),
                           show=True, figsize=(12, 4), axes=None, regions_df=None,
                           region_track=True, region_track_height=0.08):
    """
    Draw one flagged low-correlation window: the flagged group against the control.

    Each axis plots one line per unit_col (e.g. biorep_celltype), coloured by whether that unit
    belongs to the flagged group or to the control, with the flagged window shaded in gold. When the
    frame carries coverage_col, the same units' read depth is drawn as grey dashed lines on a
    log-scaled right hand axis, so a window can be read against the depth it was called on.

    Parameters:
    -----------
    reactivity_df : pandas.DataFrame
        Per-position reactivity, with gene_col, pos_col, group_col, unit_col and y_col. May be the
        full frame or a slice already cut down to this gene and the two groups - the gene/group mask
        below is applied either way and is a no-op on a slice. Callers drawing many windows should
        pass the slice: masking a ~1e7-row frame once per window dominates the runtime.
    low_corr_df : pandas.DataFrame
        Flagged windows, with gene_col, group_col, 'window_start' and 'window_end'.
    row_idx : label
        Index label of the row of low_corr_df to draw.
    flank : int, default 50
        Positions kept either side of the window in the zoom view, when zoom_range is not given.
    group_col : str, default 'cell_type'
        Column holding the flagged group and the control.
    control_label : str, default 'bulk'
        Value of group_col the flagged group is drawn against.
    method_label : str, default ''
        Free-text tag prefixed to the panel title, e.g. the correlation method or the panel name.
    y_col, unit_col, pos_col, gene_col : str
        Column names for the plotted value, the line unit, the position and the gene.
    coverage_col : str or None, default 'coverage_signal'
        Column holding the read depth, drawn as a grey dashed line per unit on a log-scaled right
        hand axis. None turns the overlay off, and so does a frame that does not carry the column -
        a caller reading the reactivity CSV with usecols has to ask for the coverage column to get
        the line. Positions at zero depth break the line. The right hand axis is scaled per panel,
        so panels sharing a reactivity y axis do NOT share the coverage one.
    zoom_range : (int, int) or None, default None
        Explicit (start, end) for the zoom view, overriding flank. Several windows drawn on one
        shared range stay directly comparable.
    views : tuple, default ('whole', 'zoom')
        Which of the two views to draw, always in whole -> zoom order.
    show : bool, default True
        Call plt.show() at the end. Ignored when axes is given.
    figsize : tuple, default (12, 4)
        Size of each figure created. Ignored when axes is given.
    regions_df : pandas.DataFrame or None, default None
        Transcript region table for lookup_window_region. When given, every panel title gains the
        5'UTR / CDS / 3'UTR the window sits in and that region's extent, and - unless region_track
        is turned off - every panel also gains the coloured region track described below.
    region_track : bool, default True
        Draw a coloured 5'UTR / CDS / 3'UTR bar along the top of each panel, with a dashed line at
        every region boundary that falls inside the view, so the gold window band can be read
        against the transcript structure and not only against its own coordinates. Needs regions_df;
        without it nothing is drawn. Regions are clipped to the view, so a region running past the
        edge gets a bar and no boundary line on that side.
    region_track_height : float, default 0.08
        Height of that track as a fraction of the axis, measured down from the top. The track is
        drawn over the data rather than beside it, so a taller track hides more of the curves.
    axes : list of matplotlib axes or None, default None
        One ready-made axis per entry of views, in that order, so the caller can lay several panels
        of the same window out on one figure. No figure is created in that case and the caller owns
        the layout.

    Returns:
    --------
    list of (ax, coverage_ax) tuples, one per drawn view, in whole -> zoom order. coverage_ax is the
    right hand axis carrying the depth lines, or None where the overlay was not drawn.
    """
    from matplotlib import pyplot as plt

    row = low_corr_df.loc[row_idx]
    gene = row[gene_col]
    win_start, win_end = row['window_start'], row['window_end']
    flagged_group = row[group_col]
    title_prefix = f'[{method_label}] ' if method_label else ''

    gene_df = reactivity_df[(reactivity_df[gene_col] == gene) &
                            (reactivity_df[group_col].isin([flagged_group, control_label]))]
    if gene_df.empty:
        raise ValueError(f"No reactivity data found for gene {gene!r}")

    palette = {flagged_group: 'tab:red', control_label: 'tab:grey'}

    ## Which transcript region the window sits in and how far that region runs, resolved once and
    ## appended to every panel title so the window can be read in context.
    region_label = ''
    if regions_df is not None:
        window_region = lookup_window_region(regions_df, gene, win_start, win_end)
        region_label = f" | {window_region['region']}"
        if window_region['region_start'] is not None:
            region_label += f" {window_region['region_start']}-{window_region['region_end']}"
        if window_region['is_boundary_window']:
            region_label += f" (straddles {window_region['n_regions_touched']} regions)"

    ## The whole 5' -> 3' structure, resolved once for every panel: the track below clips it to the
    ## view, so the same list serves the whole-transcript and the zoom panel.
    region_extents = (transcript_region_extents(regions_df, gene)
                      if regions_df is not None and region_track else [])

    ## The two views differ only in which positions they keep, what the title calls them and whether
    ## the x axis is pinned, so they are specced here and drawn by one loop below.
    if zoom_range is None:
        zoom_start, zoom_end = win_start - flank, win_end + flank
        zoom_label = f'zoom {win_start}-{win_end} ±{flank}'
    else:
        zoom_start, zoom_end = zoom_range
        zoom_label = f'window {win_start}-{win_end}, view {zoom_start}-{zoom_end}'
    zoom_df = gene_df[(gene_df[pos_col] >= zoom_start) & (gene_df[pos_col] <= zoom_end)]

    ## The x limits are set from the range and not from the data, so windows drawn on a shared range
    ## stay comparable even where a unit has no data out to the edge.
    ## The whole view has no pinned x limits, so its visible range is taken from the data - that is
    ## what matplotlib will autoscale to, give or take the margins, and it is what the region track
    ## has to be clipped and labelled against.
    view_specs = {
        'whole': (gene_df, 'whole transcript', None,
                  (gene_df[pos_col].min(), gene_df[pos_col].max())),
        'zoom': (zoom_df, zoom_label, (zoom_start, zoom_end), (zoom_start, zoom_end)),
    }

    ## Drawn only when the frame actually carries the column: callers cutting the reactivity CSV down
    ## with usecols routinely leave coverage out, and a missing overlay is the honest answer there.
    draw_coverage = coverage_col is not None and coverage_col in gene_df.columns

    def draw_panel(ax, sub_df, view_label, view_range):
        ## Depth first and on its own axis, so the reactivity curves the panel is about stay on top of
        ## it. zorder is set on the axis itself because a twin axis is drawn as a whole, above or
        ## below the host, and not per artist.
        coverage_ax = None
        if draw_coverage:
            coverage_ax = ax.twinx()
            for _, unit_df in sub_df.groupby([unit_col, group_col]):
                unit_df = unit_df.sort_values(pos_col)
                ## Depth spans orders of magnitude within a transcript, so the axis is log and not the
                ## values: the ticks stay in real reads and a shallow unit is still readable next to a
                ## deep one. A zero-depth position has no log, and is left as a gap in the line rather
                ## than dropped, which would join the two sides across it.
                coverage_ax.plot(unit_df[pos_col], unit_df[coverage_col].where(unit_df[coverage_col] > 0),
                                 color='dimgrey', linestyle='--', alpha=0.5, linewidth=1)
            coverage_ax.set_yscale('log')
            coverage_ax.set_ylabel(f'{coverage_col} (log)', color='dimgrey')
            coverage_ax.tick_params(axis='y', labelcolor='dimgrey', labelsize=8)
            ax.set_zorder(coverage_ax.get_zorder() + 1)
            ax.patch.set_visible(False)

        for (unit, group), unit_df in sub_df.groupby([unit_col, group_col]):
            unit_df = unit_df.sort_values(pos_col)
            ax.plot(unit_df[pos_col], unit_df[y_col],
                    color=palette.get(group, 'tab:blue'), alpha=0.7, linewidth=1)
        ax.axvspan(win_start, win_end, color='gold', alpha=0.25, zorder=0)

        ## Region track along the top: one bar per region that reaches into the view, clipped to it,
        ## plus a dashed line at each boundary inside the view. Drawn in axes fraction so it keeps
        ## its height whatever the reactivity y range is, and on ax rather than on the coverage twin
        ## so it stays above the depth lines.
        view_start, view_end = view_range
        track_bottom = 1.0 - region_track_height
        drawn_regions = []
        for region_name, region_start, region_end in region_extents:
            bar_start, bar_end = max(region_start, view_start), min(region_end, view_end)
            if bar_start >= bar_end:
                continue
            region_color = REGION_TRACK_COLORS.get(region_name, 'tab:blue')
            ax.axvspan(bar_start, bar_end, ymin=track_bottom, ymax=1.0,
                       color=region_color, alpha=0.35, zorder=2)
            drawn_regions.append((region_name, region_color))
            ## A sliver of a region has no room for its own name; the legend still carries it.
            if (bar_end - bar_start) > 0.08 * (view_end - view_start):
                ax.annotate(region_name, xy=((bar_start + bar_end) / 2, 1.0 - region_track_height / 2),
                            xycoords=('data', 'axes fraction'), ha='center', va='center',
                            fontsize=8, color=region_color, zorder=3)
            for boundary in (region_start, region_end):
                if view_start <= boundary <= view_end:
                    ax.axvline(boundary, color=region_color, linestyle='--',
                               alpha=0.6, linewidth=1, zorder=1)

        ax.set_xlabel('Position')
        ax.set_ylabel(y_col)
        ax.set_title(f'{title_prefix}{gene} — {flagged_group} vs {control_label} '
                     f'({view_label}){region_label}')
        legend_handles = [plt.Line2D([0], [0], color=color, lw=2, label=group)
                          for group, color in palette.items()
                          if group in sub_df[group_col].unique()]
        if draw_coverage:
            legend_handles.append(plt.Line2D([0], [0], color='dimgrey', lw=1, linestyle='--',
                                             label=f'{coverage_col} (log)'))
        legend_handles.extend(plt.Rectangle((0, 0), 1, 1, color=region_color, alpha=0.35,
                                            label=region_name)
                              for region_name, region_color in drawn_regions)
        if legend_handles:
            ax.legend(handles=legend_handles, fontsize=8, loc='best')
        return coverage_ax

    supplied_axes = list(axes) if axes is not None else None
    axes_list = []
    for view in ('whole', 'zoom'):        # fixed order, so supplied axes line up with `views`
        if view not in views:
            continue
        view_df, view_label, xlim, view_range = view_specs[view]
        if supplied_axes is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig, ax = None, supplied_axes.pop(0)
        coverage_ax = draw_panel(ax, view_df, view_label, view_range)
        if xlim is not None:
            ax.set_xlim(*xlim)
        if fig is not None:
            fig.tight_layout()
        axes_list.append((ax, coverage_ax))

    if show and axes is None:
        plt.show()

    return axes_list
