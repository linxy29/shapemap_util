## This file contains helper functions that are used when analyzing structure data

def filter_data_with_cellinfo(coverage_sparse, mutrate_sparse, genepos, cells_df,
                               cell_filters=None):
    """
    Filter sparse matrices based on cell metadata and identify all-NA cells/windows.

    Parameters:
    -----------
    coverage_sparse : scipy.sparse.csr_matrix
        Coverage matrix (positions x cells)
    mutrate_sparse : scipy.sparse.csr_matrix
        Mutation rate matrix (positions x cells)
    genepos : pd.DataFrame or list
        Position identifiers
    cells_df : pd.DataFrame
        Cell information dataframe with cell names as index.
        Example columns: n_windows, leiden, Total
    cell_filters : dict, optional
        Dictionary of column filters. Keys are column names, values can be:
        - single value: exact match (e.g., {'leiden': 'K562'})
        - list: membership (e.g., {'leiden': ['K562', 'HEK293T']})
        - tuple of (min, max): range filter (e.g., {'n_windows': (1000, 50000)})
          Use None for open-ended ranges: (1000, None) for >= 1000

    Returns:
    --------
    tuple: (filtered_coverage, filtered_mutrate, filtered_genepos, filtered_cells_df)
    """
    import numpy as np
    import pandas as pd
    from scipy import sparse

    # Ensure cells_df has cell names as index
    if cells_df.index.name != 'cell' and 'cell' in cells_df.columns:
        cells_df = cells_df.set_index('cell')

    print(f"Original shape: {coverage_sparse.shape}")
    print(f"Original cells: {len(cells_df)}")

    # Store original cell list for index mapping
    original_cell_list = cells_df.index.tolist()

    # Step 1: Apply cell metadata filters
    if cell_filters:
        print(f"\nApplying cell filters: {cell_filters}")
        cell_mask = pd.Series(True, index=cells_df.index)

        for col, condition in cell_filters.items():
            if col not in cells_df.columns:
                print(f"  Warning: column '{col}' not found in cells_df, skipping")
                continue

            if isinstance(condition, tuple) and len(condition) == 2:
                # Range filter (min, max)
                min_val, max_val = condition
                if min_val is not None:
                    cell_mask &= (cells_df[col] >= min_val)
                if max_val is not None:
                    cell_mask &= (cells_df[col] <= max_val)
                print(f"  {col}: range [{min_val}, {max_val}]")
            elif isinstance(condition, list):
                # Membership filter
                cell_mask &= cells_df[col].isin(condition)
                print(f"  {col}: in {condition}")
            else:
                # Exact match
                cell_mask &= (cells_df[col] == condition)
                print(f"  {col}: == {condition}")

        # Get indices of cells that pass metadata filters
        cells_passing_filter = cells_df[cell_mask].index.tolist()
        valid_cell_indices = [original_cell_list.index(c) for c in cells_passing_filter]

        print(f"Cells after metadata filter: {len(cells_passing_filter)}")

        # Handle empty result
        if len(valid_cell_indices) == 0:
            print("Warning: No cells pass the metadata filters!")
            empty_genepos = pd.DataFrame() if isinstance(genepos, pd.DataFrame) else []
            return (sparse.csr_matrix((0, 0)), sparse.csr_matrix((0, 0)),
                    empty_genepos, cells_df.iloc[0:0])

        # Subset matrices to only include filtered cells
        filtered_coverage = coverage_sparse[:, valid_cell_indices]
        filtered_mutrate = mutrate_sparse[:, valid_cell_indices]
        filtered_cells_df = cells_df.loc[cells_passing_filter].copy()
    else:
        # No cell filters, use all data
        filtered_coverage = coverage_sparse.copy()
        filtered_mutrate = mutrate_sparse.copy()
        filtered_cells_df = cells_df.copy()

    # Filter genepos to match (keep all positions for now)
    if isinstance(genepos, pd.DataFrame):
        filtered_genepos = genepos.copy()
    elif isinstance(genepos, (pd.Series, pd.Index)):
        filtered_genepos = genepos.copy()
    else:
        filtered_genepos = list(genepos)

    # Step 2: Calculate and print all-NA cells and windows
    cov_csr = filtered_coverage.tocsr()
    row_nnz = np.diff(cov_csr.indptr)  # number of non-zeros per row
    n_all_na_windows = np.sum(row_nnz == 0)

    cov_csc = filtered_coverage.tocsc()
    col_nnz = np.diff(cov_csc.indptr)  # number of non-zeros per column
    n_all_na_cells = np.sum(col_nnz == 0)

    print(f"All-NA cells: {n_all_na_cells} / {len(filtered_cells_df)}")
    print(f"All-NA windows: {n_all_na_windows} / {filtered_coverage.shape[0]}")
    print(f"Filtered shape: {filtered_coverage.shape}")

    return filtered_coverage, filtered_mutrate, filtered_genepos, filtered_cells_df


def filter_sparse_data(coverage_sparse, mutrate_sparse, genepos, cells, 
                       coverage_threshold, mutrate_threshold):
    """
    Filter sparse matrices based on coverage and mutation rate thresholds.
    Removes cells and positions with no data points left after filtering.
    
    Parameters:
    -----------
    coverage_sparse : scipy.sparse.csr_matrix
        Coverage matrix (positions x cells)
    mutrate_sparse : scipy.sparse.csr_matrix
        Mutation rate matrix (positions x cells)
    genepos : list
        Position identifiers (chrom.pos)
    cells : list
        Cell/sample names
    coverage_threshold : float
        Minimum coverage threshold
    mutrate_threshold : float
        Maximum mutation rate threshold
    
    Returns:
    --------
    tuple: (filtered_coverage, filtered_mutrate, filtered_genepos, filtered_cells)
    """
    import numpy as np
    import pandas as pd
    from scipy import sparse
    
    print(f"Original shape: {coverage_sparse.shape}")
    print(f"Filtering: coverage >= {coverage_threshold}, mutrate <= {mutrate_threshold}")
    
    # Convert to COO format for easier filtering
    cov_coo = coverage_sparse.tocoo()
    mut_coo = mutrate_sparse.tocoo()
    
    # Create mask for values that pass both thresholds
    mask = (cov_coo.data >= coverage_threshold) & (mut_coo.data <= mutrate_threshold)
    
    print(f"Data points: {len(cov_coo.data):,} -> {mask.sum():,} (kept {100*mask.sum()/len(mask):.2f}%)")
    
    # Apply mask
    filtered_rows = cov_coo.row[mask]
    filtered_cols = cov_coo.col[mask]
    filtered_cov_data = cov_coo.data[mask]
    filtered_mut_data = mut_coo.data[mask]
    
    # Find positions and cells that have at least one data point
    valid_positions = np.unique(filtered_rows)
    valid_cells = np.unique(filtered_cols)
    
    print(f"windows: {len(genepos):,} -> {len(valid_positions):,}")
    print(f"Cells: {len(cells):,} -> {len(valid_cells):,}")
    
    # Create mapping from old to new indices
    pos_mapping = {old_idx: new_idx for new_idx, old_idx in enumerate(valid_positions)}
    cell_mapping = {old_idx: new_idx for new_idx, old_idx in enumerate(valid_cells)}
    
    # Remap indices
    new_rows = np.array([pos_mapping[r] for r in filtered_rows])
    new_cols = np.array([cell_mapping[c] for c in filtered_cols])
    
    # Build new sparse matrices
    new_shape = (len(valid_positions), len(valid_cells))
    filtered_coverage = sparse.csr_matrix(
        (filtered_cov_data, (new_rows, new_cols)),
        shape=new_shape, dtype=np.int32
    )
    filtered_mutrate = sparse.csr_matrix(
        (filtered_mut_data, (new_rows, new_cols)),
        shape=new_shape, dtype=np.float32
    )
    
    # Filter position and cell lists
    # Handle genepos based on its type (list, array, or DataFrame)
    if isinstance(genepos, pd.DataFrame):
        # Keep all columns and reset index
        filtered_genepos = genepos.iloc[valid_positions].reset_index(drop=True)
    elif isinstance(genepos, (pd.Series, pd.Index)):
        filtered_genepos = genepos.iloc[valid_positions].tolist()
    else:
        filtered_genepos = [genepos[i] for i in valid_positions]
    
    # Handle cells based on its type
    if isinstance(cells, pd.DataFrame):
        filtered_cells = cells.iloc[valid_cells].reset_index(drop=True)
    elif isinstance(cells, (pd.Series, pd.Index)):
        filtered_cells = cells.iloc[valid_cells].tolist()
    else:
        filtered_cells = [cells[i] for i in valid_cells]
    
    print(f"Filtered shape: {filtered_coverage.shape}")
    
    return filtered_coverage, filtered_mutrate, filtered_genepos, filtered_cells

def generate_cluster_dataframes(clustering_df, coverage_sparse, mutrate_sparse, 
                                genepos, cells, cluster_column='leiden'):
    """
    Generate per-cluster dataframes with aggregated mutation data.
    
    Parameters:
    -----------
    clustering_df : pd.DataFrame
        DataFrame with 'cell' and cluster columns
    coverage_sparse : scipy.sparse.csr_matrix
        Coverage matrix (positions x cells)
    mutrate_sparse : scipy.sparse.csr_matrix
        Mutation rate matrix (positions x cells)
    genepos : pd.DataFrame or list
        Position information. If DataFrame, all columns are maintained in output.
    cells : list
        Cell/sample names
    cluster_column : str
        Name of the clustering column (default: 'leiden')
    
    Returns:
    --------
    list of pd.DataFrame
        List of dataframes, one per cluster
    """
    import numpy as np
    import pandas as pd

    # Create cell to index mapping
    cell_to_idx = {cell: idx for idx, cell in enumerate(cells)}
    
    # Get clustering info for cells that exist in our data
    clustering_subset = clustering_df[clustering_df['cell'].isin(cells)].copy()
    print(f"Cells in clustering data: {len(clustering_df)}")
    print(f"Cells in both datasets: {len(clustering_subset)}")
    
    # Get unique clusters
    clusters = clustering_subset[cluster_column].unique()
    print(f"Number of clusters: {len(clusters)}")
    
    # Check if genepos is a DataFrame
    genepos_is_df = isinstance(genepos, pd.DataFrame)
    
    cluster_dfs = []
    
    for cluster in clusters:
        print(f"\nProcessing cluster: {cluster}")
        
        # Get cells in this cluster
        cluster_cells = clustering_subset[
            clustering_subset[cluster_column] == cluster
        ]['cell'].tolist()
        
        # Get indices of these cells
        cluster_cell_indices = [cell_to_idx[cell] for cell in cluster_cells if cell in cell_to_idx]
        
        if len(cluster_cell_indices) == 0:
            print(f"  No cells found in data for cluster {cluster}")
            continue
        
        print(f"  Cells in cluster: {len(cluster_cell_indices)}")
        
        # Extract submatrices for this cluster
        cluster_cov = coverage_sparse[:, cluster_cell_indices]
        cluster_mut = mutrate_sparse[:, cluster_cell_indices]
        
        # Convert to COO for processing
        cluster_cov_coo = cluster_cov.tocoo()
        cluster_mut_coo = cluster_mut.tocoo()
        
        # Group by position to get min coverage and average mutrate
        pos_data = {}
        for i in range(len(cluster_cov_coo.data)):
            pos_idx = cluster_cov_coo.row[i]
            cov = cluster_cov_coo.data[i]
            mut = cluster_mut_coo.data[i]
            
            if pos_idx not in pos_data:
                pos_data[pos_idx] = {'coverages': [], 'mutrates': []}
            
            pos_data[pos_idx]['coverages'].append(cov)
            pos_data[pos_idx]['mutrates'].append(mut)
        
        # Build dataframe
        if genepos_is_df:
            # If genepos is DataFrame, merge with position info
            records = []
            for pos_idx, data in pos_data.items():
                record = {
                    'sample': f'dmso_{cluster}',
                    'mutrate': np.mean(data['mutrates']),
                    'coverage': int(np.min(data['coverages'])),
                    'mutant': np.nan
                }
                records.append(record)
            
            # Create dataframe with aggregated stats
            stats_df = pd.DataFrame(records)
            
            # Get the position info for positions that have data
            pos_indices = list(pos_data.keys())
            genepos_subset = genepos.iloc[pos_indices].reset_index(drop=True)
            
            # Concatenate position info with stats
            cluster_df = pd.concat([stats_df, genepos_subset], axis=1)
            
        else:
            # If genepos is list, parse gene.pos format
            records = []
            for pos_idx, data in pos_data.items():
                gene_pos = str(genepos[pos_idx])
                
                # Split gene.pos format
                parts = gene_pos.split('.')
                gene = parts[0] if len(parts) > 0 else gene_pos
                pos = parts[1] if len(parts) > 1 else '0'
                
                records.append({
                    'sample': f'dmso_{cluster}',
                    'gene': gene,
                    'pos': int(pos),
                    'mutrate': np.mean(data['mutrates']),
                    'coverage': int(np.min(data['coverages'])),
                    'mutant': np.nan
                })
            
            cluster_df = pd.DataFrame(records)
        
        print(f"  Generated {len(cluster_df)} records")
        cluster_dfs.append(cluster_df)
    
    return cluster_dfs


# Main pipeline
def process_clustering_pipeline(clustering_df, coverage_sparse, mutrate_sparse, 
                                genepos, cells, coverage_threshold, mutrate_threshold,
                                cluster_column='leiden'):
    """
    Complete pipeline: filter data and generate cluster dataframes.
    
    Returns:
    --------
    pd.DataFrame
        Concatenated dataframe of all clusters
    """
    import pandas as pd

    # Step 1: Filter data
    print("=" * 60)
    print("STEP 1: Filtering data")
    print("=" * 60)
    filtered_cov, filtered_mut, filtered_genepos, filtered_cells = filter_sparse_data(
        coverage_sparse, mutrate_sparse, genepos, cells,
        coverage_threshold, mutrate_threshold
    )
    
    # Step 2: Generate cluster dataframes
    print("\n" + "=" * 60)
    print("STEP 2: Generating cluster dataframes")
    print("=" * 60)
    cluster_dfs = generate_cluster_dataframes(
        clustering_df, filtered_cov, filtered_mut,
        filtered_genepos, filtered_cells, cluster_column
    )
    
    # Step 3: Concatenate all dataframes
    print("\n" + "=" * 60)
    print("STEP 3: Concatenating results")
    print("=" * 60)
    if len(cluster_dfs) == 0:
        print("Warning: No cluster dataframes generated!")
        return pd.DataFrame()
    
    final_df = pd.concat(cluster_dfs, ignore_index=True)
    print(f"Final dataframe shape: {final_df.shape}")
    print(f"\nSample distribution:")
    print(final_df['sample'].value_counts())
    
    return final_df