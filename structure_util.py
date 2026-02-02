## This file contains helper functions that are used when analyzing structure data

def filter_data_with_cellinfo(coverage_sparse, mutrate_sparse, genepos, cells,
                               cell_filters=None, cell_indices=None):
    """
    Filter sparse matrices based on cell metadata, barcode list, or indices.

    Parameters:
    -----------
    coverage_sparse : scipy.sparse.csr_matrix
        Coverage matrix (positions x cells)
    mutrate_sparse : scipy.sparse.csr_matrix
        Mutation rate matrix (positions x cells)
    genepos : pd.DataFrame or list
        Position identifiers
    cells : pd.DataFrame, list, or array-like
        Cell information. Can be:
        - pd.DataFrame: Cell metadata with cell names as index (supports cell_filters)
        - list/array: Simple list of cell barcodes (use cell_indices to filter)
    cell_filters : dict, optional
        Dictionary of column filters (only works when cells is a DataFrame).
        Keys are column names, values can be:
        - single value: exact match (e.g., {'leiden': 'K562'})
        - list: membership (e.g., {'leiden': ['K562', 'HEK293T']})
        - tuple of (min, max): range filter (e.g., {'n_windows': (1000, 50000)})
          Use None for open-ended ranges: (1000, None) for >= 1000
    cell_indices : list, array, or callable, optional
        Direct index-based filtering (works with both DataFrame and list cells).
        Can be:
        - list/array of integers: indices of cells to keep
        - list/array of booleans: boolean mask (True = keep)
        - list of strings: cell barcodes to keep
        - callable: function that takes cells and returns indices or mask

    Returns:
    --------
    tuple: (filtered_coverage, filtered_mutrate, filtered_genepos, filtered_cells)
        - filtered_cells: DataFrame if input was DataFrame, list if input was list

    Examples:
    ---------
    # With DataFrame and metadata filters:
    filter_data_with_cellinfo(cov, mut, genepos, cells_df,
                               cell_filters={'leiden': 'K562'})

    # With list of barcodes and index filter:
    filter_data_with_cellinfo(cov, mut, genepos, barcode_list,
                               cell_indices=[0, 1, 5, 10])  # keep these indices

    # With list of barcodes and barcode filter:
    filter_data_with_cellinfo(cov, mut, genepos, barcode_list,
                               cell_indices=['AAACCTGCAGTC', 'AAACCTGCAGTT'])

    # With boolean mask:
    mask = [True, False, True, ...]  # same length as cells
    filter_data_with_cellinfo(cov, mut, genepos, barcode_list, cell_indices=mask)

    # With callable:
    filter_data_with_cellinfo(cov, mut, genepos, barcode_list,
                               cell_indices=lambda x: [i for i, c in enumerate(x) if 'K562' in c])
    """
    import numpy as np
    import pandas as pd
    from scipy import sparse

    # Determine if cells is a DataFrame or a list
    cells_is_df = isinstance(cells, pd.DataFrame)

    if cells_is_df:
        cells_df = cells.copy()
        # Ensure cells_df has cell names as index
        if cells_df.index.name != 'cell' and 'cell' in cells_df.columns:
            cells_df = cells_df.set_index('cell')
        original_cell_list = cells_df.index.tolist()
        n_original_cells = len(cells_df)
    else:
        # Convert list to simple list if needed
        original_cell_list = list(cells)
        n_original_cells = len(original_cell_list)
        cells_df = None  # Will create later if needed

    print(f"Original shape: {coverage_sparse.shape}")
    print(f"Original cells: {n_original_cells}")

    # Determine valid cell indices based on filters
    valid_cell_indices = None

    # Step 1: Apply cell_indices filter (works for both DataFrame and list)
    if cell_indices is not None:
        print(f"\nApplying cell_indices filter...")

        if callable(cell_indices):
            # callable: call function with cells to get indices
            cell_indices = cell_indices(original_cell_list)

        cell_indices = np.asarray(cell_indices)

        if cell_indices.dtype == bool:
            # Boolean mask
            if len(cell_indices) != n_original_cells:
                raise ValueError(f"Boolean mask length ({len(cell_indices)}) must match "
                               f"number of cells ({n_original_cells})")
            valid_cell_indices = np.where(cell_indices)[0].tolist()
            print(f"  Boolean mask: {sum(cell_indices)} cells selected")
        elif cell_indices.dtype.kind in ['U', 'S', 'O']:
            # String array (barcodes)
            barcode_set = set(cell_indices)
            valid_cell_indices = [i for i, c in enumerate(original_cell_list)
                                  if c in barcode_set]
            print(f"  Barcode filter: {len(valid_cell_indices)} cells matched")
        else:
            # Integer indices
            valid_cell_indices = cell_indices.tolist()
            print(f"  Index filter: {len(valid_cell_indices)} cells selected")

    # Step 2: Apply cell_filters (only works for DataFrame)
    if cell_filters and cells_is_df:
        print(f"\nApplying cell metadata filters: {cell_filters}")
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
        metadata_indices = [original_cell_list.index(c) for c in cells_passing_filter]
        print(f"  Cells after metadata filter: {len(metadata_indices)}")

        # Combine with cell_indices if both are specified
        if valid_cell_indices is not None:
            # Intersection of both filters
            valid_cell_indices = sorted(set(valid_cell_indices) & set(metadata_indices))
            print(f"  Cells after combining filters: {len(valid_cell_indices)}")
        else:
            valid_cell_indices = metadata_indices

    elif cell_filters and not cells_is_df:
        print("Warning: cell_filters ignored because cells is not a DataFrame")

    # Step 3: Apply filtering or keep all
    if valid_cell_indices is not None:
        # Handle empty result
        if len(valid_cell_indices) == 0:
            print("Warning: No cells pass the filters!")
            empty_genepos = pd.DataFrame() if isinstance(genepos, pd.DataFrame) else []
            if cells_is_df:
                return (sparse.csr_matrix((0, 0)), sparse.csr_matrix((0, 0)),
                        empty_genepos, cells_df.iloc[0:0])
            else:
                return (sparse.csr_matrix((0, 0)), sparse.csr_matrix((0, 0)),
                        empty_genepos, [])

        # Subset matrices to only include filtered cells
        filtered_coverage = coverage_sparse[:, valid_cell_indices]
        filtered_mutrate = mutrate_sparse[:, valid_cell_indices]

        if cells_is_df:
            filtered_cells = cells_df.iloc[valid_cell_indices].copy()
        else:
            filtered_cells = [original_cell_list[i] for i in valid_cell_indices]
    else:
        # No filters, use all data
        filtered_coverage = coverage_sparse.copy()
        filtered_mutrate = mutrate_sparse.copy()
        if cells_is_df:
            filtered_cells = cells_df.copy()
        else:
            filtered_cells = original_cell_list.copy()

    # Filter genepos to match (keep all positions for now)
    if isinstance(genepos, pd.DataFrame):
        filtered_genepos = genepos.copy()
    elif isinstance(genepos, (pd.Series, pd.Index)):
        filtered_genepos = genepos.copy()
    else:
        filtered_genepos = list(genepos)

    # Step 4: Calculate and print all-NA cells and windows
    cov_csr = filtered_coverage.tocsr()
    row_nnz = np.diff(cov_csr.indptr)  # number of non-zeros per row
    n_all_na_windows = np.sum(row_nnz == 0)

    cov_csc = filtered_coverage.tocsc()
    col_nnz = np.diff(cov_csc.indptr)  # number of non-zeros per column
    n_all_na_cells = np.sum(col_nnz == 0)

    n_filtered_cells = len(filtered_cells) if isinstance(filtered_cells, list) else len(filtered_cells)
    print(f"\nAll-NA cells: {n_all_na_cells} / {n_filtered_cells}")
    print(f"All-NA windows: {n_all_na_windows} / {filtered_coverage.shape[0]}")
    print(f"Filtered shape: {filtered_coverage.shape}")

    return filtered_coverage, filtered_mutrate, filtered_genepos, filtered_cells


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