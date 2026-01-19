## This file contains helper functions that are used when analyzing structure data

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