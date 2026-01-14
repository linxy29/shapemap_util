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