"""
Filtering functions for ShapeData.

This module provides functions for filtering ShapeData objects by
coverage/mutation rate thresholds and cell metadata.
"""

import numpy as np
import pandas as pd
from scipy import sparse
from typing import Union, Optional, List, Dict, Any, Callable, TYPE_CHECKING

if TYPE_CHECKING:
    from .core import ShapeData


def filter_thresholds(
    data: 'ShapeData',
    coverage_threshold: float,
    mutrate_threshold: float
) -> 'ShapeData':
    """
    Filter data by coverage and mutation rate thresholds.
    Removes cells and positions with no data points after filtering.

    Also adds per-position cell statistics to genepos metadata:
    - n_cells: number of cells with valid data at each position
    - prop_cells: proportion of cells with valid data (n_cells / total_cells)

    Parameters:
    -----------
    data : ShapeData
        The ShapeData object to filter
    coverage_threshold : float
        Minimum coverage threshold
    mutrate_threshold : float
        Maximum mutation rate threshold

    Returns:
    --------
    ShapeData
        Filtered data with updated genepos containing 'n_cells' and 'prop_cells' columns
    """
    from .core import ShapeData

    print(f"Original shape: {data.shape}")
    print(f"Filtering: coverage >= {coverage_threshold}, mutrate <= {mutrate_threshold}")

    # Convert to COO format for easier filtering
    cov_coo = data.coverage.tocoo()
    mut_coo = data.mutrate.tocoo()

    # Verify same sparsity pattern (coverage and mutrate should be paired)
    if len(cov_coo.data) != len(mut_coo.data):
        raise ValueError(
            f"Coverage and mutrate matrices have different number of stored values: "
            f"{len(cov_coo.data)} vs {len(mut_coo.data)}. They must have identical sparsity patterns."
        )
    if not (np.array_equal(cov_coo.row, mut_coo.row) and np.array_equal(cov_coo.col, mut_coo.col)):
        raise ValueError(
            "Coverage and mutrate matrices have different sparsity patterns. "
            "They must have identical (row, col) positions."
        )

    # Create mask for values that pass both thresholds
    # Note: NaN in mutrate will fail the <= comparison (NaN <= x is False)
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

    genepos_len = len(data.genepos) if isinstance(data.genepos, list) else len(data.genepos)
    print(f"windows: {genepos_len:,} -> {len(valid_positions):,}")
    print(f"Cells: {data.n_cells:,} -> {len(valid_cells):,}")

    # Create mapping from old to new indices
    pos_mapping = {old_idx: new_idx for new_idx, old_idx in enumerate(valid_positions)}
    cell_mapping = {old_idx: new_idx for new_idx, old_idx in enumerate(valid_cells)}

    # Remap indices
    new_rows = np.array([pos_mapping[r] for r in filtered_rows])
    new_cols = np.array([cell_mapping[c] for c in filtered_cols])

    # Build new sparse matrices
    new_shape = (len(valid_positions), len(valid_cells))
    new_coverage = sparse.csr_matrix(
        (filtered_cov_data, (new_rows, new_cols)),
        shape=new_shape, dtype=np.int32
    )
    new_mutrate = sparse.csr_matrix(
        (filtered_mut_data, (new_rows, new_cols)),
        shape=new_shape, dtype=np.float32
    )

    # Filter genepos
    if isinstance(data.genepos, pd.DataFrame):
        new_genepos = data.genepos.iloc[valid_positions].reset_index(drop=True)
    else:
        new_genepos = [data.genepos[i] for i in valid_positions]

    # Filter cells
    if data.cells_is_df:
        new_cells = data.cells.iloc[valid_cells].copy()
    else:
        new_cells = [data.cells[i] for i in valid_cells]

    # Calculate per-position statistics: number and proportion of cells with valid data
    # In sparse matrix, stored values = valid data (coverage > 0 and passed thresholds)
    n_total_cells = len(valid_cells)

    # Count non-zero entries per row (position) - these are cells with valid data
    # Using CSR format for efficient row operations
    new_cov_csr = new_coverage.tocsr()
    n_cells_per_position = np.diff(new_cov_csr.indptr)  # Number of stored values per row
    prop_cells_per_position = n_cells_per_position / n_total_cells

    # Add to genepos metadata
    if isinstance(new_genepos, pd.DataFrame):
        new_genepos['n_cells'] = n_cells_per_position
        new_genepos['prop_cells'] = prop_cells_per_position

        # Calculate quantiles
        quantiles = [0, 0.25, 0.5, 0.75, 1.0]
        n_cells_quantiles = np.quantile(n_cells_per_position, quantiles)
        prop_cells_quantiles = np.quantile(prop_cells_per_position, quantiles)

        print(f"\nPer-position cell coverage statistics:")
        print(f"  n_cells   - min: {n_cells_quantiles[0]:.0f}, Q1: {n_cells_quantiles[1]:.0f}, "
              f"median: {n_cells_quantiles[2]:.0f}, Q3: {n_cells_quantiles[3]:.0f}, max: {n_cells_quantiles[4]:.0f}")
        print(f"  prop_cells - min: {prop_cells_quantiles[0]:.2%}, Q1: {prop_cells_quantiles[1]:.2%}, "
              f"median: {prop_cells_quantiles[2]:.2%}, Q3: {prop_cells_quantiles[3]:.2%}, max: {prop_cells_quantiles[4]:.2%}")
    else:
        print("Warning: genepos is not a DataFrame, skipping n_cells/prop_cells metadata")

    print(f"Filtered shape: {new_shape}")

    # Handle reactivity if it exists (dense array)
    new_reactivity = None
    if data.reactivity is not None:
        # For dense reactivity, we need to handle differently
        # Create a new dense array and fill with valid values
        new_reactivity = np.full(new_shape, np.nan, dtype=np.float32)
        for i, (old_pos, new_pos) in enumerate(pos_mapping.items()):
            for j, (old_cell, new_cell) in enumerate(cell_mapping.items()):
                new_reactivity[new_pos, new_cell] = data.reactivity[old_pos, old_cell]

    # Handle normalized_reactivity if it exists (dense array)
    new_normalized_reactivity = None
    if hasattr(data, 'normalized_reactivity') and data.normalized_reactivity is not None:
        new_normalized_reactivity = np.full(new_shape, np.nan, dtype=np.float32)
        for i, (old_pos, new_pos) in enumerate(pos_mapping.items()):
            for j, (old_cell, new_cell) in enumerate(cell_mapping.items()):
                new_normalized_reactivity[new_pos, new_cell] = data.normalized_reactivity[old_pos, old_cell]

    return ShapeData(new_coverage, new_mutrate, new_genepos, new_cells, new_reactivity, new_normalized_reactivity)


def filter_cells(
    data: 'ShapeData',
    cell_filters: Optional[Dict[str, Any]] = None,
    cell_indices: Optional[Union[List, np.ndarray, Callable]] = None
) -> 'ShapeData':
    """
    Filter data based on cell metadata or indices.

    Parameters:
    -----------
    data : ShapeData
        The ShapeData object to filter
    cell_filters : dict, optional
        Dictionary of column filters (only works when cells is a DataFrame).
        Keys are column names, values can be:
        - single value: exact match (e.g., {'leiden': 'K562'})
        - list: membership (e.g., {'leiden': ['K562', 'HEK293T']})
        - tuple of (min, max): range filter (e.g., {'n_genepos': (1000, 50000)})
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
    ShapeData
        Filtered data
    """
    from .core import ShapeData

    original_cell_list = data.cell_names
    n_original_cells = len(original_cell_list)

    print(f"Original shape: {data.shape}")
    print(f"Original cells: {n_original_cells}")

    valid_cell_indices = None

    # Step 1: Apply cell_indices filter
    if cell_indices is not None:
        print(f"\nApplying cell_indices filter...")

        if callable(cell_indices):
            cell_indices = cell_indices(original_cell_list)

        cell_indices = np.asarray(cell_indices)

        if cell_indices.dtype == bool:
            if len(cell_indices) != n_original_cells:
                raise ValueError(
                    f"Boolean mask length ({len(cell_indices)}) must match "
                    f"number of cells ({n_original_cells})"
                )
            valid_cell_indices = np.where(cell_indices)[0].tolist()
            print(f"  Boolean mask: {sum(cell_indices)} cells selected")
        elif cell_indices.dtype.kind in ['U', 'S', 'O']:
            barcode_set = set(cell_indices)
            valid_cell_indices = [i for i, c in enumerate(original_cell_list) if c in barcode_set]
            print(f"  Barcode filter: {len(valid_cell_indices)} cells matched")
        else:
            valid_cell_indices = cell_indices.tolist()
            print(f"  Index filter: {len(valid_cell_indices)} cells selected")

    # Step 2: Apply cell_filters (only works for DataFrame)
    if cell_filters and data.cells_is_df:
        # Truncate long lists in filter display
        filters_display = {}
        for k, v in cell_filters.items():
            if isinstance(v, list) and len(v) > 5:
                filters_display[k] = f"[{v[0]}, {v[1]}, ... ({len(v)} values)]"
            else:
                filters_display[k] = v
        print(f"\nApplying cell metadata filters: {filters_display}")
        cells_df = data.cells
        cell_mask = pd.Series(True, index=cells_df.index)

        for col, condition in cell_filters.items():
            if col not in cells_df.columns:
                print(f"  Warning: column '{col}' not found in cells_df, skipping")
                continue

            if isinstance(condition, tuple) and len(condition) == 2:
                min_val, max_val = condition
                if min_val is not None:
                    cell_mask &= (cells_df[col] >= min_val)
                if max_val is not None:
                    cell_mask &= (cells_df[col] <= max_val)
                print(f"  {col}: range [{min_val}, {max_val}]")
            elif isinstance(condition, list):
                cell_mask &= cells_df[col].isin(condition)
                if len(condition) > 5:
                    print(f"  {col}: in [{condition[0]}, {condition[1]}, ... ({len(condition)} values)]")
                else:
                    print(f"  {col}: in {condition}")
            else:
                cell_mask &= (cells_df[col] == condition)
                print(f"  {col}: == {condition}")

        cells_passing_filter = cells_df[cell_mask].index.tolist()
        metadata_indices = [original_cell_list.index(c) for c in cells_passing_filter]
        print(f"  Cells after metadata filter: {len(metadata_indices)}")

        if valid_cell_indices is not None:
            valid_cell_indices = sorted(set(valid_cell_indices) & set(metadata_indices))
            print(f"  Cells after combining filters: {len(valid_cell_indices)}")
        else:
            valid_cell_indices = metadata_indices

    elif cell_filters and not data.cells_is_df:
        print("Warning: cell_filters ignored because cells is not a DataFrame")

    # Step 3: Apply filtering or keep all
    if valid_cell_indices is not None:
        if len(valid_cell_indices) == 0:
            print("Warning: No cells pass the filters!")
            empty_genepos = pd.DataFrame() if isinstance(data.genepos, pd.DataFrame) else []
            if data.cells_is_df:
                return ShapeData(
                    sparse.csr_matrix((0, 0)),
                    sparse.csr_matrix((0, 0)),
                    empty_genepos,
                    data.cells.iloc[0:0],
                    None
                )
            else:
                return ShapeData(
                    sparse.csr_matrix((0, 0)),
                    sparse.csr_matrix((0, 0)),
                    empty_genepos,
                    [],
                    None
                )

        new_coverage = data.coverage[:, valid_cell_indices]
        new_mutrate = data.mutrate[:, valid_cell_indices] if data.mutrate is not None else None
        new_reactivity = data.reactivity[:, valid_cell_indices] if data.reactivity is not None else None
        new_normalized_reactivity = data.normalized_reactivity[:, valid_cell_indices] if hasattr(data, 'normalized_reactivity') and data.normalized_reactivity is not None else None

        if data.cells_is_df:
            new_cells = data.cells.iloc[valid_cell_indices].copy()
        else:
            new_cells = [data.cells[i] for i in valid_cell_indices]
    else:
        new_coverage = data.coverage.copy()
        new_mutrate = data.mutrate.copy() if data.mutrate is not None else None
        new_reactivity = data.reactivity.copy() if data.reactivity is not None else None
        new_normalized_reactivity = data.normalized_reactivity.copy() if hasattr(data, 'normalized_reactivity') and data.normalized_reactivity is not None else None
        new_cells = data.cells.copy() if data.cells_is_df else list(data.cells)

    # Keep all positions
    if isinstance(data.genepos, pd.DataFrame):
        new_genepos = data.genepos.copy()
    else:
        new_genepos = list(data.genepos)

    # Calculate stats
    cov_csr = new_coverage.tocsr()
    row_nnz = np.diff(cov_csr.indptr)
    n_all_na_windows = np.sum(row_nnz == 0)

    cov_csc = new_coverage.tocsc()
    col_nnz = np.diff(cov_csc.indptr)
    n_all_na_cells = np.sum(col_nnz == 0)

    n_filtered_cells = len(new_cells)
    print(f"\nAll-NA cells: {n_all_na_cells} / {n_filtered_cells}")
    print(f"All-NA windows: {n_all_na_windows} / {new_coverage.shape[0]}")
    print(f"Filtered shape: {new_coverage.shape}")

    return ShapeData(new_coverage, new_mutrate, new_genepos, new_cells, new_reactivity, new_normalized_reactivity)


def subset_by_cells(data: 'ShapeData', cell_list: List[str]) -> 'ShapeData':
    """
    Subset data to include only cells in the given list.
    Convenience function equivalent to filter_cells(data, cell_indices=cell_list).

    Parameters:
    -----------
    data : ShapeData
        The ShapeData object to subset
    cell_list : list of str
        Cell barcodes to keep

    Returns:
    --------
    ShapeData
        Subsetted data
    """
    return filter_cells(data, cell_indices=cell_list)


def filter_genepos(
    data: 'ShapeData',
    genepos_filters: Optional[Dict[str, Any]] = None,
    genepos_indices: Optional[Union[List, np.ndarray, Callable]] = None
) -> 'ShapeData':
    """
    Filter data based on genepos metadata or indices.

    Parameters:
    -----------
    data : ShapeData
        The ShapeData object to filter
    genepos_filters : dict, optional
        Dictionary of column filters (only works when genepos is a DataFrame).
        Keys are column names, values can be:
        - single value: exact match (e.g., {'gene': '18S'})
        - list: membership (e.g., {'gene': ['18S', '28S']})
        - tuple of (min, max): range filter (e.g., {'pos': (100, 500)})
          Use None for open-ended ranges: (100, None) for >= 100
    genepos_indices : list, array, or callable, optional
        Direct index-based filtering (works with both DataFrame and list genepos).
        Can be:
        - list/array of integers: indices of positions to keep
        - list/array of booleans: boolean mask (True = keep)
        - callable: function that takes genepos and returns indices or mask

    Returns:
    --------
    ShapeData
        Filtered data
    """
    from .core import ShapeData

    n_original_positions = data.n_positions

    print(f"Original shape: {data.shape}")
    print(f"Original positions: {n_original_positions}")

    valid_pos_indices = None

    # Step 1: Apply genepos_indices filter
    if genepos_indices is not None:
        print(f"\nApplying genepos_indices filter...")

        if callable(genepos_indices):
            genepos_indices = genepos_indices(data.genepos)

        genepos_indices = np.asarray(genepos_indices)

        if genepos_indices.dtype == bool:
            if len(genepos_indices) != n_original_positions:
                raise ValueError(
                    f"Boolean mask length ({len(genepos_indices)}) must match "
                    f"number of positions ({n_original_positions})"
                )
            valid_pos_indices = np.where(genepos_indices)[0].tolist()
            print(f"  Boolean mask: {sum(genepos_indices)} positions selected")
        else:
            valid_pos_indices = genepos_indices.tolist()
            print(f"  Index filter: {len(valid_pos_indices)} positions selected")

    # Step 2: Apply genepos_filters (only works for DataFrame)
    genepos_is_df = isinstance(data.genepos, pd.DataFrame)

    if genepos_filters and genepos_is_df:
        # Truncate long lists in filter display
        filters_display = {}
        for k, v in genepos_filters.items():
            if isinstance(v, list) and len(v) > 5:
                filters_display[k] = f"[{v[0]}, {v[1]}, ... ({len(v)} values)]"
            else:
                filters_display[k] = v
        print(f"\nApplying genepos metadata filters: {filters_display}")
        genepos_df = data.genepos
        pos_mask = pd.Series(True, index=range(len(genepos_df)))

        for col, condition in genepos_filters.items():
            if col not in genepos_df.columns:
                print(f"  Warning: column '{col}' not found in genepos, skipping")
                continue

            if isinstance(condition, tuple) and len(condition) == 2:
                min_val, max_val = condition
                if min_val is not None:
                    pos_mask &= (genepos_df[col].values >= min_val)
                if max_val is not None:
                    pos_mask &= (genepos_df[col].values <= max_val)
                print(f"  {col}: range [{min_val}, {max_val}]")
            elif isinstance(condition, list):
                pos_mask &= genepos_df[col].isin(condition).values
                if len(condition) > 5:
                    print(f"  {col}: in [{condition[0]}, {condition[1]}, ... ({len(condition)} values)]")
                else:
                    print(f"  {col}: in {condition}")
            else:
                pos_mask &= (genepos_df[col].values == condition)
                print(f"  {col}: == {condition}")

        metadata_indices = np.where(pos_mask)[0].tolist()
        print(f"  Positions after metadata filter: {len(metadata_indices)}")

        if valid_pos_indices is not None:
            valid_pos_indices = sorted(set(valid_pos_indices) & set(metadata_indices))
            print(f"  Positions after combining filters: {len(valid_pos_indices)}")
        else:
            valid_pos_indices = metadata_indices

    elif genepos_filters and not genepos_is_df:
        print("Warning: genepos_filters ignored because genepos is not a DataFrame")

    # Step 3: Apply filtering or keep all
    if valid_pos_indices is not None:
        if len(valid_pos_indices) == 0:
            print("Warning: No positions pass the filters!")
            empty_cells = data.cells.iloc[0:0] if data.cells_is_df else []
            return ShapeData(
                sparse.csr_matrix((0, 0)),
                sparse.csr_matrix((0, 0)),
                pd.DataFrame() if genepos_is_df else [],
                empty_cells,
                None
            )

        new_coverage = data.coverage[valid_pos_indices, :]
        new_mutrate = data.mutrate[valid_pos_indices, :] if data.mutrate is not None else None
        new_reactivity = data.reactivity[valid_pos_indices, :] if data.reactivity is not None else None
        new_normalized_reactivity = data.normalized_reactivity[valid_pos_indices, :] if hasattr(data, 'normalized_reactivity') and data.normalized_reactivity is not None else None

        if genepos_is_df:
            new_genepos = data.genepos.iloc[valid_pos_indices].reset_index(drop=True)
        else:
            new_genepos = [data.genepos[i] for i in valid_pos_indices]
    else:
        new_coverage = data.coverage.copy()
        new_mutrate = data.mutrate.copy() if data.mutrate is not None else None
        new_reactivity = data.reactivity.copy() if data.reactivity is not None else None
        new_normalized_reactivity = data.normalized_reactivity.copy() if hasattr(data, 'normalized_reactivity') and data.normalized_reactivity is not None else None
        new_genepos = data.genepos.copy() if genepos_is_df else list(data.genepos)

    # Keep all cells
    new_cells = data.cells.copy() if data.cells_is_df else list(data.cells)

    # Calculate stats
    cov_csr = new_coverage.tocsr()
    row_nnz = np.diff(cov_csr.indptr)
    n_all_na_positions = np.sum(row_nnz == 0)

    cov_csc = new_coverage.tocsc()
    col_nnz = np.diff(cov_csc.indptr)
    n_all_na_cells = np.sum(col_nnz == 0)

    n_filtered_positions = new_coverage.shape[0]
    print(f"\nAll-NA positions: {n_all_na_positions} / {n_filtered_positions}")
    print(f"All-NA cells: {n_all_na_cells} / {new_coverage.shape[1]}")
    print(f"Filtered shape: {new_coverage.shape}")

    return ShapeData(new_coverage, new_mutrate, new_genepos, new_cells, new_reactivity, new_normalized_reactivity)


def subset_by_genepos(data: 'ShapeData', genepos_indices: List[int]) -> 'ShapeData':
    """
    Subset data to include only positions in the given index list.
    Convenience function equivalent to filter_genepos(data, genepos_indices=genepos_indices).

    Parameters:
    -----------
    data : ShapeData
        The ShapeData object to subset
    genepos_indices : list of int
        Position indices to keep

    Returns:
    --------
    ShapeData
        Subsetted data
    """
    return filter_genepos(data, genepos_indices=genepos_indices)
