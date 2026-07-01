"""
Analysis functions for ShapeData.

This module provides functions for calculating reactivity, cell correlations,
AUROC metrics, and statistics from ShapeData objects.
"""

import numpy as np
import pandas as pd
from typing import Optional, List, Union, TYPE_CHECKING
from sklearn import metrics

if TYPE_CHECKING:
    from .core import ShapeData


def calculate_reactivity(
    data: 'ShapeData',
    cluster_column: str = 'leiden',
    control_prefix: str = 'mutrate_control_',
    store_as: str = 'reactivity',
    normalize: bool = True,
    normalize_method: str = 'box',
    store_normalized_as: str = 'normalized_reactivity',
    drop_mutrate: bool = False,
    filter_all_na: bool = False,
    verbose: bool = True
) -> np.ndarray:
    """
    Calculate reactivity matrix by subtracting cluster-specific control mutrate.

    Reactivity = mutrate - mutrate_control_{cluster}

    For each cell, the control mutrate from its assigned cluster is subtracted
    from the mutation rate. Values are set to NaN where there is no coverage
    or no control data for that position/cluster combination.

    Optionally, also computes normalized reactivity per cell per gene and stores
    it in a separate matrix.

    Parameters:
    -----------
    data : ShapeData
        The ShapeData object to calculate reactivity for
    cluster_column : str
        Column name in cells DataFrame containing cluster assignments (default: 'leiden')
    control_prefix : str
        Prefix for control mutrate columns in genepos (default: 'mutrate_control_')
        The full column name is {control_prefix}{cluster}
    store_as : str
        Store the reactivity matrix as an attribute with this name
        (default: 'reactivity'). Access via data.reactivity after calling.
    normalize : bool
        If True, also compute normalized reactivity per cell per gene
        (default: True)
    normalize_method : str
        Normalization method (only used if normalize=True):
        - 'box': box normalization (divide by mean of 90th-95th percentile values,
                 then clip negatives to 0). Standard SHAPE-MaP normalization.
        - 'zscore': z-score normalization (x - mean) / std
        - 'minmax': min-max normalization to [0, 1]
        (default: 'box')
    store_normalized_as : str
        Store the normalized reactivity matrix as an attribute with this name
        (default: 'normalized_reactivity'). Access via data.normalized_reactivity.
    drop_mutrate : bool
        If True, set data.mutrate to None after calculation to free memory
        (default: False)
    filter_all_na : bool
        If True, filter out gene.pos rows that have all NA values in the
        reactivity matrix. This also filters corresponding rows in genepos,
        coverage, and mutrate matrices. (default: False)
    verbose : bool
        Whether to print progress information (default: True)

    Returns:
    --------
    None
        The reactivity matrix is stored as data.{store_as} (default: data.reactivity)
        The normalized reactivity matrix is stored as data.{store_normalized_as}
        (default: data.normalized_reactivity) if normalize=True

    Raises:
    -------
    ValueError
        If cells is not a DataFrame or cluster_column not found
        If genepos is not a DataFrame

    Example:
    --------
    >>> # First add control data
    >>> data.add_control_data(dmso_df)
    >>> # Then calculate reactivity (with box normalization by default)
    >>> calculate_reactivity(data, cluster_column='leiden')
    >>> # Access stored matrices
    >>> print(data.reactivity.shape)           # raw reactivity
    >>> print(data.normalized_reactivity.shape)  # normalized per cell per gene
    """
    # Validate inputs
    if not data.cells_is_df:
        raise ValueError(
            "cells must be a DataFrame with cluster assignments. "
            "Use join_cell_metadata() to add cluster information first."
        )

    if cluster_column not in data.cells.columns:
        raise ValueError(
            f"Column '{cluster_column}' not found in cells DataFrame. "
            f"Available columns: {list(data.cells.columns)}"
        )

    if not isinstance(data.genepos, pd.DataFrame):
        raise ValueError(
            "genepos must be a DataFrame with control mutrate columns. "
            "Use add_control_data() to add control data first."
        )

    # Initialize reactivity matrix as dense array filled with NaN
    n_positions, n_cells = data.mutrate.shape
    reactivity_matrix = np.full((n_positions, n_cells), np.nan, dtype=np.float32)

    if verbose:
        print(f"Calculating reactivity matrix...")
        print(f"Output matrix shape: {reactivity_matrix.shape}")

    # Get unique clusters
    clusters = data.cells[cluster_column].dropna().unique()
    if verbose:
        print(f"\nClusters found: {list(clusters)}")

    for cluster in clusters:
        if verbose:
            print(f"\n{'='*60}")
            print(f"Processing cluster: {cluster}")
            print(f"{'='*60}")

        # 1) Get cell column indices for this cluster (positional index)
        cluster_mask = data.cells[cluster_column] == cluster
        cluster_cell_indices = np.where(cluster_mask)[0]
        if verbose:
            print(f"Number of cells in cluster: {len(cluster_cell_indices)}")

        # 2) Get mutrate_control for this cluster and find valid positions
        control_col = f'{control_prefix}{cluster}'
        if control_col not in data.genepos.columns:
            if verbose:
                print(f"Warning: {control_col} not found in genepos, skipping cluster")
            continue

        mutrate_control = data.genepos[control_col].values
        valid_pos_mask = ~np.isnan(mutrate_control)
        valid_pos_indices = np.where(valid_pos_mask)[0]
        if verbose:
            print(f"Positions with valid control data: {len(valid_pos_indices)} / {n_positions}")

        # 3) Extract sub-matrices for this cluster's cells and valid positions
        sub_mutrate = data.mutrate[valid_pos_indices, :][:, cluster_cell_indices].toarray()
        sub_coverage = data.coverage[valid_pos_indices, :][:, cluster_cell_indices].toarray()
        control_valid = mutrate_control[valid_pos_mask]

        # 4) Calculate reactivity: mutrate - mutrate_control
        sub_reactivity = sub_mutrate - control_valid[:, np.newaxis]

        # 5) Set to NaN where coverage is 0 (no data)
        # Note: mutrate == 0 with coverage > 0 is valid (no mutations)
        sub_reactivity[sub_coverage == 0] = np.nan

        # 6) Fill into the output matrix at correct positions
        for i, cell_idx in enumerate(cluster_cell_indices):
            reactivity_matrix[valid_pos_indices, cell_idx] = sub_reactivity[:, i]

        # Print statistics
        if verbose:
            valid_values = sub_reactivity[~np.isnan(sub_reactivity)]
            if len(valid_values) > 0:
                print(f"\nReactivity statistics for {cluster}:")
                print(f"  Non-NaN values: {len(valid_values):,}")
                print(f"  Min: {valid_values.min():.6f}")
                print(f"  Max: {valid_values.max():.6f}")
                print(f"  Mean: {valid_values.mean():.6f}")

    if verbose:
        print(f"\n{'='*60}")
        print("Final reactivity matrix:")
        print(f"  Shape: {reactivity_matrix.shape}")
        total_non_nan = np.sum(~np.isnan(reactivity_matrix))
        print(f"  Non-NaN values: {total_non_nan:,}")
        print(f"  Sparsity (NaN): {100 * (1 - total_non_nan / reactivity_matrix.size):.2f}%")

    # Store as attribute if requested
    if store_as:
        setattr(data, store_as, reactivity_matrix)
        if verbose:
            print(f"\nStored as data.{store_as}")

    # Normalize reactivity per cell per gene
    if normalize:
        if verbose:
            print(f"\n{'='*60}")
            print(f"Normalizing reactivity per cell per gene (method: {normalize_method})")

        # Get unique genes
        genes = data.genepos['gene'].unique()
        normalized_matrix = np.full_like(reactivity_matrix, np.nan)

        for gene in genes:
            # Get position indices for this gene
            gene_mask = data.genepos['gene'] == gene
            gene_pos_indices = np.where(gene_mask)[0]

            if len(gene_pos_indices) == 0:
                continue

            # Extract reactivity for this gene (positions × cells)
            gene_reactivity = reactivity_matrix[gene_pos_indices, :]

            # Normalize each cell's values for this gene
            for cell_idx in range(gene_reactivity.shape[1]):
                cell_values = gene_reactivity[:, cell_idx]
                valid_mask = ~np.isnan(cell_values)

                if valid_mask.sum() < 2:  # Need at least 2 values
                    continue

                valid_values = cell_values[valid_mask]

                if normalize_method == 'box':
                    # Box normalization (standard SHAPE-MaP):
                    # 1. Find 90th-95th percentile values
                    # 2. Take mean of those values as normalization factor
                    # 3. Divide all values by this factor
                    # 4. Clip negative values to 0
                    p90 = np.percentile(valid_values, 90)
                    p95 = np.percentile(valid_values, 95)
                    box_values = valid_values[(valid_values >= p90) & (valid_values <= p95)]

                    if len(box_values) > 0:
                        norm_factor = np.mean(box_values)
                        if norm_factor > 0:
                            normalized = cell_values / norm_factor
                            normalized = np.clip(normalized, 0, None)  # Clip negatives to 0
                            normalized_matrix[gene_pos_indices, cell_idx] = normalized
                        else:
                            # If norm factor is 0 or negative, use raw values clipped
                            normalized_matrix[gene_pos_indices, cell_idx] = np.clip(cell_values, 0, None)
                    else:
                        # If no values in 90-95 percentile range, use 95th percentile
                        if p95 > 0:
                            normalized = cell_values / p95
                            normalized = np.clip(normalized, 0, None)
                            normalized_matrix[gene_pos_indices, cell_idx] = normalized

                elif normalize_method == 'zscore':
                    # Z-score normalization: (x - mean) / std
                    mean_val = np.nanmean(valid_values)
                    std_val = np.nanstd(valid_values)
                    if std_val > 0:
                        normalized = (cell_values - mean_val) / std_val
                        normalized_matrix[gene_pos_indices, cell_idx] = normalized

                elif normalize_method == 'minmax':
                    # Min-max normalization to [0, 1]
                    min_val = np.nanmin(valid_values)
                    max_val = np.nanmax(valid_values)
                    if max_val > min_val:
                        normalized = (cell_values - min_val) / (max_val - min_val)
                        normalized_matrix[gene_pos_indices, cell_idx] = normalized

                else:
                    raise ValueError(f"Unknown normalize_method: {normalize_method}. "
                                   f"Use 'box', 'zscore', or 'minmax'.")

        # Store normalized matrix
        if store_normalized_as:
            setattr(data, store_normalized_as, normalized_matrix)

        if verbose:
            total_non_nan = np.sum(~np.isnan(normalized_matrix))
            print(f"Normalized reactivity matrix:")
            print(f"  Shape: {normalized_matrix.shape}")
            print(f"  Non-NaN values: {total_non_nan:,}")
            valid_norm = normalized_matrix[~np.isnan(normalized_matrix)]
            if len(valid_norm) > 0:
                print(f"  Min: {valid_norm.min():.4f}")
                print(f"  Max: {valid_norm.max():.4f}")
                print(f"  Mean: {valid_norm.mean():.4f}")
            print(f"\nStored as data.{store_normalized_as}")

    # Filter out gene.pos rows that have all NA values
    if filter_all_na:
        # Find rows that are NOT all NA (i.e., have at least one non-NA value)
        valid_pos_mask = ~np.all(np.isnan(reactivity_matrix), axis=1)
        n_valid = np.sum(valid_pos_mask)
        n_filtered = n_positions - n_valid

        if verbose:
            print(f"\n{'='*60}")
            print(f"Filtering all-NA positions:")
            print(f"  Positions before filtering: {n_positions:,}")
            print(f"  Positions with all NA: {n_filtered:,}")
            print(f"  Positions after filtering: {n_valid:,}")

        if n_filtered > 0:
            # Filter reactivity matrix
            reactivity_matrix = reactivity_matrix[valid_pos_mask, :]
            if store_as:
                setattr(data, store_as, reactivity_matrix)

            # Filter normalized_reactivity matrix if it exists
            if normalize and store_normalized_as:
                normalized_matrix = normalized_matrix[valid_pos_mask, :]
                setattr(data, store_normalized_as, normalized_matrix)

            # Filter genepos metadata
            if isinstance(data.genepos, pd.DataFrame):
                data.genepos = data.genepos.iloc[valid_pos_mask].reset_index(drop=True)
            else:
                data.genepos = [data.genepos[i] for i in range(len(data.genepos)) if valid_pos_mask[i]]

            # Filter coverage matrix
            data.coverage = data.coverage[valid_pos_mask, :]

            # Filter mutrate matrix if it exists
            if data.mutrate is not None:
                data.mutrate = data.mutrate[valid_pos_mask, :]

            if verbose:
                print(f"  Updated data.genepos: {len(data.genepos)} rows")
                print(f"  Updated data.coverage: {data.coverage.shape}")
                if data.mutrate is not None:
                    print(f"  Updated data.mutrate: {data.mutrate.shape}")
                if normalize and store_normalized_as:
                    print(f"  Updated data.{store_normalized_as}: {normalized_matrix.shape}")

    # Optionally drop mutrate to free memory
    if drop_mutrate:
        data.mutrate = None
        if verbose:
            print("Dropped mutrate matrix to free memory")


def calculate_reactivity_from_control(
    data: 'ShapeData',
    control_data: 'ShapeData',
    ref_column: str = 'ref_metabin',
    store_as: str = 'reactivity',
    normalize: bool = True,
    normalize_method: str = 'box',
    store_normalized_as: str = 'normalized_reactivity',
    verbose: bool = True
) -> None:
    """
    Calculate reactivity by subtracting matched control mutrate per cell.

    For each cell in data, finds the corresponding control cell via the
    ref_column mapping in control_data.cells, then computes:

        reactivity[pos, cell] = mutrate[pos, cell] - control_mutrate[pos, matched_control_cell]

    This is designed for metabin-level data where each treatment metabin
    has a matched control (e.g., DMSO) metabin.

    Parameters
    ----------
    data : ShapeData
        Treatment ShapeData (e.g., 2A3). Reactivity is stored here.
    control_data : ShapeData
        Control ShapeData (e.g., DMSO) with ref_column in cells DataFrame.
    ref_column : str
        Column in control_data.cells mapping control cells to treatment cells.
        Values should match cell names (index) in data.
    store_as : str
        Attribute name to store reactivity matrix (default: 'reactivity').
    normalize : bool
        If True, compute normalized reactivity per cell per gene (default: True).
    normalize_method : str
        Normalization method: 'box' (default), 'zscore', or 'minmax'.
    store_normalized_as : str
        Attribute name for normalized reactivity (default: 'normalized_reactivity').
    verbose : bool
        Whether to print progress information (default: True).

    Returns
    -------
    None
        Results stored as data.{store_as} and data.{store_normalized_as}.
    """
    # Validate inputs
    if data.mutrate is None:
        raise ValueError("data.mutrate is None. Treatment data must have a mutrate matrix.")
    if control_data.mutrate is None:
        raise ValueError("control_data.mutrate is None. Control data must have a mutrate matrix.")
    if not control_data.cells_is_df:
        raise ValueError("control_data.cells must be a DataFrame with ref_column.")
    if ref_column not in control_data.cells.columns:
        raise ValueError(
            f"Column '{ref_column}' not found in control_data.cells. "
            f"Available columns: {list(control_data.cells.columns)}"
        )

    # Build genepos ID lists for position alignment
    def _get_genepos_ids(genepos):
        if isinstance(genepos, pd.DataFrame):
            if genepos.index.name in ('gene.pos', 'position'):
                return list(genepos.index)
            elif 'gene' in genepos.columns and 'pos' in genepos.columns:
                return [f"{g}.{p}" for g, p in zip(genepos['gene'], genepos['pos'])]
            else:
                return list(genepos.index)
        else:
            return [str(x) for x in genepos]

    treat_genepos_ids = _get_genepos_ids(data.genepos)
    ctrl_genepos_ids = _get_genepos_ids(control_data.genepos)

    # Find shared positions and their indices in each dataset
    ctrl_id_to_idx = {gid: i for i, gid in enumerate(ctrl_genepos_ids)}
    shared_treat_indices = []
    shared_ctrl_indices = []
    for treat_i, gid in enumerate(treat_genepos_ids):
        if gid in ctrl_id_to_idx:
            shared_treat_indices.append(treat_i)
            shared_ctrl_indices.append(ctrl_id_to_idx[gid])

    shared_treat_indices = np.array(shared_treat_indices)
    shared_ctrl_indices = np.array(shared_ctrl_indices)

    if len(shared_treat_indices) == 0:
        raise ValueError("No shared positions found between treatment and control genepos.")

    if verbose:
        print(f"Position alignment:")
        print(f"  Treatment positions: {len(treat_genepos_ids)}")
        print(f"  Control positions: {len(ctrl_genepos_ids)}")
        print(f"  Shared positions: {len(shared_treat_indices)}")

    # Build mapping: ref_metabin_name -> control column index
    control_cell_names = (
        control_data.cells.index.tolist() if control_data.cells_is_df
        else control_data.cell_names
    )
    ref_to_ctrl_idx = {}
    for ctrl_idx, ctrl_name in enumerate(control_cell_names):
        ref_name = str(control_data.cells.loc[ctrl_name, ref_column])
        ref_to_ctrl_idx[ref_name] = ctrl_idx

    # Get treatment cell names
    treat_cell_names = data.cell_names

    # Initialize reactivity matrix (same shape as treatment data)
    n_positions, n_cells = data.mutrate.shape
    reactivity_matrix = np.full((n_positions, n_cells), np.nan, dtype=np.float32)

    if verbose:
        print(f"\nCalculating reactivity from matched control...")
        print(f"  Treatment cells: {n_cells}")
        print(f"  Control cells: {control_data.n_cells}")

    matched = 0
    unmatched = []
    for treat_idx, treat_name in enumerate(treat_cell_names):
        if treat_name not in ref_to_ctrl_idx:
            unmatched.append(treat_name)
            continue

        ctrl_idx = ref_to_ctrl_idx[treat_name]
        matched += 1

        # Get mutrate/coverage vectors at shared positions only
        treat_mut = data.mutrate[shared_treat_indices, treat_idx].toarray().ravel()
        ctrl_mut = control_data.mutrate[shared_ctrl_indices, ctrl_idx].toarray().ravel()
        treat_cov = data.coverage[shared_treat_indices, treat_idx].toarray().ravel()
        ctrl_cov = control_data.coverage[shared_ctrl_indices, ctrl_idx].toarray().ravel()

        # reactivity = treatment - control
        react = treat_mut - ctrl_mut

        # Set NaN where treatment or control has no coverage
        react[(treat_cov == 0) | (ctrl_cov == 0)] = np.nan

        # Fill into the output matrix at the shared treatment positions
        reactivity_matrix[shared_treat_indices, treat_idx] = react

    if verbose:
        print(f"  Matched cells: {matched}")
        if unmatched:
            print(f"  Unmatched cells ({len(unmatched)}): {unmatched[:5]}{'...' if len(unmatched) > 5 else ''}")
        total_non_nan = np.sum(~np.isnan(reactivity_matrix))
        print(f"  Non-NaN values: {total_non_nan:,}")
        valid = reactivity_matrix[~np.isnan(reactivity_matrix)]
        if len(valid) > 0:
            print(f"  Reactivity range: [{valid.min():.6f}, {valid.max():.6f}], mean={valid.mean():.6f}")

    # Store reactivity
    if store_as:
        setattr(data, store_as, reactivity_matrix)
        if verbose:
            print(f"\nStored as data.{store_as}")

    # Normalize reactivity per cell per gene
    if normalize:
        if not isinstance(data.genepos, pd.DataFrame) or 'gene' not in data.genepos.columns:
            raise ValueError(
                "Normalization requires data.genepos to be a DataFrame with a 'gene' column."
            )

        if verbose:
            print(f"\nNormalizing reactivity per cell per gene (method: {normalize_method})")

        genes = data.genepos['gene'].unique()
        normalized_matrix = np.full_like(reactivity_matrix, np.nan)

        for gene in genes:
            gene_mask = data.genepos['gene'] == gene
            gene_pos_indices = np.where(gene_mask)[0]

            if len(gene_pos_indices) == 0:
                continue

            gene_reactivity = reactivity_matrix[gene_pos_indices, :]

            for cell_idx in range(gene_reactivity.shape[1]):
                cell_values = gene_reactivity[:, cell_idx]
                valid_mask = ~np.isnan(cell_values)

                if valid_mask.sum() < 2:
                    continue

                valid_values = cell_values[valid_mask]

                if normalize_method == 'box':
                    p90 = np.percentile(valid_values, 90)
                    p95 = np.percentile(valid_values, 95)
                    box_values = valid_values[(valid_values >= p90) & (valid_values <= p95)]

                    if len(box_values) > 0:
                        norm_factor = np.mean(box_values)
                        if norm_factor > 0:
                            normalized = cell_values / norm_factor
                            normalized = np.clip(normalized, 0, None)
                            normalized_matrix[gene_pos_indices, cell_idx] = normalized
                        else:
                            normalized_matrix[gene_pos_indices, cell_idx] = np.clip(cell_values, 0, None)
                    else:
                        if p95 > 0:
                            normalized = cell_values / p95
                            normalized = np.clip(normalized, 0, None)
                            normalized_matrix[gene_pos_indices, cell_idx] = normalized

                elif normalize_method == 'zscore':
                    mean_val = np.nanmean(valid_values)
                    std_val = np.nanstd(valid_values)
                    if std_val > 0:
                        normalized = (cell_values - mean_val) / std_val
                        normalized_matrix[gene_pos_indices, cell_idx] = normalized

                elif normalize_method == 'minmax':
                    min_val = np.nanmin(valid_values)
                    max_val = np.nanmax(valid_values)
                    if max_val > min_val:
                        normalized = (cell_values - min_val) / (max_val - min_val)
                        normalized_matrix[gene_pos_indices, cell_idx] = normalized

                else:
                    raise ValueError(f"Unknown normalize_method: {normalize_method}. "
                                   f"Use 'box', 'zscore', or 'minmax'.")

        if store_normalized_as:
            setattr(data, store_normalized_as, normalized_matrix)

        if verbose:
            total_non_nan = np.sum(~np.isnan(normalized_matrix))
            print(f"Normalized reactivity matrix:")
            print(f"  Non-NaN values: {total_non_nan:,}")
            valid_norm = normalized_matrix[~np.isnan(normalized_matrix)]
            if len(valid_norm) > 0:
                print(f"  Range: [{valid_norm.min():.4f}, {valid_norm.max():.4f}], mean={valid_norm.mean():.4f}")
            print(f"\nStored as data.{store_normalized_as}")


def calculate_cell_correlation(
    data: 'ShapeData',
    gene: str,
    cluster_column: Optional[str] = None,
    method: str = 'pearson',
    min_shared_positions: int = 100,
    save_mean_as: Optional[str] = None,
    verbose: bool = True
) -> pd.DataFrame:
    """
    Calculate pairwise correlation between cells for a given gene.

    When ``cluster_column`` is provided, correlations are computed only between
    cells that share the same cluster. When ``cluster_column`` is ``None``,
    correlations are computed across all pairs of cells.

    Parameters:
    -----------
    data : ShapeData
        The ShapeData object to analyze
    gene : str
        Gene name to calculate correlations for
    cluster_column : str or None
        Column name in cells DataFrame containing cluster assignments. If None,
        correlations are computed across all pairs of cells (default: None).
    method : str
        Correlation method: 'pearson', 'spearman', or 'kendall' (default: 'pearson')
    min_shared_positions : int
        Minimum number of shared non-NaN positions required to calculate correlation
        (default: 100). Pairs with fewer shared positions get NaN correlation.
    save_mean_as : str or None
        Column name to save mean correlation for each cell in cells metadata.
        If None, uses f'mean_corr_{gene}' (default: None)
    verbose : bool
        Whether to print progress information (default: True)

    Returns:
    --------
    pd.DataFrame
        Correlation dataframe with columns:
        - cell1: first cell barcode
        - cell2: second cell barcode
        - cluster: cluster name (or 'all' when cluster_column is None)
        - correlation: pairwise correlation value
        - n_shared: number of shared non-NaN positions

    Raises:
    -------
    ValueError
        If reactivity matrix not found (call calculate_reactivity first)
        If gene not found in genepos
        If cluster_column is provided but cells is not a DataFrame or column missing

    Example:
    --------
    >>> # Calculate correlations for 18S gene within clusters
    >>> corr_df = calculate_cell_correlation(data, '18S', cluster_column='leiden')
    >>> print(corr_df.head())
    >>> # Mean correlation is saved to data.cells['mean_corr_18S']
    >>> print(data.cells[['leiden', 'mean_corr_18S']].head())
    >>>
    >>> # Calculate correlations across all cell pairs (no cluster grouping)
    >>> corr_df = calculate_cell_correlation(data, '18S')
    """
    from itertools import combinations

    # Validate inputs
    if not hasattr(data, 'reactivity') or data.reactivity is None:
        raise ValueError(
            "Reactivity matrix not found. Call calculate_reactivity() first."
        )

    if cluster_column is not None:
        if not data.cells_is_df:
            raise ValueError(
                "cells must be a DataFrame when cluster_column is specified."
            )

        if cluster_column not in data.cells.columns:
            raise ValueError(
                f"Column '{cluster_column}' not found in cells DataFrame. "
                f"Available columns: {list(data.cells.columns)}"
            )

    if not isinstance(data.genepos, pd.DataFrame):
        raise ValueError("genepos must be a DataFrame with 'gene' column.")

    if 'gene' not in data.genepos.columns:
        raise ValueError("genepos DataFrame must have a 'gene' column.")

    # Get positions for the specified gene
    gene_mask = data.genepos['gene'] == gene
    gene_positions = np.where(gene_mask)[0]

    if len(gene_positions) == 0:
        raise ValueError(
            f"Gene '{gene}' not found in genepos. "
            f"Available genes: {data.genepos['gene'].unique().tolist()}"
        )

    if verbose:
        print(f"Calculating cell correlations for gene: {gene}")
        print(f"Number of positions for {gene}: {len(gene_positions)}")

    # Extract reactivity submatrix for this gene (positions x cells)
    gene_reactivity = data.reactivity[gene_positions, :]

    # Get cell names and build cluster -> cell-indices mapping
    cell_names = data.cell_names

    if cluster_column is None:
        # Treat all cells as one group
        cluster_groups = [('all', np.arange(len(cell_names)))]
        if verbose:
            print("No cluster_column specified: computing correlations across all cell pairs")
    else:
        clusters = data.cells[cluster_column].dropna().unique()
        cluster_groups = [
            (cluster, np.where(data.cells[cluster_column] == cluster)[0])
            for cluster in clusters
        ]
        if verbose:
            print(f"Clusters: {list(clusters)}")

    # Store correlation results
    corr_results = []

    for cluster, cluster_cell_indices in cluster_groups:
        if verbose:
            print(f"\nProcessing cluster: {cluster}")

        cluster_cell_names = [cell_names[i] for i in cluster_cell_indices]
        n_cells = len(cluster_cell_indices)

        if verbose:
            print(f"  Number of cells: {n_cells}")

        if n_cells < 2:
            if verbose:
                print(f"  Skipping: need at least 2 cells for correlation")
            continue

        # Extract reactivity for this cluster's cells
        cluster_reactivity = gene_reactivity[:, cluster_cell_indices]  # (positions x cells)

        # Calculate pairwise correlations
        n_pairs = 0
        for i, j in combinations(range(n_cells), 2):
            cell1_vals = cluster_reactivity[:, i]
            cell2_vals = cluster_reactivity[:, j]

            # Find shared non-NaN positions
            valid_mask = ~np.isnan(cell1_vals) & ~np.isnan(cell2_vals)
            n_shared = np.sum(valid_mask)

            if n_shared >= min_shared_positions:
                x = cell1_vals[valid_mask]
                y = cell2_vals[valid_mask]

                if method == 'pearson':
                    # Pearson correlation
                    corr = np.corrcoef(x, y)[0, 1]
                elif method == 'spearman':
                    from scipy.stats import spearmanr
                    corr, _ = spearmanr(x, y)
                elif method == 'kendall':
                    from scipy.stats import kendalltau
                    corr, _ = kendalltau(x, y)
                else:
                    raise ValueError(f"Unknown method: {method}. Use 'pearson', 'spearman', or 'kendall'.")
            else:
                corr = np.nan

            corr_results.append({
                'cell1': cluster_cell_names[i],
                'cell2': cluster_cell_names[j],
                'cluster': cluster,
                'correlation': corr,
                'n_shared': n_shared
            })
            n_pairs += 1

        if verbose:
            valid_corrs = [r['correlation'] for r in corr_results[-n_pairs:] if not np.isnan(r['correlation'])]
            if valid_corrs:
                print(f"  Pairs calculated: {n_pairs}")
                print(f"  Valid correlations: {len(valid_corrs)}")
                print(f"  Mean correlation: {np.mean(valid_corrs):.4f}")

    # Create correlation dataframe
    corr_df = pd.DataFrame(corr_results)

    if verbose:
        print(f"\n{'='*60}")
        print(f"Total pairs: {len(corr_df)}")
        valid_corrs = corr_df['correlation'].dropna()
        print(f"Valid correlations: {len(valid_corrs)}")
        if len(valid_corrs) > 0:
            print(f"Overall mean correlation: {valid_corrs.mean():.4f}")
            print(f"Correlation range: [{valid_corrs.min():.4f}, {valid_corrs.max():.4f}]")

    # Calculate mean correlation for each cell and save to metadata
    if save_mean_as is None:
        save_mean_as = f'mean_corr_{gene}'

    # For each cell, calculate mean of all its pairwise correlations
    cell_mean_corr = {}
    for cell in cell_names:
        # Get correlations where this cell is either cell1 or cell2
        cell_corrs = corr_df[
            (corr_df['cell1'] == cell) | (corr_df['cell2'] == cell)
        ]['correlation'].dropna()

        if len(cell_corrs) > 0:
            cell_mean_corr[cell] = cell_corrs.mean()
        else:
            cell_mean_corr[cell] = np.nan

    # Add to cells metadata
    if data.cells_is_df:
        data.cells[save_mean_as] = data.cells.index.map(cell_mean_corr)

        if verbose:
            print(f"\nMean correlation saved to cells['{save_mean_as}']")
            valid_means = data.cells[save_mean_as].dropna()
            print(f"Cells with valid mean correlation: {len(valid_means)} / {len(data.cells)}")
    else:
        if verbose:
            print(f"\ncells is not a DataFrame; skipped saving mean correlation as '{save_mean_as}'")

    return corr_df


def get_cell_stats(
    data: 'ShapeData',
    verbose: bool = True
) -> None:
    """
    Compute statistics per cell and add to cells metadata.

    Adds the following columns to data.cells:
        - n_genepos: number of non-zero entries
        - mean_coverage: mean coverage across non-zero entries
        - mean_mutrate: mean mutation rate across non-zero entries

    Parameters:
    -----------
    data : ShapeData
        The ShapeData object to compute statistics for
    verbose : bool
        Whether to print progress information (default: True)

    Returns:
    --------
    None
        Statistics are added directly to data.cells
    """
    # Number of positions per cell
    n_genepos = np.array((data.coverage > 0).sum(axis=0)).flatten()

    # Mean coverage per cell (only for non-zero entries)
    cov_sum = np.array(data.coverage.sum(axis=0)).flatten()
    mean_coverage = np.divide(cov_sum, n_genepos, where=n_genepos > 0, out=np.zeros_like(cov_sum, dtype=float))

    # Mean mutation rate per cell
    mut_sum = np.array(data.mutrate.sum(axis=0)).flatten()
    mean_mutrate = np.divide(mut_sum, n_genepos, where=n_genepos > 0, out=np.zeros_like(mut_sum, dtype=float))

    # Add to cells metadata
    if data.cells_is_df:
        data.cells['n_genepos'] = n_genepos
        data.cells['mean_coverage'] = mean_coverage
        data.cells['mean_mutrate'] = mean_mutrate
        if verbose:
            print(f"Added columns to cells metadata: n_genepos, mean_coverage, mean_mutrate")
    else:
        # Convert to DataFrame first
        data.cells = pd.DataFrame({
            'n_genepos': n_genepos,
            'mean_coverage': mean_coverage,
            'mean_mutrate': mean_mutrate
        }, index=data.cell_names)
        data.cells.index.name = 'cell'
        if verbose:
            print(f"Converted cells to DataFrame and added: n_genepos, mean_coverage, mean_mutrate")


def get_position_stats(
    data: 'ShapeData',
    verbose: bool = True
) -> None:
    """
    Compute statistics per position and add to genepos metadata.

    Adds the following columns to data.genepos:
        - n_cells: number of cells with data
        - mean_coverage: mean coverage across cells
        - min_coverage: minimum coverage across cells with data
        - mean_mutrate: mean mutation rate across cells

    Parameters:
    -----------
    data : ShapeData
        The ShapeData object to compute statistics for
    verbose : bool
        Whether to print progress information (default: True)

    Returns:
    --------
    None
        Statistics are added directly to data.genepos
    """
    # Number of cells per position
    n_cells = np.array((data.coverage > 0).sum(axis=1)).flatten()

    # Mean coverage per position
    cov_sum = np.array(data.coverage.sum(axis=1)).flatten()
    mean_coverage = np.divide(cov_sum, n_cells, where=n_cells > 0, out=np.zeros_like(cov_sum, dtype=float))

    # Min coverage per position (minimum across cells that have data, i.e. coverage > 0)
    from scipy import sparse
    cov = data.coverage
    if sparse.issparse(cov):
        # Replace 0 with NaN to ignore cells without data, then take row min
        cov_dense = cov.toarray().astype(float)
    else:
        cov_dense = np.array(cov, dtype=float)
    cov_dense[cov_dense == 0] = np.nan
    with np.errstate(all='ignore'):
        min_coverage = np.nanmin(cov_dense, axis=1)
    min_coverage = np.where(np.isnan(min_coverage), 0, min_coverage)

    # Mean mutation rate per position
    mut_sum = np.array(data.mutrate.sum(axis=1)).flatten()
    mean_mutrate = np.divide(mut_sum, n_cells, where=n_cells > 0, out=np.zeros_like(mut_sum, dtype=float))

    # Add to genepos metadata
    if isinstance(data.genepos, pd.DataFrame):
        data.genepos['n_cells'] = n_cells
        data.genepos['mean_coverage'] = mean_coverage
        data.genepos['min_coverage'] = min_coverage
        data.genepos['mean_mutrate'] = mean_mutrate
        if verbose:
            print(f"Added columns to genepos metadata: n_cells, mean_coverage, min_coverage, mean_mutrate")
    else:
        # Convert to DataFrame first
        data.genepos = pd.DataFrame({
            'position': data.genepos,
            'n_cells': n_cells,
            'mean_coverage': mean_coverage,
            'min_coverage': min_coverage,
            'mean_mutrate': mean_mutrate
        })
        if verbose:
            print(f"Converted genepos to DataFrame and added: n_cells, mean_coverage, min_coverage, mean_mutrate")


def calculate_auroc(
    data: 'ShapeData',
    gene: str,
    dot_bracket: str,
    skip_positions: Optional[List[int]] = None,
    min_positions: Optional[int] = None,
    save_as: Optional[str] = None,
    verbose: bool = True
) -> None:
    """
    Calculate AUROC for each cell based on reactivity vs RNA secondary structure.

    For each cell, AUROC is computed using reactivity values as predictions and
    dot-bracket structure labels as ground truth (1.0 for unpaired '.', 0.0 for
    paired '(' or ')').

    Parameters:
    -----------
    data : ShapeData
        The ShapeData object with reactivity matrix
    gene : str
        Gene name to calculate AUROC for
    dot_bracket : str
        Dot-bracket notation string for RNA secondary structure.
        Can also be a file path to a file containing the dot-bracket string.
    skip_positions : list of int, optional
        1-based positions to exclude from AUROC calculation (e.g., primer binding sites)
    min_positions : int, optional
        Minimum number of valid positions required to calculate AUROC.
        If None, defaults to int(0.9 * len(dot_bracket)).
        Cells with fewer valid positions get NaN.
    save_as : str, optional
        Column name to save AUROC values in cells metadata.
        If None, uses f'AUROC_{gene}' (default: None)
    verbose : bool
        Whether to print progress information (default: True)

    Returns:
    --------
    None
        AUROC values are saved to cells metadata.

    Raises:
    -------
    ValueError
        If reactivity matrix not found (call calculate_reactivity first)
        If gene not found in genepos
        If dot_bracket length doesn't match gene positions

    Example:
    --------
    >>> # Calculate AUROC for 18S rRNA
    >>> dot_bracket_18S = "(((...)))..."  # or path to file
    >>> calculate_auroc(data, '18S', dot_bracket_18S)
    >>> print(data.cells[['leiden', 'AUROC_18S']].head())
    """
    import os

    # Validate inputs
    if not hasattr(data, 'reactivity') or data.reactivity is None:
        raise ValueError(
            "Reactivity matrix not found. Call calculate_reactivity() first."
        )

    if not isinstance(data.genepos, pd.DataFrame):
        raise ValueError("genepos must be a DataFrame with 'gene' column.")

    if 'gene' not in data.genepos.columns:
        raise ValueError("genepos DataFrame must have a 'gene' column.")

    # Load dot_bracket from file if it's a path
    if os.path.isfile(dot_bracket):
        with open(dot_bracket, 'r') as f:
            dot_bracket = f.read().strip()
        if verbose:
            print(f"Loaded dot-bracket structure from file")

    # Set default min_positions based on dot_bracket length
    if min_positions is None:
        min_positions = int(0.9 * len(dot_bracket))
        if verbose:
            print(f"Using min_positions = {min_positions} (90% of {len(dot_bracket)} positions)")

    # Get positions for the specified gene
    gene_mask = data.genepos['gene'] == gene
    gene_positions = np.where(gene_mask)[0]

    if len(gene_positions) == 0:
        raise ValueError(
            f"Gene '{gene}' not found in genepos. "
            f"Available genes: {data.genepos['gene'].unique().tolist()}"
        )

    # Get position numbers from genepos
    gene_genepos = data.genepos[gene_mask].copy()

    if verbose:
        print(f"Calculating AUROC for gene: {gene}")
        print(f"Number of positions for {gene}: {len(gene_positions)}")
        print(f"Dot-bracket length: {len(dot_bracket)}")

    # Create dot-bracket dataframe with position labels
    dotbracket_df = pd.DataFrame({
        "pos": range(1, len(dot_bracket) + 1),
        "dbrt": list(dot_bracket),
    })

    # Assign labels: 1.0 for '.', 0.0 for '(' or ')'
    dotbracket_df["dbrt_label"] = dotbracket_df["dbrt"].apply(
        lambda x: 1.0 if x == '.' else 0.0
    )

    # Merge gene positions with dot-bracket labels
    gene_genepos = gene_genepos.reset_index(drop=True)
    gene_genepos['_row_idx'] = gene_positions  # Store original row indices

    merged_df = gene_genepos.merge(
        dotbracket_df,
        left_on='pos',
        right_on='pos',
        how='inner'
    )

    if len(merged_df) == 0:
        raise ValueError(
            f"No matching positions found between genepos and dot-bracket. "
            f"Check that position numbering is consistent."
        )

    if verbose:
        print(f"Matched positions: {len(merged_df)}")

    # Apply skip_positions filter
    if skip_positions is not None:
        skip_set = set(skip_positions)
        merged_df = merged_df[~merged_df['pos'].isin(skip_set)]
        if verbose:
            print(f"Positions after skipping: {len(merged_df)}")

    if len(merged_df) < min_positions:
        raise ValueError(
            f"Only {len(merged_df)} positions available after filtering. "
            f"Need at least {min_positions}."
        )

    # Get the row indices in the reactivity matrix
    valid_row_indices = merged_df['_row_idx'].values
    dbrt_labels = merged_df['dbrt_label'].values

    # Extract reactivity submatrix for valid positions
    gene_reactivity = data.reactivity[valid_row_indices, :]

    # Get cell names
    cell_names = data.cell_names
    n_cells = len(cell_names)

    if verbose:
        print(f"Calculating AUROC for {n_cells} cells...")

    # Calculate AUROC for each cell
    auroc_values = []
    valid_count = 0

    for cell_idx in range(n_cells):
        cell_reactivity = gene_reactivity[:, cell_idx]

        # Get non-NaN positions
        valid_mask = ~np.isnan(cell_reactivity)
        n_valid = np.sum(valid_mask)

        if n_valid >= min_positions:
            pred = cell_reactivity[valid_mask]
            label = dbrt_labels[valid_mask]

            # Check if we have both classes
            unique_labels = np.unique(label)
            if len(unique_labels) < 2:
                auroc_values.append(np.nan)
            else:
                try:
                    fpr, tpr, _ = metrics.roc_curve(label, pred)
                    auc = metrics.auc(fpr, tpr)
                    auroc_values.append(round(auc, 3))
                    valid_count += 1
                except Exception:
                    auroc_values.append(np.nan)
        else:
            auroc_values.append(np.nan)

    # Create series with cell names as index
    auroc_series = pd.Series(auroc_values, index=cell_names, name=f'AUROC_{gene}')

    if verbose:
        print(f"\n{'='*60}")
        print(f"AUROC calculation complete:")
        print(f"  Cells with valid AUROC: {valid_count} / {n_cells}")
        valid_aurocs = auroc_series.dropna()
        if len(valid_aurocs) > 0:
            print(f"  Mean AUROC: {valid_aurocs.mean():.3f}")
            print(f"  Median AUROC: {valid_aurocs.median():.3f}")
            print(f"  AUROC range: [{valid_aurocs.min():.3f}, {valid_aurocs.max():.3f}]")

    # Save to cells metadata
    if save_as is None:
        save_as = f'AUROC_{gene}'

    if data.cells_is_df:
        data.cells[save_as] = data.cells.index.map(auroc_series)
        if verbose:
            print(f"\nSaved AUROC values to cells['{save_as}']")
    else:
        if verbose:
            print(f"\nNote: cells is a list, not DataFrame. AUROC values not saved to metadata.")
            print(f"Use data.to_cells_df() to convert cells to DataFrame first.")


def differential_reactivity_lm(
    data: 'ShapeData',
    predictors: Union[str, List[str]],
    gene: Optional[Union[str, List[str]]] = None,
    use_normalized: bool = False,
    cluster_column: Optional[str] = None,
    cluster_value: Optional[str] = None,
    min_cells: int = 10,
    correct_multiple: bool = True,
    correction_method: str = 'fdr_bh',
    verbose: bool = True
) -> pd.DataFrame:
    """
    Perform linear model analysis for differential reactivity across conditions.

    Fits a linear model: reactivity ~ b_0 + b_1*X1 + b_2*X2 + ...
    for each position. By default tests all positions; optionally filter by gene.

    Parameters:
    -----------
    data : ShapeData
        The ShapeData object with reactivity matrix and cell metadata
    predictors : str or list of str
        Column name(s) in cells metadata to use as predictors.
        Can be a single string or a list of strings.
        Examples:
        - 'diet' for reactivity ~ diet
        - ['diet', 'region'] for reactivity ~ diet + region
        - ['treatment', 'batch', 'age'] for multiple predictors
        Categorical variables are automatically dummy-coded.
    gene : str, list of str, or None
        Gene name(s) to analyze. If None (default), tests all positions.
        Can be a single gene name or a list of gene names.
    use_normalized : bool
        If True, use normalized_reactivity instead of raw reactivity (default: False)
    cluster_column : str, optional
        Column name containing cluster assignments. If provided with cluster_value,
        only cells belonging to that cluster are analyzed.
    cluster_value : str, optional
        Cluster value to filter on. Requires cluster_column to be set.
    min_cells : int
        Minimum number of cells with valid reactivity required per position
        (default: 10)
    correct_multiple : bool
        Whether to apply multiple testing correction (default: True)
    correction_method : str
        Method for multiple testing correction. Options: 'bonferroni', 'fdr_bh',
        'fdr_by', 'holm', 'sidak', 'holm-sidak' (default: 'fdr_bh')
    verbose : bool
        Whether to print progress information (default: True)

    Returns:
    --------
    pd.DataFrame
        Results dataframe with columns:
        - gene: gene name
        - pos: position number
        - variable: predictor variable name
        - coefficient: estimated coefficient (b)
        - std_err: standard error of coefficient
        - t_statistic: t-statistic
        - pvalue: p-value
        - pvalue_adj: adjusted p-value (if correct_multiple=True)
        - n_cells: number of cells used in the model
        - r_squared: R-squared of the model

    Example:
    --------
    >>> # Test all positions with reactivity ~ diet
    >>> results = differential_reactivity_lm(data, predictors='diet')

    >>> # Use normalized reactivity
    >>> results = differential_reactivity_lm(data, predictors='diet', use_normalized=True)

    >>> # Test specific gene
    >>> results = differential_reactivity_lm(data, predictors='diet', gene='18S')

    >>> # Test multiple genes
    >>> results = differential_reactivity_lm(data, predictors=['diet', 'region'], gene=['18S', '28S'])

    >>> # Test within a specific cluster
    >>> results = differential_reactivity_lm(
    ...     data, predictors='treatment',
    ...     cluster_column='leiden', cluster_value='0'
    ... )

    >>> # Filter significant results
    >>> sig_results = results[results['pvalue_adj'] < 0.05]
    """
    import statsmodels.api as sm
    from statsmodels.stats.multitest import multipletests

    # Convert single predictor to list
    if isinstance(predictors, str):
        predictors = [predictors]

    # Select which reactivity matrix to use
    if use_normalized:
        if not hasattr(data, 'normalized_reactivity') or data.normalized_reactivity is None:
            raise ValueError(
                "Normalized reactivity matrix not found. "
                "Call calculate_reactivity() with normalize=True first."
            )
        reactivity_matrix = data.normalized_reactivity
        reactivity_type = "normalized_reactivity"
    else:
        if not hasattr(data, 'reactivity') or data.reactivity is None:
            raise ValueError(
                "Reactivity matrix not found. Call calculate_reactivity() first."
            )
        reactivity_matrix = data.reactivity
        reactivity_type = "reactivity"

    # Validate inputs
    if not data.cells_is_df:
        raise ValueError("cells must be a DataFrame with metadata columns.")

    if not isinstance(data.genepos, pd.DataFrame):
        raise ValueError("genepos must be a DataFrame with 'gene' and 'pos' columns.")

    # Validate predictor variables exist in cells metadata
    for var in predictors:
        if var not in data.cells.columns:
            raise ValueError(
                f"Predictor '{var}' not found in cells metadata. "
                f"Available columns: {list(data.cells.columns)}"
            )

    # Get positions to test
    if gene is None:
        # Test all positions
        test_positions = np.arange(data.n_positions)
        test_genepos = data.genepos.copy()
        gene_label = "all genes"
    else:
        # Convert single gene to list
        if isinstance(gene, str):
            gene = [gene]
        # Filter to specified genes
        gene_mask = data.genepos['gene'].isin(gene)
        test_positions = np.where(gene_mask)[0]

        if len(test_positions) == 0:
            raise ValueError(
                f"Gene(s) {gene} not found in genepos. "
                f"Available genes: {data.genepos['gene'].unique().tolist()}"
            )

        test_genepos = data.genepos[gene_mask].copy()
        gene_label = ", ".join(gene) if len(gene) <= 3 else f"{len(gene)} genes"

    if verbose:
        print(f"Differential reactivity analysis (Linear Model)")
        print(f"{'='*60}")
        print(f"Using: {reactivity_type}")
        print(f"Genes: {gene_label}")
        print(f"Number of positions: {len(test_positions)}")
        print(f"Formula: {reactivity_type} ~ {' + '.join(predictors)}")

    # Filter cells by cluster if specified
    cell_mask = pd.Series(True, index=data.cells.index)

    if cluster_column is not None and cluster_value is not None:
        if cluster_column not in data.cells.columns:
            raise ValueError(f"Cluster column '{cluster_column}' not found in cells metadata.")
        if isinstance(cluster_value, list):
            cell_mask = data.cells[cluster_column].isin(cluster_value)
        else:
            cell_mask = data.cells[cluster_column] == cluster_value
        if verbose:
            print(f"Filtering to cluster: {cluster_column} in {cluster_value}" if isinstance(cluster_value, list)
                  else f"Filtering to cluster: {cluster_column} == {cluster_value}")
    elif cluster_column is not None:
        # Filter to cells that have non-null cluster values
        cell_mask = data.cells[cluster_column].notna()
        if verbose:
            print(f"Filtering to cells with valid {cluster_column}")

    valid_cell_indices = np.where(cell_mask)[0]
    cells_df = data.cells.iloc[valid_cell_indices].copy()

    if verbose:
        print(f"Cells in analysis: {len(valid_cell_indices)}")

    # Prepare design matrix (X) from cell metadata
    # Handle categorical variables using pd.get_dummies
    X_parts = []
    var_names = []

    for var in predictors:
        col = cells_df[var]
        if col.dtype == 'object' or col.dtype.name == 'category':
            # Categorical: create dummy variables (drop first level as reference)
            dummies = pd.get_dummies(col, prefix=var, drop_first=True, dtype=float)
            X_parts.append(dummies)
            var_names.extend(dummies.columns.tolist())
        else:
            # Numeric: use as is
            X_parts.append(col.to_frame())
            var_names.append(var)

    X = pd.concat(X_parts, axis=1)
    X = sm.add_constant(X)  # Add intercept

    if verbose:
        print(f"Design matrix shape: {X.shape}")
        print(f"Variables in model: {var_names}")

    # Extract reactivity for test positions and selected cells
    test_reactivity = reactivity_matrix[test_positions, :][:, valid_cell_indices]

    # Run linear model for each position
    results = []
    n_tested = 0

    for i, pos_idx in enumerate(test_positions):
        row = test_genepos.iloc[i]
        gene_name = row['gene']
        pos = row['pos']
        y = test_reactivity[i, :]

        # Filter out NaN values
        valid_mask = ~np.isnan(y)
        n_valid = np.sum(valid_mask)

        if n_valid < min_cells:
            continue

        y_valid = y[valid_mask]
        X_valid = X.iloc[valid_mask]

        # Check for NaN in X and drop
        complete_mask = ~X_valid.isna().any(axis=1)
        if complete_mask.sum() < min_cells:
            continue

        y_final = y_valid[complete_mask.values]
        X_final = X_valid[complete_mask]

        try:
            model = sm.OLS(y_final, X_final).fit()
            n_tested += 1

            # Extract coefficients for each variable (skip intercept)
            for var_name in var_names:
                if var_name in model.params.index:
                    results.append({
                        'gene': gene_name,
                        'pos': pos,
                        'variable': var_name,
                        'coefficient': model.params[var_name],
                        'std_err': model.bse[var_name],
                        't_statistic': model.tvalues[var_name],
                        'pvalue': model.pvalues[var_name],
                        'n_cells': len(y_final),
                        'r_squared': model.rsquared
                    })
        except Exception as e:
            if verbose:
                print(f"  Warning: Model failed at {gene_name}:{pos}: {e}")
            continue

    if verbose:
        print(f"\nPositions tested: {n_tested} / {len(test_positions)}")

    if len(results) == 0:
        if verbose:
            print("No results to return.")
        return pd.DataFrame()

    results_df = pd.DataFrame(results)

    # Apply multiple testing correction
    if correct_multiple and len(results_df) > 0:
        results_df['pvalue_adj'] = np.nan
        # Correct within each variable, skipping NaN p-values
        for var in results_df['variable'].unique():
            var_mask = results_df['variable'] == var
            pvals = results_df.loc[var_mask, 'pvalue'].values
            valid = ~np.isnan(pvals)

            if valid.sum() > 0:
                _, pvals_adj, _, _ = multipletests(pvals[valid], method=correction_method)
                adj_full = np.full(len(pvals), np.nan)
                adj_full[valid] = pvals_adj
                results_df.loc[var_mask, 'pvalue_adj'] = adj_full

        if verbose:
            for var in results_df['variable'].unique():
                var_results = results_df[results_df['variable'] == var]
                n_sig = (var_results['pvalue_adj'] < 0.05).sum()
                print(f"  {var}: {n_sig} significant positions (adj. p < 0.05)")

    return results_df


def differential_reactivity_wilcoxon(
    data: 'ShapeData',
    group_column: str,
    group1: Union[str, List[str]],
    group2: Union[str, List[str]],
    gene: Optional[Union[str, List[str]]] = None,
    use_normalized: bool = False,
    min_cells_per_group: int = 5,
    correct_multiple: bool = True,
    correction_method: str = 'fdr_bh',
    alternative: str = 'two-sided',
    verbose: bool = True
) -> pd.DataFrame:
    """
    Perform Wilcoxon rank-sum test (Mann-Whitney U) for differential reactivity.

    Tests whether reactivity differs between two groups at each position.
    By default tests all positions; optionally filter by gene.

    Parameters:
    -----------
    data : ShapeData
        The ShapeData object with reactivity matrix and cell metadata
    group_column : str
        Column name in cells metadata defining groups (e.g., 'diet', 'region')
    group1 : str or list of str
        Value(s) in group_column for first group
    group2 : str or list of str
        Value(s) in group_column for second group
    gene : str, list of str, or None
        Gene name(s) to analyze. If None (default), tests all positions.
        Can be a single gene name or a list of gene names.
    use_normalized : bool
        If True, use normalized_reactivity instead of raw reactivity (default: False)
    min_cells_per_group : int
        Minimum number of cells with valid reactivity required per group
        (default: 5)
    correct_multiple : bool
        Whether to apply multiple testing correction (default: True)
    correction_method : str
        Method for multiple testing correction (default: 'fdr_bh')
    alternative : str
        Alternative hypothesis: 'two-sided', 'less', or 'greater' (default: 'two-sided')
    verbose : bool
        Whether to print progress information (default: True)

    Returns:
    --------
    pd.DataFrame
        Results dataframe with columns:
        - gene: gene name
        - pos: position number
        - statistic: Wilcoxon test statistic (U)
        - pvalue: p-value
        - pvalue_adj: adjusted p-value (if correct_multiple=True)
        - effect_size: rank-biserial correlation (effect size)
        - mean_group1: mean reactivity in group1
        - mean_group2: mean reactivity in group2
        - diff_mean: difference in means (group1 - group2)
        - log2_fold_change: log2(mean_group1 / mean_group2)
        - n_group1: number of cells in group1
        - n_group2: number of cells in group2

    Example:
    --------
    >>> # Compare HFD vs CTRL across all positions
    >>> results = differential_reactivity_wilcoxon(
    ...     data, group_column='diet', group1='HFD', group2='CTRL'
    ... )

    >>> # Use normalized reactivity
    >>> results = differential_reactivity_wilcoxon(
    ...     data, group_column='diet', group1='HFD', group2='CTRL', use_normalized=True
    ... )

    >>> # Compare for specific gene
    >>> results = differential_reactivity_wilcoxon(
    ...     data, group_column='diet', group1='HFD', group2='CTRL', gene='18S'
    ... )

    >>> # Compare for multiple genes
    >>> results = differential_reactivity_wilcoxon(
    ...     data, group_column='region', group1='hypothalamus', group2='cortex',
    ...     gene=['18S', '28S']
    ... )

    >>> # Filter significant positions
    >>> sig_results = results[results['pvalue_adj'] < 0.05]
    """
    from scipy.stats import mannwhitneyu
    from statsmodels.stats.multitest import multipletests

    # Select which reactivity matrix to use
    if use_normalized:
        if not hasattr(data, 'normalized_reactivity') or data.normalized_reactivity is None:
            raise ValueError(
                "Normalized reactivity matrix not found. "
                "Call calculate_reactivity() with normalize=True first."
            )
        reactivity_matrix = data.normalized_reactivity
        reactivity_type = "normalized_reactivity"
    else:
        if not hasattr(data, 'reactivity') or data.reactivity is None:
            raise ValueError(
                "Reactivity matrix not found. Call calculate_reactivity() first."
            )
        reactivity_matrix = data.reactivity
        reactivity_type = "reactivity"

    # Validate inputs
    if not data.cells_is_df:
        raise ValueError("cells must be a DataFrame with metadata columns.")

    if not isinstance(data.genepos, pd.DataFrame):
        raise ValueError("genepos must be a DataFrame with 'gene' and 'pos' columns.")

    if group_column not in data.cells.columns:
        raise ValueError(
            f"Column '{group_column}' not found in cells metadata. "
            f"Available columns: {list(data.cells.columns)}"
        )

    # Convert group1/group2 to lists if single values
    if isinstance(group1, str):
        group1 = [group1]
    if isinstance(group2, str):
        group2 = [group2]

    # Get positions to test
    if gene is None:
        # Test all positions
        test_positions = np.arange(data.n_positions)
        test_genepos = data.genepos.copy()
        gene_label = "all genes"
    else:
        # Convert single gene to list
        if isinstance(gene, str):
            gene = [gene]
        # Filter to specified genes
        gene_mask = data.genepos['gene'].isin(gene)
        test_positions = np.where(gene_mask)[0]

        if len(test_positions) == 0:
            raise ValueError(
                f"Gene(s) {gene} not found in genepos. "
                f"Available genes: {data.genepos['gene'].unique().tolist()}"
            )

        test_genepos = data.genepos[gene_mask].copy()
        gene_label = ", ".join(gene) if len(gene) <= 3 else f"{len(gene)} genes"

    if verbose:
        print(f"Differential reactivity analysis (Wilcoxon Rank-Sum Test)")
        print(f"{'='*60}")
        print(f"Using: {reactivity_type}")
        print(f"Genes: {gene_label}")
        print(f"Number of positions: {len(test_positions)}")
        print(f"Group column: {group_column}")
        print(f"Group 1: {group1}")
        print(f"Group 2: {group2}")

    # Get cell indices for each group
    group1_mask = data.cells[group_column].isin(group1)
    group2_mask = data.cells[group_column].isin(group2)

    group1_indices = np.where(group1_mask)[0]
    group2_indices = np.where(group2_mask)[0]

    if verbose:
        print(f"Cells in group 1: {len(group1_indices)}")
        print(f"Cells in group 2: {len(group2_indices)}")

    if len(group1_indices) < min_cells_per_group:
        raise ValueError(f"Group 1 has fewer than {min_cells_per_group} cells.")
    if len(group2_indices) < min_cells_per_group:
        raise ValueError(f"Group 2 has fewer than {min_cells_per_group} cells.")

    # Extract reactivity for test positions
    test_reactivity = reactivity_matrix[test_positions, :]

    # Run Wilcoxon test for each position
    results = []
    n_tested = 0

    for i, pos_idx in enumerate(test_positions):
        row = test_genepos.iloc[i]
        gene_name = row['gene']
        pos = row['pos']

        # Get reactivity values for each group
        vals1 = test_reactivity[i, group1_indices]
        vals2 = test_reactivity[i, group2_indices]

        # Remove NaN values
        vals1_valid = vals1[~np.isnan(vals1)]
        vals2_valid = vals2[~np.isnan(vals2)]

        # Check minimum cells
        if len(vals1_valid) < min_cells_per_group or len(vals2_valid) < min_cells_per_group:
            continue

        n_tested += 1

        try:
            # Perform Wilcoxon rank-sum test (Mann-Whitney U)
            stat, pval = mannwhitneyu(vals1_valid, vals2_valid, alternative=alternative)

            # Calculate effect size (rank-biserial correlation)
            n1, n2 = len(vals1_valid), len(vals2_valid)
            effect_size = 1 - (2 * stat) / (n1 * n2)

            # Calculate means and log2 fold change
            mean1 = np.mean(vals1_valid)
            mean2 = np.mean(vals2_valid)

            # Handle log2 fold change (avoid division by zero or negative values)
            if mean2 != 0 and mean1 > 0 and mean2 > 0:
                log2fc = np.log2(mean1 / mean2)
            elif mean1 > 0 and mean2 <= 0:
                log2fc = np.inf
            elif mean1 <= 0 and mean2 > 0:
                log2fc = -np.inf
            else:
                log2fc = np.nan

            results.append({
                'gene': gene_name,
                'pos': pos,
                'statistic': stat,
                'pvalue': pval,
                'effect_size': effect_size,
                'mean_group1': mean1,
                'mean_group2': mean2,
                'diff_mean': mean1 - mean2,
                'log2_fold_change': log2fc,
                'n_group1': n1,
                'n_group2': n2
            })
        except Exception as e:
            if verbose:
                print(f"  Warning: Test failed at {gene_name}:{pos}: {e}")
            continue

    if verbose:
        print(f"\nPositions tested: {n_tested} / {len(test_positions)}")

    if len(results) == 0:
        if verbose:
            print("No results to return.")
        return pd.DataFrame()

    results_df = pd.DataFrame(results)

    # Apply multiple testing correction
    if correct_multiple and len(results_df) > 0:
        pvals = results_df['pvalue'].values
        _, pvals_adj, _, _ = multipletests(pvals, method=correction_method)
        results_df['pvalue_adj'] = pvals_adj

        if verbose:
            n_sig = (results_df['pvalue_adj'] < 0.05).sum()
            print(f"Significant positions (adj. p < 0.05): {n_sig}")

    # Sort by p-value
    results_df = results_df.sort_values('pvalue').reset_index(drop=True)

    return results_df
