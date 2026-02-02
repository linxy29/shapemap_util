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

    Raises:
    -------
    ValueError
        If cells is not a DataFrame or cluster_column not found
        If genepos is not a DataFrame

    Example:
    --------
    >>> # First add control data
    >>> data.add_control_data(dmso_df)
    >>> # Then calculate reactivity
    >>> calculate_reactivity(data, cluster_column='leiden')
    >>> # Access stored reactivity
    >>> print(data.reactivity.shape)
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

    # Optionally drop mutrate to free memory
    if drop_mutrate:
        data.mutrate = None
        if verbose:
            print("Dropped mutrate matrix to free memory")


def calculate_cell_correlation(
    data: 'ShapeData',
    gene: str,
    cluster_column: str = 'leiden',
    method: str = 'pearson',
    min_shared_positions: int = 100,
    save_mean_as: Optional[str] = None,
    verbose: bool = True
) -> pd.DataFrame:
    """
    Calculate pairwise correlation between cells within the same cluster for a given gene.

    For each pair of cells in the same cluster, computes correlation of their
    reactivity values across all positions of the specified gene.

    Parameters:
    -----------
    data : ShapeData
        The ShapeData object to analyze
    gene : str
        Gene name to calculate correlations for
    cluster_column : str
        Column name in cells DataFrame containing cluster assignments (default: 'leiden')
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
        - cluster: cluster name
        - correlation: pairwise correlation value
        - n_shared: number of shared non-NaN positions

    Raises:
    -------
    ValueError
        If reactivity matrix not found (call calculate_reactivity first)
        If gene not found in genepos
        If cells is not a DataFrame

    Example:
    --------
    >>> # Calculate correlations for 18S gene
    >>> corr_df = calculate_cell_correlation(data, '18S', cluster_column='leiden')
    >>> print(corr_df.head())
    >>> # Mean correlation is saved to data.cells['mean_corr_18S']
    >>> print(data.cells[['leiden', 'mean_corr_18S']].head())
    """
    from itertools import combinations

    # Validate inputs
    if not hasattr(data, 'reactivity') or data.reactivity is None:
        raise ValueError(
            "Reactivity matrix not found. Call calculate_reactivity() first."
        )

    if not data.cells_is_df:
        raise ValueError(
            "cells must be a DataFrame with cluster assignments."
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

    # Get cell names and cluster assignments
    cell_names = data.cell_names
    clusters = data.cells[cluster_column].dropna().unique()

    if verbose:
        print(f"Clusters: {list(clusters)}")

    # Store correlation results
    corr_results = []

    for cluster in clusters:
        if verbose:
            print(f"\nProcessing cluster: {cluster}")

        # Get cell indices for this cluster
        cluster_mask = data.cells[cluster_column] == cluster
        cluster_cell_indices = np.where(cluster_mask)[0]
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
    data.cells[save_mean_as] = data.cells.index.map(cell_mean_corr)

    if verbose:
        print(f"\nMean correlation saved to cells['{save_mean_as}']")
        valid_means = data.cells[save_mean_as].dropna()
        print(f"Cells with valid mean correlation: {len(valid_means)} / {len(data.cells)}")

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

    # Mean mutation rate per position
    mut_sum = np.array(data.mutrate.sum(axis=1)).flatten()
    mean_mutrate = np.divide(mut_sum, n_cells, where=n_cells > 0, out=np.zeros_like(mut_sum, dtype=float))

    # Add to genepos metadata
    if isinstance(data.genepos, pd.DataFrame):
        data.genepos['n_cells'] = n_cells
        data.genepos['mean_coverage'] = mean_coverage
        data.genepos['mean_mutrate'] = mean_mutrate
        if verbose:
            print(f"Added columns to genepos metadata: n_cells, mean_coverage, mean_mutrate")
    else:
        # Convert to DataFrame first
        data.genepos = pd.DataFrame({
            'position': data.genepos,
            'n_cells': n_cells,
            'mean_coverage': mean_coverage,
            'mean_mutrate': mean_mutrate
        })
        if verbose:
            print(f"Converted genepos to DataFrame and added: n_cells, mean_coverage, mean_mutrate")


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
    gene: str,
    predictors: Union[str, List[str]],
    cluster_column: Optional[str] = None,
    cluster_value: Optional[str] = None,
    min_cells: int = 10,
    correct_multiple: bool = True,
    correction_method: str = 'fdr_bh',
    verbose: bool = True
) -> pd.DataFrame:
    """
    Perform linear model analysis for differential reactivity across conditions.

    Fits a linear model: normalized_reactivity ~ b_0 + b_1*X1 + b_2*X2 + ...
    for each position in the specified gene. Predictors can be any column(s)
    from the cells metadata.

    Parameters:
    -----------
    data : ShapeData
        The ShapeData object with reactivity matrix and cell metadata
    gene : str
        Gene name to analyze
    predictors : str or list of str
        Column name(s) in cells metadata to use as predictors.
        Can be a single string or a list of strings.
        Examples:
        - 'diet' for reactivity ~ diet
        - ['diet', 'region'] for reactivity ~ diet + region
        - ['treatment', 'batch', 'age'] for multiple predictors
        Categorical variables are automatically dummy-coded.
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
    >>> # Test reactivity ~ diet for 18S
    >>> results = differential_reactivity_lm(data, gene='18S', predictors='diet')

    >>> # Test reactivity ~ diet + region
    >>> results = differential_reactivity_lm(
    ...     data, gene='18S', predictors=['diet', 'region']
    ... )

    >>> # Test within a specific cluster
    >>> results = differential_reactivity_lm(
    ...     data, gene='18S', predictors='treatment',
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

    # Validate inputs
    if not hasattr(data, 'reactivity') or data.reactivity is None:
        raise ValueError(
            "Reactivity matrix not found. Call calculate_reactivity() first."
        )

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

    # Get positions for the specified gene
    gene_mask = data.genepos['gene'] == gene
    gene_positions = np.where(gene_mask)[0]

    if len(gene_positions) == 0:
        raise ValueError(
            f"Gene '{gene}' not found in genepos. "
            f"Available genes: {data.genepos['gene'].unique().tolist()}"
        )

    gene_genepos = data.genepos[gene_mask].copy()

    if verbose:
        print(f"Differential reactivity analysis (Linear Model)")
        print(f"{'='*60}")
        print(f"Gene: {gene}")
        print(f"Number of positions: {len(gene_positions)}")
        print(f"Formula: reactivity ~ {' + '.join(predictors)}")

    # Filter cells by cluster if specified
    cell_mask = pd.Series(True, index=data.cells.index)

    if cluster_column is not None and cluster_value is not None:
        if cluster_column not in data.cells.columns:
            raise ValueError(f"Cluster column '{cluster_column}' not found in cells metadata.")
        cell_mask = data.cells[cluster_column] == cluster_value
        if verbose:
            print(f"Filtering to cluster: {cluster_column} == {cluster_value}")
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
            dummies = pd.get_dummies(col, prefix=var, drop_first=True)
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

    # Extract reactivity for gene positions and selected cells
    gene_reactivity = data.reactivity[gene_positions, :][:, valid_cell_indices]

    # Run linear model for each position
    results = []
    n_tested = 0

    for i, pos_idx in enumerate(gene_positions):
        pos = gene_genepos.iloc[i]['pos']
        y = gene_reactivity[i, :]

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
                        'gene': gene,
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
                print(f"  Warning: Model failed at pos {pos}: {e}")
            continue

    if verbose:
        print(f"\nPositions tested: {n_tested} / {len(gene_positions)}")

    if len(results) == 0:
        if verbose:
            print("No results to return.")
        return pd.DataFrame()

    results_df = pd.DataFrame(results)

    # Apply multiple testing correction
    if correct_multiple and len(results_df) > 0:
        # Correct within each variable
        for var in results_df['variable'].unique():
            var_mask = results_df['variable'] == var
            pvals = results_df.loc[var_mask, 'pvalue'].values

            _, pvals_adj, _, _ = multipletests(pvals, method=correction_method)
            results_df.loc[var_mask, 'pvalue_adj'] = pvals_adj

        if verbose:
            for var in results_df['variable'].unique():
                var_results = results_df[results_df['variable'] == var]
                n_sig = (var_results['pvalue_adj'] < 0.05).sum()
                print(f"  {var}: {n_sig} significant positions (adj. p < 0.05)")

    return results_df


def differential_reactivity_wilcoxon(
    data: 'ShapeData',
    gene: str,
    group_column: str,
    group1: Union[str, List[str]],
    group2: Union[str, List[str]],
    min_cells_per_group: int = 5,
    correct_multiple: bool = True,
    correction_method: str = 'fdr_bh',
    alternative: str = 'two-sided',
    verbose: bool = True
) -> pd.DataFrame:
    """
    Perform Wilcoxon rank-sum test (Mann-Whitney U) for differential reactivity.

    Tests whether reactivity differs between two groups at each position.

    Parameters:
    -----------
    data : ShapeData
        The ShapeData object with reactivity matrix and cell metadata
    gene : str
        Gene name to analyze
    group_column : str
        Column name in cells metadata defining groups (e.g., 'diet', 'region')
    group1 : str or list of str
        Value(s) in group_column for first group
    group2 : str or list of str
        Value(s) in group_column for second group
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
        - log2_fold_change: log2(mean_group1 / mean_group2)
        - n_group1: number of cells in group1
        - n_group2: number of cells in group2

    Example:
    --------
    >>> # Compare HFD vs CTRL diet
    >>> results = differential_reactivity_wilcoxon(
    ...     data, gene='18S', group_column='diet',
    ...     group1='HFD', group2='CTRL'
    ... )
    >>> # Filter significant positions
    >>> sig_results = results[results['pvalue_adj'] < 0.05]
    """
    from scipy.stats import mannwhitneyu
    from statsmodels.stats.multitest import multipletests

    # Validate inputs
    if not hasattr(data, 'reactivity') or data.reactivity is None:
        raise ValueError(
            "Reactivity matrix not found. Call calculate_reactivity() first."
        )

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

    # Get positions for the specified gene
    gene_mask = data.genepos['gene'] == gene
    gene_positions = np.where(gene_mask)[0]

    if len(gene_positions) == 0:
        raise ValueError(
            f"Gene '{gene}' not found in genepos. "
            f"Available genes: {data.genepos['gene'].unique().tolist()}"
        )

    gene_genepos = data.genepos[gene_mask].copy()

    if verbose:
        print(f"Differential reactivity analysis (Wilcoxon Rank-Sum Test)")
        print(f"{'='*60}")
        print(f"Gene: {gene}")
        print(f"Number of positions: {len(gene_positions)}")
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

    # Extract reactivity for gene positions
    gene_reactivity = data.reactivity[gene_positions, :]

    # Run Wilcoxon test for each position
    results = []
    n_tested = 0

    for i, pos_idx in enumerate(gene_positions):
        pos = gene_genepos.iloc[i]['pos']

        # Get reactivity values for each group
        vals1 = gene_reactivity[i, group1_indices]
        vals2 = gene_reactivity[i, group2_indices]

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
                'gene': gene,
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
                print(f"  Warning: Test failed at pos {pos}: {e}")
            continue

    if verbose:
        print(f"\nPositions tested: {n_tested} / {len(gene_positions)}")

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
