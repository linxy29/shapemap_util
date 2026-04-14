"""
Metadata functions for ShapeData.

This module provides functions for adding and joining metadata
to ShapeData objects.
"""

import numpy as np
import pandas as pd
from typing import Union, Optional, List, TYPE_CHECKING

if TYPE_CHECKING:
    from .core import ShapeData


def add_genepos_counts(data: 'ShapeData', inplace: bool = True) -> Optional['ShapeData']:
    """
    Add 'n_genepos' column to cells metadata, counting non-zero entries per cell.

    Parameters:
    -----------
    data : ShapeData
        The ShapeData object to modify
    inplace : bool
        If True, modify in place. If False, return a new ShapeData.

    Returns:
    --------
    ShapeData or None
        Returns new object if inplace=False, otherwise None.
    """
    # Count non-zero entries per cell (column)
    n_genepos_per_cell = np.array((data.coverage > 0).sum(axis=0)).flatten()

    if inplace:
        # Convert to DataFrame if needed
        if not data.cells_is_df:
            data.cells = pd.DataFrame({'n_genepos': n_genepos_per_cell}, index=data.cells)
            data.cells.index.name = 'cell'
        else:
            data.cells['n_genepos'] = n_genepos_per_cell
        return None
    else:
        new_data = data.copy()
        add_genepos_counts(new_data, inplace=True)
        return new_data


def join_cell_metadata(
    data: 'ShapeData',
    metadata_df: pd.DataFrame,
    columns: Union[str, List[str]],
    on: str = 'cell',
    how: str = 'left',
    inplace: bool = True
) -> Optional['ShapeData']:
    """
    Join cell metadata from an external DataFrame.

    Parameters:
    -----------
    data : ShapeData
        The ShapeData object to modify
    metadata_df : pd.DataFrame
        DataFrame with metadata to join
    columns : str or list of str
        Column(s) to join from metadata_df
    on : str
        Column in metadata_df that contains cell identifiers (default: 'cell')
    how : str
        Join type (default: 'left')
    inplace : bool
        If True, modify in place

    Returns:
    --------
    ShapeData or None
    """
    if isinstance(columns, str):
        columns = [columns]

    # Ensure cells is a DataFrame
    if inplace:
        if not data.cells_is_df:
            data.cells = pd.DataFrame(index=data.cells)
            data.cells.index.name = 'cell'

        # Join the columns
        data.cells = data.cells.join(
            metadata_df.set_index(on)[columns],
            how=how
        )
        return None
    else:
        new_data = data.to_cells_df()
        join_cell_metadata(new_data, metadata_df, columns, on, how, inplace=True)
        return new_data


def join_genepos_metadata(
    data: 'ShapeData',
    metadata_df: pd.DataFrame,
    columns: Union[str, List[str]],
    on: Union[str, List[str]] = None,
    how: str = 'left',
    inplace: bool = True
) -> Optional['ShapeData']:
    """
    Join position metadata from an external DataFrame.

    Parameters:
    -----------
    data : ShapeData
        The ShapeData object to modify
    metadata_df : pd.DataFrame
        DataFrame with metadata to join
    columns : str or list of str
        Column(s) to join from metadata_df
    on : str or list of str, optional
        Column(s) to join on (must exist in both genepos and metadata_df).
        Default: ['gene', 'pos']
    how : str
        Join type (default: 'left')
    inplace : bool
        If True, modify in place

    Returns:
    --------
    ShapeData or None

    Raises:
    -------
    ValueError
        If genepos is not a DataFrame
    """
    if on is None:
        on = ['gene', 'pos']

    if isinstance(columns, str):
        columns = [columns]

    if not isinstance(data.genepos, pd.DataFrame):
        raise ValueError(
            "genepos must be a DataFrame to join metadata. "
            "Use to_genepos_df() first or provide genepos as DataFrame during initialization."
        )

    if inplace:
        # Merge the columns
        data.genepos = data.genepos.merge(
            metadata_df[on + columns] if isinstance(on, list) else metadata_df[[on] + columns],
            on=on,
            how=how
        )
        return None
    else:
        new_data = data.copy()
        join_genepos_metadata(new_data, metadata_df, columns, on, how, inplace=True)
        return new_data


def add_control_data(
    data: 'ShapeData',
    control_df: pd.DataFrame,
    cluster_col: str = 'sample',
    cluster_prefix: str = 'dmso_',
    coverage_col: str = 'coverage',
    mutrate_col: str = 'mutrate',
    on: List[str] = None,
    control_suffix: str = 'control',
    inplace: bool = True,
    verbose: bool = True
) -> Optional['ShapeData']:
    """
    Add control (e.g., DMSO) coverage and mutation rate data to genepos metadata.

    This function processes control data where each row has a cluster identifier
    (extracted from cluster_col by removing cluster_prefix), and adds columns
    named {coverage_col}_{control_suffix}_{cluster} and {mutrate_col}_{control_suffix}_{cluster}
    for each cluster found in the control data.

    Parameters:
    -----------
    data : ShapeData
        The ShapeData object to modify
    control_df : pd.DataFrame
        DataFrame with control data containing columns for cluster identification,
        coverage, mutrate, and position identifiers (gene, pos)
    cluster_col : str
        Column name containing cluster identifiers (default: 'sample')
    cluster_prefix : str
        Prefix to remove from cluster_col values to get cluster names (default: 'dmso_')
    coverage_col : str
        Column name for coverage values (default: 'coverage')
    mutrate_col : str
        Column name for mutation rate values (default: 'mutrate')
    on : list of str, optional
        Columns to merge on (default: ['gene', 'pos'])
    control_suffix : str
        Suffix for control columns (default: 'control')
    inplace : bool
        If True, modify in place (default: True)
    verbose : bool
        If True, print progress information (default: True)

    Returns:
    --------
    ShapeData or None
        Returns new object if inplace=False, otherwise None.

    Raises:
    -------
    ValueError
        If genepos is not a DataFrame

    Example:
    --------
    >>> dmso_df = pd.read_csv("dmso_mutrate.csv")
    >>> # dmso_df has columns: sample, gene, pos, coverage, mutrate
    >>> # sample values like "dmso_K562", "dmso_HEK293T"
    >>> add_control_data(data, dmso_df)
    >>> # Now data.genepos has columns: coverage_control_K562, mutrate_control_K562, etc.
    """
    from .core import ShapeData

    if on is None:
        on = ['gene', 'pos']

    if not isinstance(data.genepos, pd.DataFrame):
        raise ValueError(
            "genepos must be a DataFrame to add control data. "
            "Use to_genepos_df() first or provide genepos as DataFrame during initialization."
        )

    # Work on a copy of control_df to avoid modifying the original
    control_df = control_df.copy()

    # Extract cluster names by removing prefix
    control_df['_cluster'] = control_df[cluster_col].str.replace(cluster_prefix, "", regex=False)
    clusters = control_df['_cluster'].unique()

    if verbose:
        print(f"Found {len(clusters)} clusters in control data: {list(clusters)}")

    if inplace:
        genepos_df = data.genepos
    else:
        genepos_df = data.genepos.copy()

    # Merge control data for each cluster
    for cluster in clusters:
        if verbose:
            print(f"\nProcessing cluster: {cluster}")

        # Filter control data for this cluster
        cluster_data = control_df[control_df['_cluster'] == cluster].copy()

        # Create renamed columns
        cov_col_new = f'{coverage_col}_{control_suffix}_{cluster}'
        mut_col_new = f'{mutrate_col}_{control_suffix}_{cluster}'

        cluster_data_renamed = cluster_data.rename(columns={
            coverage_col: cov_col_new,
            mutrate_col: mut_col_new
        })[on + [cov_col_new, mut_col_new]]

        # Merge with genepos
        genepos_df = genepos_df.merge(
            cluster_data_renamed,
            on=on,
            how='left'
        )

        # Report merge results
        if verbose:
            n_matched = genepos_df[mut_col_new].notna().sum()
            n_total = len(genepos_df)
            print(f"  Matched positions: {n_matched:,} / {n_total:,} ({100*n_matched/n_total:.2f}%)")

    if inplace:
        data.genepos = genepos_df
        return None
    else:
        return ShapeData(
            data.coverage.copy(),
            data.mutrate.copy() if data.mutrate is not None else None,
            genepos_df,
            data.cells.copy() if data.cells_is_df else list(data.cells),
            data.reactivity.copy() if data.reactivity is not None else None,
            data.normalized_reactivity.copy() if hasattr(data, 'normalized_reactivity') and data.normalized_reactivity is not None else None
        )
