"""
Metabin functions for ShapeData.

This module provides functions for grouping spatial bins into metabins
based on cluster membership and spatial proximity using recursive balanced bisection.
"""

import numpy as np
import pandas as pd
from scipy import sparse
from typing import Optional, Dict, List, TYPE_CHECKING

if TYPE_CHECKING:
    from .core import ShapeData


def _normalize_id(x):
    """Normalize bin ID to string, converting float-like values to int first."""
    try:
        return str(int(float(x)))
    except (ValueError, TypeError):
        return str(x)


def _balanced_bisection(coords, n_bins):
    """
    Recursively partition bins along the longer spatial axis into balanced groups.

    At each step, the bins are split along the axis with the larger range into
    two halves, with the split ratio matching the target group ratio so that
    each final group has approximately `n_bins` bins.

    Parameters:
    -----------
    coords : np.ndarray, shape (n, 2)
        (x, y) coordinates of bins.
    n_bins : int
        Target number of bins per group.

    Returns:
    --------
    np.ndarray
        Integer label array of length n, one label per bin.
    """
    n = len(coords)
    n_groups = max(1, round(n / n_bins))
    labels = np.zeros(n, dtype=int)

    def _recurse(indices, n_target_groups, label_offset):
        if n_target_groups <= 1:
            labels[indices] = label_offset
            return label_offset + 1

        sub_coords = coords[indices]
        # Pick the longer axis
        ranges = sub_coords.max(axis=0) - sub_coords.min(axis=0)
        axis = int(np.argmax(ranges))

        # Split proportionally to target group counts
        left_groups = n_target_groups // 2
        right_groups = n_target_groups - left_groups
        left_n = round(len(indices) * left_groups / n_target_groups)
        left_n = max(1, min(left_n, len(indices) - 1))

        order = np.argsort(sub_coords[:, axis])
        left_idx = indices[order[:left_n]]
        right_idx = indices[order[left_n:]]

        next_label = _recurse(left_idx, left_groups, label_offset)
        next_label = _recurse(right_idx, right_groups, next_label)
        return next_label

    _recurse(np.arange(n), n_groups, 0)
    return labels


def create_metabins_from_dict(
    target_data: 'ShapeData',
    metabin_dict: Dict[str, List],
    cell_id_col: str = 'bin100_cellID',
    extra_meta: Optional[Dict[str, dict]] = None,
    verbose: bool = True,
    ids_col_name: str = 'bin100_ids',
    count_col_name: str = 'n_bins',
    index_prefix: str = 'metabin',
) -> 'ShapeData':
    """
    Create metabins by aggregating bins according to a dictionary of groupings.

    For each metabin, the specified bins are aggregated: coverage is summed,
    mutation rate is recomputed as sum(mutcount)/sum(coverage) if mutcount is
    available, otherwise as a coverage-weighted average.

    Parameters:
    -----------
    target_data : ShapeData
        The ShapeData whose bins will be aggregated. Must have cells as DataFrame.
    metabin_dict : dict
        Dictionary mapping metabin name (str) to list of bin IDs.
        Bin IDs should match values in target_data.cells[cell_id_col] (or index).
    cell_id_col : str
        Column name for bin IDs in target_data cells (default: 'bin100_cellID').
        If this column is absent, the cells index is used.
    extra_meta : dict, optional
        Dictionary mapping metabin name to dict of extra metadata columns.
        These are merged into the output cells DataFrame for each metabin.
        Example: {'metabin_A': {'cluster': 'K562', 'x': 10.5, 'y': 20.3}}
    verbose : bool
        Print progress information (default: True).

    Returns:
    --------
    ShapeData
        New ShapeData with metabins as cells. The cells DataFrame always contains:
        - 'bin100_ids': list of matched bin IDs belonging to this metabin
        - 'n_bins': number of bins in this metabin
        Plus any columns from extra_meta.

    Example:
    --------
    >>> metabin_dict = {
    ...     'group_A': ['bin_1', 'bin_2', 'bin_3'],
    ...     'group_B': ['bin_4', 'bin_5'],
    ... }
    >>> result = data.create_metabins_from_dict(metabin_dict)
    """
    from .core import ShapeData

    if not target_data.cells_is_df:
        raise ValueError("target_data cells must be a DataFrame. Use to_cells_df() first.")

    # Build target bin ID -> column index lookup
    target_cells_df = target_data.cells
    if cell_id_col in target_cells_df.columns:
        target_bin_ids = target_cells_df[cell_id_col].apply(_normalize_id).values
    else:
        target_bin_ids = np.array([_normalize_id(x) for x in target_cells_df.index])
    target_id_to_idx = {bid: i for i, bid in enumerate(target_bin_ids)}

    has_mutcount = target_data.mutcount is not None
    has_mutrate = target_data.mutrate is not None
    n_positions = target_data.n_positions

    if verbose:
        print(f"Creating metabins from dict ({len(metabin_dict)} entries)")
        print(f"  Target bins: {len(target_cells_df)}")
        print(f"  Has mutcount: {has_mutcount}")

    metabin_cov_cols = []
    metabin_mut_cols = []
    metabin_mc_cols = []
    metabin_meta = []
    skipped = 0

    for metabin_name, bin_id_list in metabin_dict.items():
        # Normalize and look up column indices
        target_indices = []
        matched_ids = []
        for bid in bin_id_list:
            norm_bid = _normalize_id(bid)
            tidx = target_id_to_idx.get(norm_bid)
            if tidx is not None:
                target_indices.append(tidx)
                matched_ids.append(norm_bid)

        if len(target_indices) == 0:
            skipped += 1
            continue

        target_indices = np.array(target_indices)

        # Aggregate coverage: sum across bins
        group_cov = target_data.coverage[:, target_indices].sum(axis=1)
        group_cov = np.asarray(group_cov).flatten()

        # Aggregate mutcount if available
        if has_mutcount:
            group_mc = target_data.mutcount[:, target_indices].sum(axis=1)
            group_mc = np.asarray(group_mc).flatten()

        # Aggregate mutrate
        if has_mutrate:
            if has_mutcount:
                group_mr = np.zeros(n_positions, dtype=np.float32)
                nonzero = group_cov > 0
                group_mr[nonzero] = group_mc[nonzero] / group_cov[nonzero]
            else:
                group_mut_cov = target_data.mutrate[:, target_indices].multiply(
                    target_data.coverage[:, target_indices]
                ).sum(axis=1)
                group_mut_cov = np.asarray(group_mut_cov).flatten()
                group_mr = np.zeros(n_positions, dtype=np.float32)
                nonzero = group_cov > 0
                group_mr[nonzero] = group_mut_cov[nonzero] / group_cov[nonzero]

        # Store as sparse columns, ensuring mutrate has the same sparsity as coverage
        cov_col = sparse.csc_matrix(group_cov.reshape(-1, 1), dtype=np.int32)
        metabin_cov_cols.append(cov_col)
        if has_mutrate:
            cov_nz = cov_col.nonzero()[0]
            mr_col = sparse.csc_matrix(
                (group_mr[cov_nz], (cov_nz, np.zeros(len(cov_nz), dtype=int))),
                shape=(n_positions, 1), dtype=np.float32
            )
            metabin_mut_cols.append(mr_col)
        if has_mutcount:
            cov_nz = cov_col.nonzero()[0] if not has_mutrate else cov_nz
            mc_col = sparse.csc_matrix(
                (group_mc[cov_nz].astype(np.int32), (cov_nz, np.zeros(len(cov_nz), dtype=int))),
                shape=(n_positions, 1), dtype=np.int32
            )
            metabin_mc_cols.append(mc_col)

        # Build metadata
        meta = {
            ids_col_name: matched_ids,
            count_col_name: len(target_indices),
        }
        if extra_meta is not None and metabin_name in extra_meta:
            meta.update(extra_meta[metabin_name])
        metabin_meta.append(meta)

    if verbose:
        print(f"\nMetabins created: {len(metabin_meta)}")
        if skipped > 0:
            print(f"Skipped (no matched bins): {skipped}")

    if len(metabin_meta) == 0:
        raise ValueError("No metabins could be created - no bin ID matches found.")

    # Build output matrices
    metabin_cov = sparse.hstack(metabin_cov_cols, format='csr')
    metabin_mut = sparse.hstack(metabin_mut_cols, format='csr') if has_mutrate else None
    metabin_mc = sparse.hstack(metabin_mc_cols, format='csr') if has_mutcount else None

    # Build cells DataFrame
    metabin_cells = pd.DataFrame(metabin_meta)
    metabin_cells.index = [f"{index_prefix}_{i}" for i in range(len(metabin_cells))]
    metabin_cells.index.name = 'cell'

    # Copy genepos
    genepos_copy = target_data.genepos.copy() if isinstance(target_data.genepos, pd.DataFrame) else list(target_data.genepos)

    if verbose:
        print(f"Output shape: ({n_positions}, {len(metabin_cells)})")
        bins_per_metabin = metabin_cells[count_col_name]
        print(f"Bins per metabin: min={bins_per_metabin.min()}, "
              f"median={bins_per_metabin.median():.0f}, "
              f"max={bins_per_metabin.max()}")

    return ShapeData(
        metabin_cov,
        metabin_mut,
        genepos_copy,
        metabin_cells,
        mutcount_sparse=metabin_mc
    )


def create_metabins(
    data: 'ShapeData',
    n_bins: int,
    cluster_col: str = 'leiden',
    x_col: str = 'x',
    y_col: str = 'y',
    cell_id_col: str = 'bin100_cellID',
    verbose: bool = True
) -> 'ShapeData':
    """
    Group spatial bins from the same cluster into metabins using recursive balanced bisection.

    For each cluster, bins are recursively partitioned along the longer spatial axis
    into groups of approximately `n_bins` bins each. Coverage is summed, mutation rate
    is recomputed as a coverage-weighted average (or sum(mutcount)/sum(coverage) if
    mutcount is available).

    Parameters:
    -----------
    data : ShapeData
        The ShapeData object. Must have cells as DataFrame with cluster, x, y columns.
    n_bins : int
        Target number of bins per metabin.
    cluster_col : str
        Column name in cells metadata for cluster labels (default: 'leiden').
    x_col : str
        Column name for x coordinates (default: 'x').
    y_col : str
        Column name for y coordinates (default: 'y').
    cell_id_col : str
        Column name for bin100 cell IDs (default: 'bin100_cellID').
        If this column exists in cells metadata, it is used as the bin identifier.
        Otherwise, the cells index is used.
    verbose : bool
        Print progress information (default: True).

    Returns:
    --------
    ShapeData
        New ShapeData with metabins as cells. The cells DataFrame contains:
        - cluster_col: cluster label
        - x_col: x coordinate of the medoid bin
        - y_col: y coordinate of the medoid bin
        - 'bin100_ids': list of bin100 IDs belonging to this metabin
        - 'n_bins': number of bins in this metabin
    """
    from .core import ShapeData

    # Validate inputs
    if not data.cells_is_df:
        raise ValueError("cells must be a DataFrame. Use to_cells_df() first.")

    cells_df = data.cells
    for col in [cluster_col, x_col, y_col]:
        if col not in cells_df.columns:
            raise ValueError(f"Column '{col}' not found in cells metadata. "
                             f"Available columns: {list(cells_df.columns)}")

    has_cell_id_col = cell_id_col in cells_df.columns
    has_mutcount = data.mutcount is not None
    has_mutrate = data.mutrate is not None

    if verbose:
        print(f"Creating metabins with target size {n_bins}")
        print(f"  Cluster column: {cluster_col}")
        print(f"  Coordinates: ({x_col}, {y_col})")
        print(f"  Has mutcount: {has_mutcount}")

    clusters = cells_df[cluster_col].unique()
    n_positions = data.n_positions

    # Collect results per metabin
    metabin_cov_cols = []      # list of sparse column vectors
    metabin_mut_cols = []      # list of sparse column vectors (mutrate)
    metabin_mc_cols = []       # list of sparse column vectors (mutcount)
    metabin_meta = []          # list of dicts for cells DataFrame

    metabin_idx = 0

    for cluster in sorted(clusters):
        cluster_mask = cells_df[cluster_col] == cluster
        cluster_cell_indices = np.where(cluster_mask)[0]
        n_cluster_bins = len(cluster_cell_indices)

        if n_cluster_bins == 0:
            continue

        # Extract coordinates
        cluster_cells = cells_df.iloc[cluster_cell_indices]
        coords = cluster_cells[[x_col, y_col]].values

        # Assign bins to metabin groups via balanced bisection
        labels = _balanced_bisection(coords, n_bins)
        n_metabins = len(np.unique(labels))

        if verbose:
            print(f"  Cluster '{cluster}': {n_cluster_bins} bins -> {n_metabins} metabins")

        # Extract submatrices for this cluster's cells
        cov_sub = data.coverage[:, cluster_cell_indices]  # (n_positions, n_cluster_bins)
        if has_mutrate:
            mut_sub = data.mutrate[:, cluster_cell_indices]
        if has_mutcount:
            mc_sub = data.mutcount[:, cluster_cell_indices]

        # Get bin IDs
        if has_cell_id_col:
            bin_ids = cluster_cells[cell_id_col].values
        else:
            bin_ids = np.array(cluster_cells.index)

        # Aggregate each metabin group
        for group_label in range(n_metabins):
            group_mask = labels == group_label
            group_local_indices = np.where(group_mask)[0]

            if len(group_local_indices) == 0:
                continue

            # Aggregate coverage: sum across bins
            group_cov = cov_sub[:, group_local_indices].sum(axis=1)  # (n_positions, 1) dense matrix
            group_cov = np.asarray(group_cov).flatten()

            # Aggregate mutcount if available
            if has_mutcount:
                group_mc = mc_sub[:, group_local_indices].sum(axis=1)
                group_mc = np.asarray(group_mc).flatten()

            # Aggregate mutrate
            if has_mutrate:
                if has_mutcount:
                    # mutrate = sum(mutcount) / sum(coverage)
                    group_mr = np.zeros(n_positions, dtype=np.float32)
                    nonzero = group_cov > 0
                    group_mr[nonzero] = group_mc[nonzero] / group_cov[nonzero]
                else:
                    # Coverage-weighted average: sum(mutrate * coverage) / sum(coverage)
                    group_mut_cov = mut_sub[:, group_local_indices].multiply(
                        cov_sub[:, group_local_indices]
                    ).sum(axis=1)
                    group_mut_cov = np.asarray(group_mut_cov).flatten()
                    group_mr = np.zeros(n_positions, dtype=np.float32)
                    nonzero = group_cov > 0
                    group_mr[nonzero] = group_mut_cov[nonzero] / group_cov[nonzero]

            # Store as sparse columns, ensuring mutrate has the same sparsity as coverage
            cov_col = sparse.csc_matrix(group_cov.reshape(-1, 1), dtype=np.int32)
            metabin_cov_cols.append(cov_col)
            if has_mutrate:
                # Build mutrate column with same sparsity pattern as coverage
                # (preserves explicit zeros where coverage>0 but mutrate=0)
                cov_nz = cov_col.nonzero()[0]
                mr_col = sparse.csc_matrix(
                    (group_mr[cov_nz], (cov_nz, np.zeros(len(cov_nz), dtype=int))),
                    shape=(n_positions, 1), dtype=np.float32
                )
                metabin_mut_cols.append(mr_col)
            if has_mutcount:
                # Build mutcount column with same sparsity pattern as coverage
                cov_nz = cov_col.nonzero()[0] if not has_mutrate else cov_nz
                mc_col = sparse.csc_matrix(
                    (group_mc[cov_nz].astype(np.int32), (cov_nz, np.zeros(len(cov_nz), dtype=int))),
                    shape=(n_positions, 1), dtype=np.int32
                )
                metabin_mc_cols.append(mc_col)

            # Find medoid: bin closest to group centroid
            group_coords = coords[group_local_indices]
            centroid = group_coords.mean(axis=0)
            distances = np.sqrt(((group_coords - centroid) ** 2).sum(axis=1))
            medoid_local = group_local_indices[np.argmin(distances)]

            # Build metadata for this metabin
            metabin_meta.append({
                cluster_col: cluster,
                x_col: coords[medoid_local, 0],
                y_col: coords[medoid_local, 1],
                'bin100_ids': list(bin_ids[group_local_indices]),
                'n_bins': len(group_local_indices),
            })

            metabin_idx += 1

    if verbose:
        print(f"\nTotal metabins created: {metabin_idx}")

    # Build output matrices
    metabin_cov = sparse.hstack(metabin_cov_cols, format='csr')
    metabin_mut = sparse.hstack(metabin_mut_cols, format='csr') if has_mutrate else None
    metabin_mc = sparse.hstack(metabin_mc_cols, format='csr') if has_mutcount else None

    # Build cells DataFrame
    metabin_cells = pd.DataFrame(metabin_meta)
    metabin_cells.index = [f"metabin_{i}" for i in range(len(metabin_cells))]
    metabin_cells.index.name = 'cell'

    # Copy genepos
    genepos_copy = data.genepos.copy() if isinstance(data.genepos, pd.DataFrame) else list(data.genepos)

    if verbose:
        print(f"Output shape: ({n_positions}, {metabin_idx})")
        bins_per_metabin = metabin_cells['n_bins']
        print(f"Bins per metabin: min={bins_per_metabin.min()}, "
              f"median={bins_per_metabin.median():.0f}, "
              f"max={bins_per_metabin.max()}")

    return ShapeData(
        metabin_cov,
        metabin_mut,
        genepos_copy,
        metabin_cells,
        mutcount_sparse=metabin_mc
    )


def create_metabins_from_mapping(
    target_data: 'ShapeData',
    ref_metabin_data: 'ShapeData',
    mapping_df: pd.DataFrame,
    ref_bin_col: str = 'f2a3_bin',
    target_bin_col: str = 'mapped_dmso_bin',
    cell_id_col: str = 'bin100_cellID',
    ref_metabin_col: str = 'ref_metabin',
    verbose: bool = True
) -> 'ShapeData':
    """
    Create metabins for a target sample using groupings from a reference metabin ShapeData.

    Each metabin in the reference has a list of bin IDs (stored in 'bin100_ids').
    For each reference metabin, the corresponding target bins are looked up via
    mapping_df and aggregated into a target metabin.

    Parameters:
    -----------
    target_data : ShapeData
        The target ShapeData (e.g., DMSO). Must have cells as DataFrame.
    ref_metabin_data : ShapeData
        The reference metabin ShapeData (e.g., f2a3 metabin output from create_metabins).
        Must have 'bin100_ids' column in its cells DataFrame.
    mapping_df : pd.DataFrame
        DataFrame mapping reference bin IDs to target bin IDs.
    ref_bin_col : str
        Column in mapping_df for reference bin IDs (default: 'f2a3_bin').
    target_bin_col : str
        Column in mapping_df for target bin IDs (default: 'mapped_dmso_bin').
    cell_id_col : str
        Column name for bin IDs in target_data cells (default: 'bin100_cellID').
    ref_metabin_col : str
        Column name in the output cells DataFrame for storing the reference metabin ID
        (default: 'ref_metabin').
    verbose : bool
        Print progress information (default: True).

    Returns:
    --------
    ShapeData
        New ShapeData with metabins matching the reference grouping structure.
    """
    ref_cells = ref_metabin_data.cells
    if 'bin100_ids' not in ref_cells.columns:
        raise ValueError("ref_metabin_data must have 'bin100_ids' column (output of create_metabins).")

    # Build mapping: ref_bin -> target_bin
    ref_to_target = dict(zip(
        mapping_df[ref_bin_col].apply(_normalize_id),
        mapping_df[target_bin_col].apply(_normalize_id)
    ))

    if verbose:
        print(f"Creating target metabins from reference grouping")
        print(f"  Reference metabins: {len(ref_cells)}")
        print(f"  Mapping entries: {len(mapping_df)}")

    # Build metabin_dict and extra_meta from reference groupings
    metabin_dict = {}
    extra_meta = {}

    for metabin_name, row in ref_cells.iterrows():
        ref_bin_ids = [_normalize_id(b) for b in row['bin100_ids']]

        # Map ref bins -> target bins
        target_bin_ids = []
        for ref_bid in ref_bin_ids:
            target_bid = ref_to_target.get(ref_bid)
            if target_bid is not None:
                target_bin_ids.append(target_bid)

        metabin_dict[metabin_name] = target_bin_ids

        # Copy metadata from reference metabin (except bin100_ids and n_bins which are rebuilt)
        meta = {}
        for col in ref_cells.columns:
            if col not in ('bin100_ids', 'n_bins'):
                meta[col] = row[col]
        meta[ref_metabin_col] = metabin_name
        extra_meta[metabin_name] = meta

    return create_metabins_from_dict(
        target_data, metabin_dict, cell_id_col, extra_meta, verbose
    )


def create_metacells(
    target_data: 'ShapeData',
    group_col: str,
    cell_id_col: str = 'bin100_cellID',
    verbose: bool = True,
) -> 'ShapeData':
    """
    Create metacells by grouping cells that share the same value in `group_col`.

    Each unique non-NaN value in target_data.cells[group_col] becomes one metacell.
    Coverage / mutcount / mutrate aggregation is delegated to create_metabins_from_dict.

    Parameters
    ----------
    target_data : ShapeData
        The ShapeData whose cells will be aggregated. Must have cells as DataFrame.
    group_col : str
        Column name in target_data.cells used for groupby.
    cell_id_col : str
        Column name for cell IDs in target_data.cells (default: 'bin100_cellID').
        If absent, the cells index is used.
    verbose : bool
        Print progress information (default: True).

    Returns
    -------
    ShapeData
        New ShapeData with one metacell per unique value of `group_col`.
        Output cells DataFrame contains:
        - `group_col`: the group value
        - 'cellbarcodes': list of original cell IDs in this metacell
        - 'n_cells': number of cells in this metacell

        Cells with NaN in `group_col` are dropped (count logged when verbose).
    """
    if not target_data.cells_is_df:
        raise ValueError("target_data cells must be a DataFrame. Use to_cells_df() first.")

    cells_df = target_data.cells
    if group_col not in cells_df.columns:
        raise ValueError(f"Column '{group_col}' not found in cells metadata. "
                         f"Available columns: {list(cells_df.columns)}")

    # Resolve cell IDs (column or index)
    if cell_id_col in cells_df.columns:
        bin_ids = cells_df[cell_id_col].values
    else:
        bin_ids = np.array(cells_df.index)

    # Drop NaN in group_col
    n_total = len(cells_df)
    working_df = cells_df.dropna(subset=[group_col])
    n_dropped = n_total - len(working_df)
    if verbose and n_dropped > 0:
        print(f"Dropped {n_dropped} cells with NaN in '{group_col}'")

    if len(working_df) == 0:
        raise ValueError(f"All cells have NaN in '{group_col}'; nothing to group.")

    # Build metabin_dict and extra_meta from groupby
    metabin_dict = {}
    extra_meta = {}
    for group_value, sub_df in working_df.groupby(group_col, sort=True, dropna=True):
        positional_idx = cells_df.index.get_indexer(sub_df.index)
        key = str(group_value)
        metabin_dict[key] = list(bin_ids[positional_idx])
        extra_meta[key] = {group_col: group_value}

    if verbose:
        print(f"Grouping {len(working_df)} cells into {len(metabin_dict)} metacells by '{group_col}'")

    return create_metabins_from_dict(
        target_data,
        metabin_dict,
        cell_id_col=cell_id_col,
        extra_meta=extra_meta,
        verbose=verbose,
        ids_col_name='cellbarcodes',
        count_col_name='n_cells',
        index_prefix='metacell',
    )
