"""
Metabin functions for ShapeData.

This module provides functions for grouping spatial bins into metabins
based on cluster membership and spatial proximity.
"""

import numpy as np
import pandas as pd
from scipy import sparse
from typing import Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from .core import ShapeData


def create_metabins(
    data: 'ShapeData',
    n_bins: int,
    cluster_col: str = 'leiden',
    x_col: str = 'x',
    y_col: str = 'y',
    cell_id_col: str = 'bin100_cellID',
    n_neighbors: int = 10,
    verbose: bool = True
) -> 'ShapeData':
    """
    Group spatial bins from the same cluster into metabins using agglomerative clustering.

    For each cluster, bins are grouped into spatially compact metabins of approximately
    `n_bins` bins each. Coverage is summed, mutation rate is recomputed as a
    coverage-weighted average (or sum(mutcount)/sum(coverage) if mutcount is available).

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
    n_neighbors : int
        Number of neighbors for KNN connectivity graph (default: 10).
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
    from sklearn.cluster import AgglomerativeClustering
    from sklearn.neighbors import kneighbors_graph

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

        # Determine number of metabins for this cluster
        n_metabins = max(1, round(n_cluster_bins / n_bins))

        if verbose:
            print(f"  Cluster '{cluster}': {n_cluster_bins} bins -> {n_metabins} metabins")

        # Assign bins to metabin groups
        if n_metabins >= n_cluster_bins:
            # Each bin is its own metabin
            labels = np.arange(n_cluster_bins)
            n_metabins = n_cluster_bins
        elif n_metabins == 1:
            labels = np.zeros(n_cluster_bins, dtype=int)
        else:
            # Build KNN connectivity graph
            k = min(n_neighbors, n_cluster_bins - 1)
            connectivity = kneighbors_graph(coords, n_neighbors=k, mode='connectivity',
                                            include_self=False)
            agg = AgglomerativeClustering(
                n_clusters=n_metabins,
                connectivity=connectivity,
                linkage='ward'
            )
            labels = agg.fit_predict(coords)

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

            # Store as sparse columns
            metabin_cov_cols.append(sparse.csc_matrix(group_cov.reshape(-1, 1), dtype=np.int32))
            if has_mutrate:
                metabin_mut_cols.append(sparse.csc_matrix(group_mr.reshape(-1, 1), dtype=np.float32))
            if has_mutcount:
                metabin_mc_cols.append(sparse.csc_matrix(group_mc.reshape(-1, 1), dtype=np.int32))

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
