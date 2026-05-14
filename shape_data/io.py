"""
I/O functions for ShapeData.

This module provides functions for creating ShapeData objects from VCF files,
and combining multiple ShapeData objects.
"""

import numpy as np
import pandas as pd
from scipy import sparse
from typing import TYPE_CHECKING, Dict

if TYPE_CHECKING:
    from .core import ShapeData


def combine(
    data_dict: Dict[str, 'ShapeData'],
    sample_column: str = 'sample',
    cell_prefix: bool = True,
    join: str = 'inner',
    verbose: bool = True
) -> 'ShapeData':
    """
    Combine multiple ShapeData objects into one by horizontally stacking matrices.

    Cells (columns) are concatenated, and a sample column is added to cells
    metadata. Positions (rows) are aligned across samples according to ``join``.

    Parameters:
    -----------
    data_dict : dict
        Dictionary of {sample_name: ShapeData}.
    sample_column : str
        Column name added to cells DataFrame to track sample origin (default: 'sample').
    cell_prefix : bool
        If True, prefix cell barcodes with sample name to avoid collisions (default: True).
    join : str
        How to align positions across samples:
        - 'inner' (default): keep only positions present in ALL samples (intersection).
        - 'outer': keep positions present in ANY sample (union). For samples that
          are missing a given position, that row is filled with sparse zeros in
          coverage, mutrate, and mutcount.
    verbose : bool
        Print progress information (default: True).

    Returns:
    --------
    ShapeData

    Example:
    --------
    >>> sd1 = ShapeData(cov1, mut1, pos1, cells1)
    >>> sd2 = ShapeData(cov2, mut2, pos2, cells2)
    >>> combined = ShapeData.combine({'nain3_1': sd1, 'dmso_1': sd2})
    >>> # Keep all positions across samples instead of intersecting:
    >>> combined = ShapeData.combine({'nain3_1': sd1, 'dmso_1': sd2}, join='outer')
    """
    from .core import ShapeData

    if not data_dict:
        raise ValueError("data_dict is empty")
    if join not in ('inner', 'outer'):
        raise ValueError(f"join must be 'inner' or 'outer', got {join!r}")

    names = list(data_dict.keys())
    first = data_dict[names[0]]
    ref_pos = first.genepos

    # Build position index for each sample
    def _get_idx(sd):
        if isinstance(sd.genepos, pd.DataFrame):
            return sd.genepos.index
        return pd.Index(sd.genepos)

    # Combine position indices across ALL samples (intersection or union)
    combined_idx = _get_idx(first)
    all_match = True
    for sd in data_dict.values():
        cur_idx = _get_idx(sd)
        if not cur_idx.equals(combined_idx):
            all_match = False
            if join == 'inner':
                combined_idx = combined_idx.intersection(cur_idx)
            else:  # outer
                combined_idx = combined_idx.union(cur_idx, sort=False)

    if not all_match and verbose:
        kind = 'common' if join == 'inner' else 'union'
        print(f"WARNING: Positions differ across samples. "
              f"Using {len(combined_idx)} {kind} positions.")

    cov_list = []
    mut_list = []
    mc_list = []
    all_cells = []
    sample_labels = []
    has_mutcount = all(hasattr(sd, 'mutcount') and sd.mutcount is not None for sd in data_dict.values())

    n_target = len(combined_idx)

    def _reindex_to_outer(mat, cur_idx):
        """Reindex sparse matrix rows from cur_idx to combined_idx ordering.
        Rows present in combined_idx but missing from cur_idx are sparse zeros.
        """
        if mat is None:
            return None
        # positions[i] = row index in cur_idx for combined_idx[i], or -1 if absent
        positions = cur_idx.get_indexer(combined_idx)
        valid_target = np.where(positions >= 0)[0]
        valid_sample = positions[valid_target]
        sub = mat[valid_sample]  # rows from this sample that map into the union
        coo = sub.tocoo()
        new_rows = valid_target[coo.row]
        return sparse.coo_matrix(
            (coo.data, (new_rows, coo.col)),
            shape=(n_target, mat.shape[1]),
            dtype=mat.dtype,
        ).tocsr()

    for name, sd in data_dict.items():
        cur_idx = _get_idx(sd)

        if all_match:
            cov_list.append(sd.coverage)
            mut_list.append(sd.mutrate)
            if has_mutcount:
                mc_list.append(sd.mutcount)
        elif join == 'inner':
            mask = np.isin(cur_idx, combined_idx)
            cov_list.append(sd.coverage[mask])
            mut_list.append(sd.mutrate[mask] if sd.mutrate is not None else None)
            if has_mutcount:
                mc_list.append(sd.mutcount[mask])
        else:  # outer
            cov_list.append(_reindex_to_outer(sd.coverage, cur_idx))
            mut_list.append(_reindex_to_outer(sd.mutrate, cur_idx))
            if has_mutcount:
                mc_list.append(_reindex_to_outer(sd.mutcount, cur_idx))

        # Cell barcodes
        cell_names = sd.cell_names
        if cell_prefix:
            cell_names = [f"{name}_{c}" for c in cell_names]
        all_cells.extend(cell_names)
        sample_labels.extend([name] * sd.n_cells)

    # Stack matrices
    combined_cov = sparse.hstack(cov_list, format='csr')
    if all(m is not None for m in mut_list):
        combined_mut = sparse.hstack(mut_list, format='csr')
    else:
        combined_mut = None
    combined_mc = sparse.hstack(mc_list, format='csr') if has_mutcount else None

    # Build genepos for combined positions
    if all_match:
        genepos_copy = ref_pos.copy() if isinstance(ref_pos, pd.DataFrame) else list(ref_pos)
    elif isinstance(ref_pos, pd.DataFrame):
        if join == 'inner':
            genepos_copy = ref_pos.loc[ref_pos.index.isin(combined_idx)].copy()
        else:  # outer: gather genepos rows from all samples, keep first occurrence
            all_genepos = pd.concat(
                [sd.genepos for sd in data_dict.values()
                 if isinstance(sd.genepos, pd.DataFrame)]
            )
            all_genepos = all_genepos[~all_genepos.index.duplicated(keep='first')]
            genepos_copy = all_genepos.loc[combined_idx].copy()
    else:
        if join == 'inner':
            genepos_copy = [p for p in ref_pos if p in combined_idx]
        else:
            genepos_copy = list(combined_idx)

    result = ShapeData(combined_cov, combined_mut, genepos_copy, all_cells, mutcount_sparse=combined_mc)
    result = result.to_cells_df()
    result.cells[sample_column] = sample_labels

    if verbose:
        print(f"Combined {len(names)} samples ({join} join): "
              f"{result.n_positions:,} positions x {result.n_cells:,} cells")
        for name in names:
            n = data_dict[name].n_cells
            print(f"  {name}: {n:,} cells")

    return result


def from_vcf_pair(
    rRNA_path: str,
    steptwo_path: str,
    progress_interval: int = 100000,
    include_mutant_counts: bool = False
) -> 'ShapeData':
    """
    Create ShapeData from paired rRNA and steptwo VCF files.

    Parameters:
    -----------
    rRNA_path : str
        Path to rRNA VCF file
    steptwo_path : str
        Path to steptwo VCF file
    progress_interval : int
        Print progress every N variants
    include_mutant_counts : bool
        If True, also store mutant read count matrix (default: False)

    Returns:
    --------
    ShapeData
    """
    from .core import ShapeData
    from cellsnp_util import read_paired_vcf_to_sparse

    result = read_paired_vcf_to_sparse(
        rRNA_path, steptwo_path, progress_interval,
        include_mutant_counts=include_mutant_counts
    )
    if include_mutant_counts:
        cov, mut, ad, pos, cells = result
        return ShapeData(cov, mut, pos, cells, mutcount_sparse=ad)
    else:
        cov, mut, pos, cells = result
        return ShapeData(cov, mut, pos, cells)


def from_cellsnp_vcf(
    vcf_path: str,
    progress_interval: int = 100000,
    target_genes: list = None,
    include_mutant_counts: bool = False
) -> 'ShapeData':
    """
    Create ShapeData from cellSNP VCF file.

    Parameters:
    -----------
    vcf_path : str
        Path to cellSNP.cells.vcf.bgz.gz file
    progress_interval : int
        Print progress every N variants
    target_genes : list of str, optional
        If provided, only load positions belonging to these genes.
        Uses indexed random access (vcf.fetch) when available.
    include_mutant_counts : bool
        If True, also store mutant read count matrix (default: False)

    Returns:
    --------
    ShapeData
    """
    from .core import ShapeData
    from cellsnp_util import read_cellsnp_vcf_to_matrices_pysam_sparse

    print(f"Loading cellSNP VCF: {vcf_path}")
    if target_genes is not None:
        print(f"  target_genes: {target_genes}")
    print(f"  include_mutant_counts: {include_mutant_counts}")

    result = read_cellsnp_vcf_to_matrices_pysam_sparse(
        vcf_path, progress_interval, target_genes=target_genes,
        include_mutant_counts=include_mutant_counts
    )
    if include_mutant_counts:
        cov, mut, ad, pos, cells = result
        print(f"Creating ShapeData: {cov.shape[0]:,} positions x {cov.shape[1]:,} cells (with mutant counts)")
        return ShapeData(cov, mut, pos, cells, mutcount_sparse=ad)
    else:
        cov, mut, pos, cells = result
        print(f"Creating ShapeData: {cov.shape[0]:,} positions x {cov.shape[1]:,} cells")
        return ShapeData(cov, mut, pos, cells)
