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
    verbose: bool = True
) -> 'ShapeData':
    """
    Combine multiple ShapeData objects into one by horizontally stacking matrices.

    All objects must share the same positions (rows). Cells (columns) are
    concatenated, and a sample column is added to cells metadata.

    Parameters:
    -----------
    data_dict : dict
        Dictionary of {sample_name: ShapeData}.
    sample_column : str
        Column name added to cells DataFrame to track sample origin (default: 'sample').
    cell_prefix : bool
        If True, prefix cell barcodes with sample name to avoid collisions (default: True).
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
    """
    from .core import ShapeData

    if not data_dict:
        raise ValueError("data_dict is empty")

    names = list(data_dict.keys())
    first = data_dict[names[0]]
    ref_pos = first.genepos

    # Build position index for each sample
    def _get_idx(sd):
        if isinstance(sd.genepos, pd.DataFrame):
            return sd.genepos.index
        return pd.Index(sd.genepos)

    # Find common positions across ALL samples
    common_idx = _get_idx(first)
    all_match = True
    for sd in data_dict.values():
        cur_idx = _get_idx(sd)
        if not cur_idx.equals(common_idx):
            all_match = False
            common_idx = common_idx.intersection(cur_idx)

    if not all_match and verbose:
        print(f"WARNING: Positions differ across samples. "
              f"Using {len(common_idx)} common positions.")

    cov_list = []
    mut_list = []
    mc_list = []
    all_cells = []
    sample_labels = []
    has_mutcount = all(hasattr(sd, 'mutcount') and sd.mutcount is not None for sd in data_dict.values())

    for name, sd in data_dict.items():
        cur_idx = _get_idx(sd)

        if all_match:
            cov_list.append(sd.coverage)
            mut_list.append(sd.mutrate)
            if has_mutcount:
                mc_list.append(sd.mutcount)
        else:
            mask = np.isin(cur_idx, common_idx)
            cov_list.append(sd.coverage[mask])
            mut_list.append(sd.mutrate[mask] if sd.mutrate is not None else None)
            if has_mutcount:
                mc_list.append(sd.mutcount[mask])

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

    # Build genepos for common positions
    if all_match:
        genepos_copy = ref_pos.copy() if isinstance(ref_pos, pd.DataFrame) else list(ref_pos)
    elif isinstance(ref_pos, pd.DataFrame):
        genepos_copy = ref_pos.loc[ref_pos.index.isin(common_idx)].copy()
    else:
        genepos_copy = [p for p in ref_pos if p in common_idx]

    result = ShapeData(combined_cov, combined_mut, genepos_copy, all_cells, mutcount_sparse=combined_mc)
    result = result.to_cells_df()
    result.cells[sample_column] = sample_labels

    if verbose:
        print(f"Combined {len(names)} samples: "
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
