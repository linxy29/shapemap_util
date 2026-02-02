"""
I/O functions for ShapeData.

This module provides functions for creating ShapeData objects from VCF files.
"""

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .core import ShapeData


def from_vcf_pair(
    rRNA_path: str,
    steptwo_path: str,
    progress_interval: int = 100000
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

    Returns:
    --------
    ShapeData
    """
    from .core import ShapeData
    from cellsnp_util import read_paired_vcf_to_sparse

    cov, mut, pos, cells = read_paired_vcf_to_sparse(
        rRNA_path, steptwo_path, progress_interval
    )
    return ShapeData(cov, mut, pos, cells)


def from_cellsnp_vcf(
    vcf_path: str,
    progress_interval: int = 100000
) -> 'ShapeData':
    """
    Create ShapeData from cellSNP VCF file.

    Parameters:
    -----------
    vcf_path : str
        Path to cellSNP.cells.vcf.bgz.gz file
    progress_interval : int
        Print progress every N variants

    Returns:
    --------
    ShapeData
    """
    from .core import ShapeData
    from cellsnp_util import read_cellsnp_vcf_to_matrices_pysam_sparse

    cov, mut, pos, cells = read_cellsnp_vcf_to_matrices_pysam_sparse(
        vcf_path, progress_interval
    )
    return ShapeData(cov, mut, pos, cells)
