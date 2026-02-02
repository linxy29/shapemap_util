"""
ShapeData: A container class for single-cell SHAPE-MaP data.

This package provides the ShapeData class that encapsulates sparse coverage
and mutation rate matrices along with position and cell metadata.

Usage
-----
>>> from shape_data import ShapeData
>>>
>>> # Create from VCF files
>>> data = ShapeData.from_vcf_pair(rRNA_path, steptwo_path)
>>>
>>> # Or create directly from matrices
>>> data = ShapeData(coverage, mutrate, genepos, cells)
>>>
>>> # Filter and analyze
>>> data = data.filter_thresholds(coverage_threshold=1500, mutrate_threshold=0.2)
>>> data.add_genepos_counts()
>>> data.join_cell_metadata(clustering_df, 'leiden', on='cell')

Submodules
----------
- core: Main ShapeData class
- filter: Filtering functions (filter_thresholds, filter_cells, filter_genepos, subset_by_cells, subset_by_genepos)
- analysis: Analysis functions (calculate_reactivity, calculate_cell_correlation, get_cell_stats, get_position_stats, differential_reactivity_lm, differential_reactivity_wilcoxon)
- metadata: Metadata functions (add_genepos_counts, join_cell_metadata, join_genepos_metadata, add_control_data)
- plot: Plotting functions (plot_violin, plot_violin_multi)
- io: I/O functions (from_vcf_pair, from_cellsnp_vcf)
"""

from .core import ShapeData

# Import commonly used functions for convenience
from .filter import filter_thresholds, filter_cells, filter_genepos, subset_by_cells, subset_by_genepos
from .analysis import (
    calculate_reactivity, calculate_cell_correlation, get_cell_stats, get_position_stats,
    differential_reactivity_lm, differential_reactivity_wilcoxon
)
from .metadata import add_genepos_counts, join_cell_metadata, join_genepos_metadata, add_control_data
from .plot import plot_violin, plot_violin_multi
from .io import from_vcf_pair, from_cellsnp_vcf

__all__ = [
    # Main class
    'ShapeData',
    # Filter functions
    'filter_thresholds',
    'filter_cells',
    'filter_genepos',
    'subset_by_cells',
    'subset_by_genepos',
    # Analysis functions
    'calculate_reactivity',
    'calculate_cell_correlation',
    'get_cell_stats',
    'get_position_stats',
    'differential_reactivity_lm',
    'differential_reactivity_wilcoxon',
    # Metadata functions
    'add_genepos_counts',
    'join_cell_metadata',
    'join_genepos_metadata',
    'add_control_data',
    # Plot functions
    'plot_violin',
    'plot_violin_multi',
    # I/O functions
    'from_vcf_pair',
    'from_cellsnp_vcf',
]
