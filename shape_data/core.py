"""
Core ShapeData class.

This module provides the main ShapeData class for single-cell SHAPE-MaP data.
"""

import numpy as np
import pandas as pd
from scipy import sparse
from typing import Union, Optional, List, Dict, Any, Tuple, Callable

# Import helper functions from submodules
from .filter import (
    filter_thresholds as _filter_thresholds,
    filter_cells as _filter_cells,
    filter_genepos as _filter_genepos,
    subset_by_cells as _subset_by_cells,
    subset_by_genepos as _subset_by_genepos
)
from .analysis import (
    calculate_reactivity as _calculate_reactivity,
    calculate_cell_correlation as _calculate_cell_correlation,
    calculate_auroc as _calculate_auroc,
    get_cell_stats as _get_cell_stats,
    get_position_stats as _get_position_stats,
    differential_reactivity_lm as _differential_reactivity_lm,
    differential_reactivity_wilcoxon as _differential_reactivity_wilcoxon
)
from .metadata import (
    add_genepos_counts as _add_genepos_counts,
    join_cell_metadata as _join_cell_metadata,
    join_genepos_metadata as _join_genepos_metadata,
    add_control_data as _add_control_data
)
from .plot import (
    plot_violin as _plot_violin,
    plot_violin_multi as _plot_violin_multi
)
from .io import (
    from_vcf_pair as _from_vcf_pair,
    from_cellsnp_vcf as _from_cellsnp_vcf
)


class ShapeData:
    """
    Container for sparse single-cell SHAPE-MaP data.

    This class holds coverage and mutation rate matrices in sparse format
    along with position and cell metadata. It provides methods for filtering,
    subsetting, and joining metadata.

    Attributes:
    -----------
    coverage : scipy.sparse.csr_matrix
        Coverage matrix (positions x cells)
    mutrate : scipy.sparse.csr_matrix or None
        Mutation rate matrix (positions x cells)
    reactivity : scipy.sparse.csr_matrix or None
        Reactivity matrix (positions x cells), computed via calculate_reactivity()
    genepos : pd.DataFrame or list
        Position identifiers. If DataFrame, should have columns like 'gene', 'pos', 'ref'.
    cells : pd.DataFrame or list
        Cell information. If DataFrame, cell barcodes should be the index.

    Examples:
    ---------
    >>> # Create from VCF reading functions
    >>> cov, mut, pos, cells = read_paired_vcf_to_sparse(rRNA_path, steptwo_path)
    >>> data = ShapeData(cov, mut, pos, cells)

    >>> # Filter by coverage and mutation rate
    >>> data = data.filter_thresholds(coverage_threshold=1500, mutrate_threshold=0.2)

    >>> # Add cell metadata
    >>> data.add_genepos_counts()
    >>> data.join_cell_metadata(clustering_df, 'leiden', on='cell')

    >>> # Filter by cell metadata
    >>> data = data.filter_cells({'Total': (100000, None)})

    >>> # Filter by cell indices
    >>> data = data.filter_cells(cell_indices=['barcode1', 'barcode2'])
    """

    def __init__(
        self,
        coverage_sparse: sparse.csr_matrix,
        mutrate_sparse: sparse.csr_matrix,
        genepos: Union[pd.DataFrame, List, pd.Series, pd.Index],
        cells: Union[pd.DataFrame, List, np.ndarray],
        reactivity_sparse: Optional[sparse.csr_matrix] = None,
        normalized_reactivity: Optional[np.ndarray] = None
    ):
        """
        Initialize ShapeData.

        Parameters:
        -----------
        coverage_sparse : scipy.sparse.csr_matrix
            Coverage matrix (positions x cells)
        mutrate_sparse : scipy.sparse.csr_matrix
            Mutation rate matrix (positions x cells)
        genepos : pd.DataFrame, list, pd.Series, or pd.Index
            Position identifiers
        cells : pd.DataFrame, list, or array-like
            Cell information. Can be:
            - pd.DataFrame: Cell metadata with cell names as index
            - list/array: Simple list of cell barcodes
        reactivity_sparse : scipy.sparse.csr_matrix or None, optional
            Reactivity matrix (positions x cells), can be None (default: None)
        normalized_reactivity : np.ndarray or None, optional
            Normalized reactivity matrix (positions x cells), can be None (default: None)
        """
        # Ensure sparse matrices are CSR format (handle sparse, dense arrays, and None)
        self.coverage = self._to_csr(coverage_sparse)
        self.mutrate = self._to_csr(mutrate_sparse)
        # Keep reactivity as dense array (not sparse) since it contains NaN values
        # that don't benefit from sparse storage and analysis functions expect dense
        self.reactivity = self._keep_dense(reactivity_sparse)
        self.normalized_reactivity = self._keep_dense(normalized_reactivity)

        # Store genepos
        if isinstance(genepos, (pd.Series, pd.Index)):
            self.genepos = list(genepos)
        else:
            self.genepos = genepos

        # Store cells
        if isinstance(cells, pd.DataFrame):
            self.cells = cells.copy()
            # Ensure the index is named
            if self.cells.index.name is None and 'cell' in self.cells.columns:
                self.cells = self.cells.set_index('cell')
        elif isinstance(cells, list):
            self.cells = cells.copy() if isinstance(cells, list) else list(cells)
        else:
            self.cells = list(cells)

        # Validate dimensions
        self._validate()

    @staticmethod
    def _to_csr(matrix):
        """Convert matrix to CSR format, handling None, sparse, and dense arrays."""
        if matrix is None:
            return None
        if isinstance(matrix, sparse.csr_matrix):
            return matrix
        if sparse.issparse(matrix):
            return matrix.tocsr()
        # Handle numpy arrays
        return sparse.csr_matrix(matrix)

    @staticmethod
    def _keep_dense(matrix):
        """Keep matrix as dense numpy array, converting sparse to dense if needed.

        Used for reactivity matrices which contain NaN values and should remain
        dense for proper handling by analysis functions like np.isnan().
        """
        if matrix is None:
            return None
        if sparse.issparse(matrix):
            return matrix.toarray()
        # Already a dense numpy array
        return np.asarray(matrix)

    def _validate(self):
        """Validate that all dimensions are consistent."""
        n_positions, n_cells = self.coverage.shape

        if self.mutrate is not None and self.mutrate.shape != self.coverage.shape:
            raise ValueError(
                f"Coverage shape {self.coverage.shape} != mutrate shape {self.mutrate.shape}"
            )

        if self.reactivity is not None and self.reactivity.shape != self.coverage.shape:
            raise ValueError(
                f"Coverage shape {self.coverage.shape} != reactivity shape {self.reactivity.shape}"
            )

        if self.normalized_reactivity is not None and self.normalized_reactivity.shape != self.coverage.shape:
            raise ValueError(
                f"Coverage shape {self.coverage.shape} != normalized_reactivity shape {self.normalized_reactivity.shape}"
            )

        # Check genepos length
        genepos_len = len(self.genepos) if isinstance(self.genepos, list) else len(self.genepos)
        if genepos_len != n_positions:
            raise ValueError(
                f"Number of positions in genepos ({genepos_len}) != matrix rows ({n_positions})"
            )

        # Check cells length
        cells_len = len(self.cells)
        if cells_len != n_cells:
            raise ValueError(
                f"Number of cells ({cells_len}) != matrix columns ({n_cells})"
            )

    @property
    def shape(self) -> Tuple[int, int]:
        """Return shape as (n_positions, n_cells)."""
        return self.coverage.shape

    @property
    def n_positions(self) -> int:
        """Number of positions/windows."""
        return self.coverage.shape[0]

    @property
    def n_cells(self) -> int:
        """Number of cells."""
        return self.coverage.shape[1]

    @property
    def cells_is_df(self) -> bool:
        """Check if cells is a DataFrame."""
        return isinstance(self.cells, pd.DataFrame)

    @property
    def cell_names(self) -> List[str]:
        """Get list of cell names/barcodes."""
        if self.cells_is_df:
            return self.cells.index.tolist()
        else:
            return list(self.cells)

    def __repr__(self) -> str:
        cells_type = "DataFrame" if self.cells_is_df else "list"
        genepos_type = "DataFrame" if isinstance(self.genepos, pd.DataFrame) else "list"
        nnz = self.coverage.nnz
        sparsity = 100 * (1 - nnz / (self.n_positions * self.n_cells)) if self.n_positions * self.n_cells > 0 else 0
        has_mutrate = self.mutrate is not None
        has_reactivity = self.reactivity is not None
        has_normalized = self.normalized_reactivity is not None
        return (
            f"ShapeData(\n"
            f"  shape=({self.n_positions:,} positions x {self.n_cells:,} cells),\n"
            f"  non-zeros={nnz:,}, sparsity={sparsity:.2f}%,\n"
            f"  mutrate={has_mutrate}, reactivity={has_reactivity}, normalized_reactivity={has_normalized},\n"
            f"  genepos={genepos_type}, cells={cells_type}\n"
            f")"
        )

    def copy(self) -> 'ShapeData':
        """Create a deep copy of the data."""
        genepos_copy = self.genepos.copy() if isinstance(self.genepos, pd.DataFrame) else list(self.genepos)
        cells_copy = self.cells.copy() if self.cells_is_df else list(self.cells)
        reactivity_copy = self.reactivity.copy() if self.reactivity is not None else None
        normalized_copy = self.normalized_reactivity.copy() if self.normalized_reactivity is not None else None
        return ShapeData(
            self.coverage.copy(),
            self.mutrate.copy() if self.mutrate is not None else None,
            genepos_copy,
            cells_copy,
            reactivity_copy,
            normalized_copy
        )

    def to_cells_df(self) -> 'ShapeData':
        """
        Convert cells to DataFrame if it's currently a list.
        Returns a new ShapeData with cells as DataFrame.
        """
        if self.cells_is_df:
            return self

        cells_df = pd.DataFrame(index=self.cells)
        cells_df.index.name = 'cell'

        return ShapeData(
            self.coverage.copy(),
            self.mutrate.copy() if self.mutrate is not None else None,
            self.genepos.copy() if isinstance(self.genepos, pd.DataFrame) else list(self.genepos),
            cells_df,
            self.reactivity.copy() if self.reactivity is not None else None,
            self.normalized_reactivity.copy() if self.normalized_reactivity is not None else None
        )

    def to_genepos_df(self) -> 'ShapeData':
        """
        Convert genepos to DataFrame if it's currently a list.
        Returns a new ShapeData with genepos as DataFrame.

        If genepos is a list of 'gene.pos' format strings (e.g., '18S.100'),
        it will be parsed into separate 'gene' and 'pos' columns with
        'gene.pos' as the index. Otherwise, a simple 'position' column is created.
        """
        if isinstance(self.genepos, pd.DataFrame):
            return self

        # Check if genepos contains 'gene.pos' format strings
        sample_pos = str(self.genepos[0]) if len(self.genepos) > 0 else ''
        if '.' in sample_pos:
            # Parse 'gene.pos' format (e.g., '18S.100' -> gene='18S', pos=100)
            genes = []
            positions = []
            for pos_str in self.genepos:
                parts = str(pos_str).rsplit('.', 1)  # rsplit to handle genes with dots
                if len(parts) == 2:
                    genes.append(parts[0])
                    try:
                        positions.append(int(parts[1]))
                    except ValueError:
                        positions.append(parts[1])
                else:
                    genes.append(pos_str)
                    positions.append(None)

            genepos_df = pd.DataFrame({
                'gene': genes,
                'pos': positions
            }, index=self.genepos)
            genepos_df.index.name = 'gene.pos'
        else:
            # Simple position list without gene.pos format
            genepos_df = pd.DataFrame({'position': self.genepos}, index=self.genepos)
            genepos_df.index.name = 'position'

        return ShapeData(
            self.coverage.copy(),
            self.mutrate.copy() if self.mutrate is not None else None,
            genepos_df,
            self.cells.copy() if self.cells_is_df else list(self.cells),
            self.reactivity.copy() if self.reactivity is not None else None,
            self.normalized_reactivity.copy() if self.normalized_reactivity is not None else None
        )

    # =========================================================================
    # Metadata methods (delegating to metadata module)
    # =========================================================================

    def add_genepos_counts(self, inplace: bool = True) -> Optional['ShapeData']:
        """
        Add 'n_genepos' column to cells metadata, counting non-zero entries per cell.

        Parameters:
        -----------
        inplace : bool
            If True, modify in place. If False, return a new ShapeData.

        Returns:
        --------
        ShapeData or None
            Returns new object if inplace=False, otherwise None.
        """
        return _add_genepos_counts(self, inplace=inplace)

    def join_cell_metadata(
        self,
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
        return _join_cell_metadata(self, metadata_df, columns, on, how, inplace)

    def join_genepos_metadata(
        self,
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
        """
        return _join_genepos_metadata(self, metadata_df, columns, on, how, inplace)

    def add_control_data(
        self,
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

        Parameters:
        -----------
        control_df : pd.DataFrame
            DataFrame with control data
        cluster_col : str
            Column name containing cluster identifiers (default: 'sample')
        cluster_prefix : str
            Prefix to remove from cluster_col values (default: 'dmso_')
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
        """
        return _add_control_data(
            self, control_df, cluster_col, cluster_prefix,
            coverage_col, mutrate_col, on, control_suffix, inplace, verbose
        )

    # =========================================================================
    # Filtering methods (delegating to filter module)
    # =========================================================================

    def filter_thresholds(
        self,
        coverage_threshold: float,
        mutrate_threshold: float,
        inplace: bool = False
    ) -> 'ShapeData':
        """
        Filter data by coverage and mutation rate thresholds.

        Parameters:
        -----------
        coverage_threshold : float
            Minimum coverage threshold
        mutrate_threshold : float
            Maximum mutation rate threshold
        inplace : bool
            If True, modify in place (not recommended for sparse operations)

        Returns:
        --------
        ShapeData
            Filtered data
        """
        result = _filter_thresholds(self, coverage_threshold, mutrate_threshold)
        if inplace:
            self.coverage = result.coverage
            self.mutrate = result.mutrate
            self.genepos = result.genepos
            self.cells = result.cells
            self.reactivity = result.reactivity
            return self
        return result

    def filter_cells(
        self,
        cell_filters: Optional[Dict[str, Any]] = None,
        cell_indices: Optional[Union[List, np.ndarray, Callable]] = None,
        inplace: bool = False
    ) -> 'ShapeData':
        """
        Filter data based on cell metadata or indices.

        Parameters:
        -----------
        cell_filters : dict, optional
            Dictionary of column filters (only works when cells is a DataFrame)
        cell_indices : list, array, or callable, optional
            Direct index-based filtering
        inplace : bool
            If True, modify in place (not recommended)

        Returns:
        --------
        ShapeData
            Filtered data
        """
        result = _filter_cells(self, cell_filters, cell_indices)
        if inplace:
            self.coverage = result.coverage
            self.mutrate = result.mutrate
            self.genepos = result.genepos
            self.cells = result.cells
            self.reactivity = result.reactivity
            return self
        return result

    def subset_by_cells(self, cell_list: List[str]) -> 'ShapeData':
        """
        Subset data to include only cells in the given list.

        Parameters:
        -----------
        cell_list : list of str
            Cell barcodes to keep

        Returns:
        --------
        ShapeData
            Subsetted data
        """
        return _subset_by_cells(self, cell_list)

    def filter_genepos(
        self,
        genepos_filters: Optional[Dict[str, Any]] = None,
        genepos_indices: Optional[Union[List, np.ndarray, Callable]] = None,
        inplace: bool = False
    ) -> 'ShapeData':
        """
        Filter data based on genepos metadata or indices.

        Parameters:
        -----------
        genepos_filters : dict, optional
            Dictionary of column filters (only works when genepos is a DataFrame).
            Keys are column names, values can be:
            - single value: exact match (e.g., {'gene': '18S'})
            - list: membership (e.g., {'gene': ['18S', '28S']})
            - tuple of (min, max): range filter (e.g., {'pos': (100, 500)})
              Use None for open-ended ranges: (100, None) for >= 100
        genepos_indices : list, array, or callable, optional
            Direct index-based filtering
        inplace : bool
            If True, modify in place (not recommended)

        Returns:
        --------
        ShapeData
            Filtered data
        """
        result = _filter_genepos(self, genepos_filters, genepos_indices)
        if inplace:
            self.coverage = result.coverage
            self.mutrate = result.mutrate
            self.genepos = result.genepos
            self.cells = result.cells
            self.reactivity = result.reactivity
            return self
        return result

    def subset_by_genepos(self, genepos_indices: List[int]) -> 'ShapeData':
        """
        Subset data to include only positions in the given index list.

        Parameters:
        -----------
        genepos_indices : list of int
            Position indices to keep

        Returns:
        --------
        ShapeData
            Subsetted data
        """
        return _subset_by_genepos(self, genepos_indices)

    # =========================================================================
    # Analysis methods (delegating to analysis module)
    # =========================================================================

    def calculate_reactivity(
        self,
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

        Also computes normalized reactivity per cell per gene (default: box normalization).

        Parameters:
        -----------
        cluster_column : str
            Column name in cells DataFrame containing cluster assignments
        control_prefix : str
            Prefix for control mutrate columns in genepos
        store_as : str
            Store the reactivity matrix as an attribute with this name
        normalize : bool
            If True, also compute normalized reactivity per cell per gene (default: True)
        normalize_method : str
            Normalization method: 'box' (default), 'zscore', or 'minmax'
            - 'box': divide by mean of 90th-95th percentile, clip negatives to 0
            - 'zscore': (x - mean) / std
            - 'minmax': scale to [0, 1]
        store_normalized_as : str
            Store the normalized reactivity matrix as this attribute (default: 'normalized_reactivity')
        drop_mutrate : bool
            If True, set self.mutrate to None after calculation
        filter_all_na : bool
            If True, filter out gene.pos rows that have all NA values in the
            reactivity matrix. This also filters corresponding rows in genepos,
            coverage, and mutrate matrices. (default: False)
        verbose : bool
            Whether to print progress information

        Returns:
        --------
        None
            The reactivity matrix is stored as self.{store_as}
            The normalized reactivity is stored as self.{store_normalized_as}
        """
        return _calculate_reactivity(
            self, cluster_column, control_prefix, store_as,
            normalize, normalize_method, store_normalized_as,
            drop_mutrate, filter_all_na, verbose
        )

    def calculate_cell_correlation(
        self,
        gene: str,
        cluster_column: str = 'leiden',
        method: str = 'pearson',
        min_shared_positions: int = 100,
        save_mean_as: Optional[str] = None,
        verbose: bool = True
    ) -> pd.DataFrame:
        """
        Calculate pairwise correlation between cells within the same cluster.

        Parameters:
        -----------
        gene : str
            Gene name to calculate correlations for
        cluster_column : str
            Column name in cells DataFrame containing cluster assignments
        method : str
            Correlation method: 'pearson', 'spearman', or 'kendall'
        min_shared_positions : int
            Minimum number of shared non-NaN positions
        save_mean_as : str or None
            Column name to save mean correlation for each cell
        verbose : bool
            Whether to print progress information

        Returns:
        --------
        pd.DataFrame
            Correlation dataframe
        """
        return _calculate_cell_correlation(
            self, gene, cluster_column, method, min_shared_positions, save_mean_as, verbose
        )

    def calculate_auroc(
        self,
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

        Example:
        --------
        >>> # Calculate AUROC for 18S rRNA
        >>> dot_bracket_18S = "(((...)))..."  # or path to file
        >>> data.calculate_auroc('18S', dot_bracket_18S)
        >>> print(data.cells[['leiden', 'AUROC_18S']].head())
        """
        _calculate_auroc(
            self, gene, dot_bracket, skip_positions, min_positions, save_as, verbose
        )

    def get_cell_stats(self, verbose: bool = True) -> None:
        """
        Compute statistics per cell and add to cells metadata.

        Adds the following columns to self.cells:
            - n_genepos: number of non-zero entries
            - mean_coverage: mean coverage across non-zero entries
            - mean_mutrate: mean mutation rate across non-zero entries

        Parameters:
        -----------
        verbose : bool
            Whether to print progress information (default: True)

        Returns:
        --------
        None
            Statistics are added directly to self.cells
        """
        _get_cell_stats(self, verbose)

    def get_position_stats(self, verbose: bool = True) -> None:
        """
        Compute statistics per position and add to genepos metadata.

        Adds the following columns to self.genepos:
            - n_cells: number of cells with data
            - mean_coverage: mean coverage across cells
            - mean_mutrate: mean mutation rate across cells

        Parameters:
        -----------
        verbose : bool
            Whether to print progress information (default: True)

        Returns:
        --------
        None
            Statistics are added directly to self.genepos
        """
        _get_position_stats(self, verbose)

    def differential_reactivity_lm(
        self,
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
        Perform linear model analysis for differential reactivity.

        Fits: reactivity ~ b_0 + b_1*X1 + b_2*X2 + ...
        for each position. By default tests all positions.

        Parameters:
        -----------
        predictors : str or list of str
            Column name(s) in cells metadata to use as predictors.
            Examples: 'diet', ['diet', 'region'], ['treatment', 'batch']
        gene : str, list of str, or None
            Gene name(s) to analyze. If None (default), tests all positions.
        use_normalized : bool
            If True, use normalized_reactivity instead of raw reactivity (default: False)
        cluster_column : str, optional
            Column for cluster filtering
        cluster_value : str, optional
            Specific cluster value to analyze
        min_cells : int
            Minimum cells required per position (default: 10)
        correct_multiple : bool
            Apply multiple testing correction (default: True)
        correction_method : str
            Correction method (default: 'fdr_bh')
        verbose : bool
            Print progress (default: True)

        Returns:
        --------
        pd.DataFrame
            Results with coefficient, std_err, t_statistic, pvalue, pvalue_adj

        Example:
        --------
        >>> # Test all positions (raw reactivity)
        >>> results = data.differential_reactivity_lm(['diet', 'region'])
        >>> # Use normalized reactivity
        >>> results = data.differential_reactivity_lm('diet', use_normalized=True)
        >>> # Test specific gene
        >>> results = data.differential_reactivity_lm('diet', gene='18S')
        >>> sig = results[results['pvalue_adj'] < 0.05]
        """
        return _differential_reactivity_lm(
            self, predictors, gene, use_normalized, cluster_column, cluster_value,
            min_cells, correct_multiple, correction_method, verbose
        )

    def differential_reactivity_wilcoxon(
        self,
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
        Perform Wilcoxon rank-sum test for differential reactivity.

        Tests whether reactivity differs between two groups at each position.
        By default tests all positions.

        Parameters:
        -----------
        group_column : str
            Column defining groups (e.g., 'diet', 'region')
        group1 : str or list of str
            Value(s) for first group
        group2 : str or list of str
            Value(s) for second group
        gene : str, list of str, or None
            Gene name(s) to analyze. If None (default), tests all positions.
        use_normalized : bool
            If True, use normalized_reactivity instead of raw reactivity (default: False)
        min_cells_per_group : int
            Minimum cells per group (default: 5)
        correct_multiple : bool
            Apply multiple testing correction (default: True)
        correction_method : str
            Correction method (default: 'fdr_bh')
        alternative : str
            'two-sided', 'less', or 'greater' (default: 'two-sided')
        verbose : bool
            Print progress (default: True)

        Returns:
        --------
        pd.DataFrame
            Results with statistic, pvalue, effect_size, mean_group1/2, log2_fold_change

        Example:
        --------
        >>> # Test all positions (raw reactivity)
        >>> results = data.differential_reactivity_wilcoxon(
        ...     group_column='diet', group1='HFD', group2='CTRL'
        ... )
        >>> # Use normalized reactivity
        >>> results = data.differential_reactivity_wilcoxon(
        ...     group_column='diet', group1='HFD', group2='CTRL', use_normalized=True
        ... )
        >>> # Test specific gene
        >>> results = data.differential_reactivity_wilcoxon(
        ...     group_column='diet', group1='HFD', group2='CTRL', gene='18S'
        ... )
        >>> sig = results[results['pvalue_adj'] < 0.05]
        """
        return _differential_reactivity_wilcoxon(
            self, group_column, group1, group2, gene, use_normalized,
            min_cells_per_group, correct_multiple, correction_method,
            alternative, verbose
        )

    # =========================================================================
    # Plotting methods (delegating to plot module)
    # =========================================================================

    def plot_violin(
        self,
        column: Union[str, List[str]],
        source: str = 'auto',
        groupby: Optional[str] = None,
        title: Optional[str] = None,
        xlabel: Optional[str] = None,
        ylabel: Optional[str] = None,
        figsize: Optional[Tuple[int, int]] = None,
        dot_size: Union[float, str] = 'auto',
        dot_alpha: float = 0.5,
        dot_color: str = 'black',
        violin_alpha: float = 0.7,
        palette: Optional[str] = 'Set2',
        jitter: float = 0.15,
        show_means: bool = True,
        ax: Optional[Any] = None,
        ncols: int = 2,
        **kwargs
    ) -> Any:
        """
        Plot violin plot with dot plot overlay for metadata column(s).

        Parameters:
        -----------
        column : str or list of str
            Column name(s) to plot
        source : str
            Data source: 'genepos', 'cells', or 'auto'
        groupby : str, optional
            Column name to group by
        title : str, optional
            Plot title
        xlabel : str, optional
            X-axis label
        ylabel : str, optional
            Y-axis label
        figsize : tuple, optional
            Figure size
        ncols : int
            Number of columns in subplot grid
        dot_size : float or 'auto'
            Size of dots. If 'auto', calculates based on number of data points
            (smaller dots for more points). Default: 'auto'
        dot_alpha : float
            Transparency of dots
        dot_color : str
            Color of dots
        violin_alpha : float
            Transparency of violin plot
        palette : str, optional
            Color palette for groups
        jitter : float
            Amount of jitter for dots
        show_means : bool
            Whether to show mean markers
        ax : matplotlib.axes.Axes, optional
            Existing axes to plot on
        **kwargs
            Additional arguments passed to seaborn.violinplot

        Returns:
        --------
        matplotlib.axes.Axes
            The axes object containing the plot
        """
        return _plot_violin(
            self, column, source, groupby, title, xlabel, ylabel,
            figsize, dot_size, dot_alpha, dot_color, violin_alpha,
            palette, jitter, show_means, ax, ncols, **kwargs
        )

    def plot_violin_multi(
        self,
        columns: List[str],
        source: str = 'auto',
        groupby: Optional[str] = None,
        ncols: int = 2,
        title: Optional[str] = None,
        figsize: Optional[Tuple[int, int]] = None,
        dot_size: Union[float, str] = 'auto',
        dot_alpha: float = 0.5,
        dot_color: str = 'black',
        violin_alpha: float = 0.7,
        palette: Optional[str] = 'Set2',
        jitter: float = 0.15,
        show_means: bool = True,
        **kwargs
    ) -> Any:
        """
        Plot violin plots for multiple columns in subplots.

        Parameters:
        -----------
        columns : list of str
            Column names to plot
        source : str
            Data source: 'genepos', 'cells', or 'auto'
        groupby : str, optional
            Column name to group by
        ncols : int
            Number of columns in subplot grid
        title : str, optional
            Overall figure title
        figsize : tuple, optional
            Figure size
        dot_size : float or 'auto'
            Size of dots. If 'auto', calculates based on number of data points
            (smaller dots for more points). Default: 'auto'
        dot_alpha : float
            Transparency of dots
        dot_color : str
            Color of dots
        violin_alpha : float
            Transparency of violin
        palette : str, optional
            Color palette for groups
        jitter : float
            Amount of jitter for dots
        show_means : bool
            Whether to show mean markers
        **kwargs
            Additional arguments passed to seaborn.violinplot

        Returns:
        --------
        numpy.ndarray of matplotlib.axes.Axes
            Array of axes objects
        """
        return _plot_violin_multi(
            self, columns, source, groupby, ncols, title, figsize,
            dot_size, dot_alpha, dot_color, violin_alpha, palette,
            jitter, show_means, **kwargs
        )

    # =========================================================================
    # Class methods (delegating to io module)
    # =========================================================================

    @classmethod
    def from_vcf_pair(
        cls,
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
        return _from_vcf_pair(rRNA_path, steptwo_path, progress_interval)

    @classmethod
    def from_cellsnp_vcf(
        cls,
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
        return _from_cellsnp_vcf(vcf_path, progress_interval)

    # =========================================================================
    # Save/Load methods
    # =========================================================================

    def save(self, path: str, compress: bool = True) -> None:
        """
        Save ShapeData to a directory.

        Creates a directory containing:
        - coverage.npz: Sparse coverage matrix
        - mutrate.npz: Sparse mutation rate matrix (if exists)
        - reactivity.npy: Dense reactivity matrix (if exists)
        - normalized_reactivity.npy: Dense normalized reactivity matrix (if exists)
        - genepos.parquet or genepos.pkl: Position metadata
        - cells.parquet or cells.pkl: Cell metadata

        Parameters:
        -----------
        path : str
            Directory path to save to (will be created if doesn't exist)
        compress : bool
            Whether to use compression (default: True)

        Example:
        --------
        >>> data.save('./my_shapedata')
        >>> # Later...
        >>> data = ShapeData.load('./my_shapedata')
        """
        import os

        os.makedirs(path, exist_ok=True)

        # Save sparse matrices
        sparse.save_npz(os.path.join(path, 'coverage.npz'), self.coverage, compressed=compress)

        if self.mutrate is not None:
            sparse.save_npz(os.path.join(path, 'mutrate.npz'), self.mutrate, compressed=compress)

        # Save reactivity (dense numpy array)
        if self.reactivity is not None:
            if compress:
                np.savez_compressed(os.path.join(path, 'reactivity.npz'), reactivity=self.reactivity)
            else:
                np.save(os.path.join(path, 'reactivity.npy'), self.reactivity)

        # Save normalized_reactivity (dense numpy array)
        if self.normalized_reactivity is not None:
            if compress:
                np.savez_compressed(os.path.join(path, 'normalized_reactivity.npz'),
                                   normalized_reactivity=self.normalized_reactivity)
            else:
                np.save(os.path.join(path, 'normalized_reactivity.npy'), self.normalized_reactivity)

        # Save genepos
        if isinstance(self.genepos, pd.DataFrame):
            self.genepos.to_parquet(os.path.join(path, 'genepos.parquet'))
        else:
            pd.to_pickle(self.genepos, os.path.join(path, 'genepos.pkl'))

        # Save cells
        if self.cells_is_df:
            self.cells.to_parquet(os.path.join(path, 'cells.parquet'))
        else:
            pd.to_pickle(self.cells, os.path.join(path, 'cells.pkl'))

        print(f"Saved ShapeData to {path}")

    @classmethod
    def load(cls, path: str) -> 'ShapeData':
        """
        Load ShapeData from a directory.

        Parameters:
        -----------
        path : str
            Directory path to load from

        Returns:
        --------
        ShapeData

        Example:
        --------
        >>> data = ShapeData.load('./my_shapedata')
        """
        import os

        # Load coverage (required)
        coverage = sparse.load_npz(os.path.join(path, 'coverage.npz'))

        # Load mutrate (optional)
        mutrate_path = os.path.join(path, 'mutrate.npz')
        mutrate = sparse.load_npz(mutrate_path) if os.path.exists(mutrate_path) else None

        # Load reactivity (optional, check both formats)
        reactivity = None
        reactivity_npz = os.path.join(path, 'reactivity.npz')
        reactivity_npy = os.path.join(path, 'reactivity.npy')
        if os.path.exists(reactivity_npz):
            reactivity = np.load(reactivity_npz)['reactivity']
        elif os.path.exists(reactivity_npy):
            reactivity = np.load(reactivity_npy)

        # Load normalized_reactivity (optional, check both formats)
        normalized_reactivity = None
        norm_react_npz = os.path.join(path, 'normalized_reactivity.npz')
        norm_react_npy = os.path.join(path, 'normalized_reactivity.npy')
        if os.path.exists(norm_react_npz):
            normalized_reactivity = np.load(norm_react_npz)['normalized_reactivity']
        elif os.path.exists(norm_react_npy):
            normalized_reactivity = np.load(norm_react_npy)

        # Load genepos
        genepos_parquet = os.path.join(path, 'genepos.parquet')
        genepos_pkl = os.path.join(path, 'genepos.pkl')
        if os.path.exists(genepos_parquet):
            genepos = pd.read_parquet(genepos_parquet)
        elif os.path.exists(genepos_pkl):
            genepos = pd.read_pickle(genepos_pkl)
        else:
            raise FileNotFoundError(f"genepos file not found in {path}")

        # Load cells
        cells_parquet = os.path.join(path, 'cells.parquet')
        cells_pkl = os.path.join(path, 'cells.pkl')
        if os.path.exists(cells_parquet):
            cells = pd.read_parquet(cells_parquet)
        elif os.path.exists(cells_pkl):
            cells = pd.read_pickle(cells_pkl)
        else:
            raise FileNotFoundError(f"cells file not found in {path}")

        print(f"Loaded ShapeData from {path}")

        return cls(coverage, mutrate, genepos, cells, reactivity, normalized_reactivity)

    def to_pickle(self, path: str) -> None:
        """
        Save ShapeData to a single pickle file.

        Parameters:
        -----------
        path : str
            File path (e.g., 'data.pkl' or 'data.shapedata')

        Example:
        --------
        >>> data.to_pickle('my_data.pkl')
        >>> data = ShapeData.from_pickle('my_data.pkl')
        """
        import pickle

        with open(path, 'wb') as f:
            pickle.dump({
                'coverage': self.coverage,
                'mutrate': self.mutrate,
                'genepos': self.genepos,
                'cells': self.cells,
                'reactivity': self.reactivity,
                'normalized_reactivity': self.normalized_reactivity,
            }, f, protocol=pickle.HIGHEST_PROTOCOL)

        print(f"Saved ShapeData to {path}")

    @classmethod
    def from_pickle(cls, path: str) -> 'ShapeData':
        """
        Load ShapeData from a pickle file.

        Parameters:
        -----------
        path : str
            File path

        Returns:
        --------
        ShapeData

        Example:
        --------
        >>> data = ShapeData.from_pickle('my_data.pkl')
        """
        import pickle

        with open(path, 'rb') as f:
            d = pickle.load(f)

        print(f"Loaded ShapeData from {path}")

        return cls(
            d['coverage'],
            d['mutrate'],
            d['genepos'],
            d['cells'],
            d.get('reactivity'),
            d.get('normalized_reactivity')
        )
