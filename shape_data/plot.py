"""
Plotting functions for ShapeData.

This module provides functions for creating visualizations from ShapeData objects.
"""

import numpy as np
import pandas as pd
from typing import Union, Optional, List, Tuple, Any, TYPE_CHECKING

if TYPE_CHECKING:
    from .core import ShapeData


def _calculate_dot_size(n: int) -> float:
    """
    Calculate appropriate dot size based on number of data points.

    Uses a logarithmic scale to determine dot size:
    - n < 50: size = 5
    - n = 100: size ~ 3
    - n = 500: size ~ 1.5
    - n = 1000: size ~ 1
    - n = 5000: size ~ 0.5
    - n > 10000: size = 0.3

    Parameters:
    -----------
    n : int
        Number of data points

    Returns:
    --------
    float
        Recommended dot size
    """
    if n < 50:
        return 5.0
    elif n > 10000:
        return 0.3
    else:
        # Logarithmic scaling: size decreases as n increases
        # Formula: size = a - b * log10(n)
        # Calibrated so that n=100 -> ~3, n=1000 -> ~1
        size = 7.0 - 2.0 * np.log10(n)
        return max(0.3, min(5.0, size))


def plot_violin(
    data: 'ShapeData',
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

    If a single column is provided, creates a single violin plot.
    If multiple columns are provided, creates subplots for each column.

    Parameters:
    -----------
    data : ShapeData
        The ShapeData object to plot from
    column : str or list of str
        Column name(s) to plot (from genepos or cells metadata).
        If a list, creates subplots for each column.
    source : str
        Data source: 'genepos', 'cells', or 'auto' (default: 'auto')
        'auto' will search in both and use whichever contains the column
    groupby : str, optional
        Column name to group by (e.g., 'leiden' for cluster comparison)
        If None, creates a single violin plot
    title : str, optional
        Plot title. If None, uses column name
    xlabel : str, optional
        X-axis label
    ylabel : str, optional
        Y-axis label. If None, uses column name
    figsize : tuple, optional
        Figure size (width, height) in inches. If None, auto-calculated.
    ncols : int
        Number of columns in subplot grid when plotting multiple columns (default: 2)
    dot_size : float or 'auto'
        Size of dots in the dot plot. If 'auto', calculates based on number of
        data points (smaller dots for more points). Default: 'auto'
    dot_alpha : float
        Transparency of dots (default: 0.5)
    dot_color : str
        Color of dots (default: 'black'). Ignored if groupby is specified.
    violin_alpha : float
        Transparency of violin plot (default: 0.7)
    palette : str, optional
        Color palette for groups (default: 'Set2')
    jitter : float
        Amount of jitter for dots (default: 0.15)
    show_means : bool
        Whether to show mean markers on violin (default: True)
    ax : matplotlib.axes.Axes, optional
        Existing axes to plot on. If None, creates new figure.
    **kwargs
        Additional arguments passed to seaborn.violinplot

    Returns:
    --------
    matplotlib.axes.Axes
        The axes object containing the plot

    Example:
    --------
    >>> # Plot cell correlation distribution
    >>> plot_violin(data, 'mean_corr_18S', groupby='leiden')
    >>>
    >>> # Plot proportion of cells per position
    >>> plot_violin(data, 'prop_cells', source='genepos', groupby='gene')
    >>>
    >>> # Plot multiple columns in subplots
    >>> plot_violin(data, ['n_cells', 'prop_cells'], source='genepos')
    """
    import matplotlib.pyplot as plt
    try:
        import seaborn as sns
    except ImportError:
        raise ImportError("seaborn is required for plot_violin. Install with: pip install seaborn")

    # If multiple columns provided, use plot_violin_multi
    if isinstance(column, list):
        return plot_violin_multi(
            data,
            columns=column,
            source=source,
            groupby=groupby,
            ncols=ncols,
            title=title,
            figsize=figsize,
            dot_size=dot_size,
            dot_alpha=dot_alpha,
            dot_color=dot_color,
            violin_alpha=violin_alpha,
            palette=palette,
            jitter=jitter,
            show_means=show_means,
            **kwargs
        )

    # Single column - use default figsize if not provided
    if figsize is None:
        figsize = (8, 6)

    # Determine data source
    if source == 'auto':
        in_genepos = isinstance(data.genepos, pd.DataFrame) and column in data.genepos.columns
        in_cells = data.cells_is_df and column in data.cells.columns

        if in_genepos and in_cells:
            print(f"Warning: '{column}' found in both genepos and cells. Using genepos. "
                  f"Specify source='cells' to use cells instead.")
            source = 'genepos'
        elif in_genepos:
            source = 'genepos'
        elif in_cells:
            source = 'cells'
        else:
            available_genepos = list(data.genepos.columns) if isinstance(data.genepos, pd.DataFrame) else []
            available_cells = list(data.cells.columns) if data.cells_is_df else []
            raise ValueError(
                f"Column '{column}' not found.\n"
                f"Available in genepos: {available_genepos}\n"
                f"Available in cells: {available_cells}"
            )

    # Get the data
    if source == 'genepos':
        if not isinstance(data.genepos, pd.DataFrame):
            raise ValueError("genepos must be a DataFrame to plot from it")
        if column not in data.genepos.columns:
            raise ValueError(f"Column '{column}' not found in genepos. Available: {list(data.genepos.columns)}")
        data_df = data.genepos.copy()
    elif source == 'cells':
        if not data.cells_is_df:
            raise ValueError("cells must be a DataFrame to plot from it")
        if column not in data.cells.columns:
            raise ValueError(f"Column '{column}' not found in cells. Available: {list(data.cells.columns)}")
        data_df = data.cells.copy()
    else:
        raise ValueError(f"source must be 'genepos', 'cells', or 'auto', got '{source}'")

    # Validate groupby column
    if groupby is not None and groupby not in data_df.columns:
        raise ValueError(f"groupby column '{groupby}' not found in {source}. Available: {list(data_df.columns)}")

    # Drop NaN values for the plot column
    plot_data = data_df[[column] + ([groupby] if groupby else [])].dropna(subset=[column])

    if len(plot_data) == 0:
        raise ValueError(f"No valid (non-NaN) data in column '{column}'")

    # Calculate dot_size if 'auto'
    if dot_size == 'auto':
        dot_size = _calculate_dot_size(len(plot_data))

    # Create figure if no axes provided
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Create violin plot
    if groupby:
        # Grouped violin plot
        sns.violinplot(
            data=plot_data,
            x=groupby,
            y=column,
            ax=ax,
            alpha=violin_alpha,
            palette=palette,
            inner=None,  # Don't show inner box plot
            **kwargs
        )

        # Add strip plot (dots) on top
        sns.stripplot(
            data=plot_data,
            x=groupby,
            y=column,
            ax=ax,
            size=dot_size,
            alpha=dot_alpha,
            color=dot_color,
            jitter=jitter
        )

        # Add mean markers if requested
        if show_means:
            means = plot_data.groupby(groupby)[column].mean()
            for i, (group, mean_val) in enumerate(means.items()):
                ax.scatter([i], [mean_val], color='red', s=50, zorder=5, marker='D',
                          edgecolors='white', linewidths=1, label='Mean' if i == 0 else '')
    else:
        # Single violin plot
        sns.violinplot(
            data=plot_data,
            y=column,
            ax=ax,
            alpha=violin_alpha,
            color=palette if isinstance(palette, str) and palette.startswith('#') else 'steelblue',
            inner=None,
            **kwargs
        )

        # Add strip plot (dots) on top
        sns.stripplot(
            data=plot_data,
            y=column,
            ax=ax,
            size=dot_size,
            alpha=dot_alpha,
            color=dot_color,
            jitter=jitter
        )

        # Add mean marker if requested
        if show_means:
            mean_val = plot_data[column].mean()
            ax.scatter([0], [mean_val], color='red', s=50, zorder=5, marker='D',
                      edgecolors='white', linewidths=1, label='Mean')

    # Set labels and title
    if title is None:
        title = f"Distribution of {column}" + (f" by {groupby}" if groupby else "")
    ax.set_title(title)

    if ylabel is None:
        ylabel = column
    ax.set_ylabel(ylabel)

    if xlabel is not None:
        ax.set_xlabel(xlabel)
    elif groupby:
        ax.set_xlabel(groupby)

    # Add legend for mean marker
    if show_means:
        ax.legend(loc='upper right')

    plt.tight_layout()

    return ax


def plot_violin_multi(
    data: 'ShapeData',
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
    data : ShapeData
        The ShapeData object to plot from
    columns : list of str
        Column names to plot
    source : str
        Data source: 'genepos', 'cells', or 'auto' (default: 'auto')
    groupby : str, optional
        Column name to group by within each subplot
    ncols : int
        Number of columns in subplot grid (default: 2)
    title : str, optional
        Overall figure title
    figsize : tuple, optional
        Figure size. If None, auto-calculated
    dot_size : float or 'auto'
        Size of dots. If 'auto', calculates based on number of data points
        (smaller dots for more points). Default: 'auto'
    dot_alpha : float
        Transparency of dots (default: 0.5)
    dot_color : str
        Color of dots (default: 'black')
    violin_alpha : float
        Transparency of violin (default: 0.7)
    palette : str, optional
        Color palette for groups (default: 'Set2')
    jitter : float
        Amount of jitter for dots (default: 0.15)
    show_means : bool
        Whether to show mean markers (default: True)
    **kwargs
        Additional arguments passed to seaborn.violinplot

    Returns:
    --------
    numpy.ndarray of matplotlib.axes.Axes
        Array of axes objects

    Example:
    --------
    >>> # Plot two columns in subplots
    >>> plot_violin_multi(data, ['n_cells', 'prop_cells'], source='genepos')
    >>>
    >>> # Plot with grouping
    >>> plot_violin_multi(data, ['mean_corr_18S', 'mean_corr_28S'], groupby='leiden')
    """
    import matplotlib.pyplot as plt
    try:
        import seaborn as sns
    except ImportError:
        raise ImportError("seaborn is required. Install with: pip install seaborn")

    if len(columns) < 1:
        raise ValueError("Must provide at least one column")

    # Determine data source
    if source == 'auto':
        in_genepos = isinstance(data.genepos, pd.DataFrame) and columns[0] in data.genepos.columns
        in_cells = data.cells_is_df and columns[0] in data.cells.columns

        if in_genepos:
            source = 'genepos'
        elif in_cells:
            source = 'cells'
        else:
            available_genepos = list(data.genepos.columns) if isinstance(data.genepos, pd.DataFrame) else []
            available_cells = list(data.cells.columns) if data.cells_is_df else []
            raise ValueError(
                f"Column '{columns[0]}' not found.\n"
                f"Available in genepos: {available_genepos}\n"
                f"Available in cells: {available_cells}"
            )

    # Get the data
    if source == 'genepos':
        if not isinstance(data.genepos, pd.DataFrame):
            raise ValueError("genepos must be a DataFrame")
        data_df = data.genepos.copy()
    elif source == 'cells':
        if not data.cells_is_df:
            raise ValueError("cells must be a DataFrame")
        data_df = data.cells.copy()
    else:
        raise ValueError(f"source must be 'genepos', 'cells', or 'auto'")

    # Validate columns exist
    missing_cols = [c for c in columns if c not in data_df.columns]
    if missing_cols:
        raise ValueError(f"Columns not found in {source}: {missing_cols}")

    # Validate groupby
    if groupby is not None and groupby not in data_df.columns:
        raise ValueError(f"groupby column '{groupby}' not found in {source}")

    # Calculate grid dimensions
    n_plots = len(columns)
    nrows = (n_plots + ncols - 1) // ncols

    if figsize is None:
        figsize = (5 * ncols, 4 * nrows)

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
    axes = axes.flatten()

    for i, col in enumerate(columns):
        ax = axes[i]
        plot_data = data_df[[col] + ([groupby] if groupby else [])].dropna(subset=[col])

        # Calculate dot_size for this subplot if 'auto'
        current_dot_size = _calculate_dot_size(len(plot_data)) if dot_size == 'auto' else dot_size

        if groupby:
            sns.violinplot(
                data=plot_data,
                x=groupby,
                y=col,
                ax=ax,
                alpha=violin_alpha,
                palette=palette,
                inner=None,
                **kwargs
            )
            sns.stripplot(
                data=plot_data,
                x=groupby,
                y=col,
                ax=ax,
                size=current_dot_size,
                alpha=dot_alpha,
                color=dot_color,
                jitter=jitter
            )
            if show_means:
                means = plot_data.groupby(groupby)[col].mean()
                for j, (group, mean_val) in enumerate(means.items()):
                    ax.scatter([j], [mean_val], color='red', s=50, zorder=5, marker='D',
                              edgecolors='white', linewidths=1)
            ax.set_xlabel(groupby)
        else:
            sns.violinplot(
                data=plot_data,
                y=col,
                ax=ax,
                alpha=violin_alpha,
                color='steelblue',
                inner=None,
                **kwargs
            )
            sns.stripplot(
                data=plot_data,
                y=col,
                ax=ax,
                size=current_dot_size,
                alpha=dot_alpha,
                color=dot_color,
                jitter=jitter
            )
            if show_means:
                mean_val = plot_data[col].mean()
                ax.scatter([0], [mean_val], color='red', s=50, zorder=5, marker='D',
                          edgecolors='white', linewidths=1)

        ax.set_title(col)
        ax.set_ylabel(col)

    # Hide unused subplots
    for i in range(n_plots, len(axes)):
        axes[i].set_visible(False)

    if title:
        fig.suptitle(title, y=1.02)

    plt.tight_layout()
    return axes[:n_plots]


def plot_reactivity(
    data: 'ShapeData',
    pos_range: Optional[Tuple[int, int]] = None,
    gene: Optional[str] = None,
    cluster_col: str = 'leiden',
    clusters: Optional[List] = None,
    title: Optional[str] = None,
    figsize: Optional[Tuple[int, int]] = None,
    palette: Optional[str] = 'Set2',
    line_alpha: float = 0.8,
    fill_alpha: float = 0.2,
    variance_style: str = 'fill',
    spread: str = 'se',
    use_normalized: bool = True,
    ax: Optional[Any] = None,
    smoothing: int = 1,
    log_coverage: bool = False,
) -> Any:
    """
    Plot coverage and reactivity profiles across positions, grouped by cluster.

    Creates two vertically stacked subplots sharing the same x-axis:
    - Top: mean coverage per position per cluster (with variance)
    - Bottom: mean reactivity per position per cluster (with variance)

    Parameters:
    -----------
    data : ShapeData
        The ShapeData object to plot from. Must have coverage matrix and
        reactivity matrix (call calculate_reactivity first if needed).
        cells must be a DataFrame with cluster_col column.
    pos_range : tuple of (int, int), optional
        (start, end) row indices into genepos to select positions (inclusive).
        If None, all positions are used.
    gene : str, optional
        Gene name to filter positions. Only used when genepos is a DataFrame
        with a 'gene' column. Overrides pos_range if both given.
    cluster_col : str
        Column name in cells metadata for cluster labels (default: 'leiden').
    clusters : list, optional
        Specific cluster values to plot. If None, plots all clusters.
    title : str, optional
        Overall figure title.
    figsize : tuple, optional
        Figure size (width, height). Default: (12, 6).
    palette : str, optional
        Matplotlib colormap or seaborn palette name (default: 'Set2').
    line_alpha : float
        Transparency of mean lines (default: 0.8).
    fill_alpha : float
        Transparency of variance shading (default: 0.2).
    variance_style : str
        How to display variance: 'fill' for shaded region (mean +/- std),
        'errorbar' for error bars (default: 'fill').
    ax : numpy.ndarray of matplotlib.axes.Axes, optional
        Existing axes array of length 2 to plot on. If None, creates new figure.
    smoothing : int
        Window size for sliding-window smoothing along positions (default: 1,
        no smoothing). When smoothing=3, each position's coverage and reactivity
        become the mean of itself and its two neighbors. Edge positions use
        reflect mode.
    log_coverage : bool
        If True, plot the coverage subplot on a logarithmic y-axis (default:
        False). Zero/non-positive coverage positions are already shown as gaps,
        so they are simply omitted from the log-scaled axis.

    Returns:
    --------
    numpy.ndarray of matplotlib.axes.Axes
        Array of two axes [coverage_ax, reactivity_ax].

    Example:
    --------
    >>> # Plot full gene profile by cluster
    >>> data.plot_profile(gene='18S', cluster_col='leiden')
    >>>
    >>> # Plot a position range with specific clusters
    >>> data.plot_profile(pos_range=(100, 300), clusters=[0, 1, 2])
    """
    import matplotlib.pyplot as plt
    from scipy import sparse

    if not data.cells_is_df:
        raise ValueError("cells must be a DataFrame with cluster information")
    if cluster_col not in data.cells.columns:
        raise ValueError(f"'{cluster_col}' not found in cells. Available: {list(data.cells.columns)}")
    if use_normalized:
        if data.normalized_reactivity is None:
            raise ValueError("normalized_reactivity is None. Call calculate_reactivity() with normalize=True first.")
        reactivity_matrix = data.normalized_reactivity
    else:
        if data.reactivity is None:
            raise ValueError("reactivity matrix is None. Call calculate_reactivity() first.")
        reactivity_matrix = data.reactivity

    # Determine position mask
    n_pos = data.coverage.shape[0]
    pos_mask = np.ones(n_pos, dtype=bool)

    if gene is not None and isinstance(data.genepos, pd.DataFrame) and 'gene' in data.genepos.columns:
        pos_mask &= data.genepos['gene'].values == gene

    if pos_range is not None and isinstance(data.genepos, pd.DataFrame) and 'pos' in data.genepos.columns:
        pos_vals = data.genepos['pos'].values
        start, end = pos_range
        pos_mask &= (pos_vals >= start) & (pos_vals <= end)
    elif pos_range is not None:
        start, end = pos_range
        idx_mask = np.zeros(n_pos, dtype=bool)
        idx_mask[start:min(end + 1, n_pos)] = True
        pos_mask &= idx_mask

    pos_indices = np.where(pos_mask)[0]

    if len(pos_indices) == 0:
        raise ValueError("No positions selected")

    # Build x-axis labels from genepos
    if isinstance(data.genepos, pd.DataFrame) and 'pos' in data.genepos.columns:
        x_positions = data.genepos['pos'].values[pos_indices]
    else:
        x_positions = pos_indices

    # Get cluster assignments
    cluster_labels = data.cells[cluster_col]
    if clusters is None:
        clusters = sorted(cluster_labels.unique())

    # Get color palette
    try:
        import seaborn as sns
        colors = sns.color_palette(palette, n_colors=len(clusters))
    except ImportError:
        cmap = plt.get_cmap(palette)
        colors = [cmap(i / max(len(clusters) - 1, 1)) for i in range(len(clusters))]

    # Slice matrices to selected positions
    cov_sub = data.coverage[pos_indices, :]
    react_sub = reactivity_matrix[pos_indices, :]

    # Apply smoothing along positions if requested
    if smoothing > 1:
        from scipy.ndimage import uniform_filter1d
        # Coverage: no NaN, uniform_filter1d is fine
        if sparse.issparse(cov_sub):
            cov_sub = sparse.csr_matrix(uniform_filter1d(cov_sub.toarray().astype(float), size=smoothing, axis=0))
        else:
            cov_sub = uniform_filter1d(np.asarray(cov_sub, dtype=float), size=smoothing, axis=0)
        # Reactivity: contains NaN, use nan-aware rolling mean with reflect padding
        react_float = react_sub.astype(float)
        half = smoothing // 2
        n_pos = react_float.shape[0]
        react_smoothed = np.empty_like(react_float)
        for i in range(n_pos):
            indices = []
            for j in range(i - half, i + half + 1):
                if j < 0:
                    indices.append(-j)  # reflect
                elif j >= n_pos:
                    indices.append(2 * (n_pos - 1) - j)  # reflect
                else:
                    indices.append(j)
            react_smoothed[i] = np.nanmean(react_float[indices], axis=0)
        react_sub = react_smoothed

    # Create figure
    if figsize is None:
        figsize = (12, 6)

    if ax is None:
        fig, axes = plt.subplots(2, 1, figsize=figsize, sharex=True)
    else:
        axes = ax
        fig = axes[0].get_figure()

    ax_cov, ax_react = axes[0], axes[1]

    for ci, cluster in enumerate(clusters):
        cell_mask = (cluster_labels == cluster).values
        n_cells = cell_mask.sum()
        if n_cells == 0:
            continue

        # Coverage: sparse matrix, extract columns for this cluster
        cov_cluster = cov_sub[:, cell_mask]
        if sparse.issparse(cov_cluster):
            cov_dense = cov_cluster.toarray()
        else:
            cov_dense = np.asarray(cov_cluster)
        cov_mean = np.mean(cov_dense, axis=1)

        # Reactivity: dense array with NaN
        react_cluster = react_sub[:, cell_mask]
        react_mean = np.nanmean(react_cluster, axis=1)

        # Compute spread bands
        if spread == 'std':
            cov_lo = cov_mean - np.std(cov_dense, axis=1)
            cov_hi = cov_mean + np.std(cov_dense, axis=1)
            react_lo = react_mean - np.nanstd(react_cluster, axis=1)
            react_hi = react_mean + np.nanstd(react_cluster, axis=1)
        elif spread == 'se':
            cov_se = np.std(cov_dense, axis=1) / np.sqrt(n_cells)
            cov_lo = cov_mean - cov_se
            cov_hi = cov_mean + cov_se
            n_valid = np.sum(~np.isnan(react_cluster), axis=1).astype(float)
            n_valid[n_valid == 0] = np.nan
            react_se = np.nanstd(react_cluster, axis=1) / np.sqrt(n_valid)
            react_lo = react_mean - react_se
            react_hi = react_mean + react_se
        elif spread == 'ci95':
            cov_se = np.std(cov_dense, axis=1) / np.sqrt(n_cells)
            cov_lo = cov_mean - 1.96 * cov_se
            cov_hi = cov_mean + 1.96 * cov_se
            n_valid = np.sum(~np.isnan(react_cluster), axis=1).astype(float)
            n_valid[n_valid == 0] = np.nan
            react_se = np.nanstd(react_cluster, axis=1) / np.sqrt(n_valid)
            react_lo = react_mean - 1.96 * react_se
            react_hi = react_mean + 1.96 * react_se
        elif spread == 'iqr':
            cov_lo = np.percentile(cov_dense, 25, axis=1)
            cov_hi = np.percentile(cov_dense, 75, axis=1)
            react_lo = np.nanpercentile(react_cluster, 25, axis=1)
            react_hi = np.nanpercentile(react_cluster, 75, axis=1)
        else:
            raise ValueError(f"spread must be 'std', 'se', 'ci95', or 'iqr', got '{spread}'")

        # Skip zero-coverage positions: NaN breaks the line/fill so they appear as gaps
        zero_cov = cov_mean == 0
        cov_mean = cov_mean.astype(float).copy()
        cov_lo = np.asarray(cov_lo, dtype=float).copy()
        cov_hi = np.asarray(cov_hi, dtype=float).copy()
        cov_mean[zero_cov] = np.nan
        cov_lo[zero_cov] = np.nan
        cov_hi[zero_cov] = np.nan

        label = f"{cluster} (n={n_cells})"
        color = colors[ci]

        # Plot coverage
        ax_cov.plot(x_positions, cov_mean, label=label, color=color, alpha=line_alpha)
        if variance_style == 'fill':
            ax_cov.fill_between(x_positions, cov_lo, cov_hi,
                                color=color, alpha=fill_alpha)
        else:
            ax_cov.errorbar(x_positions, cov_mean, yerr=[cov_mean - cov_lo, cov_hi - cov_mean],
                            color=color, alpha=line_alpha, capsize=2, fmt='none')

        # Plot reactivity
        ax_react.plot(x_positions, react_mean, label=label, color=color, alpha=line_alpha)
        if variance_style == 'fill':
            ax_react.fill_between(x_positions, react_lo, react_hi,
                                  color=color, alpha=fill_alpha)
        else:
            ax_react.errorbar(x_positions, react_mean, yerr=[react_mean - react_lo, react_hi - react_mean],
                              color=color, alpha=line_alpha, capsize=2, fmt='none')

    if log_coverage:
        ax_cov.set_yscale('log')
        ax_cov.set_ylabel('Coverage (log scale)')
    else:
        ax_cov.set_ylabel('Coverage')
    ax_cov.legend(title=cluster_col, fontsize='small')
    ax_react.set_ylabel('Normalized reactivity' if use_normalized else 'Reactivity')
    ax_react.set_xlabel('Position')
    ax_react.legend(title=cluster_col, fontsize='small')

    if title:
        fig.suptitle(title)
    elif gene is not None:
        fig.suptitle(f'{gene} profile by {cluster_col}')

    plt.tight_layout()
    return axes


def plot_reactivity_distribution(
    data: 'ShapeData',
    color_by: Optional[str] = None,
    use_normalized: bool = False,
    cells: Optional[List[str]] = None,
    gene: Optional[str] = None,
    pos_range: Optional[Tuple[int, int]] = None,
    ncols: int = 4,
    figsize: Optional[Tuple[int, int]] = None,
    palette: Optional[str] = 'Set2',
    fill: bool = True,
    sharex: bool = True,
    title: Optional[str] = None,
    **kwargs
) -> Any:
    """
    Plot the distribution of reactivity values for each cell as a grid of KDE plots.

    Creates one subplot per cell, each showing a kernel density estimate of that
    cell's reactivity values across positions (NaN values are dropped). When
    ``color_by`` is given, each cell's curve is colored by the cell's value in that
    metadata column, with a shared figure legend.

    Parameters:
    -----------
    data : ShapeData
        The ShapeData object to plot from. Must have a reactivity matrix
        (call calculate_reactivity() first if needed).
    color_by : str, optional
        Column name in cells metadata used to color each cell's KDE. If None,
        all curves use a single default color.
    use_normalized : bool
        If True, use normalized_reactivity; otherwise use raw reactivity
        (default: False).
    cells : list of str, optional
        Specific cell names to plot. If None, plots all cells.
    gene : str, optional
        Gene name to filter positions. Only used when genepos is a DataFrame
        with a 'gene' column.
    pos_range : tuple of (int, int), optional
        (start, end) position filter (inclusive). Uses genepos 'pos' values when
        available, otherwise row-index slicing. Overridden by ``gene`` where both apply.
    ncols : int
        Number of columns in the subplot grid (default: 4).
    figsize : tuple, optional
        Figure size (width, height). If None, auto-calculated as (4*ncols, 3*nrows).
    palette : str, optional
        Color palette name used to map ``color_by`` groups to colors (default: 'Set2').
    fill : bool
        Whether to fill the area under each KDE curve (default: True).
    sharex : bool
        Whether subplots share the same x-axis (default: True).
    title : str, optional
        Overall figure title.
    **kwargs
        Additional arguments passed to seaborn.kdeplot.

    Returns:
    --------
    numpy.ndarray of matplotlib.axes.Axes
        Array of the axes objects used (one per plotted cell).

    Example:
    --------
    >>> # One KDE per cell, colored by RBP
    >>> data.plot_reactivity_distribution(color_by='RBP')
    >>>
    >>> # Restrict to a gene and a subset of cells
    >>> data.plot_reactivity_distribution(color_by='RBP', gene='18S',
    ...                                   cells=['EWSR1', 'HNRNPA1'])
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    try:
        import seaborn as sns
    except ImportError:
        raise ImportError("seaborn is required for plot_reactivity_distribution. "
                          "Install with: pip install seaborn")

    # Select reactivity matrix
    if use_normalized:
        if data.normalized_reactivity is None:
            raise ValueError("normalized_reactivity is None. Call calculate_reactivity() "
                             "with normalize=True first.")
        reactivity_matrix = data.normalized_reactivity
        value_label = 'Normalized reactivity'
    else:
        if data.reactivity is None:
            raise ValueError("reactivity matrix is None. Call calculate_reactivity() first.")
        reactivity_matrix = data.reactivity
        value_label = 'Reactivity'

    # Determine position mask (same logic as plot_reactivity)
    n_pos = reactivity_matrix.shape[0]
    pos_mask = np.ones(n_pos, dtype=bool)

    if gene is not None and isinstance(data.genepos, pd.DataFrame) and 'gene' in data.genepos.columns:
        pos_mask &= data.genepos['gene'].values == gene

    if pos_range is not None and isinstance(data.genepos, pd.DataFrame) and 'pos' in data.genepos.columns:
        pos_vals = data.genepos['pos'].values
        start, end = pos_range
        pos_mask &= (pos_vals >= start) & (pos_vals <= end)
    elif pos_range is not None:
        start, end = pos_range
        idx_mask = np.zeros(n_pos, dtype=bool)
        idx_mask[start:min(end + 1, n_pos)] = True
        pos_mask &= idx_mask

    pos_indices = np.where(pos_mask)[0]
    if len(pos_indices) == 0:
        raise ValueError("No positions selected")

    react_sub = np.asarray(reactivity_matrix[pos_indices, :], dtype=float)

    # Resolve cell column indices
    cell_names = data.cell_names
    if cells is None:
        cell_indices = list(range(len(cell_names)))
    else:
        name_to_idx = {name: i for i, name in enumerate(cell_names)}
        missing = [c for c in cells if c not in name_to_idx]
        if missing:
            raise ValueError(f"Cells not found: {missing}")
        cell_indices = [name_to_idx[c] for c in cells]

    if len(cell_indices) == 0:
        raise ValueError("No cells selected")

    # Build color mapping from color_by metadata column
    cell_colors = None
    group_color = None
    if color_by is not None:
        if not data.cells_is_df:
            raise ValueError("cells must be a DataFrame to use color_by")
        if color_by not in data.cells.columns:
            raise ValueError(f"color_by column '{color_by}' not found in cells. "
                             f"Available: {list(data.cells.columns)}")
        group_values = data.cells[color_by].values
        unique_groups = list(pd.unique(group_values))
        colors = sns.color_palette(palette, n_colors=len(unique_groups))
        group_color = {g: colors[i] for i, g in enumerate(unique_groups)}
        cell_colors = {idx: group_color[group_values[idx]] for idx in cell_indices}

    # Grid layout
    n_plots = len(cell_indices)
    nrows = (n_plots + ncols - 1) // ncols
    if figsize is None:
        figsize = (4 * ncols, 3 * nrows)

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False, sharex=sharex)
    axes = axes.flatten()

    for i, col in enumerate(cell_indices):
        ax = axes[i]
        values = react_sub[:, col]
        values = values[np.isfinite(values)]

        color = cell_colors[col] if cell_colors is not None else 'steelblue'

        cell_name = cell_names[col]
        if color_by is not None:
            subtitle = f"{cell_name} ({data.cells[color_by].values[col]})"
        else:
            subtitle = str(cell_name)

        # KDE needs at least 2 distinct finite values
        if values.size < 2 or np.unique(values).size < 2:
            ax.text(0.5, 0.5, 'insufficient data', ha='center', va='center',
                    transform=ax.transAxes, fontsize='small', color='gray')
        else:
            sns.kdeplot(x=values, ax=ax, color=color, fill=fill, **kwargs)

        ax.set_title(subtitle, fontsize='small')
        ax.set_xlabel(value_label)

    # Hide unused subplots
    for i in range(n_plots, len(axes)):
        axes[i].set_visible(False)

    # Shared legend for color groups
    if color_by is not None and group_color is not None:
        handles = [mpatches.Patch(color=c, label=str(g)) for g, c in group_color.items()]
        fig.legend(handles=handles, title=color_by, loc='center left',
                   bbox_to_anchor=(1.0, 0.5), fontsize='small')

    if title:
        fig.suptitle(title, y=1.02)

    plt.tight_layout()
    return axes[:n_plots]


def plot_mutrate_distribution(
    data: 'ShapeData',
    control_cols: Optional[List[str]] = None,
    control_prefix: str = 'mean_mutrate_control_',
    color_by: Optional[str] = None,
    cells: Optional[List[str]] = None,
    gene: Optional[str] = None,
    pos_range: Optional[Tuple[int, int]] = None,
    ncols: int = 4,
    figsize: Optional[Tuple[int, int]] = None,
    palette: Optional[str] = 'Set2',
    fill: bool = True,
    sharex: bool = True,
    title: Optional[str] = None,
    **kwargs
) -> Any:
    """
    Plot the distribution of mutation rate for each cell as a grid of KDE plots.

    Creates one subplot per cell. Each subplot overlays:
    - the KDE of that cell's raw mutation rate (``data.mutrate``) across positions,
      using only positions where the cell's coverage > 0 (zero-coverage positions
      are dropped); and
    - one KDE per control mutation-rate vector stored as a column in ``data.genepos``
      (per-position control profiles, e.g. added via ``add_control_data``). The same
      control curves appear in every subplot as a shared baseline.

    Parameters:
    -----------
    data : ShapeData
        The ShapeData object to plot from. Must have a mutrate matrix.
    control_cols : list of str, optional
        Names of control columns in ``data.genepos`` to overlay. If None, all genepos
        columns whose names start with ``control_prefix`` are used.
    control_prefix : str
        Prefix used to auto-discover control columns when ``control_cols`` is None
        (default: 'mean_mutrate_control_').
    color_by : str, optional
        Column in cells metadata used to color each cell's mutrate KDE. If None, a
        single default color is used.
    cells : list of str, optional
        Specific cell names to plot. If None, plots all cells.
    gene : str, optional
        Gene name to filter positions (requires genepos DataFrame with 'gene' column).
    pos_range : tuple of (int, int), optional
        (start, end) position filter (inclusive). Uses genepos 'pos' values when
        available, otherwise row-index slicing. Overridden by ``gene`` where both apply.
    ncols : int
        Number of columns in the subplot grid (default: 4).
    figsize : tuple, optional
        Figure size. If None, auto-calculated as (4*ncols, 3*nrows).
    palette : str, optional
        Color palette name used to map ``color_by`` groups to colors (default: 'Set2').
    fill : bool
        Whether to fill the area under each cell's mutrate KDE (default: True). Control
        curves are always drawn unfilled (dashed) for comparison.
    sharex : bool
        Whether subplots share the same x-axis (default: True).
    title : str, optional
        Overall figure title.
    **kwargs
        Additional arguments passed to seaborn.kdeplot for the cell mutrate curve.

    Returns:
    --------
    numpy.ndarray of matplotlib.axes.Axes
        Array of the axes objects used (one per plotted cell).

    Example:
    --------
    >>> # One mutrate KDE per cell, plus every mean_mutrate_control_* baseline
    >>> data.plot_mutrate_distribution()
    >>>
    >>> # Restrict to a gene and a single control column, colored by RBP
    >>> data.plot_mutrate_distribution(control_cols=['mean_mutrate_control_rbp'],
    ...                                gene='18S', color_by='RBP')
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import matplotlib.lines as mlines
    from scipy import sparse
    try:
        import seaborn as sns
    except ImportError:
        raise ImportError("seaborn is required for plot_mutrate_distribution. "
                          "Install with: pip install seaborn")

    if data.mutrate is None:
        raise ValueError("mutrate matrix is None. ShapeData must be created with a "
                         "mutrate matrix.")

    # Determine position mask (same logic as plot_reactivity_distribution)
    n_pos = data.mutrate.shape[0]
    pos_mask = np.ones(n_pos, dtype=bool)

    if gene is not None and isinstance(data.genepos, pd.DataFrame) and 'gene' in data.genepos.columns:
        pos_mask &= data.genepos['gene'].values == gene

    if pos_range is not None and isinstance(data.genepos, pd.DataFrame) and 'pos' in data.genepos.columns:
        pos_vals = data.genepos['pos'].values
        start, end = pos_range
        pos_mask &= (pos_vals >= start) & (pos_vals <= end)
    elif pos_range is not None:
        start, end = pos_range
        idx_mask = np.zeros(n_pos, dtype=bool)
        idx_mask[start:min(end + 1, n_pos)] = True
        pos_mask &= idx_mask

    pos_indices = np.where(pos_mask)[0]
    if len(pos_indices) == 0:
        raise ValueError("No positions selected")

    # Slice mutrate and coverage to selected positions (keep sparse, densify per-column)
    mut_sub = data.mutrate[pos_indices, :]
    cov_sub = data.coverage[pos_indices, :]

    # Resolve control columns from genepos
    if not isinstance(data.genepos, pd.DataFrame):
        raise ValueError("genepos must be a DataFrame to plot control columns")
    if control_cols is None:
        control_cols = [c for c in data.genepos.columns if c.startswith(control_prefix)]
    else:
        missing = [c for c in control_cols if c not in data.genepos.columns]
        if missing:
            available = [c for c in data.genepos.columns if c.startswith(control_prefix)]
            raise ValueError(f"Control columns not found in genepos: {missing}. "
                             f"Available control columns: {available}")

    # Extract each control as a finite-filtered per-position vector + assign a color
    control_curves = []  # list of (name, values, color)
    control_palette = sns.color_palette('dark', n_colors=max(len(control_cols), 1))
    for ci, ccol in enumerate(control_cols):
        cvals = np.asarray(data.genepos[ccol].values, dtype=float)[pos_indices]
        cvals = cvals[np.isfinite(cvals)]
        control_curves.append((ccol, cvals, control_palette[ci]))

    # Resolve cell column indices
    cell_names = data.cell_names
    if cells is None:
        cell_indices = list(range(len(cell_names)))
    else:
        name_to_idx = {name: i for i, name in enumerate(cell_names)}
        missing = [c for c in cells if c not in name_to_idx]
        if missing:
            raise ValueError(f"Cells not found: {missing}")
        cell_indices = [name_to_idx[c] for c in cells]

    if len(cell_indices) == 0:
        raise ValueError("No cells selected")

    # Build color mapping from color_by metadata column
    cell_colors = None
    group_color = None
    if color_by is not None:
        if not data.cells_is_df:
            raise ValueError("cells must be a DataFrame to use color_by")
        if color_by not in data.cells.columns:
            raise ValueError(f"color_by column '{color_by}' not found in cells. "
                             f"Available: {list(data.cells.columns)}")
        group_values = data.cells[color_by].values
        unique_groups = list(pd.unique(group_values))
        colors = sns.color_palette(palette, n_colors=len(unique_groups))
        group_color = {g: colors[i] for i, g in enumerate(unique_groups)}
        cell_colors = {idx: group_color[group_values[idx]] for idx in cell_indices}

    # Grid layout
    n_plots = len(cell_indices)
    nrows = (n_plots + ncols - 1) // ncols
    if figsize is None:
        figsize = (4 * ncols, 3 * nrows)

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False, sharex=sharex)
    axes = axes.flatten()

    for i, col in enumerate(cell_indices):
        ax = axes[i]

        mut_col = mut_sub[:, col]
        cov_col = cov_sub[:, col]
        mut_vals = mut_col.toarray().ravel() if sparse.issparse(mut_col) else np.asarray(mut_col).ravel()
        cov_vals = cov_col.toarray().ravel() if sparse.issparse(cov_col) else np.asarray(cov_col).ravel()

        # Only positions with coverage > 0 and finite mutrate
        values = mut_vals[(cov_vals > 0) & np.isfinite(mut_vals)]

        color = cell_colors[col] if cell_colors is not None else 'steelblue'

        cell_name = cell_names[col]
        if color_by is not None:
            subtitle = f"{cell_name} ({data.cells[color_by].values[col]})"
        else:
            subtitle = str(cell_name)

        # KDE needs at least 2 distinct finite values
        if values.size < 2 or np.unique(values).size < 2:
            ax.text(0.5, 0.5, 'insufficient data', ha='center', va='center',
                    transform=ax.transAxes, fontsize='small', color='gray')
        else:
            sns.kdeplot(x=values, ax=ax, color=color, fill=fill, **kwargs)

        # Overlay control distributions (shared baseline, dashed, unfilled)
        for ccol, cvals, ccolor in control_curves:
            if cvals.size >= 2 and np.unique(cvals).size >= 2:
                sns.kdeplot(x=cvals, ax=ax, color=ccolor, fill=False,
                            linestyle='--', label=ccol)

        ax.set_title(subtitle, fontsize='small')
        ax.set_xlabel('Mutation rate')
        ax.set_xlim(-0.01, 0.07)
        # Suppress per-axes legend; a shared figure legend is built below
        legend = ax.get_legend()
        if legend is not None:
            legend.remove()

    # Hide unused subplots
    for i in range(n_plots, len(axes)):
        axes[i].set_visible(False)

    # Shared legend: color_by groups (if any) + control curves
    handles = []
    if color_by is not None and group_color is not None:
        handles.extend([mpatches.Patch(color=c, label=str(g))
                        for g, c in group_color.items()])
    handles.extend([mlines.Line2D([], [], color=ccolor, linestyle='--', label=ccol)
                    for ccol, _, ccolor in control_curves])
    if handles:
        legend_title = color_by if color_by is not None else None
        fig.legend(handles=handles, title=legend_title, loc='center left',
                   bbox_to_anchor=(1.0, 0.5), fontsize='small')

    if title:
        fig.suptitle(title, y=1.02)

    plt.tight_layout()
    return axes[:n_plots]
