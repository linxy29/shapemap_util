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
