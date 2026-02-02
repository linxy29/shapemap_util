# ShapeData Package

A Python package for handling single-cell SHAPE-MaP data with sparse matrices and metadata management.

## Installation

The package is located in the `shape_data/` directory. Import it directly:

```python
from shape_data import ShapeData
```

## Quick Start

```python
from shape_data import ShapeData
import pandas as pd

# Load data from VCF files
data = ShapeData.from_vcf_pair("rRNA.vcf.gz", "steptwo.vcf.gz")

# Filter by thresholds
data = data.filter_thresholds(coverage_threshold=1500, mutrate_threshold=0.2)

# Add metadata
data.add_genepos_counts()
data.join_cell_metadata(clustering_df, columns='leiden', on='cell')

print(data)
```

## Creating ShapeData

### From VCF Files

```python
# From paired rRNA and steptwo VCF files
data = ShapeData.from_vcf_pair(
    rRNA_path="path/to/rRNA.vcf.gz",
    steptwo_path="path/to/steptwo.vcf.gz",
    progress_interval=100000
)

# From cellSNP VCF file
data = ShapeData.from_cellsnp_vcf(
    vcf_path="path/to/cellSNP.cells.vcf.bgz.gz",
    progress_interval=100000
)
```

### From Matrices Directly

```python
from scipy import sparse
import pandas as pd

# Create sparse matrices (positions x cells)
coverage_sparse = sparse.csr_matrix(coverage_array)
mutrate_sparse = sparse.csr_matrix(mutrate_array)

# Position identifiers (list or DataFrame)
genepos = ['18S.100', '18S.101', '18S.102', ...]  # or DataFrame

# Cell identifiers (list or DataFrame)
cells = ['AAACCTGCAGTAACGG-1', 'AAACCTGCATCACGTA-1', ...]  # or DataFrame

data = ShapeData(coverage_sparse, mutrate_sparse, genepos, cells)
```

## Data Structure

```python
data.coverage    # scipy.sparse.csr_matrix (positions x cells)
data.mutrate     # scipy.sparse.csr_matrix (positions x cells)
data.genepos     # list or pd.DataFrame - position metadata
data.cells       # list or pd.DataFrame - cell metadata

# Properties
data.shape       # (n_positions, n_cells)
data.n_positions # Number of positions
data.n_cells     # Number of cells
data.cell_names  # List of cell barcodes
data.cells_is_df # True if cells is a DataFrame
```

## Metadata Management

### Converting to DataFrames

```python
# Convert cells to DataFrame (for adding metadata)
data = data.to_cells_df()

# Convert genepos to DataFrame (parses 'gene.pos' format)
data = data.to_genepos_df()
# Result: DataFrame with 'gene' and 'pos' columns, 'gene.pos' as index
```

### Cell Metadata

```python
# Add n_genepos column (count of non-zero positions per cell)
data.add_genepos_counts()

# Join external cell metadata
clustering_df = pd.read_csv("clustering.csv")
# clustering_df must have a column matching cell barcodes (default: 'cell')

data.join_cell_metadata(
    metadata_df=clustering_df,
    columns=['leiden', 'sample', 'Total'],  # columns to join
    on='cell',                               # join key in metadata_df
    how='left',                              # join type
    inplace=True                             # modify in place
)

# Access cell metadata
print(data.cells.head())
print(data.cells.columns.tolist())
```

### Position (genepos) Metadata

```python
# First ensure genepos is a DataFrame
data = data.to_genepos_df()

# Join external position metadata
pos_meta = pd.read_csv("position_annotation.csv")
# pos_meta must have 'gene' and 'pos' columns

data.join_genepos_metadata(
    metadata_df=pos_meta,
    columns=['ref', 'annotation', 'structure'],  # columns to join
    on=['gene', 'pos'],                          # join keys
    how='left',
    inplace=True
)

# Access position metadata
print(data.genepos.head())
```

### Adding Control Data (e.g., DMSO)

```python
# Control DataFrame format:
# | sample     | gene | pos | coverage | mutrate |
# |------------|------|-----|----------|---------|
# | dmso_K562  | 18S  | 100 | 5000     | 0.01    |
# | dmso_HEK   | 18S  | 100 | 4500     | 0.02    |

dmso_df = pd.read_csv("dmso_mutrate.csv")

data.add_control_data(
    control_df=dmso_df,
    cluster_col='sample',        # column with cluster identifiers
    cluster_prefix='dmso_',      # prefix to remove (dmso_K562 -> K562)
    coverage_col='coverage',
    mutrate_col='mutrate',
    on=['gene', 'pos'],
    control_suffix='control',    # suffix for new columns
    verbose=True
)

# Result: adds columns like 'coverage_control_K562', 'mutrate_control_K562'
```

## Filtering

### Filter by Thresholds

```python
# Set values to 0 where coverage < threshold or mutrate > threshold
data = data.filter_thresholds(
    coverage_threshold=1500,
    mutrate_threshold=0.2
)
```

### Filter by Cell Metadata

```python
# Filter using metadata column ranges
data = data.filter_cells(
    cell_filters={
        'Total': (100000, None),      # Total >= 100000
        'n_genepos': (50, 500),       # 50 <= n_genepos <= 500
        'leiden': ['0', '1', '2']     # leiden in ['0', '1', '2']
    }
)

# Filter by specific cell barcodes
data = data.filter_cells(
    cell_indices=['AAACCTGCAGTAACGG-1', 'AAACCTGCATCACGTA-1']
)

# Filter using a callable
data = data.filter_cells(
    cell_indices=lambda df: df['leiden'] == '0'
)
```

### Subset by Cell List

```python
cell_list = ['barcode1', 'barcode2', 'barcode3']
data = data.subset_by_cells(cell_list)
```

## Analysis

### Calculate Reactivity

```python
# Subtract cluster-specific control mutrate from each cell
data.calculate_reactivity(
    cluster_column='leiden',              # cluster assignment column
    control_prefix='mutrate_control_',    # prefix for control columns
    store_as='reactivity',                # attribute name for result
    drop_mutrate=False,
    verbose=True
)

# Access reactivity matrix
print(data.reactivity)  # sparse matrix
```

### Calculate Cell Correlation

```python
# Pairwise correlation between cells within clusters
corr_df = data.calculate_cell_correlation(
    gene='18S',
    cluster_column='leiden',
    method='pearson',              # 'pearson', 'spearman', or 'kendall'
    min_shared_positions=100,
    save_mean_as='mean_corr',      # save mean correlation to cells metadata
    verbose=True
)
```

### Get Statistics

```python
# Per-cell statistics
cell_stats = data.get_cell_stats()
# Returns: n_genepos, mean_coverage, mean_mutrate per cell

# Per-position statistics
pos_stats = data.get_position_stats()
# Returns: n_cells, mean_coverage, mean_mutrate per position
```

## Plotting

### Single Violin Plot

```python
# Plot single column
ax = data.plot_violin(
    column='n_genepos',
    source='cells',           # 'cells', 'genepos', or 'auto'
    groupby='leiden',         # optional grouping
    title='Positions per Cell',
    figsize=(10, 6)
)

# Plot multiple columns
ax = data.plot_violin(
    column=['n_genepos', 'mean_corr'],
    groupby='leiden',
    ncols=2
)
```

### Multiple Violin Plots

```python
fig = data.plot_violin_multi(
    columns=['coverage_control_K562', 'coverage_control_HEK'],
    source='genepos',
    groupby='gene',
    ncols=2,
    figsize=(12, 8),
    title='Control Coverage by Gene'
)
```

## Complete Workflow Example

```python
from shape_data import ShapeData
import pandas as pd

# 1. Load data
data = ShapeData.from_vcf_pair(
    "data/rRNA.vcf.gz",
    "data/steptwo.vcf.gz"
)
print(f"Loaded: {data}")

# 2. Filter by quality thresholds
data = data.filter_thresholds(
    coverage_threshold=1500,
    mutrate_threshold=0.2
)

# 3. Convert to DataFrames for metadata
data = data.to_genepos_df()

# 4. Add cell metadata
data.add_genepos_counts()

clustering = pd.read_csv("clustering.csv")
data.join_cell_metadata(clustering, columns=['leiden', 'sample', 'Total'], on='cell')

# 5. Filter cells by metadata
data = data.filter_cells(cell_filters={
    'Total': (100000, None),
    'n_genepos': (50, None)
})

# 6. Add control data
dmso = pd.read_csv("dmso_mutrate.csv")
data.add_control_data(dmso, cluster_prefix='dmso_')

# 7. Calculate reactivity
data.calculate_reactivity(cluster_column='leiden')

# 8. Analyze correlations
corr_df = data.calculate_cell_correlation(
    gene='18S',
    cluster_column='leiden',
    save_mean_as='mean_corr'
)

# 9. Filter by correlation
data = data.filter_cells(cell_filters={'mean_corr': (0.5, None)})

# 10. Visualize
data.plot_violin('mean_corr', groupby='leiden', title='Cell Correlation by Cluster')

print(f"Final dataset: {data}")
```

## API Reference

### ShapeData Class

| Method | Description |
|--------|-------------|
| `from_vcf_pair()` | Create from paired VCF files |
| `from_cellsnp_vcf()` | Create from cellSNP VCF |
| `copy()` | Deep copy |
| `to_cells_df()` | Convert cells to DataFrame |
| `to_genepos_df()` | Convert genepos to DataFrame |

### Metadata Functions

| Method | Description |
|--------|-------------|
| `add_genepos_counts()` | Add n_genepos column to cells |
| `join_cell_metadata()` | Join external cell metadata |
| `join_genepos_metadata()` | Join external position metadata |
| `add_control_data()` | Add control (DMSO) data |

### Filter Functions

| Method | Description |
|--------|-------------|
| `filter_thresholds()` | Filter by coverage/mutrate thresholds |
| `filter_cells()` | Filter by cell metadata or indices |
| `subset_by_cells()` | Subset to specific cell list |

### Analysis Functions

| Method | Description |
|--------|-------------|
| `calculate_reactivity()` | Calculate reactivity matrix |
| `calculate_cell_correlation()` | Pairwise cell correlation |
| `get_cell_stats()` | Statistics per cell |
| `get_position_stats()` | Statistics per position |

### Plot Functions

| Method | Description |
|--------|-------------|
| `plot_violin()` | Violin plot with dot overlay |
| `plot_violin_multi()` | Multiple violin plots in grid |
