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

### From Saved Files

```python
# Load from pickle file (single file)
data = ShapeData.from_pickle('my_data.pkl')

# Load from directory (multiple files)
data = ShapeData.load('./my_shapedata/')
```

### Combining Multiple ShapeData Objects

Horizontally concatenate cells from several ShapeData objects that share the same
positions. A `sample` column is added to track origin; cell barcodes can be prefixed
with the sample name to avoid collisions.

```python
combined = ShapeData.combine(
    {'nain3_1': sd1, 'dmso_1': sd2, 'dmso_2': sd3},
    sample_column='sample',
    cell_prefix=True,
)
print(combined.cells['sample'].value_counts())
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

### Filter by Position (genepos) Metadata

```python
# Filter using metadata column ranges
data = data.filter_genepos(
    genepos_filters={
        'gene': '18S',                    # exact match
        'pos': (100, 500),                # 100 <= pos <= 500
        'n_cells': (10, None),            # n_cells >= 10
    }
)

# Filter by multiple genes
data = data.filter_genepos(
    genepos_filters={
        'gene': ['18S', '28S']            # gene in ['18S', '28S']
    }
)

# Filter by specific position indices
data = data.filter_genepos(
    genepos_indices=[0, 1, 2, 5, 10]
)

# Filter using a boolean mask
mask = data.genepos['n_cells'] >= 5
data = data.filter_genepos(genepos_indices=mask.values)

# Filter using a callable
data = data.filter_genepos(
    genepos_indices=lambda gp: gp['gene'] == '18S'
)

# Combine filters and indices
data = data.filter_genepos(
    genepos_filters={'gene': '18S'},
    genepos_indices=lambda gp: gp['n_cells'] >= 5
)
```

### Subset by Position Indices

```python
pos_indices = [0, 5, 10, 15, 20]
data = data.subset_by_genepos(pos_indices)
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

### Calculate Reactivity from a Matched Control ShapeData

When treatment and control are stored as separate ShapeData objects (e.g., 2A3 vs DMSO
metabins) with a one-to-one mapping between them, subtract per-cell control mutrate:

```python
# control_data.cells must have a column (default 'ref_metabin') whose values
# match cell names (index) in `data`.
data.calculate_reactivity_from_control(
    control_data=dmso_metabins,
    ref_column='ref_metabin',
    store_as='reactivity',
    normalize=True,                  # also produce normalized_reactivity
    normalize_method='box',          # 'box', 'zscore', or 'minmax'
)

print(data.reactivity)
print(data.normalized_reactivity)
```

### Calculate AUROC vs Known Secondary Structure

For each cell, score reactivity as a predictor of unpaired vs paired bases given a
dot-bracket structure. Result is saved as a column in cells metadata.

```python
# Dot-bracket string (or path to a file containing it)
dot_bracket_18S = "(((...)))..."

data.calculate_auroc(
    gene='18S',
    dot_bracket=dot_bracket_18S,
    skip_positions=[1, 2, 3],       # 1-based positions to exclude
    save_as='AUROC_18S',            # default: f'AUROC_{gene}'
)

print(data.cells[['leiden', 'AUROC_18S']].head())
```

### Differential Reactivity Tests

Two complementary tests across positions, both returning a tidy DataFrame with
multiple-testing-corrected p-values.

**Linear model** (multiple predictors, continuous or categorical):

```python
results = data.differential_reactivity_lm(
    predictors=['diet', 'region'],   # cells columns; categoricals dummy-coded
    gene='18S',                       # or list of genes, or None for all
    use_normalized=False,
    min_cells=10,
    correction_method='fdr_bh',
)
# Columns: gene, pos, variable, coefficient, std_err, t_statistic,
#          pvalue, pvalue_adj, n_cells, r_squared
```

**Wilcoxon rank-sum** (two-group comparison):

```python
results = data.differential_reactivity_wilcoxon(
    group_column='diet',
    group1='HFD', group2='CTRL',
    gene=['18S', '28S'],
    alternative='two-sided',
)
# Columns: gene, pos, statistic, pvalue, pvalue_adj, effect_size,
#          mean_group1, mean_group2, diff_mean, log2_fold_change,
#          n_group1, n_group2
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

## Aggregating Cells into Metabins / Metacells

Aggregate groups of cells into a smaller number of "meta" cells. Coverage and mutcount
are summed across the group; mutrate is recomputed as `sum(mutcount) / sum(coverage)`
(or as a coverage-weighted average when mutcount is unavailable). Each method returns a
new `ShapeData` whose cells are the aggregated groups.

### Group Cells by a Metadata Column (`create_metacells`)

The simplest grouping: every cell sharing the same value in `group_col` becomes one
metacell. Cells with `NaN` in `group_col` are dropped (count is logged when verbose).

```python
# One metacell per unique value of cells['leiden']
mc = data.create_metacells(group_col='leiden')

# Output cells DataFrame columns:
#   - 'leiden'        : the group value
#   - 'cellbarcodes'  : list of original cell IDs in this metacell
#   - 'n_cells'       : number of cells aggregated
print(mc.cells.head())
print(mc.shape)  # (n_positions, n_unique_leiden_values)

# If your cell IDs live in a column other than 'bin100_cellID' (default),
# specify cell_id_col; if the column is absent, cells.index is used.
mc = data.create_metacells(group_col='sample', cell_id_col='cell_barcode')
```

### Group Spatial Bins by Cluster + Spatial Bisection (`create_metabins`)

For spatial data: within each cluster, recursively partition bins along the longer
spatial axis into balanced groups of approximately `n_bins` bins each.

```python
mb = data.create_metabins(
    n_bins=20,                  # target bins per metabin
    cluster_col='leiden',
    x_col='x', y_col='y',
    cell_id_col='bin100_cellID',
)
# Output cells contain: leiden, x, y (medoid), bin100_ids, n_bins
```

### Group Bins from an Explicit Dictionary (`create_metabins_from_dict`)

When you already know which bins belong to each metabin:

```python
metabin_dict = {
    'group_A': ['bin_1', 'bin_2', 'bin_3'],
    'group_B': ['bin_4', 'bin_5'],
}
mb = data.create_metabins_from_dict(metabin_dict)

# Optionally attach extra metadata per metabin:
mb = data.create_metabins_from_dict(
    metabin_dict,
    extra_meta={'group_A': {'cluster': 'K562'}, 'group_B': {'cluster': 'HEK'}},
)
```

### Reuse a Reference Grouping on a Target Sample (`create_metabins_from_mapping`)

Apply the metabin structure of a reference ShapeData (e.g., f2a3) to a target sample
(e.g., DMSO) via a bin-to-bin mapping table:

```python
target_mb = dmso_data.create_metabins_from_mapping(
    ref_metabin_data=f2a3_metabins,
    mapping_df=bin_mapping_df,
    ref_bin_col='f2a3_bin',
    target_bin_col='mapped_dmso_bin',
)
```

## Save and Load

### Single File (Pickle)

The simplest way to save and load ShapeData:

```python
# Save to a single pickle file
data.to_pickle('my_data.pkl')

# Load from pickle file
data = ShapeData.from_pickle('my_data.pkl')
```

### Directory (Multiple Files)

Save to a directory with separate files for each component (allows inspection):

```python
# Save to directory
data.save('./my_shapedata/')
# Creates:
#   my_shapedata/
#   ├── coverage.npz      (sparse matrix)
#   ├── mutrate.npz       (sparse matrix)
#   ├── reactivity.npz    (dense array, if exists)
#   ├── genepos.parquet   (DataFrame)
#   └── cells.parquet     (DataFrame)

# Load from directory
data = ShapeData.load('./my_shapedata/')
```

| Method | Format | Best for |
|--------|--------|----------|
| `to_pickle()` | Single `.pkl` file | Quick save/load, sharing |
| `save()` | Directory with multiple files | Inspectable, version control friendly |

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

### Coverage / Reactivity Profile by Cluster

Two stacked subplots: mean coverage (top) and mean reactivity (bottom) along positions,
one line per cluster, with shaded variance bands. Requires `data.reactivity` to exist
(call `calculate_reactivity` or `calculate_reactivity_from_control` first).

```python
# Plot full gene profile by cluster
axes = data.plot_reactivity(gene='18S', cluster_col='leiden')

# Plot a position range with specific clusters and smoothing
axes = data.plot_reactivity(
    pos_range=(100, 300),
    clusters=['0', '1', '2'],
    use_normalized=True,
    spread='se',                  # 'se', 'std', or 'ci'
    smoothing=3,                  # rolling mean window
    palette='Set2',
    figsize=(12, 6),
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
| `from_pickle()` | Load from pickle file |
| `load()` | Load from directory |
| `combine()` | Concatenate multiple ShapeData objects (shared positions) |
| `to_pickle()` | Save to pickle file |
| `save()` | Save to directory |
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
| `filter_genepos()` | Filter by position metadata or indices |
| `subset_by_cells()` | Subset to specific cell list |
| `subset_by_genepos()` | Subset to specific position indices |

### Analysis Functions

| Method | Description |
|--------|-------------|
| `calculate_reactivity()` | Reactivity from per-cluster control columns |
| `calculate_reactivity_from_control()` | Reactivity from a matched control ShapeData |
| `calculate_cell_correlation()` | Pairwise cell correlation |
| `calculate_auroc()` | AUROC of reactivity vs dot-bracket structure |
| `differential_reactivity_lm()` | Linear-model differential reactivity test |
| `differential_reactivity_wilcoxon()` | Wilcoxon rank-sum differential reactivity test |
| `get_cell_stats()` | Statistics per cell |
| `get_position_stats()` | Statistics per position |

### Plot Functions

| Method | Description |
|--------|-------------|
| `plot_violin()` | Violin plot with dot overlay |
| `plot_violin_multi()` | Multiple violin plots in grid |
| `plot_reactivity()` | Coverage + reactivity profiles by cluster |

### Metabin / Metacell Functions

| Method | Description |
|--------|-------------|
| `create_metacells()` | Group cells sharing a metadata column value into metacells |
| `create_metabins()` | Group spatial bins per cluster via balanced bisection |
| `create_metabins_from_dict()` | Aggregate bins from an explicit `name -> [bin_ids]` dict |
| `create_metabins_from_mapping()` | Reuse a reference metabin grouping on a target sample |
