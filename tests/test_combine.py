"""Tests for ShapeData.combine()."""

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from shape_data import ShapeData


# ── Helpers ──────────────────────────────────────────────────────────

def make_shapedata(n_positions=5, n_cells=3, seed=0, genepos_type='df'):
    """Create a small ShapeData with random sparse data."""
    rng = np.random.default_rng(seed)

    cov = sparse.random(n_positions, n_cells, density=0.6, format='csr',
                        dtype=np.int32, random_state=rng)
    cov.data = rng.integers(10, 500, size=cov.nnz).astype(np.int32)

    mut = cov.copy().astype(np.float32)
    mut.data = rng.uniform(0, 0.1, size=mut.nnz).astype(np.float32)

    cells = [f"CELL{i}" for i in range(n_cells)]

    if genepos_type == 'df':
        genepos = pd.DataFrame({
            'gene': ['18S'] * n_positions,
            'pos':  list(range(1, n_positions + 1)),
        })
        genepos.index = [f"18S.{p}" for p in genepos['pos']]
    else:
        genepos = [f"18S.{p}" for p in range(1, n_positions + 1)]

    return ShapeData(cov, mut, genepos, cells)


# ── Tests ────────────────────────────────────────────────────────────

class TestCombineBasic:
    """Basic combine functionality."""

    def test_combine_two_samples(self):
        sd1 = make_shapedata(n_positions=5, n_cells=3, seed=0)
        sd2 = make_shapedata(n_positions=5, n_cells=4, seed=1)

        combined = ShapeData.combine({'s1': sd1, 's2': sd2}, verbose=False)

        assert combined.n_positions == 5
        assert combined.n_cells == 3 + 4
        assert combined.coverage.shape == (5, 7)
        assert combined.mutrate.shape == (5, 7)

    def test_combine_three_samples(self):
        sd1 = make_shapedata(n_positions=5, n_cells=2, seed=0)
        sd2 = make_shapedata(n_positions=5, n_cells=3, seed=1)
        sd3 = make_shapedata(n_positions=5, n_cells=4, seed=2)

        combined = ShapeData.combine(
            {'a': sd1, 'b': sd2, 'c': sd3}, verbose=False
        )

        assert combined.n_cells == 2 + 3 + 4
        assert combined.n_positions == 5

    def test_combine_single_sample(self):
        sd = make_shapedata(n_positions=5, n_cells=3, seed=0)
        combined = ShapeData.combine({'only': sd}, verbose=False)

        assert combined.n_cells == 3
        assert combined.n_positions == 5


class TestCombineSampleColumn:
    """Sample column is correctly added."""

    def test_sample_labels(self):
        sd1 = make_shapedata(n_positions=5, n_cells=3, seed=0)
        sd2 = make_shapedata(n_positions=5, n_cells=4, seed=1)

        combined = ShapeData.combine({'nain3': sd1, 'dmso': sd2}, verbose=False)

        assert 'sample' in combined.cells.columns
        labels = combined.cells['sample'].tolist()
        assert labels == ['nain3'] * 3 + ['dmso'] * 4

    def test_custom_sample_column_name(self):
        sd1 = make_shapedata(n_positions=5, n_cells=2, seed=0)
        sd2 = make_shapedata(n_positions=5, n_cells=2, seed=1)

        combined = ShapeData.combine(
            {'a': sd1, 'b': sd2}, sample_column='condition', verbose=False
        )

        assert 'condition' in combined.cells.columns
        assert 'sample' not in combined.cells.columns


class TestCombineCellPrefix:
    """Cell barcode prefixing."""

    def test_prefix_enabled(self):
        sd1 = make_shapedata(n_positions=5, n_cells=2, seed=0)
        sd2 = make_shapedata(n_positions=5, n_cells=2, seed=1)

        combined = ShapeData.combine(
            {'s1': sd1, 's2': sd2}, cell_prefix=True, verbose=False
        )

        names = combined.cell_names
        assert all(n.startswith('s1_') for n in names[:2])
        assert all(n.startswith('s2_') for n in names[2:])

    def test_prefix_disabled(self):
        sd1 = make_shapedata(n_positions=5, n_cells=2, seed=0)
        sd2 = make_shapedata(n_positions=5, n_cells=2, seed=1)

        combined = ShapeData.combine(
            {'s1': sd1, 's2': sd2}, cell_prefix=False, verbose=False
        )

        names = combined.cell_names
        assert not any(n.startswith('s1_') or n.startswith('s2_') for n in names)


class TestCombineMatrixValues:
    """Matrix data is preserved correctly after combine."""

    def test_coverage_values_preserved(self):
        sd1 = make_shapedata(n_positions=5, n_cells=2, seed=10)
        sd2 = make_shapedata(n_positions=5, n_cells=3, seed=20)

        combined = ShapeData.combine({'a': sd1, 'b': sd2}, verbose=False)

        # First 2 columns should match sd1
        np.testing.assert_array_equal(
            combined.coverage[:, :2].toarray(),
            sd1.coverage.toarray()
        )
        # Last 3 columns should match sd2
        np.testing.assert_array_equal(
            combined.coverage[:, 2:].toarray(),
            sd2.coverage.toarray()
        )

    def test_mutrate_values_preserved(self):
        sd1 = make_shapedata(n_positions=5, n_cells=2, seed=10)
        sd2 = make_shapedata(n_positions=5, n_cells=3, seed=20)

        combined = ShapeData.combine({'a': sd1, 'b': sd2}, verbose=False)

        np.testing.assert_array_almost_equal(
            combined.mutrate[:, :2].toarray(),
            sd1.mutrate.toarray()
        )
        np.testing.assert_array_almost_equal(
            combined.mutrate[:, 2:].toarray(),
            sd2.mutrate.toarray()
        )


class TestCombineGenepos:
    """Genepos handling."""

    def test_genepos_dataframe_preserved(self):
        sd1 = make_shapedata(n_positions=5, n_cells=2, seed=0, genepos_type='df')
        sd2 = make_shapedata(n_positions=5, n_cells=3, seed=1, genepos_type='df')

        combined = ShapeData.combine({'a': sd1, 'b': sd2}, verbose=False)

        assert isinstance(combined.genepos, pd.DataFrame)
        assert list(combined.genepos.index) == list(sd1.genepos.index)

    def test_genepos_list_preserved(self):
        sd1 = make_shapedata(n_positions=5, n_cells=2, seed=0, genepos_type='list')
        sd2 = make_shapedata(n_positions=5, n_cells=3, seed=1, genepos_type='list')

        combined = ShapeData.combine({'a': sd1, 'b': sd2}, verbose=False)

        assert isinstance(combined.genepos, list)
        assert combined.genepos == sd1.genepos


class TestCombineMismatchedPositions:
    """Handling of samples with different positions."""

    def test_different_positions_uses_common(self):
        sd1 = make_shapedata(n_positions=5, n_cells=2, seed=0)

        # Create sd2 with different positions (only 3 overlap with sd1)
        cov2 = sparse.random(4, 3, density=0.6, format='csr', dtype=np.int32)
        cov2.data = np.ones(cov2.nnz, dtype=np.int32) * 100
        mut2 = cov2.copy().astype(np.float32)
        mut2.data = np.ones(mut2.nnz, dtype=np.float32) * 0.05
        genepos2 = pd.DataFrame({
            'gene': ['18S'] * 4,
            'pos':  [2, 3, 4, 99],  # positions 2,3,4 overlap; 99 is extra
        })
        genepos2.index = [f"18S.{p}" for p in genepos2['pos']]
        sd2 = ShapeData(cov2, mut2, genepos2, ['X1', 'X2', 'X3'])

        # sd1 has positions 1-5, sd2 has 2,3,4,99 → common = 2,3,4
        combined = ShapeData.combine({'a': sd1, 'b': sd2}, verbose=False)

        assert combined.n_positions == 3  # only common positions
        assert combined.n_cells == 2 + 3
        assert list(combined.genepos.index) == ['18S.2', '18S.3', '18S.4']


class TestCombineEdgeCases:
    """Edge cases."""

    def test_empty_dict_raises(self):
        with pytest.raises(ValueError, match="empty"):
            ShapeData.combine({})

    def test_cells_is_dataframe(self):
        sd1 = make_shapedata(n_positions=5, n_cells=2, seed=0)
        sd2 = make_shapedata(n_positions=5, n_cells=3, seed=1)

        combined = ShapeData.combine({'a': sd1, 'b': sd2}, verbose=False)

        assert combined.cells_is_df
        assert combined.cells.index.name == 'cell'

    def test_mutrate_none_handled(self):
        """If one sample has mutrate=None, combined mutrate should be None."""
        sd1 = make_shapedata(n_positions=5, n_cells=2, seed=0)

        # Create sd2 with mutrate=None
        cov2 = sparse.random(5, 3, density=0.6, format='csr', dtype=np.int32)
        cov2.data = np.ones(cov2.nnz, dtype=np.int32) * 100
        genepos2 = sd1.genepos.copy()
        sd2 = ShapeData(cov2, None, genepos2, ['X1', 'X2', 'X3'])

        combined = ShapeData.combine({'a': sd1, 'b': sd2}, verbose=False)
        assert combined.mutrate is None


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
