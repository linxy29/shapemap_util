"""Tabix-backed bulk access for *.mutrate.txt.bgz files.

Each *.mutrate.txt.bgz file (built by ``bgzip_tabix_mutrate.sh``) is a
13-column TSV with no header. Column 1 is the gene name, columns 2/3 are the
start/end positions of the per-base record. The accompanying .tbi lets
``pysam.TabixFile.fetch`` pull rows for one gene in O(log N).
"""

from __future__ import annotations

import os
from typing import Iterable, Sequence, Union

import pandas as pd
import pysam

MUTRATE_COLUMNS = ["gene", "pos", "pos+1", "gene.position", "mutrate", "strand", "coverage", "coverage_withIndel", "mutant", "normalized_cov", "g_readcount", "refnt", "detail"]

def _as_gene_list(genes: Union[str, Iterable[str], None]) -> list[str] | None:
    if genes is None:
        return None
    if isinstance(genes, str):
        return [genes]
    return list(genes)


def read_genes(
    bgz_path: str,
    genes: Union[str, Iterable[str], None],
    columns: Sequence[str] | None = None,
    min_coverage: float | None = None,
    coverage_col: str = "coverage",
) -> pd.DataFrame:
    """Fetch all rows for ``genes`` from a tabix-indexed mutrate file.

    Parameters
    ----------
    bgz_path : path to a ``*.mutrate.txt.bgz`` file (``.tbi`` must sit next to it).
    genes    : single gene name, iterable of gene names, or ``None`` to read
               every gene in the file.
    columns  : optional subset of ``MUTRATE_COLUMNS`` to return; default is all.
    min_coverage : if set, drop rows where ``coverage_col`` < this value.
    coverage_col : which column to threshold on (``coverage`` or
                   ``coverage_withIndel``).

    Returns
    -------
    DataFrame with one row per genomic position, ordered as in the file.
    Genes that aren't present produce no rows (no error).
    """
    tbi = bgz_path + ".tbi"
    if not os.path.exists(bgz_path):
        raise FileNotFoundError(bgz_path)
    if not os.path.exists(tbi):
        raise FileNotFoundError(f"missing tabix index: {tbi}")

    gene_list = _as_gene_list(genes)
    rows: list[list[str]] = []
    with pysam.TabixFile(bgz_path) as tbx:
        if gene_list is None:
            for line in tbx.fetch():
                rows.append(line.split("\t"))
        else:
            available = set(tbx.contigs)
            for g in gene_list:
                if g not in available:
                    continue
                for line in tbx.fetch(g):
                    rows.append(line.split("\t"))

    df = pd.DataFrame(rows, columns=MUTRATE_COLUMNS)
    if df.empty:
        return df if columns is None else df[list(columns)]

    numeric_cols = ["pos", "pos+1", "mutrate", "coverage",
                    "coverage_withIndel", "mutant", "normalized_cov", "g_readcount"]
    for c in numeric_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    if min_coverage is not None:
        if coverage_col not in df.columns:
            raise KeyError(f"coverage_col {coverage_col!r} not in {list(df.columns)}")
        df = df[df[coverage_col] >= min_coverage].reset_index(drop=True)

    return df if columns is None else df[list(columns)]


def read_genes_multi(
    bgz_paths: dict[str, str] | Sequence[str],
    genes: Union[str, Iterable[str], None],
    columns: Sequence[str] | None = None,
    source_col: str = "source",
    min_coverage: float | None = None,
    coverage_col: str = "coverage",
) -> pd.DataFrame:
    """Fetch ``genes`` from many bgz files and concatenate the results.

    Parameters
    ----------
    bgz_paths : either a {label: path} dict (label written into ``source_col``)
                or a sequence of paths (file basename used as the label).
    genes     : single gene name, iterable of gene names, or ``None`` to read
                every gene in each file.
    columns   : optional subset of ``MUTRATE_COLUMNS`` to return (the
                ``source_col`` is always added).
    source_col: name of the column that records which file a row came from.
    min_coverage : if set, drop rows where ``coverage_col`` < this value.
    coverage_col : which column to threshold on (``coverage`` or
                   ``coverage_withIndel``).

    Returns
    -------
    Long-format DataFrame with all matching rows, tagged by ``source_col``.
    """
    if isinstance(bgz_paths, dict):
        items = list(bgz_paths.items())
    else:
        items = [(os.path.basename(p), p) for p in bgz_paths]

    parts = []
    for label, path in items:
        sub = read_genes(
            path, genes,
            columns=columns,
            min_coverage=min_coverage,
            coverage_col=coverage_col,
        )
        if sub.empty:
            continue
        sub.insert(0, source_col, label)
        parts.append(sub)

    if not parts:
        out_cols = (list(columns) if columns is not None else MUTRATE_COLUMNS)
        return pd.DataFrame(columns=[source_col] + out_cols)
    return pd.concat(parts, ignore_index=True)
