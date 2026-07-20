"""Tests for oncocartograph.preprocessing.rna_seq.

Fixtures mirror the real GDC STAR-Counts file format (leading #-commented
gene-model line, header row, four N_* alignment-summary rows, then gene
rows) confirmed against downloaded files. The VST test uses a small
synthetic Poisson count matrix -- large enough for pydeseq2's dispersion
trend fitting to run, not asserting exact DESeq2 internals.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from oncocartograph.preprocessing.rna_seq import (
    build_counts_matrix,
    filter_low_expression,
    normalize_and_vst,
    read_star_counts,
)

_STAR_COUNTS_CONTENT = (
    "# gene-model: GENCODE v36\n"
    "gene_id\tgene_name\tgene_type\tunstranded\tstranded_first\tstranded_second"
    "\ttpm_unstranded\tfpkm_unstranded\tfpkm_uq_unstranded\n"
    "N_unmapped\t\t\t100\t100\t100\t\t\t\n"
    "N_multimapping\t\t\t200\t200\t200\t\t\t\n"
    "N_noFeature\t\t\t300\t300\t300\t\t\t\n"
    "N_ambiguous\t\t\t400\t400\t400\t\t\t\n"
    "ENSG00000000001.1\tGENE1\tprotein_coding\t1000\t500\t500\t10.0\t9.0\t9.5\n"
    "ENSG00000000002.1\tGENE2\tprotein_coding\t5\t2\t3\t0.1\t0.05\t0.06\n"
)


@pytest.fixture
def star_counts_file(tmp_path: Path) -> Path:
    """Write a synthetic STAR-Counts file and return its path."""
    path = tmp_path / "sample.rna_seq.augmented_star_gene_counts.tsv"
    path.write_text(_STAR_COUNTS_CONTENT)
    return path


def test_read_star_counts_drops_star_summary_rows(star_counts_file: Path) -> None:
    """N_unmapped/N_multimapping/N_noFeature/N_ambiguous must not appear as genes."""
    counts = read_star_counts(star_counts_file)

    assert set(counts.index) == {"ENSG00000000001.1", "ENSG00000000002.1"}
    assert counts["ENSG00000000001.1"] == 1000
    assert counts["ENSG00000000002.1"] == 5


def test_build_counts_matrix_combines_patients(tmp_path: Path) -> None:
    """The matrix must be gene-indexed with one raw-count column per patient."""
    path_a = tmp_path / "a.tsv"
    path_a.write_text(_STAR_COUNTS_CONTENT)
    path_b = tmp_path / "b.tsv"
    path_b.write_text(_STAR_COUNTS_CONTENT)

    matrix = build_counts_matrix({"case-a": path_a, "case-b": path_b})

    assert list(matrix.columns) == ["case-a", "case-b"]
    assert matrix.loc["ENSG00000000001.1", "case-a"] == 1000


def test_filter_low_expression_drops_low_count_genes() -> None:
    """A gene below the CPM threshold in most samples must be dropped."""
    # Library sizes are ~10,000,000 (dominated by high_expr), so low_expr's
    # raw counts of 0-2 correspond to CPM well under 1.0 in every sample.
    counts = pd.DataFrame(
        {
            "case-a": [10_000_000, 1, 10_000_000],
            "case-b": [10_000_000, 0, 10_000_000],
            "case-c": [10_000_000, 2, 0],
        },
        index=["high_expr", "low_expr", "sparse_high_expr"],
    )

    filtered = filter_low_expression(counts, min_cpm=1.0, min_fraction_samples=0.5)

    assert "low_expr" not in filtered.index
    assert "high_expr" in filtered.index


def test_filter_low_expression_keeps_genes_expressed_in_enough_samples() -> None:
    """A gene meeting the CPM threshold in enough samples must be retained
    even if it fails in a minority of samples."""
    counts = pd.DataFrame(
        {
            "case-a": [1000, 0],
            "case-b": [1000, 0],
            "case-c": [1000, 1000],
        },
        index=["always_high", "sometimes_high"],
    )

    filtered = filter_low_expression(counts, min_cpm=1.0, min_fraction_samples=0.3)

    assert "sometimes_high" in filtered.index


def test_normalize_and_vst_preserves_shape_and_labels_with_no_nans() -> None:
    """VST output must match input shape/labels and contain no NaNs."""
    rng = np.random.default_rng(42)
    n_genes, n_samples = 60, 8
    counts = pd.DataFrame(
        rng.poisson(lam=200, size=(n_genes, n_samples)),
        index=[f"gene{i}" for i in range(n_genes)],
        columns=[f"case-{i}" for i in range(n_samples)],
    )

    vst = normalize_and_vst(counts)

    assert vst.shape == counts.shape
    assert list(vst.index) == list(counts.index)
    assert list(vst.columns) == list(counts.columns)
    assert not vst.isna().to_numpy().any()


def test_normalize_and_vst_is_deterministic() -> None:
    """Given identical input, VST must produce identical output (no hidden randomness)."""
    rng = np.random.default_rng(7)
    n_genes, n_samples = 60, 8
    counts = pd.DataFrame(
        rng.poisson(lam=150, size=(n_genes, n_samples)),
        index=[f"gene{i}" for i in range(n_genes)],
        columns=[f"case-{i}" for i in range(n_samples)],
    )

    vst_first = normalize_and_vst(counts)
    vst_second = normalize_and_vst(counts)

    pd.testing.assert_frame_equal(vst_first, vst_second)
