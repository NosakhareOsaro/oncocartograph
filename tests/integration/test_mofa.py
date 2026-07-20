"""Tests for oncocartograph.integration.mofa.

Uses small synthetic multi-omic matrices and a fast/low-iteration MOFA+
training configuration so the end-to-end training test runs quickly --
this is a real (if tiny) MOFA+ training run, not a mock.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from oncocartograph.integration.mofa import (
    build_mofa_input,
    build_view_long_df,
    get_factor_values,
    get_variance_explained,
    load_mofa_model,
    train_mofa_model,
)


def test_build_view_long_df_melts_and_labels_view() -> None:
    """The long df must have one row per (feature, sample) pair, labelled with the view."""
    matrix = pd.DataFrame({"case-a": [1.0, 2.0], "case-b": [3.0, 4.0]}, index=["gene1", "gene2"])

    long_df = build_view_long_df("rna_seq", matrix)

    assert set(long_df.columns) == {"feature", "sample", "value", "view", "group"}
    assert len(long_df) == 4
    assert (long_df["view"] == "rna_seq").all()
    assert (long_df["group"] == "all").all()


def test_build_view_long_df_drops_nan_rows() -> None:
    """A missing value must produce no row at all, not a NaN-valued row."""
    matrix = pd.DataFrame({"case-a": [1.0, float("nan")]}, index=["gene1", "gene2"])

    long_df = build_view_long_df("rna_seq", matrix)

    assert len(long_df) == 1
    assert long_df.iloc[0]["feature"] == "gene1"


def test_build_mofa_input_combines_views_with_partial_overlap() -> None:
    """Views with different patient sets must combine without error, samples missing
    from a view simply absent from that view's rows."""
    view_a = pd.DataFrame({"case-1": [1.0], "case-2": [2.0]}, index=["f1"])
    view_b = pd.DataFrame({"case-1": [3.0]}, index=["f1"])

    combined = build_mofa_input({"view_a": view_a, "view_b": view_b})

    assert set(combined["view"]) == {"view_a", "view_b"}
    assert set(combined.loc[combined["view"] == "view_b", "sample"]) == {"case-1"}


def test_train_mofa_model_raises_on_mismatched_view_and_likelihood_keys(tmp_path: Path) -> None:
    """A likelihoods dict not matching views' keys must fail loudly before touching mofapy2."""
    views = {"rna_seq": pd.DataFrame({"case-1": [1.0]}, index=["f1"])}
    likelihoods = {"wrong_name": "gaussian"}

    with pytest.raises(ValueError, match="same keys"):
        train_mofa_model(
            views,
            likelihoods,  # type: ignore[arg-type]
            n_factors=2,
            seed=42,
            outfile=tmp_path / "model.hdf5",
        )


@pytest.fixture(scope="module")
def small_trained_model_path(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Train a small real MOFA+ model once and reuse it across result-extraction tests."""
    rng = np.random.default_rng(0)
    n_samples, n_features = 20, 15
    samples = [f"case-{i}" for i in range(n_samples)]

    gaussian_view = pd.DataFrame(
        rng.normal(size=(n_features, n_samples)),
        index=[f"gene{i}" for i in range(n_features)],
        columns=samples,
    )
    bernoulli_view = pd.DataFrame(
        rng.integers(0, 2, size=(n_features, n_samples)),
        index=[f"mut{i}" for i in range(n_features)],
        columns=samples,
    ).astype(float)

    outfile = tmp_path_factory.mktemp("mofa") / "test_model.hdf5"
    train_mofa_model(
        {"gaussian_view": gaussian_view, "bernoulli_view": bernoulli_view},
        {"gaussian_view": "gaussian", "bernoulli_view": "bernoulli"},
        n_factors=3,
        seed=42,
        outfile=outfile,
        convergence_mode="fast",
        max_iterations=20,
    )
    return outfile


def test_train_mofa_model_creates_output_file(small_trained_model_path: Path) -> None:
    """Training must produce a readable HDF5 file at the requested path."""
    assert small_trained_model_path.exists()


def test_get_factor_values_returns_sample_by_factor_matrix(small_trained_model_path: Path) -> None:
    """Factor values must be indexed by sample, one column per factor."""
    model = load_mofa_model(small_trained_model_path)

    factors = get_factor_values(model)

    assert factors.shape == (20, 3)


def test_get_variance_explained_covers_every_view_and_factor_combination(
    small_trained_model_path: Path,
) -> None:
    """Variance explained must have one row per (factor, view) pair, both views present."""
    model = load_mofa_model(small_trained_model_path)

    r2 = get_variance_explained(model)

    assert set(r2["View"]) == {"gaussian_view", "bernoulli_view"}
    assert r2["Factor"].nunique() == 3
    assert len(r2) == 3 * 2
    assert (r2["R2"] >= 0).all()
