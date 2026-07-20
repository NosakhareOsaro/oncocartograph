"""Tests for oncocartograph.preprocessing.methylation.

Beta value fixtures mirror the real SeSAMe level3betas.txt file format
(two-column, no header) confirmed against downloaded files.
"""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd
import pytest

from oncocartograph.preprocessing.methylation import (
    beta_to_m_value,
    build_methylation_matrix,
    filter_by_missingness,
    read_beta_values,
    select_top_variable_probes,
)


def test_read_beta_values_parses_headerless_two_column_file(tmp_path: Path) -> None:
    """Real SeSAMe files have no header row -- parsing must not eat the first data row."""
    path = tmp_path / "sample.methylation_array.sesame.level3betas.txt"
    path.write_text("cg00000029\t0.0671066940271816\ncg00000108\t0.976245089775103\n")

    beta = read_beta_values(path)

    assert len(beta) == 2
    assert beta["cg00000029"] == pytest.approx(0.0671066940271816)
    assert beta["cg00000108"] == pytest.approx(0.976245089775103)


def test_read_beta_values_treats_na_as_missing(tmp_path: Path) -> None:
    """SeSAMe's own masking marks unreliable probes as NA; must become NaN, not a string."""
    path = tmp_path / "sample.txt"
    path.write_text("cg00000029\tNA\ncg00000108\t0.5\n")

    beta = read_beta_values(path)

    assert math.isnan(beta["cg00000029"])
    assert beta["cg00000108"] == pytest.approx(0.5)


def test_beta_to_m_value_zero_point_five_maps_to_zero() -> None:
    """beta=0.5 is the logit midpoint: M = log2(0.5/0.5) = log2(1) = 0."""
    result = beta_to_m_value(pd.Series([0.5]))
    assert result.iloc[0] == pytest.approx(0.0)


def test_beta_to_m_value_is_finite_at_boundaries() -> None:
    """beta=0 and beta=1 must produce finite values (clipped), not +/-inf."""
    result = beta_to_m_value(pd.Series([0.0, 1.0]))
    assert math.isfinite(result.iloc[0])
    assert math.isfinite(result.iloc[1])
    assert result.iloc[0] < 0
    assert result.iloc[1] > 0


def test_beta_to_m_value_is_monotonic_increasing() -> None:
    """Higher beta must always produce a higher M-value."""
    result = beta_to_m_value(pd.Series([0.1, 0.3, 0.5, 0.7, 0.9]))
    assert list(result) == sorted(result)


def test_build_methylation_matrix_combines_patients(tmp_path: Path) -> None:
    """The matrix must be probe-indexed with one M-value column per patient."""
    path_a = tmp_path / "a.txt"
    path_a.write_text("cg00000029\t0.5\n")
    path_b = tmp_path / "b.txt"
    path_b.write_text("cg00000029\t0.9\n")

    matrix = build_methylation_matrix({"case-a": path_a, "case-b": path_b})

    assert list(matrix.columns) == ["case-a", "case-b"]
    assert matrix.loc["cg00000029", "case-a"] == pytest.approx(0.0)
    assert matrix.loc["cg00000029", "case-b"] == pytest.approx(
        beta_to_m_value(pd.Series([0.9])).iloc[0]
    )


def test_filter_by_missingness_drops_high_missingness_probes() -> None:
    """A probe missing in more than the allowed fraction of patients must be dropped."""
    matrix = pd.DataFrame(
        {
            "case-a": [1.0, float("nan"), float("nan")],
            "case-b": [1.0, 1.0, float("nan")],
            "case-c": [1.0, 1.0, float("nan")],
        },
        index=["clean", "one_missing", "all_missing"],
    )

    filtered = filter_by_missingness(matrix, max_missing_fraction=0.5)

    assert set(filtered.index) == {"clean", "one_missing"}


def test_select_top_variable_probes_keeps_highest_variance() -> None:
    """The n most variable probes (by variance across patients) must be kept, in variance order."""
    matrix = pd.DataFrame(
        {
            "case-a": [0.0, 0.0, 0.0],
            "case-b": [10.0, 1.0, 0.0],
            "case-c": [-10.0, -1.0, 0.0],
        },
        index=["high_var", "mid_var", "zero_var"],
    )

    top = select_top_variable_probes(matrix, n=2)

    assert list(top.index) == ["high_var", "mid_var"]


def test_select_top_variable_probes_keeps_all_when_n_exceeds_row_count() -> None:
    """Requesting more probes than exist must return all of them, not error."""
    matrix = pd.DataFrame({"case-a": [1.0, 2.0]}, index=["p1", "p2"])

    top = select_top_variable_probes(matrix, n=100)

    assert len(top) == 2
