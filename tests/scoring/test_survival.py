"""Tests for oncocartograph.scoring.survival.

Synthetic survival data is generated from a known planted hazard
relationship (covariate -> exponential hazard), so tests check that
fitting recovers the correct *direction* and *significance* of a known
effect, not exact coefficient values (which have real sampling
variability even at n=200) -- the appropriate precision level for a
statistical-fitting test.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from oncocartograph.scoring.survival import (
    SurvivalEvidence,
    fit_univariate_cox,
    screen_survival_associations,
)


def _simulate_survival_data(
    n: int, true_log_hr: float, seed: int
) -> tuple[pd.Series, pd.Series, pd.Series]:
    """Simulate (covariate, duration, event) with a known planted hazard ratio."""
    rng = np.random.default_rng(seed)
    index = [f"patient-{i}" for i in range(n)]
    covariate = pd.Series(rng.normal(size=n), index=index)
    baseline_hazard = 0.01
    hazard = baseline_hazard * np.exp(true_log_hr * covariate.to_numpy())
    true_duration = rng.exponential(1 / hazard)
    censoring_time = rng.exponential(200, size=n)
    duration = pd.Series(np.minimum(true_duration, censoring_time), index=index)
    event = pd.Series((true_duration <= censoring_time).astype(int), index=index)
    return covariate, duration, event


def test_fit_univariate_cox_recovers_harmful_direction() -> None:
    """A covariate planted with a positive log-hazard must fit HR > 1 and be significant."""
    covariate, duration, event = _simulate_survival_data(n=200, true_log_hr=1.0, seed=1)

    result = fit_univariate_cox(covariate, duration, event)

    assert result is not None
    assert result.hazard_ratio > 1.0
    assert result.p_value < 0.05
    assert result.hazard_ratio_ci_low < result.hazard_ratio < result.hazard_ratio_ci_high


def test_fit_univariate_cox_recovers_protective_direction() -> None:
    """A covariate planted with a negative log-hazard must fit HR < 1 and be significant."""
    covariate, duration, event = _simulate_survival_data(n=200, true_log_hr=-1.0, seed=2)

    result = fit_univariate_cox(covariate, duration, event)

    assert result is not None
    assert result.hazard_ratio < 1.0
    assert result.p_value < 0.05


def test_fit_univariate_cox_reports_sample_and_event_counts() -> None:
    """n_samples/n_events must reflect the data actually used, after dropping NaNs."""
    covariate, duration, event = _simulate_survival_data(n=50, true_log_hr=0.5, seed=3)

    result = fit_univariate_cox(covariate, duration, event)

    assert result is not None
    assert result.n_samples == 50
    assert result.n_events == int(event.sum())


def test_fit_univariate_cox_drops_missing_values() -> None:
    """Patients missing a biomarker value must be excluded, not crash the fit."""
    covariate, duration, event = _simulate_survival_data(n=100, true_log_hr=1.0, seed=4)
    covariate_with_gaps = covariate.copy()
    covariate_with_gaps.iloc[:10] = float("nan")

    result = fit_univariate_cox(covariate_with_gaps, duration, event)

    assert result is not None
    assert result.n_samples == 90


def test_fit_univariate_cox_returns_none_for_zero_variance_biomarker() -> None:
    """A biomarker with the same value in every patient cannot be fit and must not raise."""
    n = 50
    index = [f"patient-{i}" for i in range(n)]
    constant_values = pd.Series([1.0] * n, index=index)
    duration = pd.Series(np.random.default_rng(5).exponential(100, n), index=index)
    event = pd.Series([1] * n, index=index)

    result = fit_univariate_cox(constant_values, duration, event)

    assert result is None


def test_fit_univariate_cox_returns_none_for_degenerate_sparse_binary_covariate() -> None:
    """A rare binary marker with zero events in its mutated subgroup must return None,
    not a SurvivalEvidence with NaN statistics.

    Reproduces a real bug found running against the actual TCGA-BRCA TNBC
    mutation data: lifelines does not raise on this degenerate case, it
    returns a "successful" fit with NaN summary statistics, which later
    broke FDR correction (scipy raises on out-of-range p-values). ~11% of
    real recurrence-filtered mutation genes hit this in practice.
    """
    n = 122
    # Marker is 1 for exactly 3 patients, all of whom are censored (event=0)
    # -- zero events in the mutated subgroup, the actual real-world pattern.
    marker = pd.Series([1] * 3 + [0] * (n - 3))
    duration = pd.Series(np.random.default_rng(10).exponential(500, n))
    event = pd.Series([0, 0, 0] + [1] * 10 + [0] * (n - 13))

    result = fit_univariate_cox(marker, duration, event)

    assert result is None


def test_fit_univariate_cox_returns_none_for_too_few_samples() -> None:
    """A single-patient input cannot support a Cox fit."""
    result = fit_univariate_cox(
        pd.Series([1.0], index=["p1"]),
        pd.Series([100.0], index=["p1"]),
        pd.Series([1], index=["p1"]),
    )

    assert result is None


def test_screen_survival_associations_screens_multiple_features() -> None:
    """Screening must return one row per fittable feature, with an FDR-adjusted p-value column.

    duration/event are simulated from a single hidden covariate, so that
    "harmful_gene" (the hidden covariate itself) and "protective_gene"
    (its negation) have a genuine, known relationship to survival --
    not independently-simulated data compared against unrelated
    survival times.
    """
    hidden_covariate, duration, event = _simulate_survival_data(n=150, true_log_hr=1.2, seed=6)
    rng = np.random.default_rng(7)

    harmful = hidden_covariate
    protective = -hidden_covariate
    noise = pd.Series(rng.normal(size=150), index=duration.index)
    constant = pd.Series([1.0] * 150, index=duration.index)
    # Degenerate sparse binary marker (zero events in its mutated subgroup)
    # must not reach FDR correction with a NaN p-value and crash the batch.
    # Deterministically pick 3 censored (event=0) patients as "mutated" so
    # the mutated subgroup is guaranteed to have zero observed events.
    censored_patients = event[event == 0].index[:3]
    degenerate = pd.Series(0, index=duration.index)
    degenerate.loc[censored_patients] = 1

    matrix = pd.DataFrame(
        [harmful, protective, noise, constant, degenerate],
        index=["harmful_gene", "protective_gene", "noise_gene", "constant_gene", "degenerate_gene"],
    )

    results = screen_survival_associations(matrix, duration, event)

    assert set(results.index) == {"harmful_gene", "protective_gene", "noise_gene"}
    assert "constant_gene" not in results.index
    assert "degenerate_gene" not in results.index
    assert "p_adj" in results.columns
    assert results.loc["harmful_gene", "hazard_ratio"] > 1.0
    assert results.loc["protective_gene", "hazard_ratio"] < 1.0
    assert (results["p_adj"] >= results["p_value"]).all()
    assert results["p_value"].between(0, 1).all()


def test_screen_survival_associations_empty_matrix_returns_empty_frame() -> None:
    """An empty feature matrix must return an empty DataFrame, not error."""
    duration = pd.Series([100.0, 200.0], index=["p1", "p2"])
    event = pd.Series([1, 0], index=["p1", "p2"])
    matrix = pd.DataFrame(index=[], columns=["p1", "p2"])

    results = screen_survival_associations(matrix, duration, event)

    assert results.empty


def test_survival_evidence_is_frozen() -> None:
    """SurvivalEvidence must be immutable (frozen dataclass)."""
    evidence = SurvivalEvidence(
        hazard_ratio=1.5,
        hazard_ratio_ci_low=1.1,
        hazard_ratio_ci_high=2.0,
        p_value=0.01,
        n_samples=100,
        n_events=20,
    )
    with pytest.raises(AttributeError):
        evidence.hazard_ratio = 2.0  # type: ignore[misc]
