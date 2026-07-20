"""Tests for oncocartograph.validation.burstein_check."""

from __future__ import annotations

from oncocartograph.scoring.survival import SurvivalEvidence
from oncocartograph.validation.burstein_check import (
    KNOWN_BIOLOGY_MARKERS,
    check_known_biology_markers,
)


def _evidence(hazard_ratio: float) -> SurvivalEvidence:
    return SurvivalEvidence(
        hazard_ratio=hazard_ratio,
        hazard_ratio_ci_low=hazard_ratio * 0.5,
        hazard_ratio_ci_high=hazard_ratio * 1.5,
        p_value=0.2,
        n_samples=100,
        n_events=20,
    )


def test_check_known_biology_markers_covers_all_markers_in_order() -> None:
    """Output must have one result per marker, in the declared order, even with no evidence."""
    results = check_known_biology_markers({})
    assert [r.gene_symbol for r in results] == [m.gene_symbol for m in KNOWN_BIOLOGY_MARKERS]
    assert all(r.plausible is None for r in results)


def test_check_known_biology_markers_protective_expectation_met() -> None:
    """AR (expected protective) with an observed HR<1 must be marked plausible."""
    results = check_known_biology_markers({"AR": _evidence(0.5)})
    ar_result = next(r for r in results if r.gene_symbol == "AR")
    assert ar_result.observed_direction == "protective"
    assert ar_result.plausible is True


def test_check_known_biology_markers_protective_expectation_violated() -> None:
    """AR (expected protective) with an observed HR>1 must be marked implausible, not hidden."""
    results = check_known_biology_markers({"AR": _evidence(1.8)})
    ar_result = next(r for r in results if r.gene_symbol == "AR")
    assert ar_result.observed_direction == "harmful"
    assert ar_result.plausible is False


def test_check_known_biology_markers_missing_gene_reports_unavailable() -> None:
    """A marker gene absent from the evidence dict must report None, not raise."""
    results = check_known_biology_markers({"PTEN": _evidence(1.2)})
    cd274_result = next(r for r in results if r.gene_symbol == "CD274")
    assert cd274_result.observed_hazard_ratio is None
    assert cd274_result.plausible is None
