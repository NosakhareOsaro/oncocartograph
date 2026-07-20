"""Tests for oncocartograph.scoring.association."""

from __future__ import annotations

import pandas as pd
import pytest
from scipy.stats import fisher_exact

from oncocartograph.scoring.association import AssociationEvidence, fisher_exact_association


def test_fisher_exact_association_matches_scipy_reference() -> None:
    """Our wrapper's odds ratio/p-value must match scipy's own fisher_exact directly."""
    marker = pd.Series([1, 1, 1, 1, 0, 0, 0, 0, 0, 0])
    outcome = pd.Series(
        ["Dead", "Dead", "Dead", "Alive", "Alive", "Alive", "Alive", "Alive", "Alive", "Alive"]
    )

    result = fisher_exact_association(marker, outcome)

    assert result is not None
    # Rows sorted by marker (0, 1), columns sorted by outcome label ("Alive" < "Dead"):
    # marker=0 (6 patients): all Alive -> [6, 0]
    # marker=1 (4 patients): 1 Alive, 3 Dead -> [1, 3]
    expected_table = [[6, 0], [1, 3]]
    expected_odds_ratio, expected_p = fisher_exact(expected_table)
    assert result.odds_ratio == pytest.approx(expected_odds_ratio)
    assert result.p_value == pytest.approx(expected_p)
    assert result.contingency_table == ((6, 0), (1, 3))


def test_fisher_exact_association_contingency_table_matches_counts() -> None:
    """The retained contingency table must reflect the actual marker/outcome counts."""
    marker = pd.Series([1, 1, 0, 0])
    outcome = pd.Series(["A", "B", "A", "A"])

    result = fisher_exact_association(marker, outcome)

    assert result is not None
    # rows=marker (0,1), cols=outcome (A,B): marker=0->(A:2,B:0); marker=1->(A:1,B:1)
    assert result.contingency_table == ((2, 0), (1, 1))
    assert result.n_samples == 4


def test_fisher_exact_association_drops_missing_values() -> None:
    """Patients missing either value must be excluded before counting."""
    marker = pd.Series([1, 0, None, 1])
    outcome = pd.Series(["A", "A", "B", "B"])

    result = fisher_exact_association(marker, outcome)

    assert result is not None
    assert result.n_samples == 3


def test_fisher_exact_association_returns_none_for_single_category_marker() -> None:
    """A marker with only one distinct value cannot support the test."""
    marker = pd.Series([1, 1, 1, 1])
    outcome = pd.Series(["A", "A", "B", "B"])

    result = fisher_exact_association(marker, outcome)

    assert result is None


def test_fisher_exact_association_returns_none_for_non_binary_outcome() -> None:
    """An outcome with other than exactly two categories cannot support this test."""
    marker = pd.Series([1, 0, 1, 0])
    outcome = pd.Series(["A", "B", "C", "A"])

    result = fisher_exact_association(marker, outcome)

    assert result is None


def test_association_evidence_is_frozen() -> None:
    """AssociationEvidence must be immutable."""
    evidence = AssociationEvidence(
        odds_ratio=2.0, p_value=0.05, contingency_table=((1, 2), (3, 4)), n_samples=10
    )
    with pytest.raises(AttributeError):
        evidence.odds_ratio = 3.0  # type: ignore[misc]
