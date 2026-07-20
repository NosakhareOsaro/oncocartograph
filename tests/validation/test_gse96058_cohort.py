"""Tests for oncocartograph.validation.gse96058_cohort.

Values are synthetic, shaped after the real GSE96058 er/pgr/her2 status
encoding (0=negative, 1=positive, "NA"=missing) confirmed against the
real downloaded data.
"""

from __future__ import annotations

import pandas as pd
import pytest

from oncocartograph.validation.gse96058_cohort import (
    EXCLUSION_REASON_COLUMN,
    IS_TNBC_COLUMN,
    Gse96058CohortDecision,
    build_gse96058_cohort_audit,
    classify_gse96058_sample,
    select_gse96058_tnbc_cohort,
)


def test_classify_all_negative_is_tnbc() -> None:
    """ER-/PgR-/HER2- must be included with no exclusion reason."""
    decision = classify_gse96058_sample("0", "0", "0")
    assert decision == Gse96058CohortDecision(is_tnbc=True, exclusion_reason=None)


def test_classify_er_positive_excludes() -> None:
    """An ER-positive sample must be excluded, reason naming ER."""
    decision = classify_gse96058_sample("1", "0", "0")
    assert decision.is_tnbc is False
    assert decision.exclusion_reason is not None
    assert "ER" in decision.exclusion_reason
    assert "Receptor positive" in decision.exclusion_reason


def test_classify_missing_value_excludes_as_indeterminate() -> None:
    """ "NA" must be treated as indeterminate, not negative."""
    decision = classify_gse96058_sample("NA", "0", "0")
    assert decision.is_tnbc is False
    assert decision.exclusion_reason is not None
    assert "ER status missing or indeterminate" in decision.exclusion_reason


def test_classify_indeterminate_takes_precedence_over_positive() -> None:
    """A sample with one indeterminate and one positive marker must report
    the indeterminacy, not a receptor-positive verdict from incomplete data."""
    decision = classify_gse96058_sample("NA", "1", "0")
    assert decision.is_tnbc is False
    assert "indeterminate" in decision.exclusion_reason.lower()
    assert "Receptor positive" not in decision.exclusion_reason


def _synthetic_clinical() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"title": "F1", "er status": "0", "pgr status": "0", "her2 status": "0"},
            {"title": "F1repl", "er status": "0", "pgr status": "0", "her2 status": "0"},
            {"title": "F2", "er status": "1", "pgr status": "0", "her2 status": "0"},
            {"title": "F3", "er status": "0", "pgr status": "0", "her2 status": "NA"},
        ]
    )


def test_build_gse96058_cohort_audit_drops_replicates_before_classifying() -> None:
    """A technical replicate (title ending in 'repl') must be excluded before
    classification, not double-counted alongside its primary sample."""
    audit = build_gse96058_cohort_audit(_synthetic_clinical())

    assert "F1repl" not in set(audit["title"])
    assert len(audit) == 3


def test_build_gse96058_cohort_audit_classifies_correctly() -> None:
    """Audit output must match expected per-sample classification."""
    audit = build_gse96058_cohort_audit(_synthetic_clinical())

    by_title = audit.set_index("title")
    assert by_title.loc["F1", IS_TNBC_COLUMN]
    assert not by_title.loc["F2", IS_TNBC_COLUMN]
    assert not by_title.loc["F3", IS_TNBC_COLUMN]
    assert by_title.loc["F2", EXCLUSION_REASON_COLUMN] is not None


def test_build_gse96058_cohort_audit_raises_on_missing_column() -> None:
    """A clinical DataFrame missing a required column must fail loudly."""
    clinical = _synthetic_clinical().drop(columns=["her2 status"])

    with pytest.raises(KeyError):
        build_gse96058_cohort_audit(clinical)


def test_select_gse96058_tnbc_cohort_filters_and_drops_audit_columns() -> None:
    """select_gse96058_tnbc_cohort must return only included samples, without audit columns."""
    audit = build_gse96058_cohort_audit(_synthetic_clinical())

    cohort = select_gse96058_tnbc_cohort(audit)

    assert set(cohort["title"]) == {"F1"}
    assert IS_TNBC_COLUMN not in cohort.columns
    assert EXCLUSION_REASON_COLUMN not in cohort.columns
