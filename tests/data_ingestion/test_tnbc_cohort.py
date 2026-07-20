"""Tests for oncocartograph.data_ingestion.tnbc_cohort.

All clinical values used here are synthetic fixtures constructed to
exercise the ADR 0001 classification rules -- they are not real TCGA
patient data.
"""

from __future__ import annotations

import math

import pandas as pd
import pytest

from oncocartograph.data_ingestion.tnbc_cohort import (
    EXCLUSION_REASON_COLUMN,
    IS_TNBC_COLUMN,
    CohortDecision,
    build_tnbc_cohort_audit,
    classify_patient,
    select_tnbc_cohort,
)


def test_all_negative_is_classified_as_tnbc() -> None:
    """ER-/PR-/HER2- (clean case) must be included with no exclusion reason."""
    decision = classify_patient("Negative", "Negative", "Negative", None)
    assert decision == CohortDecision(is_tnbc=True, exclusion_reason=None)


def test_er_positive_excludes_with_reason() -> None:
    """An ER-positive patient is not TNBC and the reason must name ER."""
    decision = classify_patient("Positive", "Negative", "Negative", None)
    assert decision.is_tnbc is False
    assert decision.exclusion_reason is not None
    assert "ER" in decision.exclusion_reason
    assert "Receptor positive" in decision.exclusion_reason


def test_pr_positive_excludes_with_reason() -> None:
    """A PR-positive patient is not TNBC and the reason must name PR."""
    decision = classify_patient("Negative", "Positive", "Negative", None)
    assert decision.is_tnbc is False
    assert decision.exclusion_reason is not None
    assert "PR" in decision.exclusion_reason


def test_her2_ihc_positive_excludes_with_reason() -> None:
    """A HER2 IHC-positive patient is not TNBC."""
    decision = classify_patient("Negative", "Negative", "Positive", None)
    assert decision.is_tnbc is False
    assert decision.exclusion_reason is not None
    assert "HER2" in decision.exclusion_reason


def test_her2_equivocal_with_fish_negative_is_included() -> None:
    """IHC 2+ resolved negative by FISH counts as HER2-negative (ADR 0001)."""
    decision = classify_patient("Negative", "Negative", "Equivocal", "Negative")
    assert decision == CohortDecision(is_tnbc=True, exclusion_reason=None)


def test_her2_equivocal_with_fish_positive_excludes_as_positive() -> None:
    """IHC 2+ resolved positive by FISH is a real HER2-positive result, not missing data."""
    decision = classify_patient("Negative", "Negative", "Equivocal", "Positive")
    assert decision.is_tnbc is False
    assert decision.exclusion_reason is not None
    assert "Receptor positive" in decision.exclusion_reason
    assert "HER2" in decision.exclusion_reason


def test_her2_equivocal_with_no_fish_is_indeterminate_not_negative() -> None:
    """The key ADR 0001 judgment call: equivocal-without-FISH must be excluded,
    not silently treated as negative."""
    decision = classify_patient("Negative", "Negative", "Equivocal", None)
    assert decision.is_tnbc is False
    assert decision.exclusion_reason is not None
    assert "no FISH follow-up" in decision.exclusion_reason


@pytest.mark.parametrize("missing_value", [None, math.nan, "", "[Not Evaluated]", "Indeterminate"])
def test_er_missing_or_indeterminate_excludes(missing_value: object) -> None:
    """Any non-Positive/Negative ER value must be treated as indeterminate, not negative."""
    decision = classify_patient(missing_value, "Negative", "Negative", None)
    assert decision.is_tnbc is False
    assert decision.exclusion_reason is not None
    assert "ER status missing or indeterminate" in decision.exclusion_reason


def test_her2_ihc_missing_is_indeterminate() -> None:
    """A HER2 IHC value that is neither Positive/Negative/Equivocal must be indeterminate."""
    decision = classify_patient("Negative", "Negative", None, None)
    assert decision.is_tnbc is False
    assert decision.exclusion_reason is not None
    assert "HER2 IHC status missing or indeterminate" in decision.exclusion_reason


def test_multiple_indeterminate_markers_combine_into_one_reason() -> None:
    """When both ER and HER2 are indeterminate, both must appear in the reason."""
    decision = classify_patient(None, "Negative", "Equivocal", None)
    assert decision.is_tnbc is False
    assert decision.exclusion_reason is not None
    assert "ER status missing or indeterminate" in decision.exclusion_reason
    assert "no FISH follow-up" in decision.exclusion_reason


def test_indeterminate_takes_precedence_over_positive() -> None:
    """If ER is indeterminate AND PR is positive, the reason must report the
    indeterminacy, not a receptor-positive verdict based on incomplete data."""
    decision = classify_patient(None, "Positive", "Negative", None)
    assert decision.is_tnbc is False
    assert decision.exclusion_reason is not None
    assert "indeterminate" in decision.exclusion_reason.lower()
    assert "Receptor positive" not in decision.exclusion_reason


def _synthetic_clinical_fixture() -> pd.DataFrame:
    """A small synthetic clinical table exercising every classification path."""
    return pd.DataFrame(
        [
            {
                "bcr_patient_barcode": "SYNTH-0001",
                "er_status_by_ihc": "Negative",
                "pr_status_by_ihc": "Negative",
                "her2_status_by_ihc": "Negative",
                "her2_fish_status": None,
            },
            {
                "bcr_patient_barcode": "SYNTH-0002",
                "er_status_by_ihc": "Positive",
                "pr_status_by_ihc": "Negative",
                "her2_status_by_ihc": "Negative",
                "her2_fish_status": None,
            },
            {
                "bcr_patient_barcode": "SYNTH-0003",
                "er_status_by_ihc": "Negative",
                "pr_status_by_ihc": "Negative",
                "her2_status_by_ihc": "Equivocal",
                "her2_fish_status": "Negative",
            },
            {
                "bcr_patient_barcode": "SYNTH-0004",
                "er_status_by_ihc": "Negative",
                "pr_status_by_ihc": "Negative",
                "her2_status_by_ihc": "Equivocal",
                "her2_fish_status": None,
            },
        ]
    )


def test_build_tnbc_cohort_audit_adds_expected_columns() -> None:
    """The audit table must contain one is_tnbc/exclusion_reason pair per input row."""
    clinical = _synthetic_clinical_fixture()

    audit = build_tnbc_cohort_audit(clinical)

    assert len(audit) == len(clinical)
    assert list(audit[IS_TNBC_COLUMN]) == [True, False, True, False]
    assert (
        audit.loc[audit["bcr_patient_barcode"] == "SYNTH-0002", EXCLUSION_REASON_COLUMN].iloc[0]
        is not None
    )


def test_build_tnbc_cohort_audit_raises_on_missing_column() -> None:
    """A clinical DataFrame missing a required column must fail loudly, not silently."""
    clinical = _synthetic_clinical_fixture().drop(columns=["her2_fish_status"])

    with pytest.raises(KeyError):
        build_tnbc_cohort_audit(clinical)


def test_select_tnbc_cohort_filters_and_drops_audit_columns() -> None:
    """select_tnbc_cohort must return only included patients, without audit columns."""
    audit = build_tnbc_cohort_audit(_synthetic_clinical_fixture())

    cohort = select_tnbc_cohort(audit)

    assert set(cohort["bcr_patient_barcode"]) == {"SYNTH-0001", "SYNTH-0003"}
    assert IS_TNBC_COLUMN not in cohort.columns
    assert EXCLUSION_REASON_COLUMN not in cohort.columns
