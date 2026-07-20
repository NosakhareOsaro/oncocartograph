"""Tests for oncocartograph.scoring.evidence."""

from __future__ import annotations

import pytest

from oncocartograph.scoring.evidence import (
    BiomarkerEvidence,
    DruggabilityEvidence,
    IntegrationEvidence,
    RecurrenceEvidence,
)


def test_recurrence_evidence_mutation_fraction() -> None:
    """mutation_fraction must be n_patients_mutated / cohort_size."""
    evidence = RecurrenceEvidence(n_patients_mutated=12, cohort_size=122)

    assert evidence.mutation_fraction == pytest.approx(12 / 122)


def test_recurrence_evidence_association_defaults_to_none() -> None:
    """association must be optional, defaulting to None when no categorical test was run."""
    evidence = RecurrenceEvidence(n_patients_mutated=5, cohort_size=100)

    assert evidence.association is None


def test_biomarker_evidence_all_axes_optional() -> None:
    """A BiomarkerEvidence with only candidate_id must be constructible (evidence all optional)."""
    evidence = BiomarkerEvidence(candidate_id="TP53")

    assert evidence.survival is None
    assert evidence.integration is None
    assert evidence.recurrence is None
    assert evidence.druggability is None


def test_integration_evidence_is_frozen() -> None:
    """IntegrationEvidence must be immutable."""
    evidence = IntegrationEvidence(factor="Factor1", weight=0.5, view_variance_explained=10.0)
    with pytest.raises(AttributeError):
        evidence.weight = 1.0  # type: ignore[misc]


def test_druggability_evidence_chembl_max_phase_optional() -> None:
    """chembl_max_phase must default to None (e.g. no known compound for this target)."""
    evidence = DruggabilityEvidence(tractability_score=0.7)

    assert evidence.chembl_max_phase is None
