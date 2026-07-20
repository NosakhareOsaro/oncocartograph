"""Tests for oncocartograph.scoring.composite."""

from __future__ import annotations

import pytest

from oncocartograph.scoring.composite import ScoreWeights, composite_biomarker_score
from oncocartograph.scoring.evidence import (
    BiomarkerEvidence,
    DruggabilityEvidence,
    IntegrationEvidence,
    RecurrenceEvidence,
)
from oncocartograph.scoring.survival import SurvivalEvidence


def _survival(
    hazard_ratio: float, p_value: float, p_value_adjusted: float | None = None
) -> SurvivalEvidence:
    return SurvivalEvidence(
        hazard_ratio=hazard_ratio,
        hazard_ratio_ci_low=hazard_ratio * 0.8,
        hazard_ratio_ci_high=hazard_ratio * 1.2,
        p_value=p_value,
        p_value_adjusted=p_value_adjusted,
        n_samples=100,
        n_events=16,
    )


def test_composite_score_with_all_axes_present() -> None:
    """All three axes present: weighted average using the default weights."""
    evidence = BiomarkerEvidence(
        candidate_id="GENE1",
        survival=_survival(hazard_ratio=2.0, p_value=0.04, p_value_adjusted=0.1),
        druggability=DruggabilityEvidence(tractability_score=0.8),
        integration=IntegrationEvidence(factor="Factor2", weight=1.5, view_variance_explained=50.0),
    )

    score = composite_biomarker_score(evidence)

    # survival=1-0.1=0.9; druggability=0.8; selection_pathway=50/100=0.5
    # composite = 0.9*0.5 + 0.8*0.35 + 0.5*0.15 = 0.45 + 0.28 + 0.075 = 0.805
    assert score == pytest.approx(0.805)


def test_composite_score_renormalizes_when_druggability_missing() -> None:
    """Missing druggability must renormalize weights over the remaining two axes, not score as 0."""
    evidence = BiomarkerEvidence(
        candidate_id="GENE2",
        survival=_survival(hazard_ratio=2.0, p_value=0.04, p_value_adjusted=0.1),
        integration=IntegrationEvidence(factor="Factor2", weight=1.5, view_variance_explained=50.0),
    )

    score = composite_biomarker_score(evidence)

    # available weights: survival=0.5, selection_pathway=0.15, total=0.65
    # composite = (0.9*0.5 + 0.5*0.15) / 0.65
    expected = (0.9 * 0.5 + 0.5 * 0.15) / 0.65
    assert score == pytest.approx(expected)


def test_composite_score_floors_protective_survival_to_zero() -> None:
    """A protective association (HR<=1) must score 0 on the survival axis, not negative."""
    evidence = BiomarkerEvidence(
        candidate_id="GENE3",
        survival=_survival(hazard_ratio=0.5, p_value=0.001, p_value_adjusted=0.01),
        druggability=DruggabilityEvidence(tractability_score=0.8),
        integration=IntegrationEvidence(factor="Factor2", weight=1.5, view_variance_explained=50.0),
    )

    score = composite_biomarker_score(evidence)

    # survival_score = 0 (floored); composite = (0*0.5 + 0.8*0.35 + 0.5*0.15) / 1.0
    assert score == pytest.approx(0.8 * 0.35 + 0.5 * 0.15)


def test_composite_score_hazard_ratio_exactly_one_is_not_harmful() -> None:
    """HR == 1.0 (no effect) must also floor to 0, not be treated as harmful."""
    evidence = BiomarkerEvidence(
        candidate_id="GENE_NULL",
        survival=_survival(hazard_ratio=1.0, p_value=0.5, p_value_adjusted=0.9),
    )

    score = composite_biomarker_score(evidence)

    assert score == pytest.approx(0.0)


def test_composite_score_for_mutation_candidate_uses_recurrence_not_integration() -> None:
    """A mutation-derived candidate (recurrence, no integration) must still score correctly."""
    evidence = BiomarkerEvidence(
        candidate_id="TP53",
        survival=_survival(hazard_ratio=3.0, p_value=0.01, p_value_adjusted=0.05),
        recurrence=RecurrenceEvidence(n_patients_mutated=10, cohort_size=100),
    )

    score = composite_biomarker_score(evidence)

    # survival_score = 1 - 0.05 = 0.95; selection_pathway_score = 10/100 = 0.1
    # available weights: survival=0.5, selection_pathway=0.15, total=0.65
    expected = (0.95 * 0.5 + 0.1 * 0.15) / 0.65
    assert score == pytest.approx(expected)


def test_composite_score_prefers_adjusted_p_value_over_raw() -> None:
    """When p_value_adjusted is set, it must be used instead of the raw p_value."""
    evidence_with_adj = BiomarkerEvidence(
        candidate_id="G1", survival=_survival(hazard_ratio=2.0, p_value=0.001, p_value_adjusted=0.5)
    )
    evidence_without_adj = BiomarkerEvidence(
        candidate_id="G2",
        survival=_survival(hazard_ratio=2.0, p_value=0.001, p_value_adjusted=None),
    )

    score_with_adj = composite_biomarker_score(evidence_with_adj)
    score_without_adj = composite_biomarker_score(evidence_without_adj)

    assert score_with_adj == pytest.approx(1 - 0.5)
    assert score_without_adj == pytest.approx(1 - 0.001)
    assert score_with_adj < score_without_adj


def test_composite_score_respects_custom_weights() -> None:
    """Custom ScoreWeights must be used verbatim (before renormalization), not the defaults."""
    evidence = BiomarkerEvidence(
        candidate_id="GENE4",
        survival=_survival(hazard_ratio=2.0, p_value=0.04, p_value_adjusted=0.1),
        druggability=DruggabilityEvidence(tractability_score=0.8),
    )
    weights = ScoreWeights(survival=1.0, druggability=1.0, selection_pathway=0.0)

    score = composite_biomarker_score(evidence, weights=weights)

    # both weights equal 1.0 -> simple average of the two available axis scores
    assert score == pytest.approx((0.9 + 0.8) / 2)


def test_composite_score_raises_when_no_evidence_at_all() -> None:
    """A candidate with zero evidence axes must fail loudly, not silently return 0."""
    evidence = BiomarkerEvidence(candidate_id="EMPTY")

    with pytest.raises(ValueError, match="EMPTY"):
        composite_biomarker_score(evidence)
