"""Composite biomarker priority score.

Combines survival, druggability, and selection-pathway (MOFA+
integration or mutation recurrence) evidence into one ranked score per
candidate. Missing evidence axes -- e.g. a mutation candidate has no
``IntegrationEvidence`` at all, since MOFA+ loadings do not apply to it
(see ``docs/methods.md`` §3.4) -- are excluded from the weighted average
and the remaining weights renormalized, rather than penalizing a
candidate for lacking an evidence type that structurally cannot apply
to it.
"""

from __future__ import annotations

from dataclasses import dataclass

from oncocartograph.scoring.evidence import BiomarkerEvidence

#: Upper bound (percent) mofax's variance_explained convention uses;
#: divide by this to normalise IntegrationEvidence into a [0, 1] score.
_VARIANCE_EXPLAINED_MAX_PERCENT = 100.0


@dataclass(frozen=True)
class ScoreWeights:
    """Composite score weights across evidence axes.

    All weights are explicit, named, and overridable per call -- never
    hidden constants. They are renormalized over whichever axes are
    actually present for a given candidate (see
    :func:`composite_biomarker_score`), so the values here represent
    *relative* importance, not an assumption that all three axes are
    always available.

    Attributes:
        survival: Weight for survival association evidence -- the core
            statistical claim.
        druggability: Weight for druggability evidence (actionability).
        selection_pathway: Weight for the selection-pathway bonus: MOFA+
            integration confidence, or mutation recurrence fraction,
            whichever applies to a given candidate.
    """

    survival: float = 0.5
    druggability: float = 0.35
    selection_pathway: float = 0.15


def _survival_score(evidence: BiomarkerEvidence) -> float | None:
    """Score in [0, 1]: rewards harmful (HR>1), significant associations only.

    Protective associations (hazard_ratio <= 1) score 0 rather than a
    negative value -- this project's explicit choice to treat "druggable
    biomarker" as "something harmful to disrupt," not "anything
    survival-associated regardless of direction." Significance uses
    ``1 - p`` (preferring the FDR-adjusted p-value when available) as a
    simple, bounded transform requiring no arbitrary cap constant.
    """
    survival = evidence.survival
    if survival is None:
        return None
    if survival.hazard_ratio <= 1.0:
        return 0.0
    p = survival.p_value_adjusted if survival.p_value_adjusted is not None else survival.p_value
    return 1.0 - p


def _selection_pathway_score(evidence: BiomarkerEvidence) -> float | None:
    """Score in [0, 1] from whichever selection-pathway evidence is present.

    For MOFA+-derived candidates, uses the supporting factor's
    view-variance-explained (how well-supported the factor itself is),
    not the candidate's raw loading weight, which has no natural [0, 1]
    bound. For mutation-derived candidates, uses the mutation's
    prevalence fraction, which is naturally in [0, 1] already.
    """
    if evidence.integration is not None:
        return min(
            evidence.integration.view_variance_explained / _VARIANCE_EXPLAINED_MAX_PERCENT, 1.0
        )
    if evidence.recurrence is not None:
        return evidence.recurrence.mutation_fraction
    return None


def _druggability_score(evidence: BiomarkerEvidence) -> float | None:
    if evidence.druggability is None:
        return None
    return evidence.druggability.tractability_score


def composite_biomarker_score(
    evidence: BiomarkerEvidence, weights: ScoreWeights | None = None
) -> float:
    """Combine available evidence axes into one composite priority score.

    Args:
        evidence: All evidence collected for one candidate.
        weights: Per-axis weights; renormalized over whichever axes are
            actually present for this candidate. Defaults to
            :class:`ScoreWeights`'s defaults if not supplied.

    Returns:
        A composite score in [0, 1]. Higher means higher priority.

    Raises:
        ValueError: If no evidence axis at all is available for this
            candidate (nothing to score).
    """
    weights = weights if weights is not None else ScoreWeights()
    axis_scores = {
        "survival": (_survival_score(evidence), weights.survival),
        "druggability": (_druggability_score(evidence), weights.druggability),
        "selection_pathway": (_selection_pathway_score(evidence), weights.selection_pathway),
    }
    available = {
        name: (score, weight) for name, (score, weight) in axis_scores.items() if score is not None
    }
    if not available:
        raise ValueError(f"No evidence available at all for candidate {evidence.candidate_id!r}")

    total_weight = sum(weight for _, weight in available.values())
    return sum(score * weight for score, weight in available.values()) / total_weight
