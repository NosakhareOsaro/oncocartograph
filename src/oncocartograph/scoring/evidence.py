"""Evidence data model: the distinct input types the composite score combines.

MOFA+-derived evidence (``IntegrationEvidence``) and mutation-recurrence
evidence (``RecurrenceEvidence``) are deliberately separate, explicitly
documented types feeding the same composite score -- not folded into a
single MOFA+-only pipeline. This matters because MOFA+ factor loadings do
not apply to mutation-derived candidates: in this project's actual MOFA+
run, the mutation view contributed essentially no variance to any factor
(<=0.003% everywhere, see ``docs/methods.md`` §3.4), so a scoring design
that only understood "MOFA+ evidence" would have no way to score mutation
candidates at all.
"""

from __future__ import annotations

from dataclasses import dataclass

from oncocartograph.scoring.association import AssociationEvidence
from oncocartograph.scoring.survival import SurvivalEvidence


@dataclass(frozen=True)
class IntegrationEvidence:
    """MOFA+ provenance for a candidate selected via factor loading.

    Attributes:
        factor: Which MOFA+ factor this candidate was selected from
            (e.g. ``"Factor1"``).
        weight: The feature's loading (weight) on that factor, as
            reported by ``mofax.get_weights`` -- retained for
            auditability, not itself used as the composite score's
            selection-pathway term (see ``composite.py``).
        view_variance_explained: The percentage (0-100, matching
            ``mofax.get_variance_explained``'s convention) of that
            view's variance the factor explains -- i.e. how
            well-supported the factor itself is.
    """

    factor: str
    weight: float
    view_variance_explained: float


@dataclass(frozen=True)
class RecurrenceEvidence:
    """Mutation-recurrence provenance for a candidate selected by recurrence filtering.

    Attributes:
        n_patients_mutated: Number of patients with a qualifying
            (non-synonymous) variant in this gene.
        cohort_size: Total number of patients with mutation data
            available -- the denominator for ``mutation_fraction``.
        association: An optional categorical association test result
            (e.g. vs. ``vital_status`` or a subtype label), if one was
            run for this candidate.
    """

    n_patients_mutated: int
    cohort_size: int
    association: AssociationEvidence | None = None

    @property
    def mutation_fraction(self) -> float:
        """Fraction of the mutation-data cohort carrying this mutation."""
        return self.n_patients_mutated / self.cohort_size


@dataclass(frozen=True)
class DruggabilityEvidence:
    """Druggability evidence schema.

    This package defines the shape only. Populating real values from
    Open Targets/ChEMBL is ``feat/drug-target-scoring``'s responsibility,
    not this package's -- until then, composite scoring simply treats
    this evidence axis as absent for every candidate.

    Attributes:
        tractability_score: A [0, 1] score reflecting how tractable this
            target is judged to be (source-defined upstream).
        chembl_max_phase: Highest clinical trial phase reached by any
            ChEMBL compound against this target (0-4), or ``None`` if no
            such compound is known.
    """

    tractability_score: float
    chembl_max_phase: float | None = None


@dataclass(frozen=True)
class BiomarkerEvidence:
    """All evidence collected for one candidate biomarker.

    Attributes:
        candidate_id: Gene or probe identifier.
        survival: Survival association evidence, present for any
            candidate whose Cox model could be fit.
        integration: MOFA+ provenance, present only if this candidate
            was selected via factor loading.
        recurrence: Mutation-recurrence provenance, present only if this
            candidate was selected via recurrence filtering.
        druggability: Druggability evidence, present only once
            ``feat/drug-target-scoring`` has populated it.
    """

    candidate_id: str
    survival: SurvivalEvidence | None = None
    integration: IntegrationEvidence | None = None
    recurrence: RecurrenceEvidence | None = None
    druggability: DruggabilityEvidence | None = None
