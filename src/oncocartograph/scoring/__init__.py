"""Composite biomarker scoring: survival validation + druggability prioritisation.

This is the standalone, independently versioned and unit-tested contribution
of OncoCartograph: a scoring function that combines statistical survival
association evidence with druggability evidence (Open Targets / ChEMBL) into
a single, reproducible composite biomarker priority score. See
``src/oncocartograph/scoring/README.md`` for package-level documentation.

This package has zero dependency on the rest of ``oncocartograph`` --
every module here operates on plain pandas Series/DataFrames and
dataclasses, never importing ``oncocartograph.data_ingestion``,
``oncocartograph.preprocessing``, or ``oncocartograph.integration``. This
is deliberate: it is meant to be extractable to a standalone,
independently publishable package with no more than an import-path
change.
"""

from oncocartograph.scoring.association import AssociationEvidence, fisher_exact_association
from oncocartograph.scoring.composite import ScoreWeights, composite_biomarker_score
from oncocartograph.scoring.evidence import (
    BiomarkerEvidence,
    DruggabilityEvidence,
    IntegrationEvidence,
    RecurrenceEvidence,
)
from oncocartograph.scoring.survival import (
    SurvivalEvidence,
    fit_univariate_cox,
    screen_survival_associations,
)

__all__ = [
    "AssociationEvidence",
    "BiomarkerEvidence",
    "DruggabilityEvidence",
    "IntegrationEvidence",
    "RecurrenceEvidence",
    "ScoreWeights",
    "SurvivalEvidence",
    "composite_biomarker_score",
    "fisher_exact_association",
    "fit_univariate_cox",
    "screen_survival_associations",
]
