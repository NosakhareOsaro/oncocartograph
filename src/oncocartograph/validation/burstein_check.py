"""Lightweight known-biology plausibility check (Burstein et al. 2015).

Scope was explicitly reduced, with the project owner's sign-off, from
reproducing Burstein et al.'s full LAR/MES/BLIS/BLIA transcriptomic
subtyping to a much smaller descriptive check: whether a handful of
genes with well-documented TNBC biology behave in the direction the
literature would suggest, using real GSE96058 survival evidence. This
is not a pass/fail gate on the pipeline as a whole -- it is reported as
a qualitative, per-gene plausibility note alongside the quantitative
direction-concordance result in :mod:`oncocartograph.validation.replication`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from oncocartograph.scoring.survival import SurvivalEvidence

Direction = Literal["harmful", "protective"]


@dataclass(frozen=True)
class KnownBiologyMarker:
    """A gene whose expected survival-association direction is documented in the literature.

    Attributes:
        gene_symbol: HGNC gene symbol.
        rationale: Why this gene and this expected direction were chosen.
        expected_direction: ``"harmful"`` (HR > 1 expected) or
            ``"protective"`` (HR < 1 expected).
    """

    gene_symbol: str
    rationale: str
    expected_direction: Direction


#: Genes and their literature-expected direction. Immune-checkpoint genes
#: (CD274/PD-L1, PDCD1/PD-1, CTLA4) are included as markers of
#: tumor-infiltrating lymphocytes, not as therapeutic targets here -- TIL
#: infiltration is repeatedly reported as favorable-prognosis in TNBC.
KNOWN_BIOLOGY_MARKERS: tuple[KnownBiologyMarker, ...] = (
    KnownBiologyMarker(
        gene_symbol="AR",
        rationale=(
            "Defines Burstein et al.'s luminal androgen receptor (LAR) subtype; "
            "AR-driven TNBC is repeatedly reported as a less proliferative, "
            "relatively less aggressive phenotype."
        ),
        expected_direction="protective",
    ),
    KnownBiologyMarker(
        gene_symbol="PTEN",
        rationale=(
            "Tumor suppressor with recurrent loss across TNBC subtypes "
            "(Burstein et al. 2015; Shah et al. 2012); loss of a tumor "
            "suppressor is expected to associate with worse outcome, i.e. "
            "lower expression with higher hazard -- so higher PTEN expression "
            "is expected to be protective."
        ),
        expected_direction="protective",
    ),
    KnownBiologyMarker(
        gene_symbol="CD274",
        rationale="PD-L1; marker of tumor-infiltrating lymphocytes, repeatedly favorable in TNBC.",
        expected_direction="protective",
    ),
    KnownBiologyMarker(
        gene_symbol="PDCD1",
        rationale="PD-1; marker of tumor-infiltrating lymphocytes, repeatedly favorable in TNBC.",
        expected_direction="protective",
    ),
    KnownBiologyMarker(
        gene_symbol="CTLA4",
        rationale="CTLA-4; marker of tumor-infiltrating lymphocytes, repeatedly favorable in TNBC.",
        expected_direction="protective",
    ),
)


@dataclass(frozen=True)
class KnownBiologyCheckResult:
    """The plausibility outcome for one known-biology marker.

    Attributes:
        gene_symbol: HGNC gene symbol.
        expected_direction: The literature-expected direction.
        observed_hazard_ratio: The gene's real GSE96058 hazard ratio, or
            ``None`` if it could not be fit.
        observed_direction: ``"harmful"``/``"protective"`` derived from
            the observed hazard ratio, or ``None`` if unavailable.
        plausible: Whether observed and expected direction agree, or
            ``None`` if unavailable. This is descriptive, not a
            statistical test -- no significance claim is attached.
    """

    gene_symbol: str
    expected_direction: Direction
    observed_hazard_ratio: float | None
    observed_direction: Direction | None
    plausible: bool | None


def _observed_direction(hazard_ratio: float) -> Direction:
    return "harmful" if hazard_ratio > 1 else "protective"


def check_known_biology_markers(
    evidence: dict[str, SurvivalEvidence | None],
) -> list[KnownBiologyCheckResult]:
    """Compare each known-biology marker's real GSE96058 direction to its literature expectation.

    Args:
        evidence: GSE96058 :class:`SurvivalEvidence` keyed by gene
            symbol (as produced by re-fitting
            :func:`oncocartograph.scoring.survival.fit_univariate_cox`
            on the genes in :data:`KNOWN_BIOLOGY_MARKERS`); missing or
            ``None`` entries are reported as unavailable rather than
            omitted.

    Returns:
        One :class:`KnownBiologyCheckResult` per marker in
        :data:`KNOWN_BIOLOGY_MARKERS`, in that order.
    """
    results: list[KnownBiologyCheckResult] = []
    for marker in KNOWN_BIOLOGY_MARKERS:
        marker_evidence = evidence.get(marker.gene_symbol)
        if marker_evidence is None:
            results.append(
                KnownBiologyCheckResult(
                    gene_symbol=marker.gene_symbol,
                    expected_direction=marker.expected_direction,
                    observed_hazard_ratio=None,
                    observed_direction=None,
                    plausible=None,
                )
            )
            continue

        observed = _observed_direction(marker_evidence.hazard_ratio)
        results.append(
            KnownBiologyCheckResult(
                gene_symbol=marker.gene_symbol,
                expected_direction=marker.expected_direction,
                observed_hazard_ratio=marker_evidence.hazard_ratio,
                observed_direction=observed,
                plausible=observed == marker.expected_direction,
            )
        )
    return results
