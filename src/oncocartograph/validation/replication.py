"""External replication of TCGA survival-association direction in GSE96058.

Given that the TCGA-BRCA screen itself produced 0/709 FDR-significant
biomarkers (see docs/methods.md), demanding that GSE96058 replicate
nominal statistical significance would be incoherent -- there is nothing
significant to replicate. The pre-registered primary criterion is
therefore direction concordance: does the sign of the log hazard ratio
in GSE96058 agree with TCGA's sign, more often than the 50% expected by
chance, tested with a one-sided exact binomial test at alpha=0.05.
Nominal p-value replication is reported alongside as secondary,
informational-only context.

This criterion (and alpha) were fixed before the analysis was run, in
this module's docstring and the accompanying ADR, precisely so the
pass/fail bar cannot be quietly redefined after seeing the result.
"""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd
from scipy.stats import binomtest

from oncocartograph.scoring.survival import SurvivalEvidence

#: The pre-registered significance threshold for the direction-concordance test.
PRE_REGISTERED_ALPHA = 0.05

#: The pre-registered null concordance rate (chance agreement between two
#: independent, arbitrarily-signed hazard ratios).
CHANCE_CONCORDANCE_RATE = 0.5


@dataclass(frozen=True)
class CandidateReplication:
    """One TCGA candidate's replication status in GSE96058.

    Attributes:
        candidate_id: The original TCGA candidate identifier (e.g.
            ``"rna_seq:ENSG00000104267.10"``).
        gene_symbol: The gene symbol used to look it up in GSE96058, or
            ``None`` if it could not be resolved.
        tcga_hazard_ratio: The candidate's hazard ratio in TCGA.
        gse96058_evidence: The candidate's re-fit :class:`SurvivalEvidence`
            in GSE96058, or ``None`` if it could not be resolved (no gene
            symbol) or could not be fit (absent from the expression
            matrix, or a degenerate Cox fit).
        concordant: ``True``/``False`` if both directions are known,
            ``None`` if GSE96058 evidence is unavailable.
    """

    candidate_id: str
    gene_symbol: str | None
    tcga_hazard_ratio: float
    gse96058_evidence: SurvivalEvidence | None
    concordant: bool | None


def _direction_sign(hazard_ratio: float) -> int:
    """Sign of the log hazard ratio: +1 (harmful), -1 (protective), 0 (HR==1)."""
    if hazard_ratio > 1:
        return 1
    if hazard_ratio < 1:
        return -1
    return 0


def build_replication_table(
    tcga_hazard_ratios: dict[str, float],
    gene_symbols: dict[str, str | None],
    gse96058_evidence: dict[str, SurvivalEvidence | None],
) -> list[CandidateReplication]:
    """Assemble per-candidate replication records.

    Args:
        tcga_hazard_ratios: TCGA hazard ratio per candidate ID (e.g.
            ``"rna_seq:ENSG00000104267.10"``).
        gene_symbols: Gene symbol per candidate ID, as resolved via
            :meth:`OpenTargetsClient.fetch_targets`; ``None`` for
            candidates whose Ensembl ID could not be resolved.
        gse96058_evidence: Re-fit GSE96058 :class:`SurvivalEvidence` per
            gene symbol (not candidate ID); ``None`` or absent for genes
            not found in the expression matrix or whose Cox model did
            not converge.

    Returns:
        One :class:`CandidateReplication` per candidate in
        ``tcga_hazard_ratios``, in the same order.
    """
    records: list[CandidateReplication] = []
    for candidate_id, tcga_hr in tcga_hazard_ratios.items():
        symbol = gene_symbols.get(candidate_id)
        evidence = gse96058_evidence.get(symbol) if symbol is not None else None

        concordant: bool | None = None
        if evidence is not None:
            tcga_sign = _direction_sign(tcga_hr)
            gse_sign = _direction_sign(evidence.hazard_ratio)
            concordant = tcga_sign != 0 and tcga_sign == gse_sign

        records.append(
            CandidateReplication(
                candidate_id=candidate_id,
                gene_symbol=symbol,
                tcga_hazard_ratio=tcga_hr,
                gse96058_evidence=evidence,
                concordant=concordant,
            )
        )
    return records


@dataclass(frozen=True)
class ConcordanceTestResult:
    """Outcome of the pre-registered direction-concordance test.

    Attributes:
        n_total_candidates: Total TCGA candidates considered.
        n_fittable: Candidates with usable GSE96058 evidence (gene
            symbol resolved, present in the expression matrix, and a
            convergent Cox fit) -- the denominator of the test.
        n_concordant: Of ``n_fittable``, how many agreed in direction.
        concordance_rate: ``n_concordant / n_fittable``.
        p_value: One-sided exact binomial test p-value against the
            :data:`CHANCE_CONCORDANCE_RATE` null.
        alpha: The pre-registered significance threshold.
        success: ``p_value < alpha`` -- the single pre-registered
            pass/fail verdict for this validation.
    """

    n_total_candidates: int
    n_fittable: int
    n_concordant: int
    concordance_rate: float
    p_value: float
    alpha: float
    success: bool


def run_direction_concordance_test(
    replications: list[CandidateReplication], *, alpha: float = PRE_REGISTERED_ALPHA
) -> ConcordanceTestResult:
    """Run the pre-registered one-sided binomial concordance test.

    Args:
        replications: Output of :func:`build_replication_table`.
        alpha: Significance threshold (default: the pre-registered 0.05).

    Returns:
        A :class:`ConcordanceTestResult`. Raises no error and reframes
        nothing based on the outcome -- callers must report ``success``
        as-is, including a ``False`` result.
    """
    fittable = [r for r in replications if r.concordant is not None]
    n_fittable = len(fittable)
    n_concordant = sum(1 for r in fittable if r.concordant)

    test = binomtest(n_concordant, n_fittable, p=CHANCE_CONCORDANCE_RATE, alternative="greater")

    return ConcordanceTestResult(
        n_total_candidates=len(replications),
        n_fittable=n_fittable,
        n_concordant=n_concordant,
        concordance_rate=n_concordant / n_fittable,
        p_value=float(test.pvalue),
        alpha=alpha,
        success=float(test.pvalue) < alpha,
    )


def replication_table_to_frame(replications: list[CandidateReplication]) -> pd.DataFrame:
    """Flatten replication records into a DataFrame for reporting/export."""
    rows = [
        {
            "candidate_id": r.candidate_id,
            "gene_symbol": r.gene_symbol,
            "tcga_hazard_ratio": r.tcga_hazard_ratio,
            "gse96058_hazard_ratio": (
                r.gse96058_evidence.hazard_ratio if r.gse96058_evidence else None
            ),
            "gse96058_p_value": (r.gse96058_evidence.p_value if r.gse96058_evidence else None),
            "gse96058_n_samples": (r.gse96058_evidence.n_samples if r.gse96058_evidence else None),
            "concordant": r.concordant,
        }
        for r in replications
    ]
    return pd.DataFrame(rows)
