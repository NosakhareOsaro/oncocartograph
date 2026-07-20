"""TNBC sub-cohort classification from TCGA-BRCA clinical receptor status.

Implements the exact rules recorded in
``docs/adr/0001-tnbc-cohort-definition.md`` and ``docs/methods.md`` §1.2:
ER/PR negative by the ASCO/CAP <1% IHC staining convention (Hammond et al.
2010), and HER2 negative by IHC 0/1+ or FISH-confirmed-negative IHC-2+
(Wolff et al. 2013/2018). Any patient with an indeterminate call on any of
the three markers -- including an IHC-2+ HER2 result with no recorded FISH
follow-up -- is excluded and the reason logged, never silently imputed.

This module is deliberately pure and side-effect free (no I/O) so the
classification rules themselves can be unit tested against small synthetic
clinical fixtures independent of how the clinical data was retrieved.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import pandas as pd

ReceptorCall = Literal["Positive", "Negative", "Indeterminate"]

_POSITIVE = "Positive"
_NEGATIVE = "Negative"
_EQUIVOCAL = "Equivocal"

#: Column names expected in the input clinical DataFrame.
ER_STATUS_COLUMN = "er_status_by_ihc"
PR_STATUS_COLUMN = "pr_status_by_ihc"
HER2_IHC_STATUS_COLUMN = "her2_status_by_ihc"
HER2_FISH_STATUS_COLUMN = "her2_fish_status"

#: Column names added to the audit table output.
IS_TNBC_COLUMN = "is_tnbc"
EXCLUSION_REASON_COLUMN = "exclusion_reason"


@dataclass(frozen=True)
class CohortDecision:
    """The classification outcome for a single patient.

    Attributes:
        is_tnbc: True if and only if ER, PR, and HER2 were all
            unambiguously negative per the ASCO/CAP conventions cited in
            ``docs/adr/0001-tnbc-cohort-definition.md``.
        exclusion_reason: Human-readable reason the patient was excluded
            from the TNBC cohort, or ``None`` if ``is_tnbc`` is True.
    """

    is_tnbc: bool
    exclusion_reason: str | None


def _normalise(value: object) -> str | None:
    """Normalise a raw clinical field value to a plain string or None.

    Args:
        value: A raw value from the clinical DataFrame, which may be a
            string, ``NaN`` (from pandas), or ``None``.

    Returns:
        The stripped string value, or ``None`` if the value is missing.
    """
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    text = str(value).strip()
    return text or None


def _call_er_or_pr(raw_value: object, marker_name: str) -> tuple[ReceptorCall, str | None]:
    """Classify a raw ER or PR IHC status value.

    Args:
        raw_value: Raw ``er_status_by_ihc``/``pr_status_by_ihc`` value.
        marker_name: "ER" or "PR", used only in the indeterminacy message.

    Returns:
        A tuple of (call, indeterminacy detail). The detail is ``None``
        unless the call is "Indeterminate".
    """
    value = _normalise(raw_value)
    if value == _POSITIVE:
        return "Positive", None
    if value == _NEGATIVE:
        return "Negative", None
    detail = f"{marker_name} status missing or indeterminate ({value!r})"
    return "Indeterminate", detail


def _call_her2(raw_ihc_value: object, raw_fish_value: object) -> tuple[ReceptorCall, str | None]:
    """Classify HER2 status from IHC (and, if equivocal, reflex FISH).

    Args:
        raw_ihc_value: Raw ``her2_status_by_ihc`` value.
        raw_fish_value: Raw ``her2_fish_status`` value, consulted only
            when the IHC call is Equivocal.

    Returns:
        A tuple of (call, indeterminacy detail). The detail is ``None``
        unless the call is "Indeterminate".
    """
    ihc_value = _normalise(raw_ihc_value)
    if ihc_value == _POSITIVE:
        return "Positive", None
    if ihc_value == _NEGATIVE:
        return "Negative", None
    if ihc_value == _EQUIVOCAL:
        fish_value = _normalise(raw_fish_value)
        if fish_value == _NEGATIVE:
            return "Negative", None
        if fish_value == _POSITIVE:
            return "Positive", None
        return (
            "Indeterminate",
            "HER2 IHC equivocal with no FISH follow-up recorded "
            f"(her2_fish_status={fish_value!r})",
        )
    return "Indeterminate", f"HER2 IHC status missing or indeterminate ({ihc_value!r})"


def classify_patient(
    er_status_by_ihc: object,
    pr_status_by_ihc: object,
    her2_status_by_ihc: object,
    her2_fish_status: object,
) -> CohortDecision:
    """Apply the TNBC cohort definition rules to one patient's receptor calls.

    Args:
        er_status_by_ihc: Raw ER IHC status field.
        pr_status_by_ihc: Raw PR IHC status field.
        her2_status_by_ihc: Raw HER2 IHC status field.
        her2_fish_status: Raw HER2 FISH status field (used only if the IHC
            call is Equivocal).

    Returns:
        The resulting :class:`CohortDecision`.
    """
    er_call, er_detail = _call_er_or_pr(er_status_by_ihc, "ER")
    pr_call, pr_detail = _call_er_or_pr(pr_status_by_ihc, "PR")
    her2_call, her2_detail = _call_her2(her2_status_by_ihc, her2_fish_status)

    indeterminate_details = [d for d in (er_detail, pr_detail, her2_detail) if d is not None]
    if indeterminate_details:
        return CohortDecision(is_tnbc=False, exclusion_reason="; ".join(indeterminate_details))

    positive_markers = [
        name
        for name, call in (("ER", er_call), ("PR", pr_call), ("HER2", her2_call))
        if call == _POSITIVE
    ]
    if positive_markers:
        return CohortDecision(
            is_tnbc=False,
            exclusion_reason=f"Receptor positive: {', '.join(positive_markers)}",
        )

    return CohortDecision(is_tnbc=True, exclusion_reason=None)


def build_tnbc_cohort_audit(clinical: pd.DataFrame) -> pd.DataFrame:
    """Classify every patient in a clinical DataFrame and return an audit table.

    Args:
        clinical: A DataFrame with at least the columns
            ``er_status_by_ihc``, ``pr_status_by_ihc``,
            ``her2_status_by_ihc``, and ``her2_fish_status``.

    Returns:
        A copy of ``clinical`` with two additional columns: ``is_tnbc``
        (bool) and ``exclusion_reason`` (str, or ``None`` for included
        patients). This audit table -- not just the filtered cohort -- is
        the artifact that makes the cohort definition reproducible: every
        patient's raw field values and the resulting decision are visible
        together.

    Raises:
        KeyError: If any required column is missing from ``clinical``.
    """
    required_columns = {
        ER_STATUS_COLUMN,
        PR_STATUS_COLUMN,
        HER2_IHC_STATUS_COLUMN,
        HER2_FISH_STATUS_COLUMN,
    }
    missing_columns = required_columns - set(clinical.columns)
    if missing_columns:
        raise KeyError(f"clinical DataFrame is missing required columns: {sorted(missing_columns)}")

    is_tnbc_flags: list[bool] = []
    exclusion_reasons: list[str | None] = []
    for row in clinical.itertuples(index=False):
        decision = classify_patient(
            getattr(row, ER_STATUS_COLUMN),
            getattr(row, PR_STATUS_COLUMN),
            getattr(row, HER2_IHC_STATUS_COLUMN),
            getattr(row, HER2_FISH_STATUS_COLUMN),
        )
        is_tnbc_flags.append(decision.is_tnbc)
        exclusion_reasons.append(decision.exclusion_reason)

    audit = clinical.copy()
    audit[IS_TNBC_COLUMN] = is_tnbc_flags
    audit[EXCLUSION_REASON_COLUMN] = exclusion_reasons
    return audit


def select_tnbc_cohort(audit: pd.DataFrame) -> pd.DataFrame:
    """Filter a cohort audit table down to the included TNBC patients.

    Args:
        audit: The output of :func:`build_tnbc_cohort_audit`.

    Returns:
        The subset of rows where ``is_tnbc`` is True, with the audit
        columns dropped (callers that need the audit trail should keep a
        reference to the full audit table separately).
    """
    included = audit[audit[IS_TNBC_COLUMN]].copy()
    return included.drop(columns=[IS_TNBC_COLUMN, EXCLUSION_REASON_COLUMN])
