"""GSE96058 TNBC sub-cohort classification.

Applies the same rule as ``oncocartograph.data_ingestion.tnbc_cohort``
(ER/PR/HER2 all negative, excluding any missing/indeterminate marker
rather than imputing it) to GSE96058's simpler encoding: real GSE96058
``er status``/``pgr status``/``her2 status`` fields are already
histopathology-resolved to ``"0"`` (negative), ``"1"`` (positive), or
``"NA"`` (missing) -- confirmed against the real downloaded data before
this module was written, which is why this is a separate, simpler
function rather than reusing ``classify_patient`` (TCGA's IHC strings
and HER2-equivocal/FISH nuance do not apply here: HER2 status here is
already fully resolved to a single binary/NA field).

Technical replicates (a second sequencing library/run for the same
underlying tumor sample, present in GSE96058 as ``<title>repl``) must be
excluded before classification so a patient is not double-counted.
"""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

#: Value indicating a negative histopathology call for ER/PgR/HER2.
_NEGATIVE = "0"

#: Value indicating a positive histopathology call for ER/PgR/HER2.
_POSITIVE = "1"

ER_STATUS_COLUMN = "er status"
PGR_STATUS_COLUMN = "pgr status"
HER2_STATUS_COLUMN = "her2 status"

IS_TNBC_COLUMN = "is_tnbc"
EXCLUSION_REASON_COLUMN = "exclusion_reason"


@dataclass(frozen=True)
class Gse96058CohortDecision:
    """The TNBC classification outcome for a single GSE96058 sample.

    Attributes:
        is_tnbc: True if and only if ER, PgR, and HER2 status are all
            unambiguously negative.
        exclusion_reason: Human-readable reason the sample was excluded,
            or ``None`` if ``is_tnbc`` is True.
    """

    is_tnbc: bool
    exclusion_reason: str | None


def classify_gse96058_sample(
    er_status: object, pgr_status: object, her2_status: object
) -> Gse96058CohortDecision:
    """Apply the TNBC rule to one GSE96058 sample's receptor status fields.

    Args:
        er_status: Raw ``er status`` value ("0", "1", "NA", or missing).
        pgr_status: Raw ``pgr status`` value.
        her2_status: Raw ``her2 status`` value.

    Returns:
        The resulting :class:`Gse96058CohortDecision`.
    """
    calls = {"ER": er_status, "PgR": pgr_status, "HER2": her2_status}
    indeterminate = [
        name for name, value in calls.items() if str(value).strip() not in (_NEGATIVE, _POSITIVE)
    ]
    if indeterminate:
        details = "; ".join(f"{name} status missing or indeterminate" for name in indeterminate)
        return Gse96058CohortDecision(is_tnbc=False, exclusion_reason=details)

    positive = [name for name, value in calls.items() if str(value).strip() == _POSITIVE]
    if positive:
        return Gse96058CohortDecision(
            is_tnbc=False, exclusion_reason=f"Receptor positive: {', '.join(positive)}"
        )

    return Gse96058CohortDecision(is_tnbc=True, exclusion_reason=None)


def build_gse96058_cohort_audit(clinical: pd.DataFrame) -> pd.DataFrame:
    """Classify every non-replicate GSE96058 sample and return an audit table.

    Args:
        clinical: A DataFrame with at least ``title``, ``er status``,
            ``pgr status``, and ``her2 status`` columns (e.g. from
            :func:`oncocartograph.validation.gse96058_clinical.read_combined_clinical`).

    Returns:
        A copy of ``clinical``, with technical replicates (``title``
        ending in "repl", case-insensitive) dropped first, and two
        additional columns: ``is_tnbc`` (bool) and ``exclusion_reason``
        (str, or ``None`` for included samples).
    """
    required_columns = {"title", ER_STATUS_COLUMN, PGR_STATUS_COLUMN, HER2_STATUS_COLUMN}
    missing_columns = required_columns - set(clinical.columns)
    if missing_columns:
        raise KeyError(f"clinical DataFrame is missing required columns: {sorted(missing_columns)}")

    is_replicate = clinical["title"].str.contains("repl", case=False, na=False)
    primary = clinical.loc[~is_replicate].copy()

    # Column names contain spaces ("er status"), so itertuples (which
    # requires valid Python identifiers for named access) is not a good
    # fit here -- iterrows is used instead, consistent with this being a
    # small (~3,400 row) table, not a performance-sensitive hot path.
    is_tnbc_flags: list[bool] = []
    exclusion_reasons: list[str | None] = []
    for _, row in primary.iterrows():
        decision = classify_gse96058_sample(
            row[ER_STATUS_COLUMN], row[PGR_STATUS_COLUMN], row[HER2_STATUS_COLUMN]
        )
        is_tnbc_flags.append(decision.is_tnbc)
        exclusion_reasons.append(decision.exclusion_reason)

    primary[IS_TNBC_COLUMN] = is_tnbc_flags
    primary[EXCLUSION_REASON_COLUMN] = exclusion_reasons
    return primary


def select_gse96058_tnbc_cohort(audit: pd.DataFrame) -> pd.DataFrame:
    """Filter a GSE96058 cohort audit table down to the included TNBC samples.

    Args:
        audit: The output of :func:`build_gse96058_cohort_audit`.

    Returns:
        The subset of rows where ``is_tnbc`` is True, with the audit
        columns dropped.
    """
    included = audit[audit[IS_TNBC_COLUMN]].copy()
    return included.drop(columns=[IS_TNBC_COLUMN, EXCLUSION_REASON_COLUMN])
