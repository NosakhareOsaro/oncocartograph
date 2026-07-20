"""Parsing for the TCGA-BRCA clinical supplement (BCR Biotab format).

The GDC "Clinical Supplement" data type for TCGA projects is distributed as
tab-delimited BCR Biotab files with a fixed-layout header: the first row is
the real column name, followed by a small number of metadata rows (a
human-readable description row and a controlled-vocabulary/CDE-ID row)
before the actual patient data begins. This module isolates that
file-format quirk from the rest of the ingestion pipeline.

The exact column names asserted here (``er_status_by_ihc``, etc.) match
the biotab layout as documented for TCGA-BRCA; per
``docs/adr/0004-gdc-rest-client-over-tcgabiolinks.md``, this is flagged as
a known item to validate against a real downloaded file during the first
live pull, and this module will be updated if the real file's column names
differ.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

#: Number of non-header metadata rows following the column-name row in a
#: BCR Biotab file, before actual patient data begins.
BIOTAB_METADATA_ROWS = 2

#: Clinical columns required for TNBC cohort classification
#: (see oncocartograph.data_ingestion.tnbc_cohort).
RECEPTOR_STATUS_COLUMNS = (
    "bcr_patient_barcode",
    "er_status_by_ihc",
    "pr_status_by_ihc",
    "her2_status_by_ihc",
    "her2_fish_status",
)

#: Clinical columns required to derive overall-survival duration/event
#: (see :func:`derive_survival_outcome`). Confirmed against the real
#: downloaded clinical file (column names match exactly).
SURVIVAL_COLUMNS = (
    "bcr_patient_uuid",
    "vital_status",
    "death_days_to",
    "last_contact_days_to",
)

#: The ``vital_status`` value indicating the death event was observed.
_DEAD = "Dead"


def read_biotab(path: Path, *, metadata_rows: int = BIOTAB_METADATA_ROWS) -> pd.DataFrame:
    """Read a BCR Biotab tab-delimited clinical file into a DataFrame.

    Args:
        path: Path to the biotab file (already downloaded).
        metadata_rows: Number of metadata rows to discard immediately
            after the header row (defaults to the standard TCGA biotab
            layout of 2: a description row and a CDE-ID row).

    Returns:
        A DataFrame with one row per patient and columns named from the
        biotab file's header row, with metadata rows removed.
    """
    table = pd.read_csv(path, sep="\t", header=0, dtype=str)
    return table.iloc[metadata_rows:].reset_index(drop=True)


def extract_receptor_status(clinical: pd.DataFrame) -> pd.DataFrame:
    """Select and validate the columns needed for TNBC cohort classification.

    Args:
        clinical: A parsed clinical DataFrame, e.g. from :func:`read_biotab`.

    Returns:
        A DataFrame restricted to :data:`RECEPTOR_STATUS_COLUMNS`, in that
        column order.

    Raises:
        KeyError: If any required column is absent, naming exactly which
            ones -- this fails loudly rather than silently proceeding with
            a partial cohort definition.
    """
    missing = [c for c in RECEPTOR_STATUS_COLUMNS if c not in clinical.columns]
    if missing:
        raise KeyError(
            f"Clinical data is missing required receptor status columns: {missing}. "
            "See docs/adr/0004-gdc-rest-client-over-tcgabiolinks.md -- the biotab "
            "column layout may need updating against the real downloaded file."
        )
    return clinical[list(RECEPTOR_STATUS_COLUMNS)].copy()


def _to_days_or_none(value: object) -> float | None:
    """Coerce a raw clinical day-count field to float.

    Args:
        value: A raw ``death_days_to``/``last_contact_days_to`` value,
            which may already be numeric-looking text or a TCGA sentinel
            string (e.g. ``"[Not Available]"``).

    Returns:
        The value as a float, or ``None`` if it cannot be parsed as one.
    """
    try:
        return float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None


def derive_survival_outcome(clinical: pd.DataFrame) -> pd.DataFrame:
    """Derive overall-survival duration and event indicators from raw clinical fields.

    Cox PH on overall survival is the only model TCGA-BRCA's clinical file
    supports (see ``docs/adr/0007-survival-methodology-and-composite-score.md``);
    this function produces the ``duration``/``event`` inputs
    :func:`oncocartograph.scoring.survival.fit_univariate_cox` and
    :func:`~oncocartograph.scoring.survival.screen_survival_associations`
    expect.

    Args:
        clinical: A parsed clinical DataFrame (e.g. from :func:`read_biotab`)
            with at least :data:`SURVIVAL_COLUMNS`.

    Returns:
        A DataFrame indexed by ``bcr_patient_uuid`` with ``duration``
        (float days) and ``event`` (int 0/1) columns. ``duration`` is
        ``death_days_to`` when the patient died (an observed event), else
        ``last_contact_days_to`` (censored at last contact). Patients
        missing both day fields are dropped, since neither a duration nor
        a censoring time is available for them.

    Raises:
        KeyError: If any of :data:`SURVIVAL_COLUMNS` is missing.
    """
    missing = [c for c in SURVIVAL_COLUMNS if c not in clinical.columns]
    if missing:
        raise KeyError(f"clinical DataFrame is missing required survival columns: {missing}")

    death_days = clinical["death_days_to"].apply(_to_days_or_none)
    last_contact_days = clinical["last_contact_days_to"].apply(_to_days_or_none)
    duration = death_days.fillna(last_contact_days)
    event = (clinical["vital_status"] == _DEAD).astype(int)

    # Built positionally (via .to_numpy()), not by passing `index=` to the
    # DataFrame constructor alongside Series that still carry the original
    # RangeIndex -- doing so silently reindexes duration/event against the
    # bcr_patient_uuid labels (which don't match integer positions),
    # producing an all-NaN result. Caught while writing this function's
    # first test.
    outcome = pd.DataFrame(
        {
            "bcr_patient_uuid": clinical["bcr_patient_uuid"].to_numpy(),
            "duration": duration.to_numpy(),
            "event": event.to_numpy(),
        }
    ).set_index("bcr_patient_uuid")
    return outcome.dropna(subset=["duration"])
