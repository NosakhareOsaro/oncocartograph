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
