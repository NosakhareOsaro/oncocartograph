"""Tests for oncocartograph.data_ingestion.clinical.

The biotab fixture file written here is a synthetic stand-in for the real
TCGA-BRCA clinical supplement layout, not real patient data.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from oncocartograph.data_ingestion.clinical import (
    RECEPTOR_STATUS_COLUMNS,
    extract_receptor_status,
    read_biotab,
)

_SYNTHETIC_BIOTAB_CONTENT = (
    "bcr_patient_barcode\ter_status_by_ihc\tpr_status_by_ihc\ther2_status_by_ihc\ther2_fish_status\textra_column\n"
    "Patient identifier\tER status\tPR status\tHER2 IHC status\tHER2 FISH status\tExtra\n"
    "CDE_ID:2003301\tCDE_ID:2957503\tCDE_ID:2957504\tCDE_ID:2957505\tCDE_ID:2957506\tCDE_ID:0000000\n"
    "SYNTH-0001\tNegative\tNegative\tNegative\t[Not Applicable]\tfoo\n"
    "SYNTH-0002\tPositive\tNegative\tNegative\t[Not Applicable]\tbar\n"
)


@pytest.fixture
def synthetic_biotab_file(tmp_path: Path) -> Path:
    """Write a small synthetic biotab-format file and return its path."""
    path = tmp_path / "synthetic_clinical_patient_brca.txt"
    path.write_text(_SYNTHETIC_BIOTAB_CONTENT)
    return path


def test_read_biotab_drops_metadata_rows_and_keeps_data(synthetic_biotab_file: Path) -> None:
    """The two metadata rows (description, CDE_ID) must be dropped, data rows kept."""
    table = read_biotab(synthetic_biotab_file)

    assert list(table.columns) == [
        "bcr_patient_barcode",
        "er_status_by_ihc",
        "pr_status_by_ihc",
        "her2_status_by_ihc",
        "her2_fish_status",
        "extra_column",
    ]
    assert len(table) == 2
    assert table.iloc[0]["bcr_patient_barcode"] == "SYNTH-0001"
    assert table.iloc[1]["er_status_by_ihc"] == "Positive"


def test_extract_receptor_status_selects_expected_columns(synthetic_biotab_file: Path) -> None:
    """extract_receptor_status must return exactly RECEPTOR_STATUS_COLUMNS, in order."""
    table = read_biotab(synthetic_biotab_file)

    receptor_status = extract_receptor_status(table)

    assert list(receptor_status.columns) == list(RECEPTOR_STATUS_COLUMNS)
    assert len(receptor_status) == 2


def test_extract_receptor_status_raises_on_missing_column(synthetic_biotab_file: Path) -> None:
    """A clinical table missing a required column must fail with a clear message."""
    table = read_biotab(synthetic_biotab_file).drop(columns=["her2_fish_status"])

    with pytest.raises(KeyError, match="her2_fish_status"):
        extract_receptor_status(table)
