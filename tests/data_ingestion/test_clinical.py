"""Tests for oncocartograph.data_ingestion.clinical.

The biotab fixture file written here is a synthetic stand-in for the real
TCGA-BRCA clinical supplement layout, not real patient data.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from oncocartograph.data_ingestion.clinical import (
    RECEPTOR_STATUS_COLUMNS,
    derive_survival_outcome,
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


def _survival_clinical() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "bcr_patient_uuid": "u1",
                "vital_status": "Dead",
                "death_days_to": "967",
                "last_contact_days_to": "[Not Applicable]",
            },
            {
                "bcr_patient_uuid": "u2",
                "vital_status": "Alive",
                "death_days_to": "[Not Applicable]",
                "last_contact_days_to": "852",
            },
            {
                "bcr_patient_uuid": "u3",
                "vital_status": "Alive",
                "death_days_to": "[Not Applicable]",
                "last_contact_days_to": "[Not Available]",
            },
        ]
    )


def test_derive_survival_outcome_uses_death_days_for_events() -> None:
    """A dead patient's duration must come from death_days_to, event=1."""
    outcome = derive_survival_outcome(_survival_clinical())

    assert outcome.loc["u1", "duration"] == 967.0
    assert outcome.loc["u1", "event"] == 1


def test_derive_survival_outcome_uses_last_contact_days_for_censored() -> None:
    """A living patient's duration must come from last_contact_days_to, event=0."""
    outcome = derive_survival_outcome(_survival_clinical())

    assert outcome.loc["u2", "duration"] == 852.0
    assert outcome.loc["u2", "event"] == 0


def test_derive_survival_outcome_drops_patients_missing_both_day_fields() -> None:
    """A patient with neither day field parseable must be dropped, not kept with NaN duration."""
    outcome = derive_survival_outcome(_survival_clinical())

    assert "u3" not in outcome.index
    assert len(outcome) == 2


def test_derive_survival_outcome_raises_on_missing_column() -> None:
    """A clinical DataFrame missing a required survival column must fail loudly."""
    clinical = _survival_clinical().drop(columns=["vital_status"])

    with pytest.raises(KeyError, match="vital_status"):
        derive_survival_outcome(clinical)
