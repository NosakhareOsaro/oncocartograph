"""Tests for oncocartograph.validation.gse96058_clinical.

The fixture format mirrors the real GSE96058 series-matrix files
(confirmed by downloading and inspecting them while planning this work
package), not synthetic guesswork about GEO's format.
"""

from __future__ import annotations

import gzip
from pathlib import Path

from oncocartograph.validation.gse96058_clinical import (
    read_combined_clinical,
    read_series_matrix,
)


def _series_matrix_row(prefix: str, *values: str) -> str:
    """Build one series-matrix row (e.g. a !Sample_characteristics_ch1 line)."""
    quoted = "\t".join(f'"{v}"' for v in values)
    return f"{prefix}\t{quoted}\n"


_SYNTHETIC_SERIES_MATRIX = (
    _series_matrix_row("!Sample_title", "F1", "F2", "F3")
    + _series_matrix_row("!Sample_geo_accession", "GSM1", "GSM2", "GSM3")
    + _series_matrix_row(
        "!Sample_characteristics_ch1", "er status: 0", "er status: 1", "er status: NA"
    )
    + _series_matrix_row(
        "!Sample_characteristics_ch1", "pgr status: 0", "pgr status: 1", "pgr status: NA"
    )
    + _series_matrix_row(
        "!Sample_characteristics_ch1", "her2 status: 0", "her2 status: 0", "her2 status: NA"
    )
    + _series_matrix_row(
        "!Sample_characteristics_ch1",
        "overall survival days: 1000",
        "overall survival days: 2000",
        "overall survival days: 500",
    )
    + _series_matrix_row(
        "!Sample_characteristics_ch1",
        "overall survival event: 0",
        "overall survival event: 1",
        "overall survival event: 0",
    )
    + _series_matrix_row(
        "!Sample_characteristics_ch1",
        "pam50 subtype: Basal",
        "pam50 subtype: LumA",
        "pam50 subtype: Basal",
    )
)


def test_read_series_matrix_parses_plain_text(tmp_path: Path) -> None:
    """Fields must be extracted by name (before the colon), one row per sample."""
    path = tmp_path / "series_matrix.txt"
    path.write_text(_SYNTHETIC_SERIES_MATRIX)

    table = read_series_matrix(path)

    assert len(table) == 3
    assert list(table["title"]) == ["F1", "F2", "F3"]
    assert list(table["geo_accession"]) == ["GSM1", "GSM2", "GSM3"]
    assert list(table["er status"]) == ["0", "1", "NA"]
    assert list(table["overall survival event"]) == ["0", "1", "0"]
    assert list(table["pam50 subtype"]) == ["Basal", "LumA", "Basal"]


def test_read_series_matrix_parses_gzip(tmp_path: Path) -> None:
    """A .gz series matrix (the real GEO download format) must parse identically."""
    path = tmp_path / "series_matrix.txt.gz"
    with gzip.open(path, "wt", encoding="utf-8") as f:
        f.write(_SYNTHETIC_SERIES_MATRIX)

    table = read_series_matrix(path)

    assert len(table) == 3
    assert list(table["er status"]) == ["0", "1", "NA"]


def test_read_combined_clinical_concatenates_multiple_files(tmp_path: Path) -> None:
    """Two platform-specific files (GSE96058's real split) must combine into one table."""
    path_a = tmp_path / "gpl_a.txt"
    path_a.write_text(_SYNTHETIC_SERIES_MATRIX)
    path_b = tmp_path / "gpl_b.txt"
    path_b.write_text(
        '!Sample_title\t"F4"\n'
        '!Sample_geo_accession\t"GSM4"\n'
        '!Sample_characteristics_ch1\t"er status: 1"\n'
        '!Sample_characteristics_ch1\t"pgr status: 1"\n'
        '!Sample_characteristics_ch1\t"her2 status: 0"\n'
        '!Sample_characteristics_ch1\t"overall survival days: 300"\n'
        '!Sample_characteristics_ch1\t"overall survival event: 1"\n'
        '!Sample_characteristics_ch1\t"pam50 subtype: LumB"\n'
    )

    combined = read_combined_clinical([path_a, path_b])

    assert len(combined) == 4
    assert list(combined["title"]) == ["F1", "F2", "F3", "F4"]


def test_read_series_matrix_ignores_unrelated_lines(tmp_path: Path) -> None:
    """Lines that are neither title/geo_accession/characteristics (e.g. series-level
    metadata) must be ignored, not crash or corrupt parsing."""
    path = tmp_path / "series_matrix.txt"
    path.write_text(
        "^SERIES = GSE96058\n"
        "!Series_title = Some title\n"
        + _series_matrix_row("!Sample_title", "F1")
        + _series_matrix_row("!Sample_geo_accession", "GSM1")
        + _series_matrix_row("!Sample_characteristics_ch1", "er status: 0")
    )

    table = read_series_matrix(path)

    assert list(table["title"]) == ["F1"]
    assert list(table["er status"]) == ["0"]


def test_read_series_matrix_skips_malformed_characteristics_row(tmp_path: Path) -> None:
    """A characteristics row with no 'field: value' colon must be skipped, not crash."""
    path = tmp_path / "series_matrix.txt"
    path.write_text(
        _series_matrix_row("!Sample_title", "F1")
        + _series_matrix_row("!Sample_geo_accession", "GSM1")
        + _series_matrix_row("!Sample_characteristics_ch1", "no colon here")
        + _series_matrix_row("!Sample_characteristics_ch1", "er status: 0")
    )

    table = read_series_matrix(path)

    assert list(table["er status"]) == ["0"]
