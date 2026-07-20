"""Parsing for the GSE96058 (SCAN-B) series-matrix clinical/survival metadata.

GEO series-matrix files store one field per row and one sample per column
(``!Sample_title``, ``!Sample_geo_accession``, and repeated
``!Sample_characteristics_ch1`` rows in ``"field: value"`` format) --
confirmed against the real downloaded files before this module was
written. GSE96058 splits its 3,409 samples (3,273 primary + 136 technical
replicates) across two platform-specific series-matrix files
(GPL11154/HiSeq 2000, GPL18573/NextSeq 500); both must be parsed and
combined to get the full cohort.
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import pandas as pd

_CHARACTERISTIC_PATTERN = re.compile(r"^([^:]+):\s*(.*)$")


def _open_text(path: Path):  # type: ignore[no-untyped-def]
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open(encoding="utf-8", errors="replace")


def read_series_matrix(path: Path) -> pd.DataFrame:
    """Parse one GEO series-matrix file into a per-sample clinical DataFrame.

    Args:
        path: Path to a ``*_series_matrix.txt`` or ``*_series_matrix.txt.gz``
            file.

    Returns:
        A DataFrame with one row per sample and columns for ``title``,
        ``geo_accession``, and every ``!Sample_characteristics_ch1``
        field found (column names taken from the field name before the
        first colon, e.g. ``"er status"``).
    """
    fields: dict[str, list[str | None]] = {}
    geo_accessions: list[str] = []
    with _open_text(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("!Sample_geo_accession"):
                geo_accessions = [p.strip('"') for p in line.split("\t")[1:]]
            elif line.startswith("!Sample_title"):
                fields["title"] = [p.strip('"') for p in line.split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                values = [p.strip('"') for p in line.split("\t")[1:]]
                match = _CHARACTERISTIC_PATTERN.match(values[0])
                if not match:
                    continue
                field_name = match.group(1).strip()
                value_matches = [_CHARACTERISTIC_PATTERN.match(v) for v in values]
                fields[field_name] = [m.group(2).strip() if m else None for m in value_matches]

    table = pd.DataFrame(fields)
    table["geo_accession"] = geo_accessions
    return table


def read_combined_clinical(paths: list[Path]) -> pd.DataFrame:
    """Parse and combine multiple platform-specific series-matrix files.

    Args:
        paths: Paths to each platform's series-matrix file (e.g. GSE96058's
            GPL11154 and GPL18573 files).

    Returns:
        The row-wise concatenation of :func:`read_series_matrix` applied to
        each path, with a fresh integer index.
    """
    tables = [read_series_matrix(path) for path in paths]
    return pd.concat(tables, ignore_index=True)
