"""Tests for oncocartograph.validation.gse96058_expression.

The fixture format (quoted header, gene symbol first column, log2 values)
mirrors the real GSE96058 expression file, confirmed by downloading and
inspecting a real chunk of it while planning this work package.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from oncocartograph.validation.gse96058_expression import read_selected_gene_expression

_SYNTHETIC_EXPRESSION_CSV = (
    '"","F1","F2","F3"\n'
    '"5_8S_rRNA",-3.32192809488736,-3.32192809488736,-3.32192809488736\n'
    '"TP53",5.1,4.9,6.2\n'
    '"AR",2.1,7.8,3.0\n'
    '"GAPDH",10.1,10.2,9.9\n'
)


@pytest.fixture
def expression_file(tmp_path: Path) -> Path:
    path = tmp_path / "expression.csv.gz"
    with gzip.open(path, "wt", encoding="utf-8") as f:
        f.write(_SYNTHETIC_EXPRESSION_CSV)
    return path


def test_read_selected_gene_expression_extracts_requested_genes(expression_file: Path) -> None:
    """Only the requested genes must appear in the result, with correct sample columns."""
    result = read_selected_gene_expression(expression_file, ["TP53", "AR"])

    assert set(result.index) == {"TP53", "AR"}
    assert list(result.columns) == ["F1", "F2", "F3"]
    assert result.loc["TP53", "F1"] == pytest.approx(5.1)
    assert result.loc["AR", "F2"] == pytest.approx(7.8)


def test_read_selected_gene_expression_ignores_unrequested_genes(expression_file: Path) -> None:
    """Genes present in the file but not requested must not appear in the result."""
    result = read_selected_gene_expression(expression_file, ["TP53"])

    assert "GAPDH" not in result.index
    assert "5_8S_rRNA" not in result.index


def test_read_selected_gene_expression_missing_requested_gene_is_simply_absent(
    expression_file: Path,
) -> None:
    """A requested gene not present in the file must be silently absent, not an error."""
    result = read_selected_gene_expression(expression_file, ["TP53", "NOTAREALGENE"])

    assert set(result.index) == {"TP53"}


def test_read_selected_gene_expression_empty_request_returns_empty_frame(
    expression_file: Path,
) -> None:
    """Requesting zero genes must return an empty DataFrame, not error."""
    result = read_selected_gene_expression(expression_file, [])

    assert result.empty
