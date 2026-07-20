"""Tests for oncocartograph.preprocessing.copy_number.

Fixture values are synthetic; workflow_type strings mirror the real GDC
values confirmed via a live metadata query (ASCAT3/ASCAT2/ABSOLUTE
LiftOver/AscatNGS).
"""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd
import pytest

from oncocartograph.preprocessing.copy_number import (
    DIPLOID_COPY_NUMBER,
    build_copy_number_matrix,
    read_gene_level_copy_number,
    relative_log2_copy_number,
    resolve_copy_number_files,
    workflow_priority_tie_break,
)


def _hit(
    file_id: str, case_id: str, workflow_type: str, sample_type: str = "Primary Tumor"
) -> dict:
    return {
        "file_id": file_id,
        "file_name": f"{file_id}.tsv",
        "analysis": {"workflow_type": workflow_type},
        "cases": [{"case_id": case_id, "samples": [{"sample_type": sample_type}]}],
    }


def test_workflow_priority_tie_break_prefers_ascat3() -> None:
    """Given a choice, ASCAT3 must win over ASCAT2/ABSOLUTE LiftOver/AscatNGS."""
    hits = [
        _hit("f-absolute", "case-1", "ABSOLUTE LiftOver"),
        _hit("f-ascat2", "case-1", "ASCAT2"),
        _hit("f-ascat3", "case-1", "ASCAT3"),
    ]
    chosen = workflow_priority_tie_break()(hits)
    assert chosen["file_id"] == "f-ascat3"


def test_workflow_priority_tie_break_falls_back_through_priority_order() -> None:
    """Without ASCAT3, ASCAT2 must be preferred over ABSOLUTE LiftOver/AscatNGS."""
    hits = [
        _hit("f-ngs", "case-1", "AscatNGS"),
        _hit("f-absolute", "case-1", "ABSOLUTE LiftOver"),
        _hit("f-ascat2", "case-1", "ASCAT2"),
    ]
    chosen = workflow_priority_tie_break()(hits)
    assert chosen["file_id"] == "f-ascat2"


def test_workflow_priority_tie_break_falls_back_to_default_for_unknown_workflow() -> None:
    """A workflow_type not in the priority list must fall back to the deterministic default."""
    hits = [
        _hit("zzz", "case-1", "SomeFutureWorkflow"),
        _hit("aaa", "case-1", "SomeFutureWorkflow"),
    ]
    chosen = workflow_priority_tie_break()(hits)
    assert chosen["file_id"] == "aaa"


def test_resolve_copy_number_files_picks_ascat3_when_available() -> None:
    """End-to-end: resolve_copy_number_files must select ASCAT3 for a patient
    with multiple workflow variants available."""
    hits = [
        _hit("f-ascat2", "case-1", "ASCAT2"),
        _hit("f-ascat3", "case-1", "ASCAT3"),
        _hit("f-single", "case-2", "AscatNGS"),
    ]

    resolved = resolve_copy_number_files(hits)

    assert resolved["case-1"].file_id == "f-ascat3"
    assert resolved["case-2"].file_id == "f-single"


def test_read_gene_level_copy_number_parses_values_and_missing(tmp_path: Path) -> None:
    """Parsing must index by gene_id and leave uncalled genes as NaN, not 0."""
    content = (
        "gene_id\tgene_name\tchromosome\tstart\tend\tcopy_number\tmin_copy_number\tmax_copy_number\n"
        "ENSG00000000001.1\tGENE1\tchr1\t100\t200\t2\t2\t2\n"
        "ENSG00000000002.1\tGENE2\tchr1\t300\t400\t\t\t\n"
        "ENSG00000000003.1\tGENE3\tchr1\t500\t600\t4\t4\t4\n"
    )
    path = tmp_path / "sample.gene_level_copy_number.v36.tsv"
    path.write_text(content)

    series = read_gene_level_copy_number(path)

    assert series["ENSG00000000001.1"] == 2
    assert math.isnan(series["ENSG00000000002.1"])
    assert series["ENSG00000000003.1"] == 4


def test_relative_log2_copy_number_zeroes_diploid_baseline() -> None:
    """Diploid (2) must map to exactly 0; the naive (CN+1)/2 formula would not."""
    result = relative_log2_copy_number(pd.Series([DIPLOID_COPY_NUMBER]))
    assert result.iloc[0] == pytest.approx(0.0)


def test_relative_log2_copy_number_is_finite_at_homozygous_deletion() -> None:
    """CN=0 (homozygous deletion) must produce a finite negative value, not -inf."""
    result = relative_log2_copy_number(pd.Series([0]))
    assert math.isfinite(result.iloc[0])
    assert result.iloc[0] < 0


def test_relative_log2_copy_number_is_monotonic_increasing() -> None:
    """Higher absolute copy number must always produce a higher relative value."""
    result = relative_log2_copy_number(pd.Series([0, 1, 2, 3, 4, 6]))
    assert list(result) == sorted(result)


def test_build_copy_number_matrix_combines_patients_into_gene_by_patient(tmp_path: Path) -> None:
    """The final matrix must be gene-indexed with one column per patient."""
    content_a = (
        "gene_id\tgene_name\tchromosome\tstart\tend\tcopy_number\tmin_copy_number\tmax_copy_number\n"
        "ENSG00000000001.1\tGENE1\tchr1\t100\t200\t2\t2\t2\n"
    )
    content_b = (
        "gene_id\tgene_name\tchromosome\tstart\tend\tcopy_number\tmin_copy_number\tmax_copy_number\n"
        "ENSG00000000001.1\tGENE1\tchr1\t100\t200\t4\t4\t4\n"
    )
    path_a = tmp_path / "a.tsv"
    path_a.write_text(content_a)
    path_b = tmp_path / "b.tsv"
    path_b.write_text(content_b)

    matrix = build_copy_number_matrix({"case-a": path_a, "case-b": path_b})

    assert list(matrix.columns) == ["case-a", "case-b"]
    assert matrix.loc["ENSG00000000001.1", "case-a"] == pytest.approx(0.0)
    assert matrix.loc["ENSG00000000001.1", "case-b"] == pytest.approx(
        relative_log2_copy_number(pd.Series([4])).iloc[0]
    )
