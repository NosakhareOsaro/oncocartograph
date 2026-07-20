"""Tests for oncocartograph.preprocessing.sample_manifest.

Fixture shapes mirror real GDC file-hit responses inspected during the
live pull (single-sample files for RNA-seq/methylation, paired
tumor+normal sample lists for copy number/mutation callers) -- these are
synthetic stand-ins, not real patient data.
"""

from __future__ import annotations

import pytest

from oncocartograph.preprocessing.sample_manifest import (
    PRIMARY_TUMOR_SAMPLE_TYPE,
    ResolvedFile,
    default_tie_break,
    file_case_id,
    file_sample_types,
    is_primary_tumor,
    resolve_primary_tumor_files,
)


def _hit(file_id: str, case_id: str, sample_types: list[str], file_name: str | None = None) -> dict:
    return {
        "file_id": file_id,
        "file_name": file_name or f"{file_id}.tsv",
        "cases": [{"case_id": case_id, "samples": [{"sample_type": st} for st in sample_types]}],
    }


def test_file_sample_types_extracts_single_sample() -> None:
    """A single-sample file (e.g. RNA-seq) must report exactly that one type."""
    hit = _hit("f1", "case-1", ["Primary Tumor"])
    assert file_sample_types(hit) == {"Primary Tumor"}


def test_file_sample_types_extracts_paired_tumor_normal() -> None:
    """A paired-caller file (e.g. CNV/mutation) lists both tumor and normal."""
    hit = _hit("f1", "case-1", ["Primary Tumor", "Blood Derived Normal"])
    assert file_sample_types(hit) == {"Primary Tumor", "Blood Derived Normal"}


def test_file_sample_types_handles_no_cases() -> None:
    """A file with no cases metadata must return an empty set, not error."""
    assert file_sample_types({"file_id": "f1"}) == set()


def test_file_sample_types_skips_samples_with_no_sample_type() -> None:
    """A sample entry with a missing/empty sample_type must be skipped, not added as None."""
    hit = {
        "file_id": "f1",
        "cases": [
            {
                "case_id": "case-1",
                "samples": [{"sample_type": None}, {"sample_type": "Primary Tumor"}],
            }
        ],
    }
    assert file_sample_types(hit) == {"Primary Tumor"}


def test_file_case_id_extracts_case() -> None:
    """file_case_id must return the case UUID from the first listed case."""
    hit = _hit("f1", "case-42", ["Primary Tumor"])
    assert file_case_id(hit) == "case-42"


def test_file_case_id_raises_when_no_case() -> None:
    """A file with no associated case must fail loudly, not silently skip."""
    with pytest.raises(KeyError):
        file_case_id({"file_id": "f1", "cases": []})


def test_is_primary_tumor_true_for_single_sample() -> None:
    """A single Primary Tumor sample must be recognised as such."""
    assert is_primary_tumor(_hit("f1", "case-1", ["Primary Tumor"])) is True


def test_is_primary_tumor_true_for_paired_file() -> None:
    """A paired tumor+normal file must still count as Primary Tumor
    (membership check, not exclusivity -- CNV/mutation files always list
    the germline reference sample alongside the tumor sample)."""
    assert is_primary_tumor(_hit("f1", "case-1", ["Primary Tumor", "Blood Derived Normal"])) is True


def test_is_primary_tumor_false_for_normal_only() -> None:
    """A Solid Tissue Normal-only file must not be classified as Primary Tumor."""
    assert is_primary_tumor(_hit("f1", "case-1", ["Solid Tissue Normal"])) is False


def test_default_tie_break_picks_smallest_file_id() -> None:
    """default_tie_break must be deterministic: smallest file_id wins."""
    hits = [_hit("zzz", "case-1", ["Primary Tumor"]), _hit("aaa", "case-1", ["Primary Tumor"])]
    assert default_tie_break(hits)["file_id"] == "aaa"


def test_resolve_primary_tumor_files_basic_one_per_case() -> None:
    """The common case: one Primary Tumor file per patient, all included."""
    hits = [
        _hit("f1", "case-1", ["Primary Tumor"]),
        _hit("f2", "case-2", ["Primary Tumor"]),
    ]

    resolved = resolve_primary_tumor_files(hits)

    assert resolved == {
        "case-1": ResolvedFile(case_id="case-1", file_id="f1", file_name="f1.tsv"),
        "case-2": ResolvedFile(case_id="case-2", file_id="f2", file_name="f2.tsv"),
    }


def test_resolve_primary_tumor_files_excludes_normal_and_metastatic() -> None:
    """Matched normal / metastatic files for a case must not appear in the result
    unless that case also has a genuine Primary Tumor file."""
    hits = [
        _hit("f1", "case-1", ["Primary Tumor"]),
        _hit("f2", "case-1", ["Solid Tissue Normal"]),
        _hit("f3", "case-2", ["Metastatic"]),
    ]

    resolved = resolve_primary_tumor_files(hits)

    assert set(resolved) == {"case-1"}
    assert resolved["case-1"].file_id == "f1"


def test_resolve_primary_tumor_files_includes_paired_caller_files() -> None:
    """CNV/mutation-style paired files (tumor + germline reference) must resolve normally."""
    hits = [_hit("f1", "case-1", ["Primary Tumor", "Blood Derived Normal"])]

    resolved = resolve_primary_tumor_files(hits)

    assert resolved["case-1"].file_id == "f1"


def test_resolve_primary_tumor_files_applies_tie_break_for_duplicates() -> None:
    """A case with two Primary Tumor files must be resolved via the tie-break rule."""
    hits = [
        _hit("zzz", "case-1", ["Primary Tumor"]),
        _hit("aaa", "case-1", ["Primary Tumor"]),
    ]

    resolved = resolve_primary_tumor_files(hits)

    assert len(resolved) == 1
    assert resolved["case-1"].file_id == "aaa"


def test_resolve_primary_tumor_files_respects_custom_tie_break() -> None:
    """A caller-supplied tie_break must override the default."""
    hits = [
        _hit("aaa", "case-1", ["Primary Tumor"]),
        _hit("zzz", "case-1", ["Primary Tumor"]),
    ]

    resolved = resolve_primary_tumor_files(hits, tie_break=lambda hs: hs[-1])

    assert resolved["case-1"].file_id == "zzz"


def test_primary_tumor_sample_type_constant_matches_gdc_value() -> None:
    """Guard against an accidental typo in the constant used throughout this module."""
    assert PRIMARY_TUMOR_SAMPLE_TYPE == "Primary Tumor"
