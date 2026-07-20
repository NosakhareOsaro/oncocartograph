"""Tests for oncocartograph.preprocessing.mutation.

MAF fixtures are synthetic, shaped after the real column layout and
Variant_Classification vocabulary confirmed against downloaded MAF files
during this work package.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd

from oncocartograph.preprocessing.mutation import (
    NON_SYNONYMOUS_CLASSIFICATIONS,
    build_mutation_matrix,
    filter_by_recurrence,
    read_non_synonymous_genes,
)

_MAF_HEADER = "Hugo_Symbol\tVariant_Classification\tOther_Column\n"


def _write_maf(path: Path, rows: list[tuple[str, str]], compress: bool = False) -> Path:
    content = "#version 2.4\n#filter comment line\n" + _MAF_HEADER
    content += "".join(f"{gene}\t{classification}\tignored\n" for gene, classification in rows)
    if compress:
        path.write_bytes(gzip.compress(content.encode()))
    else:
        path.write_text(content)
    return path


def test_read_non_synonymous_genes_keeps_only_allowed_classifications(tmp_path: Path) -> None:
    """Silent/Intron/etc. must be excluded; Missense/Nonsense/etc. kept."""
    path = _write_maf(
        tmp_path / "sample.maf",
        [
            ("TP53", "Missense_Mutation"),
            ("PIK3CA", "Silent"),
            ("BRCA1", "Nonsense_Mutation"),
            ("PTEN", "Intron"),
        ],
    )

    genes = read_non_synonymous_genes(path)

    assert genes == {"TP53", "BRCA1"}


def test_read_non_synonymous_genes_handles_gzip_compressed_maf(tmp_path: Path) -> None:
    """Real MC3 MAFs are gzip-compressed; parsing must work transparently."""
    path = _write_maf(tmp_path / "sample.maf.gz", [("TP53", "Missense_Mutation")], compress=True)

    genes = read_non_synonymous_genes(path)

    assert genes == {"TP53"}


def test_read_non_synonymous_genes_deduplicates_multiple_variants_same_gene(tmp_path: Path) -> None:
    """Two non-synonymous variants in the same gene must count as one gene, not two."""
    path = _write_maf(
        tmp_path / "sample.maf",
        [("TP53", "Missense_Mutation"), ("TP53", "Nonsense_Mutation")],
    )

    genes = read_non_synonymous_genes(path)

    assert genes == {"TP53"}


def test_non_synonymous_classifications_excludes_silent_and_intron() -> None:
    """Guard against accidentally including passenger-like categories in the allow-list."""
    assert "Silent" not in NON_SYNONYMOUS_CLASSIFICATIONS
    assert "Intron" not in NON_SYNONYMOUS_CLASSIFICATIONS
    assert "Missense_Mutation" in NON_SYNONYMOUS_CLASSIFICATIONS


def test_build_mutation_matrix_produces_binary_zero_not_nan(tmp_path: Path) -> None:
    """A patient with no variant in a gene must get 0, not NaN (unlike other omic views)."""
    path_a = _write_maf(tmp_path / "a.maf", [("TP53", "Missense_Mutation")])
    path_b = _write_maf(tmp_path / "b.maf", [("BRCA1", "Nonsense_Mutation")])

    matrix = build_mutation_matrix({"case-a": path_a, "case-b": path_b})

    assert matrix.loc["TP53", "case-a"] == 1
    assert matrix.loc["TP53", "case-b"] == 0
    assert matrix.loc["BRCA1", "case-a"] == 0
    assert matrix.loc["BRCA1", "case-b"] == 1
    assert not matrix.isna().to_numpy().any()


def test_build_mutation_matrix_handles_empty_input() -> None:
    """No resolved files must produce an empty matrix, not an error."""
    matrix = build_mutation_matrix({})
    assert matrix.empty


def test_filter_by_recurrence_keeps_genes_meeting_threshold() -> None:
    """Genes mutated in fewer than min_patients patients must be dropped."""
    matrix = pd.DataFrame(
        {"case-a": [1, 1, 1], "case-b": [1, 1, 0], "case-c": [1, 0, 0]},
        index=["RECURRENT", "TWO_ONLY", "ONE_ONLY"],
    )

    filtered = filter_by_recurrence(matrix, min_patients=3)

    assert list(filtered.index) == ["RECURRENT"]


def test_filter_by_recurrence_with_threshold_two_keeps_two_genes() -> None:
    """A lower threshold must retain more genes, exercising a different cutoff."""
    matrix = pd.DataFrame(
        {"case-a": [1, 1, 1], "case-b": [1, 1, 0], "case-c": [1, 0, 0]},
        index=["RECURRENT", "TWO_ONLY", "ONE_ONLY"],
    )

    filtered = filter_by_recurrence(matrix, min_patients=2)

    assert set(filtered.index) == {"RECURRENT", "TWO_ONLY"}
