"""Tests for oncocartograph.validation.replication."""

from __future__ import annotations

import pytest

from oncocartograph.scoring.survival import SurvivalEvidence
from oncocartograph.validation.replication import (
    CandidateReplication,
    build_replication_table,
    replication_table_to_frame,
    run_direction_concordance_test,
)


def _evidence(hazard_ratio: float) -> SurvivalEvidence:
    return SurvivalEvidence(
        hazard_ratio=hazard_ratio,
        hazard_ratio_ci_low=hazard_ratio * 0.5,
        hazard_ratio_ci_high=hazard_ratio * 1.5,
        p_value=0.1,
        n_samples=100,
        n_events=20,
    )


def test_build_replication_table_marks_concordant_when_directions_agree() -> None:
    """A TCGA HR>1 and a GSE96058 HR>1 for the same gene must be marked concordant."""
    records = build_replication_table(
        tcga_hazard_ratios={"rna_seq:ENSG1": 1.5},
        gene_symbols={"rna_seq:ENSG1": "TP53"},
        gse96058_evidence={"TP53": _evidence(2.0)},
    )
    assert records == [
        CandidateReplication(
            candidate_id="rna_seq:ENSG1",
            gene_symbol="TP53",
            tcga_hazard_ratio=1.5,
            gse96058_evidence=_evidence(2.0),
            concordant=True,
        )
    ]


def test_build_replication_table_marks_discordant_when_directions_disagree() -> None:
    """A TCGA HR>1 (harmful) and a GSE96058 HR<1 (protective) must be marked discordant."""
    records = build_replication_table(
        tcga_hazard_ratios={"rna_seq:ENSG1": 1.5},
        gene_symbols={"rna_seq:ENSG1": "TP53"},
        gse96058_evidence={"TP53": _evidence(0.5)},
    )
    assert records[0].concordant is False


def test_build_replication_table_unresolved_symbol_is_not_fittable() -> None:
    """A candidate whose Ensembl ID could not be mapped to a symbol must have concordant=None."""
    records = build_replication_table(
        tcga_hazard_ratios={"rna_seq:ENSG1": 1.5},
        gene_symbols={"rna_seq:ENSG1": None},
        gse96058_evidence={},
    )
    assert records[0].gene_symbol is None
    assert records[0].concordant is None


def test_build_replication_table_tcga_hr_of_exactly_one_is_never_concordant() -> None:
    """An exact HR of 1 (zero-signed log HR) must never count as concordant, even on both sides."""
    records = build_replication_table(
        tcga_hazard_ratios={"rna_seq:ENSG1": 1.0},
        gene_symbols={"rna_seq:ENSG1": "TP53"},
        gse96058_evidence={"TP53": _evidence(1.0)},
    )
    assert records[0].concordant is False


def test_build_replication_table_gene_absent_from_gse96058_is_not_fittable() -> None:
    """A resolved gene symbol absent from the GSE96058 evidence dict must have concordant=None."""
    records = build_replication_table(
        tcga_hazard_ratios={"rna_seq:ENSG1": 1.5},
        gene_symbols={"rna_seq:ENSG1": "TP53"},
        gse96058_evidence={},
    )
    assert records[0].concordant is None


def test_direction_concordance_uses_only_fittable_candidates_as_denominator() -> None:
    """Candidates with no GSE96058 evidence are excluded from n_fittable, not scored discordant."""
    records = build_replication_table(
        tcga_hazard_ratios={"rna_seq:ENSG1": 1.5, "rna_seq:ENSG2": 1.2},
        gene_symbols={"rna_seq:ENSG1": "TP53", "rna_seq:ENSG2": "AR"},
        gse96058_evidence={"TP53": _evidence(2.0)},
    )
    result = run_direction_concordance_test(records)
    assert result.n_total_candidates == 2
    assert result.n_fittable == 1
    assert result.n_concordant == 1


def test_direction_concordance_all_concordant_is_significant() -> None:
    """A perfect concordance run with enough candidates must clear the pre-registered alpha=0.05."""
    records = [
        CandidateReplication(f"c{i}", f"g{i}", 1.5, _evidence(2.0), concordant=True)
        for i in range(10)
    ]
    result = run_direction_concordance_test(records)
    assert result.n_fittable == 10
    assert result.concordance_rate == 1.0
    assert result.success is True
    assert result.p_value < 0.05


def test_direction_concordance_chance_level_is_not_significant() -> None:
    """A 50/50 concordance split must not clear the pre-registered alpha=0.05 bar."""
    records = [
        CandidateReplication(f"c{i}", f"g{i}", 1.5, _evidence(2.0), concordant=(i % 2 == 0))
        for i in range(20)
    ]
    result = run_direction_concordance_test(records)
    assert result.concordance_rate == pytest.approx(0.5)
    assert result.success is False


def test_replication_table_to_frame_has_expected_columns() -> None:
    """The exported DataFrame must carry both TCGA and GSE96058 hazard ratios for review."""
    records = build_replication_table(
        tcga_hazard_ratios={"rna_seq:ENSG1": 1.5},
        gene_symbols={"rna_seq:ENSG1": "TP53"},
        gse96058_evidence={"TP53": _evidence(2.0)},
    )
    frame = replication_table_to_frame(records)
    assert list(frame.columns) == [
        "candidate_id",
        "gene_symbol",
        "tcga_hazard_ratio",
        "gse96058_hazard_ratio",
        "gse96058_p_value",
        "gse96058_n_samples",
        "concordant",
    ]
    assert frame.loc[0, "gse96058_hazard_ratio"] == 2.0
