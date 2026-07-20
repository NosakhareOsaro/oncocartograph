"""Tests for oncocartograph.drug_targets.tractability.

The TP53 fixture is real Open Targets tractability data (captured via a
live GraphQL query against ENSG00000141510 on 2026-07-20, while planning
this work package), not synthetic -- included as a real-world regression
check alongside the synthetic edge cases below.
"""

from __future__ import annotations

from oncocartograph.drug_targets.tractability import (
    TIER_ANY_EVIDENCE,
    TIER_APPROVED,
    TIER_CLINICAL,
    TIER_NONE,
    score_tractability,
)

# Real Open Targets tractability response for TP53 (ENSG00000141510),
# captured 2026-07-20. Advanced Clinical is true (SM and OC modalities);
# Approved Drug is false everywhere.
_REAL_TP53_TRACTABILITY = [
    {"label": "Approved Drug", "modality": "SM", "value": False},
    {"label": "Advanced Clinical", "modality": "SM", "value": True},
    {"label": "Phase 1 Clinical", "modality": "SM", "value": False},
    {"label": "Structure with Ligand", "modality": "SM", "value": True},
    {"label": "High-Quality Ligand", "modality": "SM", "value": True},
    {"label": "High-Quality Pocket", "modality": "SM", "value": False},
    {"label": "Med-Quality Pocket", "modality": "SM", "value": True},
    {"label": "Druggable Family", "modality": "SM", "value": True},
    {"label": "Approved Drug", "modality": "AB", "value": False},
    {"label": "Advanced Clinical", "modality": "AB", "value": False},
    {"label": "Phase 1 Clinical", "modality": "AB", "value": False},
    {"label": "Approved Drug", "modality": "PR", "value": False},
    {"label": "UniProt Ubiquitination", "modality": "PR", "value": True},
    {"label": "Database Ubiquitination", "modality": "PR", "value": True},
    {"label": "Small Molecule Binder", "modality": "PR", "value": True},
    {"label": "Approved Drug", "modality": "OC", "value": False},
    {"label": "Advanced Clinical", "modality": "OC", "value": True},
    {"label": "Phase 1 Clinical", "modality": "OC", "value": False},
]

# Real Open Targets tractability response for GTF3C1 (ENSG00000077235),
# captured 2026-07-20: only weak structural/biological evidence, no
# clinical-stage or approved-drug evidence in any modality.
_REAL_GTF3C1_TRACTABILITY = [
    {"label": "Approved Drug", "modality": "SM", "value": False},
    {"label": "Advanced Clinical", "modality": "SM", "value": False},
    {"label": "Phase 1 Clinical", "modality": "SM", "value": False},
    {"label": "Structure with Ligand", "modality": "SM", "value": False},
    {"label": "High-Quality Ligand", "modality": "SM", "value": False},
    {"label": "High-Quality Pocket", "modality": "SM", "value": False},
    {"label": "Med-Quality Pocket", "modality": "SM", "value": False},
    {"label": "Druggable Family", "modality": "SM", "value": False},
    {"label": "Approved Drug", "modality": "AB", "value": False},
    {"label": "Advanced Clinical", "modality": "AB", "value": False},
    {"label": "Phase 1 Clinical", "modality": "AB", "value": False},
    {"label": "UniProt loc high conf", "modality": "AB", "value": True},
    {"label": "Approved Drug", "modality": "PR", "value": False},
    {"label": "Small Molecule Binder", "modality": "PR", "value": True},
    {"label": "Approved Drug", "modality": "OC", "value": False},
    {"label": "Advanced Clinical", "modality": "OC", "value": False},
    {"label": "Phase 1 Clinical", "modality": "OC", "value": False},
]


def test_score_tractability_on_real_tp53_data_is_clinical_tier() -> None:
    """TP53 has Advanced Clinical evidence (SM, OC) but no Approved Drug -> 0.66."""
    assert score_tractability(_REAL_TP53_TRACTABILITY) == TIER_CLINICAL


def test_score_tractability_on_real_gtf3c1_data_is_any_evidence_tier() -> None:
    """GTF3C1 has only weak structural evidence, no clinical-stage evidence -> 0.33."""
    assert score_tractability(_REAL_GTF3C1_TRACTABILITY) == TIER_ANY_EVIDENCE


def test_score_tractability_approved_drug_outranks_clinical() -> None:
    """A target with both Approved Drug and Advanced Clinical true must score as Approved."""
    buckets = [
        {"label": "Approved Drug", "modality": "SM", "value": True},
        {"label": "Advanced Clinical", "modality": "SM", "value": True},
    ]
    assert score_tractability(buckets) == TIER_APPROVED


def test_score_tractability_no_evidence_at_all() -> None:
    """A target with every bucket false must score 0."""
    buckets = [
        {"label": "Approved Drug", "modality": "SM", "value": False},
        {"label": "Druggable Family", "modality": "SM", "value": False},
    ]
    assert score_tractability(buckets) == TIER_NONE


def test_score_tractability_empty_buckets_is_none_tier() -> None:
    """No tractability data at all (empty list) must score 0, not error."""
    assert score_tractability([]) == TIER_NONE


def test_score_tractability_phase_1_clinical_alone_is_clinical_tier() -> None:
    """Phase 1 Clinical (not just Advanced Clinical) must also count as clinical-stage."""
    buckets = [{"label": "Phase 1 Clinical", "modality": "AB", "value": True}]
    assert score_tractability(buckets) == TIER_CLINICAL
