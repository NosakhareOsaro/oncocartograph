"""Tests for oncocartograph.drug_targets.druggability.

Uses stub clients (not HTTP mocking -- the individual clients are already
unit tested elsewhere) to isolate the orchestration/combination logic:
how Open Targets tractability and ChEMBL max phase combine into one
DruggabilityEvidence per gene.
"""

from __future__ import annotations

from typing import Any

from oncocartograph.drug_targets.druggability import build_druggability_evidence


class _StubOpenTargetsClient:
    def __init__(self, targets: dict[str, dict[str, Any]]) -> None:
        self._targets = targets

    def fetch_targets(self, ensembl_ids: list[str]) -> dict[str, dict[str, Any]]:
        return {eid: self._targets[eid] for eid in ensembl_ids if eid in self._targets}


class _StubChEMBLClient:
    def __init__(
        self,
        target_ids_by_accession: dict[str, str | None],
        max_phase_by_target_id: dict[str, float | None],
    ) -> None:
        self._target_ids_by_accession = target_ids_by_accession
        self._max_phase_by_target_id = max_phase_by_target_id

    def resolve_accessions_to_target_ids(self, accessions: list[str]) -> dict[str, str | None]:
        return {a: self._target_ids_by_accession.get(a) for a in accessions}

    def fetch_max_phase(self, target_chembl_ids: list[str]) -> dict[str, float | None]:
        return {t: self._max_phase_by_target_id.get(t) for t in target_chembl_ids}


def _bucket(label: str, value: bool, modality: str = "SM") -> dict[str, Any]:
    return {"label": label, "modality": modality, "value": value}


def test_build_druggability_evidence_combines_both_signals_via_max() -> None:
    """When ChEMBL's max_phase/4 exceeds the tractability tier score, the combined score
    must use the ChEMBL-derived value."""
    ot_client = _StubOpenTargetsClient(
        {
            "ENSG1": {
                "approved_symbol": "GENE1",
                "tractability": [_bucket("Druggable Family", True)],  # tier = 0.33
                "uniprot_accession": "P00001",
            }
        }
    )
    chembl_client = _StubChEMBLClient(
        target_ids_by_accession={"P00001": "CHEMBL1"},
        max_phase_by_target_id={"CHEMBL1": 4.0},  # 4/4 = 1.0, exceeds 0.33
    )

    result = build_druggability_evidence(
        ["ENSG1"], open_targets_client=ot_client, chembl_client=chembl_client
    )

    assert result["ENSG1"].tractability_score == 1.0
    assert result["ENSG1"].chembl_max_phase == 4.0


def test_build_druggability_evidence_keeps_stronger_tractability_signal() -> None:
    """When the OT tractability tier score exceeds max_phase/4, the combined score
    must keep the tractability-derived value, not be dragged down by a weak ChEMBL phase."""
    ot_client = _StubOpenTargetsClient(
        {
            "ENSG1": {
                "approved_symbol": "GENE1",
                "tractability": [_bucket("Approved Drug", True)],  # tier = 1.0
                "uniprot_accession": "P00001",
            }
        }
    )
    chembl_client = _StubChEMBLClient(
        target_ids_by_accession={"P00001": "CHEMBL1"},
        max_phase_by_target_id={"CHEMBL1": 1.0},  # 1/4 = 0.25, weaker than 1.0
    )

    result = build_druggability_evidence(
        ["ENSG1"], open_targets_client=ot_client, chembl_client=chembl_client
    )

    assert result["ENSG1"].tractability_score == 1.0
    assert result["ENSG1"].chembl_max_phase == 1.0


def test_build_druggability_evidence_handles_missing_uniprot_accession() -> None:
    """A gene with no UniProt accession must still get tractability-only evidence,
    with chembl_max_phase left None, not crash on the missing lookup key."""
    ot_client = _StubOpenTargetsClient(
        {
            "ENSG1": {
                "approved_symbol": "GENE1",
                "tractability": [_bucket("Druggable Family", True)],
                "uniprot_accession": None,
            }
        }
    )
    chembl_client = _StubChEMBLClient(target_ids_by_accession={}, max_phase_by_target_id={})

    result = build_druggability_evidence(
        ["ENSG1"], open_targets_client=ot_client, chembl_client=chembl_client
    )

    assert result["ENSG1"].tractability_score == 0.33
    assert result["ENSG1"].chembl_max_phase is None


def test_build_druggability_evidence_handles_unresolvable_chembl_target() -> None:
    """A UniProt accession with no matching ChEMBL target must leave chembl_max_phase None."""
    ot_client = _StubOpenTargetsClient(
        {
            "ENSG1": {
                "approved_symbol": "GENE1",
                "tractability": [_bucket("Druggable Family", True)],
                "uniprot_accession": "P00001",
            }
        }
    )
    chembl_client = _StubChEMBLClient(
        target_ids_by_accession={"P00001": None}, max_phase_by_target_id={}
    )

    result = build_druggability_evidence(
        ["ENSG1"], open_targets_client=ot_client, chembl_client=chembl_client
    )

    assert result["ENSG1"].chembl_max_phase is None
    assert result["ENSG1"].tractability_score == 0.33


def test_build_druggability_evidence_omits_genes_not_found_in_open_targets() -> None:
    """An Ensembl ID Open Targets has no record for must be absent, not a default entry."""
    ot_client = _StubOpenTargetsClient({})
    chembl_client = _StubChEMBLClient(target_ids_by_accession={}, max_phase_by_target_id={})

    result = build_druggability_evidence(
        ["ENSG_UNKNOWN"], open_targets_client=ot_client, chembl_client=chembl_client
    )

    assert result == {}


def test_build_druggability_evidence_empty_input() -> None:
    """No genes to look up must return {} without error."""
    ot_client = _StubOpenTargetsClient({})
    chembl_client = _StubChEMBLClient(target_ids_by_accession={}, max_phase_by_target_id={})

    result = build_druggability_evidence(
        [], open_targets_client=ot_client, chembl_client=chembl_client
    )

    assert result == {}
