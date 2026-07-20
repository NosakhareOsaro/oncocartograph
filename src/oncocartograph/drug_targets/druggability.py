"""Orchestrates Open Targets + ChEMBL lookups into DruggabilityEvidence.

Combines Open Targets tractability evidence and ChEMBL clinical trial
phase evidence into the single ``tractability_score`` the composite
scoring formula actually reads (see
``oncocartograph.scoring.composite._druggability_score``), via
``max(bucket_tier_score, chembl_max_phase / 4)`` -- either signal alone
can indicate a real drug exists or is in development, so this project
takes the stronger of the two rather than averaging them down.
``chembl_max_phase`` is retained on the evidence object separately for
transparent reporting (e.g. manuscript tables), even though it is not a
second independent term in the composite score.
"""

from __future__ import annotations

from collections.abc import Sequence

from oncocartograph.drug_targets.chembl_client import ChEMBLClient
from oncocartograph.drug_targets.open_targets_client import OpenTargetsClient
from oncocartograph.drug_targets.tractability import score_tractability
from oncocartograph.scoring.evidence import DruggabilityEvidence

#: ChEMBL's max_phase scale is 0-4; dividing by this puts it on the same
#: [0, 1] scale as the tractability tier score.
_CHEMBL_MAX_PHASE_CEILING = 4.0


def build_druggability_evidence(
    ensembl_ids: Sequence[str],
    *,
    open_targets_client: OpenTargetsClient,
    chembl_client: ChEMBLClient,
) -> dict[str, DruggabilityEvidence]:
    """Build DruggabilityEvidence for a batch of genes from live Open Targets + ChEMBL data.

    Args:
        ensembl_ids: Unversioned Ensembl gene IDs to look up (e.g.
            ``["ENSG00000141510"]``).
        open_targets_client: A configured :class:`OpenTargetsClient`.
        chembl_client: A configured :class:`ChEMBLClient`.

    Returns:
        A dict mapping each Ensembl ID Open Targets has a record for to
        its :class:`DruggabilityEvidence`. IDs Open Targets has no record
        for are simply absent, not mapped to a zero/default evidence
        object.
    """
    targets = open_targets_client.fetch_targets(ensembl_ids)

    accessions = [
        info["uniprot_accession"] for info in targets.values() if info["uniprot_accession"]
    ]
    target_ids_by_accession = chembl_client.resolve_accessions_to_target_ids(accessions)

    chembl_target_ids = [
        target_id for target_id in target_ids_by_accession.values() if target_id is not None
    ]
    max_phase_by_target_id = chembl_client.fetch_max_phase(chembl_target_ids)

    evidence: dict[str, DruggabilityEvidence] = {}
    for ensembl_id, info in targets.items():
        tractability_score = score_tractability(info["tractability"])

        accession = info["uniprot_accession"]
        chembl_target_id = target_ids_by_accession.get(accession) if accession else None
        max_phase = max_phase_by_target_id.get(chembl_target_id) if chembl_target_id else None

        combined_score = tractability_score
        if max_phase is not None:
            combined_score = max(tractability_score, max_phase / _CHEMBL_MAX_PHASE_CEILING)

        evidence[ensembl_id] = DruggabilityEvidence(
            tractability_score=combined_score,
            chembl_max_phase=max_phase,
        )
    return evidence
