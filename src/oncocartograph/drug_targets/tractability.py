"""Tractability scoring: collapses Open Targets' tractability buckets into [0, 1].

Open Targets reports ~28 boolean buckets per target across four modalities
(SM=Small Molecule, AB=Antibody, PR=PROTAC, OC=Other Clinical), confirmed
via a live query against real targets (TP53, GTF3C1) before this module
was written. The clinical-stage buckets ("Approved Drug", "Advanced
Clinical", "Phase 1 Clinical") use identical label text regardless of
modality, so tiering only needs to inspect the label, not a per-modality
lookup table.

This project's chosen collapse is a 3-tier scheme: clinical-stage evidence
(a real drug exists or is being tested) outranks any purely structural or
biological tractability evidence, which in turn outranks no evidence at
all.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

#: Bucket label indicating an approved drug already exists for this target.
APPROVED_DRUG_LABEL = "Approved Drug"

#: Bucket labels indicating a drug is in clinical development (any modality).
CLINICAL_STAGE_LABELS = frozenset({"Advanced Clinical", "Phase 1 Clinical"})

#: Tier scores, highest to lowest.
TIER_APPROVED = 1.0
TIER_CLINICAL = 0.66
TIER_ANY_EVIDENCE = 0.33
TIER_NONE = 0.0


def score_tractability(buckets: Sequence[Mapping[str, object]]) -> float:
    """Collapse a target's Open Targets tractability buckets into one [0, 1] score.

    Args:
        buckets: The raw ``tractability`` list from an Open Targets
            ``target``/``targets`` GraphQL response: each item a mapping
            with ``label`` (str), ``modality`` (str), and ``value``
            (bool) keys.

    Returns:
        1.0 if any modality has an approved drug; 0.66 if any modality
        has a drug in clinical development (Phase 1 or more advanced);
        0.33 if any other tractability bucket is true; 0.0 if no bucket
        is true at all (including an empty ``buckets`` sequence).
    """
    if any(bucket["value"] and bucket["label"] == APPROVED_DRUG_LABEL for bucket in buckets):
        return TIER_APPROVED
    if any(bucket["value"] and bucket["label"] in CLINICAL_STAGE_LABELS for bucket in buckets):
        return TIER_CLINICAL
    if any(bucket["value"] for bucket in buckets):
        return TIER_ANY_EVIDENCE
    return TIER_NONE
