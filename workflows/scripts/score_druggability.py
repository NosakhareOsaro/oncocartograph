#!/usr/bin/env python3
"""Attach real Open Targets/ChEMBL druggability evidence and re-score candidates.

Methylation candidates are excluded from druggability lookup (CpG probe
IDs are not gene identifiers -- see docs/adr/0008 and docs/methods.md
§6.1); they retain druggability=None, renormalized away in the composite
score, same as before this step.

Usage:
    python workflows/scripts/score_druggability.py --config workflows/config.yaml
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml

from oncocartograph.drug_targets.chembl_client import ChEMBLClient
from oncocartograph.drug_targets.druggability import build_druggability_evidence
from oncocartograph.drug_targets.open_targets_client import OpenTargetsClient
from oncocartograph.scoring.composite import composite_biomarker_score
from oncocartograph.scoring.evidence import (
    BiomarkerEvidence,
    IntegrationEvidence,
    RecurrenceEvidence,
)
from oncocartograph.scoring.survival import SurvivalEvidence

_DRUGGABLE_VIEWS = ("rna_seq", "copy_number")


def _row_to_evidence(row: pd.Series) -> BiomarkerEvidence:
    """Reconstruct a BiomarkerEvidence from one candidate_evidence.csv row."""
    survival = SurvivalEvidence(
        hazard_ratio=row["hazard_ratio"],
        hazard_ratio_ci_low=row["hazard_ratio_ci_low"],
        hazard_ratio_ci_high=row["hazard_ratio_ci_high"],
        p_value=row["p_value"],
        p_value_adjusted=row["p_value_adjusted"],
        n_samples=int(row["n_samples"]),
        n_events=int(row["n_events"]),
    )
    integration = None
    recurrence = None
    if row["pathway"] == "mofa":
        integration = IntegrationEvidence(
            factor=row["factor"],
            weight=row["weight"],
            view_variance_explained=row["view_variance_explained"],
        )
    else:
        recurrence = RecurrenceEvidence(
            n_patients_mutated=int(row["n_patients_mutated"]), cohort_size=int(row["cohort_size"])
        )
    return BiomarkerEvidence(
        candidate_id=row["candidate_id"],
        survival=survival,
        integration=integration,
        recurrence=recurrence,
    )


def main() -> None:
    """Re-score all candidates with real druggability evidence and write the final table."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("workflows/config.yaml"))
    args = parser.parse_args()

    config = yaml.safe_load(args.config.read_text())
    processed_dir = Path(config["data_dir"]) / "processed"

    table = pd.read_csv(processed_dir / "candidate_evidence.csv")
    evidence = {row["candidate_id"]: _row_to_evidence(row) for _, row in table.iterrows()}

    ot_client = OpenTargetsClient(base_url=config["open_targets_api_base_url"])
    chembl_client = ChEMBLClient(base_url=config["chembl_api_base_url"])

    ensembl_id_by_candidate: dict[str, str] = {}
    for candidate_id in evidence:
        view, feature = candidate_id.split(":", 1)
        if view in _DRUGGABLE_VIEWS:
            ensembl_id_by_candidate[candidate_id] = feature.split(".")[0]

    mutation_symbols = [cid.split(":", 1)[1] for cid in evidence if cid.startswith("mutation:")]
    symbol_to_ensembl = ot_client.map_symbols_to_ensembl_ids(mutation_symbols)

    all_ensembl_ids = set(ensembl_id_by_candidate.values()) | {
        v for v in symbol_to_ensembl.values() if v is not None
    }
    druggability_by_ensembl = build_druggability_evidence(
        list(all_ensembl_ids), open_targets_client=ot_client, chembl_client=chembl_client
    )

    rows = []
    for candidate_id, ev in evidence.items():
        # _row_to_evidence always sets `survival`; the dataclass's Optional
        # type reflects a general schema, not this script's invariant.
        assert ev.survival is not None
        if candidate_id in ensembl_id_by_candidate:
            ensembl_id: str | None = ensembl_id_by_candidate[candidate_id]
        elif candidate_id.startswith("mutation:"):
            ensembl_id = symbol_to_ensembl.get(candidate_id.split(":", 1)[1])
        else:
            ensembl_id = None
        drug = druggability_by_ensembl.get(ensembl_id) if ensembl_id else None

        ev_with_drug = BiomarkerEvidence(
            candidate_id=ev.candidate_id,
            survival=ev.survival,
            integration=ev.integration,
            recurrence=ev.recurrence,
            druggability=drug,
        )
        rows.append(
            {
                "candidate_id": candidate_id,
                "composite_score": composite_biomarker_score(ev_with_drug),
                "hazard_ratio": ev.survival.hazard_ratio,
                "p_value_adjusted": ev.survival.p_value_adjusted,
                "tractability_score": drug.tractability_score if drug else None,
                "chembl_max_phase": drug.chembl_max_phase if drug else None,
            }
        )

    result = pd.DataFrame(rows).sort_values("composite_score", ascending=False)
    result.to_csv(processed_dir / "composite_scores.csv", index=False)

    n_with_drug = result["tractability_score"].notna().sum()
    print(f"Druggability evidence populated for {n_with_drug}/{len(result)} candidates")


if __name__ == "__main__":
    main()
