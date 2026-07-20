#!/usr/bin/env python3
"""Screen MOFA+-derived and mutation-recurrence-derived candidates for survival association.

Two candidate-generation pathways feed one composite score (docs/methods.md
§4.1): MOFA+ factor loading (top 20 features per factor, per view, across
every factor clearing the variance-explained screening threshold) for
RNA-seq/methylation/copy-number, and mutation-recurrence (every gene
surviving the recurrence filter) for mutations -- bypassing MOFA+
entirely, since the mutation view contributes essentially no variance to
any factor (docs/methods.md §3.4).

This does not yet include druggability evidence; see score_druggability.py.

Usage:
    python workflows/scripts/score_candidates.py --config workflows/config.yaml
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml
from mofax.core import mofa_model

from oncocartograph.integration.mofa import get_variance_explained, load_mofa_model
from oncocartograph.scoring.composite import composite_biomarker_score
from oncocartograph.scoring.evidence import (
    BiomarkerEvidence,
    IntegrationEvidence,
    RecurrenceEvidence,
)
from oncocartograph.scoring.survival import SurvivalEvidence, screen_survival_associations

_TOP_N_FEATURES_PER_FACTOR = 20


def _mofa_derived_evidence(
    model: mofa_model,
    views: dict[str, pd.DataFrame],
    duration: pd.Series,
    event: pd.Series,
    variance_threshold: float,
) -> dict[str, BiomarkerEvidence]:
    """Screen MOFA+-derived candidates (RNA-seq/methylation/copy-number) for survival."""
    variance_explained = get_variance_explained(model)
    pivot = variance_explained.pivot(index="Factor", columns="View", values="R2")
    informative_factors = pivot[pivot.max(axis=1) >= variance_threshold].index.tolist()

    top_features = model.get_top_features(
        factors=informative_factors, n_features=_TOP_N_FEATURES_PER_FACTOR, df=True, per_view=True
    )
    non_mutation = top_features[top_features["view"] != "mutation"]
    deduped = non_mutation.loc[non_mutation.groupby(["feature", "view"])["value_abs"].idxmax()]

    evidence: dict[str, BiomarkerEvidence] = {}
    for view_name, view_matrix in views.items():
        view_candidates = deduped[deduped["view"] == view_name]
        features = [f for f in view_candidates["feature"] if f in view_matrix.index]
        if not features:
            continue
        sub_matrix = view_matrix.loc[features]
        common = sub_matrix.columns.intersection(duration.index)
        sub_matrix = sub_matrix[common]
        screen = screen_survival_associations(sub_matrix, duration.loc[common], event.loc[common])

        for feature in screen.index:
            row = view_candidates[view_candidates["feature"] == feature].iloc[0]
            factor = row["factor"]
            variance = float(pivot.loc[factor, view_name]) if factor in pivot.index else 0.0
            surv_row = screen.loc[feature]
            candidate_id = f"{view_name}:{feature}"
            evidence[candidate_id] = BiomarkerEvidence(
                candidate_id=candidate_id,
                survival=SurvivalEvidence(
                    hazard_ratio=surv_row["hazard_ratio"],
                    hazard_ratio_ci_low=surv_row["hazard_ratio_ci_low"],
                    hazard_ratio_ci_high=surv_row["hazard_ratio_ci_high"],
                    p_value=surv_row["p_value"],
                    p_value_adjusted=surv_row["p_adj"],
                    n_samples=int(surv_row["n_samples"]),
                    n_events=int(surv_row["n_events"]),
                ),
                integration=IntegrationEvidence(
                    factor=factor, weight=float(row["value"]), view_variance_explained=variance
                ),
            )
    return evidence


def _mutation_derived_evidence(
    mutation_view: pd.DataFrame, duration: pd.Series, event: pd.Series
) -> dict[str, BiomarkerEvidence]:
    """Screen mutation-recurrence-derived candidates for survival association."""
    common = mutation_view.columns.intersection(duration.index)
    sub_matrix = mutation_view[common]
    screen = screen_survival_associations(sub_matrix, duration.loc[common], event.loc[common])

    cohort_size = sub_matrix.shape[1]
    evidence: dict[str, BiomarkerEvidence] = {}
    for gene in screen.index:
        row = screen.loc[gene]
        candidate_id = f"mutation:{gene}"
        evidence[candidate_id] = BiomarkerEvidence(
            candidate_id=candidate_id,
            survival=SurvivalEvidence(
                hazard_ratio=row["hazard_ratio"],
                hazard_ratio_ci_low=row["hazard_ratio_ci_low"],
                hazard_ratio_ci_high=row["hazard_ratio_ci_high"],
                p_value=row["p_value"],
                p_value_adjusted=row["p_adj"],
                n_samples=int(row["n_samples"]),
                n_events=int(row["n_events"]),
            ),
            recurrence=RecurrenceEvidence(
                n_patients_mutated=int(sub_matrix.loc[gene].sum()), cohort_size=cohort_size
            ),
        )
    return evidence


def _evidence_to_frame(evidence: dict[str, BiomarkerEvidence]) -> pd.DataFrame:
    """Flatten evidence + composite score into one exportable table."""
    rows = []
    for candidate_id, ev in evidence.items():
        # Every entry built above always sets `survival` (both pathways
        # screen for it before an evidence record is even created); the
        # dataclass's Optional type reflects a general schema, not this
        # script's actual invariant.
        assert ev.survival is not None
        view = candidate_id.split(":", 1)[0]
        rows.append(
            {
                "candidate_id": candidate_id,
                "view": view,
                "pathway": "mofa" if ev.integration is not None else "mutation_recurrence",
                "hazard_ratio": ev.survival.hazard_ratio,
                "hazard_ratio_ci_low": ev.survival.hazard_ratio_ci_low,
                "hazard_ratio_ci_high": ev.survival.hazard_ratio_ci_high,
                "p_value": ev.survival.p_value,
                "p_value_adjusted": ev.survival.p_value_adjusted,
                "n_samples": ev.survival.n_samples,
                "n_events": ev.survival.n_events,
                "factor": ev.integration.factor if ev.integration else None,
                "weight": ev.integration.weight if ev.integration else None,
                "view_variance_explained": (
                    ev.integration.view_variance_explained if ev.integration else None
                ),
                "n_patients_mutated": ev.recurrence.n_patients_mutated if ev.recurrence else None,
                "cohort_size": ev.recurrence.cohort_size if ev.recurrence else None,
                "composite_score_no_druggability": composite_biomarker_score(ev),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    """Screen all candidates and write the pre-druggability evidence table."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("workflows/config.yaml"))
    args = parser.parse_args()

    config = yaml.safe_load(args.config.read_text())
    data_dir = Path(config["data_dir"])
    processed_dir = data_dir / "processed"

    survival = pd.read_csv(processed_dir / "survival.csv", index_col=0)
    duration, event = survival["duration"], survival["event"]

    views = {
        name: pd.read_parquet(processed_dir / f"{name}_view.parquet")
        for name in ("rna_seq", "methylation", "copy_number")
    }
    mutation_view = pd.read_parquet(processed_dir / "mutation_view.parquet")

    model = load_mofa_model(processed_dir / "mofa_model.hdf5")
    mofa_evidence = _mofa_derived_evidence(
        model, views, duration, event, config["mofa"]["variance_explained_threshold_percent"]
    )
    mutation_evidence = _mutation_derived_evidence(mutation_view, duration, event)

    all_evidence = {**mofa_evidence, **mutation_evidence}
    table = _evidence_to_frame(all_evidence)
    table.to_csv(processed_dir / "candidate_evidence.csv", index=False)

    alpha = config["scoring"]["fdr_alpha"]
    n_sig = int((table["p_value_adjusted"] < alpha).sum())
    print(f"Screened {len(table)} candidates; {n_sig} survived FDR correction (alpha={alpha})")


if __name__ == "__main__":
    main()
