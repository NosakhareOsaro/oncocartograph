#!/usr/bin/env python3
"""Run the pre-registered GSE96058 replication analysis and Burstein plausibility check.

See docs/adr/0009-external-replication-methodology-and-result.md for the
full methodology and the real result this reproduces.

Usage:
    python workflows/scripts/validate_external.py --config workflows/config.yaml
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml

from oncocartograph.drug_targets.open_targets_client import OpenTargetsClient
from oncocartograph.scoring.survival import fit_univariate_cox
from oncocartograph.validation.burstein_check import check_known_biology_markers
from oncocartograph.validation.gse96058_clinical import read_combined_clinical
from oncocartograph.validation.gse96058_cohort import (
    build_gse96058_cohort_audit,
    select_gse96058_tnbc_cohort,
)
from oncocartograph.validation.gse96058_expression import read_selected_gene_expression
from oncocartograph.validation.replication import (
    build_replication_table,
    replication_table_to_frame,
    run_direction_concordance_test,
)


def main() -> None:
    """Run the GSE96058 replication analysis end-to-end and write the results."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("workflows/config.yaml"))
    args = parser.parse_args()

    config = yaml.safe_load(args.config.read_text())
    data_dir = Path(config["data_dir"])
    processed_dir = data_dir / "processed"
    external_dir = data_dir / "external" / "gse96058"

    # --- TCGA rna_seq candidates: all fittable, per the confirmed judgment call ---
    candidates = pd.read_csv(processed_dir / "candidate_evidence.csv")
    rna_candidates = candidates[candidates["view"] == "rna_seq"]
    tcga_hazard_ratios = dict(
        zip(rna_candidates["candidate_id"], rna_candidates["hazard_ratio"], strict=True)
    )

    unversioned_to_candidate: dict[str, str] = {
        candidate_id.split(":", 1)[1].split(".")[0]: candidate_id
        for candidate_id in tcga_hazard_ratios
    }
    ot_client = OpenTargetsClient(base_url=config["open_targets_api_base_url"])
    targets = ot_client.fetch_targets(list(unversioned_to_candidate))
    gene_symbols: dict[str, str | None] = {
        candidate_id: None for candidate_id in tcga_hazard_ratios
    }
    for unversioned, candidate_id in unversioned_to_candidate.items():
        target = targets.get(unversioned)
        gene_symbols[candidate_id] = target["approved_symbol"] if target else None

    # --- GSE96058 TNBC cohort + survival ---
    clinical_paths = [external_dir / Path(url).name for url in config["gse96058"]["clinical_urls"]]
    clinical = read_combined_clinical(clinical_paths)
    audit = build_gse96058_cohort_audit(clinical)
    tnbc_cohort = select_gse96058_tnbc_cohort(audit)

    gse_duration = pd.to_numeric(tnbc_cohort["overall survival days"], errors="coerce")
    gse_event = pd.to_numeric(tnbc_cohort["overall survival event"], errors="coerce")
    gse_duration.index = tnbc_cohort["title"].to_numpy()
    gse_event.index = tnbc_cohort["title"].to_numpy()
    print(f"GSE96058 TNBC cohort: {len(tnbc_cohort)} samples, {int(gse_event.sum())} events")

    # --- Expression: stream only the genes we need ---
    burstein_symbols = config["validation"]["burstein_marker_genes"]
    wanted_symbols = {s for s in gene_symbols.values() if s is not None} | set(burstein_symbols)
    expression_path = external_dir / Path(config["gse96058"]["expression_url"]).name
    expression = read_selected_gene_expression(expression_path, wanted_symbols)

    common_samples = expression.columns.intersection(gse_duration.dropna().index).intersection(
        gse_event.dropna().index
    )

    gse96058_evidence = {
        gene: fit_univariate_cox(
            expression.loc[gene].loc[common_samples],
            gse_duration.loc[common_samples],
            gse_event.loc[common_samples],
        )
        for gene in expression.index
    }

    # --- Pre-registered direction-concordance test ---
    replications = build_replication_table(tcga_hazard_ratios, gene_symbols, gse96058_evidence)
    result = run_direction_concordance_test(
        replications, alpha=config["validation"]["concordance_alpha"]
    )

    print("\n=== PRE-REGISTERED DIRECTION-CONCORDANCE RESULT ===")
    print(f"n_total_candidates = {result.n_total_candidates}")
    print(f"n_fittable         = {result.n_fittable}")
    print(f"n_concordant       = {result.n_concordant}")
    print(f"concordance_rate   = {result.concordance_rate:.4f}")
    print(f"p_value            = {result.p_value:.6f}")
    print(f"SUCCESS (pre-registered) = {result.success}")

    replication_table_to_frame(replications).to_csv(
        processed_dir / "gse96058_replication_table.csv", index=False
    )

    # --- Burstein plausibility check ---
    burstein_results = check_known_biology_markers(gse96058_evidence)
    print("\n=== BURSTEIN ET AL. 2015 KNOWN-BIOLOGY PLAUSIBILITY CHECK ===")
    for r in burstein_results:
        print(
            f"{r.gene_symbol:8s} expected={r.expected_direction:11s} "
            f"observed_hr={r.observed_hazard_ratio} plausible={r.plausible}"
        )


if __name__ == "__main__":
    main()
