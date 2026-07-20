#!/usr/bin/env python3
"""Train the MOFA+ multi-omics model and extract factor values / variance explained.

Usage:
    python workflows/scripts/integrate_mofa.py --config workflows/config.yaml
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml

from oncocartograph.integration.mofa import (
    Likelihood,
    get_factor_values,
    get_variance_explained,
    load_mofa_model,
    train_mofa_model,
)

_LIKELIHOODS: dict[str, Likelihood] = {
    "rna_seq": "gaussian",
    "methylation": "gaussian",
    "copy_number": "gaussian",
    "mutation": "bernoulli",
}


def main() -> None:
    """Train the MOFA+ model and write factor/variance-explained tables."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("workflows/config.yaml"))
    args = parser.parse_args()

    config = yaml.safe_load(args.config.read_text())
    data_dir = Path(config["data_dir"])
    processed_dir = data_dir / "processed"

    views = {name: pd.read_parquet(processed_dir / f"{name}_view.parquet") for name in _LIKELIHOODS}

    model_path = processed_dir / "mofa_model.hdf5"
    train_mofa_model(
        views,
        _LIKELIHOODS,
        n_factors=config["mofa"]["n_factors"],
        seed=config["random_seed"],
        outfile=model_path,
        convergence_mode=config["mofa"]["convergence_mode"],
        max_iterations=config["mofa"]["max_iterations"],
    )

    model = load_mofa_model(model_path)
    get_factor_values(model).to_csv(processed_dir / "mofa_factor_values.csv")
    get_variance_explained(model).to_csv(processed_dir / "mofa_variance_explained.csv", index=False)
    print(f"Trained MOFA+ model with {config['mofa']['n_factors']} factors -> {model_path}")


if __name__ == "__main__":
    main()
