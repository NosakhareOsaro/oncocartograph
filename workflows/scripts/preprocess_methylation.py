#!/usr/bin/env python3
"""Build the filtered, M-value, top-variable-probe methylation view.

Usage:
    python workflows/scripts/preprocess_methylation.py --config workflows/config.yaml
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import yaml

from oncocartograph.preprocessing.methylation import (
    build_methylation_matrix,
    filter_by_missingness,
    select_top_variable_probes,
)


def main() -> None:
    """Build and write the methylation view as Parquet."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("workflows/config.yaml"))
    args = parser.parse_args()

    config = yaml.safe_load(args.config.read_text())
    data_dir = Path(config["data_dir"])
    processed_dir = data_dir / "processed"

    manifest = json.loads((processed_dir / "resolved_paths_methylation.json").read_text())
    resolved_files = {case_id: Path(path) for case_id, path in manifest.items()}

    matrix = build_methylation_matrix(resolved_files)
    filtered = filter_by_missingness(
        matrix, max_missing_fraction=config["preprocessing"]["methylation_max_missing_fraction"]
    )
    view = select_top_variable_probes(
        filtered, n=config["preprocessing"]["methylation_top_n_probes"]
    )

    view.to_parquet(processed_dir / "methylation_view.parquet")
    print(f"methylation view: {view.shape}")


if __name__ == "__main__":
    main()
