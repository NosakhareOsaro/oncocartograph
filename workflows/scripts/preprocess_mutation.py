#!/usr/bin/env python3
"""Build the binary, recurrence-filtered mutation view.

Usage:
    python workflows/scripts/preprocess_mutation.py --config workflows/config.yaml
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import yaml

from oncocartograph.preprocessing.mutation import build_mutation_matrix, filter_by_recurrence


def main() -> None:
    """Build and write the mutation view as Parquet."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("workflows/config.yaml"))
    args = parser.parse_args()

    config = yaml.safe_load(args.config.read_text())
    data_dir = Path(config["data_dir"])
    processed_dir = data_dir / "processed"

    manifest = json.loads((processed_dir / "resolved_paths_mutation.json").read_text())
    resolved_files = {case_id: Path(path) for case_id, path in manifest.items()}

    matrix = build_mutation_matrix(resolved_files)
    view = filter_by_recurrence(
        matrix, min_patients=config["preprocessing"]["mutation_min_recurrent_patients"]
    )

    view.to_parquet(processed_dir / "mutation_view.parquet")
    print(f"mutation view: {view.shape}")


if __name__ == "__main__":
    main()
