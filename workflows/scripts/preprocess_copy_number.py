#!/usr/bin/env python3
"""Build the relative-to-diploid log2, top-variable-gene copy number view.

Usage:
    python workflows/scripts/preprocess_copy_number.py --config workflows/config.yaml
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import yaml

from oncocartograph.preprocessing.copy_number import (
    build_copy_number_matrix,
    select_top_variable_genes,
)


def main() -> None:
    """Build and write the copy number view as Parquet."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("workflows/config.yaml"))
    args = parser.parse_args()

    config = yaml.safe_load(args.config.read_text())
    data_dir = Path(config["data_dir"])
    processed_dir = data_dir / "processed"

    manifest = json.loads((processed_dir / "resolved_paths_copy_number.json").read_text())
    resolved_files = {case_id: Path(path) for case_id, path in manifest.items()}

    matrix = build_copy_number_matrix(resolved_files)
    view = select_top_variable_genes(matrix, n=config["preprocessing"]["copy_number_top_n_genes"])

    view.to_parquet(processed_dir / "copy_number_view.parquet")
    print(f"copy_number view: {view.shape}")


if __name__ == "__main__":
    main()
