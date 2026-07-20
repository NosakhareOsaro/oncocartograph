#!/usr/bin/env python3
"""Build the filtered, VST-normalized, top-variable-gene RNA-seq view.

Usage:
    python workflows/scripts/preprocess_rna_seq.py --config workflows/config.yaml
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import yaml

from oncocartograph.preprocessing.rna_seq import (
    build_counts_matrix,
    filter_low_expression,
    normalize_and_vst,
    select_top_variable_genes,
)


def main() -> None:
    """Build and write the RNA-seq view as Parquet."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("workflows/config.yaml"))
    args = parser.parse_args()

    config = yaml.safe_load(args.config.read_text())
    data_dir = Path(config["data_dir"])
    processed_dir = data_dir / "processed"

    manifest = json.loads((processed_dir / "resolved_paths_rna_seq.json").read_text())
    resolved_files = {case_id: Path(path) for case_id, path in manifest.items()}

    raw_counts = build_counts_matrix(resolved_files)
    filtered = filter_low_expression(
        raw_counts,
        min_cpm=config["preprocessing"]["rna_seq_min_cpm"],
        min_fraction_samples=config["preprocessing"]["rna_seq_min_fraction_samples"],
    )
    vst = normalize_and_vst(filtered)
    view = select_top_variable_genes(vst, n=config["preprocessing"]["rna_seq_top_n_genes"])

    view.to_parquet(processed_dir / "rna_seq_view.parquet")
    print(f"rna_seq view: {view.shape}")


if __name__ == "__main__":
    main()
