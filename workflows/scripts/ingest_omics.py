#!/usr/bin/env python3
"""Resolve and download one omic layer's Primary Tumor files for the TNBC cohort.

Queries GDC for the requested layer's files restricted to the TNBC
cohort's case UUIDs, resolves down to one Primary Tumor file per patient
(docs/methods.md §2.1; ASCAT-workflow-priority for copy number, docs/adr/0005),
downloads only the resolved subset (skipping files already present on
disk, so a re-run doesn't re-download), and writes a case_id -> local
path manifest for the corresponding preprocessing rule to consume.

Usage:
    python workflows/scripts/ingest_omics.py --layer rna_seq --config workflows/config.yaml
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from oncocartograph.data_ingestion.gdc_client import GDCClient
from oncocartograph.data_ingestion.omics_ingestion import (
    TCGA_PROJECT_ID,
    copy_number_file_filter,
    methylation_file_filter,
    mutation_file_filter,
    rna_seq_file_filter,
)
from oncocartograph.data_ingestion.provenance import record_download
from oncocartograph.preprocessing.copy_number import CNV_MANIFEST_FIELDS, resolve_copy_number_files
from oncocartograph.preprocessing.sample_manifest import (
    MANIFEST_FIELDS,
    ResolvedFile,
    resolve_primary_tumor_files,
)

_FILTER_BUILDERS = {
    "rna_seq": rna_seq_file_filter,
    "methylation": methylation_file_filter,
    "copy_number": copy_number_file_filter,
    "mutation": mutation_file_filter,
}


def resolve_layer_files(layer: str, file_hits: list[dict[str, Any]]) -> dict[str, ResolvedFile]:
    """Resolve one layer's file hits to one Primary Tumor file per patient.

    Args:
        layer: One of ``_FILTER_BUILDERS``' keys.
        file_hits: GDC file hits for this layer.

    Returns:
        A dict mapping case_id to a :class:`ResolvedFile`, using the
        copy-number-specific workflow-priority resolver for
        ``"copy_number"`` and the generic resolver otherwise.
    """
    if layer == "copy_number":
        return resolve_copy_number_files(file_hits)
    return resolve_primary_tumor_files(file_hits)


def main() -> None:
    """Ingest one omic layer's resolved files and write its manifest."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--layer", required=True, choices=sorted(_FILTER_BUILDERS))
    parser.add_argument("--config", type=Path, default=Path("workflows/config.yaml"))
    parser.add_argument("--case-ids-csv", type=Path, default=None)
    args = parser.parse_args()

    config = yaml.safe_load(args.config.read_text())
    data_dir = Path(config["data_dir"])
    case_ids_csv = args.case_ids_csv or data_dir / "processed" / "tnbc_case_ids.csv"
    raw_dir = data_dir / "raw" / args.layer
    processed_dir = data_dir / "processed"
    processed_dir.mkdir(parents=True, exist_ok=True)

    cohort = pd.read_csv(case_ids_csv)
    case_ids = cohort["bcr_patient_uuid"].tolist()

    client = GDCClient(config["gdc_api_base_url"])
    filters = _FILTER_BUILDERS[args.layer](case_ids)
    fields = CNV_MANIFEST_FIELDS if args.layer == "copy_number" else MANIFEST_FIELDS
    file_hits = list(client.query_files(filters, fields=list(fields)))
    resolved = resolve_layer_files(args.layer, file_hits)

    resolved_paths: dict[str, str] = {}
    for case_id, resolved_file in resolved.items():
        destination = raw_dir / resolved_file.file_name
        if not destination.exists():
            client.download_file(resolved_file.file_id, destination)
            record_download(
                source="GDC",
                query_description=f"TNBC cohort {args.layer} files, TCGA project {TCGA_PROJECT_ID}",
                accession_or_file_id=resolved_file.file_id,
                artifact_path=destination,
                extra={"case_id": case_id, "layer": args.layer},
            )
        resolved_paths[case_id] = str(destination)

    manifest_path = processed_dir / f"resolved_paths_{args.layer}.json"
    manifest_path.write_text(json.dumps(resolved_paths, indent=2, sort_keys=True))
    print(f"{args.layer}: resolved {len(resolved_paths)}/{len(case_ids)} patients")


if __name__ == "__main__":
    main()
