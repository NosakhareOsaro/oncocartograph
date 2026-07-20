#!/usr/bin/env python3
"""Pull the TCGA-BRCA clinical supplement and build the TNBC cohort + survival table.

Queries GDC for the "Clinical Supplement" / "BCR Biotab" file whose name
contains "clinical_patient" (the patient-level file with receptor status
and survival fields, distinct from the follow-up/drug/nte biotab files
GDC also returns for this filter -- confirmed via a live query before
this script was written), downloads it, classifies the TNBC sub-cohort
(docs/adr/0001), and derives the overall-survival duration/event table
(docs/adr/0007).

Usage:
    python workflows/scripts/ingest_cohort.py --config workflows/config.yaml
"""

from __future__ import annotations

import argparse
from pathlib import Path

import yaml

from oncocartograph.data_ingestion.clinical import (
    derive_survival_outcome,
    extract_receptor_status,
    read_biotab,
)
from oncocartograph.data_ingestion.gdc_client import GDCClient
from oncocartograph.data_ingestion.provenance import record_download
from oncocartograph.data_ingestion.tnbc_cohort import build_tnbc_cohort_audit, select_tnbc_cohort

_TCGA_PROJECT_ID = "TCGA-BRCA"
_CLINICAL_FILENAME_MARKER = "clinical_patient"


def find_clinical_patient_file(client: GDCClient) -> dict[str, str]:
    """Find the single TCGA-BRCA patient-level clinical BCR Biotab file.

    Args:
        client: A configured GDC client.

    Returns:
        The matching file hit dict (``file_id``, ``file_name``).

    Raises:
        ValueError: If not exactly one matching file is found.
    """
    filters = {
        "op": "and",
        "content": [
            {
                "op": "in",
                "content": {"field": "cases.project.project_id", "value": [_TCGA_PROJECT_ID]},
            },
            {"op": "in", "content": {"field": "data_type", "value": ["Clinical Supplement"]}},
            {"op": "in", "content": {"field": "data_format", "value": ["BCR Biotab"]}},
        ],
    }
    hits = list(client.query_files(filters, fields=["file_id", "file_name"]))
    matches = [h for h in hits if _CLINICAL_FILENAME_MARKER in h["file_name"]]
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one patient-level clinical file, found {len(matches)}: {matches}"
        )
    return matches[0]


def main() -> None:
    """Ingest the clinical file and write the TNBC cohort + survival outcome CSVs."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("workflows/config.yaml"))
    parser.add_argument("--raw-dir", type=Path, default=None)
    parser.add_argument("--processed-dir", type=Path, default=None)
    args = parser.parse_args()

    config = yaml.safe_load(args.config.read_text())
    data_dir = Path(config["data_dir"])
    raw_dir = args.raw_dir or data_dir / "raw" / "clinical"
    processed_dir = args.processed_dir or data_dir / "processed"
    processed_dir.mkdir(parents=True, exist_ok=True)

    client = GDCClient(config["gdc_api_base_url"])
    file_hit = find_clinical_patient_file(client)
    destination = raw_dir / file_hit["file_name"]
    client.download_file(file_hit["file_id"], destination)
    record_download(
        source="GDC",
        query_description=(
            "cases.project.project_id=TCGA-BRCA AND data_type=Clinical Supplement "
            "AND data_format=BCR Biotab, filtered to the patient-level file by name"
        ),
        accession_or_file_id=file_hit["file_id"],
        artifact_path=destination,
        extra={"file_name": file_hit["file_name"]},
    )

    clinical = read_biotab(destination)
    receptor_status = extract_receptor_status(clinical)
    audit = build_tnbc_cohort_audit(receptor_status)
    cohort = select_tnbc_cohort(audit)

    audit.to_csv(processed_dir / "tnbc_cohort_audit.csv", index=False)
    cohort.to_csv(processed_dir / "tnbc_cohort.csv", index=False)

    tnbc_barcodes = set(cohort["bcr_patient_barcode"])
    # tnbc_cohort.csv (like the real cohort-definition artifact it mirrors)
    # deliberately carries only receptor-status fields, not the GDC case
    # UUID -- ingest_omics.py needs that UUID to query per-patient omics
    # files, so it's recorded separately here rather than widening the
    # cited cohort-definition CSV's schema.
    case_ids = clinical.loc[
        clinical["bcr_patient_barcode"].isin(tnbc_barcodes),
        ["bcr_patient_uuid", "bcr_patient_barcode"],
    ]
    case_ids.to_csv(processed_dir / "tnbc_case_ids.csv", index=False)
    tnbc_clinical = clinical[clinical["bcr_patient_barcode"].isin(tnbc_barcodes)]
    survival = derive_survival_outcome(tnbc_clinical)
    survival.to_csv(processed_dir / "survival.csv")

    print(
        f"TNBC cohort: {len(cohort)} patients; survival table: {len(survival)} patients, "
        f"{int(survival['event'].sum())} events"
    )


if __name__ == "__main__":
    main()
