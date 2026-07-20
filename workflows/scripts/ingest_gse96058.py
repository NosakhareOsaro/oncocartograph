#!/usr/bin/env python3
"""Download the GSE96058 (SCAN-B) clinical and expression files.

Usage:
    python workflows/scripts/ingest_gse96058.py --config workflows/config.yaml
"""

from __future__ import annotations

import argparse
from pathlib import Path
from urllib.parse import urlparse

import requests
import yaml

from oncocartograph.data_ingestion.provenance import record_download


def _download(url: str, destination: Path) -> None:
    """Stream-download one URL to disk, creating parent directories as needed."""
    if destination.exists():
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, stream=True, timeout=120) as response:
        response.raise_for_status()
        with destination.open("wb") as fh:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                fh.write(chunk)


def main() -> None:
    """Download GSE96058's clinical series-matrix files and expression matrix."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("workflows/config.yaml"))
    args = parser.parse_args()

    config = yaml.safe_load(args.config.read_text())
    external_dir = Path(config["data_dir"]) / "external" / "gse96058"

    for url in config["gse96058"]["clinical_urls"]:
        destination = external_dir / Path(urlparse(url).path).name
        _download(url, destination)
        record_download(
            source="GEO",
            query_description="GSE96058 series-matrix clinical/survival metadata",
            accession_or_file_id="GSE96058",
            artifact_path=destination,
            extra={"url": url},
        )

    expression_url = config["gse96058"]["expression_url"]
    expression_destination = external_dir / Path(urlparse(expression_url).path).name
    _download(expression_url, expression_destination)
    record_download(
        source="GEO",
        query_description="GSE96058 gene expression matrix (log2(FPKM+0.1))",
        accession_or_file_id="GSE96058",
        artifact_path=expression_destination,
        extra={"url": expression_url},
    )
    print(f"GSE96058 files downloaded to {external_dir}")


if __name__ == "__main__":
    main()
