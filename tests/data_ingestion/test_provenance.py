"""Tests for oncocartograph.data_ingestion.provenance."""

from __future__ import annotations

import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path

from oncocartograph.data_ingestion.provenance import (
    compute_sha256,
    provenance_path_for,
    record_download,
)


def test_compute_sha256_matches_hashlib_reference(tmp_path: Path) -> None:
    """compute_sha256 must match hashlib's own digest for the same content."""
    file_path = tmp_path / "artifact.bin"
    file_path.write_bytes(b"some synthetic file content")

    result = compute_sha256(file_path)

    assert result == hashlib.sha256(b"some synthetic file content").hexdigest()


def test_provenance_path_for_appends_suffix() -> None:
    """The sidecar path must be the artifact path with .provenance.json appended."""
    artifact = Path("/data/raw/sample.tsv")

    assert provenance_path_for(artifact) == Path("/data/raw/sample.tsv.provenance.json")


def test_record_download_writes_expected_json(tmp_path: Path) -> None:
    """record_download must write a JSON sidecar with all expected fields."""
    artifact_path = tmp_path / "sample.tsv"
    artifact_path.write_bytes(b"synthetic content")
    fixed_time = datetime(2026, 7, 20, 12, 0, 0, tzinfo=UTC)

    record = record_download(
        source="GDC",
        query_description="files endpoint, TCGA-BRCA STAR-Counts",
        accession_or_file_id="fake-file-uuid-0001",
        artifact_path=artifact_path,
        extra={"data_category": "Transcriptome Profiling"},
        clock=lambda: fixed_time,
    )

    assert record.source == "GDC"
    assert record.sha256_checksum == hashlib.sha256(b"synthetic content").hexdigest()
    assert record.downloaded_at == fixed_time

    sidecar_path = tmp_path / "sample.tsv.provenance.json"
    assert sidecar_path.exists()
    payload = json.loads(sidecar_path.read_text())
    assert payload["source"] == "GDC"
    assert payload["accession_or_file_id"] == "fake-file-uuid-0001"
    assert payload["downloaded_at"] == "2026-07-20T12:00:00+00:00"
    assert payload["extra"] == {"data_category": "Transcriptome Profiling"}


def test_record_download_defaults_extra_to_empty_dict(tmp_path: Path) -> None:
    """Omitting extra must produce an empty dict, not None, in the record."""
    artifact_path = tmp_path / "sample2.tsv"
    artifact_path.write_bytes(b"more synthetic content")

    record = record_download(
        source="GEO",
        query_description="GSE96058 series matrix",
        accession_or_file_id="GSE96058",
        artifact_path=artifact_path,
    )

    assert record.extra == {}
