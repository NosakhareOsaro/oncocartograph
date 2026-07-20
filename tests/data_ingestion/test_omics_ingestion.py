"""Tests for oncocartograph.data_ingestion.omics_ingestion.

Uses a stub GDC client (no real network calls) and synthetic file-hit
fixtures shaped like real GDC API responses.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from oncocartograph.data_ingestion.omics_ingestion import (
    TCGA_PROJECT_ID,
    copy_number_file_filter,
    download_files,
    find_files,
    methylation_file_filter,
    mutation_file_filter,
    rna_seq_file_filter,
)

_CASE_IDS = ["case-uuid-1", "case-uuid-2"]


def _filter_field_values(filters: dict[str, Any], field: str) -> list[Any]:
    for clause in filters["content"]:
        if clause["content"]["field"] == field:
            value: list[Any] = clause["content"]["value"]
            return value
    raise AssertionError(f"field {field!r} not found in filter content")


def test_rna_seq_filter_restricts_to_project_cases_and_star_counts() -> None:
    """The RNA-seq filter must scope to TCGA-BRCA, the given cases, and STAR-Counts."""
    filters = rna_seq_file_filter(_CASE_IDS)

    assert _filter_field_values(filters, "cases.project.project_id") == [TCGA_PROJECT_ID]
    assert _filter_field_values(filters, "cases.case_id") == _CASE_IDS
    assert _filter_field_values(filters, "data_category") == ["Transcriptome Profiling"]
    assert _filter_field_values(filters, "analysis.workflow_type") == ["STAR - Counts"]


def test_methylation_filter_restricts_to_450k_platform() -> None:
    """The methylation filter must scope to the 450K platform and processed beta values only.

    Without the data_type constraint, GDC also returns raw .idat intensity
    files under the same data_category/platform -- a real bug found by
    running the live pull (see git history), which wasted ~1.8GB pulling
    files this pipeline never uses.
    """
    filters = methylation_file_filter(_CASE_IDS)

    assert _filter_field_values(filters, "data_category") == ["DNA Methylation"]
    assert _filter_field_values(filters, "platform") == ["Illumina Human Methylation 450"]
    assert _filter_field_values(filters, "data_type") == ["Methylation Beta Value"]


def test_copy_number_filter_restricts_to_gene_level_copy_number() -> None:
    """The CNV filter must scope to gene-level copy number data."""
    filters = copy_number_file_filter(_CASE_IDS)

    assert _filter_field_values(filters, "data_category") == ["Copy Number Variation"]
    assert _filter_field_values(filters, "data_type") == ["Gene Level Copy Number"]


def test_mutation_filter_restricts_to_masked_somatic_mutation() -> None:
    """The mutation filter must scope to Masked Somatic Mutation (MC3-derived) MAFs."""
    filters = mutation_file_filter(_CASE_IDS)

    assert _filter_field_values(filters, "data_category") == ["Simple Nucleotide Variation"]
    assert _filter_field_values(filters, "data_type") == ["Masked Somatic Mutation"]


class _StubGDCClient:
    """Stub matching the subset of GDCClient's interface omics_ingestion uses."""

    def __init__(self, hits: list[dict[str, Any]]) -> None:
        self._hits = hits
        self.queried_filters: list[dict[str, Any]] = []
        self.downloaded: list[tuple[str, Path]] = []

    def query_files(self, filters: dict[str, Any], fields: Any) -> list[dict[str, Any]]:
        self.queried_filters.append(filters)
        return self._hits

    def download_file(self, file_id: str, destination: Path) -> Path:
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(f"synthetic content for {file_id}".encode())
        self.downloaded.append((file_id, destination))
        return destination


def test_find_files_returns_query_files_results() -> None:
    """find_files must pass the filter through and return the raw hit list."""
    hits = [{"file_id": "file-1", "file_name": "a.tsv"}]
    client = _StubGDCClient(hits)

    result = find_files(client, {"op": "and", "content": []})  # type: ignore[arg-type]

    assert result == hits
    assert len(client.queried_filters) == 1


def test_download_files_downloads_and_writes_provenance(tmp_path: Path) -> None:
    """Each file must be downloaded to destination_dir/file_name with a provenance sidecar."""
    hits = [
        {"file_id": "file-1", "file_name": "sample1.tsv", "data_type": "Gene Expression"},
        {"file_id": "file-2", "file_name": "sample2.tsv", "data_type": "Gene Expression"},
    ]
    client = _StubGDCClient(hits)

    records = download_files(
        client,  # type: ignore[arg-type]
        hits,
        tmp_path,
        source="GDC",
        query_description="synthetic test query",
    )

    assert len(records) == 2
    assert (tmp_path / "sample1.tsv").exists()
    assert (tmp_path / "sample2.tsv").exists()

    sidecar = json.loads((tmp_path / "sample1.tsv.provenance.json").read_text())
    assert sidecar["source"] == "GDC"
    assert sidecar["accession_or_file_id"] == "file-1"
    assert sidecar["extra"]["file_name"] == "sample1.tsv"


def test_download_files_uses_file_id_as_filename_fallback(tmp_path: Path) -> None:
    """If a file hit has no file_name, the file_id must be used as the destination name."""
    hits = [{"file_id": "file-no-name"}]
    client = _StubGDCClient(hits)

    download_files(client, hits, tmp_path, source="GDC", query_description="q")  # type: ignore[arg-type]

    assert (tmp_path / "file-no-name").exists()
