"""Per-omic GDC file discovery and download for the TCGA-BRCA TNBC cohort.

Builds the GDC filter expressions for each omic layer used in this project
(RNA-seq STAR-Counts, Illumina 450K methylation, gene-level copy number,
MC3 public MAF -- see ``docs/data_sources.md``), restricted to a
given set of case UUIDs (the TNBC cohort selected by
``oncocartograph.data_ingestion.tnbc_cohort``), and downloads the matching
files with a provenance record per file.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

from oncocartograph.data_ingestion.gdc_client import GDCClient
from oncocartograph.data_ingestion.provenance import ProvenanceRecord, record_download

#: TCGA project this pipeline is scoped to (see docs/adr/0001).
TCGA_PROJECT_ID = "TCGA-BRCA"

_FILE_FIELDS = ("file_id", "file_name", "data_category", "data_type", "cases.case_id")


def _in_filter(field: str, values: list[str]) -> dict[str, Any]:
    """Build a single GDC "field is in values" filter clause.

    Args:
        field: GDC entity field name, e.g. "data_category".
        values: Allowed values for that field.

    Returns:
        A GDC ``{"op": "in", ...}`` filter clause.
    """
    return {"op": "in", "content": {"field": field, "value": values}}


def _case_and_project_filter(
    case_ids: Sequence[str], extra_clauses: list[dict[str, Any]]
) -> dict[str, Any]:
    """Build the shared "project + case + ..." GDC filter prefix.

    Args:
        case_ids: GDC case UUIDs to restrict the query to.
        extra_clauses: Additional filter clauses specific to one omic
            layer (e.g. data_category/data_type/workflow_type constraints).

    Returns:
        A GDC ``filters`` expression combining project, case, and the
        caller-supplied clauses with ``"op": "and"``.
    """
    return {
        "op": "and",
        "content": [
            _in_filter("cases.project.project_id", [TCGA_PROJECT_ID]),
            _in_filter("cases.case_id", list(case_ids)),
            *extra_clauses,
        ],
    }


def rna_seq_file_filter(case_ids: Sequence[str]) -> dict[str, Any]:
    """GDC filter for harmonized RNA-seq STAR-Counts gene expression files.

    Args:
        case_ids: GDC case UUIDs to restrict the query to.

    Returns:
        A GDC filter expression for :func:`GDCClient.query_files`.
    """
    return _case_and_project_filter(
        case_ids,
        [
            _in_filter("data_category", ["Transcriptome Profiling"]),
            _in_filter("data_type", ["Gene Expression Quantification"]),
            _in_filter("analysis.workflow_type", ["STAR - Counts"]),
        ],
    )


def methylation_file_filter(case_ids: Sequence[str]) -> dict[str, Any]:
    """GDC filter for Illumina Infinium HumanMethylation450 (450K) files.

    Args:
        case_ids: GDC case UUIDs to restrict the query to.

    Returns:
        A GDC filter expression for :func:`GDCClient.query_files`.
    """
    return _case_and_project_filter(
        case_ids,
        [
            _in_filter("data_category", ["DNA Methylation"]),
            _in_filter("platform", ["Illumina Human Methylation 450"]),
            _in_filter("data_type", ["Methylation Beta Value"]),
        ],
    )


def copy_number_file_filter(case_ids: Sequence[str]) -> dict[str, Any]:
    """GDC filter for gene-level copy number files.

    These report absolute integer total copy number per gene, not
    GISTIC2 thresholded categorical calls as originally assumed --
    confirmed against real downloaded files, see
    ``docs/adr/0005-copy-number-workflow-and-transform.md``.

    Args:
        case_ids: GDC case UUIDs to restrict the query to.

    Returns:
        A GDC filter expression for :func:`GDCClient.query_files`.
    """
    return _case_and_project_filter(
        case_ids,
        [
            _in_filter("data_category", ["Copy Number Variation"]),
            _in_filter("data_type", ["Gene Level Copy Number"]),
        ],
    )


def mutation_file_filter(case_ids: Sequence[str]) -> dict[str, Any]:
    """GDC filter for MC3 public somatic mutation MAF files.

    Args:
        case_ids: GDC case UUIDs to restrict the query to.

    Returns:
        A GDC filter expression for :func:`GDCClient.query_files`.
    """
    return _case_and_project_filter(
        case_ids,
        [
            _in_filter("data_category", ["Simple Nucleotide Variation"]),
            _in_filter("data_type", ["Masked Somatic Mutation"]),
        ],
    )


def find_files(client: GDCClient, filters: dict[str, Any]) -> list[dict[str, Any]]:
    """Query the GDC ``files`` endpoint for a given filter expression.

    Args:
        client: A configured :class:`GDCClient`.
        filters: A GDC filter expression, e.g. from :func:`rna_seq_file_filter`.

    Returns:
        A list of matching file hit dicts (``file_id``, ``file_name``,
        ``data_category``, ``data_type``, ``cases``).
    """
    return list(client.query_files(filters, fields=_FILE_FIELDS))


def download_files(
    client: GDCClient,
    files: Sequence[dict[str, Any]],
    destination_dir: Path,
    *,
    source: str,
    query_description: str,
) -> list[ProvenanceRecord]:
    """Download a list of GDC files and record provenance for each.

    Args:
        client: A configured :class:`GDCClient`.
        files: File hit dicts as returned by :func:`find_files`.
        destination_dir: Directory to download files into (created if
            needed).
        source: Provenance ``source`` label, e.g. "GDC".
        query_description: Human-readable description of the query that
            produced ``files``, recorded in each file's provenance sidecar.

    Returns:
        One :class:`ProvenanceRecord` per downloaded file, in the same
        order as ``files``.
    """
    records = []
    for file_info in files:
        file_id = file_info["file_id"]
        file_name = file_info.get("file_name", file_id)
        destination = destination_dir / file_name
        client.download_file(file_id, destination)
        records.append(
            record_download(
                source=source,
                query_description=query_description,
                accession_or_file_id=file_id,
                artifact_path=destination,
                extra={"file_name": file_name, "data_type": file_info.get("data_type")},
            )
        )
    return records
