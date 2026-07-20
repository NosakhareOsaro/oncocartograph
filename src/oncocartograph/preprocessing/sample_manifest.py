"""Resolve raw GDC file hits to one canonical Primary Tumor file per patient.

Every omic layer ingested in ``feat/data-ingestion`` can return more files
per patient than one: RNA-seq/methylation may include matched normal or
metastatic samples alongside the primary tumor, and copy number/mutation
files (tumor-normal paired callers) list *both* the tumor and germline
reference sample against the same file. This module isolates that
resolution logic -- restrict to Primary Tumor, deterministically pick one
file when a patient has more than one -- so every per-omic preprocessing
module consumes an already-resolved one-file-per-patient mapping rather
than reimplementing this each time.

Patients with no Primary Tumor file for a given omic layer are simply
absent from the resolved mapping, not an error: MOFA+ is designed to
factorize partially-overlapping views, so this project does not force a
smaller complete-case cohort across all four omics.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Any

#: GDC sample_type value identifying a primary tumor sample.
PRIMARY_TUMOR_SAMPLE_TYPE = "Primary Tumor"

#: GDC file fields required for manifest resolution; pass to GDCClient.query_files.
MANIFEST_FIELDS = ("file_id", "file_name", "cases.case_id", "cases.samples.sample_type")


def file_sample_types(file_hit: dict[str, Any]) -> set[str]:
    """Extract the set of sample_type values associated with a GDC file hit.

    Args:
        file_hit: A GDC file hit dict including expanded ``cases.samples``.

    Returns:
        The set of distinct sample_type strings across all cases/samples
        listed for this file (usually one case, one or two samples).
    """
    types: set[str] = set()
    for case in file_hit.get("cases", []):
        for sample in case.get("samples", []):
            sample_type = sample.get("sample_type")
            if sample_type:
                types.add(sample_type)
    return types


def file_case_id(file_hit: dict[str, Any]) -> str:
    """Extract the case UUID a GDC file hit belongs to.

    Args:
        file_hit: A GDC file hit dict including expanded ``cases.case_id``.

    Returns:
        The case UUID string.

    Raises:
        KeyError: If the file has no associated case. This shouldn't
            happen for files queried via a case-scoped filter, but fails
            loudly rather than silently skipping or picking an arbitrary
            case.
    """
    cases = file_hit.get("cases", [])
    if not cases:
        raise KeyError(f"file {file_hit.get('file_id')!r} has no associated case")
    return str(cases[0]["case_id"])


def is_primary_tumor(file_hit: dict[str, Any]) -> bool:
    """Check whether a file hit is associated with a Primary Tumor sample.

    Args:
        file_hit: A GDC file hit dict including expanded ``cases.samples``.

    Returns:
        True if "Primary Tumor" appears among the file's sample types
        (paired-caller files list both tumor and normal; this checks
        membership, not exclusivity).
    """
    return PRIMARY_TUMOR_SAMPLE_TYPE in file_sample_types(file_hit)


def default_tie_break(file_hits: Sequence[dict[str, Any]]) -> dict[str, Any]:
    """Deterministically choose one file when a patient has more than one.

    Args:
        file_hits: Two or more file hits for the same case, all Primary
            Tumor.

    Returns:
        The hit with the lexicographically smallest ``file_id``. The
        choice of "smallest file_id" is arbitrary with respect to biology
        -- there is no principled reason to prefer one aliquot over
        another for this project's purposes -- but it is fully
        reproducible given the same input, which is what matters here.
    """
    return sorted(file_hits, key=lambda hit: hit["file_id"])[0]


@dataclass(frozen=True)
class ResolvedFile:
    """One patient's chosen Primary Tumor file for a given omic layer.

    Attributes:
        case_id: GDC case UUID.
        file_id: GDC file UUID of the chosen file.
        file_name: The chosen file's name, for locating it on disk under
            ``data/raw/<layer>/``.
    """

    case_id: str
    file_id: str
    file_name: str


def resolve_primary_tumor_files(
    file_hits: Sequence[dict[str, Any]],
    *,
    tie_break: Callable[[Sequence[dict[str, Any]]], dict[str, Any]] = default_tie_break,
) -> dict[str, ResolvedFile]:
    """Resolve GDC file hits down to one Primary Tumor file per case.

    Args:
        file_hits: GDC file hits, e.g. from ``GDCClient.query_files`` with
            fields including ``cases.case_id`` and
            ``cases.samples.sample_type`` (see :data:`MANIFEST_FIELDS`).
        tie_break: Callable choosing one hit when a case has more than one
            Primary Tumor file. Deterministic by default.

    Returns:
        A dict mapping case_id to :class:`ResolvedFile`, one entry per
        case with at least one Primary Tumor file.
    """
    by_case: dict[str, list[dict[str, Any]]] = {}
    for hit in file_hits:
        if not is_primary_tumor(hit):
            continue
        by_case.setdefault(file_case_id(hit), []).append(hit)

    resolved: dict[str, ResolvedFile] = {}
    for case_id, hits in by_case.items():
        chosen = hits[0] if len(hits) == 1 else tie_break(hits)
        resolved[case_id] = ResolvedFile(
            case_id=case_id,
            file_id=chosen["file_id"],
            file_name=chosen.get("file_name", chosen["file_id"]),
        )
    return resolved
