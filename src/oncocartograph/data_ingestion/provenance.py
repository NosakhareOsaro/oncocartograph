"""Provenance logging for ingested data artifacts.

Every file this pipeline downloads must be traceable back to the exact
source, query, accession, and download time that produced it (see
``docs/data_sources.md``). This module writes one small JSON sidecar file
per downloaded artifact rather than relying on humans to remember to note
this information down separately.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Callable, Mapping
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

_CHUNK_SIZE = 1024 * 1024


def compute_sha256(path: Path) -> str:
    """Compute the SHA-256 checksum of a file's contents.

    Args:
        path: Path to the file to checksum.

    Returns:
        The hex-encoded SHA-256 digest.
    """
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(_CHUNK_SIZE), b""):
            digest.update(chunk)
    return digest.hexdigest()


@dataclass(frozen=True)
class ProvenanceRecord:
    """A single artifact's provenance metadata.

    Attributes:
        source: Where the artifact came from, e.g. "GDC" or "GEO".
        query_description: Human-readable description of the query that
            produced this artifact (e.g. the GDC filter used).
        accession_or_file_id: The source system's identifier for this
            artifact (a GDC file UUID, a GEO accession, etc.).
        downloaded_at: UTC timestamp of when the artifact was downloaded.
        sha256_checksum: SHA-256 checksum of the downloaded file's
            contents, for later integrity verification.
        extra: Any additional query parameters worth recording (page size,
            filters, API version, etc.).
    """

    source: str
    query_description: str
    accession_or_file_id: str
    downloaded_at: datetime
    sha256_checksum: str
    extra: Mapping[str, Any] = field(default_factory=dict)

    def to_json_dict(self) -> dict[str, Any]:
        """Serialise this record to a JSON-compatible dict.

        Returns:
            A dict with the timestamp rendered as an ISO-8601 string.
        """
        payload = asdict(self)
        payload["downloaded_at"] = self.downloaded_at.isoformat()
        return payload


def provenance_path_for(artifact_path: Path) -> Path:
    """Return the sidecar provenance file path for a downloaded artifact.

    Args:
        artifact_path: Path to the downloaded data file.

    Returns:
        A path alongside ``artifact_path`` with a ``.provenance.json``
        suffix appended.
    """
    return artifact_path.with_name(artifact_path.name + ".provenance.json")


def record_download(
    *,
    source: str,
    query_description: str,
    accession_or_file_id: str,
    artifact_path: Path,
    extra: Mapping[str, Any] | None = None,
    clock: Callable[[], datetime] = lambda: datetime.now(UTC),
) -> ProvenanceRecord:
    """Compute and persist a provenance record for a just-downloaded artifact.

    Args:
        source: Where the artifact came from, e.g. "GDC" or "GEO".
        query_description: Human-readable description of the query used.
        accession_or_file_id: The source system's identifier for this
            artifact.
        artifact_path: Path to the already-downloaded file to checksum.
        extra: Any additional query parameters worth recording.
        clock: Callable returning the current UTC time; overridable for
            deterministic testing.

    Returns:
        The :class:`ProvenanceRecord` that was written.
    """
    record = ProvenanceRecord(
        source=source,
        query_description=query_description,
        accession_or_file_id=accession_or_file_id,
        downloaded_at=clock(),
        sha256_checksum=compute_sha256(artifact_path),
        extra=dict(extra or {}),
    )
    provenance_path = provenance_path_for(artifact_path)
    provenance_path.write_text(json.dumps(record.to_json_dict(), indent=2, sort_keys=True))
    return record
