"""Typed client for the public GDC (Genomic Data Commons) REST API.

This client is a deliberately thin wrapper around
https://api.gdc.cancer.gov (``files``, ``cases``, and ``data`` endpoints).
It exists so the rest of the ingestion pipeline never constructs raw GDC
filter JSON or handles pagination/retries itself; see
``docs/adr/0004-gdc-rest-client-over-tcgabiolinks.md`` for why this talks to
GDC directly instead of going through TCGAbiolinks/R.
"""

from __future__ import annotations

import logging
import time
from collections.abc import Iterator, Mapping, Sequence
from pathlib import Path
from typing import Any

import requests

logger = logging.getLogger(__name__)

GDCFilter = Mapping[str, Any]
"""A GDC query filter, e.g. ``{"op": "in", "content": {"field": ..., "value": [...]}}``."""


class GDCRequestError(RuntimeError):
    """Raised when a GDC API request fails after all retries are exhausted."""


class GDCClient:
    """Minimal typed client for the GDC REST API.

    Args:
        base_url: Root URL of the GDC API (see ``Settings.gdc_api_base_url``).
        page_size: Number of hits requested per page for paginated queries.
        timeout_seconds: Per-request timeout in seconds.
        max_retries: Number of retries for transient failures (connection
            errors or 5xx responses) before raising :class:`GDCRequestError`.
        session: Optional pre-configured :class:`requests.Session`, mainly
            for test injection. A new session is created if not provided.
    """

    def __init__(
        self,
        base_url: str,
        *,
        page_size: int = 100,
        timeout_seconds: float = 30.0,
        max_retries: int = 3,
        session: requests.Session | None = None,
    ) -> None:
        """Initialise the client with connection and retry parameters."""
        self._base_url = base_url.rstrip("/")
        self._page_size = page_size
        self._timeout_seconds = timeout_seconds
        self._max_retries = max_retries
        self._session = session or requests.Session()

    def _request(
        self,
        method: str,
        path: str,
        *,
        params: Mapping[str, Any] | None = None,
        json_body: Mapping[str, Any] | None = None,
        stream: bool = False,
    ) -> requests.Response:
        """Issue an HTTP request with retry on transient failure.

        Args:
            method: HTTP method, e.g. "GET" or "POST".
            path: API path relative to the base URL, e.g. "/files".
            params: Query string parameters.
            json_body: JSON request body (used for POST queries).
            stream: Whether to stream the response body (for file downloads).

        Returns:
            The successful :class:`requests.Response`.

        Raises:
            GDCRequestError: If all retry attempts fail.
        """
        url = f"{self._base_url}{path}"
        last_error: Exception | None = None
        for attempt in range(self._max_retries + 1):
            try:
                response = self._session.request(
                    method,
                    url,
                    params=params,
                    json=json_body,
                    timeout=self._timeout_seconds,
                    stream=stream,
                )
                if response.status_code >= 500:
                    raise GDCRequestError(f"GDC API returned {response.status_code} for {url}")
                response.raise_for_status()
                return response
            except (requests.RequestException, GDCRequestError) as exc:
                last_error = exc
                if attempt < self._max_retries:
                    backoff_seconds = 2**attempt
                    logger.warning(
                        "GDC request to %s failed (attempt %d/%d): %s. Retrying in %ds.",
                        url,
                        attempt + 1,
                        self._max_retries + 1,
                        exc,
                        backoff_seconds,
                    )
                    time.sleep(backoff_seconds)
        raise GDCRequestError(
            f"GDC request to {url} failed after {self._max_retries + 1} attempts"
        ) from last_error

    def query_files(
        self,
        filters: GDCFilter,
        fields: Sequence[str],
    ) -> Iterator[dict[str, Any]]:
        """Query the GDC ``files`` endpoint, transparently paginating all hits.

        Args:
            filters: A GDC filter expression (see module-level ``GDCFilter``).
            fields: File-entity fields to request, e.g.
                ``["file_id", "file_name", "cases.case_id"]``.

        Yields:
            One dict per matching file, in GDC's returned order.
        """
        yield from self._paginated_query("/files", filters, fields)

    def query_cases(
        self,
        filters: GDCFilter,
        fields: Sequence[str],
    ) -> Iterator[dict[str, Any]]:
        """Query the GDC ``cases`` endpoint, transparently paginating all hits.

        Args:
            filters: A GDC filter expression (see module-level ``GDCFilter``).
            fields: Case-entity fields to request, e.g.
                ``["case_id", "submitter_id", "diagnoses.primary_diagnosis"]``.

        Yields:
            One dict per matching case, in GDC's returned order.
        """
        yield from self._paginated_query("/cases", filters, fields)

    def _paginated_query(
        self,
        path: str,
        filters: GDCFilter,
        fields: Sequence[str],
    ) -> Iterator[dict[str, Any]]:
        offset = 0
        total: int | None = None
        while total is None or offset < total:
            body = {
                "filters": dict(filters),
                "fields": ",".join(fields),
                "format": "JSON",
                "size": self._page_size,
                "from": offset,
            }
            response = self._request("POST", path, json_body=body)
            payload = response.json()["data"]
            hits = payload["hits"]
            total = payload["pagination"]["total"]
            yield from hits
            if not hits:
                break
            offset += len(hits)

    def download_file(self, file_id: str, destination: Path) -> Path:
        """Download a single GDC file by its ``file_id`` (UUID).

        Args:
            file_id: The GDC file UUID (as returned by :meth:`query_files`).
            destination: Path to write the downloaded file to. Parent
                directories are created if needed.

        Returns:
            The ``destination`` path, for chaining.
        """
        destination.parent.mkdir(parents=True, exist_ok=True)
        response = self._request("GET", f"/data/{file_id}", stream=True)
        with destination.open("wb") as fh:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                fh.write(chunk)
        return destination
