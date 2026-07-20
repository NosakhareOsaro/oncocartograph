"""Typed client for the public ChEMBL REST API.

UniProt-accession-to-target resolution and max clinical trial phase
lookup, confirmed against the real API before this client was written --
including that ChEMBL's free-text target search is unreliable for exact
gene matching (searching "TP53" returned "TP53-binding protein 1" as the
top hit, not TP53 itself), so this client resolves targets by exact
UniProt accession match instead, restricted to ``target_type=SINGLE
PROTEIN`` to avoid multi-component complexes.
"""

from __future__ import annotations

import logging
import time
from collections.abc import Sequence
from typing import Any

import requests

logger = logging.getLogger(__name__)

#: Restrict target resolution to single-protein targets, not complexes.
_SINGLE_PROTEIN_TARGET_TYPE = "SINGLE PROTEIN"


class ChEMBLRequestError(RuntimeError):
    """Raised when a ChEMBL API request fails after all retries are exhausted."""


class ChEMBLClient:
    """Minimal typed client for the ChEMBL REST API.

    Args:
        base_url: The ChEMBL REST API base URL (see
            ``Settings.chembl_api_base_url``).
        timeout_seconds: Per-request timeout in seconds.
        max_retries: Number of retries for transient failures before
            raising :class:`ChEMBLRequestError`.
        session: Optional pre-configured :class:`requests.Session`,
            mainly for test injection.
    """

    def __init__(
        self,
        base_url: str,
        *,
        timeout_seconds: float = 30.0,
        max_retries: int = 3,
        session: requests.Session | None = None,
    ) -> None:
        """Initialise the client with connection and retry parameters."""
        self._base_url = base_url.rstrip("/")
        self._timeout_seconds = timeout_seconds
        self._max_retries = max_retries
        self._session = session or requests.Session()

    def _get(self, path: str, params: dict[str, Any]) -> dict[str, Any]:
        """Issue a GET request with retry on transient failure.

        Args:
            path: API path relative to the base URL, e.g. "/target.json".
            params: Query string parameters.

        Returns:
            The parsed JSON response body.

        Raises:
            ChEMBLRequestError: If all retry attempts fail.
        """
        url = f"{self._base_url}{path}"
        last_error: Exception | None = None
        for attempt in range(self._max_retries + 1):
            try:
                response = self._session.get(url, params=params, timeout=self._timeout_seconds)
                if response.status_code >= 500:
                    raise ChEMBLRequestError(
                        f"ChEMBL API returned {response.status_code} for {url}"
                    )
                response.raise_for_status()
                return dict(response.json())
            except (requests.RequestException, ChEMBLRequestError) as exc:
                last_error = exc
                if attempt < self._max_retries:
                    backoff_seconds = 2**attempt
                    logger.warning(
                        "ChEMBL request to %s failed (attempt %d/%d): %s. Retrying in %ds.",
                        url,
                        attempt + 1,
                        self._max_retries + 1,
                        exc,
                        backoff_seconds,
                    )
                    time.sleep(backoff_seconds)
        raise ChEMBLRequestError(
            f"ChEMBL request to {url} failed after {self._max_retries + 1} attempts"
        ) from last_error

    def resolve_accessions_to_target_ids(
        self, accessions: Sequence[str], *, batch_size: int = 50
    ) -> dict[str, str | None]:
        """Batch-resolve UniProt accessions to ChEMBL single-protein target IDs.

        Args:
            accessions: UniProt accessions to resolve (e.g. ``["P04637"]``).
            batch_size: Maximum number of accessions per request; requests
                are chunked to stay within this limit.

        Returns:
            A dict mapping each input accession to its ``target_chembl_id``,
            or ``None`` if no single-protein ChEMBL target exists for it.
        """
        results: dict[str, str | None] = dict.fromkeys(dict.fromkeys(accessions))
        unique_accessions = list(results)
        for start in range(0, len(unique_accessions), batch_size):
            chunk = unique_accessions[start : start + batch_size]
            payload = self._get(
                "/target.json",
                {
                    "target_components__accession__in": ",".join(chunk),
                    "target_type": _SINGLE_PROTEIN_TARGET_TYPE,
                    "limit": batch_size,
                },
            )
            for target in payload["targets"]:
                for component in target["target_components"]:
                    accession = component["accession"]
                    if accession in results:
                        results[accession] = target["target_chembl_id"]
        return results

    def fetch_max_phase(
        self, target_chembl_ids: Sequence[str], *, batch_size: int = 50
    ) -> dict[str, float | None]:
        """Batch-fetch the maximum clinical trial phase reached against each target.

        Args:
            target_chembl_ids: ChEMBL target IDs (e.g. ``["CHEMBL4096"]``).
            batch_size: Maximum number of target IDs per request; requests
                are chunked to stay within this limit.

        Returns:
            A dict mapping each input target ID to the maximum
            ``max_phase`` (0-4) across all of ChEMBL's recorded
            mechanisms against it, or ``None`` if no mechanism record
            exists for that target.
        """
        results: dict[str, float | None] = dict.fromkeys(dict.fromkeys(target_chembl_ids))
        unique_ids = list(results)
        for start in range(0, len(unique_ids), batch_size):
            chunk = unique_ids[start : start + batch_size]
            payload = self._get(
                "/mechanism.json",
                {"target_chembl_id__in": ",".join(chunk), "limit": 1000},
            )
            for mechanism in payload["mechanisms"]:
                target_id = mechanism["target_chembl_id"]
                max_phase = mechanism["max_phase"]
                if max_phase is None:
                    continue
                current = results.get(target_id)
                if current is None or max_phase > current:
                    results[target_id] = float(max_phase)
        return results
