"""Typed client for the public Open Targets GraphQL API.

Base URL, gene-symbol resolution (``mapIds``), and target tractability
lookup (``targets``), confirmed against the real API before this client
was written -- including that Open Targets sorts nothing for us and
returns an empty ``hits`` list (not an error) for an unresolvable gene
symbol, and that ``proteinIds`` includes non-canonical TrEMBL entries
alongside the canonical ``uniprot_swissprot`` one, which must be
filtered explicitly.
"""

from __future__ import annotations

import logging
import time
from collections.abc import Sequence
from typing import Any

import requests

logger = logging.getLogger(__name__)

_MAP_IDS_QUERY = """
query MapIds($queryTerms: [String!]!) {
  mapIds(queryTerms: $queryTerms, entityNames: ["target"]) {
    mappings {
      term
      hits {
        id
      }
    }
  }
}
"""

_TARGETS_QUERY = """
query Targets($ensemblIds: [String!]!) {
  targets(ensemblIds: $ensemblIds) {
    id
    approvedSymbol
    tractability {
      label
      modality
      value
    }
    proteinIds {
      id
      source
    }
  }
}
"""

#: The proteinIds source identifying a target's canonical, reviewed
#: UniProt accession (as opposed to unreviewed TrEMBL entries).
_CANONICAL_UNIPROT_SOURCE = "uniprot_swissprot"


class OpenTargetsRequestError(RuntimeError):
    """Raised when an Open Targets API request fails after all retries are exhausted."""


class OpenTargetsClient:
    """Minimal typed client for the Open Targets GraphQL API.

    Args:
        base_url: The Open Targets GraphQL endpoint (see
            ``Settings.open_targets_api_base_url``).
        timeout_seconds: Per-request timeout in seconds.
        max_retries: Number of retries for transient failures before
            raising :class:`OpenTargetsRequestError`.
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
        self._base_url = base_url
        self._timeout_seconds = timeout_seconds
        self._max_retries = max_retries
        self._session = session or requests.Session()

    def _post(self, query: str, variables: dict[str, Any]) -> dict[str, Any]:
        """Issue a GraphQL POST request with retry on transient failure.

        Args:
            query: The GraphQL query document.
            variables: GraphQL query variables.

        Returns:
            The response's ``data`` object.

        Raises:
            OpenTargetsRequestError: If all retry attempts fail, or the
                response contains GraphQL errors.
        """
        last_error: Exception | None = None
        for attempt in range(self._max_retries + 1):
            try:
                response = self._session.post(
                    self._base_url,
                    json={"query": query, "variables": variables},
                    timeout=self._timeout_seconds,
                )
                if response.status_code >= 500:
                    raise OpenTargetsRequestError(
                        f"Open Targets API returned {response.status_code}"
                    )
                response.raise_for_status()
                payload = response.json()
                if payload.get("errors"):
                    raise OpenTargetsRequestError(f"Open Targets API errors: {payload['errors']}")
                return dict(payload["data"])
            except (requests.RequestException, OpenTargetsRequestError) as exc:
                last_error = exc
                if attempt < self._max_retries:
                    backoff_seconds = 2**attempt
                    logger.warning(
                        "Open Targets request failed (attempt %d/%d): %s. Retrying in %ds.",
                        attempt + 1,
                        self._max_retries + 1,
                        exc,
                        backoff_seconds,
                    )
                    time.sleep(backoff_seconds)
        raise OpenTargetsRequestError(
            f"Open Targets request failed after {self._max_retries + 1} attempts"
        ) from last_error

    def map_symbols_to_ensembl_ids(self, symbols: Sequence[str]) -> dict[str, str | None]:
        """Batch-resolve gene symbols to Ensembl gene IDs via ``mapIds``.

        Args:
            symbols: Gene symbols to resolve (e.g. ``["TP53", "GTF3C1"]``).

        Returns:
            A dict mapping each input symbol to its resolved Ensembl ID,
            or ``None`` if the symbol could not be resolved to exactly
            one target. If a symbol resolves to multiple targets, the
            first hit is used (logged as ambiguous) rather than treated
            as unresolved, since Open Targets returns hits ranked by
            relevance.
        """
        if not symbols:
            return {}
        data = self._post(_MAP_IDS_QUERY, {"queryTerms": list(symbols)})
        results: dict[str, str | None] = {}
        for mapping in data["mapIds"]["mappings"]:
            hits = mapping["hits"]
            if not hits:
                results[mapping["term"]] = None
                continue
            if len(hits) > 1:
                logger.warning(
                    "Symbol %r resolved to %d targets; using the first (%s).",
                    mapping["term"],
                    len(hits),
                    hits[0]["id"],
                )
            results[mapping["term"]] = str(hits[0]["id"])
        return results

    def fetch_targets(
        self, ensembl_ids: Sequence[str], *, batch_size: int = 50
    ) -> dict[str, dict[str, Any]]:
        """Batch-fetch tractability and canonical UniProt accession for targets.

        Args:
            ensembl_ids: Ensembl gene IDs to fetch (unversioned, e.g.
                ``"ENSG00000141510"``).
            batch_size: Maximum number of IDs per GraphQL call; requests
                are chunked to stay within this limit.

        Returns:
            A dict mapping each resolved Ensembl ID to a dict with keys
            ``approved_symbol``, ``tractability`` (the raw bucket list,
            see :func:`oncocartograph.drug_targets.tractability.score_tractability`),
            and ``uniprot_accession`` (the canonical ``uniprot_swissprot``
            accession, or ``None`` if not found). IDs Open Targets has no
            record for are simply absent from the result.
        """
        results: dict[str, dict[str, Any]] = {}
        unique_ids = list(dict.fromkeys(ensembl_ids))
        for start in range(0, len(unique_ids), batch_size):
            chunk = unique_ids[start : start + batch_size]
            data = self._post(_TARGETS_QUERY, {"ensemblIds": chunk})
            for target in data["targets"]:
                accession = next(
                    (
                        p["id"]
                        for p in target["proteinIds"]
                        if p["source"] == _CANONICAL_UNIPROT_SOURCE
                    ),
                    None,
                )
                results[target["id"]] = {
                    "approved_symbol": target["approvedSymbol"],
                    "tractability": target["tractability"],
                    "uniprot_accession": accession,
                }
        return results
