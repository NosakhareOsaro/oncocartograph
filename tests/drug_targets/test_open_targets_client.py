"""Tests for oncocartograph.drug_targets.open_targets_client.

Uses a fake requests.Session so no real network calls are made. Response
shapes mirror the real Open Targets GraphQL API, confirmed via live
queries while planning this work package (e.g. mapIds returning an empty
hits list rather than an error for an unresolvable symbol, proteinIds
including non-canonical TrEMBL entries alongside the canonical
uniprot_swissprot one).
"""

from __future__ import annotations

from typing import Any

import pytest
import requests

from oncocartograph.drug_targets.open_targets_client import (
    OpenTargetsClient,
    OpenTargetsRequestError,
)


class FakeResponse:
    def __init__(self, status_code: int = 200, json_data: dict[str, Any] | None = None) -> None:
        self.status_code = status_code
        self._json_data = json_data or {}

    def json(self) -> dict[str, Any]:
        return self._json_data

    def raise_for_status(self) -> None:
        if self.status_code >= 400:
            raise requests.HTTPError(f"{self.status_code} error")


class FakeSession:
    def __init__(self, responses: list[FakeResponse | Exception]) -> None:
        self._responses = list(responses)
        self.calls: list[dict[str, Any]] = []

    def post(self, url: str, **kwargs: Any) -> FakeResponse:
        self.calls.append({"url": url, **kwargs})
        item = self._responses.pop(0)
        if isinstance(item, Exception):
            raise item
        return item


def _make_client(session: FakeSession, **kwargs: Any) -> OpenTargetsClient:
    return OpenTargetsClient(
        base_url="https://api.platform.opentargets.org/api/v4/graphql",
        session=session,  # type: ignore[arg-type]
        max_retries=kwargs.pop("max_retries", 3),
        **kwargs,
    )


def test_map_symbols_to_ensembl_ids_resolves_known_symbols() -> None:
    """A resolvable symbol must map to its single hit's Ensembl ID."""
    session = FakeSession(
        [
            FakeResponse(
                json_data={
                    "data": {
                        "mapIds": {
                            "mappings": [
                                {"term": "TP53", "hits": [{"id": "ENSG00000141510"}]},
                                {"term": "GTF3C1", "hits": [{"id": "ENSG00000077235"}]},
                            ]
                        }
                    }
                }
            )
        ]
    )
    client = _make_client(session)

    result = client.map_symbols_to_ensembl_ids(["TP53", "GTF3C1"])

    assert result == {"TP53": "ENSG00000141510", "GTF3C1": "ENSG00000077235"}


def test_map_symbols_to_ensembl_ids_returns_none_for_unresolvable_symbol() -> None:
    """An empty hits list (real API behaviour, not an error) must map to None."""
    session = FakeSession(
        [
            FakeResponse(
                json_data={"data": {"mapIds": {"mappings": [{"term": "NOTAGENE", "hits": []}]}}}
            )
        ]
    )
    client = _make_client(session)

    result = client.map_symbols_to_ensembl_ids(["NOTAGENE"])

    assert result == {"NOTAGENE": None}


def test_map_symbols_to_ensembl_ids_uses_first_hit_when_ambiguous() -> None:
    """A symbol resolving to multiple targets must use the first hit, not error."""
    session = FakeSession(
        [
            FakeResponse(
                json_data={
                    "data": {
                        "mapIds": {
                            "mappings": [
                                {
                                    "term": "AMBIGUOUS",
                                    "hits": [{"id": "ENSG00000000001"}, {"id": "ENSG00000000002"}],
                                }
                            ]
                        }
                    }
                }
            )
        ]
    )
    client = _make_client(session)

    result = client.map_symbols_to_ensembl_ids(["AMBIGUOUS"])

    assert result == {"AMBIGUOUS": "ENSG00000000001"}


def test_map_symbols_to_ensembl_ids_empty_input_returns_empty_dict() -> None:
    """No symbols to resolve must return {} without making a request."""
    session = FakeSession([])
    client = _make_client(session)

    result = client.map_symbols_to_ensembl_ids([])

    assert result == {}
    assert len(session.calls) == 0


def test_fetch_targets_extracts_tractability_and_canonical_accession() -> None:
    """Canonical uniprot_swissprot accession must be selected over TrEMBL entries."""
    session = FakeSession(
        [
            FakeResponse(
                json_data={
                    "data": {
                        "targets": [
                            {
                                "id": "ENSG00000141510",
                                "approvedSymbol": "TP53",
                                "tractability": [
                                    {"label": "Approved Drug", "modality": "SM", "value": False}
                                ],
                                "proteinIds": [
                                    {"id": "A0A087WT22", "source": "uniprot_trembl"},
                                    {"id": "P04637", "source": "uniprot_swissprot"},
                                ],
                            }
                        ]
                    }
                }
            )
        ]
    )
    client = _make_client(session)

    result = client.fetch_targets(["ENSG00000141510"])

    assert result["ENSG00000141510"]["approved_symbol"] == "TP53"
    assert result["ENSG00000141510"]["uniprot_accession"] == "P04637"
    assert len(result["ENSG00000141510"]["tractability"]) == 1


def test_fetch_targets_handles_missing_canonical_accession() -> None:
    """A target with no uniprot_swissprot entry must get None, not crash."""
    session = FakeSession(
        [
            FakeResponse(
                json_data={
                    "data": {
                        "targets": [
                            {
                                "id": "ENSG00000000001",
                                "approvedSymbol": "FAKE1",
                                "tractability": [],
                                "proteinIds": [{"id": "X1", "source": "uniprot_trembl"}],
                            }
                        ]
                    }
                }
            )
        ]
    )
    client = _make_client(session)

    result = client.fetch_targets(["ENSG00000000001"])

    assert result["ENSG00000000001"]["uniprot_accession"] is None


def test_fetch_targets_chunks_across_batch_size() -> None:
    """More IDs than batch_size must trigger multiple requests, results merged."""
    session = FakeSession(
        [
            FakeResponse(
                json_data={
                    "data": {
                        "targets": [
                            {
                                "id": "ENSG1",
                                "approvedSymbol": "G1",
                                "tractability": [],
                                "proteinIds": [],
                            }
                        ]
                    }
                }
            ),
            FakeResponse(
                json_data={
                    "data": {
                        "targets": [
                            {
                                "id": "ENSG2",
                                "approvedSymbol": "G2",
                                "tractability": [],
                                "proteinIds": [],
                            }
                        ]
                    }
                }
            ),
        ]
    )
    client = _make_client(session)

    result = client.fetch_targets(["ENSG1", "ENSG2"], batch_size=1)

    assert len(session.calls) == 2
    assert set(result) == {"ENSG1", "ENSG2"}


def test_fetch_targets_deduplicates_repeated_ids() -> None:
    """A repeated Ensembl ID in the input must not trigger a duplicate lookup."""
    session = FakeSession(
        [
            FakeResponse(
                json_data={
                    "data": {
                        "targets": [
                            {
                                "id": "ENSG1",
                                "approvedSymbol": "G1",
                                "tractability": [],
                                "proteinIds": [],
                            }
                        ]
                    }
                }
            )
        ]
    )
    client = _make_client(session)

    result = client.fetch_targets(["ENSG1", "ENSG1"], batch_size=50)

    assert len(session.calls) == 1
    assert set(result) == {"ENSG1"}


def test_post_retries_on_connection_error_then_succeeds() -> None:
    """A transient connection error must be retried, not raised immediately."""
    session = FakeSession(
        [
            requests.ConnectionError("boom"),
            FakeResponse(json_data={"data": {"mapIds": {"mappings": []}}}),
        ]
    )
    client = _make_client(session, max_retries=1)

    result = client.map_symbols_to_ensembl_ids(["X"])

    assert result == {}
    assert len(session.calls) == 2


def test_post_raises_after_exhausting_retries() -> None:
    """Once retries are exhausted, OpenTargetsRequestError must be raised."""
    session = FakeSession(
        [requests.ConnectionError("boom"), requests.ConnectionError("boom again")]
    )
    client = _make_client(session, max_retries=1)

    with pytest.raises(OpenTargetsRequestError):
        client.map_symbols_to_ensembl_ids(["X"])


def test_post_raises_on_graphql_errors_field() -> None:
    """A 200 response containing a GraphQL 'errors' field must raise, not return partial data.

    The final raised error carries a generic "failed after N attempts"
    message (matching the GDCClient retry-wrapping convention elsewhere
    in this project), with the original GraphQL error chained as its
    cause rather than surfaced in the message text directly.
    """
    session = FakeSession([FakeResponse(json_data={"errors": [{"message": "Cannot query field"}]})])
    client = _make_client(session, max_retries=0)

    with pytest.raises(OpenTargetsRequestError) as exc_info:
        client.map_symbols_to_ensembl_ids(["X"])

    assert "Cannot query field" in str(exc_info.value.__cause__)


def test_post_treats_5xx_as_retryable() -> None:
    """A 5xx response should be retried like a connection error."""
    session = FakeSession(
        [
            FakeResponse(status_code=503),
            FakeResponse(json_data={"data": {"mapIds": {"mappings": []}}}),
        ]
    )
    client = _make_client(session, max_retries=1)

    result = client.map_symbols_to_ensembl_ids(["X"])

    assert result == {}
