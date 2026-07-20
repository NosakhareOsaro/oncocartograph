"""Tests for oncocartograph.drug_targets.chembl_client.

Uses a fake requests.Session so no real network calls are made. Response
shapes mirror the real ChEMBL REST API, confirmed via live queries while
planning this work package (e.g. mechanism.json returning multiple
records per target, each with its own max_phase, requiring a max()
aggregation rather than assuming one record per target).
"""

from __future__ import annotations

from typing import Any

import pytest
import requests

from oncocartograph.drug_targets.chembl_client import ChEMBLClient, ChEMBLRequestError


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

    def get(self, url: str, **kwargs: Any) -> FakeResponse:
        self.calls.append({"url": url, **kwargs})
        item = self._responses.pop(0)
        if isinstance(item, Exception):
            raise item
        return item


def _make_client(session: FakeSession, **kwargs: Any) -> ChEMBLClient:
    return ChEMBLClient(
        base_url="https://www.ebi.ac.uk/chembl/api/data",
        session=session,  # type: ignore[arg-type]
        max_retries=kwargs.pop("max_retries", 3),
        **kwargs,
    )


def test_resolve_accessions_to_target_ids_exact_match() -> None:
    """An accession must resolve to the target whose components list contains it."""
    session = FakeSession(
        [
            FakeResponse(
                json_data={
                    "targets": [
                        {
                            "target_chembl_id": "CHEMBL4096",
                            "target_components": [{"accession": "P04637"}],
                        }
                    ]
                }
            )
        ]
    )
    client = _make_client(session)

    result = client.resolve_accessions_to_target_ids(["P04637"])

    assert result == {"P04637": "CHEMBL4096"}


def test_resolve_accessions_to_target_ids_missing_accession_is_none() -> None:
    """An accession with no matching target must map to None, not be absent."""
    session = FakeSession([FakeResponse(json_data={"targets": []})])
    client = _make_client(session)

    result = client.resolve_accessions_to_target_ids(["UNKNOWN"])

    assert result == {"UNKNOWN": None}


def test_resolve_accessions_to_target_ids_only_assigns_matching_component() -> None:
    """A response with multiple targets must assign each accession to its own target only."""
    session = FakeSession(
        [
            FakeResponse(
                json_data={
                    "targets": [
                        {
                            "target_chembl_id": "CHEMBL4096",
                            "target_components": [{"accession": "P04637"}],
                        },
                        {
                            "target_chembl_id": "CHEMBL1293249",
                            "target_components": [{"accession": "Q13887"}],
                        },
                    ]
                }
            )
        ]
    )
    client = _make_client(session)

    result = client.resolve_accessions_to_target_ids(["P04637", "Q13887"])

    assert result == {"P04637": "CHEMBL4096", "Q13887": "CHEMBL1293249"}


def test_resolve_accessions_to_target_ids_ignores_unrequested_accessions() -> None:
    """A target component accession we did not ask about must be ignored, not added."""
    session = FakeSession(
        [
            FakeResponse(
                json_data={
                    "targets": [
                        {
                            "target_chembl_id": "CHEMBL4096",
                            "target_components": [
                                {"accession": "P04637"},
                                {"accession": "UNASKED_ACCESSION"},
                            ],
                        }
                    ]
                }
            )
        ]
    )
    client = _make_client(session)

    result = client.resolve_accessions_to_target_ids(["P04637"])

    assert result == {"P04637": "CHEMBL4096"}
    assert "UNASKED_ACCESSION" not in result


def test_resolve_accessions_to_target_ids_chunks_across_batch_size() -> None:
    """More accessions than batch_size must trigger multiple requests."""
    session = FakeSession(
        [
            FakeResponse(
                json_data={
                    "targets": [
                        {"target_chembl_id": "CHEMBL1", "target_components": [{"accession": "A1"}]}
                    ]
                }
            ),
            FakeResponse(
                json_data={
                    "targets": [
                        {"target_chembl_id": "CHEMBL2", "target_components": [{"accession": "A2"}]}
                    ]
                }
            ),
        ]
    )
    client = _make_client(session)

    result = client.resolve_accessions_to_target_ids(["A1", "A2"], batch_size=1)

    assert len(session.calls) == 2
    assert result == {"A1": "CHEMBL1", "A2": "CHEMBL2"}


def test_fetch_max_phase_takes_maximum_across_multiple_mechanisms() -> None:
    """Real ChEMBL data has multiple mechanism records per target (e.g. TP53: phases 2 and 3);
    the result must be the max, not the first or last record."""
    session = FakeSession(
        [
            FakeResponse(
                json_data={
                    "mechanisms": [
                        {"target_chembl_id": "CHEMBL4096", "max_phase": 2},
                        {"target_chembl_id": "CHEMBL4096", "max_phase": 3},
                    ]
                }
            )
        ]
    )
    client = _make_client(session)

    result = client.fetch_max_phase(["CHEMBL4096"])

    assert result == {"CHEMBL4096": 3.0}


def test_fetch_max_phase_does_not_decrease_when_later_record_is_smaller() -> None:
    """A later mechanism record with a smaller max_phase must not overwrite the running max."""
    session = FakeSession(
        [
            FakeResponse(
                json_data={
                    "mechanisms": [
                        {"target_chembl_id": "CHEMBL1", "max_phase": 4},
                        {"target_chembl_id": "CHEMBL1", "max_phase": 1},
                    ]
                }
            )
        ]
    )
    client = _make_client(session)

    result = client.fetch_max_phase(["CHEMBL1"])

    assert result == {"CHEMBL1": 4.0}


def test_fetch_max_phase_missing_target_is_none() -> None:
    """A target with no mechanism records at all must map to None."""
    session = FakeSession([FakeResponse(json_data={"mechanisms": []})])
    client = _make_client(session)

    result = client.fetch_max_phase(["CHEMBL_NO_DRUGS"])

    assert result == {"CHEMBL_NO_DRUGS": None}


def test_fetch_max_phase_ignores_null_max_phase_records() -> None:
    """A mechanism record with max_phase=null must not overwrite a real value with None."""
    session = FakeSession(
        [
            FakeResponse(
                json_data={
                    "mechanisms": [
                        {"target_chembl_id": "CHEMBL1", "max_phase": 4},
                        {"target_chembl_id": "CHEMBL1", "max_phase": None},
                    ]
                }
            )
        ]
    )
    client = _make_client(session)

    result = client.fetch_max_phase(["CHEMBL1"])

    assert result == {"CHEMBL1": 4.0}


def test_fetch_max_phase_chunks_across_batch_size() -> None:
    """More target IDs than batch_size must trigger multiple requests."""
    session = FakeSession(
        [
            FakeResponse(json_data={"mechanisms": [{"target_chembl_id": "C1", "max_phase": 1}]}),
            FakeResponse(json_data={"mechanisms": [{"target_chembl_id": "C2", "max_phase": 2}]}),
        ]
    )
    client = _make_client(session)

    result = client.fetch_max_phase(["C1", "C2"], batch_size=1)

    assert len(session.calls) == 2
    assert result == {"C1": 1.0, "C2": 2.0}


def test_get_retries_on_connection_error_then_succeeds() -> None:
    """A transient connection error must be retried, not raised immediately."""
    session = FakeSession(
        [requests.ConnectionError("boom"), FakeResponse(json_data={"mechanisms": []})]
    )
    client = _make_client(session, max_retries=1)

    result = client.fetch_max_phase(["C1"])

    assert result == {"C1": None}
    assert len(session.calls) == 2


def test_get_raises_after_exhausting_retries() -> None:
    """Once retries are exhausted, ChEMBLRequestError must be raised."""
    session = FakeSession(
        [requests.ConnectionError("boom"), requests.ConnectionError("boom again")]
    )
    client = _make_client(session, max_retries=1)

    with pytest.raises(ChEMBLRequestError):
        client.fetch_max_phase(["C1"])


def test_get_treats_5xx_as_retryable() -> None:
    """A 5xx response should be retried like a connection error."""
    session = FakeSession(
        [FakeResponse(status_code=503), FakeResponse(json_data={"mechanisms": []})]
    )
    client = _make_client(session, max_retries=1)

    result = client.fetch_max_phase(["C1"])

    assert result == {"C1": None}


def test_resolve_accessions_empty_input_makes_no_request() -> None:
    """No accessions to resolve must return {} without making a request."""
    session = FakeSession([])
    client = _make_client(session)

    result = client.resolve_accessions_to_target_ids([])

    assert result == {}
    assert len(session.calls) == 0
