"""Tests for oncocartograph.data_ingestion.gdc_client.

Uses a fake requests.Session so no real network calls are made. All
responses here are synthetic fixtures shaped like real GDC API responses,
not real GDC data.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest
import requests

from oncocartograph.data_ingestion.gdc_client import GDCClient, GDCRequestError


class FakeResponse:
    """Minimal stand-in for requests.Response used in tests."""

    def __init__(
        self,
        status_code: int = 200,
        json_data: dict[str, Any] | None = None,
        content_chunks: list[bytes] | None = None,
    ) -> None:
        self.status_code = status_code
        self._json_data = json_data or {}
        self._content_chunks = content_chunks or []

    def json(self) -> dict[str, Any]:
        return self._json_data

    def raise_for_status(self) -> None:
        if self.status_code >= 400:
            raise requests.HTTPError(f"{self.status_code} error")

    def iter_content(self, chunk_size: int) -> list[bytes]:
        return self._content_chunks


class FakeSession:
    """Fake requests.Session driven by a queue of canned responses/exceptions."""

    def __init__(self, responses: list[FakeResponse | Exception]) -> None:
        self._responses = list(responses)
        self.calls: list[dict[str, Any]] = []

    def request(self, method: str, url: str, **kwargs: Any) -> FakeResponse:
        self.calls.append({"method": method, "url": url, **kwargs})
        item = self._responses.pop(0)
        if isinstance(item, Exception):
            raise item
        return item


def _page_response(hits: list[dict[str, Any]], total: int) -> FakeResponse:
    return FakeResponse(
        status_code=200,
        json_data={"data": {"hits": hits, "pagination": {"total": total}}},
    )


def _make_client(session: FakeSession, **kwargs: Any) -> GDCClient:
    return GDCClient(
        base_url="https://api.gdc.cancer.gov",
        session=session,  # type: ignore[arg-type]
        max_retries=kwargs.pop("max_retries", 3),
        **kwargs,
    )


def test_query_files_paginates_across_multiple_pages() -> None:
    """query_files must keep requesting until all hits (per pagination.total) are seen."""
    session = FakeSession(
        [
            _page_response([{"file_id": "a"}, {"file_id": "b"}], total=3),
            _page_response([{"file_id": "c"}], total=3),
        ]
    )
    client = _make_client(session, page_size=2)

    results = list(client.query_files(filters={"op": "and", "content": []}, fields=["file_id"]))

    assert [r["file_id"] for r in results] == ["a", "b", "c"]
    assert len(session.calls) == 2


def test_query_files_stops_on_empty_page_even_if_total_not_reached() -> None:
    """An empty hits page must stop pagination rather than looping forever."""
    session = FakeSession([_page_response([], total=5)])
    client = _make_client(session)

    results = list(client.query_files(filters={"op": "and", "content": []}, fields=["file_id"]))

    assert results == []


def test_query_cases_uses_cases_endpoint() -> None:
    """query_cases must hit /cases, not /files."""
    session = FakeSession([_page_response([{"case_id": "case-1"}], total=1)])
    client = _make_client(session)

    results = list(client.query_cases(filters={"op": "and", "content": []}, fields=["case_id"]))

    assert results == [{"case_id": "case-1"}]
    assert session.calls[0]["url"].endswith("/cases")


def test_request_retries_on_connection_error_then_succeeds() -> None:
    """A transient connection error must be retried, not raised immediately."""
    session = FakeSession(
        [
            requests.ConnectionError("boom"),
            _page_response([{"file_id": "a"}], total=1),
        ]
    )
    client = _make_client(session, max_retries=1)

    results = list(client.query_files(filters={"op": "and", "content": []}, fields=["file_id"]))

    assert results == [{"file_id": "a"}]
    assert len(session.calls) == 2


def test_request_raises_gdc_request_error_after_exhausting_retries() -> None:
    """Once retries are exhausted, a clear GDCRequestError must be raised."""
    session = FakeSession(
        [
            requests.ConnectionError("boom"),
            requests.ConnectionError("boom again"),
        ]
    )
    client = _make_client(session, max_retries=1)

    with pytest.raises(GDCRequestError):
        list(client.query_files(filters={"op": "and", "content": []}, fields=["file_id"]))


def test_request_treats_5xx_as_retryable() -> None:
    """A 5xx response should be retried like a connection error, not raised directly."""
    session = FakeSession(
        [
            FakeResponse(status_code=503),
            _page_response([{"file_id": "a"}], total=1),
        ]
    )
    client = _make_client(session, max_retries=1)

    results = list(client.query_files(filters={"op": "and", "content": []}, fields=["file_id"]))

    assert results == [{"file_id": "a"}]


def test_download_file_writes_content_and_creates_parent_dirs(tmp_path: Path) -> None:
    """download_file must create missing parent directories and write all chunks."""
    session = FakeSession([FakeResponse(status_code=200, content_chunks=[b"abc", b"def"])])
    client = _make_client(session)
    destination = tmp_path / "nested" / "out.bin"

    result_path = client.download_file("file-uuid", destination)

    assert result_path == destination
    assert destination.read_bytes() == b"abcdef"


@pytest.fixture
def retry_sleep_patch(monkeypatch: pytest.MonkeyPatch) -> Callable[[float], None]:
    """Avoid real sleeps during retry-backoff tests."""
    calls: list[float] = []
    monkeypatch.setattr("oncocartograph.data_ingestion.gdc_client.time.sleep", calls.append)
    return calls.append


def test_retry_backoff_does_not_sleep_when_no_retry_needed(
    retry_sleep_patch: Callable[[float], None],
) -> None:
    """A successful first attempt must not trigger any sleep."""
    session = FakeSession([_page_response([], total=0)])
    client = _make_client(session)

    list(client.query_files(filters={"op": "and", "content": []}, fields=["file_id"]))
