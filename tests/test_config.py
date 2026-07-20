"""Tests for oncocartograph.config."""

from pathlib import Path

import pytest
from pydantic import ValidationError

from oncocartograph.config import Settings, get_settings


def test_settings_defaults() -> None:
    """Default settings should be deterministic and use a relative data dir."""
    settings = Settings()
    assert settings.data_dir == Path("data")
    assert settings.random_seed == 42
    assert settings.log_level == "INFO"
    assert settings.gdc_api_base_url == "https://api.gdc.cancer.gov"
    assert settings.gdc_page_size == 100
    assert settings.gdc_request_timeout_seconds == 30.0
    assert settings.gdc_max_retries == 3
    assert (
        settings.open_targets_api_base_url == "https://api.platform.opentargets.org/api/v4/graphql"
    )
    assert settings.open_targets_batch_size == 50
    assert settings.chembl_api_base_url == "https://www.ebi.ac.uk/chembl/api/data"
    assert settings.chembl_batch_size == 50


def test_gdc_page_size_must_be_positive() -> None:
    """A non-positive page size would make GDC pagination loop forever or fail."""
    with pytest.raises(ValidationError):
        Settings(gdc_page_size=0)


def test_derived_data_dirs_are_relative_to_data_dir() -> None:
    """raw/processed/external dirs must always nest under data_dir."""
    settings = Settings(data_dir=Path("/tmp/oncocartograph-data"))
    assert settings.raw_data_dir == Path("/tmp/oncocartograph-data/raw")
    assert settings.processed_data_dir == Path("/tmp/oncocartograph-data/processed")
    assert settings.external_data_dir == Path("/tmp/oncocartograph-data/external")


def test_get_settings_returns_settings_instance() -> None:
    """get_settings() is the public entry point and must return Settings."""
    assert isinstance(get_settings(), Settings)


def test_env_prefix_overrides_default(monkeypatch) -> None:  # type: ignore[no-untyped-def]
    """ONCOCARTOGRAPH_-prefixed env vars must override defaults."""
    monkeypatch.setenv("ONCOCARTOGRAPH_RANDOM_SEED", "7")
    settings = Settings()
    assert settings.random_seed == 7
