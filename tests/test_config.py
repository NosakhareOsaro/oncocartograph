"""Tests for oncocartograph.config."""

from pathlib import Path

from oncocartograph.config import Settings, get_settings


def test_settings_defaults() -> None:
    """Default settings should be deterministic and use a relative data dir."""
    settings = Settings()
    assert settings.data_dir == Path("data")
    assert settings.random_seed == 42
    assert settings.log_level == "INFO"


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
