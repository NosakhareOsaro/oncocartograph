"""Central pipeline configuration.

All paths, cache locations, and default stochastic seeds used across
OncoCartograph are defined here via ``pydantic-settings`` rather than as
constants scattered through individual modules. Values can be overridden by
environment variables prefixed ``ONCOCARTOGRAPH_`` or by a ``.env`` file, so
the same codebase runs unmodified on a laptop, CI, or an HPC cluster.
"""

from pathlib import Path

from pydantic import Field
from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    """Runtime configuration for the OncoCartograph pipeline.

    Attributes:
        data_dir: Root directory under which raw, processed, and external
            data subdirectories are created.
        random_seed: Default random seed for stochastic pipeline steps
            (e.g. MOFA+ initialisation, train/test splits). Individual
            pipeline stages may accept an explicit override, but must log
            whichever seed is actually used.
        log_level: Python ``logging`` level name applied to all
            OncoCartograph loggers.
    """

    model_config = SettingsConfigDict(
        env_prefix="ONCOCARTOGRAPH_",
        env_file=".env",
        extra="ignore",
    )

    data_dir: Path = Field(default=Path("data"))
    random_seed: int = Field(default=42)
    log_level: str = Field(default="INFO")

    @property
    def raw_data_dir(self) -> Path:
        """Directory for unmodified, as-downloaded source data."""
        return self.data_dir / "raw"

    @property
    def processed_data_dir(self) -> Path:
        """Directory for pipeline-generated intermediate and final data."""
        return self.data_dir / "processed"

    @property
    def external_data_dir(self) -> Path:
        """Directory for external validation cohort data (e.g. GSE96058)."""
        return self.data_dir / "external"


def get_settings() -> Settings:
    """Return a fresh :class:`Settings` instance loaded from env/.env.

    Returns:
        A ``Settings`` instance reflecting current environment variables
        and any ``.env`` file in the working directory.
    """
    return Settings()
