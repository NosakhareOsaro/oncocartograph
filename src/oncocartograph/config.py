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
        gdc_api_base_url: Base URL for the GDC REST API (see
            https://docs.gdc.cancer.gov/API/Users_Guide/Getting_Started/).
        gdc_page_size: Number of hits requested per page when paginating
            GDC ``files``/``cases`` queries.
        gdc_request_timeout_seconds: Per-request timeout for GDC API calls.
        gdc_max_retries: Number of retries for transient GDC API failures
            (connection errors, 5xx responses) before giving up.
        open_targets_api_base_url: Base URL for the Open Targets GraphQL
            API (https://api.platform.opentargets.org/api/v4/graphql).
        open_targets_batch_size: Number of targets requested per
            ``targets(ensemblIds: ...)``/``mapIds`` GraphQL call.
        open_targets_request_timeout_seconds: Per-request timeout for
            Open Targets API calls.
        open_targets_max_retries: Number of retries for transient Open
            Targets API failures before giving up.
        chembl_api_base_url: Base URL for the ChEMBL REST API
            (https://www.ebi.ac.uk/chembl/api/data).
        chembl_batch_size: Number of identifiers requested per batched
            ChEMBL ``__in`` filter call.
        chembl_request_timeout_seconds: Per-request timeout for ChEMBL
            API calls.
        chembl_max_retries: Number of retries for transient ChEMBL API
            failures before giving up.
    """

    model_config = SettingsConfigDict(
        env_prefix="ONCOCARTOGRAPH_",
        env_file=".env",
        extra="ignore",
    )

    data_dir: Path = Field(default=Path("data"))
    random_seed: int = Field(default=42)
    log_level: str = Field(default="INFO")

    gdc_api_base_url: str = Field(default="https://api.gdc.cancer.gov")
    gdc_page_size: int = Field(default=100, gt=0, le=2000)
    gdc_request_timeout_seconds: float = Field(default=30.0, gt=0)
    gdc_max_retries: int = Field(default=3, ge=0)

    open_targets_api_base_url: str = Field(
        default="https://api.platform.opentargets.org/api/v4/graphql"
    )
    open_targets_batch_size: int = Field(default=50, gt=0, le=200)
    open_targets_request_timeout_seconds: float = Field(default=30.0, gt=0)
    open_targets_max_retries: int = Field(default=3, ge=0)

    chembl_api_base_url: str = Field(default="https://www.ebi.ac.uk/chembl/api/data")
    chembl_batch_size: int = Field(default=50, gt=0, le=200)
    chembl_request_timeout_seconds: float = Field(default=30.0, gt=0)
    chembl_max_retries: int = Field(default=3, ge=0)

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
