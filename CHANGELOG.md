# Changelog

All notable changes to this project are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions
follow the milestone scheme in the README (v0.1.0 = ingestion +
preprocessing working, v0.2.0 = MOFA+ integration working, v0.3.0 =
scoring package + tests passing, v0.4.0 = external validation complete,
v1.0.0 = full pipeline + report + docs complete and reproducible
end-to-end).

## [Unreleased]

### Added

- Repository scaffolding: full `src/oncocartograph` package layout
  (data_ingestion, preprocessing, integration, scoring, validation,
  drug_targets, reporting), each with a scoped docstring.
- `oncocartograph.config.Settings` — single pydantic-settings source of
  truth for paths, random seed, and log level.
- pytest test suite (currently covering `config`), ruff/black/mypy
  pre-commit hooks, and GitHub Actions CI (lint, type-check, test, docs
  presence check).
- MIT license, `.gitignore` tuned for omics data and workflow-engine run
  artifacts.
- `docs/methods.md`, `docs/data_sources.md`, and ADRs 0001-0003 recording
  the TNBC cohort definition thresholds, the Snakemake workflow engine
  choice, and the GSE96058 external validation cohort choice.
- `CITATION.cff` for citability.

### Pending

- `feat/data-ingestion`: TCGA-BRCA TNBC sub-cohort script and multi-omics
  ingestion.
- `feat/preprocessing`, `feat/mofa-integration`, `feat/scoring-package`,
  `feat/validation`, `feat/drug-target-scoring`, `feat/reporting`.
