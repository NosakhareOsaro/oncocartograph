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
- ADR 0004: direct GDC REST API client instead of TCGAbiolinks/R, keeping
  the ingestion path pure Python and testable without an R runtime.
- `oncocartograph.data_ingestion.gdc_client.GDCClient` — typed GDC REST API
  client (files/cases query with pagination, single-file download,
  retry-with-backoff on transient failures).
- `oncocartograph.data_ingestion.clinical` — BCR Biotab clinical supplement
  parser and receptor-status column extraction.
- `oncocartograph.data_ingestion.tnbc_cohort` — TNBC cohort classification
  implementing ADR 0001's exact rules, producing a full per-patient audit
  table (raw values in, include/exclude decision + reason out).
- `oncocartograph.data_ingestion.provenance` — SHA-256-checksummed
  provenance sidecar records for every downloaded artifact.
- `oncocartograph.data_ingestion.omics_ingestion` — per-omic GDC file
  filters (RNA-seq STAR-Counts, 450K methylation, GISTIC2 gene-level copy
  number, MC3 Masked Somatic Mutation MAF) and download orchestration.
- 45 tests across the ingestion module, 100% coverage, verified against
  the actual Python 3.11 target via the project Docker image.

### Deferred (by design, this iteration)

- Live pull against the real GDC API (real cohort N, real downloaded
  files) — the ingestion code above is built and unit-tested against
  mocked/synthetic fixtures only. Running it live against GDC, and
  updating `docs/data_sources.md`/`docs/methods.md` §1 with the resulting
  real numbers, is a deliberate next step, not an oversight.
- `feat/preprocessing`, `feat/mofa-integration`, `feat/scoring-package`,
  `feat/validation`, `feat/drug-target-scoring`, `feat/reporting`.
