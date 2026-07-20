# ADR 0004: Direct GDC REST API client instead of TCGAbiolinks/R

## Status

Accepted (2026-07-20)

## Context

TCGAbiolinks (the tool named in the original data plan) is an R/Bioconductor
package. This project's engineering standard is a typed, tested Python
3.11 codebase (ruff/black/mypy, >85% coverage target), and the rest of the
pipeline besides the MOFA+ step (ADR 0002 already accepts an R dependency
there, scoped narrowly to that one stage) is pure Python. Bridging
TCGAbiolinks into the ingestion path would mean either an R subprocess
step or an in-process rpy2 bridge for functionality that is, underneath,
just calls to the public GDC REST API.

## Decision

Query the [GDC REST API](https://api.gdc.cancer.gov) directly from Python
(`files`, `cases`, `data` endpoints) via a typed `GDCClient`, rather than
going through TCGAbiolinks.

## Rationale

- **No R dependency in the ingestion path.** Keeps ingestion testable with
  ordinary mocked HTTP responses in pytest, with no subprocess or R runtime
  needed in CI for this stage.
- **The GDC REST API is public, stable, and well-documented** — TCGAbiolinks
  itself is a wrapper around exactly these endpoints, so nothing is lost in
  terms of what data is reachable.
- **Consistent with the "typed Python everywhere except MOFA+" architecture**
  already established in ADR 0002, rather than introducing a second R
  dependency surface for a different reason.

## Alternatives considered

- **R subprocess bridge** (TCGAbiolinks script invoked from Python, writing
  flat files Python then reads): rejected — adds a second language and an
  untestable-in-pytest subprocess boundary to the ingestion path for
  functionality reachable directly over HTTP.
- **rpy2 in-process bridge**: rejected — most fragile option
  (R/Python object marshalling), hardest to keep working across
  environments (local venv, CI, Docker) long-term.

## Consequences

- Clinical ER/PR/HER2/FISH status fields are sourced from the GDC
  "Clinical Supplement" BCR Biotab file (tab-delimited, fixed 3-row header)
  rather than via `TCGAbiolinks::GDCquery_clinic()`. The exact column
  layout of this file will be validated against a real downloaded file
  during the first live pull (tracked as a known open item, not assumed
  correct in advance) and `docs/data_sources.md` updated accordingly.
- If a future need arises for a GDC dataset not cleanly reachable via the
  REST API, that would be a new ADR, not a silent reintroduction of
  TCGAbiolinks.
