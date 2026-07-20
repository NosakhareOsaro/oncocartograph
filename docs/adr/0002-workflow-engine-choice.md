# ADR 0002: Workflow engine — Snakemake over Nextflow

## Status

Accepted (2026-07-20)

## Context

The pipeline needs a workflow engine to orchestrate ingestion,
preprocessing, MOFA+ integration, scoring, validation, and reporting as a
reproducible DAG, with per-step environment isolation (MOFA+ itself is an
R/Python hybrid tool, so at least one step needs an R environment alongside
otherwise Python-native code).

The two realistic candidates were **Snakemake** and **Nextflow (DSL2)**.

## Decision

Use **Snakemake**.

## Rationale

- **Python-native rules.** The rest of the codebase is Python 3.11 with
  type hints and pydantic-settings config; Snakemake rules can call Python
  functions and scripts directly and share the same `config.yaml`/Settings
  conventions, whereas Nextflow's Groovy-based DSL2 would introduce a
  second language and a second configuration surface for a solo-maintained
  project.
- **Per-rule conda/mamba environments** are first-class in Snakemake,
  which is exactly what's needed to isolate the R-based MOFA+ step from the
  Python steps without a heavier container-per-step setup.
- **Scale profile matches the project.** Nextflow's main advantages —
  cloud/HPC executor abstraction, nf-core ecosystem conventions — matter
  most at production/cloud scale. This project runs on a single workstation
  or a single HPC node against a cohort of ~150-200 TNBC samples; Nextflow's
  additional complexity would not be earning its keep here.
- **DAG visualization** (`snakemake --dag`) gives the pipeline DAG diagram
  referenced in the README directly from the workflow definition, rather
  than a hand-maintained mermaid diagram that can drift from the real
  pipeline.

## Alternatives considered

**Nextflow DSL2.** Rejected primarily due to the added Groovy/DSL2 learning
surface and configuration duplication for a solo maintainer, and because
its cloud-scale executor strengths are not needed at this project's scale.
It remains a reasonable choice if this pipeline is ever productionised
across many cohorts/institutions — noted here so the tradeoff isn't
forgotten if that need arises.

## Consequences

- `workflows/Snakefile` is the single source of truth for pipeline
  orchestration; nothing load-bearing lives only in a notebook.
- The MOFA+ integration rule will declare its own conda environment
  (R + MOFA2) separate from the Python steps' environment.
