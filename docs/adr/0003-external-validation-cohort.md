# ADR 0003: External validation cohort — GSE96058 (SCAN-B)

## Status

Accepted (2026-07-20)

## Context

The project requires an independent external cohort to validate biomarker
scores derived from the TCGA-BRCA TNBC pipeline, per the "reproducible
before novel" design goal. Candidates considered: GSE96058 (SCAN-B) and
GSE58812 (Jézéquel et al. 2015).

## Decision

Use **GSE96058 (SCAN-B)** as the primary external validation cohort.

## Rationale

| | GSE96058 (chosen) | GSE58812 |
|---|---|---|
| Platform | RNA-seq | Affymetrix U133 Plus 2.0 microarray |
| N (TNBC subset) | ~140-150 (from ~3,273 total, subset by receptor annotation) | 107 (TNBC-exclusive design) |
| Survival endpoint | Overall survival, 52mo median follow-up | Event-free survival only |
| Institution(s) | Multi-site (Sweden, SCAN-B network) | Single institution (Nantes, France) |

The decisive factor is **platform match**: TCGA-BRCA RNA-seq data is
processed via the GDC harmonized STAR-Counts workflow, and validating
expression-derived MOFA+ factors and composite scores against a microarray
cohort (GSE58812) would introduce a cross-platform normalization confound
that is hard to distinguish from genuine biological non-replication.
GSE96058's RNA-seq platform avoids that confound. Its larger N and
multi-site design also give more statistical power for the TNBC subset
than GSE58812's 107 samples.

## Alternatives considered

**GSE58812.** Attractive because it is TNBC-exclusive by design (no
subsetting needed) and single-institution (less batch heterogeneity), but
rejected as the *primary* cohort due to the microarray/RNA-seq platform
mismatch and EFS-only (not OS) survival endpoint. Retained as a candidate
**secondary/sensitivity cohort** if time permits after `feat/validation` —
not committed to as required scope, to avoid inflating this work package.

## Consequences

- `feat/validation` must map GSE96058's GEO-provided receptor status
  metadata fields to the same ER/PR/HER2 negativity rules as
  `docs/adr/0001-tnbc-cohort-definition.md`; exact field names differ from
  GDC's and will be documented in `docs/data_sources.md` once ingested.
- Any expression-based score comparison between TCGA-BRCA and GSE96058
  must account for the two cohorts' independent RNA-seq processing
  pipelines (this project does not attempt to re-process GSE96058 raw
  reads through the GDC pipeline).
