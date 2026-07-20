# ADR 0005: Copy number workflow priority and relative-to-diploid transform

## Status

Accepted (2026-07-20)

## Context

The original data plan (and `docs/data_sources.md` as first written) assumed
TCGA-BRCA copy number would be **GISTIC2 thresholded calls** (categorical
-2/-1/0/1/2). Inspecting the actual files downloaded during the live pull
(`feat/data-ingestion`) showed this was wrong: GDC's current harmonized
pipeline for TCGA-BRCA gene-level copy number reports **absolute integer
total copy number per gene** (2 = diploid, 0 = homozygous deletion, etc.),
produced by up to four different calling workflows per patient. A live GDC
metadata query (`analysis.workflow_type`) confirmed the real distribution
across the 143-patient TNBC cohort: ASCAT3 (138 files), ASCAT2 (141 files
— the *older* version of the same algorithm, not a duplicate of ASCAT3),
ABSOLUTE LiftOver (131 files, a different algorithm entirely), and
AscatNGS (16 files, a WGS-specific ASCAT variant for cases without SNP
array data).

Two decisions were needed: which workflow to prefer per patient, and how
to transform absolute copy number into a value suitable for a Gaussian
MOFA+ view.

## Decision

**Workflow priority:** ASCAT3 > ASCAT2 > ABSOLUTE LiftOver > AscatNGS.
Implemented as `oncocartograph.preprocessing.copy_number.CNV_WORKFLOW_PRIORITY`,
applied via a `workflow_priority_tie_break` factory plugged into
`sample_manifest.resolve_primary_tumor_files`'s generic tie-break
mechanism, rather than a copy-number-specific resolution path.

**Transform:** relative-to-diploid log2,
`log2((copy_number + 1) / (DIPLOID_COPY_NUMBER + 1))` with
`DIPLOID_COPY_NUMBER = 2`. Adding 1 to *both* the numerator and the
diploid reference (not just the numerator) is deliberate: it makes
diploid (CN=2) map to exactly 0 while keeping CN=0 (homozygous deletion)
finite instead of `log2(0) = -inf`.

## Rationale

- Preferring ASCAT3 over ASCAT2 is a same-algorithm-family newest-version
  preference — more defensible than jumping straight to a different
  algorithm (ABSOLUTE) when a newer version of the same caller is
  available.
- AscatNGS is ranked last because it's WGS-specific and only exists for
  the small subset of cases lacking the SNP-array-based pipelines; it's a
  necessary fallback, not a preferred source.
- The naive `log2((CN+1)/2)` transform (tried first, caught before
  committing) maps diploid to 0.585, not 0 — verified numerically. The
  corrected formula was checked against `DIPLOID_COPY_NUMBER=2 -> 0` and
  `CN=0 -> finite` before being used anywhere.

## Alternatives considered

**Attempt to reconstruct GISTIC2-style categorical calls** by thresholding
the absolute copy number against an estimated tumor ploidy. Rejected for
this iteration: correcting for tumor ploidy properly requires a
purity/ploidy estimate per sample (which ABSOLUTE's own output could
provide, but only for the subset of patients with an ABSOLUTE LiftOver
file), adding real scope for a benefit (categorical vs. continuous
Gaussian view) that MOFA+ doesn't strictly need. Documented as a
limitation in `docs/methods.md` rather than engineered around.

## Consequences

- `docs/data_sources.md`'s original "GISTIC2 thresholded calls" wording is
  corrected to describe what's actually used.
- The CNV Gaussian view in MOFA+ does not account for whole-genome
  doubling or tumor purity; genes in samples with unusual ploidy will have
  systematically shifted relative copy number. Noted in
  `docs/methods.md` §8 (Limitations).
