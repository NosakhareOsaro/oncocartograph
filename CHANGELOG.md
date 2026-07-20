# Changelog

All notable changes to this project are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions
follow the milestone scheme in the README (v0.1.0 = ingestion +
preprocessing working, v0.2.0 = MOFA+ integration working, v0.3.0 =
scoring package + tests passing, v0.4.0 = external validation complete,
v1.0.0 = full pipeline + report + docs complete and reproducible
end-to-end).

## [Unreleased]

### Deferred (by design)

- `feat/mofa-integration`, `feat/scoring-package`, `feat/validation`,
  `feat/drug-target-scoring`, `feat/reporting`.

## [0.1.0] - 2026-07-20

Ingestion and preprocessing working end-to-end against real TCGA-BRCA
data, not just synthetic fixtures.

### Data ingestion

- ADR 0004: direct GDC REST API client instead of TCGAbiolinks/R.
- `oncocartograph.data_ingestion.gdc_client.GDCClient` — typed GDC REST
  API client (files/cases query with pagination, single-file download,
  retry-with-backoff on transient failures).
- `oncocartograph.data_ingestion.clinical` — BCR Biotab clinical
  supplement parser and receptor-status column extraction.
- `oncocartograph.data_ingestion.tnbc_cohort` — TNBC cohort
  classification implementing ADR 0001's exact rules, producing a full
  per-patient audit table (raw values in, include/exclude decision +
  reason out).
- `oncocartograph.data_ingestion.provenance` — SHA-256-checksummed
  provenance sidecar records for every downloaded artifact.
- `oncocartograph.data_ingestion.omics_ingestion` — per-omic GDC file
  filters (RNA-seq STAR-Counts, 450K methylation beta values, gene-level
  copy number, Masked Somatic Mutation MAF) and download orchestration.
- **Live pull (2026-07-20, GDC Data Release 45.0):** 143/1,097 TCGA-BRCA
  patients (13.0%) classify as TNBC. Ingested 158 RNA-seq + 348
  methylation + 426 copy number + 126 mutation files for the cohort
  (~5.4 GB), each with a checksummed provenance sidecar.
- Fixed a real bug found by the live pull: 3 patients have
  `her2_fish_status="Equivocal"` (FISH performed but inconclusive), but
  the exclusion message wrongly implied no FISH was attempted.
- Fixed the methylation filter to exclude raw `.idat` files (wasted
  ~1.8GB on the first live pull before the fix).

### Preprocessing

- `oncocartograph.preprocessing.sample_manifest` — resolves every omic
  layer to one Primary Tumor sample per patient, handling both
  single-sample files (RNA-seq/methylation) and paired tumor+normal
  files (copy number/mutation callers) correctly.
- `oncocartograph.preprocessing.rna_seq` — CPM-based low-expression
  filtering, DESeq2-style size-factor normalization + VST via
  `pydeseq2`, top-2,000-variable-gene selection.
- `oncocartograph.preprocessing.copy_number` — workflow priority
  resolution (ASCAT3 > ASCAT2 > ABSOLUTE LiftOver > AscatNGS, confirmed
  via live GDC metadata, not assumed) and relative-to-diploid log2
  transform. ADR 0005 corrects the original "GISTIC2 thresholded calls"
  assumption — real data is absolute integer copy number.
- `oncocartograph.preprocessing.methylation` — beta value loading,
  missingness filtering, beta→M-value transform, top-5,000-variable-probe
  selection.
- `oncocartograph.preprocessing.mutation` — non-synonymous variant
  allow-list filtering, binary gene×patient matrix, ≥3-patient recurrence
  filter.
- Real-data sanity checks (not just synthetic-fixture tests) for every
  module: real cohort N, real workflow distribution, real TP53 mutation
  frequency (7/10 patients, matching TNBC's known ~80% rate), real 450K
  probe count (486,427).
- 96 tests total across data ingestion + preprocessing, 100% coverage,
  verified against the actual Python 3.11 target via the project Docker
  image throughout (not just the local dev environment).
