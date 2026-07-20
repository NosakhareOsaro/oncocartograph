# Changelog

All notable changes to this project are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions
follow the milestone scheme in the README (v0.1.0 = ingestion +
preprocessing working, v0.2.0 = MOFA+ integration working, v0.3.0 =
scoring package + tests passing, v0.4.0 = external validation complete,
v1.0.0 = full pipeline + report + docs complete and reproducible
end-to-end). Patch versions (e.g. v0.3.1) mark real completed work
between named milestones that the original scheme didn't assign its
own number.

## [Unreleased]

## [0.4.0] - 2026-07-20

External validation against GSE96058 (SCAN-B), a real independent RNA-seq
TNBC cohort.

- `oncocartograph.validation` — GSE96058 series-matrix clinical parser,
  TNBC sub-cohort classifier (real result: N=143, 26 events), a
  streaming gene-symbol-filtered reader for the ~564MB expression
  matrix, a pre-registered direction-concordance replication test, and
  a Burstein et al. (2015) known-biology plausibility check
  (scope-reduced from full LAR/MES/BLIS/BLIA subtyping, confirmed with
  sign-off).
- Pre-registered the primary success/failure criterion — one-sided
  exact binomial test on direction concordance vs. chance, alpha=0.05 —
  before running the real analysis, specifically so it could not be
  redefined after seeing the result (ADR 0009).
- **Real result, reported honestly: the primary criterion FAILED.**
  45/109 (41.3%) TCGA candidates were direction-concordant in GSE96058,
  below the 50% chance rate (p=0.973). Consistent with
  `feat/scoring-package`'s already-documented 0/709 FDR-significant
  screen. The Burstein plausibility check passed 5/5 (AR, PTEN, CD274,
  PDCD1, CTLA4 all directionally consistent with the literature) — a
  real but separate result that does not offset the primary null
  finding. Documented as a limitation, not reframed.
- Per-candidate replication table persisted at
  `data/processed/gse96058_replication_table.csv`.
- 30 new tests, 100% coverage on all new modules, verified against the
  Python 3.11 target and the real live Open Targets API and real
  downloaded GEO data before each commit.

## [0.3.1] - 2026-07-20

Drug-target/druggability evidence working end-to-end against the real
709-candidate set.

- `oncocartograph.drug_targets` — Open Targets GraphQL client
  (`mapIds` gene-symbol resolution, `targets` tractability + UniProt
  accession lookup) and ChEMBL REST client (exact UniProt-accession
  target resolution — free-text search confirmed unreliable — and
  max clinical trial phase via the `mechanism` endpoint).
- 3-tier tractability scoring (Approved Drug=1.0, clinical-stage=0.66,
  any other evidence=0.33, none=0.0), combined with ChEMBL max phase via
  `max(tier_score, max_phase/4)` into the `tractability_score` the
  composite formula reads.
- Methylation candidates (CpG probe IDs) explicitly deferred — need an
  Illumina 450K manifest not in scope this iteration (ADR 0008).
- **Real re-scoring run:** 480/480 comparable candidates got real
  druggability evidence. Ranking changed substantially (Spearman
  ρ=0.656 vs. the druggability-absent ranking, only 4/20 previous top
  candidates remained in the top 20) — the mutation-derived `TP53`
  candidate ranked 7th overall with real ChEMBL evidence
  (`max_phase=3.0`), a candidate MOFA+ could never have surfaced.
- 38 new tests (173 total), 100% coverage on all new modules, verified
  against the Python 3.11 target and against the real live APIs (not
  just mocked HTTP) before each commit.

### Deferred (by design)

- `feat/validation`, `feat/reporting`.

## [0.3.0] - 2026-07-20

Composite biomarker scoring package working end-to-end against the real
preprocessed cohort -- this project's core novel, independently citable
contribution.

- `oncocartograph.scoring` — standalone package with zero cross-imports
  from the rest of `oncocartograph`, mechanically enforced by
  `tests/scoring/test_decoupling.py`.
- Two distinct, explicitly typed candidate-generation pathways feed the
  same composite score, per the v0.2.0 finding that MOFA+ contributes
  ~0% variance to the mutation view: `IntegrationEvidence` (MOFA+ factor
  loading, for RNA-seq/methylation/CNV) and `RecurrenceEvidence`
  (mutation prevalence + optional Fisher's-exact categorical
  association), both scored via identical `SurvivalEvidence` (univariate
  Cox PH, continuous or binary covariate).
- Cox PH on overall survival, not Fine-Gray competing-risks — confirmed
  as a data-availability constraint (no cause-of-death/recurrence coding
  in the real TCGA-BRCA clinical file), not a stylistic choice.
- `composite_biomarker_score` renormalizes weights over whichever
  evidence axes are present, so a mutation candidate is never penalized
  for lacking `IntegrationEvidence`. `DruggabilityEvidence` is schema-only
  here; `feat/drug-target-scoring` populates real values.
- Fixed a real bug caught by running the real 143-patient survival
  screen: lifelines can return a "successful" Cox fit with non-finite
  summary statistics (NaN, or an infinite confidence bound) for sparse
  binary covariates whose mutated subgroup has too few of the cohort's
  16 events. An initial NaN-only check caught 92/845 (~11%) of
  recurrence-filtered mutation genes; broadening to `np.isfinite` caught
  the true rate: 712/845 (84%).
- **Real scoring run (2026-07-20):** 709 candidates screened (576
  MOFA+-derived + 133 mutation-derived). 0 survived FDR correction — the
  expected honest outcome of screening hundreds of candidates against
  only 16 events, flagged as a limitation before this work package began.
  Documented as hypothesis-generating, not confirmatory.
- 30 new tests (135 total), 99% coverage, verified against the Python
  3.11 target.

## [0.2.0] - 2026-07-20

MOFA+ integration working end-to-end against the real preprocessed
cohort, not just synthetic fixtures.

- ADR 0006: `mofapy2`/`mofax` (pure Python) for MOFA+ training and
  interpretation, retracting an earlier assumption that this stage would
  need an R + MOFA2 environment — verified directly before writing any
  code.
- `oncocartograph.integration.mofa` — long-format view construction
  (verified to correctly handle both whole-patient absence from a view
  and scattered missing values within a view), training with confirmed
  hyperparameters (K=15, `scale_views=True`, `convergence_mode="slow"`,
  seed from config), and result extraction (factor values, variance
  explained).
- Fixed a real bug caught by the first end-to-end test: mofapy2 sorts
  view names alphabetically internally regardless of input order, so a
  likelihoods list built in dict-insertion order can silently misassign
  a likelihood to the wrong view. Fixed by sorting view names before
  building the likelihoods list.
- Fixed a real gap found while preparing this work package: the copy
  number view had ~60,623 genes with no feature-selection step (a 30x
  imbalance against the other views' 2,000-5,000). Added the same
  top-variable-gene selection already used for RNA-seq/methylation.
- **Real training run (2026-07-20):** copy number (2,000×142), RNA-seq
  (2,000×142), methylation (5,000×104), mutation (845×122); 143 patients
  total via view union, no forced complete-case subset. 12/15 factors
  clear a ≥2%-variance-explained screening threshold. Factor1 is
  overwhelmingly copy-number-driven (56.5% CNV variance, <1.1%
  elsewhere); the mutation view contributes essentially no variance to
  any factor (≤0.003%) — both reported as real findings in
  `docs/methods.md` §3.4/§8, with implications for the scoring work
  package (mutation biomarkers will need direct statistics, not MOFA+
  loadings).
- 105 tests total, 100% coverage, verified against the Python 3.11
  target throughout.

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
