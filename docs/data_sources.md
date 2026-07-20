# Data sources

Every dataset used in OncoCartograph is listed here with its exact
accession/version, access method, download date, and license. Entries are
added at the point a work package actually ingests the data — this file
never lists a dataset as "in use" before the corresponding ingestion code
exists and has been run.

## Primary cohort: TCGA-BRCA

| Field | Value |
|---|---|
| Project | TCGA-BRCA (Breast Invasive Carcinoma) |
| Access | NCI Genomic Data Commons (GDC) REST API, open-tier harmonized data, queried directly (see ADR 0004 — no TCGAbiolinks/R dependency) |
| GDC data release at time of pull | Data Release 45.0 (2025-12-04), API commit `8f7c2a51ab0084b216ad1b62a3fae8b945439c53` |
| Clinical / receptor status source file | `nationwidechildrens.org_clinical_patient_brca.txt` (BCR Biotab, `data_type=Clinical Supplement`), GDC file UUID `8162d394-8b64-4da2-9f5b-d164c54b9608`, 1,097 patient rows |
| **TNBC cohort N (live pull, 2026-07-20)** | **143** of 1,097 TCGA-BRCA patients — see `docs/methods.md` §1.2 for the full exclusion breakdown and `data/processed/tnbc_cohort_audit.csv` for the per-patient audit table (gitignored; regenerate via `oncocartograph.data_ingestion.tnbc_cohort`) |
| RNA-seq | GDC harmonized STAR-Counts workflow — 158 files, 763 MB, downloaded 2026-07-20 |
| DNA methylation | Illumina Infinium HumanMethylation450 (450K) — 348 files, 3.2 GB, downloaded 2026-07-20 |
| Copy number | GDC Gene Level Copy Number — 426 files, 1.4 GB, downloaded 2026-07-20 |
| Somatic mutation | GDC Masked Somatic Mutation (MC3-derived) — 126 files, 6.0 MB, downloaded 2026-07-20 |
| File counts exceeding cohort N | Expected: some patients have multiple samples/aliquots per omic layer (e.g. tumor + matched normal, multiple variant-calling pipelines for mutation data). Resolved during `feat/preprocessing` dedup, not during ingestion. |
| Per-file provenance | One `<filename>.provenance.json` sidecar per downloaded file under `data/raw/{rna_seq,methylation,copy_number,mutation,clinical}/` (gitignored), each with a SHA-256 checksum, the exact GDC filter used, and download timestamp — this table is the human-readable summary, the sidecars are authoritative. |
| License | GDC open-tier data usage policy (no additional restriction beyond standard TCGA open access terms) |

## External validation cohort: GEO GSE96058 (SCAN-B)

| Field | Value |
|---|---|
| Accession | [GSE96058](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96058) |
| Description | Sweden Cancerome Analysis Network—Breast (SCAN-B) RNA-seq cohort |
| Platform | RNA-seq (Illumina HiSeq 2000 / NextSeq 500) |
| N | 3,409 total samples (3,273 primary + 136 technical replicates) across two platform-specific series-matrix files; TNBC subset via histopathology `er status`/`pgr status`/`her2 status` fields (all-negative, `"NA"` excluded), same rule as Section 1.2 of `methods.md` — real result: **N=143, 26 events** |
| Survival | Overall survival (`overall survival days`/`overall survival event` fields), median follow-up ~52 months |
| Associated publications | Saal et al. 2015, *Genome Med* 7:20 (cohort description); Brueffer et al. 2018, *JCO Precis Oncol* (biomarker classifiers) |
| Why chosen over alternatives | See [`docs/adr/0003-external-validation-cohort.md`](adr/0003-external-validation-cohort.md) |
| License | GEO public data, no restriction beyond standard NCBI terms |
| Files used | Series-matrix clinical files: `GSE96058-GPL11154_series_matrix.txt.gz` (HiSeq 2000), `GSE96058-GPL18573_series_matrix.txt.gz` (NextSeq 500). Expression matrix: `GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz` (~564MB; gene-symbol-indexed, log2(FPKM+0.1)-transformed), from `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96058/suppl/` |
| Download date | 2026-07-20 |

## Published-hit reproduction benchmark

| Field | Value |
|---|---|
| Study | Burstein et al. 2015, *Clin Cancer Res* 21(7):1688-1698, PMID 25208879 |
| Reported hits used for comparison | Subtype-specific druggable targets: androgen receptor (AR) in the LAR subtype; recurrent PTEN/PI3K-pathway alterations across subtypes; immune-checkpoint gene expression in the BLIA subtype |
| Cohort relationship to TCGA-BRCA | Independent (Baylor College of Medicine discovery/validation cohorts plus seven public TNBC datasets) — not a subset of TCGA-BRCA, so reproducing these hits from our TCGA-BRCA-derived pipeline is a genuine cross-cohort replication check |
| License / access | Published paper; no raw data redistribution needed since we compare our own TCGA-derived findings against their *reported* gene/pathway-level conclusions |
| Scope actually implemented | Reduced, with sign-off, from full LAR/MES/BLIS/BLIA subtype reproduction to a five-gene (AR, PTEN, CD274, PDCD1, CTLA4) directional-plausibility check against real GSE96058 survival evidence — see [`docs/adr/0009-external-replication-methodology-and-result.md`](adr/0009-external-replication-methodology-and-result.md). Real result: 5/5 markers plausible. |

## Drug/target evidence

| Field | Value |
|---|---|
| Open Targets | Public GraphQL API, `https://api.platform.opentargets.org/api/v4/graphql` — no API key required for the request volumes this project needs |
| Open Targets queries used | `mapIds` (gene symbol → Ensembl ID, mutation candidates only); `targets(ensemblIds: ...)` for `tractability` (28 boolean buckets, 4 modalities) and `proteinIds` (canonical `uniprot_swissprot` accession) |
| ChEMBL | Public REST API, `https://www.ebi.ac.uk/chembl/api/data` |
| ChEMBL queries used | `target.json?target_components__accession__in=...&target_type=SINGLE PROTEIN` (exact UniProt accession → ChEMBL target ID — free-text search confirmed unreliable, e.g. "TP53" search returns "TP53-binding protein 1" first); `mechanism.json?target_chembl_id__in=...` (max `max_phase` across all mechanism records per target) |
| Coverage | RNA-seq/copy-number (Ensembl ID) and mutation (gene symbol) candidates only. Methylation (CpG probe ID) candidates are out of scope this iteration — see [`docs/adr/0008-druggability-evidence-sources.md`](adr/0008-druggability-evidence-sources.md) |
| Query date | 2026-07-20 (`feat/drug-target-scoring`). API responses can change over time; re-running the pipeline will reflect the Open Targets/ChEMBL data current at that time, not this snapshot. |

## Provenance logging convention

Every ingestion script writes a provenance record (source, query
parameters, accession/version, download timestamp, and — where
applicable — file checksums) alongside its output, so any intermediate
file in this pipeline can be traced back to the exact external query that
produced it. This table is the human-readable summary; the machine-readable
provenance logs are the authoritative source.
