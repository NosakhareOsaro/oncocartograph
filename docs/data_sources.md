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
| Access | NCI Genomic Data Commons (GDC), open-tier harmonized data, via TCGAbiolinks |
| RNA-seq | GDC harmonized STAR-Counts workflow |
| DNA methylation | Illumina Infinium HumanMethylation450 (450K) |
| Copy number | GISTIC2 thresholded calls (source: Broad Firehose or GDC — final choice logged at ingestion time) |
| Somatic mutation | MC3 public MAF |
| Clinical / receptor status | TCGA-BRCA clinical supplement + TCGA-CDR (Liu et al. 2018, PMID 29625055) |
| License | GDC open-tier data usage policy (no additional restriction beyond standard TCGA open access terms) |
| Exact file UUIDs / query parameters / download date | _To be recorded here by the ingestion script's provenance log once `feat/data-ingestion` runs — this table will be replaced with the actual GDC query manifest, not summarised by hand._ |

## External validation cohort: GEO GSE96058 (SCAN-B)

| Field | Value |
|---|---|
| Accession | [GSE96058](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96058) |
| Description | Sweden Cancerome Analysis Network—Breast (SCAN-B) RNA-seq cohort |
| Platform | RNA-seq (Illumina HiSeq 2000 / NextSeq 500) |
| N | ~3,273 samples total; TNBC subset defined by the same rules as Section 1.2 of `methods.md`, applied to whatever ER/PR/HER2 annotation fields GEO provides for this series (exact field mapping to be documented at ingestion time, since GEO metadata field names differ from GDC's) |
| Survival | Overall survival, median follow-up ~52 months |
| Associated publications | Saal et al. 2015, *Genome Med* 7:20 (cohort description); Brueffer et al. 2018, *JCO Precis Oncol* (biomarker classifiers) |
| Why chosen over alternatives | See [`docs/adr/0003-external-validation-cohort.md`](adr/0003-external-validation-cohort.md) |
| License | GEO public data, no restriction beyond standard NCBI terms |
| Download date | _To be recorded once `feat/validation` ingests this series._ |

## Published-hit reproduction benchmark

| Field | Value |
|---|---|
| Study | Burstein et al. 2015, *Clin Cancer Res* 21(7):1688-1698, PMID 25208879 |
| Reported hits used for comparison | Subtype-specific druggable targets: androgen receptor (AR) in the LAR subtype; recurrent PTEN/PI3K-pathway alterations across subtypes; immune-checkpoint gene expression in the BLIA subtype |
| Cohort relationship to TCGA-BRCA | Independent (Baylor College of Medicine discovery/validation cohorts plus seven public TNBC datasets) — not a subset of TCGA-BRCA, so reproducing these hits from our TCGA-BRCA-derived pipeline is a genuine cross-cohort replication check |
| License / access | Published paper; no raw data redistribution needed since we compare our own TCGA-derived findings against their *reported* gene/pathway-level conclusions |

## Drug/target evidence

| Field | Value |
|---|---|
| Open Targets | Public GraphQL API, `https://api.platform.opentargets.org/api/v4/graphql` — no API key required for the request volumes this project needs |
| ChEMBL | Public REST API, `https://www.ebi.ac.uk/chembl/api/data` |
| Exact query parameters and response snapshot dates | _To be recorded once `feat/drug-target-scoring` lands, since API responses can change over time and every score must be traceable to a specific query date._ |

## Provenance logging convention

Every ingestion script writes a provenance record (source, query
parameters, accession/version, download timestamp, and — where
applicable — file checksums) alongside its output, so any intermediate
file in this pipeline can be traced back to the exact external query that
produced it. This table is the human-readable summary; the machine-readable
provenance logs are the authoritative source.
