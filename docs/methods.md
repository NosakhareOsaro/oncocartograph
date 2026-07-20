# Methods

> **Status:** this document is written incrementally as each pipeline stage
> is implemented. Sections for stages not yet built are marked `_Pending_`
> rather than left silently absent, so the document's completeness matches
> the pipeline's actual state at every commit.

## 1. Cohort definition

### 1.1 Base cohort

All samples are drawn from **TCGA-BRCA** (The Cancer Genome Atlas Breast
Invasive Carcinoma cohort), accessed via the NCI Genomic Data Commons (GDC)
harmonized data and the TCGA Pan-Cancer Clinical Data Resource (TCGA-CDR;
Liu et al. 2018, *Cell* 173(2):400-416.e11, PMID 29625055). Exact accession
numbers, file UUIDs, and download dates are logged in
[`data_sources.md`](data_sources.md) and, per-run, in the ingestion script's
provenance output.

### 1.2 Triple-negative sub-cohort definition

A sample is classified **TNBC** if and only if all three of the following
hold, using the IHC/FISH fields recorded in the TCGA-BRCA clinical
supplement:

| Marker | Negative definition | Field(s) used | Guideline cited |
|---|---|---|---|
| ER | IHC nuclear staining < 1% of tumour cells | `er_status_by_ihc == "Negative"` | Hammond et al. 2010, *Arch Pathol Lab Med* 134:907-922 |
| PR | IHC nuclear staining < 1% of tumour cells | `pr_status_by_ihc == "Negative"` | Hammond et al. 2010 (as above); reaffirmed in the 2020 update, Allison et al., *J Clin Oncol* 38:1346-1366 |
| HER2 | IHC 0/1+ (negative), **or** IHC 2+ (equivocal) with reflex FISH HER2/CEP17 ratio < 2.0 | `her2_status_by_ihc`, `her2_fish_status` | Wolff et al. 2013/2018 ASCO/CAP HER2 testing guideline update, *J Clin Oncol* / *Arch Pathol Lab Med* |

**Exclusion rule (not imputation):** any patient with a missing status call
for ER, PR, or HER2, or an IHC-2+ (equivocal) HER2 call with **no** recorded
FISH follow-up, is excluded from the TNBC cohort and logged, with reason, in
the cohort-definition script's audit table output. This is a deliberate
choice: many published TNBC analyses do not state their exact ER/PR/HER2
cutoffs or how they handled equivocal/missing calls, which makes their
cohort definitions hard to audit or reproduce. OncoCartograph's ingestion
script instead emits, for every TCGA-BRCA patient, the raw field values used
and the resulting include/exclude decision, so the final cohort N is a
verifiable consequence of documented rules rather than an opaque filter.

**Sample-level rules:** primary tumour samples only (GDC sample type code
`01`), one sample per patient (classification operates at the patient
level on the clinical supplement; per-omic sample deduplication happens in
`feat/preprocessing`).

**Empirical result (live pull, 2026-07-20, GDC Data Release 45.0):** of
1,097 TCGA-BRCA patients in `nationwidechildrens.org_clinical_patient_brca.txt`,
**143 (13.0%) were classified as TNBC.** Exclusion breakdown for the
remaining 954:

| Reason category | N |
|---|---|
| Receptor positive (ER and/or PR and/or HER2) | 737 |
| ER and/or PR and/or HER2 IHC status missing/indeterminate (`[Not Evaluated]`, `[Not Available]`, `Indeterminate`) | 196 |
| HER2 IHC equivocal, not resolved by FISH | 21 |
| **Total excluded** | **954** |

These three categories are mutually exclusive in this dataset (they sum
exactly to 954) because indeterminacy always takes precedence over a
receptor-positive verdict when both would otherwise apply (see
`classify_patient`'s precedence rule, tested in
`tests/data_ingestion/test_tnbc_cohort.py`). The authoritative per-patient
breakdown, including the full raw field values behind every decision, is
`data/processed/tnbc_cohort_audit.csv` (gitignored; regenerable from the
source clinical file via
`oncocartograph.data_ingestion.tnbc_cohort.build_tnbc_cohort_audit`).

This confirmed a real accuracy issue in the initial implementation: 3
patients have `her2_fish_status="Equivocal"` (a FISH result *was* recorded,
it just didn't resolve the case) rather than a missing value, and the
original exclusion message wrongly implied no FISH had been attempted.
Fixed before this pull's numbers were finalized (see git history for
`oncocartograph.data_ingestion.tnbc_cohort`) — the classification outcome
was unaffected, only the audit message's accuracy.

### 1.3 Rationale for this specific rule set

See [`docs/adr/0001-tnbc-cohort-definition.md`](adr/0001-tnbc-cohort-definition.md)
for the full decision record, including the alternative (treating
IHC-equivocal-without-FISH as negative) that was considered and rejected.

## 2. Data preprocessing

_Pending — will document per-omic QC thresholds (RNA-seq low-count
filtering, methylation probe filtering/masking, CNV segmentation handling,
MAF variant filtering) and normalisation methods once
`feat/preprocessing` lands._

## 3. Multi-omics integration (MOFA+)

_Pending — will document factor count justification, convergence
criteria, random seed handling, and factor-to-biology interpretation
approach once `feat/mofa-integration` lands. See
`docs/adr/0002-workflow-engine-choice.md` for the orchestration decision
made ahead of this stage._

## 4. Composite biomarker scoring

_Pending — will document the survival model (Cox proportional hazards vs.
Fine-Gray competing-risks, chosen based on observed competing-event
structure in the TNBC cohort), multiple-testing correction method,
p-value/effect-size thresholds, and the druggability scoring formula once
`feat/scoring-package` lands._

## 5. External validation

_Pending — will document the GSE96058 (SCAN-B) cohort's TNBC subsetting
criteria (matched to Section 1.2 where the available fields allow), and the
comparison methodology against Burstein et al. 2015 (PMID 25208879)
subtype-specific druggable target calls, once `feat/validation` lands._

## 6. Drug-target evidence and prioritisation

_Pending — will document the Open Targets and ChEMBL query strategy and
tractability/bioactivity evidence weighting once
`feat/drug-target-scoring` lands._

## 7. Software and versions

Recorded per-release in `CHANGELOG.md` and pinned in `pyproject.toml` /
lockfiles. Full version table will be added once the pipeline has a first
end-to-end run.

## 8. Limitations

_Pending final version, but noted here early so it is not forgotten:_
TCGA-BRCA IHC calls were made across many contributing institutions without
fully centralised scoring, so some inter-site variability in ER/PR/HER2
calls is expected and cannot be corrected for retrospectively. The GSE96058
external cohort is Swedish and RNA-seq-only (no methylation/CNV/mutation
validation), so external validation in this project is necessarily
expression-centric. The Burstein et al. 2015 benchmark cohort is not
population-matched to TCGA-BRCA. These limitations will be expanded with
concrete impact assessments once results exist to discuss.
