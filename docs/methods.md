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

### 2.1 Sample resolution (all omic layers)

Files ingested in `feat/data-ingestion` can include more than one file per
patient per layer: matched normal or metastatic samples alongside the
primary tumor (RNA-seq, methylation), and — for copy number and mutation,
which are produced by tumor-normal paired callers — a germline reference
sample listed alongside the tumor sample against the *same* file, not a
separate one. `oncocartograph.preprocessing.sample_manifest` resolves
every layer to one **Primary Tumor** sample per patient (membership check
on the file's sample types, not exclusivity), with a deterministic
tie-break (smallest file UUID) for the rare patient with more than one
qualifying file. Patients with no Primary Tumor file for a given layer are
simply absent from that layer's matrix — this project does not force a
smaller complete-case cohort across all four omics, relying instead on
MOFA+'s native support for partially-overlapping views.

Confirmed patient counts after resolution (live pull, 2026-07-20, of the
143-patient TNBC cohort):

| Layer | Resolved patients | Note |
|---|---|---|
| RNA-seq | 142 | 1 patient has no RNA-seq file at all |
| Methylation | ~106 | Largest gap — many TNBC patients lack a 450K profile in TCGA-BRCA |
| Copy number | 142 | See §2.3 for the additional workflow-resolution step |
| Mutation | ~122 | Some patients lack WXS/mutation calls |

### 2.2 RNA-seq (`oncocartograph.preprocessing.rna_seq`)

- **Input:** GDC STAR-Counts `augmented_star_gene_counts.tsv`, raw
  (`unstranded`) counts column. STAR's four alignment-summary rows
  (`N_unmapped`, `N_multimapping`, `N_noFeature`, `N_ambiguous`) are
  removed before use.
- **Low-expression filter:** genes must reach ≥1 CPM in ≥10% of samples
  (edgeR `filterByExpr`-style convention) to be retained. On a 10-patient
  real-data check this reduced 60,660 to 25,186 genes.
- **Normalization:** DESeq2-style median-of-ratios size-factor
  normalization and variance-stabilizing transformation (VST), via
  `pydeseq2` (`design="~1"`, no differential-expression contrast — MOFA+
  discovers structure directly from the resulting continuous values).
  Chosen over a simpler log2(TPM+1) approach because size factors correct
  for library-composition differences more robustly than TPM alone.
- **Feature selection:** top 2,000 genes by variance across patients, for
  MOFA+ tractability.

### 2.3 Copy number (`oncocartograph.preprocessing.copy_number`)

- **Input:** GDC gene-level copy number files. **Correction to the
  original data plan:** these report absolute integer total copy number
  per gene (2 = diploid), not GISTIC2 thresholded categorical calls as
  first assumed — see
  [`docs/adr/0005-copy-number-workflow-and-transform.md`](adr/0005-copy-number-workflow-and-transform.md).
- **Workflow resolution:** up to four calling workflows exist per
  patient (ASCAT3, ASCAT2, ABSOLUTE LiftOver, AscatNGS — confirmed via
  live GDC metadata query, not assumed from filenames). Preference order:
  ASCAT3 > ASCAT2 > ABSOLUTE LiftOver > AscatNGS. On the real 143-patient
  cohort this resolved 142 patients: 138 via ASCAT3 directly, 3 falling
  back to ASCAT2, 1 to ABSOLUTE LiftOver.
- **Transform:** relative-to-diploid log2,
  `log2((copy_number + 1) / (2 + 1))` — zeroes the diploid baseline
  exactly while keeping homozygous deletion (CN=0) finite. Does not
  correct for tumor purity/ploidy (§8 Limitations).

### 2.4 DNA methylation (`oncocartograph.preprocessing.methylation`)

- **Input:** GDC SeSAMe `level3betas.txt` processed beta values only —
  the ingestion filter was corrected during this work package to exclude
  raw `.idat` intensity files it had also been pulling (see git history
  for `feat(data-ingestion)`). 486,427 probes per patient, matching the
  450K array; ~14% already `NA` from SeSAMe's own probe masking, so no
  additional cross-reactive-probe blacklist is applied here.
- **Missingness filter:** probes missing in more than an operator-set
  fraction of profiled patients are dropped.
- **Transform:** beta → M-value via the logit transform (Du et al. 2010,
  *BMC Bioinformatics* 11:587), with beta clipped to
  `[1e-6, 1 - 1e-6]` first so 0/1 stay finite. M-values are the
  recommended representation for regression/factor-analysis use, since
  beta values are bounded and heteroscedastic near the extremes.
- **Feature selection:** top 5,000 probes by variance across patients,
  kept at probe (not gene-aggregated) resolution — gene-level
  interpretation happens later, during biomarker scoring, by annotating
  top-loading probes back to genes.

### 2.5 Somatic mutation (`oncocartograph.preprocessing.mutation`)

- **Input:** GDC Masked Somatic Mutation MAFs (MC3-derived ensemble
  calls), already one tumor-normal-paired file per patient.
- **Variant filter:** an explicit allow-list of non-synonymous
  `Variant_Classification` values (`Missense_Mutation`,
  `Nonsense_Mutation`, `Nonstop_Mutation`, `Frame_Shift_Del`,
  `Frame_Shift_Ins`, `In_Frame_Del`, `In_Frame_Ins`, `Splice_Site`,
  `Translation_Start_Site`), confirmed against the real vocabulary
  present in downloaded files. `Silent`/`Intron`/UTR/`RNA` variants are
  excluded.
- **Representation:** binary gene × patient matrix (1 = ≥1 non-synonymous
  variant in that gene for that patient), for a Bernoulli MOFA+ view.
- **Recurrence filter:** genes mutated in fewer than 3 patients (~2% of
  the mutation-data subset) are dropped. No restriction to a curated
  cancer gene panel — deliberately whole-exome, to preserve discovery
  potential for non-obvious biomarkers rather than only recovering
  already-known genes.
- **Real sanity check:** on a 10-patient subset, TP53 was mutated in 7/10
  patients, consistent with its well-established ~80% mutation frequency
  in TNBC.

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

**Added during `feat/preprocessing`:** the copy number relative-log2
transform (§2.3) does not correct for tumor purity or ploidy — a real
tumor with whole-genome doubling will show systematically shifted
relative copy number values that this pipeline cannot distinguish from
true focal amplification/deletion without a purity/ploidy estimate per
sample (see ADR 0005). Methylation coverage in the TNBC cohort is
substantially incomplete (~106/143 patients have a usable 450K profile,
§2.1) — MOFA+'s partial-view handling means this doesn't block the
pipeline, but it does mean the methylation view's contribution to any
downstream factor is supported by a smaller effective N than the other
three omics.
