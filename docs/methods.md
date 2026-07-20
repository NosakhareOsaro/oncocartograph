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

### 3.1 Implementation

MOFA+ models are trained with `mofapy2` and read/interpreted with `mofax`
— both pure Python (see
[`docs/adr/0006-mofa-plus-implementation-and-training.md`](adr/0006-mofa-plus-implementation-and-training.md)
for why this needed no R environment, correcting an earlier assumption).

### 3.2 View construction and likelihoods

Each of the four preprocessed omic matrices (§2) is melted into MOFA+'s
long-format input (`sample, feature, view, group, value`) and combined
into a single training input. Missing values — whether a patient
entirely absent from a view, or a scattered missing value within an
otherwise-present sample — are represented by the row's absence, not by
an explicit `NaN`; both cases were verified directly to train correctly
before this was relied on.

| View | Likelihood | Rationale |
|---|---|---|
| RNA-seq (VST, top-2,000 genes) | Gaussian | Continuous, approximately normal after VST |
| Methylation (M-value, top-5,000 probes) | Gaussian | Continuous, unbounded (unlike beta values) |
| Copy number (relative log2, top-2,000 genes) | Gaussian | Continuous, log-ratio scale |
| Mutation (binary, recurrence-filtered) | Bernoulli | Discrete presence/absence |

`scale_views=True` scales each view to unit variance before training, so
no view dominates purely through differences in value magnitude.

### 3.3 Training configuration

- **Factors:** initialised with K=15, a middle-of-the-road choice for a
  ~100-150 sample cohort matching common practice in published MOFA+
  multi-omics analyses at this scale. ARD priors on feature weights
  (mofapy2 default) provide within-factor sparsity; factors explaining
  <2% of variance in every view are excluded from downstream
  interpretation (not training) — this is a post-hoc screening step, not
  a change to K itself.
- **Convergence:** `convergence_mode="slow"`, the tightest ELBO tolerance
  mofapy2 offers, chosen for a final citable analysis over faster/looser
  exploratory settings.
- **Seed:** `Settings.random_seed` (project default 42), logged by
  mofapy2 at training time.

### 3.4 Real training result (2026-07-20)

Trained on the full real preprocessed cohort — copy number (2,000 x 142),
RNA-seq (2,000 x 142), methylation (5,000 x 104), mutation (845 x 122) —
with no forced complete-case restriction; 143 total patients are
represented in the resulting factor matrix (the union across views).
Training took 77.7s and ran the full 1,000-iteration cap without
formally reaching `"slow"` mode's ELBO tolerance, though by the final
iterations ELBO was changing by <0.0002% per step — practically
converged for interpretation purposes. A longer `max_iterations` cap is
a candidate refinement for a final canonical run, not applied here.

**Variance explained per factor** (%, screening threshold ≥2% in at
least one view):

| Factor | Copy number | Methylation | Mutation | RNA-seq |
|---|---|---|---|---|
| Factor1 | 56.50 | 0.01 | 0.00 | 1.09 |
| Factor2 | 0.31 | 17.67 | 0.00 | 2.77 |
| Factor3 | 0.20 | 5.77 | 0.00 | 7.93 |
| Factor4 | 1.00 | 4.78 | 0.00 | 5.14 |
| Factor6 | 4.91 | 0.05 | 0.00 | 0.09 |
| Factor5 | 0.11 | 4.55 | 0.00 | 2.32 |
| Factor7 | 0.48 | 4.40 | 0.00 | 0.11 |
| Factor8 | 3.26 | 0.32 | 0.00 | 0.24 |
| Factor10 | 2.95 | 0.04 | 0.00 | 0.23 |
| Factor9 | -1.20 | 2.82 | 0.00 | 1.87 |
| Factor13 | 2.42 | 0.00 | 0.00 | 0.01 |
| Factor14 | 2.07 | 0.02 | 0.00 | 0.01 |
| ~~Factor15~~ | 1.81 | 0.01 | 0.00 | 0.16 |
| ~~Factor12~~ | 1.79 | 0.31 | 0.00 | 0.57 |
| ~~Factor11~~ | 1.61 | 0.75 | 0.00 | 0.47 |

**12 of 15 factors clear the ≥2% screening threshold** in at least one
view (Factors 11, 12, 15 struck through above do not and are excluded
from downstream interpretation). Two things worth noting honestly:

- **Factor1 is overwhelmingly copy-number-driven** (56.5% of CNV
  variance, <1.1% everywhere else) — likely reflects a broad genomic
  instability / aneuploidy axis rather than a multi-omic signal, and
  should be interpreted as such rather than assumed to reflect shared
  cross-omic biology.
- **The mutation view contributes essentially nothing to any factor**
  (≤0.003% variance explained everywhere). This is a real result, not a
  bug: after the ≥3-patient recurrence filter, 845 genes remain but
  mutation events are inherently sparse and largely private to individual
  patients, giving MOFA+ little shared structure to extract via a
  Bernoulli likelihood at this cohort size. Mutation-derived biomarker
  candidates in this project will need to rely more on direct
  recurrence/association statistics (`feat/scoring-package`) than on
  MOFA+ factor loadings.
- Most substantial shared cross-omic structure appears in
  Factors 2-5 and 9 (methylation + RNA-seq, and to a lesser extent
  copy number), consistent with expression and methylation both reading
  out overlapping regulatory biology.

## 4. Composite biomarker scoring

Implemented in `oncocartograph.scoring`, a standalone package with zero
dependency on the rest of `oncocartograph` (mechanically enforced by
`tests/scoring/test_decoupling.py`), designed to be extractable to a
PyPI package independently.

### 4.1 Two candidate-generation pathways, one composite score

Per the v0.2.0 MOFA+ finding that the mutation view contributes
essentially no variance to any factor (§3.4), candidates enter the
scoring pipeline via one of two distinct, explicitly typed pathways:

- **MOFA+-derived** (RNA-seq/methylation/CNV): selected by factor
  loading on one of the 12 factors passing the ≥2% variance-explained
  screen (`mofax.get_top_features`). Recorded as `IntegrationEvidence`
  (factor, loading weight, the factor's view-variance-explained).
- **Mutation-derived**: selected by the recurrence filter already applied
  in preprocessing (≥3 patients), bypassing MOFA+ entirely. Recorded as
  `RecurrenceEvidence` (mutation count/fraction, plus an optional
  Fisher's-exact categorical association test against any supplied
  categorical outcome -- e.g. `vital_status` now, or a subtype label once
  `feat/validation` defines one).

Both pathways receive **identical survival-association treatment** (§4.2)
and feed the **same composite formula** (§4.3) -- the pathway only
determines which "selection pathway" evidence is attached, not a
different scoring formula.

### 4.2 Survival association: Cox PH on overall survival

Univariate Cox proportional hazards (`lifelines`), continuous covariate
for RNA-seq/methylation/CNV, binary covariate for mutations -- the same
statistical treatment regardless of candidate origin. Cox PH on overall
survival is used because it is the only model the real data supports:
the TCGA-BRCA clinical file has no cause-of-death or recurrence coding to
define a Fine-Gray competing event against, and no PFS/recurrence fields
at all are populated for this cohort (see
[`docs/adr/0007-survival-methodology-and-composite-score.md`](adr/0007-survival-methodology-and-composite-score.md)).

**A fit is discarded (excluded, not included with invalid statistics) if
any of the hazard ratio, either 95% CI bound, or the p-value is
non-finite.** This matters more than it might sound: on the real
recurrence-filtered mutation data (845 genes, 122 patients with mutation
data, only 16 total events across the cohort), **712/845 genes (84%)
produced a degenerate fit** -- lifelines does not always raise on this;
it can return `NaN` statistics or a finite-looking but meaningless
estimate (hazard ratio near zero, infinite upper CI bound) when a rare
mutation's carrier subgroup contains too few of the 16 events to support
a stable estimate. An initial NaN-only check caught just 92/845 (~11%);
broadening to `np.isfinite` across all four statistics caught the true
84%. This was found and fixed by running the real screen, not anticipated
in advance.

Multiple testing: Benjamini-Hochberg FDR correction
(`scipy.stats.false_discovery_control`), applied across every candidate
screened together in one call.

### 4.3 Composite score

`composite_biomarker_score` combines up to three evidence axes as a
weighted average, **renormalized over whichever axes are actually
present** for a given candidate:

| Axis | Default weight | Source |
|---|---|---|
| Survival | 0.5 | §4.2, every candidate |
| Druggability | 0.35 | Schema-defined here; populated by `feat/drug-target-scoring` |
| Selection pathway | 0.15 | MOFA+ view-variance-explained, or mutation prevalence fraction |

Renormalization matters specifically for mutation candidates: they
structurally have no `IntegrationEvidence` (MOFA+ does not apply to
them), and without renormalization they would be unfairly penalized for
lacking an axis that could never apply. Survival scoring uses `1 - p`
(preferring the FDR-adjusted p-value when available) and **floors
protective associations (hazard ratio ≤ 1) to 0** rather than a negative
value -- this project's explicit choice to treat "druggable biomarker" as
"something harmful to disrupt," not "anything survival-associated
regardless of direction."

### 4.4 Real scoring run (2026-07-20)

Screened 709 total candidates (576 MOFA+-derived across RNA-seq/
methylation/copy number + 133 mutation-derived, of the 845 mutation genes
attempted) against the real 143-patient survival table (16 events), with
no druggability evidence yet (§4.3's druggability axis renormalizes away
until `feat/drug-target-scoring` populates it).

**0 of 709 candidates survived FDR correction (p_adj < 0.05)** -- the
expected honest outcome of screening hundreds of candidates against only
16 events, flagged as a real limitation before this work package began,
not a surprise after the fact. 349/709 candidates showed a harmful
direction (HR>1).

The top-ranked candidate overall was an RNA-seq gene
(`ENSG00000172551.11`, composite score 0.61, HR=1.25, p_adj=0.21) selected
via MOFA+ factor loading. A mutation-derived candidate (`GTF3C1`, HR≈11.9,
p_adj=0.43) ranked 6th, demonstrating the two-pathway architecture works
as designed -- both selection pathways produced real evidence that ranks
together in one composite ordering, not two separate pipelines. `GTF3C1`'s
extreme hazard ratio should be read cautiously: it reflects a very small
number of mutation carriers and events, exactly the kind of estimate
instability §4.2 describes, not a confidently large effect.

**Given no candidate reached significance after correction, this
project's rankings should be treated as hypothesis-generating, not
confirmatory** -- consistent with what a 16-event, 143-patient TNBC
sub-cohort of TCGA-BRCA can realistically support.

## 5. External validation

_Pending — will document the GSE96058 (SCAN-B) cohort's TNBC subsetting
criteria (matched to Section 1.2 where the available fields allow), and the
comparison methodology against Burstein et al. 2015 (PMID 25208879)
subtype-specific druggable target calls, once `feat/validation` lands._

## 6. Drug-target evidence and prioritisation

Implemented in `oncocartograph.drug_targets`, which imports from (but is
not imported by) `oncocartograph.scoring` — the scoring package's
zero-cross-import guarantee (§4) is one-directional.

### 6.1 Identifier resolution

Three candidate identifier types exist in this project's real candidate
set, requiring different resolution:

- **RNA-seq/copy number** (versioned Ensembl IDs, e.g.
  `ENSG00000172551.11`): version suffix stripped, used directly.
- **Mutation** (gene symbols, e.g. `GTF3C1`): batch-resolved to Ensembl
  IDs via Open Targets' `mapIds` query, confirmed to return an empty hit
  list (not an error) for an unresolvable symbol.
- **Methylation** (CpG probe IDs, e.g. `cg00000029`): **out of scope for
  this work package.** Probe IDs are not gene identifiers; mapping to a
  nearest/associated gene needs the Illumina 450K manifest, a new
  reference dataset not part of any prior work package (see
  [`docs/adr/0008-druggability-evidence-sources.md`](adr/0008-druggability-evidence-sources.md)).
  Methylation candidates retain `druggability=None` (renormalized away)
  until this is addressed.

### 6.2 Open Targets tractability

Open Targets' GraphQL `targets(ensemblIds: ...)` query returns ~28
boolean tractability buckets per target across four modalities (SM=Small
Molecule, AB=Antibody, PR=PROTAC, OC=Other Clinical) — confirmed via live
queries against real targets (TP53, GTF3C1) before any code was written.
Collapsed to a single [0, 1] score via a 3-tier scheme (identical
bucket-label text across modalities means this only needs to inspect the
label, not a per-modality table):

| Condition | Score |
|---|---|
| Approved Drug true, any modality | 1.0 |
| Advanced Clinical or Phase 1 Clinical true, any modality | 0.66 |
| Any other bucket true | 0.33 |
| No bucket true | 0.0 |

### 6.3 ChEMBL max clinical phase

ChEMBL's free-text target search is unreliable for exact gene matching
(confirmed: searching "TP53" returns "TP53-binding protein 1" as the top
hit). Targets are instead resolved by exact UniProt accession match
(`target_components__accession__in`, restricted to
`target_type=SINGLE PROTEIN`), using the canonical `uniprot_swissprot`
accession from Open Targets' `proteinIds` (filtering out non-canonical
TrEMBL entries, which are returned alongside it). The maximum `max_phase`
(0-4) across all of ChEMBL's `mechanism` records for that target is used
— real data has multiple mechanism records per target with different
phases (confirmed: TP53 shows both phase 2 and phase 3 records), so a
max, not first-or-last, aggregation is required.

### 6.4 Combining into `tractability_score`

`max(tractability_tier_score, chembl_max_phase / 4)` — either signal
alone can indicate a real drug exists or is in development, so this
project takes the stronger rather than averaging them down.
`chembl_max_phase` is retained on `DruggabilityEvidence` separately for
transparent reporting even though the composite formula (§4.3) reads
only the combined `tractability_score`.

### 6.5 Re-scoring result with all three evidence axes (2026-07-20)

Re-ran composite scoring on the 480 RNA-seq/copy-number/mutation
candidates from §4.4 (methylation's 229 candidates excluded per §6.1's
scope decision) with real Open Targets + ChEMBL druggability evidence
populated for **480/480 (100%)** — a large improvement over §4.4's run,
where druggability was absent and renormalized away for every single
candidate.

**The ranking changed substantially, not marginally:** Spearman
correlation between the before- and after-druggability rankings is
ρ=0.656 (480 candidates, p=1.8×10⁻⁶⁰) — positive and significant, but far
from 1.0. Only **4 of the previous top 20 candidates remained in the top
20** once druggability was added (16/20 displaced); at the top 50, 22/50
(44%) were displaced. Mean absolute rank change across all 480 candidates
was ~78 positions (out of 480, i.e. ~16% of the full list).

**Two clear, explainable patterns drove the reshuffling:**

1. **Candidates with weak survival evidence but strong real druggability
   jumped sharply.** The single biggest mover, `rna_seq:ENSG00000050555.19`,
   went from rank 469.5→53 (Δ=+416.5) — its survival signal alone was
   negligible (composite score 0.0002 before), but real druggability
   evidence (an existing approved/late-phase drug) pulled it into the top
   quarter once that 0.35-weight axis was no longer renormalized away.
   Several other candidates showed the same pattern (Δ>390 positions).
2. **Candidates with moderate survival/selection-pathway evidence but no
   real druggability dropped sharply.** A cluster of copy-number
   candidates sharing one genomic segment (the same shared-loading
   pattern noted in §3.4) dropped from rank 64→351.5 (Δ=−287.5) once
   their actual (weak) druggability evidence was factored in at full
   weight rather than excluded.

**A real, biologically sensible result:** the mutation-derived `TP53`
candidate — bypassing MOFA+ entirely per this project's two-pathway
architecture (§4.1) — ranked **7th overall** (composite score 0.466),
with `tractability_score=0.75` and `chembl_max_phase=3.0`, both
reflecting TP53's genuine, well-documented drug development history
(confirmed directly against the real ChEMBL API in §6.3). This is exactly
the kind of candidate the mutation-recurrence pathway was designed to
surface, since MOFA+ factor loadings could never have found it (§3.4:
mutation view contributed ≤0.003% variance to every factor).

**Interpretation:** this reshuffling is evidence the composite score is
doing real work, not just re-deriving the survival ranking with cosmetic
reweighting — a substantial fraction of top candidates are only
"top" once actionability (does a real drug exist or is one in
development) is considered alongside statistical association. Given
§4.4's honest finding that 0 candidates reached FDR-corrected
significance, this ranking should still be read as hypothesis-generating:
druggability evidence changes *which* statistically-modest associations
are worth following up on, not the underlying statistical confidence in
any of them.

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

**Added during `feat/mofa-integration`:** the real trained model's
mutation view contributes essentially no variance to any of the 15
factors (≤0.003% everywhere, §3.4) — mutation-derived biomarker
candidates will need direct recurrence/association statistics rather
than MOFA+ factor loadings. Factor1 (56.5% of copy number variance,
<1.1% elsewhere) most likely reflects broad genomic instability rather
than a genuine shared multi-omic axis, and should not be over-interpreted
as such during biomarker scoring. Training reached the 1,000-iteration
cap without formally satisfying `"slow"` mode's ELBO tolerance (though
the remaining change per iteration was <0.0002%); a longer iteration cap
is a candidate refinement for a final canonical run.

**Added during `feat/scoring-package`:** with only 16 observed events in
the 143-patient cohort, 0 of 709 screened candidates survived FDR
correction in the real scoring run (§4.4) -- an expected consequence of
the event count, not a pipeline defect, but it means current rankings are
hypothesis-generating rather than confirmatory. More strikingly, 84% of
recurrence-filtered mutation genes (712/845) could not produce a finite
Cox estimate at all, because their mutated subgroup contained too few of
the 16 total events to support a stable fit -- this is a statistical
ceiling imposed by TNBC's reduced representation within TCGA-BRCA, not
something a different modeling choice within this project could fix.
Any future work drawing conclusions from the mutation-recurrence pathway
should treat surviving candidates' effect sizes (hazard ratios) with real
caution, given how few carriers/events typically support them.
