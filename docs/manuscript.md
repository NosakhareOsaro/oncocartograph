# OncoCartograph: a reproducible multi-omics pipeline for TNBC biomarker prioritisation, with a pre-registered external validation that failed

**Nosakhare Osaro**

*Preprint / project write-up. Not yet submitted to a journal. Code, data provenance, and full methods: https://github.com/NosakhareOsaro/oncocartograph.*

## Abstract

OncoCartograph is a reproducible, auditable pipeline for TNBC biomarker prioritisation from TCGA-BRCA (1,097 patients → 143 TNBC by explicit, cited ER/PR/HER2 rules), integrating RNA-seq, methylation, copy number, and mutation data via MOFA+ and scoring candidates with a standalone, independently-citable composite package combining survival association (Cox PH) and druggability evidence (Open Targets/ChEMBL); one honest finding from integration is that the mutation view contributes essentially no variance to any MOFA+ factor, so mutation-derived candidates rely on direct recurrence statistics instead. Of 709 screened candidates, none survived FDR correction, reflecting the TNBC sub-cohort's limited event count — rankings are hypothesis-generating, not confirmatory. We pre-registered a falsifiable external-validation criterion (direction concordance vs. chance) before testing against an independent RNA-seq TNBC cohort (GSE96058/SCAN-B, N=143); **this primary criterion failed** (41.3% concordance, below the 50% chance rate, p=0.97), a pattern consistent with — though not proof of — the discovery screen's limited statistical power, while a secondary five-gene known-biology plausibility check (Burstein et al. 2015) was directionally consistent in 5/5 cases. We present this as evidence that the pipeline's pre-registered, falsifiable design works as intended — it detected its own lack of replicable signal rather than obscuring it — and argue the reproducible cohort definition, standalone scoring package, and validation methodology are re-usable contributions independent of whether TCGA-BRCA's TNBC subset yields significant biomarkers on its own.

## 1. Introduction

Triple-negative breast cancer (TNBC) — ER-negative, PR-negative, HER2-negative by IHC/FISH — has the fewest approved targeted therapies and the poorest 5-year survival of any breast cancer subtype. Many public multi-omics integration demonstrations run on TCGA pan-cancer or generic BRCA cohorts without a rigorous, reproducible subtype definition, which limits their clinical relevance and makes their biomarker calls hard to audit or trust: readers cannot tell whether a paper's exact ER/PR/HER2 cutoffs, or its handling of equivocal/missing calls, would reproduce the same cohort from the same source data.

This project set out to build three things, each intended to stand on its own as re-usable evidence of rigor, independent of whether any specific biomarker candidate holds up:

1. A **reproducible, auditable TNBC sub-cohort definition** from TCGA-BRCA, with explicit thresholds cited to the ASCO/CAP conventions (Hammond et al. 2010; Wolff et al. 2013/2018) and a full per-patient audit trail, rather than an undocumented filter.
2. A **standalone, independently-citable composite biomarker scoring package** (`oncocartograph.scoring`), combining survival-association evidence with druggability evidence, with a zero-cross-import architecture enforced by an AST-based test so it can be extracted to its own package without modification.
3. A **pre-registered, falsifiable external validation** against an independent cohort, with the primary success/failure criterion fixed *before* the analysis was run, specifically so it could not be redefined after seeing the result.

The third of these is the paper's central methodological argument: a validation step is only evidence of quality if it was capable of failing, and if a negative result is reported as such. We built the criterion in `feat/validation` before running the real analysis (docs/adr/0009), and it did fail. We report that failure in the same terms as every positive finding in this pipeline, in this abstract, in this section, and in the conclusion — not softened, hedged, or confined to a limitations appendix.

## 2. Methods (condensed)

Full methodology, exact real-data numbers, and every judgment call's rationale are in [`docs/methods.md`](methods.md) and the architecture decision records in [`docs/adr/`](adr/); this section summarizes only what is needed to follow the Results.

**2.1 Cohort.** TCGA-BRCA patients are classified TNBC if and only if ER, PR, and HER2 are all unambiguously negative by IHC (HER2 IHC-2+ resolved by reflex FISH); any missing or indeterminate call is excluded and logged, never imputed (methods.md §1, ADR 0001).

**2.2 Preprocessing.** RNA-seq: GDC STAR-Counts, low-expression filtered, DESeq2 median-of-ratios normalization + VST (`pydeseq2`), top 2,000 variable genes. Methylation: SeSAMe beta values, logit-transformed to M-values, top 5,000 variable probes. Copy number: gene-level absolute integer copy number (not GISTIC2 categorical calls, see §5.1), relative-to-diploid log2 transform, top 2,000 variable genes. Mutation: MC3 MAFs, non-synonymous variant allow-list, binary gene×patient matrix, ≥3-patient recurrence filter (methods.md §2).

**2.3 Integration.** MOFA+ (`mofapy2`/`mofax`, pure Python), K=15 factors, `scale_views=True`, `convergence_mode="slow"`, seed 42, factors screened post-hoc at a ≥2%-variance-explained threshold (methods.md §3).

**2.4 Scoring.** Two candidate-generation pathways feed one composite formula: MOFA+ factor loading (RNA-seq/methylation/copy-number) and mutation-recurrence (bypassing MOFA+ entirely, since integration contributes no usable signal for mutations — see Results 3.2). Both receive identical univariate Cox PH survival-association treatment on overall survival (the only model TCGA-BRCA's clinical file supports), Benjamini-Hochberg FDR correction, and a weighted composite score (survival 0.5, druggability 0.35, selection-pathway 0.15, renormalized over whichever axes apply to a given candidate) (methods.md §4).

**2.5 Druggability.** Open Targets GraphQL tractability buckets (3-tier scoring) combined with ChEMBL maximum clinical trial phase via `max(tier_score, max_phase/4)` (methods.md §6).

**2.6 External validation.** GSE96058 (SCAN-B), an independent RNA-seq TNBC cohort (chosen over the microarray alternative GSE58812 specifically to avoid a cross-platform confound, ADR 0003), N=143 TNBC by histopathology status, 26 events. The **pre-registered primary criterion**: direction concordance (sign of log hazard ratio) between TCGA and GSE96058, tested against the 50% chance rate with a one-sided exact binomial test at α=0.05 — chosen instead of demanding replicated significance, because the discovery screen itself produced zero FDR-significant hits (§2.4; ADR 0009). A secondary, scope-reduced check of five genes with documented TNBC biology (AR, PTEN, CD274, PDCD1, CTLA4) against Burstein et al. (2015)'s subtype findings is reported as a qualitative plausibility note, not a statistical test.

**2.7 Reproducibility.** Every stage is implemented as tested library code (`src/oncocartograph/`) and wired into one real Snakemake DAG (`workflows/Snakefile`) callable end-to-end from raw GDC/GEO/Open Targets/ChEMBL data; every stochastic step logs its seed; every downloaded artifact carries a checksummed provenance sidecar (methods.md §7, `docs/data_sources.md`).

## 3. Results

### 3.1 Cohort definition

Of 1,097 TCGA-BRCA patients, **143 (13.0%) were classified TNBC**. Of the 954 excluded: 737 receptor-positive, 196 with a missing/indeterminate IHC call, 21 HER2-equivocal-unresolved-by-FISH — mutually exclusive categories summing exactly to 954, because indeterminacy takes precedence over a receptor-positive verdict when both apply. 158 RNA-seq, 348 methylation, 426 copy number, and 126 mutation files were ingested for this cohort (file counts exceed 143 because of multiple samples/aliquots per patient; resolved to one Primary Tumor sample per patient per omic layer, yielding 142/104–106/142/122 patients respectively).

### 3.2 Multi-omics integration

Trained on copy number (2,000×142), RNA-seq (2,000×142), methylation (5,000×104), and mutation (845×122) — 143 patients total via the union of views, no forced complete-case restriction. 12 of 15 factors cleared a ≥2%-variance-explained screening threshold. Two findings are load-bearing, not incidental:

- **Factor1 is overwhelmingly copy-number-driven** (56.5% of CNV variance, <1.1% everywhere else), most likely reflecting a broad genomic-instability axis rather than shared cross-omic biology.
- **The mutation view contributes essentially no variance to any factor** (≤0.003% everywhere). This is a real result of TNBC's sparse, largely-private mutation landscape at this cohort size, not a bug — and it directly shaped the scoring architecture (§2.4): mutation-derived candidates use direct recurrence statistics, because MOFA+ factor loadings give them nothing to select on.

### 3.3 Composite scoring

709 candidates screened (576 MOFA+-derived + 133 mutation-derived of 845 attempted) against real survival data (16 events across 143 patients). **0 of 709 candidates survived FDR correction** — an expected consequence of the event count, flagged before the screen was run, not a surprise after the fact. The top-ranked candidate (RNA-seq, `ENSG00000172551.11`, composite score 0.61, HR=1.25, p_adj=0.21) was MOFA+-derived; a mutation-derived candidate (`GTF3C1`, HR≈11.9, p_adj=0.43) ranked 6th, demonstrating both selection pathways contribute to one ranked list, not two disconnected pipelines — though `GTF3C1`'s extreme hazard ratio reflects a very small number of carriers/events and should be read with corresponding caution. Screening the real recurrence-filtered mutation data also surfaced a genuine bug: 84% of genes (712/845) produced a degenerate Cox fit (too few events in the carrier subgroup) that an initial NaN-only check caught only 11% of; broadening to check all four fit statistics for finiteness caught the true rate.

### 3.4 Druggability re-scoring

Real Open Targets + ChEMBL evidence was obtained for 480/480 (100%) of the RNA-seq/copy-number/mutation candidates (methylation's 229 candidates excluded — CpG probe IDs are not gene identifiers without a 450K manifest, out of scope here). Adding druggability substantially reshuffled the ranking (Spearman ρ=0.656 vs. the druggability-absent ranking; only 4/20 previous top-20 candidates remained in the top 20). The mutation-derived `TP53` candidate rose to 7th overall on the strength of real, well-documented drug-development evidence (`tractability_score=0.75`, `chembl_max_phase=3.0`) that MOFA+ factor loadings could never have surfaced on their own (§3.2). This is evidence the composite score does real integrative work, not merely re-deriving the survival ranking with cosmetic reweighting — though, given §3.3's null significance result, it should still be read as changing which statistically-modest associations are worth following up on, not as adding statistical confidence to any of them.

### 3.5 External validation: the primary criterion failed

GSE96058's TNBC sub-cohort: N=143, 26 events. All 152 fittable TCGA RNA-seq candidates were carried forward (not a top-N subset); 152/152 resolved to a gene symbol via the live Open Targets API; 109 had usable GSE96058 Cox evidence.

**Primary, pre-registered criterion: failed.** 45/109 (41.3%) candidates were direction-concordant between TCGA and GSE96058 — *below* the 50% rate expected by chance. One-sided exact binomial test: p=0.973, nowhere near α=0.05. This is not a marginal miss.

**Secondary, informational: 0/109 replicated at nominal p<0.05** in GSE96058 — consistent with, not a separate surprise beyond, TCGA's own null screen (§3.3).

**Burstein et al. (2015) plausibility check: 5/5.** AR, PTEN, CD274, PDCD1, and CTLA4 all showed the literature-expected protective (HR<1) direction in real GSE96058 data. This is a real, separate, reassuring signal about the validation pipeline's mechanics (correct cohort definition, correct expression parsing, correct Cox fitting) — it does not offset or soften the primary null result, and is reported alongside it rather than in its place.

**Independent reproduction.** Re-running the entire pipeline a second time through the real Snakemake DAG built for this project's final work package — a freshly re-derived RNA-seq view from raw counts and a freshly-retrained MOFA+ model, not a cached artifact — reproduced the same qualitative conclusion with different exact numbers, as expected given MOFA+'s stochastic training: 39.8% concordance (113 fittable candidates, 45 concordant), p=0.988, Burstein check still 5/5. The primary criterion's failure is a property of the underlying data and design, not an artifact of one specific run.

## 4. Discussion

The central methodological question this project asked was not "did TCGA-BRCA's TNBC subset yield validated biomarkers" — with 16 survival events across 143 patients, that was always a long shot, flagged as a real statistical-power constraint before any screening was run (§3.3). The question was whether a validation step designed to be capable of failing would, in fact, report a failure honestly when one occurred. It did.

A pattern consistent with — though not proof of — the explanation is straightforward: a discovery screen that produced zero FDR-significant hits selected candidates whose survival associations were, individually, indistinguishable from noise; there was no principled reason to expect noise-level associations to point in a consistent direction in an independent cohort, and on real data they did not (41.3%, then 39.8% on independent reproduction — both below chance, not scattered around it). We cannot rule out other contributing factors (residual technical heterogeneity between the two RNA-seq processing pipelines, TCGA-BRCA's TNBC subset being non-representative of SCAN-B's) without a second, adequately-powered discovery cohort to test them against — a limitation we state plainly rather than resolve by assertion.

The Burstein et al. (2015) check's 5/5 result is worth discussing on its own terms, separate from the primary finding: it demonstrates the GSE96058 ingestion, cohort classification, and Cox-fitting machinery behaves correctly on genes with genuine, independently-documented TNBC biology, even though the TCGA-derived candidate set built by the same machinery did not replicate. This distinguishes a "the validation pipeline is broken" explanation (which the Burstein result argues against) from "TCGA-BRCA's TNBC subset did not contain replicable signal at this sample size" (which the primary result is most consistent with).

## 5. Limitations

**5.1 Cohort and data.** TCGA-BRCA IHC calls were made across many contributing institutions without fully centralized scoring; some inter-site variability in ER/PR/HER2 calls is expected and cannot be corrected for retrospectively. GSE96058 is Swedish and RNA-seq-only (no methylation/CNV/mutation validation), so external validation here is necessarily expression-centric. The Burstein et al. (2015) benchmark cohort is not population-matched to TCGA-BRCA.

**5.2 Copy number.** The relative-to-diploid log2 transform does not correct for tumor purity or ploidy; a tumor with whole-genome doubling will show systematically shifted relative copy number this pipeline cannot distinguish from true focal amplification/deletion.

**5.3 Methylation coverage.** Only ~104–106/143 TNBC patients have a usable 450K profile; MOFA+'s partial-view handling means this does not block the pipeline, but the methylation view's contribution to any factor is supported by a smaller effective N than the other omics.

**5.4 MOFA+ training.** Training reached the 1,000-iteration cap without formally satisfying `"slow"` mode's ELBO tolerance (though the per-iteration change was <0.0002% at the end) — a longer iteration cap is a candidate refinement, not applied here.

**5.5 Statistical power.** With only 16 observed events in the 143-patient TCGA discovery cohort, 0/709 candidates survived FDR correction, and 84% of recurrence-filtered mutation genes could not produce a finite Cox estimate at all. Every composite-scored ranking in this project should be read as hypothesis-generating, not confirmatory.

**5.6 External validation.** The pre-registered primary criterion failed (§3.5); this is the paper's headline limitation, not a footnote. Composite-scored candidates from this pipeline should not be treated as externally validated biomarkers pending a better-powered discovery cohort.

**5.7 Reproducibility scope.** The Snakemake pipeline (§2.7) reproduces every result from raw source data end-to-end, but a from-scratch run involves multi-gigabyte downloads and live third-party API calls (GDC, GEO, Open Targets, ChEMBL) whose availability and exact response content this project does not control going forward.

## 6. Conclusion

This project's contribution is not a validated TNBC biomarker panel — the pre-registered external validation failed on its primary criterion (41.3%, then 39.8% on independent reproduction, both below the 50% chance rate), and we report that plainly rather than around it. The contribution is a reproducible, auditable cohort definition; a standalone, independently-citable composite scoring package with an enforced zero-cross-import architecture; and, most centrally, a pre-registered, falsifiable validation methodology that was capable of failing and did — surfacing that TCGA-BRCA's TNBC subset does not, on its own, support externally-replicable biomarker discovery at this event count, rather than obscuring that fact behind a positive-sounding but unfalsifiable design. We argue these are re-usable contributions independent of whether TCGA-BRCA's TNBC subset yields significant biomarkers on its own, and that reporting a negative pre-registered result this plainly is itself evidence of the methodology's rigor, not a weakness to be minimized.

## References

- Hammond MEH, Hayes DF, Dowsett M, et al. American Society of Clinical Oncology/College of American Pathologists guideline recommendations for immunohistochemical testing of estrogen and progesterone receptors in breast cancer. *Arch Pathol Lab Med*. 2010;134(7):907-922.
- Wolff AC, Hammond MEH, Hicks DG, et al. Recommendations for human epidermal growth factor receptor 2 testing in breast cancer: American Society of Clinical Oncology/College of American Pathologists clinical practice guideline update. *J Clin Oncol*. 2013;31(31):3997-4013 (updated 2018, Allison KH et al., *J Clin Oncol*. 2020;38(12):1346-1366).
- Liu J, Lichtenberg T, Hoadley KA, et al. An integrated TCGA pan-cancer clinical data resource to drive high-quality survival outcome analytics. *Cell*. 2018;173(2):400-416.e11.
- Argelaguet R, Velten B, Arnol D, et al. Multi-Omics Factor Analysis—a framework for unsupervised integration of multi-omics data sets. *Mol Syst Biol*. 2018;14(6):e8124.
- Argelaguet R, Arnol D, Bredikhin D, Deloro Y, Velten B, Marioni JC, Stegle O. MOFA+: a statistical framework for comprehensive integration of multi-modal single-cell data. *Genome Biol*. 2020;21:111.
- Du P, Zhang X, Huang CC, et al. Comparison of Beta-value and M-value methods for quantifying methylation levels by microarray analysis. *BMC Bioinformatics*. 2010;11:587.
- Saal LH, Vallon-Christersson J, Häkkinen J, et al. The Sweden Cancerome Analysis Network - Breast (SCAN-B) Initiative: a large-scale multicenter infrastructure towards implementation of breast cancer genomic analyses in the clinical routine. *Genome Med*. 2015;7:20.
- Brueffer C, Vallon-Christersson J, Grabau D, et al. Clinical value of RNA sequencing-based classifiers for prediction of the five conventional breast cancer biomarkers: a report from the population-based multicenter Sweden Cancerome Analysis Network-Breast Initiative. *JCO Precis Oncol*. 2018;2:1-18.
- Burstein MD, Tsimelzon A, Poage GM, et al. Comprehensive genomic analysis identifies novel subtypes and targets of triple-negative breast cancer. *Clin Cancer Res*. 2015;21(7):1688-1698.
