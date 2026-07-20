# ADR 0009: External replication methodology, pre-registered criterion, and result

## Status

Accepted (2026-07-20)

## Context

`feat/scoring-package` (ADR 0007) produced 0/709 FDR-significant TCGA-BRCA
TNBC biomarker candidates -- an honest, underpowered-but-real screening
result, not a bug. `feat/validation` needed a way to test these candidates
against an independent cohort (GSE96058/SCAN-B, ADR 0003) that did not
implicitly assume TCGA's associations were real signal worth replicating
at nominal significance, since there is no significant TCGA result to
replicate in that sense.

A second, narrower question was how much of Burstein et al. (2015)'s TNBC
subtyping work to reproduce: the full LAR/MES/BLIS/BLIA transcriptomic
classifier is out of scope for this work package's time budget.

Both the statistical criterion and the Burstein scope reduction were
fixed *before* the real analysis was run, specifically so neither could be
quietly redefined after seeing the result.

## Decision

**Primary, pre-registered criterion: direction concordance.** For each
TCGA candidate with usable GSE96058 evidence, compare the sign of
`log(hazard_ratio)` between TCGA and GSE96058. Test the observed
concordance rate against the 50% chance rate with a one-sided exact
binomial test (`scipy.stats.binomtest`, `alternative="greater"`) at
`alpha=0.05`. `success = p_value < 0.05`, nothing else.

**Secondary, informational-only: nominal p-value replication.** Reported
alongside, explicitly not part of the pass/fail bar -- with 0/709
FDR-significant TCGA hits, demanding significance-replication would be
statistically incoherent.

**Burstein check, scope-reduced:** rather than reproducing the LAR/MES/
BLIS/BLIA classifier, check whether five genes with well-documented TNBC
biology (AR, PTEN, CD274/PD-L1, PDCD1/PD-1, CTLA4) show the
literature-expected hazard direction in real GSE96058 data. This is a
descriptive plausibility note (`oncocartograph.validation.burstein_check`),
not a statistical test and not a pipeline-wide pass/fail gate.

Implementation: `oncocartograph.validation.replication` and
`oncocartograph.validation.burstein_check`.

## Rationale

Direction concordance is the correct question to ask of a discovery
cohort that produced no FDR-significant hits: it asks "do these
candidates' effect *directions* generalize better than chance," which is
answerable and meaningful even when none of the underlying estimates were
individually significant. Testing for replicated significance instead
would guarantee failure by construction (TCGA already failed its own
significance bar) and would not distinguish "the pipeline is unsound"
from "TCGA-BRCA's TNBC subset is underpowered," which is the actual,
already-documented (ADR 0007) constraint.

## Alternatives considered

**Requiring nominal p-value replication as the primary bar.** Rejected:
incoherent given the source screen's 0/709 FDR-significant result, and
would conflate two different questions (direction vs. significance).

**Full Burstein LAR/MES/BLIS/BLIA subtype reproduction.** Rejected as
out of scope for this work package's time budget; confirmed with the
project owner as an explicit scope reduction, not a silent one.

## Consequences: the real, run result

Running the pre-registered analysis on real data (152 TCGA rna_seq
candidates re-derived from raw STAR-Counts files via the committed
preprocessing pipeline; 152/152 Ensembl IDs resolved to gene symbols via
the live Open Targets API; GSE96058 TNBC cohort N=143, 26 events; 109
candidates with usable GSE96058 Cox evidence):

- **Primary criterion: FAILED.** 45/109 = 41.3% direction concordance --
  *below* the 50% chance rate. One-sided binomial p=0.973, nowhere near
  `alpha=0.05`. `success=False`.
- **Secondary (informational): 0/109 candidates replicated at nominal
  p<0.05 in GSE96058** (consistent with TCGA's own null screen; not a
  separate surprise).
- **Burstein plausibility check: 5/5 markers plausible.** AR, PTEN,
  CD274, PDCD1, and CTLA4 all showed the literature-expected protective
  (HR<1) direction in real GSE96058 data.

This is reported and documented (`docs/methods.md` §7, `README.md`) as a
genuine limitation: the TCGA-BRCA TNBC discovery screen's candidates do
not externally replicate above chance in an independent RNA-seq TNBC
cohort. The Burstein check's 5/5 result is a real, separate, reassuring
data point about the GSE96058 ingestion/scoring machinery itself, but it
does not offset or override the primary null result -- both are reported
as-is, per the explicit instruction not to reframe the pipeline's claims
to fit whatever the analysis returned.

The per-candidate replication table (TCGA vs. GSE96058 hazard ratios and
concordance flags for all 152 candidates) is persisted at
`data/processed/gse96058_replication_table.csv` as durable evidence of
the real analysis.
