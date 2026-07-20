# ADR 0001: TNBC sub-cohort definition thresholds

## Status

Accepted (2026-07-20)

## Context

TCGA-BRCA is a pan-breast-cancer cohort; TNBC is a clinically defined
subset (ER-negative, PR-negative, HER2-negative) that must be derived from
receptor status fields recorded in the clinical supplement. Many published
analyses that claim to work on "TNBC" do not state their exact IHC
positivity cutoffs or how they handled IHC-equivocal HER2 calls, which
makes their cohort definitions — and therefore their downstream biomarker
claims — difficult to audit or reproduce.

Two decisions had to be made explicitly:

1. What positivity/negativity cutoff to use for ER and PR IHC scoring.
2. How to classify HER2 IHC-equivocal (2+) samples that have no recorded
   reflex FISH result.

## Decision

- **ER/PR:** negative defined as <1% nuclear staining, per the ASCO/CAP
  guideline (Hammond et al. 2010, *Arch Pathol Lab Med* 134:907-922;
  reaffirmed 2020, Allison et al., *J Clin Oncol* 38:1346-1366). This is
  the prevailing clinical convention and the one TCGA's own IHC scoring
  protocol was designed against.
- **HER2:** negative defined as IHC 0/1+, or IHC 2+ (equivocal) **with**
  FISH HER2/CEP17 ratio <2.0, per the ASCO/CAP HER2 testing guideline
  update (Wolff et al. 2013/2018).
- **IHC-equivocal HER2 with no FISH follow-up recorded: excluded and
  logged**, not imputed as negative.

## Alternatives considered

**Treat HER2-equivocal-without-FISH as negative.** This would increase the
final cohort N, which is attractive given TNBC is already a minority
subset of TCGA-BRCA. Rejected because it is a real methodological
assumption — IHC 2+ is, by definition, indeterminate without FISH — and
silently folding indeterminate calls into "negative" would make the
resulting cohort's HER2-negativity claim weaker than it appears, undermining
the project's core goal of an auditable, defensible cohort definition. If
the empirical loss of samples to this exclusion turns out to be large
enough to threaten downstream statistical power, that will be reported
transparently in `docs/methods.md` §8 (Limitations) rather than resolved by
loosening this rule after the fact.

## Consequences

- The ingestion script must emit a full per-patient audit table (raw field
  values in, include/exclude decision out), not just a filtered sample
  list — this is a deliverable, not incidental logging.
- Final cohort N is not fixed in advance; it is an empirical output of
  these rules applied to current GDC data.
