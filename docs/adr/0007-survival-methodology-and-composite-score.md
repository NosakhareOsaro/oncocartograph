# ADR 0007: Survival methodology (Cox PH) and composite biomarker score design

## Status

Accepted (2026-07-20)

## Context

The scoring package (`oncocartograph.scoring`) is this project's core novel
contribution and is designed to be extractable to a standalone PyPI
package. Two things had to be settled before writing any code: which
survival model to use, and how to combine survival evidence with
druggability evidence and each candidate's selection pathway (MOFA+
factor loading, for RNA-seq/methylation/CNV; mutation recurrence, for
mutations) into one ranked score.

### Survival model: Cox PH is a data constraint, not a preference

The original project plan asked whether Cox proportional hazards or
Fine-Gray competing-risks regression should be used. Inspecting the real
downloaded TCGA-BRCA clinical file settled this: only overall survival
fields (`vital_status`, `death_days_to`, `last_contact_days_to`) are
populated for the 143-patient TNBC cohort (16 events, 11.2%);
`days_to_patient_progression_free` and `days_to_tumor_progression` are
0/143 populated. There is no cause-of-death or recurrence coding to
define a competing event against, so Fine-Gray is not available as an
option, not merely unnecessary. Cox PH on overall survival is used.

**16 events is a real statistical power constraint**, not a footnote: per
the standard ~10-events-per-covariate rule of thumb, this is right at the
floor for stable *univariate* estimates and precludes any multivariate
adjustment. This was flagged before any code was written, and the actual
screening run confirmed the expected consequence: 0 of 709 screened
candidates survived Benjamini-Hochberg FDR correction.

### A more severe consequence than anticipated: degenerate mutation fits

Running the real screen also surfaced something more serious than
"underpowered": **712 of 845 (84%) recurrence-filtered mutation genes
could not produce a usable Cox estimate at all** -- not just a
non-significant one. lifelines does not always raise an exception for
this; it can return a "successful" fit with either `NaN` summary
statistics, or a finite-looking but meaningless one (hazard ratio near
zero with an infinite upper confidence bound), when a rare binary
covariate's mutated subgroup contains zero or too few of the 16 total
events to support a stable partial-likelihood estimate. An initial fix
that checked only for `NaN` caught 92/845 (~11%) but missed the
finite-but-infinite-CI pattern entirely; broadening the check to
`np.isfinite` across all four summary statistics (hazard ratio and both
CI bounds, plus p-value) caught the true 84%.

## Decision

**Survival:** univariate Cox PH (`lifelines`) on overall survival,
applied identically whether the candidate's value is continuous
(RNA-seq/methylation/CNV) or binary (mutation), via
`oncocartograph.scoring.survival.fit_univariate_cox`. A fit is discarded
(returns `None`, excluded from downstream screening and FDR correction)
if any of hazard ratio, both CI bounds, or the p-value is non-finite.

**Composite score:** a weighted average over three evidence axes --
survival (0.5), druggability (0.35), selection-pathway (0.15) -- computed
per-candidate via `oncocartograph.scoring.composite.composite_biomarker_score`,
with weights renormalized over whichever axes are actually present for a
given candidate (see ADR discussion in `evidence.py`/`composite.py`
docstrings for why this matters specifically for mutation candidates,
which structurally have no `IntegrationEvidence`). Survival scoring
floors protective associations (HR<=1) to 0, reflecting the project's
choice to treat "druggable biomarker" as "something harmful to disrupt."

## Alternatives considered

**Treating the mutation view's degenerate fits as merely non-significant
rather than excluding them.** Rejected -- an infinite confidence interval
or NaN standard error means the model did not converge to an informative
estimate at all; including such a "result" in FDR correction is not just
uninformative, it is mathematically invalid input (scipy's
`false_discovery_control` raises on out-of-range p-values, and even where
it wouldn't, an unbounded HR estimate should never be presented as
evidence).

**Dichotomizing continuous biomarkers by median split** and using a
log-rank test uniformly instead of continuous Cox regression. Rejected --
loses information and introduces an arbitrary cutpoint, which matters
more, not less, given the small number of events available.

## Consequences

- The real screening run's honest result (0 candidates significant after
  FDR correction, 84% of mutation candidates un-fittable) is documented
  in `docs/methods.md` §4 rather than smoothed over. This project's
  biomarker rankings should be read as hypothesis-generating, not
  confirmatory, given TNBC's reduced sample size within TCGA-BRCA.
- Any future re-run with a larger event count (e.g. if `feat/validation`'s
  external cohort has more events) should re-examine whether the 84%
  mutation attrition rate persists -- it may be substantially an artifact
  of this specific cohort's event scarcity rather than a general property
  of mutation-survival testing.
