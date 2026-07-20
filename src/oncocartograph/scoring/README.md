# oncocartograph.scoring

Standalone composite biomarker scoring package for prioritising druggable
biomarkers from multi-omics survival-association evidence.

> **Status:** implemented (`feat/scoring-package`, 2026-07-20). Zero
> dependency on the rest of `oncocartograph` -- mechanically enforced by
> `tests/scoring/test_decoupling.py`, not just documented convention.
> Extracting this directory to a standalone PyPI package should require no
> more than an import-path change.

## Purpose

Takes per-biomarker statistical evidence (survival association,
selection-pathway provenance) and druggability evidence, and produces a
single, documented composite priority score per candidate biomarker.

Two distinct, explicitly typed evidence pathways feed the same composite
score -- deliberately not a single MOFA+-only pipeline. In this project's
own MOFA+ run, the mutation view contributed essentially no variance to
any factor (see `docs/methods.md` §3.4), so mutation-derived candidates
need a genuinely different selection mechanism (recurrence filtering) even
though they're scored by the same formula as everything else:

- **`IntegrationEvidence`** -- for candidates selected via MOFA+ factor
  loading (RNA-seq/methylation/CNV).
- **`RecurrenceEvidence`** -- for candidates selected via mutation
  recurrence filtering, with an optional generic Fisher's-exact
  categorical association test against any supplied outcome.
- **`SurvivalEvidence`** -- computed identically for *every* candidate via
  univariate Cox PH, regardless of which pathway selected it.
- **`DruggabilityEvidence`** -- schema only here; populated by
  `feat/drug-target-scoring`.

## Design goals

- No hidden magic numbers: every weight in the composite score is an
  explicit, documented, user-overridable parameter (`ScoreWeights`).
- Deterministic given identical inputs (no unseeded randomness).
- Missing evidence axes are handled by weight renormalization, not
  penalization -- a mutation candidate structurally has no
  `IntegrationEvidence` and must not be scored as if that's a deficiency.
- Zero dependency on the rest of the OncoCartograph pipeline's
  internals -- every module here operates on plain pandas
  Series/DataFrames and dataclasses.

## API

```python
from oncocartograph.scoring import (
    fit_univariate_cox,
    screen_survival_associations,
    fisher_exact_association,
    BiomarkerEvidence,
    IntegrationEvidence,
    RecurrenceEvidence,
    DruggabilityEvidence,
    ScoreWeights,
    composite_biomarker_score,
)

# Survival association -- identical treatment for continuous or binary values.
screen = screen_survival_associations(expression_matrix, duration, event)

# Assemble evidence for one candidate (any subset of these may be present).
evidence = BiomarkerEvidence(
    candidate_id="TP53",
    survival=fit_univariate_cox(values, duration, event),
    recurrence=RecurrenceEvidence(n_patients_mutated=10, cohort_size=122),
)

score = composite_biomarker_score(evidence, weights=ScoreWeights())
```

See `docs/methods.md` §4 for the full statistical justification (why Cox
PH rather than Fine-Gray, the FDR correction approach, and the composite
formula's weighting rationale) and §4.4 for this project's real scoring
run result, including a real bug this package's tests caught (a
significant fraction of sparse mutation fits produce non-finite Cox
estimates that must be excluded, not passed through as evidence).
