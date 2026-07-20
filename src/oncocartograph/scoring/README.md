# oncocartograph.scoring

Standalone composite biomarker scoring package for prioritising druggable
biomarkers from multi-omics survival-association evidence.

> **Status:** scaffolded, not yet implemented. This README documents the
> intended design so the package's scope is clear from commit 1. It will be
> filled in with real usage examples, API docs, and validation results as
> part of the `feat/scoring-package` work package.

## Purpose

This package is designed to be usable, tested, and citable independently of
the rest of the OncoCartograph pipeline. It takes per-biomarker statistical
evidence (e.g. survival hazard ratios and multiple-testing-corrected
p-values from a Cox or Fine-Gray model) and druggability evidence (Open
Targets tractability data, ChEMBL bioactivity data) and produces a single,
documented composite priority score per candidate biomarker.

## Design goals

- No hidden magic numbers: every weight and threshold in the composite score
  is an explicit, documented, user-overridable parameter.
- Deterministic given identical inputs (no unseeded randomness).
- >90% unit test coverage, including edge cases (missing druggability
  evidence, ties, all-negative survival evidence).
- Zero dependency on the rest of the OncoCartograph pipeline's internals —
  it operates on plain pandas DataFrames with a documented schema, so it can
  be extracted and published to PyPI independently if desired.

## Planned API (subject to change until `feat/scoring-package` lands)

```python
from oncocartograph.scoring import composite_biomarker_score

result = composite_biomarker_score(
    survival_evidence=survival_df,   # gene, hazard_ratio, p_adj, ...
    druggability_evidence=drug_df,   # gene, tractability_score, chembl_max_phase, ...
    weights=ScoreWeights(survival=0.6, druggability=0.4),
)
```

Full statistical justification for the scoring formula will be documented in
`docs/methods.md` and in this README once implemented.
