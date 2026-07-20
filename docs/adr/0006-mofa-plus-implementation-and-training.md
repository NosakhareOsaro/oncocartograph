# ADR 0006: MOFA+ implementation (mofapy2/mofax, no R) and training configuration

## Status

Accepted (2026-07-20)

## Context

The project architecture names "MOFA+ wrapper" as the integration stage.
MOFA+ is most commonly used via its R interface (the `MOFA2` R package),
and an earlier note in this repo's `Dockerfile` assumed an R + MOFA2
environment would be needed for this stage, extending the Python-only
image built for ingestion/preprocessing.

Before writing any integration code, this was checked directly rather
than assumed: `mofapy2` (the official Python package for training MOFA+
models) and `mofax` (the companion Python package, same MOFA+ team, for
reading a trained model's factor values, feature weights, and variance
explained) were installed and tested. Both work standalone — training and
interpretation are both fully achievable in pure Python, with no R
dependency at all for this stage.

Two further things were verified directly (not assumed) before finalizing
the design:
1. `mofapy2`'s long-format DataFrame input (`sample, feature, view, group,
   value`) correctly handles patients **entirely absent** from a view
   (tested: a 10-sample view A alongside a 7-sample view B trained
   successfully, `N=10` reported overall).
2. It also correctly handles **scattered missing values** within an
   otherwise-present sample/view (tested: dropping specific
   sample-feature rows, not just whole samples, from the long-format
   input trains successfully) — needed since methylation and copy number
   both have residual missingness after their preprocessing filters.

## Decision

**Implementation:** `mofapy2` for training, `mofax` for reading/
interpreting the trained model. No R dependency anywhere in this stage.
The `Dockerfile`'s earlier note about needing an R + MOFA2 environment is
retracted.

**View → likelihood mapping:**

| View | Likelihood | Features (post-selection) | Patients |
|---|---|---|---|
| RNA-seq (VST) | Gaussian | 2,000 | 142 |
| Methylation (M-value) | Gaussian | 5,000 | 104 |
| Copy number (relative log2) | Gaussian | 2,000 (see below) | 142 |
| Mutation (binary) | Bernoulli | 845 | 122 |

**Sample handling:** union of all patients across views (not a forced
complete-case subset), built as a long-format DataFrame with missing
sample/view combinations simply absent as rows, per the verification
above. `scale_views=True` (unit-variance scaling per view) to prevent a
naturally higher-magnitude view from dominating the model regardless of
biological signal.

**Copy number feature cap:** the CNV view initially had ~60,623 genes
with no selection step (a real gap found while preparing this ADR, by
running full preprocessing against the entire real cohort rather than
small subsets) — a 30x imbalance against the other Gaussian views. Capped
to the top 2,000 most-variable genes, matching RNA-seq's count.

**Training hyperparameters:**
- `factors=15`: middle-of-the-road starting count for a ~100-150 sample
  cohort, per common practice in published MOFA+ multi-omics analyses at
  this scale. ARD priors on feature weights (`ard_weights=True`, the
  mofapy2 default) provide sparsity; factors explaining <2% variance in
  every view are excluded during interpretation (not training) — see
  `docs/methods.md` §3.
- `convergence_mode="slow"`: tightest ELBO convergence tolerance mofapy2
  offers. Chosen over the faster/looser default because this is a final,
  citable analysis, not exploratory work.
- `seed`: taken from `Settings.random_seed` (project-wide default 42),
  logged at training time — consistent with the "every stochastic step
  logs its seed" engineering standard.

## Alternatives considered

**R + MOFA2** for training and/or interpretation. Rejected once `mofapy2`
+ `mofax` were confirmed sufficient — avoids introducing a second
language/runtime into a project that has otherwise stayed pure Python
except where genuinely unavoidable (this stage turned out not to require
that trade-off at all).

**Forcing a complete-case cohort** across all four views (~intersection
of 142/104/142/122 patients, likely well under 100). Rejected — MOFA+'s
core value proposition is factorizing partially-overlapping views, and
this was directly verified to work as expected.

## Consequences

- `Dockerfile` no longer needs an R layer for MOFA+; it stays Python-only
  end-to-end through this stage.
- The methylation view's smaller N (104) means the mutation and
  methylation views individually support less statistical weight in any
  given factor than RNA-seq/CNV (142 patients each) — noted in
  `docs/methods.md` §8 Limitations, not silently absorbed.
