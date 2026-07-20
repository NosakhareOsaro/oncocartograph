# Global Talent visa evidence mapping

This document maps OncoCartograph's actual, real project history — not a
hypothetical one — to the evidence checklist established at the start of
this project. Every claim below is checked against the real repository
state (commit hashes, file paths, live API/CI results) rather than
recalled from memory; where something has not happened yet, that is
stated plainly rather than implied.

Repository: https://github.com/NosakhareOsaro/oncocartograph (public, MIT license)

## 1. Personal, significant technical contribution

**The standalone composite scoring package** (`src/oncocartograph/scoring/`)
is the strongest evidence here. It combines three independently-sourced
evidence types — survival association (Cox PH), druggability (Open
Targets/ChEMBL), and selection pathway (MOFA+ factor loading or mutation
recurrence) — into one composite score, with:

- A **mechanically enforced zero-cross-import architecture**:
  `tests/scoring/test_decoupling.py` parses the package's AST and fails
  the test suite if `oncocartograph.scoring` imports anything from the
  rest of the project, so the package is verifiably extractable to a
  standalone PyPI package, not just described as modular.
- **100% test coverage** on every module in the package (`survival.py`,
  `association.py`, `evidence.py`, `composite.py`), including a
  real bug found and fixed by running the actual screen against real
  data: `lifelines` can return non-finite Cox statistics for sparse
  binary covariates, which an initial NaN-only check caught in only
  ~11% of affected genes; broadening the check to all four fit
  statistics caught the true 84% (`src/oncocartograph/scoring/survival.py`,
  documented in `docs/methods.md` §4.2).
- A **documented, cited design rationale** for every non-obvious choice
  (Cox PH over Fine-Gray as a data-availability constraint, not a
  preference; the renormalized-weighted-average composite formula;
  floors on protective associations) in
  [`docs/adr/0007-survival-methodology-and-composite-score.md`](adr/0007-survival-methodology-and-composite-score.md).

**The pre-registered, falsifiable external validation** (`src/oncocartograph/validation/`,
[`docs/adr/0009-external-replication-methodology-and-result.md`](adr/0009-external-replication-methodology-and-result.md))
is the second pillar: a validation design that fixed its
success/failure criterion *before* running the real analysis, and
reported the resulting failure honestly rather than redefining the
criterion afterward. See §3 below for why this specifically demonstrates
quality, not just effort.

**Two real mistakes, caught and corrected through direct verification
against real data rather than assumption**, are evidence of engineering
judgment, not just code volume:

- **The GISTIC2 assumption.** The original project plan assumed
  TCGA-BRCA copy number would ship as GISTIC2 thresholded categorical
  calls. Inspecting the actual downloaded files during the live GDC pull
  showed this was wrong — GDC's current harmonized pipeline reports
  absolute integer total copy number from up to four different calling
  workflows per patient. Corrected in commit `8ac2bc8`
  (`feat(preprocessing): add copy number workflow resolution and CN
  transform`) and documented in
  [`docs/adr/0005-copy-number-workflow-and-transform.md`](adr/0005-copy-number-workflow-and-transform.md).
- **The MOFA+ view-ordering bug.** `mofapy2` silently re-sorts view
  names alphabetically internally, independent of input order; the
  likelihoods list was initially built in dict-insertion order, which
  would have silently assigned the wrong likelihood model (e.g.
  Bernoulli) to the wrong view. This surfaced immediately as an
  `AssertionError` ("Data must be binary") while writing the first
  end-to-end integration test — before it ever reached a real training
  run — and was fixed by explicitly sorting view names to match
  mofapy2's internal order, with a regression test whose fixture
  deliberately has mismatched insertion order
  (`src/oncocartograph/integration/mofa.py:106-110`, commit `fc23118`,
  `tests/integration/test_mofa.py`).

Both are evidence that this project's stated engineering discipline
("verify against real data/APIs before writing code") was actually
practiced, not just asserted: both bugs were caught by testing against
real data or a real first integration test, not found later in
production.

## 2. Evidence of work adopted/usable by others

Not a personal script: this is a packaged, installable, tested,
documented, continuously-integrated project.

- **Installable**: `pip install -e ".[dev,workflow]"` (`pyproject.toml`);
  MIT-licensed, public GitHub repository.
- **Tested**: 207 tests, 99% overall coverage, 100% on every module
  added since `feat/scoring-package` (`pyproject.toml`'s `pytest`/
  `coverage` configuration; verified in the Python 3.11 Docker image
  that mirrors CI, not just the local dev environment).
- **Continuously integrated**: `.github/workflows/ci.yml` runs lint
  (ruff+black), type checking (mypy --strict), tests, a docs-presence
  check, and — as of this work package — a Snakemake DAG dry-run
  validation, on every push/PR to `main`. As of this writing, CI is
  genuinely green on the latest `main` commit (verified via the GitHub
  Actions API, not assumed from a local exit code).
- **Documented**: `README.md`, `docs/methods.md` (full methodology),
  `docs/data_sources.md` (every dataset's exact accession/query/download
  date/license), nine architecture decision records in `docs/adr/`, and
  a real, runnable Snakemake pipeline (`workflows/Snakefile`) that
  reproduces every result end-to-end from raw source data — not just
  described in prose, but wired into one DAG and verified against real
  data during this work package (dry-run resolves cleanly from a
  completely clean `data/` state; every fast/already-cached rule was
  actually executed and reproduced the already-published real numbers).
- **Reviewable**: every commit is small, atomic, and states *why* a
  change was made, not just what changed; every milestone
  (`v0.1.0`–`v0.4.0`) is an annotated, pushed, independently-verified git
  tag.

## 3. Evidence of external validation of quality

This is the criterion most worth being precise about, because the real
result is not a straightforward success.

**What succeeded**: the validation *design* meets a high bar for
rigor — a real independent cohort (GSE96058/SCAN-B, not a synthetic or
held-out TCGA split), chosen specifically to avoid a cross-platform
confound (ADR 0003), with a primary success/failure criterion
(direction-concordance vs. chance, one-sided exact binomial test,
α=0.05) implemented and committed (`src/oncocartograph/validation/replication.py`,
commit `b70d037`) *before* the real analysis was ever run against real
GSE96058 data — the commit that fixed the test came first, and the
commit that reports the result (`b517483`, which also adds the prose
ADR) came strictly after it, on the same branch. A secondary, scope-reduced check against
Burstein et al. (2015)'s published TNBC subtype biology — five genes
with well-documented literature directions (AR, PTEN, CD274, PDCD1,
CTLA4) — **passed 5/5**, independently confirmed by a second full
pipeline run during this work package.

**What did not succeed**: the primary, pre-registered quantitative
replication criterion **failed** — 41.3% direction concordance (below
the 50% chance rate), p=0.973, reproduced a second time at 39.8%/p=0.988
via the newly-built Snakemake pipeline. Composite-scored candidates from
this pipeline are not, on the current evidence, externally validated
biomarkers.

**Why this is still evidence of quality, not a weakness to minimize**: a
validation step is only meaningful evidence if it was capable of
reporting failure, and this one did — honestly, in the abstract and
conclusion of `docs/manuscript.md`, not buried in a limitations
appendix. The alternative — a validation designed so it could only ever
report success, or a negative result quietly reframed after the fact —
would be worth far less as evidence of rigor, however much better it
would read. The Burstein check's 5/5 result demonstrates the validation
*machinery* (cohort definition, expression parsing, Cox fitting) works
correctly; the primary result demonstrates the *candidates* it was
testing did not hold up outside their discovery cohort. Reporting both,
plainly, in the same documents, is the actual evidence of quality on
offer here.

## 4. A citable artifact

- **Preprint**: [`docs/manuscript.md`](manuscript.md) — a real write-up
  of the full pipeline and its results, including the null validation
  result stated in the abstract and conclusion, not just the methods
  section.
- **`CITATION.cff`**: exists at the repository root, kept in sync with
  the current release version.
- **Zenodo DOI**: **does not exist yet.** Minting one requires linking
  this GitHub repository to a Zenodo account and creating a GitHub
  Release (Zenodo archives the release and assigns a DOI automatically)
  — a real action for you to take, not something fabricated here. See
  the recognition-potential section below for the concrete next step.

## 5. Recognition potential

Realistic next steps, none of which have happened yet — listed as
concrete options for you to pursue, not claimed as already achieved:

1. **Mint a real Zenodo DOI for the `v1.0.0` release.** Connect
   https://zenodo.org to your GitHub account (Zenodo → GitHub → toggle
   this repository on), then create a GitHub Release for the `v1.0.0`
   tag; Zenodo archives it and assigns a DOI automatically within
   minutes. This turns `docs/manuscript.md` and the codebase into a
   permanently citable artifact with a real identifier, closing the gap
   noted in §4.
2. **Extract `oncocartograph.scoring` to a standalone PyPI package.**
   Its zero-cross-import architecture is already mechanically verified
   (`tests/scoring/test_decoupling.py`); publishing it separately (e.g.
   as `biomarker-composite-score` or similar) would let it be `pip
   install`-able and cited independently of the full TNBC pipeline,
   directly relevant to the "adopted/usable by others" criterion.
3. **Submit the manuscript for external technical review or a
   preprint server.** `docs/manuscript.md` is already structured as a
   preprint; posting it to bioRxiv (or a similar venue) or requesting a
   review from a named bioinformatics researcher working in
   multi-omics integration or TNBC would convert the current
   self-assessed rigor into third-party-recognized rigor — the kind of
   recognition this checklist's other criteria cannot substitute for.
