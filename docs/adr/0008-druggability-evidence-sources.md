# ADR 0008: Druggability evidence sources and combination

## Status

Accepted (2026-07-20)

## Context

`oncocartograph.scoring.evidence.DruggabilityEvidence` was defined in
`feat/scoring-package` as a schema-only placeholder (`tractability_score`,
`chembl_max_phase`), deferred to this work package to populate with real
Open Targets and ChEMBL data. Three things needed real verification
before writing any code, all confirmed via live queries against the
actual APIs rather than assumed from documentation:

1. Open Targets' `tractability` field is ~28 boolean buckets across four
   modalities (SM/AB/PR/OC), not a pre-computed [0, 1] score -- this
   project needed its own collapse.
2. ChEMBL's free-text target search is unreliable for exact gene
   matching (searching "TP53" returns "TP53-binding protein 1" as the
   top hit). Exact resolution requires going via UniProt accession
   (obtained from Open Targets) and ChEMBL's
   `target_components__accession__in` filter restricted to
   `target_type=SINGLE PROTEIN`.
3. Our real candidate set has three different identifier types: versioned
   Ensembl IDs (RNA-seq/CNV), gene symbols (mutation), and CpG probe IDs
   (methylation) -- only the first two map to genes without an additional
   reference dataset.

## Decision

**Identifier resolution:** strip version suffixes for RNA-seq/CNV
candidates (already Ensembl IDs); resolve mutation candidates' gene
symbols via Open Targets' `mapIds` query (batch, returns an empty hit
list rather than erroring for unresolvable symbols).

**Methylation candidates are out of scope for this work package.**
CpG probe IDs (e.g. `cg00000029`) are not gene identifiers; mapping a
probe to its nearest/associated gene requires the Illumina 450K manifest
(~450K rows), a new reference data dependency not part of any prior work
package, with its own non-trivial judgment call (nearest TSS? within gene
body? multiple candidate genes per probe?). Deferred rather than
silently worked around with an approximate mapping.

**Tractability scoring:** a 3-tier collapse of the raw buckets --
Approved Drug (any modality) = 1.0, clinical-stage evidence (Advanced
Clinical or Phase 1 Clinical, any modality) = 0.66, any other tractability
bucket = 0.33, no evidence = 0.0. Implemented in
`oncocartograph.drug_targets.tractability.score_tractability`.

**Combining Open Targets tractability with ChEMBL max phase:**
`max(tractability_tier_score, chembl_max_phase / 4)`. Either signal alone
can indicate a real drug exists or is in development; taking the max
rather than an average means a target strong on one axis but merely
average on the other still gets credit for the strong signal.
`chembl_max_phase` remains on `DruggabilityEvidence` separately for
transparent reporting (e.g. manuscript tables) even though it is folded
into `tractability_score` for the actual composite score, which reads
only that one field.

## Alternatives considered

**Finer-grained tractability tiers** distinguishing structural evidence
quality (e.g. ranking "Structure with Ligand" above "Druggable Family").
Rejected for this iteration -- would require classifying ~20 individual
bucket labels into sub-tiers, each choice more arbitrary than the
clinical-stage/approved-drug/any-evidence/none split, which follows
Open Targets' own natural "how close to being a real drug" ordering.

**Keeping `chembl_max_phase` fully separate from the composite score**
(informational only). Rejected -- would mean a target with strong ChEMBL
clinical evidence but average Open Targets tractability buckets gets no
credit for it in the actual ranking, understating real druggability
evidence the project has in hand.

**Adding the Illumina 450K manifest now** to unblock methylation
candidates. Rejected for this iteration given the added scope (new data
source, its own provenance/licensing considerations, a probe-to-gene
assignment method needing its own justification) -- a candidate follow-up
if methylation-derived biomarkers turn out to matter significantly.

## Consequences

- Methylation-derived candidates (229/576 of the real MOFA+-derived
  candidate set, the largest single view) will retain
  `druggability=None` and be scored on survival + selection-pathway
  evidence only (renormalized weights) until a probe-to-gene mapping is
  added.
- `docs/data_sources.md` records the exact Open Targets/ChEMBL query
  patterns used, per the project's provenance-logging convention.
