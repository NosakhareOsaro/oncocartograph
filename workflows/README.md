# workflows/

Pipeline orchestration via Snakemake (see
[`docs/adr/0002-workflow-engine-choice.md`](../docs/adr/0002-workflow-engine-choice.md)
for why Snakemake over Nextflow).

`Snakefile` and per-stage conda environment files will be added
incrementally as each work package lands, starting with
`feat/data-ingestion`. Nothing is stubbed out here in advance of the rules
it would actually run.
