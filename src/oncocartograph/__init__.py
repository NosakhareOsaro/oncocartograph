"""OncoCartograph: multi-omics integration and druggable biomarker prioritisation for TNBC.

This package implements a reproducible pipeline that integrates RNA-seq,
DNA methylation, copy number, and somatic mutation data from TCGA-BRCA
samples restricted to a rigorously defined triple-negative breast cancer
(TNBC) sub-cohort, using MOFA+ for multi-omics factor analysis and a
standalone composite scoring package to prioritise druggable biomarkers.

See ``docs/methods.md`` for the full methodology and ``docs/data_sources.md``
for exact dataset accessions and versions used.
"""

__version__ = "0.3.0"
