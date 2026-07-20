"""Data ingestion: TCGA-BRCA TNBC sub-cohort definition and multi-omics data retrieval.

Responsible for querying GDC/TCGAbiolinks and GEO for RNA-seq, methylation,
copy number, mutation, and clinical data, and for producing the auditable
TNBC sub-cohort sample list described in ``docs/methods.md``.
"""
