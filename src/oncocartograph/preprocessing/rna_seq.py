"""RNA-seq preprocessing: low-expression filtering and DESeq2-style normalization.

Consumes GDC's STAR-Counts ``augmented_star_gene_counts.tsv`` files
(confirmed format: a leading ``#``-commented gene-model line, a header
row, four ``N_*`` alignment-summary rows to discard, then one row per
gene with raw counts in the ``unstranded`` column). Normalization uses
``pydeseq2`` (a Python port of DESeq2) for median-of-ratios size-factor
normalization and variance-stabilizing transformation (VST), per
ADR-confirmed choice over a simpler log2(TPM+1) approach -- size factors
correct for library-composition differences more robustly than TPM alone.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference

#: gene_id prefix marking STAR's alignment-summary rows, not genes.
_STAR_SUMMARY_ROW_PREFIX = "N_"


def read_star_counts(path: Path) -> pd.Series:
    """Read a STAR-Counts file into a gene_id -> raw count Series.

    Args:
        path: Path to a downloaded ``*.rna_seq.augmented_star_gene_counts.tsv``
            file.

    Returns:
        A Series of raw (``unstranded``) integer counts indexed by
        versioned Ensembl gene ID, with STAR's ``N_*`` alignment-summary
        rows removed.
    """
    table = pd.read_csv(path, sep="\t", comment="#", header=0)
    table = table[~table["gene_id"].str.startswith(_STAR_SUMMARY_ROW_PREFIX)]
    return table.set_index("gene_id")["unstranded"]


def build_counts_matrix(resolved_files: dict[str, Path]) -> pd.DataFrame:
    """Build a gene x patient raw counts matrix from per-patient STAR-Counts files.

    Args:
        resolved_files: Mapping of case_id to that patient's resolved
            STAR-Counts file path.

    Returns:
        A DataFrame indexed by gene_id, one column per case_id, with raw
        integer counts.
    """
    columns = {case_id: read_star_counts(path) for case_id, path in resolved_files.items()}
    return pd.DataFrame(columns)


def filter_low_expression(
    counts: pd.DataFrame, *, min_cpm: float = 1.0, min_fraction_samples: float = 0.1
) -> pd.DataFrame:
    """Drop genes with low expression across the cohort (edgeR filterByExpr-style).

    Args:
        counts: A gene x patient raw counts matrix.
        min_cpm: Minimum counts-per-million a gene must reach in a sample
            for that sample to count towards ``min_fraction_samples``.
        min_fraction_samples: Minimum fraction of samples that must meet
            ``min_cpm`` for a gene to be retained.

    Returns:
        The subset of rows (genes) meeting the expression threshold.
    """
    library_sizes = counts.sum(axis=0)
    cpm = counts.div(library_sizes, axis=1) * 1_000_000
    fraction_expressed = (cpm >= min_cpm).mean(axis=1)
    return counts.loc[fraction_expressed >= min_fraction_samples]


def normalize_and_vst(counts: pd.DataFrame) -> pd.DataFrame:
    """Apply DESeq2 median-of-ratios normalization and variance-stabilizing transform.

    Args:
        counts: A gene x patient raw (post-filtering) counts matrix.

    Returns:
        A gene x patient DataFrame of VST-transformed values, same shape
        and index/column labels as ``counts``. Uses ``design="~1"``
        (intercept only) since this pipeline does not test a specific
        contrast here -- MOFA+ discovers structure from the resulting
        continuous values directly, rather than from a pre-specified
        differential expression design.
    """
    samples_by_genes = counts.T
    metadata = pd.DataFrame(index=samples_by_genes.index)
    dds = DeseqDataSet(
        counts=samples_by_genes,
        metadata=metadata,
        design="~1",
        inference=DefaultInference(n_cpus=1),
        quiet=True,
    )
    dds.vst(use_design=False)
    vst_values = pd.DataFrame(
        dds.layers["vst_counts"],
        index=samples_by_genes.index,
        columns=samples_by_genes.columns,
    )
    return vst_values.T


def select_top_variable_genes(matrix: pd.DataFrame, n: int) -> pd.DataFrame:
    """Select the n most variable genes by variance across patients.

    Args:
        matrix: A gene x patient VST-transformed matrix.
        n: Number of top-variable genes to keep. If ``matrix`` has fewer
            than ``n`` rows, all rows are kept.

    Returns:
        The subset of rows with the highest variance, in descending order
        of variance.
    """
    variances = matrix.var(axis=1, skipna=True)
    top_genes = variances.sort_values(ascending=False).head(n).index
    return matrix.loc[top_genes]
