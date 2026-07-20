"""Mutation preprocessing: non-synonymous filtering and binary gene x patient matrix.

Converts per-patient MC3-derived masked somatic mutation MAFs (one file
per patient, already tumor/normal-paired by the calling pipeline -- see
``oncocartograph.preprocessing.sample_manifest``) into a single binary
gene x patient matrix suitable for a Bernoulli MOFA+ view: 1 if that gene
carries at least one non-synonymous variant in that patient, 0 otherwise.

Real ``Variant_Classification`` values were confirmed against downloaded
MAF files during this work package (Missense_Mutation, Silent,
Nonsense_Mutation, Frame_Shift_Del, Splice_Site, etc.), matching the
standard TCGA MAF vocabulary this module's filter list assumes.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

#: Variant_Classification values treated as functionally non-synonymous.
#: An explicit allow-list (not a deny-list) so an unrecognised future
#: category defaults to excluded rather than silently included.
NON_SYNONYMOUS_CLASSIFICATIONS: frozenset[str] = frozenset(
    {
        "Missense_Mutation",
        "Nonsense_Mutation",
        "Nonstop_Mutation",
        "Frame_Shift_Del",
        "Frame_Shift_Ins",
        "In_Frame_Del",
        "In_Frame_Ins",
        "Splice_Site",
        "Translation_Start_Site",
    }
)


def read_non_synonymous_genes(path: Path) -> set[str]:
    """Read a MAF file and return the set of genes with a non-synonymous variant.

    Args:
        path: Path to a (optionally gzip-compressed) MAF file. Comment
            lines beginning with ``#`` are skipped automatically.

    Returns:
        The set of distinct ``Hugo_Symbol`` values among variants whose
        ``Variant_Classification`` is in :data:`NON_SYNONYMOUS_CLASSIFICATIONS`.
    """
    maf = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        usecols=["Hugo_Symbol", "Variant_Classification"],
        dtype=str,
    )
    non_synonymous = maf[maf["Variant_Classification"].isin(NON_SYNONYMOUS_CLASSIFICATIONS)]
    return set(non_synonymous["Hugo_Symbol"].dropna())


def build_mutation_matrix(resolved_files: dict[str, Path]) -> pd.DataFrame:
    """Build a binary gene x patient mutation matrix from per-patient MAFs.

    Args:
        resolved_files: Mapping of case_id to that patient's resolved MAF
            file path.

    Returns:
        A DataFrame indexed by ``Hugo_Symbol`` gene name, one column per
        case_id, with 1 where that patient has a non-synonymous variant in
        that gene and 0 otherwise (no ``NaN`` -- absence of a mutation
        call is a genuine 0, unlike the other omic views).
    """
    genes_by_case = {
        case_id: read_non_synonymous_genes(path) for case_id, path in resolved_files.items()
    }
    all_genes = sorted(set().union(*genes_by_case.values())) if genes_by_case else []
    matrix = pd.DataFrame(
        {
            case_id: [1 if gene in genes else 0 for gene in all_genes]
            for case_id, genes in genes_by_case.items()
        },
        index=all_genes,
    )
    return matrix


def filter_by_recurrence(matrix: pd.DataFrame, min_patients: int) -> pd.DataFrame:
    """Keep only genes mutated in at least ``min_patients`` patients.

    Args:
        matrix: A binary gene x patient matrix, e.g. from
            :func:`build_mutation_matrix`.
        min_patients: Minimum number of patients a gene must be mutated in
            to be retained.

    Returns:
        The subset of rows meeting the recurrence threshold, in the same
        column order.
    """
    recurrence = matrix.sum(axis=1)
    return matrix.loc[recurrence >= min_patients]
