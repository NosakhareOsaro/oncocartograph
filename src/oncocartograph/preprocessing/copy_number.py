"""Copy number preprocessing: workflow resolution and relative-to-diploid transform.

TCGA-BRCA gene-level copy number in the current GDC harmonized pipeline is
produced by up to four distinct calling workflows per patient (confirmed
via a live GDC metadata query, not assumed): ASCAT3 (current), ASCAT2
(superseded by ASCAT3), ABSOLUTE LiftOver (a different algorithm
entirely), and AscatNGS (a WGS-specific ASCAT variant used for the small
subset of cases with whole-genome rather than SNP-array data). Per
``docs/adr/0005-copy-number-workflow-and-transform.md``, ASCAT3 is
preferred, falling back through the others in that order.

Each file reports **absolute integer total copy number per gene** (e.g. 2
= diploid, 0 = homozygous deletion), not GISTIC2's categorical thresholded
calls -- correcting an assumption in the original data plan (see the ADR).
This module converts that to a relative-to-diploid log2 value for use as a
continuous (Gaussian) MOFA+ view, without correcting for tumor
purity/ploidy (a documented limitation, not a silent gap).
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any, cast

import numpy as np
import pandas as pd

from oncocartograph.preprocessing.sample_manifest import (
    ResolvedFile,
    default_tie_break,
    resolve_primary_tumor_files,
)

#: Workflow preference order when a patient has more than one CNV calling
#: workflow available, most to least preferred. See ADR 0005.
CNV_WORKFLOW_PRIORITY: tuple[str, ...] = ("ASCAT3", "ASCAT2", "ABSOLUTE LiftOver", "AscatNGS")

#: GDC file fields required to resolve CNV files by workflow.
CNV_MANIFEST_FIELDS = (
    "file_id",
    "file_name",
    "cases.case_id",
    "cases.samples.sample_type",
    "analysis.workflow_type",
)


def workflow_priority_tie_break(
    priority: Sequence[str] = CNV_WORKFLOW_PRIORITY,
) -> Callable[[Sequence[dict[str, Any]]], dict[str, Any]]:
    """Build a tie-break function that prefers files by ``analysis.workflow_type``.

    Args:
        priority: Workflow type names in preference order, most preferred
            first.

    Returns:
        A callable suitable for :func:`resolve_primary_tumor_files`'s
        ``tie_break`` parameter: given multiple file hits for one patient,
        returns the one whose workflow type appears earliest in
        ``priority``. Falls back to :func:`default_tie_break` if none of
        the hits' workflow types appear in ``priority`` at all.
    """

    def _tie_break(hits: Sequence[dict[str, Any]]) -> dict[str, Any]:
        for workflow in priority:
            for hit in hits:
                if hit.get("analysis", {}).get("workflow_type") == workflow:
                    return hit
        return default_tie_break(hits)

    return _tie_break


def resolve_copy_number_files(file_hits: Sequence[dict[str, Any]]) -> dict[str, ResolvedFile]:
    """Resolve CNV file hits to one Primary Tumor, workflow-preferred file per patient.

    Args:
        file_hits: GDC file hits with fields including
            :data:`CNV_MANIFEST_FIELDS`.

    Returns:
        A dict mapping case_id to the chosen :class:`ResolvedFile`.
    """
    return resolve_primary_tumor_files(file_hits, tie_break=workflow_priority_tie_break())


def read_gene_level_copy_number(path: Path) -> pd.Series:
    """Parse a GDC gene-level copy number TSV into a gene_id -> copy_number Series.

    Args:
        path: Path to a downloaded ``*.gene_level_copy_number*.tsv`` file.

    Returns:
        A Series indexed by ``gene_id`` (versioned Ensembl ID) with the
        absolute integer total copy number per gene. Genes outside any
        called segment are ``NaN``, not imputed.
    """
    table = pd.read_csv(path, sep="\t", usecols=["gene_id", "copy_number"])
    return table.set_index("gene_id")["copy_number"]


#: Assumed baseline (non-tumor, non-ploidy-corrected) diploid copy number.
DIPLOID_COPY_NUMBER = 2


def relative_log2_copy_number(copy_number: pd.Series) -> pd.Series:
    """Convert absolute copy number to a relative-to-diploid log2 value.

    Args:
        copy_number: Absolute integer total copy number per gene (2 =
            diploid), as returned by :func:`read_gene_level_copy_number`.

    Returns:
        ``log2((copy_number + 1) / (DIPLOID_COPY_NUMBER + 1))``. Adding 1
        to both numerator and the diploid reference (rather than just the
        numerator) ensures diploid (2) maps to exactly 0 while CN=0
        (homozygous deletion) stays finite instead of ``-inf``. Losses are
        negative, gains are positive. Does not correct for tumor
        ploidy/purity -- see the limitations note in
        ``docs/adr/0005-copy-number-workflow-and-transform.md``.
    """
    return cast(pd.Series, np.log2((copy_number + 1) / (DIPLOID_COPY_NUMBER + 1)))


def build_copy_number_matrix(resolved_files: dict[str, Path]) -> pd.DataFrame:
    """Build a gene x patient matrix of relative log2 copy number.

    Args:
        resolved_files: Mapping of case_id to the local path of that
            patient's resolved copy number file (e.g. from combining
            :func:`resolve_copy_number_files` with the downloaded file
            locations).

    Returns:
        A DataFrame indexed by gene_id, one column per case_id, with
        relative log2 copy number values (``NaN`` where a gene had no
        called value for that patient).
    """
    columns = {}
    for case_id, path in resolved_files.items():
        raw = read_gene_level_copy_number(path)
        columns[case_id] = relative_log2_copy_number(raw)
    return pd.DataFrame(columns)
