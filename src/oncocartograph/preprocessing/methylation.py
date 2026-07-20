"""Methylation preprocessing: beta value loading, M-value transform, and probe selection.

Consumes GDC's already-processed SeSAMe beta value files
(``*.methylation_array.sesame.level3betas.txt``, per the
``feat/data-ingestion`` filter fix that excludes raw IDAT intensity
files). These are two-column, header-less, tab-delimited files
(probe_id, beta_value) confirmed against real downloaded files during
this work package -- SeSAMe's own masking already flags unreliable
probes as missing (about 14% NaN in the files inspected), so this module
does not re-derive a separate cross-reactive-probe blacklist.
"""

from __future__ import annotations

from pathlib import Path
from typing import cast

import numpy as np
import pandas as pd

#: Beta values are clipped to this range before the M-value transform to
#: avoid log(0) at the boundaries (beta is nominally in [0, 1]).
_BETA_CLIP_EPSILON = 1e-6


def read_beta_values(path: Path) -> pd.Series:
    """Read a SeSAMe level3betas.txt file into a probe_id -> beta value Series.

    Args:
        path: Path to a downloaded ``*.sesame.level3betas.txt`` file (no
            header row: probe_id, beta_value).

    Returns:
        A Series indexed by Illumina probe ID (e.g. ``cg00000029``), with
        missing/masked probes as ``NaN``.
    """
    table = pd.read_csv(path, sep="\t", header=None, names=["probe_id", "beta_value"])
    return table.set_index("probe_id")["beta_value"]


def beta_to_m_value(beta: pd.Series) -> pd.Series:
    """Convert beta values to M-values via the logit transform.

    Args:
        beta: Beta values in (nominally) [0, 1].

    Returns:
        ``log2(beta / (1 - beta))``, with beta clipped to
        ``[eps, 1 - eps]`` first to keep the result finite at the
        boundaries. M-values are the recommended representation for
        regression/factor-analysis use (Du et al. 2010, *BMC
        Bioinformatics* 11:587), since beta values are bounded and
        heteroscedastic near 0 and 1.
    """
    clipped = beta.clip(lower=_BETA_CLIP_EPSILON, upper=1 - _BETA_CLIP_EPSILON)
    return cast(pd.Series, np.log2(clipped / (1 - clipped)))


def build_methylation_matrix(resolved_files: dict[str, Path]) -> pd.DataFrame:
    """Build a probe x patient M-value matrix from per-patient beta value files.

    Args:
        resolved_files: Mapping of case_id to that patient's resolved
            beta value file path.

    Returns:
        A DataFrame indexed by probe_id, one column per case_id, with
        M-values (not raw beta values).
    """
    columns = {
        case_id: beta_to_m_value(read_beta_values(path)) for case_id, path in resolved_files.items()
    }
    return pd.DataFrame(columns)


def filter_by_missingness(matrix: pd.DataFrame, max_missing_fraction: float) -> pd.DataFrame:
    """Drop probes with too much missingness across the profiled patients.

    Args:
        matrix: A probe x patient matrix, e.g. from
            :func:`build_methylation_matrix`.
        max_missing_fraction: Maximum allowed fraction of patients with a
            missing (``NaN``) value for a probe to be retained.

    Returns:
        The subset of rows meeting the missingness threshold.
    """
    missing_fraction = matrix.isna().mean(axis=1)
    return matrix.loc[missing_fraction <= max_missing_fraction]


def select_top_variable_probes(matrix: pd.DataFrame, n: int) -> pd.DataFrame:
    """Select the n most variable probes by variance across patients.

    Args:
        matrix: A probe x patient M-value matrix.
        n: Number of top-variable probes to keep. If ``matrix`` has fewer
            than ``n`` rows, all rows are kept.

    Returns:
        The subset of rows with the highest variance, in descending order
        of variance.
    """
    variances = matrix.var(axis=1, skipna=True)
    top_probes = variances.sort_values(ascending=False).head(n).index
    return matrix.loc[top_probes]
