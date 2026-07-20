"""Generic categorical association testing (Fisher's exact test).

Not specific to any particular categorical variable -- usable to test a
binary biomarker (e.g. mutation presence) against any two-category
outcome, such as ``vital_status`` today, or a molecular subtype label
once one is defined elsewhere in the pipeline (``feat/validation``).
Zero dependency on the rest of oncocartograph.
"""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd
from scipy.stats import fisher_exact


@dataclass(frozen=True)
class AssociationEvidence:
    """Fisher's exact test result for a binary marker vs. a categorical outcome.

    Attributes:
        odds_ratio: Sample odds ratio from the 2x2 contingency table.
        p_value: Two-sided Fisher's exact test p-value.
        contingency_table: The 2x2 table actually used
            ((marker_neg/outcome_a, marker_neg/outcome_b),
            (marker_pos/outcome_a, marker_pos/outcome_b)), retained for
            auditability.
        n_samples: Total patients used, after dropping missing values.
    """

    odds_ratio: float
    p_value: float
    contingency_table: tuple[tuple[int, int], tuple[int, int]]
    n_samples: int


def fisher_exact_association(marker: pd.Series, outcome: pd.Series) -> AssociationEvidence | None:
    """Test association between a binary marker and a binary categorical outcome.

    Args:
        marker: Binary (0/1 or boolean) values indexed by patient.
        outcome: A two-category variable indexed by patient (e.g.
            ``vital_status``, or a subtype label) -- must have exactly
            two distinct non-null values.

    Returns:
        An :class:`AssociationEvidence`, or ``None`` if the data doesn't
        support the test (fewer than two categories present in either
        variable after dropping missing values).
    """
    frame = pd.DataFrame({"marker": marker, "outcome": outcome}).dropna()
    if frame["marker"].nunique() < 2 or frame["outcome"].nunique() != 2:
        return None

    marker_categories = sorted(frame["marker"].unique(), key=str)
    outcome_categories = sorted(frame["outcome"].unique(), key=str)
    table = pd.crosstab(frame["marker"], frame["outcome"]).reindex(
        index=marker_categories, columns=outcome_categories, fill_value=0
    )
    counts = table.to_numpy()
    odds_ratio, p_value = fisher_exact(counts)
    contingency_table = (
        (int(counts[0, 0]), int(counts[0, 1])),
        (int(counts[1, 0]), int(counts[1, 1])),
    )
    return AssociationEvidence(
        odds_ratio=float(odds_ratio),
        p_value=float(p_value),
        contingency_table=contingency_table,
        n_samples=len(frame),
    )
