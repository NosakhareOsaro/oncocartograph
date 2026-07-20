"""Univariate survival association testing (Cox proportional hazards).

Provides one uniform statistical treatment for a candidate biomarker's
association with overall survival, whether the biomarker's value is
continuous (e.g. gene expression, methylation M-value, relative copy
number) or binary (e.g. mutation presence/absence). Cox PH on overall
survival is used rather than Fine-Gray competing-risks regression
because the TCGA-BRCA clinical file has no cause-of-death or recurrence
coding available to define a competing event against -- this is a data
availability constraint, not a stylistic choice (see docs/methods.md).

This module has zero dependency on the rest of oncocartograph: it
operates purely on plain pandas Series/DataFrames, so it can be
extracted to a standalone package without modification.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.exceptions import ConvergenceError
from numpy.linalg import LinAlgError
from scipy.stats import false_discovery_control


@dataclass(frozen=True)
class SurvivalEvidence:
    """Univariate Cox PH association between one biomarker and overall survival.

    Attributes:
        hazard_ratio: exp(coefficient). >1 means higher biomarker values
            (or presence, for binary biomarkers) associate with worse
            survival; <1 means better survival.
        hazard_ratio_ci_low: Lower bound of the 95% CI for hazard_ratio.
        hazard_ratio_ci_high: Upper bound of the 95% CI for hazard_ratio.
        p_value: Unadjusted Wald test p-value for the biomarker's
            coefficient.
        p_value_adjusted: Benjamini-Hochberg FDR-corrected p-value,
            filled in once this evidence is placed in the context of a
            batch of simultaneously-tested candidates (see
            :func:`screen_survival_associations`). ``None`` for a single
            ad hoc test with no such batch context -- callers that build
            a composite score across many candidates must supply this
            themselves rather than use the raw ``p_value``.
        n_samples: Number of patients used to fit this model (after
            dropping missing values).
        n_events: Number of observed events (deaths) among those patients.
    """

    hazard_ratio: float
    hazard_ratio_ci_low: float
    hazard_ratio_ci_high: float
    p_value: float
    n_samples: int
    n_events: int
    p_value_adjusted: float | None = None


def fit_univariate_cox(
    values: pd.Series, duration: pd.Series, event: pd.Series
) -> SurvivalEvidence | None:
    """Fit a univariate Cox PH model for one biomarker against overall survival.

    Args:
        values: Biomarker values indexed by patient (continuous, or
            binary 0/1 for e.g. mutation presence). Patients missing a
            value are dropped.
        duration: Time-to-event or time-to-censoring, indexed by patient.
        event: 1 if the event (death) was observed, 0 if censored,
            indexed by patient.

    Returns:
        A :class:`SurvivalEvidence`, or ``None`` if the model could not
        be fit -- e.g. the biomarker has zero variance among patients
        with complete data, which happens for some genes/probes in
        small subsets and must not crash a batch screening run.
    """
    frame = pd.DataFrame({"value": values, "duration": duration, "event": event}).dropna()
    if frame["value"].nunique() < 2 or len(frame) < 2:
        return None

    cph = CoxPHFitter()
    try:
        cph.fit(frame, duration_col="duration", event_col="event")
    except (ConvergenceError, ValueError, LinAlgError):
        return None

    summary = cph.summary.loc["value"]
    key_stats = summary[["exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]].to_numpy(
        dtype=float
    )
    if not np.all(np.isfinite(key_stats)):
        # lifelines does not always raise on degenerate fits -- a very rare
        # or very sparse binary covariate (e.g. mutated in only 1-4 of 122
        # patients, with zero events in that subgroup) can converge to a
        # trivial, uninformative solution instead of raising: sometimes NaN
        # summary statistics, sometimes a finite-but-meaningless estimate
        # with an infinite confidence bound (HR~0, CI upper=inf). Both were
        # confirmed on real TCGA-BRCA mutation data (~11% of
        # recurrence-filtered genes) and both must be excluded here, not
        # passed through to FDR correction (which raises on out-of-range
        # p-values and would be misled by an unbounded CI regardless).
        return None

    return SurvivalEvidence(
        hazard_ratio=float(summary["exp(coef)"]),
        hazard_ratio_ci_low=float(summary["exp(coef) lower 95%"]),
        hazard_ratio_ci_high=float(summary["exp(coef) upper 95%"]),
        p_value=float(summary["p"]),
        n_samples=len(frame),
        n_events=int(frame["event"].sum()),
    )


def screen_survival_associations(
    matrix: pd.DataFrame, duration: pd.Series, event: pd.Series
) -> pd.DataFrame:
    """Screen many candidate biomarkers for survival association at once.

    Args:
        matrix: A feature x patient DataFrame (index=feature/gene,
            columns=patient), continuous or binary values.
        duration: Time-to-event or time-to-censoring, indexed by patient.
        event: 1 if the event was observed, 0 if censored, indexed by
            patient.

    Returns:
        A DataFrame indexed by feature, with columns ``hazard_ratio``,
        ``hazard_ratio_ci_low``, ``hazard_ratio_ci_high``, ``p_value``,
        ``n_samples``, ``n_events``, and ``p_adj`` (Benjamini-Hochberg
        FDR-corrected p-value across every feature screened in this
        call). Features whose model could not be fit are omitted
        entirely, not included with null values.
    """
    rows: dict[str, dict[str, float | int]] = {}
    for feature, row in matrix.iterrows():
        evidence = fit_univariate_cox(row, duration, event)
        if evidence is not None:
            rows[str(feature)] = {
                "hazard_ratio": evidence.hazard_ratio,
                "hazard_ratio_ci_low": evidence.hazard_ratio_ci_low,
                "hazard_ratio_ci_high": evidence.hazard_ratio_ci_high,
                "p_value": evidence.p_value,
                "n_samples": evidence.n_samples,
                "n_events": evidence.n_events,
            }

    table = pd.DataFrame.from_dict(rows, orient="index")
    if not table.empty:
        table["p_adj"] = false_discovery_control(table["p_value"], method="bh")
    return table
