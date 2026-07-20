"""MOFA+ multi-omics factor analysis: training input, training, and result extraction.

Wraps ``mofapy2`` (training) and ``mofax`` (reading a trained model) --
confirmed sufficient in pure Python before this module was written; see
``docs/adr/0006-mofa-plus-implementation-and-training.md`` for why this
needed no R dependency, contrary to an earlier assumption in this repo,
and for the verification that mofapy2's long-format input correctly
handles both patients entirely absent from a view and scattered missing
values within a view.
"""

from __future__ import annotations

from pathlib import Path
from typing import Literal, cast

import pandas as pd
from mofapy2.run.entry_point import entry_point
from mofax.core import mofa_model

#: MOFA+ likelihood models this project uses (Gaussian for continuous
#: views, Bernoulli for the binary mutation view).
Likelihood = Literal["gaussian", "bernoulli"]

#: mofapy2 ELBO convergence tolerance settings, tightest to loosest.
ConvergenceMode = Literal["slow", "medium", "fast"]


def build_view_long_df(view: str, matrix: pd.DataFrame) -> pd.DataFrame:
    """Melt one feature x patient matrix into MOFA+'s long-format input.

    Args:
        view: The view's name (e.g. "rna_seq"), recorded in the ``view``
            column.
        matrix: A feature x patient DataFrame (index=feature,
            columns=patient/case_id).

    Returns:
        A DataFrame with columns (sample, feature, view, group, value).
        Missing (``NaN``) values are dropped as rows entirely -- mofapy2
        interprets a row's absence as "not observed", never as zero.
    """
    long_df = matrix.reset_index(names="feature").melt(
        id_vars="feature", var_name="sample", value_name="value"
    )
    long_df["view"] = view
    long_df["group"] = "all"
    return long_df.dropna(subset=["value"])


def build_mofa_input(views: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Combine multiple omic views into one long-format MOFA+ input.

    Args:
        views: Mapping of view name to its feature x patient matrix.

    Returns:
        The concatenation of :func:`build_view_long_df` across all views.
        Patients missing from a given view simply have no rows for that
        view -- MOFA+ is designed to factorize partially-overlapping
        views, so this project does not force a complete-case cohort.
    """
    return pd.concat(
        [build_view_long_df(name, matrix) for name, matrix in views.items()],
        ignore_index=True,
    )


def train_mofa_model(
    views: dict[str, pd.DataFrame],
    likelihoods: dict[str, Likelihood],
    *,
    n_factors: int,
    seed: int,
    outfile: Path,
    convergence_mode: ConvergenceMode = "slow",
    max_iterations: int = 1000,
) -> None:
    """Train a MOFA+ model and save it to an HDF5 file.

    Args:
        views: Mapping of view name to its feature x patient matrix.
        likelihoods: Mapping of view name to its MOFA+ likelihood model.
            Must have exactly the same keys as ``views``.
        n_factors: Number of latent factors to initialise training with.
            ARD priors on feature weights provide sparsity; factors
            explaining negligible variance are excluded during
            interpretation, not training (see ``docs/methods.md`` §3).
        seed: Random seed for reproducible variational initialisation.
            Logged by mofapy2 itself at training time.
        outfile: Path to write the trained model to (HDF5). Parent
            directories are created if needed.
        convergence_mode: mofapy2 ELBO convergence tolerance setting.
        max_iterations: Maximum training iterations.

    Raises:
        ValueError: If ``views`` and ``likelihoods`` do not have exactly
            the same set of view names.
    """
    if set(views) != set(likelihoods):
        raise ValueError(
            f"views and likelihoods must have the same keys: {set(views)} != {set(likelihoods)}"
        )

    data = build_mofa_input(views)
    # mofapy2 sorts view names alphabetically internally
    # (np.sort(data["view"].unique())) regardless of input order, so the
    # likelihoods list must follow that same order -- passing likelihoods
    # in dict-insertion order silently misassigns them to the wrong view.
    view_names = sorted(views)

    ep = entry_point()
    ep.set_data_options(scale_views=True)
    ep.set_data_df(data, likelihoods=[likelihoods[name] for name in view_names])
    ep.set_model_options(factors=n_factors)
    ep.set_train_options(
        iter=max_iterations,
        convergence_mode=convergence_mode,
        seed=seed,
        quiet=False,
        verbose=False,
    )
    ep.build()
    ep.run()

    outfile.parent.mkdir(parents=True, exist_ok=True)
    ep.save(str(outfile))


def load_mofa_model(outfile: Path) -> mofa_model:
    """Load a trained MOFA+ model for interpretation.

    Args:
        outfile: Path to a trained model HDF5 file, e.g. from
            :func:`train_mofa_model`.

    Returns:
        A ``mofax.mofa_model`` providing factor value, feature weight,
        and variance-explained accessors.
    """
    return mofa_model(str(outfile))


def get_factor_values(model: mofa_model) -> pd.DataFrame:
    """Extract per-patient factor values from a trained model.

    Args:
        model: A loaded :func:`load_mofa_model` result.

    Returns:
        A DataFrame indexed by sample (patient), one column per factor.
    """
    return cast(pd.DataFrame, model.get_factors(df=True))


def get_variance_explained(model: mofa_model) -> pd.DataFrame:
    """Extract the fraction of variance each factor explains in each view.

    Args:
        model: A loaded :func:`load_mofa_model` result.

    Returns:
        A long-format DataFrame with one row per (factor, view) pair and
        columns ``Factor``, ``View``, ``Group``, ``R2`` (the fraction of
        that view's variance the factor explains).
    """
    return cast(pd.DataFrame, model.get_variance_explained())
