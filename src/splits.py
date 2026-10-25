"""The two validation schemes, and a guard that proves they are what they claim.

Scheme 1 (random): the number most papers report. Samples from the same hospital
appear in both train and test.
Scheme 2 (leave-one-cohort-out): train on every cohort but one, test on the one
held out. This is the number that answers "will it work somewhere new".
"""
from __future__ import annotations

from typing import Iterator

import numpy as np
import pandas as pd
from sklearn.model_selection import LeaveOneGroupOut, StratifiedKFold

Split = tuple[np.ndarray, np.ndarray]


def random_splits(y: pd.Series, n_splits: int = 5, seed: int = 0) -> Iterator[tuple[str, Split]]:
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    for i, (tr, te) in enumerate(skf.split(np.zeros(len(y)), y.to_numpy()), start=1):
        yield f"fold{i}", (tr, te)


def cohort_splits(cohort: pd.Series) -> Iterator[tuple[str, Split]]:
    logo = LeaveOneGroupOut()
    groups = cohort.to_numpy()
    for tr, te in logo.split(np.zeros(len(cohort)), groups=groups):
        held_out = str(pd.unique(groups[te])[0])
        yield held_out, (tr, te)


def assert_no_cohort_overlap(cohort: pd.Series, split: Split) -> None:
    """Raise if any cohort appears on both sides. Called on every fold at run time."""
    tr, te = split
    groups = cohort.to_numpy()
    both = set(groups[tr]) & set(groups[te])
    if both:
        raise AssertionError(f"cohort leakage: {sorted(both)} in both train and test")


def rank_normalize_within_cohort(X: pd.DataFrame, cohort: pd.Series) -> pd.DataFrame:
    """Rank each gene within its own cohort, then scale to [0, 1].

    This is what published signature work does to survive platform differences.
    It uses no labels, so it does not leak the outcome, but it does assume you
    hold a whole cohort at once rather than one patient at a time. Report results
    with and without it; the difference is itself informative.
    """
    out = X.groupby(cohort, observed=True).transform(
        lambda col: col.rank(pct=True) if col.notna().any() else col
    )
    return out.fillna(0.5)


def rank_normalize_within_sample(X: pd.DataFrame) -> pd.DataFrame:
    """Rank each sample's own genes against each other, scaled to (0, 1].

    The bedside-legal variant: row i is computed from row i alone, so no other
    patient is read, in train or test. Only gene ORDER survives, not absolute
    level. NaN becomes 0.5, the same fill as the cohort version.
    """
    return X.rank(axis=1, pct=True).fillna(0.5)
