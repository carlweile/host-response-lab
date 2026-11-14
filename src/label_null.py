"""Label-shuffle nulls for the raw random-CV headline.

Two nulls through the SAME code path as the headline (splits.random_splits,
train._fit_score, make_models()["logistic_l2"]): shuffle labels over all arrays
(global), or only inside each cohort (within-cohort). Thresholds are the pre-registered,
pre-registered in LABBOOK 2026-11-12 (PLAN commit 81252f2).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import splits, train

LEAK_MEAN = 0.55        # global mean above this = LEAK
LEAK_SINGLE = 0.65      # any single shuffle above this = LEAK
CLEAN_LOW = 0.45        # global mean below this (no LEAK) = ODD-LOW
SHORTCUT_SMALL = 0.60   # within-cohort mean below this = SMALL
SHORTCUT_LARGE = 0.70   # at or above this = LARGE, else MEDIUM
MODEL = "logistic_l2"


def shuffle_global(y: pd.Series, seed: int) -> pd.Series:
    rng = np.random.default_rng(seed)
    return pd.Series(rng.permutation(y.to_numpy()), index=y.index, name=y.name)


def shuffle_within_cohort(y: pd.Series, cohort: pd.Series, seed: int) -> pd.Series:
    """Permute labels only among arrays of the same cohort."""
    rng = np.random.default_rng(seed)
    out = y.to_numpy().copy()
    groups = cohort.to_numpy()
    for c in pd.unique(groups):
        idx = np.flatnonzero(groups == c)
        out[idx] = rng.permutation(out[idx])
    return pd.Series(out, index=y.index, name=y.name)


def check_shuffle(y: pd.Series, y_new: pd.Series, cohort: pd.Series, kind: str) -> int:
    """Guards per run; returns how many labels changed. Raises if the shuffle is wrong."""
    if not y.index.equals(y_new.index):
        raise AssertionError("shuffled labels lost the index order")
    changed = int((y.to_numpy() != y_new.to_numpy()).sum())
    if changed == 0:
        raise AssertionError(f"{kind}: shuffle changed no label")
    if int(y.sum()) != int(y_new.sum()):
        raise AssertionError(f"{kind}: class totals changed")
    if kind == "within_cohort":
        before = y.groupby(cohort.to_numpy()).sum()
        after = y_new.groupby(cohort.to_numpy()).sum()
        if not before.equals(after):
            raise AssertionError("within_cohort: per-cohort class counts changed (shuffle crossed cohorts)")
    return changed


def random_cv_auc(X: pd.DataFrame, y: pd.Series, fold_seed: int = 0, n_splits: int = 5) -> dict:
    """Mean/min/max of the 5 fold AUCs, same definition as summary_raw.csv."""
    model = train.make_models()[MODEL]
    aucs = []
    for _, split in splits.random_splits(y, n_splits=n_splits, seed=fold_seed):
        auc, note = train._fit_score(model, X, y, split)
        if auc is None:
            raise ValueError(f"undefined fold AUC: {note}")
        aucs.append(auc)
    return dict(mean_auc=float(np.mean(aucs)), min_auc=float(np.min(aucs)), max_auc=float(np.max(aucs)))


def loco_auc(X: pd.DataFrame, y: pd.Series, cohort: pd.Series) -> float:
    model = train.make_models()[MODEL]
    aucs = []
    for _, split in splits.cohort_splits(cohort):
        auc, note = train._fit_score(model, X, y, split)
        if auc is None:
            raise ValueError(f"undefined LOCO AUC: {note}")
        aucs.append(auc)
    return float(np.mean(aucs))


def verdict(global_mean: float, global_max_single: float) -> str:
    if global_mean > LEAK_MEAN or global_max_single > LEAK_SINGLE:
        return "LEAK"
    if global_mean < CLEAN_LOW:
        return "ODD-LOW"
    return "CLEAN"


def shortcut_word(within_mean: float) -> str:
    if within_mean < SHORTCUT_SMALL:
        return "SMALL"
    return "LARGE" if within_mean >= SHORTCUT_LARGE else "MEDIUM"


def cohort_only_auc(y: pd.Series, cohort: pd.Series) -> float:
    """Reference: AUC of a score that knows only the cohort's viral share (labels only, no model)."""
    from sklearn.metrics import roc_auc_score
    share = y.groupby(cohort.to_numpy()).transform("mean")
    return float(roc_auc_score(y.to_numpy(), share.to_numpy()))
