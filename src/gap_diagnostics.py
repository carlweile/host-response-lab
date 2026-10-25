"""Helpers that take the negative logistic gap apart (LABBOOK 2026-10-19 PLAN).

`pair_auc` splits one AUC into same-cohort and cross-cohort pairs, so a random-CV
number that mixes hospitals can be compared with a held-out number that does not.
`balanced_subsample` makes every cohort 50 % bacterial, or (size control) keeps
each cohort's own bacterial share at the same size. `rank_against_reference`
ranks new patients against a fixed reference only, so no other new patient is
ever read; `cohort_rank_inductive` uses it to fit cohort-wise ranking on the
training rows of a fold.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class PairAUC:
    pooled: float      # all (bacterial, viral) pairs: the number on file
    within: float      # pairs whose two patients share a cohort
    cross: float       # pairs whose two patients come from different cohorts
    n_within: int
    n_cross: int


def pair_auc(y, score, group) -> PairAUC:
    """AUC over all (bacterial, viral) pairs, and over same-group / cross-group pairs.

    AUC is the share of pairs in which the bacterial patient scores higher (ties
    count half), so `pooled` is exactly the pair-weighted mix of `within` and
    `cross`. Adding a constant to every score of one group cannot change `within`.
    """
    y = np.asarray(y).astype(int)
    s = np.asarray(score, dtype=float)
    g = np.asarray(group)
    pos, neg = y == 1, y == 0
    if not pos.any() or not neg.any():
        raise ValueError("pair_auc needs both classes")
    sp, sn = s[pos][:, None], s[neg][None, :]
    wins = (sp > sn) + 0.5 * (sp == sn)
    same = g[pos][:, None] == g[neg][None, :]

    def share(mask: np.ndarray) -> tuple[float, int]:
        n = int(mask.sum())
        return (float(wins[mask].sum() / n) if n else float("nan")), n

    pooled, _ = share(np.ones_like(same))
    within, n_within = share(same)
    cross, n_cross = share(~same)
    return PairAUC(pooled, within, cross, n_within, n_cross)


def balanced_subsample(
    y: pd.Series, cohort: pd.Series, seed: int, keep_share: bool = False
) -> np.ndarray:
    """Sorted positional indices of a random subset, 2 x (smaller class) per cohort.

    keep_share=False: half bacterial, half viral in every cohort.
    keep_share=True (the size control): the same number of patients per cohort,
    but the cohort's own bacterial share, rounded.
    Labels are read on purpose: this builds a diagnostic dataset, it is not a
    model step, and the model never sees how the subset was chosen.
    """
    rng = np.random.default_rng(seed)
    yv, cv = y.to_numpy().astype(int), cohort.to_numpy()
    keep = []
    for g in pd.unique(cv):
        idx_b = np.flatnonzero((cv == g) & (yv == 1))
        idx_v = np.flatnonzero((cv == g) & (yv == 0))
        size = 2 * min(len(idx_b), len(idx_v))
        if size == 0:
            raise ValueError(f"cohort {g} lacks one class")
        if keep_share:
            n_b = int(round(size * len(idx_b) / (len(idx_b) + len(idx_v))))
        else:
            n_b = size // 2
        keep.append(rng.choice(idx_b, n_b, replace=False))
        keep.append(rng.choice(idx_v, size - n_b, replace=False))
    return np.sort(np.concatenate(keep))


def rank_against_reference(ref: np.ndarray, new: np.ndarray) -> np.ndarray:
    """Percentile of each new value among the reference values of the same gene.

    Row i of the result is computed from `ref` and row i of `new` only: its rank
    in (reference + itself), ties averaged, divided by len(ref) + 1. That equals
    pandas rank(pct=True) of the row after it joined the reference, without
    reading any other new row.
    """
    ref = np.asarray(ref, dtype=float)
    new = np.asarray(new, dtype=float)
    if np.isnan(ref).any() or np.isnan(new).any():
        raise ValueError("rank_against_reference expects no NaN (build_matrix imputes)")
    out = np.empty_like(new)
    denom = ref.shape[0] + 1
    for i, row in enumerate(new):
        less = (ref < row).sum(axis=0)
        ties = (ref == row).sum(axis=0)
        out[i] = (1.0 + less + 0.5 * ties) / denom
    return out


def cohort_rank_inductive(
    X: pd.DataFrame, cohort: pd.Series, train: np.ndarray, test: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Cohort-wise ranking fitted on the training rows of one fold only.

    Training rows are ranked within the training rows of their own cohort (what
    `rank_normalize_within_cohort` does, minus the test rows). Each test row is
    ranked against the training rows of its own cohort, one row at a time.
    """
    Xv = X.to_numpy(dtype=float)
    cv = cohort.to_numpy()
    missing = sorted(set(cv[test]) - set(cv[train]))
    if missing:
        raise ValueError(f"test cohorts with no training rows: {missing}")
    X_tr = np.empty((len(train), Xv.shape[1]))
    X_te = np.empty((len(test), Xv.shape[1]))
    for g in pd.unique(cv[train]):
        tr_g = np.flatnonzero(cv[train] == g)
        ref = Xv[train[tr_g]]
        X_tr[tr_g] = pd.DataFrame(ref).rank(pct=True).to_numpy()
        te_g = np.flatnonzero(cv[test] == g)
        if len(te_g):
            X_te[te_g] = rank_against_reference(ref, Xv[test[te_g]])
    return X_tr, X_te
