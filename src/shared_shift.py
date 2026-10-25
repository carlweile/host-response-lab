"""Runner for the pre-registered shared-shift test (LABBOOK PLAN 2fd8865).

Two arms are new against `src/seesaw.py`:
- "I_bal": each patient is ranked against its cohort's TRAINING rows,
  re-weighted to equal bacterial and viral mass. Training labels only; no test
  row is ever read.
- "W_recentre": the whole-cohort-rank model's test scores of cohort g are moved
  by minus the class-adjusted offset of cohort g measured on TRAINING scores.

Every whole-cohort arm is a transductive diagnostic of a random-CV artefact.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .gap_diagnostics import cohort_rank_inductive, pair_auc
from .gap_experiments import _fit, _logit, _n_iter
from .seesaw import class_offset, gate_numbers, logistic, whole_matrix
from .splits import random_splits


def balanced_weights(y_ref: np.ndarray) -> np.ndarray:
    """Weights that give each class half the mass (one class only: equal weights)."""
    y_ref = np.asarray(y_ref).astype(int)
    n_b, n_v = int(y_ref.sum()), int((1 - y_ref).sum())
    if n_b == 0 or n_v == 0:
        return np.full(len(y_ref), 1.0 / len(y_ref))
    return np.where(y_ref == 1, 0.5 / n_b, 0.5 / n_v)


def rank_weighted(ref: np.ndarray, w: np.ndarray, new: np.ndarray) -> np.ndarray:
    """Weighted percentile of each new row among `ref`: mass below + half the tied mass.
    Row i of the result reads only `ref`, `w` and row i of `new`."""
    out = np.empty_like(new, dtype=float)
    for i, row in enumerate(new):
        out[i] = (w[:, None] * (ref < row)).sum(0) + 0.5 * (w[:, None] * (ref == row)).sum(0)
    return out


def rank_balanced_inductive(X, y, cohort, tr, te) -> tuple[np.ndarray, np.ndarray]:
    Xv, yv, cv = X.to_numpy(dtype=float), y.to_numpy().astype(int), cohort.to_numpy()
    missing = sorted(set(cv[te]) - set(cv[tr]))
    if missing:
        raise ValueError(f"test cohorts with no training rows: {missing}")
    X_tr, X_te = np.empty((len(tr), Xv.shape[1])), np.empty((len(te), Xv.shape[1]))
    for g in pd.unique(cv[tr]):
        tr_g, te_g = np.flatnonzero(cv[tr] == g), np.flatnonzero(cv[te] == g)
        ref, w = Xv[tr[tr_g]], balanced_weights(yv[tr[tr_g]])
        X_tr[tr_g] = rank_weighted(ref, w, ref)
        if len(te_g):
            X_te[te_g] = rank_weighted(ref, w, Xv[te[te_g]])
    return X_tr, X_te


def recentred(s_tr, y_tr, c_tr, s_te, c_te) -> np.ndarray:
    """Test scores minus the TRAINING class-adjusted offset of their cohort."""
    shift = {g: class_offset(s_tr, y_tr, c_tr, g) for g in pd.unique(c_tr)}
    return s_te - np.array([shift[g] for g in c_te])


def _row(key, arm, pa, est, yv, tr, te, extra=None):
    return dict(key, arm=arm, n_test=len(te), n_pos_test=int(yv[te].sum()),
                auc_within=pa.within, auc_cross=pa.cross, n_pairs_within=pa.n_within,
                n_pairs_cross=pa.n_cross, n_iter=_n_iter(est), **(extra or {}))


def run_k5(X, y, cohort, seed: int, tag: dict) -> tuple[list[dict], list[dict]]:
    """k=5, C=0.05: arms W, W_recentre, I, I_bal, raw on the same folds."""
    yv, cv = y.to_numpy().astype(int), cohort.to_numpy()
    model = logistic(0.05)
    Xw, Xr = whole_matrix("W", X, cohort), whole_matrix("raw", X, cohort)
    pooled = float(yv.mean())
    folds, offsets = [], []
    for fold, (tr, te) in random_splits(y, n_splits=5, seed=seed):
        key = dict(tag, k=5, C=0.05, seed=seed, fold=fold)
        mats = {"W": (Xw[tr], Xw[te]), "raw": (Xr[tr], Xr[te]),
                "I": cohort_rank_inductive(X, cohort, tr, te),
                "I_bal": rank_balanced_inductive(X, y, cohort, tr, te)}
        for arm, (X_tr, X_te) in mats.items():
            est = _fit(model, X_tr, yv[tr])
            s_tr, s_te = _logit(est, X_tr), _logit(est, X_te)
            gates = gate_numbers(s_tr, s_te, yv[tr], yv[te], cv[tr], cv[te])
            folds.append(_row(key, arm, pair_auc(yv[te], s_te, cv[te]), est, yv, tr, te, gates))
            if arm == "W":
                s_rc = recentred(s_tr, yv[tr], cv[tr], s_te, cv[te])
                folds.append(_row(key, "W_recentre", pair_auc(yv[te], s_rc, cv[te]), est,
                                  yv, tr, te, dict(algebra=np.nan, pin=np.nan)))
            for g in pd.unique(cv):
                offsets.append(dict(key, arm=arm, cohort=g, share=float(yv[cv == g].mean()),
                                    pooled_share=pooled,
                                    o_test=class_offset(s_te, yv[te], cv[te], g),
                                    o_train=class_offset(s_tr, yv[tr], cv[tr], g)))
    return folds, offsets
