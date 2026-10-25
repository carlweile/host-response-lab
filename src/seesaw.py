"""One fold runner for the pre-registered seesaw test (LABBOOK PLAN 9daaf08).

An arm is a per-cohort transform. WHOLE arms are fitted on the whole cohort
before the split, so a cohort's test rows are tied to its training rows. TRAIN
arms are their twins: fitted on the fold's training rows of each cohort only,
so no test row is ever read. Every whole-cohort arm is a transductive
diagnostic of a random-CV artefact, not a generalisation number.

Per fold the runner keeps the same-cohort / cross-cohort pair AUC (pairs only
inside that fold, one model), the class-adjusted cohort offsets of the test and
training scores, and the two numbers the P3 gates need.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .gap_diagnostics import cohort_rank_inductive, pair_auc
from .gap_experiments import _fit, _logit, _n_iter
from .splits import random_splits, rank_normalize_within_cohort
from .train import make_models

LGBM_THREADS = 8
WHOLE = ("W", "raw", "center_whole", "scale_whole", "zscore_whole")
TRAIN = ("I", "center_train", "scale_train", "zscore_train")
AFFINE = {"center": (True, False), "scale": (False, True), "zscore": (True, True)}


def c_for_k(k: int) -> float:
    """Same L2 penalty per training patient as the file value C=0.05 at k=5."""
    return 0.04 * k / (k - 1)


def logistic(C: float):
    return make_models()["logistic_l2"].set_params(clf__C=C)


def lgbm(stumps: bool = False):
    model = make_models()["lightgbm"].set_params(n_jobs=LGBM_THREADS)
    return model.set_params(num_leaves=2, max_depth=1) if stumps else model


def _stats(Xv: np.ndarray, cv: np.ndarray, rows: np.ndarray) -> dict:
    out = {}
    for g in pd.unique(cv[rows]):
        block = Xv[rows[cv[rows] == g]]
        sd = block.std(axis=0)
        out[g] = (block.mean(axis=0), np.where(sd > 0, sd, 1.0))
    return out


def _affine(Xv, cv, rows, stats, center: bool, scale: bool) -> np.ndarray:
    out = Xv[rows].copy()
    for g, (mu, sd) in stats.items():
        m = cv[rows] == g
        if center:
            out[m] -= mu
        if scale:
            out[m] /= sd
    return out


def whole_matrix(arm: str, X: pd.DataFrame, cohort: pd.Series) -> np.ndarray:
    """A WHOLE arm, computed once per cohort set before any split."""
    if arm == "W":
        return rank_normalize_within_cohort(X, cohort).to_numpy(dtype=float)
    Xv, cv = X.to_numpy(dtype=float), cohort.to_numpy()
    if arm == "raw":
        return Xv
    every = np.arange(len(cv))
    return _affine(Xv, cv, every, _stats(Xv, cv, every), *AFFINE[arm.split("_")[0]])


def train_matrices(arm, X, cohort, tr, te) -> tuple[np.ndarray, np.ndarray]:
    """A TRAIN arm: statistics from the training rows of each cohort only."""
    if arm == "I":
        return cohort_rank_inductive(X, cohort, tr, te)
    Xv, cv = X.to_numpy(dtype=float), cohort.to_numpy()
    stats = _stats(Xv, cv, tr)
    missing = sorted(set(cv[te]) - set(stats))
    if missing:
        raise ValueError(f"test cohorts with no training rows: {missing}")
    flags = AFFINE[arm.split("_")[0]]
    return _affine(Xv, cv, tr, stats, *flags), _affine(Xv, cv, te, stats, *flags)


def treatment_gate(X: pd.DataFrame, cohort: pd.Series) -> dict:
    """P2 gate, before any fit: median over genes of the SD across cohorts of the
    cohort means, raw vs after center_whole vs after scale_whole."""
    def spread(M: np.ndarray) -> float:
        means = pd.DataFrame(M).groupby(cohort.to_numpy()).mean()
        return float(np.median(means.std(axis=0, ddof=0)))
    raw = spread(X.to_numpy(dtype=float))
    center = spread(whole_matrix("center_whole", X, cohort))
    scale = spread(whole_matrix("scale_whole", X, cohort))
    return dict(spread_raw=raw, spread_center_whole=center, spread_scale_whole=scale,
                center_ok=center <= 1e-9, scale_ok=scale >= 0.1 * raw)


def class_offset(s: np.ndarray, y: np.ndarray, c: np.ndarray, g) -> float:
    """Mean score of cohort g minus the mean over all rows, per class, averaged
    over the classes cohort g has in these rows (NaN if it has none)."""
    vals = [s[(c == g) & (y == k)].mean() - s[y == k].mean()
            for k in (0, 1) if ((c == g) & (y == k)).any()]
    return float(np.mean(vals)) if vals else float("nan")


def pooled_within_sd(s: np.ndarray, y: np.ndarray) -> float:
    resid = np.concatenate([s[y == k] - s[y == k].mean() for k in (0, 1)])
    return float(np.sqrt((resid ** 2).sum() / (len(resid) - 2)))


def gate_numbers(s_tr, s_te, y_tr, y_te, c_tr, c_te) -> dict:
    """algebra: SD across cohorts of the (train + test) mean score / within-class SD.
    pin: slope of the per-cohort training-mean score on the cohort's training
    bacterial share / within-class SD of the training scores."""
    s, y, c = (np.concatenate(p) for p in ((s_tr, s_te), (y_tr, y_te), (c_tr, c_te)))
    groups = pd.unique(c)
    totals = np.array([s[c == g].mean() for g in groups])
    shares = np.array([y_tr[c_tr == g].mean() for g in pd.unique(c_tr)])
    tr_means = np.array([s_tr[c_tr == g].mean() for g in pd.unique(c_tr)])
    slope = np.polyfit(shares, tr_means, 1)[0] if len(shares) > 1 else float("nan")
    return dict(algebra=float(totals.std() / pooled_within_sd(s, y)),
                pin=float(slope / pooled_within_sd(s_tr, y_tr)))


def run_arm(X, y, cohort, arm: str, model, k: int, seed: int, tag: dict,
            only_first_fold: bool = False) -> tuple[list[dict], list[dict]]:
    """Fit `model` on every fold of one stratified k-fold split; return fold rows
    and per-cohort offset rows. `tag` (set, part, model name, C) is copied in."""
    yv, cv = y.to_numpy().astype(int), cohort.to_numpy()
    Xw = whole_matrix(arm, X, cohort) if arm in WHOLE else None
    pooled_share = float(yv.mean())
    folds, offsets = [], []
    for fold, (tr, te) in random_splits(y, n_splits=k, seed=seed):
        X_tr, X_te = (Xw[tr], Xw[te]) if Xw is not None else train_matrices(arm, X, cohort,
                                                                              tr, te)
        est = _fit(model, X_tr, yv[tr])
        pa = pair_auc(yv[te], est.predict_proba(X_te)[:, 1], cv[te])
        s_tr, s_te = _logit(est, X_tr), _logit(est, X_te)
        key = dict(tag, arm=arm, k=k, seed=seed, fold=fold)
        folds.append(dict(key, n_train=len(tr), n_test=len(te), n_pos_test=int(yv[te].sum()),
                          auc_within=pa.within, auc_cross=pa.cross,
                          n_pairs_within=pa.n_within, n_pairs_cross=pa.n_cross,
                          n_iter=_n_iter(est),
                          **gate_numbers(s_tr, s_te, yv[tr], yv[te], cv[tr], cv[te])))
        for g in pd.unique(cv):
            offsets.append(dict(key, cohort=g, share=float(yv[cv == g].mean()),
                                pooled_share=pooled_share,
                                o_test=class_offset(s_te, yv[te], cv[te], g),
                                o_train=class_offset(s_tr, yv[tr], cv[tr], g)))
        if only_first_fold:
            break
    return folds, offsets
