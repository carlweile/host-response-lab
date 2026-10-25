"""Experiments E0-E4 behind LABBOOK 2026-10-19 (the negative logistic gap).

Each function takes the raw matrix (X, y, cohort from build_matrix) and returns
tidy rows. Models and splits come unchanged from src.train / src.splits, so the
numbers are comparable with results/folds_*.csv.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.metrics import roc_auc_score

from .gap_diagnostics import (balanced_subsample, cohort_rank_inductive, pair_auc,
                              rank_against_reference)
from .splits import (assert_no_cohort_overlap, cohort_splits, random_splits,
                     rank_normalize_within_cohort, rank_normalize_within_sample)
from .train import make_models

TRANSFORMS = {
    "raw": lambda X, c: X,
    "cohort_rank": rank_normalize_within_cohort,
    "within_sample": lambda X, c: rank_normalize_within_sample(X),
}
NAN = float("nan")


def _fit(model, X: np.ndarray, y: np.ndarray):
    est = clone(model)
    est.fit(X, y)
    return est


def _logit(est, X: np.ndarray) -> np.ndarray:
    """Score on the log-odds scale (no clipping): logistic margin or LightGBM raw score."""
    if hasattr(est, "decision_function"):
        return est.decision_function(X)
    return est.predict(X, raw_score=True)


def _n_iter(est) -> float:
    clf = getattr(est, "named_steps", {}).get("clf")
    return float(np.max(clf.n_iter_)) if clf is not None and hasattr(clf, "n_iter_") else NAN


def _auc(y: np.ndarray, s: np.ndarray) -> float:
    return float(roc_auc_score(y, s)) if len(set(y)) == 2 else NAN


def both_schemes_with_scores(X, y, cohort, transform: str, models=None, n_splits: int = 5):
    """E0 + E1: rerun both schemes, keep every test score, split the random-CV AUC by pair type.

    Returns (fold rows, per-cohort rows). The per-cohort rows put the random-CV AUC
    measured inside one cohort (out-of-fold scores) next to that cohort's held-out AUC.
    """
    Xt = TRANSFORMS[transform](X, cohort)
    if not (Xt.index.equals(y.index) and Xt.index.equals(cohort.index)):
        raise AssertionError("rows of X, y and cohort are not aligned")
    Xv, yv, cv = Xt.to_numpy(dtype=float), y.to_numpy().astype(int), cohort.to_numpy()
    folds, cohorts = [], []
    for name, model in (models or make_models()).items():
        oof_p, oof_logit = np.full(len(yv), NAN), np.full(len(yv), NAN)
        for fold, (tr, te) in random_splits(y, n_splits=n_splits):
            est = _fit(model, Xv[tr], yv[tr])
            p = est.predict_proba(Xv[te])[:, 1]
            oof_p[te], oof_logit[te] = p, _logit(est, Xv[te])
            pa = pair_auc(yv[te], p, cv[te])
            if abs(pa.pooled - roc_auc_score(yv[te], p)) > 1e-9:
                raise AssertionError("pair_auc pooled disagrees with sklearn")
            folds.append(dict(transform=transform, model=name, scheme="random_cv", fold=fold,
                              n_test=len(te), n_pos_test=int(yv[te].sum()),
                              auc_pooled=pa.pooled, auc_within=pa.within, auc_cross=pa.cross,
                              n_pairs_within=pa.n_within, n_pairs_cross=pa.n_cross,
                              n_iter=_n_iter(est)))
        loco = {}
        for held_out, (tr, te) in cohort_splits(cohort):
            assert_no_cohort_overlap(cohort, (tr, te))
            est = _fit(model, Xv[tr], yv[tr])
            auc = _auc(yv[te], est.predict_proba(Xv[te])[:, 1])
            loco[held_out] = auc
            n_pos = int(yv[te].sum())
            folds.append(dict(transform=transform, model=name, scheme="leave_one_cohort_out",
                              fold=held_out, n_test=len(te), n_pos_test=n_pos,
                              auc_pooled=auc, auc_within=auc, auc_cross=NAN,
                              n_pairs_within=n_pos * (len(te) - n_pos), n_pairs_cross=0,
                              n_iter=_n_iter(est)))
        for g in pd.unique(cv):
            m = cv == g
            b, v = m & (yv == 1), m & (yv == 0)
            resid = np.concatenate([oof_logit[b] - oof_logit[b].mean(),
                                    oof_logit[v] - oof_logit[v].mean()])
            cohorts.append(dict(transform=transform, model=name, cohort=g, n=int(m.sum()),
                                bacterial_share=float(yv[m].mean()),
                                rcv_auc_in_cohort=_auc(yv[m], oof_p[m]), loco_auc=loco[g],
                                mean_logit_bacterial=float(oof_logit[b].mean()),
                                mean_logit_viral=float(oof_logit[v].mean()),
                                sd_logit_within_class=float(np.sqrt((resid ** 2).sum()
                                                                    / (len(resid) - 2)))))
    return pd.DataFrame(folds), pd.DataFrame(cohorts)


def inductive_random_cv(X, y, cohort, models=None, n_splits: int = 5) -> pd.DataFrame:
    """E3: cohort-wise ranking fitted per fold; test rows ranked against training rows only."""
    yv, cv = y.to_numpy().astype(int), cohort.to_numpy()
    rows = []
    for name, model in (models or make_models()).items():
        for fold, (tr, te) in random_splits(y, n_splits=n_splits):
            X_tr, X_te = cohort_rank_inductive(X, cohort, tr, te)
            est = _fit(model, X_tr, yv[tr])
            pa = pair_auc(yv[te], est.predict_proba(X_te)[:, 1], cv[te])
            rows.append(dict(transform="cohort_rank_inductive", model=name, scheme="random_cv",
                             fold=fold, n_test=len(te), n_pos_test=int(yv[te].sum()),
                             auc_pooled=pa.pooled, auc_within=pa.within, auc_cross=pa.cross,
                             n_pairs_within=pa.n_within, n_pairs_cross=pa.n_cross,
                             n_iter=_n_iter(est)))
    return pd.DataFrame(rows)


def split_half_heldout(X, y, cohort, models=None, n_halves: int = 20) -> pd.DataFrame:
    """E4: score half B of a held-out cohort two ways and compare the AUC on B.

    whole: B ranked within the whole held-out cohort (the version on file, it reads
    every other B patient). half_a: B ranked against half A only, one patient at a time.
    """
    Xv, yv = X.to_numpy(dtype=float), y.to_numpy().astype(int)
    X_cr = rank_normalize_within_cohort(X, cohort).to_numpy(dtype=float)
    rows = []
    for name, model in (models or make_models()).items():
        for held_out, (tr, te) in cohort_splits(cohort):
            est = _fit(model, X_cr[tr], yv[tr])
            p_full = est.predict_proba(X_cr[te])[:, 1]
            auc_full = _auc(yv[te], p_full)
            half = len(te) // 2
            for k in range(n_halves):
                perm = np.random.default_rng(k).permutation(len(te))
                a, b = te[perm[:half]], te[perm[half:]]
                p_half = est.predict_proba(rank_against_reference(Xv[a], Xv[b]))[:, 1]
                rows.append(dict(model=name, cohort=held_out, split=k, n_a=len(a), n_b=len(b),
                                 n_pos_b=int(yv[b].sum()), auc_full_cohort=auc_full,
                                 auc_b_whole=_auc(yv[b], p_full[perm[half:]]),
                                 auc_b_half_a=_auc(yv[b], p_half)))
    df = pd.DataFrame(rows)
    df["diff_whole_minus_half_a"] = df["auc_b_whole"] - df["auc_b_half_a"]
    return df


def score_mirror(X, y, cohort, models=None, n_splits: int = 5) -> pd.DataFrame:
    """Exploratory, added after E1/E3 (NOT pre-registered): per random-CV fold model,
    mean log-odds per cohort for its TRAIN and its TEST patients.

    ranking='whole_cohort' is the version on file (each cohort ranked before the
    split); 'train_only' is E3. With whole-cohort ranks, a cohort's mean features
    are almost fixed (a percentile rank averages (n+1)/2n), so for a linear model
    train + test scores of a cohort sum to an almost fixed total (exactly fixed only
    for equal-sized cohorts): pushing the training patients one way pushes the test
    patients the other.
    """
    X_cr = rank_normalize_within_cohort(X, cohort).to_numpy(dtype=float)
    yv, cv = y.to_numpy().astype(int), cohort.to_numpy()
    rows = []
    for name, model in (models or make_models()).items():
        for fold, (tr, te) in random_splits(y, n_splits=n_splits):
            for ranking in ("whole_cohort", "train_only"):
                if ranking == "whole_cohort":
                    X_tr, X_te = X_cr[tr], X_cr[te]
                else:
                    X_tr, X_te = cohort_rank_inductive(X, cohort, tr, te)
                est = _fit(model, X_tr, yv[tr])
                for part, idx, Xp in (("train", tr, X_tr), ("test", te, X_te)):
                    z, yp, cp = _logit(est, Xp), yv[idx], cv[idx]
                    for g in pd.unique(cv):
                        m = cp == g
                        b, v = m & (yp == 1), m & (yp == 0)
                        rows.append(dict(model=name, ranking=ranking, fold=fold, part=part,
                                         cohort=g, n=int(m.sum()), n_bacterial=int(b.sum()),
                                         cohort_bacterial_share=float(yv[cv == g].mean()),
                                         mean_logit=float(z[m].mean()),
                                         mean_logit_bacterial=float(z[b].mean()) if b.any() else NAN,
                                         mean_logit_viral=float(z[v].mean()) if v.any() else NAN))
    return pd.DataFrame(rows)


def balance_control(X, y, cohort, design: str, seeds=range(10), models=None) -> pd.DataFrame:
    """E2a (design='balanced') / E2b (design='size_matched'): rerun on random subsets.

    The transform is recomputed on each subset, so cohort-wise ranks see the
    subset's own class balance.
    """
    keep_share = {"balanced": False, "size_matched": True}[design]
    parts = []
    for seed in seeds:
        idx = balanced_subsample(y, cohort, seed, keep_share=keep_share)
        Xs, ys, cs = X.iloc[idx], y.iloc[idx], cohort.iloc[idx]
        for transform in TRANSFORMS:
            folds, _ = both_schemes_with_scores(Xs, ys, cs, transform, models)
            folds.insert(0, "seed", seed)
            folds.insert(0, "design", design)
            parts.append(folds)
    return pd.concat(parts, ignore_index=True)


def summarize_balance(df: pd.DataFrame) -> pd.DataFrame:
    """Per subset: gap = mean random-CV AUC (5 folds) - mean held-out AUC (5 cohorts)."""
    keys = ["design", "seed", "transform", "model"]
    per = df.groupby(keys + ["scheme"])["auc_pooled"].mean().unstack("scheme")
    per["gap"] = per["random_cv"] - per["leave_one_cohort_out"]
    pairs = df[df["scheme"] == "random_cv"].groupby(keys)[["auc_within", "auc_cross"]].mean()
    per = per.join(pairs).reset_index()
    return (per.groupby(["design", "transform", "model"])
            .agg(random_cv=("random_cv", "mean"), loco=("leave_one_cohort_out", "mean"),
                 gap_mean=("gap", "mean"), gap_sd=("gap", "std"), gap_min=("gap", "min"),
                 gap_max=("gap", "max"), within=("auc_within", "mean"),
                 cross=("auc_cross", "mean"), subsets=("gap", "count"))
            .reset_index())
