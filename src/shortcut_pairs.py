"""Within-cohort-pair AUC vs pooled AUC (LABBOOK PLAN 2026-11-13).

If pooled AUC is no higher than the AUC over same-cohort pairs, knowing the hospital
is not what lifts the headline: a same-cohort pair cannot be separated by the cohort.
All fits go through the E1 path of `gap_experiments` (the one that reproduces
`results/gapdiag_folds.csv`). Thresholds are the pre-registered, fixed in the PLAN entry.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .gap_diagnostics import PairAUC, pair_auc
from .gap_experiments import _fit, _logit
from .splits import random_splits

NOT_USED_MAX = 0.01      # mean D <= this -> NOT-USED
USED_MIN = 0.02          # mean D >= this -> USED
INSTRUMENT_LOW, INSTRUMENT_HIGH = 0.45, 0.55   # shuffled within-pair mean must sit in here
GATE_TOL = 1e-9


@dataclass(frozen=True)
class Fold:
    name: str
    te: np.ndarray      # positional test indices
    p: np.ndarray       # predicted probability of "bacterial"
    logit: np.ndarray   # log-odds score (no clipping)


def fold_scores(Xv: np.ndarray, y: pd.Series, model, fold_seed: int) -> list[Fold]:
    """Fit on each training fold, keep every test score. Xv is float64, rows aligned with y."""
    yv = y.to_numpy().astype(int)
    out = []
    for name, (tr, te) in random_splits(y, seed=fold_seed):
        est = _fit(model, Xv[tr], yv[tr])
        out.append(Fold(name, te, est.predict_proba(Xv[te])[:, 1], _logit(est, Xv[te])))
    return out


def aggregate(pairs: list[PairAUC]) -> dict:
    """Pair-weighted mix over folds: pooled by all pairs, within by same-cohort, cross by cross."""
    n_w = sum(p.n_within for p in pairs)
    n_c = sum(p.n_cross for p in pairs)
    pooled = sum(p.pooled * (p.n_within + p.n_cross) for p in pairs) / (n_w + n_c)
    within = sum(p.within * p.n_within for p in pairs if p.n_within) / n_w
    cross = sum(p.cross * p.n_cross for p in pairs if p.n_cross) / n_c
    return dict(pooled=pooled, within=within, cross=cross, D=pooled - within,
                n_within=n_w, n_cross=n_c)


def seed_row(y: pd.Series, cohort: pd.Series, folds: list[Fold], score: str = "p",
             offset: np.ndarray | None = None) -> dict:
    """One fold seed: per-fold pair_auc on `score` ("p" or "logit"), optional per-patient offset."""
    yv, cv = y.to_numpy().astype(int), cohort.to_numpy()
    pairs = []
    for f in folds:
        s = getattr(f, score)
        if offset is not None:
            s = s + offset[f.te]
        pairs.append(pair_auc(yv[f.te], s, cv[f.te]))
    return aggregate(pairs)


def verdict_word(mean_d: float) -> str:
    """The pre-registered reading rule on the plain mean of D over seeds; edges inclusive."""
    if mean_d <= NOT_USED_MAX:
        return "NOT-USED"
    if mean_d >= USED_MIN:
        return "USED"
    return "UNCLEAR"


def instrument_ok(mean_within: float) -> bool:
    """Shuffled-label within-pair mean must be about 0.50, else the metric is not blind to biology-free data."""
    return INSTRUMENT_LOW <= mean_within <= INSTRUMENT_HIGH


def lightgbm_direction(mean_d: float) -> str:
    return "pooled above within" if mean_d > 0 else "pooled at or below within"


def gate0_rows(folds_by_model: dict[str, list[Fold]], y: pd.Series, cohort: pd.Series,
               ref: pd.DataFrame) -> dict:
    """Largest abs diff to the file, per model, over every fold; pair counts must match exactly."""
    yv, cv = y.to_numpy().astype(int), cohort.to_numpy()
    ref = ref[(ref["transform"] == "raw") & (ref["scheme"] == "random_cv")]
    worst, count_bad = {}, 0
    for model, folds in folds_by_model.items():
        r = ref[ref["model"] == model].set_index("fold")
        if len(r) != len(folds):
            return dict(ok=False, reason=f"{model}: {len(r)} rows on file, {len(folds)} folds")
        d = 0.0
        for f in folds:
            pa, row = pair_auc(yv[f.te], f.p, cv[f.te]), r.loc[f.name]
            d = max(d, abs(pa.pooled - row.auc_pooled), abs(pa.within - row.auc_within),
                    abs(pa.cross - row.auc_cross))
            count_bad += int(pa.n_within != row.n_pairs_within) + int(pa.n_cross != row.n_pairs_cross)
        worst[model] = d
    ok = count_bad == 0 and all(v <= GATE_TOL for v in worst.values())  # NaN fails too
    return dict(ok=bool(ok), max_abs_diff=worst, pair_count_mismatches=count_bad)
