"""Control 2 (planted hospital offset) and the patient bootstrap for D.

Both read finished fold scores; nothing is refitted here.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .gap_diagnostics import PairAUC
from .shortcut_pairs import USED_MIN, Fold, aggregate, seed_row

K_GRID = (0.1, 0.25, 0.5, 1.0)
BOOT_B = 2000
SENS_K = 0.5
ORACLE_GRID = np.linspace(-1.0, 1.0, 41)
ORACLE_SWEEPS = 4


def planted_offset(y: pd.Series, cohort: pd.Series, folds: list[Fold], k: float) -> np.ndarray:
    """k x sd(out-of-fold logits) x z; z = patient's cohort bacterial share, standardised over patients.

    Constant inside a cohort by construction, so it cannot touch the within-cohort pairs.
    """
    logits = np.concatenate([f.logit for f in folds])
    share = y.groupby(cohort.to_numpy()).transform("mean").to_numpy(dtype=float)
    z = (share - share.mean()) / share.std()
    return k * logits.std() * z


def planted_rows(y: pd.Series, cohort: pd.Series, folds_by_seed: dict[int, list[Fold]]) -> pd.DataFrame:
    rows = []
    for seed, folds in folds_by_seed.items():
        for k in (0.0,) + K_GRID:
            off = None if k == 0.0 else planted_offset(y, cohort, folds, k)
            rows.append(dict(seed=seed, k=k, **seed_row(y, cohort, folds, "logit", off)))
    return pd.DataFrame(rows)


def planted_summary(df: pd.DataFrame) -> dict:
    """Invariance (within does not move) and sensitivity (|D(k=0.5) - D(0)| >= 0.02)."""
    m = df.groupby("k")[["pooled", "within", "D"]].mean()
    move = float(abs(m.loc[SENS_K, "D"] - m.loc[0.0, "D"]))
    return dict(mean_D_by_k={str(k): float(m.loc[k, "D"]) for k in m.index},
                within_shift=float((m["within"] - m.loc[0.0, "within"]).abs().max()),
                sensitivity_move_at_k=move, sensitivity_k=SENS_K, passed=bool(move >= USED_MIN))


def oracle_offset_D(y: pd.Series, cohort: pd.Series, folds: list[Fold]) -> dict:
    """Upper bound: best label-aware per-hospital constant offset, objective = pooled AUC.

    Start at zero, coordinate grid (ORACLE_GRID x sd of out-of-fold logits), ORACLE_SWEEPS sweeps.
    Only ever accepts a strict gain, so the result is never below the no-offset D.
    """
    cv = cohort.to_numpy()
    names = list(pd.unique(cv))
    sd = float(np.concatenate([f.logit for f in folds]).std())

    def offset(d: np.ndarray) -> np.ndarray:
        out = np.zeros(len(cv))
        for n, v in zip(names, d):
            out[cv == n] = v * sd
        return out

    d = np.zeros(len(names))
    best = seed_row(y, cohort, folds, "logit")
    for _ in range(ORACLE_SWEEPS):
        for i in range(len(names)):
            for g in ORACLE_GRID:
                e = d.copy()
                e[i] = g
                r = seed_row(y, cohort, folds, "logit", offset(e))
                if r["pooled"] > best["pooled"]:
                    d, best = e, r
    return dict(D=best["D"], pooled=best["pooled"], within=best["within"],
                offsets_in_sd=dict(zip(map(str, names), map(float, d))))


def _weighted_pair(wins: np.ndarray, same: np.ndarray, wp: np.ndarray, wn: np.ndarray) -> PairAUC:
    w = wp[:, None] * wn[None, :]

    def share(mask):
        d = float((w * mask).sum())
        return (float((w * mask * wins).sum() / d) if d else float("nan")), d

    pooled, n_all = share(np.ones_like(wins, dtype=bool))
    within, n_w = share(same)
    cross, n_c = share(~same)
    return PairAUC(pooled, within, cross, n_w, n_c)


def bootstrap_mean_d(y: pd.Series, cohort: pd.Series, folds_by_seed: dict[int, list[Fold]],
                     B: int = BOOT_B, seed: int = 0) -> dict:
    """Patient-level bootstrap of the mean-over-seeds D. One count vector per replicate, all seeds."""
    yv, cv = y.to_numpy().astype(int), cohort.to_numpy()
    blocks = []                                   # (wins, same, pos idx, neg idx) per seed x fold
    for folds in folds_by_seed.values():
        for f in folds:
            yt, gt = yv[f.te], cv[f.te]
            pos, neg = yt == 1, yt == 0
            sp, sn = f.p[pos][:, None], f.p[neg][None, :]
            blocks.append(((sp > sn) + 0.5 * (sp == sn), gt[pos][:, None] == gt[neg][None, :],
                           f.te[pos], f.te[neg]))
    n_folds = len(next(iter(folds_by_seed.values())))
    rng, ds = np.random.default_rng(seed), []
    for _ in range(B):
        cnt = rng.multinomial(len(yv), np.full(len(yv), 1.0 / len(yv))).astype(float)
        pairs = [_weighted_pair(w, s, cnt[ip], cnt[ing]) for w, s, ip, ing in blocks]
        ds.append(np.mean([aggregate(pairs[i:i + n_folds])["D"]
                           for i in range(0, len(pairs), n_folds)]))
    lo, hi = np.percentile(ds, [2.5, 97.5])
    return dict(lo=float(lo), hi=float(hi), B=B, boot_sd=float(np.std(ds, ddof=1)))
