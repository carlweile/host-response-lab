"""GSE63990 batch split (LABBOOK 2026-11-07 PLAN + 2026-11-08 clarifications, the batch-split run).

Pure helpers. Models come from `transfer.make_model` (settings equal `src/train.py`).
Every batch-out fit asserts no batch sits on both sides (`assert_no_cohort_overlap`, with
the BATCH vector as the group vector); the random control mixes batches on purpose and
has its own checks (disjoint, covers all arrays, class counts equal the targets).
"""
from __future__ import annotations

import re
from collections import Counter

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

from .splits import assert_no_cohort_overlap

ACCESSION = "GSE63990"
EXPECTED = {("0", 1): 58, ("0", 0): 88, ("1", 1): 12, ("1", 0): 27}   # (batch, bacterial) -> n
N_DRAWS = 20
N_BOOT = 1000
_BATCH = re.compile(r'^"?batch:\s*(\S+?)"?$')


def parse_batches(chars: dict[str, list[str]], samples) -> dict[str, str]:
    """sample -> batch string. Exactly one `batch:` string per array, else raise."""
    out = {}
    for s in samples:
        found = {m.group(1) for c in chars.get(s, []) if (m := _BATCH.match(c.strip()))}
        if len(found) != 1:
            raise ValueError(f"{s}: expected exactly one batch string, found {sorted(found)}")
        out[s] = found.pop()
    return out


def batch_counts(yv, bv) -> dict:
    return dict(Counter(zip(map(str, bv), map(int, yv))))


def assert_expected_counts(yv, bv, expected: dict = EXPECTED) -> None:
    got = batch_counts(yv, bv)
    if got != expected:
        raise AssertionError(f"batch x class counts {sorted(got.items())} != {sorted(expected.items())}")


def auc_of(y, s) -> float:
    return float(roc_auc_score(y, s)) if len(set(y)) == 2 else float("nan")


def fit_scores(est, Xv: np.ndarray, yv: np.ndarray, groups: np.ndarray,
               tr: np.ndarray, te: np.ndarray, mixed: bool = False) -> np.ndarray:
    """Fit on rows `tr`, return test scores. Unless `mixed`, raises if a group is on both sides."""
    if not mixed:
        assert_no_cohort_overlap(pd.Series(groups), (tr, te))
    est = clone(est)
    est.fit(Xv[tr], yv[tr])
    return est.predict_proba(Xv[te])[:, 1]


def gate_b(est, Xv: np.ndarray, bv: np.ndarray, seed: int = 0) -> dict:
    """Predict BATCH from expression (infection labels unused): mean of 5 fold AUCs + pooled."""
    t = (bv == "1").astype(int)
    folds, oof = [], np.zeros(len(t))
    for tr, te in StratifiedKFold(5, shuffle=True, random_state=seed).split(Xv, t):
        oof[te] = fit_scores(est, Xv, t, bv, tr, te, mixed=True)
        folds.append(auc_of(t[te], oof[te]))
    return dict(fold_aucs=folds, mean_auc=float(np.mean(folds)), pooled_auc=auc_of(t, oof))


def batch_rows(bv: np.ndarray) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """direction -> (train rows, test rows): b0_to_b1 trains on batch 0, b1_to_b0 the reverse."""
    b0, b1 = np.flatnonzero(bv == "0"), np.flatnonzero(bv == "1")
    return {"b0_to_b1": (b0, b1), "b1_to_b0": (b1, b0)}


def partition(yv: np.ndarray, draw: int, n_pos_small: int, n_neg_small: int) -> tuple[np.ndarray, np.ndarray]:
    """(big side, small side): the small side has the held-out batch's class counts, drawn
    with default_rng(draw); the big side is everything else. Both sort ascending."""
    rng = np.random.default_rng(draw)
    pos, neg = np.flatnonzero(yv == 1), np.flatnonzero(yv == 0)
    small = np.sort(np.concatenate([rng.choice(pos, n_pos_small, replace=False),
                                    rng.choice(neg, n_neg_small, replace=False)]))
    return np.setdiff1d(np.arange(len(yv)), small), small


def check_partition(yv, bv, big, small, n_pos_small: int, n_neg_small: int) -> None:
    """The control must match the batch split's n and class counts, be disjoint, cover all
    arrays, and mix batches on both sides."""
    if set(big) & set(small) or len(big) + len(small) != len(yv):
        raise AssertionError("control sides overlap or do not cover all arrays")
    if (int(yv[small].sum()), int((1 - yv[small]).sum())) != (n_pos_small, n_neg_small):
        raise AssertionError("small control side is not class-matched to the held-out batch")
    if int(yv[big].sum()) != int(yv.sum()) - n_pos_small or len(big) != len(yv) - n_pos_small - n_neg_small:
        raise AssertionError("big control side is not class-matched to the training batch")
    for side in (big, small):
        if len(set(bv[side])) < 2:
            raise AssertionError("a control side does not mix both batches")


def bootstrap_ci(y: np.ndarray, s: np.ndarray, n: int = N_BOOT, seed: int = 0) -> dict:
    """95 % percentile interval of the AUC over test-array resamples (model not refitted)."""
    rng = np.random.default_rng(seed)
    vals, skipped = [], 0
    for _ in range(n):
        i = rng.integers(0, len(y), len(y))
        if len(set(y[i])) < 2:
            skipped += 1
            continue
        vals.append(roc_auc_score(y[i], s[i]))
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return dict(lo=float(lo), hi=float(hi), skipped=skipped)


def run_arm_a(est, Xv, yv, bv) -> dict:
    """Both directions, batch-out. Returns direction -> dict(auc, scores, tr, te)."""
    out = {}
    for d, (tr, te) in batch_rows(bv).items():
        s = fit_scores(est, Xv, yv, bv, tr, te)
        out[d] = dict(auc=auc_of(yv[te], s), scores=s, tr=tr, te=te)
    return out


def run_control(est, Xv, yv, bv, n_pos_small: int, n_neg_small: int, n_draws: int = N_DRAWS) -> pd.DataFrame:
    """Per draw and direction: train n/classes = training batch's, test n/classes = held-out
    batch's, both sides mixing batches. Direction b0_to_b1 trains on the big side."""
    rows = []
    for d in range(n_draws):
        big, small = partition(yv, d, n_pos_small, n_neg_small)
        check_partition(yv, bv, big, small, n_pos_small, n_neg_small)
        for name, (tr, te) in (("b0_to_b1", (big, small)), ("b1_to_b0", (small, big))):
            s = fit_scores(est, Xv, yv, bv, tr, te, mixed=True)
            rows.append(dict(direction=name, draw=d, n_train=len(tr), n_test=len(te),
                             auc=auc_of(yv[te], s)))
    return pd.DataFrame(rows)
