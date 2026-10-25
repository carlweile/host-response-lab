"""Fit the same models under both validation schemes and report the gap.

No feature selection happens outside the cross-validation loop, on purpose:
picking genes on the full dataset first is the most common way this exact
experiment is accidentally rigged.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass, asdict

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from .splits import assert_no_cohort_overlap, cohort_splits, random_splits

warnings.filterwarnings("ignore", category=UserWarning)


def make_models(model_seed: int = 0) -> dict[str, object]:
    models: dict[str, object] = {
        "logistic_l2": Pipeline([
            ("scale", StandardScaler()),
            ("clf", LogisticRegression(C=0.05, max_iter=3000, class_weight="balanced")),
        ]),
    }
    try:
        from lightgbm import LGBMClassifier
        models["lightgbm"] = LGBMClassifier(
            n_estimators=300, learning_rate=0.05, num_leaves=15,
            min_child_samples=10, subsample=0.8, colsample_bytree=0.3,
            class_weight="balanced", verbose=-1, random_state=model_seed,
        )
    except ImportError:
        pass
    return models


@dataclass
class FoldResult:
    scheme: str
    model: str
    fold: str
    n_train: int
    n_test: int
    n_pos_test: int
    auc: float | None
    note: str = ""


def _fit_score(model, X: pd.DataFrame, y: pd.Series, split) -> tuple[float | None, str]:
    tr, te = split
    y_tr, y_te = y.iloc[tr], y.iloc[te]
    if y_tr.nunique() < 2:
        return None, "training fold has one class only"
    if y_te.nunique() < 2:
        return None, "test fold has one class only, AUC undefined"
    from sklearn.base import clone
    est = clone(model)
    est.fit(X.iloc[tr].to_numpy(), y_tr.to_numpy())
    scores = est.predict_proba(X.iloc[te].to_numpy())[:, 1]
    return float(roc_auc_score(y_te.to_numpy(), scores)), ""


def run_both_schemes(
    X: pd.DataFrame, y: pd.Series, cohort: pd.Series, n_splits: int = 5,
    seed: int = 0, model_seed: int = 0,
) -> list[FoldResult]:
    """`seed` shuffles the random-CV folds only; `model_seed` is LightGBM's own
    randomness. Leave-one-cohort-out has no fold randomness at all."""
    if not (X.index.equals(y.index) and X.index.equals(cohort.index)):
        raise ValueError("X, y and cohort must share one index in one order")
    results: list[FoldResult] = []
    models = make_models(model_seed)

    for name, model in models.items():
        for fold, split in random_splits(y, n_splits=n_splits, seed=seed):
            auc, note = _fit_score(model, X, y, split)
            results.append(FoldResult("random_cv", name, fold,
                                      len(split[0]), len(split[1]),
                                      int(y.iloc[split[1]].sum()), auc, note))
        for held_out, split in cohort_splits(cohort):
            assert_no_cohort_overlap(cohort, split)
            auc, note = _fit_score(model, X, y, split)
            results.append(FoldResult("leave_one_cohort_out", name, held_out,
                                      len(split[0]), len(split[1]),
                                      int(y.iloc[split[1]].sum()), auc, note))
    return results


def summarize(results: list[FoldResult]) -> pd.DataFrame:
    rows = [asdict(r) for r in results]
    df = pd.DataFrame(rows)
    scored = df[df["auc"].notna()]
    if scored.empty:
        return pd.DataFrame()
    agg = (scored.groupby(["model", "scheme"])["auc"]
           .agg(mean_auc="mean", min_auc="min", max_auc="max", folds="count")
           .reset_index())
    wide = agg.pivot(index="model", columns="scheme", values="mean_auc")
    if {"random_cv", "leave_one_cohort_out"}.issubset(wide.columns):
        wide["gap"] = wide["random_cv"] - wide["leave_one_cohort_out"]
    return agg.merge(
        wide.reset_index()[["model"] + [c for c in ("gap",) if c in wide.columns]],
        on="model", how="left",
    )
