"""Train-on-one-cohort / test-on-another transfer matrix (LABBOOK 2026-11-05 PLAN, the transfer-matrix run).

Pure helpers: models come unchanged from `train.make_models`, only `class_weight`
is switched with `set_params`. Every fit asserts the test cohort is not in the
training rows; the size-matched arm asserts the cap really was applied.
"""
from __future__ import annotations

from itertools import permutations

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.metrics import roc_auc_score

from .splits import assert_no_cohort_overlap, cohort_splits
from .train import make_models

SIZE_CAP = 30
N_DRAWS = 20
BALANCED = "balanced"


def make_model(name: str, class_weight, model_seed: int = 0):
    """The headline model with class_weight overridden (BALANCED = the on-file setting)."""
    est = clone(make_models(model_seed)[name])
    est.set_params(**({"clf__class_weight": class_weight} if name == "logistic_l2"
                      else {"class_weight": class_weight}))
    return est


def class_weight_of(est):
    """What the fitted-to-be estimator really carries (used by the tests)."""
    return est.named_steps["clf"].class_weight if hasattr(est, "named_steps") else est.class_weight


def cohort_tags(cfg: dict) -> tuple[dict, dict]:
    """(platform_family, age_group) per accession, from config/cohorts.json."""
    fam = {c["accession"]: c["platform_family"] for c in cfg["cohorts"]}
    age = {c["accession"]: c["age_group"] for c in cfg["cohorts"]}
    return fam, age


def ordered_pairs(cohorts) -> list[tuple[str, str]]:
    return list(permutations(sorted(set(cohorts)), 2))


def pair_cell(train: str, test: str, tag: dict) -> str:
    return "same" if tag[train] == tag[test] else "cross"


def age_testable(age: dict) -> bool:
    """R-age needs every cohort pinned to children or adults."""
    return all(a in ("children", "adults") for a in age.values())


def cell_counts(pairs, tag: dict) -> dict:
    out = {"same": 0, "cross": 0}
    for tr, te in pairs:
        out[pair_cell(tr, te, tag)] += 1
    return out


def _auc(y, s) -> float:
    return float(roc_auc_score(y, s)) if len(set(y)) == 2 else float("nan")


def fit_auc(est, Xv: np.ndarray, yv: np.ndarray, cv: np.ndarray,
            tr: np.ndarray, te: np.ndarray) -> float:
    """Fit on rows `tr`, score rows `te`. Raises if any cohort sits on both sides."""
    assert_no_cohort_overlap(pd.Series(cv), (tr, te))
    est = clone(est)
    est.fit(Xv[tr], yv[tr])
    return _auc(yv[te], est.predict_proba(Xv[te])[:, 1])


def pair_rows(cv: np.ndarray, train: str, test: str) -> tuple[np.ndarray, np.ndarray]:
    return np.flatnonzero(cv == train), np.flatnonzero(cv == test)


def stratified_draw(yv: np.ndarray, rows: np.ndarray, draw: int, n: int = SIZE_CAP) -> np.ndarray:
    """n rows from `rows`, bacterial count = round(n x the cohort's bacterial share), seed = draw."""
    rng = np.random.default_rng(draw)
    pos, neg = rows[yv[rows] == 1], rows[yv[rows] == 0]
    n_pos = int(round(n * len(pos) / len(rows)))
    if n_pos > len(pos) or n - n_pos > len(neg):
        raise ValueError("cohort too small for the size cap")
    got = np.concatenate([rng.choice(pos, n_pos, replace=False),
                          rng.choice(neg, n - n_pos, replace=False)])
    if len(got) != n:
        raise AssertionError("size cap not applied")
    return np.sort(got)


def run_pairs(Xv, yv, cv, model: str, class_weight=BALANCED, model_seed: int = 0) -> pd.DataFrame:
    est = make_model(model, class_weight, model_seed)
    rows = []
    for a, b in ordered_pairs(cv):
        tr, te = pair_rows(cv, a, b)
        rows.append(dict(train=a, test=b, n_train=len(tr), n_test=len(te),
                         auc=fit_auc(est, Xv, yv, cv, tr, te)))
    return pd.DataFrame(rows)


def run_sizematched(Xv, yv, cv, model: str = "logistic_l2") -> pd.DataFrame:
    est = make_model(model, BALANCED)
    rows = []
    for a, b in ordered_pairs(cv):
        full_tr, te = pair_rows(cv, a, b)
        for d in range(N_DRAWS):
            tr = stratified_draw(yv, full_tr, d)
            if len(tr) != SIZE_CAP:
                raise AssertionError("n_train != size cap")
            rows.append(dict(train=a, test=b, draw=d, n_train=len(tr),
                             auc=fit_auc(est, Xv, yv, cv, tr, te)))
    return pd.DataFrame(rows)


def run_loco(Xv, yv, cv, cohort: pd.Series, model: str, class_weight, model_seed: int = 0) -> pd.DataFrame:
    est = make_model(model, class_weight, model_seed)
    rows = []
    for held_out, (tr, te) in cohort_splits(cohort):
        rows.append(dict(model=model, class_weight=str(class_weight), model_seed=model_seed,
                         held_out=held_out, auc=fit_auc(est, Xv, yv, cv, tr, te)))
    return pd.DataFrame(rows)


def cell_mean(df: pd.DataFrame, tag: dict, cell: str, col: str = "auc") -> float:
    """Unweighted mean over the ordered pairs of a cell (size-matched: first mean over draws)."""
    per_pair = df.groupby(["train", "test"], as_index=False)[col].mean()
    keep = [pair_cell(a, b, tag) == cell for a, b in zip(per_pair["train"], per_pair["test"])]
    return float(per_pair[col][np.array(keep)].mean())


def run_balance(Xv, yv, cv, cohort: pd.Series, seeds, models=("logistic_l2", "lightgbm")) -> pd.DataFrame:
    """Arm B: LOCO with and without class weighting; logistic once (deterministic), LightGBM per seed."""
    return pd.concat([run_loco(Xv, yv, cv, cohort, m, cw, s)
                      for m in models for cw in (BALANCED, None)
                      for s in ((0,) if m == models[0] else tuple(seeds))], ignore_index=True)


def balance_summary(bal: pd.DataFrame, model_seeded: str = "lightgbm") -> dict:
    """Reading (d): per-cohort mean over seeds, then mean over cohorts; plus the per-seed spread."""
    out = {}
    for m, b in bal.groupby("model"):
        per = b.groupby(["class_weight", "held_out"]).auc.mean().unstack(0)
        d = per["None"] - per[BALANCED]
        out[m] = float(d.mean())
        out[m + "_per_cohort"] = {k: float(v) for k, v in d.items()}
    b = bal[bal.model == model_seeded]
    ps = (b[b.class_weight == "None"].groupby("model_seed").auc.mean()
          - b[b.class_weight == BALANCED].groupby("model_seed").auc.mean())
    out[model_seeded + "_per_seed"] = {int(k): float(v) for k, v in ps.items()}
    out[model_seeded + "_seed_spread"] = dict(min=float(ps.min()), max=float(ps.max()),
                                              sd=float(ps.std(ddof=1)) if len(ps) > 1 else float("nan"))
    return out
