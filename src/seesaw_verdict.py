"""Aggregate the seesaw fold rows and apply the pre-registered rules (PLAN 9daaf08).

The thresholds in `judge` are copied from the PLAN entry committed before any
real-data run. They must not be tuned after looking at results.
"""
from __future__ import annotations

import math

import numpy as np
import pandas as pd

from .seesaw import c_for_k

GROUP = ["set", "part", "model", "arm", "k", "C"]
NAN = float("nan")


def pair_weighted(folds: pd.DataFrame) -> pd.DataFrame:
    """Per seed: same- and cross-cohort AUC pair-weighted over that seed's folds
    (pairs never cross a fold), D = same - cross; then the mean over seeds."""
    f = folds.assign(w_in=folds.auc_within * folds.n_pairs_within,
                     w_cr=folds.auc_cross * folds.n_pairs_cross)
    per_seed = f.groupby(GROUP + ["seed"]).agg(w_in=("w_in", "sum"), w_cr=("w_cr", "sum"),
                                               n_in=("n_pairs_within", "sum"),
                                               n_cr=("n_pairs_cross", "sum"),
                                               algebra=("algebra", "mean"),
                                               pin=("pin", "mean"),
                                               n_iter=("n_iter", "max")).reset_index()
    per_seed["same"] = per_seed.w_in / per_seed.n_in
    per_seed["cross"] = per_seed.w_cr / per_seed.n_cr
    per_seed["D"] = per_seed.same - per_seed.cross
    return per_seed.groupby(GROUP).agg(same=("same", "mean"), cross=("cross", "mean"),
                                       D=("D", "mean"), algebra=("algebra", "mean"),
                                       pin=("pin", "mean"), n_iter=("n_iter", "max"),
                                       seeds=("seed", "nunique")).reset_index()


def offset_slope(offsets: pd.DataFrame, k: int, C: float, cohort_set: str) -> float:
    """|S(k)|: per cohort, (W minus I) test offset averaged over folds and seeds,
    regressed on (cohort bacterial share - pooled share)."""
    o = offsets[(offsets.set == cohort_set) & (offsets.part == "p1") & (offsets.k == k)
                & np.isclose(offsets.C, C) & (offsets.model == "logistic_l2")]
    per = o.groupby(["cohort", "arm"]).agg(o=("o_test", "mean"), share=("share", "first"),
                                           pooled=("pooled_share", "first")).reset_index()
    wide = per.pivot(index="cohort", columns="arm", values="o")
    x = per.drop_duplicates("cohort").set_index("cohort")
    x = (x.share - x.pooled).loc[wide.index]
    if np.ptp(x.to_numpy()) == 0:        # all cohorts at the pooled share: no slope
        return NAN
    return abs(float(np.polyfit(x.to_numpy(), (wide["W"] - wide["I"]).to_numpy(), 1)[0]))


def quantities(summary: pd.DataFrame, offsets: pd.DataFrame, gates: dict) -> dict:
    """Every judged number, per cohort set."""
    out = {}
    for cs in summary.set.unique():
        s = summary[summary.set == cs]

        def D(part, model, arm, k=5, C=None):
            C = c_for_k(k) if C is None else C
            r = s[(s.part == part) & (s.model == model) & (s.arm == arm) & (s.k == k)
                  & np.isclose(s.C, C)]
            return float(r.D.iloc[0]) if len(r) else NAN

        def row(model, arm, col):
            r = s[(s.part == "p3") & (s.model == model) & (s.arm == arm)]
            return float(r[col].iloc[0]) if len(r) else NAN

        q = {f"X{k}": D("p1", "logistic_l2", "W", k) - D("p1", "logistic_l2", "I", k)
             for k in (2, 5, 10)}
        q.update({f"S{k}": offset_slope(offsets, k, c_for_k(k), cs) for k in (2, 5, 10)})
        q["Yc"] = D("p2", "logistic_l2", "center_whole") - D("p2", "logistic_l2", "center_train")
        q["Ys"] = D("p2", "logistic_l2", "scale_whole") - D("p2", "logistic_l2", "scale_train")
        q["X_logistic"] = q["X5"]
        for m in ("lightgbm", "lgbm_stumps"):
            q[f"X_{m}"] = D("p3", m, "W") - D("p3", m, "I")
            q[f"algebra_{m}"] = row(m, "W", "algebra")
            q[f"pin_{m}"] = row(m, "W", "pin")
        p1w = s[(s.part == "p1") & (s.model == "logistic_l2") & (s.arm == "W") & (s.k == 5)
                & np.isclose(s.C, c_for_k(5))]
        q["algebra_logistic"] = float(p1w.algebra.iloc[0]) if len(p1w) else NAN
        q["max_n_iter"] = float(s.n_iter.max())
        q.update({f"gate_{key}": v for key, v in gates.get(cs, {}).items()
                  if key in ("center_ok", "scale_ok")})
        out[cs] = q
    return out


def _band(passed: bool, killed: bool) -> str:
    return "KILL" if killed else "PASS" if passed else "PARTIAL"


def judge(q: dict) -> dict:
    """Pass / kill / partial for one cohort set (PLAN 9daaf08, word for word)."""
    if any(isinstance(v, float) and math.isnan(v) for v in q.values()):
        return {p: "VOID: missing number" for p in ("P1", "P2", "P3")}
    if q["max_n_iter"] >= 3000:
        return {p: "VOID: a logistic fit hit max_iter" for p in ("P1", "P2", "P3")}
    X2, X5, X10 = q["X2"], q["X5"], q["X10"]
    p1 = _band(X5 > 0 and X2 <= 0.5 * X5 and X10 >= X5 and q["S2"] < q["S5"] < q["S10"],
               X5 <= 0 or X2 >= 0.8 * X5 or X10 < X2)
    Yc, Ys = q["Yc"], q["Ys"]
    if not (q["gate_center_ok"] and q["gate_scale_ok"]):
        p2 = "VOID: treatment gate failed"
    else:
        p2 = _band(Yc > 0 and Yc >= 0.5 * X5 and abs(Ys) <= 0.25 * Yc,
                   Yc <= 0 or Yc < 0.25 * X5 or Ys >= 0.5 * Yc)
    Xs, Xl = q["X_lgbm_stumps"], q["X_lightgbm"]
    algebra_ok = (q["algebra_lgbm_stumps"] <= 0.5 * q["algebra_lightgbm"]
                  and q["algebra_logistic"] <= 0.5 * q["algebra_lightgbm"])
    if not algebra_ok:
        p3 = "STOP: algebra gate failed (code broken, not interpretable)"
    elif q["pin_lgbm_stumps"] < 0.5 * q["pin_lightgbm"]:
        p3 = "NOT INTERPRETABLE: fit gate failed"
    else:
        p3 = _band(Xs > 0 and Xs >= 2 * max(Xl, 0) and Xs >= 0.5 * q["X_logistic"], Xs <= Xl)
    return {"P1": p1, "P2": p2, "P3": p3}


def overall(verdicts: dict) -> str:
    """Supported only if every prediction passes in every cohort set."""
    if all(v == "PASS" for per_set in verdicts.values() for v in per_set.values()):
        return "SUPPORTED"
    bad = [f"{p} {v} in set '{cs}'" for cs, per_set in verdicts.items()
           for p, v in per_set.items() if v != "PASS"]
    return "NOT SUPPORTED: " + "; ".join(bad)
