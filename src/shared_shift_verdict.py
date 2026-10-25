"""Quantities and pre-registered rules of the shared-shift test (PLAN 2fd8865).

Thresholds are copied from the PLAN entry committed before any real-data run.
One addition, decided before any run and logged in the RESULT entry: a rule
whose reference damage is zero or negative (cross_raw <= cross_W for H1,
D(I) <= 0 for H2) is VOID, because "remove half the damage" has no meaning.
"""
from __future__ import annotations

import math

import numpy as np
import pandas as pd

from .seesaw import c_for_k

NAN = float("nan")


def _slope(off: pd.DataFrame, arm: str) -> float:
    o = off[off.arm == arm].groupby("cohort").agg(o=("o_test", "mean"), share=("share", "first"),
                                                  pooled=("pooled_share", "first"))
    x = (o.share - o.pooled).to_numpy()
    return abs(float(np.polyfit(x, o.o.to_numpy(), 1)[0])) if np.ptp(x) > 0 else NAN


def quantities(summary: pd.DataFrame, offsets: pd.DataFrame) -> dict:
    """summary = seesaw_verdict.pair_weighted(folds); offsets from run_k5."""
    out = {}
    for cs in summary.set.unique():
        s, off = summary[summary.set == cs], offsets[offsets.set == cs]

        def get(arm, col="D", k=5, C=0.05):
            r = s[(s.arm == arm) & (s.k == k) & np.isclose(s.C, C)]
            return float(r[col].iloc[0]) if len(r) else NAN

        q = {f"D_{a}": get(a) for a in ("W", "I", "I_bal", "raw")}
        q.update({f"cross_{a}": get(a, "cross") for a in ("W", "W_recentre", "raw")})
        q["X2"] = get("W", k=2, C=c_for_k(2)) - get("I", k=2, C=c_for_k(2))
        q["X5"] = q["D_W"] - q["D_I"]
        q["X10"] = get("W", k=10, C=c_for_k(10)) - get("I", k=10, C=c_for_k(10))
        damage = q["cross_raw"] - q["cross_W"]
        q["R"] = (q["cross_W_recentre"] - q["cross_W"]) / damage if damage > 0 else NAN
        w = off[off.arm == "W"].groupby("cohort")[["o_test", "o_train"]].mean()
        q["signs_agree"] = bool((np.sign(w.o_test) == np.sign(w.o_train)).all())
        q["slope_I"], q["slope_I_bal"] = _slope(off, "I"), _slope(off, "I_bal")
        q["max_n_iter"] = float(s.n_iter.max())
        out[cs] = q
    return out


def _band(passed: bool, killed: bool) -> str:
    return "KILL" if killed else "PASS" if passed else "PARTIAL"


def judge(q: dict) -> dict:
    if q["max_n_iter"] >= 3000:
        return {h: "VOID: a logistic fit hit max_iter" for h in ("H1", "H2", "H3")}
    nan = lambda *v: any(isinstance(x, float) and math.isnan(x) for x in v)  # noqa: E731
    if q["cross_raw"] - q["cross_W"] <= 0 or nan(q["R"]):
        h1 = "VOID: no cross damage to remove"
    else:
        h1 = _band(q["signs_agree"] and q["R"] >= 0.5, q["R"] <= 0)
    if q["D_I"] <= 0 or nan(q["D_I_bal"], q["slope_I"], q["slope_I_bal"]):
        h2 = "VOID: no damage under I"
    else:
        h2 = _band(q["D_I_bal"] <= 0.5 * q["D_I"] and q["slope_I_bal"] <= 0.5 * q["slope_I"],
                   q["D_I_bal"] >= 0.9 * q["D_I"])
    if nan(q["X2"], q["X5"], q["X10"]):
        h3 = "VOID: missing number"
    else:
        h3 = _band(abs(q["X10"] - q["X2"]) <= 0.25 * q["X5"] and q["X5"] > 0,
                   q["X10"] >= 1.5 * q["X2"])
    return {"H1": h1, "H2": h2, "H3": h3}


def overall(verdicts: dict) -> str:
    if all(v == "PASS" for per_set in verdicts.values() for v in per_set.values()):
        return "SUPPORTED"
    return "NOT SUPPORTED: " + "; ".join(f"{h} {v} in set '{cs}'" for cs, per in verdicts.items()
                                        for h, v in per.items() if v != "PASS")
