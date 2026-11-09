"""The pre-registered rules of the 2026-11-07 PLAN (the batch-split run), as pure functions.

Thresholds are copied from the PLAN entry + the pre-run clarifications entry
(2026-11-08), both committed before any code or run. Do not tune after seeing data.
"""
from __future__ import annotations

import math

GATE_B_MIN = 0.80        # batch must be visible in the values, else the test is uninformative
COST = 0.03              # delta at or above this (both directions) = "costs AUC"
NO_COST = 0.01           # delta at or below this (both directions) = "no batch cost"
SCALE_HALF = 0.5
EPS = 1e-12              # edges are inclusive
DIRECTIONS = ("b0_to_b1", "b1_to_b0")
GATE_B_STOP = "batch not visible in GEO's values; test uninformative"


def _bad(x) -> bool:
    return x is None or (isinstance(x, float) and math.isnan(x))


def gate_b_passes(batch_auc: float | None) -> bool:
    """Inclusive at 0.80. A missing or NaN Gate-B number is an error, never a pass."""
    if _bad(batch_auc):
        raise ValueError("Gate B was not run or has no number; no verdict may be produced")
    return float(batch_auc) >= GATE_B_MIN - EPS


def below_all(auc: float, draws) -> bool:
    """Strictly below the smallest control draw; a tie is NOT below."""
    return bool(auc < min(draws) - EPS)


def delta(control_aucs, batch_out_auc: float) -> float:
    return float(sum(control_aucs) / len(control_aucs)) - batch_out_auc


def r_batch(deltas: dict, below: dict) -> str:
    """deltas / below: direction -> value (logistic, raw). Needs both directions."""
    if set(deltas) != set(DIRECTIONS) or set(below) != set(DIRECTIONS):
        raise ValueError("R-batch needs both directions")
    if all(deltas[d] >= COST - EPS for d in DIRECTIONS) and all(below[d] for d in DIRECTIONS):
        return "costs AUC"
    if all(deltas[d] <= NO_COST + EPS for d in DIRECTIONS):
        return "no batch cost"
    return "undecided"


def r_scale(batch_verdict: str, raw_deltas: dict, ws_deltas: dict) -> str:
    """Only read if R-batch fired."""
    if batch_verdict != "costs AUC":
        return "not evaluated"
    if all(ws_deltas[d] <= SCALE_HALF * raw_deltas[d] + EPS for d in DIRECTIONS):
        return "mostly a scale offset"
    return "not mostly a scale offset"


def decide(gate_b_auc: float | None, raw_deltas: dict | None = None, raw_below: dict | None = None,
           ws_deltas: dict | None = None) -> dict:
    """Gate B first. Below the bar: only the stop text, no R-batch. Missing: raises."""
    if not gate_b_passes(gate_b_auc):
        return {"gate_b": "stop", "r_batch": GATE_B_STOP, "r_scale": "not evaluated"}
    if raw_deltas is None or raw_below is None or ws_deltas is None:
        raise ValueError("Gate B passed but the arm results are missing")
    rb = r_batch(raw_deltas, raw_below)
    return {"gate_b": "pass", "r_batch": rb, "r_scale": r_scale(rb, raw_deltas, ws_deltas)}
