"""The pre-registered rules of the 2026-11-08 PLAN (the deciding batch test), as pure functions.

Thresholds are copied from the PLAN entry committed before any code or run (0e529bd).
Do not tune after seeing data. Count form, not interpolated percentiles, so every edge is exact.
"""
from __future__ import annotations

import math

from .batch_split_verdict import DIRECTIONS, EPS, GATE_B_MIN, GATE_B_STOP, gate_b_passes, r_scale  # noqa: F401

COSTS_SHARE = 0.05       # k_le / N <= this in every direction = "costs AUC"
NO_COST_SHARE = 0.25     # k_lt / N >= this in every direction = "no batch cost"
MIN_BACTERIAL = 5        # power gate, per held-out date
MIN_ARRAYS = 10
PRIMARY_N = 1000
CARD_N = 200             # the plan's literal count; descriptive only


def tail_counts(obs: float, draws) -> dict:
    """k_le: draws at or below obs (a tie counts); k_lt: draws strictly below obs."""
    d = list(draws)
    return dict(k_le=sum(1 for x in d if x <= obs + EPS), k_lt=sum(1 for x in d if x < obs - EPS), n=len(d))


def wilson(k: int, n: int, z: float = 1.959964) -> tuple[float, float]:
    """95 % Wilson interval of a tail share k/n (how much of it is Monte-Carlo noise)."""
    p = k / n
    den = 1 + z * z / n
    mid = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return mid - half, mid + half


def r_batch_n(tails: dict, keys=DIRECTIONS) -> str:
    """tails: key -> tail_counts(...). Needs every key (both directions in Arm A)."""
    if set(tails) != set(keys):
        raise ValueError(f"R-batch-1000 needs exactly {sorted(keys)}")
    if any(t["n"] <= 0 for t in tails.values()):
        raise ValueError("empty null")
    if all(t["k_le"] / t["n"] <= COSTS_SHARE + EPS for t in tails.values()):
        return "costs AUC"
    if all(t["k_lt"] / t["n"] >= NO_COST_SHARE - EPS for t in tails.values()):
        return "no batch cost"
    return "undecided"


def power_gate(counts: dict) -> dict:
    """counts: date -> (n_bacterial, n_viral). Every date must have >= 5 bacterial AND >= 10 arrays."""
    bad = {d: (b, v) for d, (b, v) in counts.items() if b < MIN_BACTERIAL or b + v < MIN_ARRAYS}
    return dict(passes=not bad and len(counts) >= 2, failing=bad)


def decide(gate_b_auc: float | None, tails: dict | None = None, raw_deltas: dict | None = None,
           ws_deltas: dict | None = None) -> dict:
    """Gate B first (missing = error, below bar = stop text only), then R-batch-1000, then R-scale."""
    if not gate_b_passes(gate_b_auc):
        return {"gate_b": "stop", "r_batch": GATE_B_STOP, "r_scale": "not evaluated"}
    if tails is None or raw_deltas is None or ws_deltas is None:
        raise ValueError("Gate B passed but the arm results are missing")
    rb = r_batch_n(tails)
    return {"gate_b": "pass", "r_batch": rb, "r_scale": r_scale(rb, raw_deltas, ws_deltas)}
