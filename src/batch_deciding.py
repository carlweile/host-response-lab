"""Deciding test of the batch cause (LABBOOK 2026-11-08 PLAN, the deciding batch test).

Arm A reuses the batch-split run's code path unchanged (`batch_split.run_control` with more draws).
Arm B (GSE40396, leave-one-date-out) is guarded by the power gate: `run_arm_b` raises before
any fit if the gate fails, so a bypass is impossible rather than merely forbidden.
"""
from __future__ import annotations

import re
from collections import Counter

import numpy as np
import pandas as pd

from . import batch_deciding_verdict as v
from . import batch_split as bs

DATE_OUT = "date_out"
_DATE = re.compile(r'^"?hybridization date:\s*(\S+?)"?$')


def parse_dates(chars: dict[str, list[str]], samples) -> dict[str, str]:
    """sample -> hybridization date. Exactly one per array, else raise."""
    out = {}
    for s in samples:
        found = {m.group(1) for c in chars.get(s, []) if (m := _DATE.match(c.strip()))}
        if len(found) != 1:
            raise ValueError(f"{s}: expected exactly one hybridization date, found {sorted(found)}")
        out[s] = found.pop()
    return out


def date_counts(labels: dict[str, str], dates: dict[str, str]) -> dict[str, tuple[int, int]]:
    """date -> (n_bacterial, n_viral) over the labelled arrays only."""
    c = Counter((dates[s], l) for s, l in labels.items() if l in ("bacterial", "viral"))
    return {d: (c[(d, "bacterial")], c[(d, "viral")]) for d in sorted(set(dates.values()))}


def run_null(est, Xv, yv, bv, n_pos: int, n_neg: int, n_draws: int = v.PRIMARY_N) -> pd.DataFrame:
    """The the batch-split run control, only longer: seeds 0..n_draws-1, same partition and per-draw checks."""
    return bs.run_control(est, Xv, yv, bv, n_pos, n_neg, n_draws=n_draws)


def gate_0b(new: pd.DataFrame, old: pd.DataFrame, tol: float = 1e-6) -> dict:
    """Draws 0-19 of the new null must reproduce the committed `batch_control.csv`."""
    key = ["transform", "direction", "draw"]
    m = old.merge(new[key + ["auc"]], on=key, how="left", suffixes=("_old", "_new"))
    if m.auc_new.isna().any():
        return dict(ok=False, max_abs_diff=float("nan"), rows=len(m))
    diff = float((m.auc_old - m.auc_new).abs().max())
    return dict(ok=diff <= tol, max_abs_diff=diff, rows=len(m))


def prefix_tails(obs: float, draws: np.ndarray, n: int) -> dict:
    return v.tail_counts(obs, draws[:n])


def run_arm_b(est, Xv, yv, dv, counts: dict, n_draws: int = v.PRIMARY_N) -> dict:
    """Leave-one-date-out vs size+class-matched mixed splits. Power gate FIRST, before any fit."""
    gate = v.power_gate(counts)
    if not gate["passes"]:
        raise RuntimeError(f"power gate failed {gate['failing']}: Arm B must not run")
    tails, aucs = {}, {}
    for d in sorted(set(dv)):
        te = np.flatnonzero(dv == d)
        tr = np.flatnonzero(dv != d)
        s = bs.fit_scores(est, Xv, yv, dv, tr, te)             # raises if a date is on both sides
        aucs[d] = bs.auc_of(yv[te], s)
        n_pos, n_neg = int(yv[te].sum()), int((1 - yv[te]).sum())
        draws = []
        for k in range(n_draws):
            big, small = bs.partition(yv, k, n_pos, n_neg)
            bs.check_partition(yv, dv, big, small, n_pos, n_neg)
            draws.append(bs.auc_of(yv[small], bs.fit_scores(est, Xv, yv, dv, big, small, mixed=True)))
        tails[d] = v.tail_counts(aucs[d], draws)
    return dict(aucs=aucs, tails=tails, verdict=v.r_batch_n(tails, keys=tuple(tails)))
