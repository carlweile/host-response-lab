"""Probe-choice control (outside review 2026-10-21, astra finding 1; LABBOOK PLAN 2026-10-21).

    python run_probe_control.py

The headline matrix keeps, per gene, the probe with the highest mean over every
patient of the series, test patients included. This re-runs the raw headline
(fold seed 0, model seed 0) with the lowest probe id per gene instead, a rule that
reads no expression value, and puts both gaps side by side.
"""
from __future__ import annotations

import json
import sys
import time
from datetime import date

import pandas as pd

import run_pipeline as rp
from run_gap_diagnostics import load_matrix
from run_seed_sweep import sweep_rows
from run_seesaw_test import _commit

MOVE_LINE = 0.02  # pre-registered: a gap moving more than this = the probe choice matters
OUT = rp.RESULTS / "probecontrol_raw.csv"
META = rp.RESULTS / "run_probecontrol.json"


def main() -> int:
    w0 = time.time()
    rows = []
    for rule in ("mean", "first"):
        X, y, cohort = load_matrix(None, collapse=rule)
        print(f"{rule}: {X.shape[0]} samples x {X.shape[1]} genes", flush=True)
        rows += [dict(collapse=rule, **r) for r in sweep_rows(X, y, cohort, seed=0)]
    df = pd.DataFrame(rows)
    wide = df.pivot(index="model", columns="collapse", values="gap")
    wide["gap_change"] = wide["first"] - wide["mean"]
    moved = bool((wide["gap_change"].abs() > MOVE_LINE).any())
    df.to_csv(OUT, index=False)
    META.write_text(json.dumps(dict(
        date=date.today().isoformat(), commit=_commit(), move_line=MOVE_LINE, moved=moved,
        gap_change=wide["gap_change"].to_dict(), wall_s=round(time.time() - w0, 1),
    ), indent=2), encoding="utf-8")
    print(df.to_string(index=False) + "\n\n" + wide.to_string())
    print(f"\nprobe choice moves a gap by more than {MOVE_LINE}: {moved}\nwrote {OUT}, {META}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
