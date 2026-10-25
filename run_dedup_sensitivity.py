"""GSE63990 duplicate-array sensitivity (LABBOOK PLAN 2026-10-21). NOT the headline.

    python run_dedup_sensitivity.py

GSE63990 holds 280 arrays of 273 subjects: "two microarray batches with seven
overlapping samples" (Tsalik 2016). GEO gives no subject ids. The last 7 arrays
(sample_274..sample_280, batch 1) carry exactly the surplus over the published
counts (3 bacterial, 2 viral, 2 non-infectious), so they are the candidate
re-runs. A patient on both sides of a random-CV split inflates random-CV only.
This re-runs the raw headline (fold seed 0, model seed 0) with and without them.

Historical since 2026-10-21: config/cohorts.json now excludes these arrays, so both
arms are the same data; the result on file was made on the 463-array matrix.
"""
from __future__ import annotations

import json
import sys
from datetime import date

import pandas as pd

import run_pipeline as rp
from run_gap_diagnostics import load_matrix
from run_seed_sweep import sweep_rows
from run_seesaw_test import _commit

CANDIDATES = [f"GSE63990:GSM{n}" for n in range(1562133, 1562140)]  # sample_274..280
OUT = rp.RESULTS / "dedup_sensitivity_raw.csv"
META = rp.RESULTS / "run_dedup_sensitivity.json"


def main() -> int:
    X, y, cohort = load_matrix(None)
    drop = [i for i in CANDIDATES if i in X.index]
    print(f"candidates in the labelled matrix: {len(drop)} of {len(CANDIDATES)}: {drop}")
    keep = ~X.index.isin(drop)
    rows = [dict(arm="all_arrays", n=int(len(X)), **r) for r in sweep_rows(X, y, cohort, seed=0)]
    rows += [dict(arm="drop_7_candidates", n=int(keep.sum()), **r)
             for r in sweep_rows(X[keep], y[keep], cohort[keep], seed=0)]
    df = pd.DataFrame(rows)
    wide = df.pivot(index="model", columns="arm", values="gap")
    wide["gap_change"] = wide["drop_7_candidates"] - wide["all_arrays"]
    df.to_csv(OUT, index=False)
    META.write_text(json.dumps(dict(date=date.today().isoformat(), commit=_commit(),
                                    dropped=drop, gap_change=wide["gap_change"].to_dict()),
                               indent=2), encoding="utf-8")
    print(df.to_string(index=False) + "\n\n" + wide.to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
