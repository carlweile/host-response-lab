"""Label-shuffle nulls of the raw random-CV headline (LABBOOK PLAN 2026-11-12).

    python run_label_null.py                 # Gate 0, then 2 nulls x 50 seeds
    python run_label_null.py --cache FILE    # reuse a pickled matrix

Gate 0: unshuffled logistic random-CV and held-out must reproduce results/summary_raw.csv
to 1e-12, else exit 1 before any shuffle and no file is written.
"""
from __future__ import annotations

import argparse
import json
import time
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

import run_pipeline as rp
from run_gap_diagnostics import load_matrix
from run_seesaw_test import _commit
from src import label_null as ln

PLAN_COMMIT = "81252f2"
SEEDS = list(range(50))
GATE_TOL = 1e-12
OUT_CSV = rp.RESULTS / "label_null_raw.csv"
OUT_JSON = rp.RESULTS / "run_label_null.json"
SUMMARY = rp.RESULTS / "summary_raw.csv"


def gate0(X, y, cohort, ref: Path = SUMMARY) -> dict:
    old = pd.read_csv(ref)
    old = old[old["model"] == ln.MODEL].set_index("scheme")["mean_auc"]
    rcv = ln.random_cv_auc(X, y)["mean_auc"]
    loco = ln.loco_auc(X, y, cohort)
    return dict(random_cv=abs(rcv - float(old["random_cv"])),
                held_out=abs(loco - float(old["leave_one_cohort_out"])),
                held_out_on_file=float(old["leave_one_cohort_out"]))


def stats(aucs: list[float]) -> dict:
    a = np.asarray(aucs)
    return dict(mean=float(a.mean()), p95=float(np.percentile(a, 95)), max=float(a.max()),
                min=float(a.min()), sd=float(a.std(ddof=1)))


def run(X, y, cohort, out_csv: Path = OUT_CSV, out_json: Path = OUT_JSON, ref: Path = SUMMARY,
        seeds: list[int] = SEEDS) -> int:
    t0 = time.time()
    g0 = gate0(X, y, cohort, ref)
    print(f"Gate 0: random_cv diff {g0['random_cv']:.2e}, held_out diff {g0['held_out']:.2e}", flush=True)
    if not (g0["random_cv"] <= GATE_TOL and g0["held_out"] <= GATE_TOL):  # NaN fails too
        print("GATE 0 FAILED: the null would not measure the number on file. Stopping, no files written.")
        return 1
    rows = []
    for kind, fn in (("global", lambda s: ln.shuffle_global(y, s)),
                     ("within_cohort", lambda s: ln.shuffle_within_cohort(y, cohort, s))):
        for seed in seeds:
            y_new = fn(seed)
            changed = ln.check_shuffle(y, y_new, cohort, kind)
            rows.append(dict(null=kind, seed=seed, **ln.random_cv_auc(X, y_new), n_labels_changed=changed))
        print(f"  {kind} done ({time.time() - t0:.0f} s)", flush=True)
    df = pd.DataFrame(rows)
    per = {k: stats(df[df["null"] == k]["mean_auc"].tolist()) for k in ("global", "within_cohort")}
    out = dict(date=date.today().isoformat(), commit=_commit(), plan_commit=PLAN_COMMIT,
               n_arrays=int(X.shape[0]), n_genes=int(X.shape[1]), seeds=[seeds[0], seeds[-1]],
               gate0=g0, nulls=per, verdict=ln.verdict(per["global"]["mean"], per["global"]["max"]),
               shortcut_word=ln.shortcut_word(per["within_cohort"]["mean"]),
               cohort_only_reference_auc=ln.cohort_only_auc(y, cohort),
               wall_s=round(time.time() - t0, 1))
    df.to_csv(out_csv, index=False)
    out_json.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2) + f"\nwrote {out_csv}, {out_json}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", type=Path, default=None)
    args = ap.parse_args()
    X, y, cohort = load_matrix(args.cache)
    return run(X, y, cohort)


if __name__ == "__main__":
    raise SystemExit(main())
