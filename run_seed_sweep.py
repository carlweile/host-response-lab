"""10-seed sweep of the raw cohort gap (, LABBOOK PLAN 2026-10-21).

    python run_seed_sweep.py                  # repro gate, then seeds 200-209
    python run_seed_sweep.py --cache FILE     # reuse a pickled matrix
    python run_seed_sweep.py --within-sample  # same sweep, within-sample scale

Raw scale by default; --within-sample ranks each patient's own genes first (no LightGBM
model-seed side check there). Raw scale: Sweep: the random-CV fold seed moves, the LightGBM model seed stays
at 0. Leave-one-cohort-out (LOCO) has no fold randomness, so with a fixed model
seed it must be identical in every row (checked). Side check: LOCO for LightGBM
with model seeds 200-209, to show how much LightGBM's own randomness moves it.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

import run_pipeline as rp
from run_gap_diagnostics import load_matrix
from run_seesaw_test import _commit
from src import splits, train

SEEDS = list(range(200, 210))
KILL_GAP = 0.01  #: gap <= 0.01 for BOTH models in most seeds = seed artefact
OUT_SWEEP = rp.RESULTS / "seedsweep_raw.csv"
OUT_MODEL = rp.RESULTS / "seedsweep_lgbm_modelseed.csv"
OUT_SUMMARY = rp.RESULTS / "seedsweep_raw_summary.csv"
META = rp.RESULTS / "run_seedsweep.json"
OUT_SWEEP_WS = rp.RESULTS / "seedsweep_withinsample.csv"
OUT_SUMMARY_WS = rp.RESULTS / "seedsweep_withinsample_summary.csv"
META_WS = rp.RESULTS / "run_seedsweep_withinsample.json"


def sweep_rows(X, y, cohort, seed: int, model_seed: int = 0) -> list[dict]:
    """One row per model: mean random-CV AUC, mean LOCO AUC, gap."""
    s = train.summarize(train.run_both_schemes(X, y, cohort, seed=seed, model_seed=model_seed))
    want = {"random_cv": 5, "leave_one_cohort_out": cohort.nunique()}
    short = s[s["folds"] != s["scheme"].map(want)]
    if not short.empty:  # an undefined-AUC fold was dropped from a mean
        raise ValueError(f"incomplete folds: {short[['model', 'scheme', 'folds']].to_dict('records')}")
    wide = s.pivot(index="model", columns="scheme", values="mean_auc")
    return [dict(fold_seed=seed, model_seed=model_seed, model=m,
                 random_cv_auc=r["random_cv"], loco_auc=r["leave_one_cohort_out"],
                 gap=r["random_cv"] - r["leave_one_cohort_out"])
            for m, r in wide.iterrows()]


def repro_gate(X, y, cohort, ref: str = "summary_raw.csv") -> float:
    """Fold seed 0 + model seed 0 must reproduce results/<ref>."""
    fresh = pd.DataFrame(sweep_rows(X, y, cohort, seed=0)).set_index("model")
    old = pd.read_csv(rp.RESULTS / ref)
    old = old.pivot(index="model", columns="scheme", values="mean_auc")
    if set(old.index) != set(fresh.index):
        return float("inf")
    diffs = [abs(fresh.loc[m, "random_cv_auc"] - old.loc[m, "random_cv"]) for m in old.index]
    diffs += [abs(fresh.loc[m, "loco_auc"] - old.loc[m, "leave_one_cohort_out"]) for m in old.index]
    return float(max(diffs)) if all(np.isfinite(diffs)) else float("inf")


def summarize_gaps(sweep: pd.DataFrame) -> pd.DataFrame:
    g = sweep.groupby("model")["gap"]
    out = g.agg(mean="mean", sd="std", min="min", max="max", n="count").reset_index()
    small = sweep.assign(k=sweep["gap"] <= KILL_GAP).groupby("model")["k"].sum()
    out["seeds_gap_le_kill"] = out["model"].map(small)
    return out


def kill_fired(sweep: pd.DataFrame) -> bool:
    """True when BOTH models have gap <= KILL_GAP in a majority of seeds."""
    per_model = sweep.assign(k=sweep["gap"] <= KILL_GAP).groupby("model")["k"].mean()
    return bool((per_model > 0.5).all())


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache", type=Path, default=None, help="pickled (X, y, cohort)")
    ap.add_argument("--within-sample", action="store_true", help="rank each patient's own genes first")
    args = ap.parse_args(argv)
    ws = args.within_sample
    ref, gate_tol = ("summary_withinsample.csv", 1e-12) if ws else ("summary_raw.csv", 1e-9)

    t0, w0 = time.process_time(), time.time()
    X, y, cohort = load_matrix(args.cache)
    print(f"matrix: {X.shape[0]} samples x {X.shape[1]} genes, {cohort.nunique()} cohorts", flush=True)

    if ws:
        X = splits.rank_normalize_within_sample(X)
    gate = repro_gate(X, y, cohort, ref)
    print(f"repro gate (seed 0 vs {ref}): max |diff| = {gate:.2e}", flush=True)
    if not np.isfinite(gate) or gate > gate_tol:
        print("REPRO GATE FAILED: the sweep would not measure the number on file. Stopping.")
        return 1

    rows = []
    for seed in SEEDS:
        rows += sweep_rows(X, y, cohort, seed)
        print(f"  fold seed {seed} done ({time.time() - w0:.0f} s)", flush=True)
    sweep = pd.DataFrame(rows)
    loco_spread = sweep.groupby("model")["loco_auc"].agg(lambda s: s.max() - s.min())
    if (loco_spread > 1e-12).any():
        print(f"LOCO moved with the fold seed (should not): {loco_spread.to_dict()}")
        return 1

    lgbm = []
    for ms in [] if ws else SEEDS:
        res = train.run_both_schemes(X, y, cohort, seed=0, model_seed=ms)
        aucs = [r.auc for r in res if r.model == "lightgbm" and r.scheme == "leave_one_cohort_out"]
        if len(aucs) != cohort.nunique() or any(a is None for a in aucs):
            raise ValueError(f"model seed {ms}: undefined LOCO AUC")
        lgbm.append(dict(model_seed=ms, model="lightgbm", loco_auc=sum(aucs) / len(aucs)))
        print(f"  lightgbm model seed {ms} done ({time.time() - w0:.0f} s)", flush=True)
    lgbm = pd.DataFrame(lgbm)
    out_sweep, out_summary, meta = (OUT_SWEEP_WS, OUT_SUMMARY_WS, META_WS) if ws else (OUT_SWEEP, OUT_SUMMARY, META)

    summary = summarize_gaps(sweep)
    killed = kill_fired(sweep)
    sweep.to_csv(out_sweep, index=False)
    if not ws:
        lgbm.to_csv(OUT_MODEL, index=False)
    summary.to_csv(out_summary, index=False)
    lgbm_meta = None if ws else dict(mean=lgbm["loco_auc"].mean(), sd=lgbm["loco_auc"].std(),
                                     min=lgbm["loco_auc"].min(), max=lgbm["loco_auc"].max())
    meta.write_text(json.dumps(dict(
        date=date.today().isoformat(), commit=_commit(), seeds=SEEDS, scale="within_sample" if ws else "raw",
        n_samples=int(X.shape[0]), n_genes=int(X.shape[1]), repro_gate_max_diff=gate,
        kill_gap=KILL_GAP, kill_fired=killed,
        lgbm_modelseed_loco=lgbm_meta,
        cpu_s=round(time.process_time() - t0, 1), wall_s=round(time.time() - w0, 1),
    ), indent=2), encoding="utf-8")
    print("\n" + sweep.to_string(index=False) + "\n\n" + summary.to_string(index=False))
    print(f"\nkill rule fired: {killed}\nwrote {out_sweep}, {out_summary}, {meta}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
