"""Within-cohort-pair AUC vs pooled AUC, fold seeds 0-9 (LABBOOK PLAN 2026-11-13).

    python run_shortcut_pairs.py --cache FILE   # Gate 0, instrument check, main read, controls

Gate 0 first: seed 0, both models, must reproduce results/gapdiag_folds.csv (raw random_cv)
to 1e-9, else exit 1 and no file is written. If the instrument check fails, only the JSON
is written, with verdict "NO VERDICT (instrument)".
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
from src import shortcut_extras as ex
from src import shortcut_pairs as sp
from src.train import make_models

PLAN_COMMIT = "65a9d0b"
ADDENDUM_COMMIT = "892e79e"
SEEDS = list(range(10))
OUT_CSV = rp.RESULTS / "shortcut_pairs_raw.csv"
OUT_PLANTED = rp.RESULTS / "shortcut_pairs_planted.csv"
OUT_JSON = rp.RESULTS / "run_shortcut_pairs.json"
REF_FOLDS = rp.RESULTS / "gapdiag_folds.csv"
REF_NULL = rp.RESULTS / "label_null_raw.csv"


def _instrument(Xv, y, cohort, model, seeds=SEEDS, ref_null: Path | None = REF_NULL) -> dict:
    rows = []
    for s in seeds:
        y_new = ln.shuffle_within_cohort(y, cohort, s)
        ln.check_shuffle(y, y_new, cohort, "within_cohort")
        rows.append(dict(seed=s, **sp.seed_row(y_new, cohort, sp.fold_scores(Xv, y_new, model, 0))))
    df = pd.DataFrame(rows)
    out = dict(within_mean=float(df["within"].mean()), within_sd=float(df["within"].std(ddof=1)),
               pooled_mean=float(df["pooled"].mean()), per_seed_within=df["within"].tolist())
    if ref_null is not None and ref_null.exists():
        old = pd.read_csv(ref_null)
        old = old[old["null"] == "within_cohort"].set_index("seed")["mean_auc"]
        out["pooled_diff_vs_label_null_raw"] = [float(p - old[s]) for s, p in zip(df["seed"], df["pooled"])]
    out["ok"] = sp.instrument_ok(out["within_mean"])
    return out


def run(X, y, cohort, models: dict | None = None, out_csv: Path = OUT_CSV,
        out_planted: Path = OUT_PLANTED, out_json: Path = OUT_JSON, ref_folds: Path = REF_FOLDS,
        ref_null: Path | None = REF_NULL, seeds: list[int] = SEEDS, boot_b: int = ex.BOOT_B) -> int:
    t0 = time.time()
    models = models or make_models()
    Xv = X.to_numpy(dtype=float)
    by_seed = {m: {s: sp.fold_scores(Xv, y, mod, s) for s in seeds[:1]} for m, mod in models.items()}
    g0 = sp.gate0_rows({m: v[seeds[0]] for m, v in by_seed.items()}, y, cohort, pd.read_csv(ref_folds))
    print(f"Gate 0: {g0}", flush=True)
    if not g0["ok"]:
        print("GATE 0 FAILED: not the numbers on file. Stopping, no files written.")
        return 1
    meta = dict(date=date.today().isoformat(), commit=_commit(), plan_commit=PLAN_COMMIT, plan_addendum_commit=ADDENDUM_COMMIT,
                n_arrays=int(X.shape[0]), n_genes=int(X.shape[1]), gate0=g0)
    inst = _instrument(Xv, y, cohort, models[ln.MODEL], seeds, ref_null)
    print(f"Instrument: within mean {inst['within_mean']:.4f}, pooled {inst['pooled_mean']:.4f} "
          f"({time.time() - t0:.0f} s)", flush=True)
    if not inst["ok"]:
        out = dict(meta, instrument=inst, verdict="NO VERDICT (instrument)", wall_s=round(time.time() - t0, 1))
        out_json.write_text(json.dumps(out, indent=2), encoding="utf-8")
        print("INSTRUMENT CHECK FAILED: no verdict, no main read.")
        return 0
    rows = []
    for m, mod in models.items():
        for s in seeds:
            folds = by_seed[m].get(s) or sp.fold_scores(Xv, y, mod, s)
            by_seed[m][s] = folds
            rows.append(dict(seed=s, model=m, **sp.seed_row(y, cohort, folds)))
        print(f"  {m} done ({time.time() - t0:.0f} s)", flush=True)
    df = pd.DataFrame(rows)
    mean_d = {m: float(df[df["model"] == m]["D"].mean()) for m in models}
    sd_d = {m: float(df[df["model"] == m]["D"].std(ddof=1)) for m in models}
    means = {m: df[df["model"] == m][["pooled", "within", "cross", "D"]].mean().to_dict() for m in models}
    boot = {m: ex.bootstrap_mean_d(y, cohort, by_seed[m], B=boot_b) for m in models}
    planted = ex.planted_rows(y, cohort, by_seed[ln.MODEL])
    oracle = {s: ex.oracle_offset_D(y, cohort, by_seed[ln.MODEL][s]) for s in seeds}
    d_max = float(np.mean([o["D"] for o in oracle.values()]))
    out = dict(meta, seeds=[seeds[0], seeds[-1]], instrument=inst, means=means, mean_D=mean_d,
               sd_D_across_seeds=sd_d, bootstrap_mean_D=boot, planted_control=ex.planted_summary(planted),
               headroom=dict(mean_D_max=d_max, usable_line_reachable=bool(d_max >= sp.USED_MIN),
                             per_seed=oracle),
               verdict=sp.verdict_word(mean_d[ln.MODEL]),
               lightgbm_direction=sp.lightgbm_direction(mean_d["lightgbm"]) if "lightgbm" in mean_d else None,
               wall_s=round(time.time() - t0, 1))
    df.to_csv(out_csv, index=False)
    planted.to_csv(out_planted, index=False)
    out_json.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2) + f"\nwrote {out_csv}, {out_planted}, {out_json}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", type=Path, default=None)
    X, y, cohort = load_matrix(ap.parse_args().cache)
    return run(X, y, cohort)


if __name__ == "__main__":
    raise SystemExit(main())
