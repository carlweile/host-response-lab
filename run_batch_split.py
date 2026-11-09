"""Does GSE63990's own array batch already cost AUC? (LABBOOK 2026-11-07 PLAN + 2026-11-08 clarifications, the batch-split run)

    python run_batch_split.py                # Gate 0 -> Gate B -> Arm A -> control -> verdicts (CPU, ~5-10 min)
    python run_batch_split.py --gates-only   # Gate 0 + Gate B only, writes nothing
    python run_batch_split.py --cache FILE   # reuse a pickled matrix

Everything is held in memory and written only at the end. Gate 0 miss = nothing written.
Gate B below 0.80 = only run_batch_split.json + the Gate-B rows of batch_split.csv.
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
from src import batch_split as bs
from src import batch_split_verdict as bv_
from src import parse_geo
from src import transfer as tf
from src.splits import rank_normalize_within_sample

LR, GBM = "logistic_l2", "lightgbm"
GATE0_TOL = 0.005
GBM_SEEDS = (0, 1, 2)
OUT_SPLIT, OUT_CTRL = rp.RESULTS / "batch_split.csv", rp.RESULTS / "batch_control.csv"
META = rp.RESULTS / "run_batch_split.json"
PLAN_COMMIT = "4531829"


def gate0(X: pd.DataFrame, y: pd.Series, cohort: pd.Series, batches: dict) -> dict:
    """Held-out-GSE63990 AUCs must equal the file's; the 185 arrays must count 58/88/12/27."""
    Xv, yv, cv = X.to_numpy(), y.to_numpy(), cohort.to_numpy()
    te, tr = np.flatnonzero(cv == bs.ACCESSION), np.flatnonzero(cv != bs.ACCESSION)
    old = pd.read_csv(rp.RESULTS / "folds_raw.csv")
    old = old[(old.scheme == "leave_one_cohort_out") & (old.fold == bs.ACCESSION)].set_index("model").auc
    diffs = {}
    for m in (LR, GBM):
        s = bs.fit_scores(tf.make_model(m, tf.BALANCED, 0), Xv, yv, cv, tr, te)
        diffs[m] = abs(bs.auc_of(yv[te], s) - float(old[m]))
    names = [i.split(":", 1)[1] for i in X.index[te]]
    b = np.array([batches[n] for n in names])
    bs.assert_expected_counts(yv[te], b)
    return dict(diffs=diffs, max_abs_diff=max(diffs.values()), ok=max(diffs.values()) <= GATE0_TOL,
                batch=b, rows=te)


def comparators() -> dict:
    """Read from the files on disk, never typed in."""
    s = pd.read_csv(rp.RESULTS / "summary_raw.csv")
    gap = float(s[s.model == LR].gap.iloc[0])
    c = pd.read_csv(rp.RESULTS / "gapdiag_cohorts.csv")
    c = c[(c["transform"] == "raw") & (c["model"] == LR) & (c["cohort"] == bs.ACCESSION)].iloc[0]
    return dict(overall_logistic_gap=gap, own_heldout_loss=float(c.rcv_auc_in_cohort - c.loco_auc))


def arm_rows(arm: dict, transform: str, model: str, seed: int, y: np.ndarray) -> list[dict]:
    rows = []
    for d, r in arm.items():
        ci = bs.bootstrap_ci(y[r["te"]], r["scores"]) if model == LR else dict(lo=None, hi=None, skipped=None)
        rows.append(dict(part="arm_a", model=model, transform=transform, direction=d, model_seed=seed,
                         n_train=len(r["tr"]), n_test=len(r["te"]), auc=r["auc"],
                         ci_lo=ci["lo"], ci_hi=ci["hi"], boot_skipped=ci["skipped"]))
    return rows


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gates-only", action="store_true")
    ap.add_argument("--cache", type=Path)
    args = ap.parse_args(argv)
    t0 = time.time()
    X, y, cohort = load_matrix(args.cache)
    on_file = json.loads((rp.RESULTS / "run_raw.json").read_text(encoding="utf-8"))
    if X.shape != (on_file["n_samples"], on_file["n_genes"]):
        raise SystemExit(f"matrix {X.shape} differs from run_raw.json; stop")
    series = parse_geo.parse_series_matrix(rp.DATA / "series" / f"{bs.ACCESSION}_series_matrix.txt.gz", bs.ACCESSION)
    batches = bs.parse_batches(series.characteristics, series.samples)

    g0 = gate0(X, y, cohort, batches)
    print(f"Gate 0 max abs diff {g0['max_abs_diff']:.6f} {g0['diffs']}", flush=True)
    if not g0["ok"]:
        print("GATE 0 MISSED: nothing written to results/")
        return 1
    rows = g0["rows"]
    yv, bvec = y.to_numpy()[rows], g0["batch"]
    Xraw = X.to_numpy()[rows]
    Xws = rank_normalize_within_sample(X.iloc[rows]).to_numpy()

    gb = bs.gate_b(tf.make_model(LR, tf.BALANCED), Xraw, bvec)
    print(f"Gate B batch-AUC (mean of 5 folds) {gb['mean_auc']:.4f}, pooled {gb['pooled_auc']:.4f}", flush=True)
    if args.gates_only:
        return 0
    split_rows = [dict(part="gate_b_fold", model=LR, transform="raw", direction="batch", model_seed=0, auc=a)
                  for a in gb["fold_aucs"]]
    split_rows += [dict(part="gate_b_mean", model=LR, transform="raw", direction="batch", model_seed=0,
                        auc=gb["mean_auc"]),
                   dict(part="gate_b_pooled", model=LR, transform="raw", direction="batch", model_seed=0,
                        auc=gb["pooled_auc"])]
    meta = dict(date=date.today().isoformat(), plan_commit=PLAN_COMMIT, card="the batch-split run",
                gate0=dict(max_abs_diff=g0["max_abs_diff"], diffs=g0["diffs"], tolerance=GATE0_TOL),
                gate_b=dict(mean_auc=gb["mean_auc"], pooled_auc=gb["pooled_auc"], fold_aucs=gb["fold_aucs"],
                            bar=bv_.GATE_B_MIN),
                batch_counts={f"batch{k[0]}_{'bacterial' if k[1] else 'viral'}": v
                              for k, v in sorted(bs.batch_counts(yv, bvec).items())})
    if not bv_.gate_b_passes(gb["mean_auc"]):
        meta["verdicts"] = bv_.decide(gb["mean_auc"])
        return finish(meta, split_rows, None, t0)

    lr = tf.make_model(LR, tf.BALANCED)
    n_pos, n_neg = (int(yv[bvec == "1"].sum()), int((1 - yv[bvec == "1"]).sum()))
    results, ctrl_frames = {}, []
    for tname, Xt in (("raw", Xraw), ("within_sample", Xws)):
        arm = bs.run_arm_a(lr, Xt, yv, bvec)
        split_rows += arm_rows(arm, tname, LR, 0, yv)
        ctrl = bs.run_control(lr, Xt, yv, bvec, n_pos, n_neg)
        ctrl.insert(0, "transform", tname)
        ctrl_frames.append(ctrl)
        results[tname] = {d: dict(auc=arm[d]["auc"],
                                  draws=ctrl[ctrl.direction == d].auc.to_numpy()) for d in bv_.DIRECTIONS}
        print(f"  {tname}: " + ", ".join(f"{d} {results[tname][d]['auc']:.4f} (control mean "
              f"{results[tname][d]['draws'].mean():.4f})" for d in bv_.DIRECTIONS), flush=True)
    gbm_ranges = {}
    for sd in GBM_SEEDS:
        arm = bs.run_arm_a(tf.make_model(GBM, tf.BALANCED, sd), Xraw, yv, bvec)
        split_rows += arm_rows(arm, "raw", GBM, sd, yv)
    g = pd.DataFrame(split_rows)
    g = g[(g.part == "arm_a") & (g.model == GBM)]
    for d in bv_.DIRECTIONS:
        a = g[g.direction == d].auc
        gbm_ranges[d] = dict(min=float(a.min()), max=float(a.max()))

    d_raw = {d: bv_.delta(results["raw"][d]["draws"], results["raw"][d]["auc"]) for d in bv_.DIRECTIONS}
    d_ws = {d: bv_.delta(results["within_sample"][d]["draws"], results["within_sample"][d]["auc"])
            for d in bv_.DIRECTIONS}
    below = {d: bv_.below_all(results["raw"][d]["auc"], results["raw"][d]["draws"]) for d in bv_.DIRECTIONS}
    meta.update(
        summary={t: {d: dict(batch_out_auc=results[t][d]["auc"], control_mean=float(results[t][d]["draws"].mean()),
                             control_min=float(results[t][d]["draws"].min()),
                             control_max=float(results[t][d]["draws"].max()),
                             control_sd=float(results[t][d]["draws"].std(ddof=1)))
                     for d in bv_.DIRECTIONS} for t in results},
        delta_raw=d_raw, delta_within_sample=d_ws, below_all_20_raw=below,
        lightgbm_raw_ranges_seeds_0_2=gbm_ranges, comparators=comparators(),
        verdicts=bv_.decide(gb["mean_auc"], d_raw, below, d_ws))
    return finish(meta, split_rows, pd.concat(ctrl_frames, ignore_index=True), t0)


def finish(meta: dict, split_rows: list, ctrl: pd.DataFrame | None, t0: float) -> int:
    df = pd.DataFrame(split_rows)
    if df.auc.isna().any() or (ctrl is not None and ctrl.auc.isna().any()):
        raise AssertionError("NaN AUC: nothing written")
    meta["runtime_s"] = round(time.time() - t0, 1)
    df.to_csv(OUT_SPLIT, index=False)
    if ctrl is not None:
        ctrl.to_csv(OUT_CTRL, index=False)
    META.write_text(json.dumps(meta, indent=2, default=float), encoding="utf-8")
    print(json.dumps({k: meta[k] for k in meta if k in ("verdicts", "delta_raw", "delta_within_sample",
                                                         "below_all_20_raw", "comparators")}, indent=2, default=float))
    return 0


if __name__ == "__main__":
    sys.exit(main())
