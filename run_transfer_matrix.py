"""Per-cohort drop analysis (LABBOOK 2026-11-05 PLAN + pre-run clarifications, the transfer-matrix run).

    python run_transfer_matrix.py                # Gate 0, then the full run (CPU, ~10-15 min)
    python run_transfer_matrix.py --gate0-only   # reproduction check only, writes nothing
    python run_transfer_matrix.py --cache FILE   # reuse a pickled matrix
    python run_transfer_matrix.py --balance-seeds 10   # the 10-seed Arm B run: Arm B only, 10 LightGBM seeds, new files

Gate 0 first: leave-one-cohort-out with the on-file setting must reproduce all ten AUCs
of results/folds_raw.csv within 0.005, else nothing is written to results/.
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
from src import transfer as tf
from src import transfer_verdict as tv
from src.splits import rank_normalize_within_sample

LR, GBM = "logistic_l2", "lightgbm"
GATE_TOL = 0.005
SEEDS = (0, 1, 2)
OUT = {k: rp.RESULTS / f"{v}.csv" for k, v in
       dict(matrix="transfer_matrix", size="transfer_sizematched", balance="balance_control").items()}
META = rp.RESULTS / "run_transfer.json"
OUT10, META10 = rp.RESULTS / "balance_control_10seeds.csv", rp.RESULTS / "run_balance_10seeds.json"


def gate0(Xv, yv, cv, cohort) -> tuple[bool, pd.DataFrame]:
    old = pd.read_csv(rp.RESULTS / "folds_raw.csv")
    old = old[old.scheme == "leave_one_cohort_out"].rename(columns={"auc": "auc_on_file", "fold": "held_out"})
    rows = [tf.run_loco(Xv, yv, cv, cohort, m, tf.BALANCED) for m in (LR, GBM)]
    new = pd.concat(rows).merge(old[["model", "held_out", "auc_on_file"]], on=["model", "held_out"])
    new["abs_diff"] = (new.auc - new.auc_on_file).abs()
    return len(new) == 10 and bool((new.abs_diff <= GATE_TOL).all()), new


def balance_only(Xv, yv, cv, cohort, n_seeds: int, g0, t0) -> int:
    """Card the 10-seed Arm B run (LABBOOK PLAN 'Arm B with 10 LightGBM seeds'): same thresholds, new files only."""
    if n_seeds < 2:
        raise SystemExit("--balance-seeds needs at least 2 for a spread")
    seeds = tuple(range(n_seeds))
    bal = tf.run_balance(Xv, yv, cv, cohort, seeds)
    summ = tf.balance_summary(bal, GBM)
    verdict = tv.r_balance({m: summ[m] for m in (LR, GBM)})
    meta = dict(date=date.today().isoformat(), plan_commit="c146725", card="the 10-seed Arm B run", model_seeds=list(seeds),
                gate0_max_abs_diff=float(g0.abs_diff.max()), balance=summ, verdict=verdict,
                thresholds=dict(matters=tv.BAL_MATTERS, not_the_drop=tv.BAL_NOT), runtime_s=round(time.time() - t0, 1))
    bal.to_csv(OUT10, index=False)
    META10.write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(json.dumps({"verdict": verdict, **summ}, indent=2))
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gate0-only", action="store_true")
    ap.add_argument("--cache", type=Path)
    ap.add_argument("--balance-seeds", type=int, default=len(SEEDS), metavar="N",
                    help="LightGBM model seeds 0..N-1 for Arm B; default 3 = the the transfer-matrix run run; "
                         "anything else runs Arm B alone and writes the *_10seeds files")
    args = ap.parse_args(argv)
    t0 = time.time()
    cfg = rp.load_config()
    fam, age = tf.cohort_tags(cfg)
    X, y, cohort = load_matrix(args.cache)
    cv, yv = cohort.to_numpy(), y.to_numpy().astype(int)
    Xv = X.to_numpy(dtype=float)
    print(f"matrix {X.shape}, cohorts {sorted(set(cv))}", flush=True)

    ok, g0 = gate0(Xv, yv, cv, cohort)
    print(g0.to_string(index=False))
    print(f"GATE 0: {'PASS' if ok else 'FAIL'} (max abs diff {g0.abs_diff.max():.2e}, tol {GATE_TOL})", flush=True)
    if not ok or args.gate0_only:
        return 0 if ok else 2

    if args.balance_seeds != len(SEEDS):
        return balance_only(Xv, yv, cv, cohort, args.balance_seeds, g0, t0)

    pairs = tf.ordered_pairs(cv)
    assert tf.cell_counts(pairs, fam) == {"same": 8, "cross": 12}, "platform cell counts"
    Xw = rank_normalize_within_sample(X).to_numpy(dtype=float)
    parts = []
    for tag, mat, model, seed in ([("raw", Xv, LR, 0), ("within_sample", Xw, LR, 0)]
                                  + [("raw", Xv, GBM, s) for s in SEEDS]):
        d = tf.run_pairs(mat, yv, cv, model, tf.BALANCED, seed)
        parts.append(d.assign(transform=tag, model=model, model_seed=seed))
        print(f"  arm A {tag} {model} seed {seed}: done", flush=True)
    matrix = pd.concat(parts, ignore_index=True)
    matrix["same_family"] = [tf.pair_cell(a, b, fam) for a, b in zip(matrix.train, matrix.test)]
    matrix["train_age"], matrix["test_age"] = matrix.train.map(age), matrix.test.map(age)
    size = tf.run_sizematched(Xv, yv, cv)
    print("  size-matched control: done", flush=True)
    bal = tf.run_balance(Xv, yv, cv, cohort, SEEDS)
    print("  arm B: done", flush=True)

    raw = matrix[(matrix["transform"] == "raw") & (matrix.model == LR)]
    ws = matrix[(matrix["transform"] == "within_sample") & (matrix.model == LR)]
    d_main = tf.cell_mean(raw, fam, "same") - tf.cell_mean(raw, fam, "cross")
    d_ctrl = tf.cell_mean(size, fam, "same") - tf.cell_mean(size, fam, "cross")
    d_ws = tf.cell_mean(ws, fam, "same") - tf.cell_mean(ws, fam, "cross")
    if any(pd.isna(v) for v in (d_main, d_ctrl, d_ws)):
        raise AssertionError("NaN in a platform difference: an arm is empty, nothing written")
    testable = tf.age_testable(age)
    a_main = a_ctrl = float("nan")
    if testable:
        a_main = tf.cell_mean(raw, age, "same") - tf.cell_mean(raw, age, "cross")
        a_ctrl = tf.cell_mean(size, age, "same") - tf.cell_mean(size, age, "cross")
    v_plat = tv.r_platform(d_main, d_ctrl)
    v_age = tv.r_age(a_main, a_ctrl, testable)
    deltas = {}
    for m in (LR, GBM):
        b = bal[bal.model == m]
        per = b.groupby(["class_weight", "held_out"]).auc.mean().unstack(0)  # mean over seeds
        deltas[m] = float((per["None"] - per[tf.BALANCED]).mean())
        deltas[m + "_per_cohort"] = {k: float(v) for k, v in (per["None"] - per[tf.BALANCED]).items()}
    verdicts = {"R-platform": v_plat, "R-age": v_age, "R-both": tv.r_both(v_plat, v_age) if testable else "skipped (R-age not testable)",
                "R-scale": tv.r_scale(v_plat, d_main, d_ws),
                "R-balance": tv.r_balance({m: deltas[m] for m in (LR, GBM)})}
    gap = pd.read_csv(rp.RESULTS / "gapdiag_cohorts.csv")
    gap = gap[(gap["transform"] == "raw") & (gap.model == LR)].set_index("cohort")
    sp = float(gap.bacterial_share.rank().corr((gap.rcv_auc_in_cohort - gap.loco_auc).rank()))
    if pd.isna(sp) or len(gap) != 5:
        raise AssertionError("class-mix table empty: nothing written")
    gbm_seeds = matrix[matrix.model == GBM]
    gbm_d = [tf.cell_mean(gbm_seeds[gbm_seeds.model_seed == s], fam, "same")
             - tf.cell_mean(gbm_seeds[gbm_seeds.model_seed == s], fam, "cross") for s in SEEDS]
    meta = dict(date=date.today().isoformat(), plan_commits=["c7202b2", "fb2b99b"], gate0_max_abs_diff=float(g0.abs_diff.max()),
                n_samples=int(X.shape[0]), n_genes=int(X.shape[1]), platform_cells=tf.cell_counts(pairs, fam),
                age_testable=testable, age_tags=age, family_tags=fam,
                platform=dict(main_diff=d_main, sizematched_diff=d_ctrl, within_sample_diff=d_ws,
                              same_mean=tf.cell_mean(raw, fam, "same"), cross_mean=tf.cell_mean(raw, fam, "cross"),
                              lightgbm_diff_by_seed=gbm_d),
                age=dict(main_diff=a_main, sizematched_diff=a_ctrl),
                balance_delta=deltas, verdicts=verdicts,
                class_mix_spearman_bacterial_share_vs_gap_n5=sp, runtime_s=round(time.time() - t0, 1),
                note="GSE40396 has exactly 30 arrays: its 20 size-matched draws are identical")
    rp.RESULTS.mkdir(exist_ok=True)
    matrix.to_csv(OUT["matrix"], index=False)
    size.to_csv(OUT["size"], index=False)
    bal.to_csv(OUT["balance"], index=False)
    META.write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(json.dumps({"verdicts": verdicts, **{k: meta[k] for k in ("platform", "age", "balance_delta")}}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
