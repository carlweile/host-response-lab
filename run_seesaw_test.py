"""Pre-registered seesaw test (LABBOOK PLAN 9daaf08, CORE §2 condition 3).

    python run_seesaw_test.py --parts synthetic         # instrument check, no real data
    python run_seesaw_test.py --parts repro             # reproduction gate vs results/gapdiag_*
    python run_seesaw_test.py --parts probe             # one-fold timing probe
    python run_seesaw_test.py                           # p1 p2 p3 on both cohort sets
    python run_seesaw_test.py --evaluate-only           # re-judge results/seesaw_* files

Every whole-cohort arm is a transductive diagnostic of a random-CV artefact. No
rank-normalised AUC from here is a generalisation number.
"""
from __future__ import annotations

import argparse
import json
import platform
import subprocess
import sys
import time
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

import run_pipeline as rp
from run_gap_diagnostics import load_matrix
from src import seesaw as ss
from src import seesaw_verdict as sv

PLAN_COMMIT = "9daaf08"
SEEDS = [100, 101, 102, 103, 104]
SETS = {"all": None, "no40396": "GSE40396"}
OUT = {k: rp.RESULTS / f"seesaw_{k}.csv" for k in ("folds", "offsets", "summary")}
META = rp.RESULTS / "run_seesaw.json"


def _subset(X, y, cohort, drop):
    keep = cohort != drop if drop else cohort == cohort
    return X[keep], y[keep], cohort[keep]


def _jobs(seeds):
    """(part, arm, model name, model, k, C, seeds) for every judged and report-only arm."""
    jobs = []
    for k in (2, 5, 10):
        C = ss.c_for_k(k)
        jobs += [("p1", arm, "logistic_l2", ss.logistic(C), k, C, seeds) for arm in ("W", "I")]
    jobs += [("p1", arm, "logistic_l2", ss.logistic(0.05), k, 0.05, seeds[:1])
             for k in (2, 10) for arm in ("W", "I")]                       # fixed-C, report
    jobs.append(("p1", "raw", "logistic_l2", ss.logistic(0.05), 5, 0.05, seeds))  # report
    jobs += [("p2", arm, "logistic_l2", ss.logistic(0.05), 5, 0.05, seeds)
             for arm in ss.WHOLE[2:] + ss.TRAIN[1:]]
    jobs += [("p3", arm, name, model, 5, ss.c_for_k(5), seeds)
             for name, model in (("lightgbm", ss.lgbm()), ("lgbm_stumps", ss.lgbm(True)))
             for arm in ("W", "I")]
    return jobs


def run_parts(X, y, cohort, parts, seeds, cpu: dict, first_fold=False, sets=SETS):
    folds, offsets = [], []
    for name, drop in sets.items():
        Xs, ys, cs = _subset(X, y, cohort, drop)
        for part, arm, mname, model, k, C, sd in _jobs(seeds):
            if part not in parts:
                continue
            t0, w0 = time.process_time(), time.time()
            for seed in sd:
                f, o = ss.run_arm(Xs, ys, cs, arm, model, k, seed,
                                  dict(set=name, part=part, model=mname, C=C), first_fold)
                folds += f
                offsets += o
            cpu[f"{name}/{part}/{mname}/{arm}/k{k}/C{C:.4f}"] = dict(
                cpu_s=round(time.process_time() - t0, 1), wall_s=round(time.time() - w0, 1),
                k=k, n_seeds=len(sd))
            print(f"  {name} {part} {mname} {arm} k={k} C={C:.4f}: "
                  f"{cpu[list(cpu)[-1]]['wall_s']} s", flush=True)
    return pd.DataFrame(folds), pd.DataFrame(offsets)


def synthetic_check() -> dict:
    """Positive control (shares 20/20/80/80) must show the effect, negative (all 50 %)
    must not. Uses the same runner and the same aggregation as the real run."""
    from tests.test_gap_diagnostics import _synthetic
    res = {}
    for label, shares in (("positive", [.2, .2, .8, .8]), ("negative", [.5] * 4)):
        X, y, c = _synthetic(shares)
        f, o = run_parts(X, y, c, {"p1", "p2"}, SEEDS, {}, sets={"all": None})
        res[label] = sv.quantities(sv.pair_weighted(f), o, {})["all"]
    pos, neg = res["positive"], res["negative"]
    ok_pos = pos["X10"] > pos["X2"] > 0 and pos["Yc"] > 0
    lim = 0.25 * pos["X5"]
    ok_neg = all(abs(neg[key]) < lim for key in ("X2", "X5", "X10", "Yc"))
    return dict(positive={k: pos[k] for k in ("X2", "X5", "X10", "Yc")},
                negative={k: neg[k] for k in ("X2", "X5", "X10", "Yc")},
                positive_ok=bool(ok_pos), negative_ok=bool(ok_neg))


def repro_gate(X, y, cohort) -> dict:
    """k=5, fold seed 0, C=0.05 must reproduce the logistic cohort_rank rows of
    gapdiag_folds.csv (W) and gapdiag_inductive.csv (I), auc_within and auc_cross."""
    old_w = pd.read_csv(rp.RESULTS / "gapdiag_folds.csv")
    old_w = old_w[(old_w["transform"] == "cohort_rank") & (old_w.model == "logistic_l2")
                  & (old_w.scheme == "random_cv")]
    old_i = pd.read_csv(rp.RESULTS / "gapdiag_inductive.csv")
    old_i = old_i[old_i.model == "logistic_l2"]
    worst = 0.0
    for arm, old in (("W", old_w), ("I", old_i)):
        new, _ = ss.run_arm(X, y, cohort, arm, ss.logistic(0.05), 5, 0, {})
        new = pd.DataFrame(new).merge(old, on="fold", suffixes=("", "_old"))
        if len(new) != 5:
            raise SystemExit(f"repro gate: {arm} arm matched {len(new)} of 5 folds")
        for col in ("auc_within", "auc_cross"):
            worst = max(worst, float((new[col] - new[col + "_old"]).abs().max()))
    return dict(max_abs_diff=worst, ok=worst <= 1e-9)


def _commit() -> str:
    head = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True,
                          text=True, cwd=rp.ROOT).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "--", "src", "run_seesaw_test.py"],
                           capture_output=True, text=True, cwd=rp.ROOT).stdout.strip()
    return head + ("-dirty" if dirty else "")


def evaluate(meta: dict) -> dict:
    folds, offsets = pd.read_csv(OUT["folds"]), pd.read_csv(OUT["offsets"])
    summary = sv.pair_weighted(folds)
    summary.to_csv(OUT["summary"], index=False)
    q = sv.quantities(summary, offsets, meta.get("treatment_gate", {}))
    verdicts = {cs: sv.judge(v) for cs, v in q.items()}
    return dict(quantities=q, verdicts=verdicts, overall=sv.overall(verdicts))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--parts", nargs="+", default=["p1", "p2", "p3"],
                    choices=["synthetic", "repro", "probe", "p1", "p2", "p3"])
    ap.add_argument("--seeds", nargs="+", type=int, default=SEEDS)
    ap.add_argument("--evaluate-only", action="store_true")
    ap.add_argument("--cache", type=Path, help="pickle file to reuse the parsed matrix")
    args = ap.parse_args(argv)

    meta = json.loads(META.read_text(encoding="utf-8")) if META.exists() else {}
    meta.update(plan_commit=PLAN_COMMIT, code_commit=_commit(), python=platform.python_version(),
                lgbm_threads=ss.LGBM_THREADS)
    if not args.evaluate_only:
        if "synthetic" in args.parts:
            meta["synthetic"] = dict(synthetic_check(), date=date.today().isoformat())
            print(json.dumps(meta["synthetic"], indent=2))
        real = [p for p in args.parts if p != "synthetic"]
        if real:
            X, y, cohort = load_matrix(args.cache)
            print(f"matrix {X.shape}, {int(y.sum())} bacterial, {cohort.nunique()} cohorts")
            meta["treatment_gate"] = {n: ss.treatment_gate(*_subset(X, y, cohort, d)[::2])
                                      for n, d in SETS.items()}
        if "repro" in real:
            meta["repro_gate"] = repro_gate(X, y, cohort)
            print(meta["repro_gate"])
        if "probe" in real:
            cpu: dict = {}
            run_parts(X, y, cohort, {"p1", "p2", "p3"}, SEEDS[:1], cpu, first_fold=True)
            full = {key: v["k"] * (len(args.seeds) if v["n_seeds"] > 1 else 1)
                    for key, v in cpu.items()}
            meta["probe"] = dict(per_job=cpu, projected_wall_min=round(
                sum(cpu[key]["wall_s"] * n for key, n in full.items()) / 60, 1),
                projected_cpu_min=round(sum(cpu[key]["cpu_s"] * n for key, n in full.items())
                                        / 60, 1))
            print(json.dumps({k: v for k, v in meta["probe"].items() if k != "per_job"}))
        todo = [p for p in real if p in ("p1", "p2", "p3")]  # not "probe" (review 2026-10-21)
        if todo:
            meta["cpu"] = {}
            f, o = run_parts(X, y, cohort, set(todo), args.seeds, meta["cpu"])
            f.to_csv(OUT["folds"], index=False)
            o.to_csv(OUT["offsets"], index=False)
            meta.update(seeds=args.seeds, run_date=date.today().isoformat())
    if OUT["folds"].exists() and (args.evaluate_only or "p1" in args.parts):
        meta["result"] = evaluate(meta)
        print(json.dumps(meta["result"], indent=2, default=str))
    META.write_text(json.dumps(meta, indent=2, default=str), encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
