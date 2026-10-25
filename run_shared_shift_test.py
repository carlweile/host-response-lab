"""Pre-registered shared-shift test (LABBOOK PLAN 2fd8865, CORE §2 condition 3).

    python run_shared_shift_test.py --parts synthetic   # instrument gate, 60 + 2,000 genes
    python run_shared_shift_test.py --parts repro       # reproduction gate vs results/gapdiag_*
    python run_shared_shift_test.py --parts run         # H1-H3 on both cohort sets
    python run_shared_shift_test.py --evaluate-only     # re-judge results/sharedshift_*

Every whole-cohort arm is a transductive diagnostic of a random-CV artefact.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from datetime import date
from pathlib import Path

import pandas as pd

import run_pipeline as rp
from run_gap_diagnostics import load_matrix
from run_seesaw_test import SEEDS, SETS, _commit, _subset
from src import seesaw as ss
from src import seesaw_verdict as sv
from src import shared_shift as sh
from src import shared_shift_verdict as hv

PLAN_COMMIT = "2fd8865"
OUT = {k: rp.RESULTS / f"sharedshift_{k}.csv" for k in ("folds", "offsets", "summary")}
META = rp.RESULTS / "run_sharedshift.json"


def run_all(X, y, cohort, seeds, sets, cpu: dict):
    folds, offsets = [], []
    for name, drop in sets.items():
        Xs, ys, cs = _subset(X, y, cohort, drop)
        t0, w0 = time.process_time(), time.time()
        for seed in seeds:
            f, o = sh.run_k5(Xs, ys, cs, seed, dict(set=name, part="k5", model="logistic_l2"))
            folds += f
            offsets += o
            for k in (2, 10):
                for arm in ("W", "I"):
                    f, _ = ss.run_arm(Xs, ys, cs, arm, ss.logistic(ss.c_for_k(k)), k, seed,
                                      dict(set=name, part="h3", model="logistic_l2",
                                           C=ss.c_for_k(k)))
                    folds += f
        cpu[name] = dict(cpu_s=round(time.process_time() - t0, 1),
                         wall_s=round(time.time() - w0, 1))
        print(f"  {name}: {cpu[name]['wall_s']} s", flush=True)
    return pd.DataFrame(folds), pd.DataFrame(offsets)


def judge_frames(folds, offsets):
    summary = sv.pair_weighted(folds)
    q = hv.quantities(summary, offsets)
    verdicts = {cs: hv.judge(v) for cs, v in q.items()}
    return summary, dict(quantities=q, verdicts=verdicts, overall=hv.overall(verdicts))


def synthetic_gate() -> dict:
    """Positive control must give H1 PASS and H2 PASS at 60 and 2,000 genes; negative
    control must stay below 0.25 x the positive |D_W - D_raw| and |D_I - D_raw|."""
    from tests.test_gap_diagnostics import _synthetic
    out = {}
    for genes in (60, 2000):
        res = {}
        for label, shares in (("positive", [.2, .2, .8, .8]), ("negative", [.5] * 4)):
            f, o = run_all(*_synthetic(shares, genes=genes), SEEDS, {"all": None}, {})
            res[label] = judge_frames(f, o)[1]
        pos, neg = res["positive"]["quantities"]["all"], res["negative"]["quantities"]["all"]
        v = res["positive"]["verdicts"]["all"]
        lim = {a: 0.25 * abs(pos[f"D_{a}"] - pos["D_raw"]) for a in ("W", "I")}
        out[f"genes_{genes}"] = dict(
            positive=pos, negative=neg, positive_verdicts=v,
            positive_ok=v["H1"] == "PASS" and v["H2"] == "PASS",
            negative_ok=all(abs(neg[f"D_{a}"] - neg["D_raw"]) < lim[a] for a in ("W", "I")))
    out["ok"] = all(g["positive_ok"] and g["negative_ok"] for g in out.values()
                    if isinstance(g, dict))
    return out


def repro_gate(X, y, cohort) -> dict:
    """run_k5 at fold seed 0 must reproduce the logistic W and I rows on file to 1e-9."""
    new = pd.DataFrame(sh.run_k5(X, y, cohort, 0, {})[0])
    old_w = pd.read_csv(rp.RESULTS / "gapdiag_folds.csv")
    old_w = old_w[(old_w["transform"] == "cohort_rank") & (old_w.model == "logistic_l2")
                  & (old_w.scheme == "random_cv")]
    old_i = pd.read_csv(rp.RESULTS / "gapdiag_inductive.csv")
    old_i = old_i[old_i.model == "logistic_l2"]
    worst = 0.0
    for arm, old in (("W", old_w), ("I", old_i)):
        m = new[new.arm == arm].merge(old, on="fold", suffixes=("", "_old"))
        if len(m) != 5:
            return dict(ok=False, error=f"{arm} matched {len(m)} of 5 folds")
        for col in ("auc_within", "auc_cross"):
            worst = max(worst, float((m[col] - m[col + "_old"]).abs().max()))
    return dict(max_abs_diff=worst, ok=worst <= 1e-9)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--parts", nargs="+", default=["run"], choices=["synthetic", "repro", "run"])
    ap.add_argument("--evaluate-only", action="store_true")
    ap.add_argument("--cache", type=Path, help="pickle file to reuse the parsed matrix")
    args = ap.parse_args(argv)

    meta = json.loads(META.read_text(encoding="utf-8")) if META.exists() else {}
    meta.update(plan_commit=PLAN_COMMIT, code_commit=_commit())
    if not args.evaluate_only:
        if "synthetic" in args.parts:
            meta["synthetic"] = dict(synthetic_gate(), date=date.today().isoformat())
            print(json.dumps({k: (v if not isinstance(v, dict) else
                                  {x: v[x] for x in ("positive_verdicts", "positive_ok",
                                                     "negative_ok")})
                              for k, v in meta["synthetic"].items()}, indent=2, default=str))
        if "repro" in args.parts or "run" in args.parts:
            X, y, cohort = load_matrix(args.cache)
            print(f"matrix {X.shape}, {int(y.sum())} bacterial, {cohort.nunique()} cohorts")
        if "repro" in args.parts:
            meta["repro_gate"] = repro_gate(X, y, cohort)
            print(meta["repro_gate"])
        if "run" in args.parts:
            gates = meta.get("synthetic", {}).get("ok") and meta.get("repro_gate", {}).get("ok")
            if not gates:
                meta["run_refused"] = "gates not passed (synthetic + repro)"
                print(meta["run_refused"])
                args.parts.remove("run")
        if "run" in args.parts:
            meta["cpu"] = {}
            f, o = run_all(X, y, cohort, SEEDS, SETS, meta["cpu"])
            f.to_csv(OUT["folds"], index=False)
            o.to_csv(OUT["offsets"], index=False)
            meta["run_date"] = date.today().isoformat()
    if OUT["folds"].exists() and (args.evaluate_only or "run" in args.parts):
        summary, meta["result"] = judge_frames(pd.read_csv(OUT["folds"]),
                                               pd.read_csv(OUT["offsets"]))
        summary.to_csv(OUT["summary"], index=False)
        print(json.dumps(meta["result"], indent=2, default=str))
    META.write_text(json.dumps(meta, indent=2, default=str), encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
