"""Take the negative logistic gap apart (LABBOOK 2026-10-19 PLAN, CORE §2 condition 3).

    python run_gap_diagnostics.py                    # E0-E4 + mirror (CPU, ~20-60 min)
    python run_gap_diagnostics.py --parts e1 e3 e4   # some parts
    python run_gap_diagnostics.py --evaluate-only    # re-judge results/gapdiag_* files
    python run_gap_diagnostics.py --cache FILE       # reuse a pickled matrix

Reads data/ the same way as `run_pipeline.py --train`. The thresholds in
evaluate() are copied from the PLAN entry committed before any run (cc7ba64);
they must not be tuned after looking at results.
"""
from __future__ import annotations

import argparse
import json
import math
import platform
import sys
import time
from datetime import date
from pathlib import Path

import pandas as pd

import run_pipeline as rp
from src import gap_experiments as gx
from src import harmonize, parse_geo
from src import labels as lab

PLAN_COMMIT = "cc7ba64"
FILE_TAG = {"raw": "raw", "cohort_rank": "ranknorm", "within_sample": "withinsample"}
OUT = {k: rp.RESULTS / f"gapdiag_{k}.csv"
       for k in ("folds", "cohorts", "inductive", "splithalf", "balance_folds", "balance_summary",
                 "mirror")}
META = rp.RESULTS / "run_gapdiag.json"
LR, GBM = "logistic_l2", "lightgbm"


def load_matrix(cache: Path | None, collapse: str = "mean"):
    """Same steps and skip rules as run_pipeline.cmd_train; E0 checks it matches."""
    if cache is not None and cache.exists():
        return pd.read_pickle(cache)
    per_cohort, label_map = {}, {}
    loaded = rp._load_all(rp.load_config())
    lab.refuse_mismatches([result for _, result in loaded.values()])
    for acc, (series, result) in loaded.items():
        annot = rp.DATA / "platforms" / f"{series.platform}.annot.gz"
        if len(result.usable) < 20 or not annot.exists():
            continue
        expr, _ = harmonize.maybe_log2(series.expression)
        genes = harmonize.to_symbols(expr, parse_geo.parse_platform_annot(annot), collapse)
        if not genes.empty:
            per_cohort[acc], label_map[acc] = genes, result.usable
    X, y, cohort, _ = harmonize.build_matrix(per_cohort, label_map)
    if cache is not None:
        pd.to_pickle((X, y, cohort), cache)
    return X, y, cohort


def _with_file_auc(folds: pd.DataFrame, transform: str) -> pd.DataFrame:
    """E0: put the AUC on file (results/folds_<tag>.csv) next to the fresh one."""
    old = pd.read_csv(rp.RESULTS / f"folds_{FILE_TAG[transform]}.csv")
    old = old.rename(columns={"auc": "auc_on_file"})[["scheme", "model", "fold", "auc_on_file"]]
    out = folds.merge(old, on=["scheme", "model", "fold"], how="left")
    out["abs_diff_vs_file"] = (out["auc_pooled"] - out["auc_on_file"]).abs()
    return out


def run_e1(X, y, cohort) -> None:
    folds, cohorts = [], []
    for t in gx.TRANSFORMS:
        f, c = gx.both_schemes_with_scores(X, y, cohort, t)
        folds.append(_with_file_auc(f, t))
        cohorts.append(c)
        print(f"  E1 {t}: done", flush=True)
    pd.concat(folds).to_csv(OUT["folds"], index=False)
    pd.concat(cohorts).to_csv(OUT["cohorts"], index=False)


def run_e3(X, y, cohort) -> None:
    gx.inductive_random_cv(X, y, cohort).to_csv(OUT["inductive"], index=False)


def run_e4(X, y, cohort) -> None:
    gx.split_half_heldout(X, y, cohort).to_csv(OUT["splithalf"], index=False)


def _run_e2(design: str):
    def runner(X, y, cohort) -> None:
        new = gx.balance_control(X, y, cohort, design)
        old = pd.read_csv(OUT["balance_folds"]) if OUT["balance_folds"].exists() else None
        if old is not None:
            new = pd.concat([old[old["design"] != design], new], ignore_index=True)
        new.to_csv(OUT["balance_folds"], index=False)
        gx.summarize_balance(new).to_csv(OUT["balance_summary"], index=False)
    return runner


def run_mirror(X, y, cohort) -> None:
    """Exploratory, added after seeing E1/E3; no prediction is judged on it."""
    gx.score_mirror(X, y, cohort).to_csv(OUT["mirror"], index=False)


RUNNERS = {"e1": run_e1, "e3": run_e3, "e4": run_e4, "e2a": _run_e2("balanced"),
           "e2b": _run_e2("size_matched"), "mirror": run_mirror}


def _band(value: float, supported: bool, failed: bool, fail_word: str = "KILLED") -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "not measured"
    return "supported" if supported else fail_word if failed else "partial"


def _read(key: str) -> pd.DataFrame | None:
    return pd.read_csv(OUT[key]) if OUT[key].exists() else None


def evaluate() -> dict:
    """Judge each pre-registered prediction (P1-P5) from the files on disk."""
    v: dict = {}
    f, c, ind, sh, bal = (_read(k) for k in
                          ("folds", "cohorts", "inductive", "splithalf", "balance_summary"))
    if f is not None:
        v["E0"] = {f"{t}/{m}": {"max_abs_diff": float(g["abs_diff_vs_file"].max()),
                               "tolerance": 0.005 if m == LR else 0.02}
                   for (t, m), g in f.groupby(["transform", "model"])}
        for row in v["E0"].values():
            row["verdict"] = "reproduced" if row["max_abs_diff"] <= row["tolerance"] else "DRIFT"
        rcv = f[(f["transform"] == "cohort_rank") & (f["scheme"] == "random_cv")]
        lr, gbm = rcv[rcv["model"] == LR], rcv[rcv["model"] == GBM]
        w = float(lr["auc_within"].mean())
        v["P1"] = {"within": w, "cross": float(lr["auc_cross"].mean()),
                   "pooled": float(lr["auc_pooled"].mean()),
                   "verdict": _band(w, w >= 0.88, w < 0.85)}
        d = float(gbm["auc_within"].mean() - gbm["auc_pooled"].mean())
        v["LightGBM_direction"] = {"within_minus_pooled": d,
                                   "verdict": "as predicted" if d >= 0 else "NOT as predicted"}
        v["iteration_cap"] = {"logistic_fits_at_3000": int(((f["model"] == LR)
                                                           & (f["n_iter"] >= 3000)).sum())}
        if ind is not None:
            ref, new = float(lr["auc_pooled"].mean()), float(
                ind.loc[ind["model"] == LR, "auc_pooled"].mean())
            v["P4"] = {"inductive_pooled": new, "transductive_pooled": ref,
                       "verdict": _band(new, abs(new - ref) <= 0.03, new >= 0.85)}
    if c is not None:
        lr = c[(c["transform"] == "cohort_rank") & (c["model"] == LR)]
        d = float((lr["rcv_auc_in_cohort"] - lr["loco_auc"]).mean())
        v["P1b"] = {"in_cohort_rcv_minus_loco": d, "verdict": _band(d, d >= -0.02, d <= -0.05)}
        v["share_vs_mean_logit_spearman"] = {
            f"{t}/{m}": {cls: float(g["bacterial_share"].corr(g[f"mean_logit_{cls}"],
                                                              method="spearman"))
                         for cls in ("bacterial", "viral")}
            for (t, m), g in c.groupby(["transform", "model"])}
    if bal is not None:
        def gap(design: str) -> float:
            row = bal[(bal["design"] == design) & (bal["transform"] == "cohort_rank")
                      & (bal["model"] == LR)]
            return float(row["gap_mean"].iloc[0]) if len(row) else float("nan")
        g2, g3 = gap("balanced"), gap("size_matched")
        v["P2"] = {"gap_mean": g2, "verdict": _band(g2, g2 >= -0.02, g2 <= -0.05)}
        v["P3"] = {"gap_mean": g3,
                   "verdict": _band(g3, g3 <= -0.05, g3 >= -0.02, "NOT SUPPORTED")}
    if sh is not None:
        per_cohort = sh[sh["model"] == LR].groupby("cohort")["diff_whole_minus_half_a"].mean()
        d = float(per_cohort.mean())
        v["P5"] = {"mean_diff": d, "per_cohort": per_cohort.round(4).to_dict(),
                   "verdict": _band(d, d <= 0.02, d > 0.02, "NOT SUPPORTED")}
    return v


def _clean(obj):
    if isinstance(obj, dict):
        return {k: _clean(x) for k, x in obj.items()}
    if isinstance(obj, float) and math.isnan(obj):
        return None
    return obj


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--parts", nargs="+", choices=list(RUNNERS),
                    default=["e1", "e3", "e4", "e2a", "e2b", "mirror"])
    ap.add_argument("--evaluate-only", action="store_true")
    ap.add_argument("--cache", type=Path, help="pickle file to reuse the parsed matrix")
    args = ap.parse_args(argv)

    meta = json.loads(META.read_text(encoding="utf-8")) if META.exists() else {}
    if not args.evaluate_only:
        X, y, cohort = load_matrix(args.cache)
        on_file = json.loads((rp.RESULTS / "run_raw.json").read_text(encoding="utf-8"))
        if X.shape != (on_file["n_samples"], on_file["n_genes"]):
            raise SystemExit(f"matrix {X.shape} differs from run_raw.json; stop")
        print(f"matrix {X.shape}, {int(y.sum())} bacterial, {cohort.nunique()} cohorts",
              flush=True)
        for part in args.parts:
            t0 = time.time()
            RUNNERS[part](X, y, cohort)
            secs = round(time.time() - t0, 1)
            meta.setdefault("parts", {})[part] = {"date": date.today().isoformat(),
                                                  "seconds": secs}
            print(f"{part}: {secs} s", flush=True)
    import lightgbm, numpy, sklearn  # noqa: E401  (versions for the record)
    meta.update(plan_commit=PLAN_COMMIT, python=platform.python_version(),
                versions={"numpy": numpy.__version__, "pandas": pd.__version__,
                          "scikit-learn": sklearn.__version__,
                          "lightgbm": lightgbm.__version__},
                verdicts=_clean(evaluate()))
    META.write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(json.dumps(meta["verdicts"], indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
