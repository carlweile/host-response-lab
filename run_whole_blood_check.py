"""Whole-blood check: does the one PBMC cohort (GSE6269) change the held-out number?

    python run_whole_blood_check.py              # Gate 0, then the 4-cohort run
    python run_whole_blood_check.py --gate0-only

Pre-registered in LABBOOK 2026-11-08 (PLAN: whole-blood check). Rules there,
not here. The gene set is FIXED to the 5-cohort shared set (12,040) so that only
the training data changes; the natural 4-cohort intersection is a secondary,
descriptive number.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from src import harmonize, splits, train  # noqa: E402
import run_pipeline as rp  # noqa: E402

RESULTS = ROOT / "results"
PBMC_COHORT = "GSE6269"
GATE0_TOL = 0.005
LGBM_SEEDS = (0, 1, 2)


def load_per_cohort(cfg: dict) -> tuple[dict, dict]:
    """Same per-cohort processing as run_pipeline.cmd_train (raw scale)."""
    loaded = rp._load_all(cfg)
    rp.lab.refuse_mismatches([r for _, r in loaded.values()])
    per_cohort, label_map = {}, {}
    for acc, (series, result) in loaded.items():
        mapping = rp.parse_geo.parse_platform_annot(
            rp.DATA / "platforms" / f"{series.platform}.annot.gz")
        expr, _ = harmonize.maybe_log2(series.expression)
        per_cohort[acc] = harmonize.to_symbols(expr, mapping)
        label_map[acc] = result.usable
    return per_cohort, label_map


def drop_cohort(X: pd.DataFrame, y: pd.Series, cohort: pd.Series, acc: str):
    keep = cohort != acc
    return X.loc[keep], y.loc[keep], cohort.loc[keep]


def assert_cohort_absent(cohort: pd.Series, acc: str) -> None:
    if (cohort == acc).any():
        raise AssertionError(f"{acc} must be absent from the whole-blood matrix")


def loco_rows(X, y, cohort, model_seed: int = 0, only_held_out=None) -> pd.DataFrame:
    rows = []
    for name, model in train.make_models(model_seed).items():
        for held_out, split in splits.cohort_splits(cohort):
            if only_held_out and held_out not in only_held_out:
                continue
            splits.assert_no_cohort_overlap(cohort, split)
            auc, note = train._fit_score(model, X, y, split)
            rows.append({"model": name, "held_out": held_out, "model_seed": model_seed,
                         "n_train": len(split[0]), "n_test": len(split[1]), "auc": auc})
    return pd.DataFrame(rows)


def gate0(X, y, cohort, on_file: pd.DataFrame) -> float:
    """Re-fit the 5-cohort LOCO rows and compare with results/folds_raw.csv."""
    fresh = loco_rows(X, y, cohort)
    ref = on_file[on_file["scheme"] == "leave_one_cohort_out"]
    merged = fresh.merge(ref, left_on=["model", "held_out"], right_on=["model", "fold"],
                         suffixes=("_new", "_ref"))
    if len(merged) != len(fresh):
        raise AssertionError("Gate 0: on-file rows missing")
    return float((merged["auc_new"] - merged["auc_ref"]).abs().max())


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gate0-only", action="store_true")
    args = ap.parse_args(argv)
    cfg = rp.load_config()

    per_cohort, label_map = load_per_cohort(cfg)
    X5, y5, c5, _ = harmonize.build_matrix(per_cohort, label_map)
    on_file = pd.read_csv(RESULTS / "folds_raw.csv")
    g0 = gate0(X5, y5, c5, on_file)
    print(f"Gate 0: max |diff| vs folds_raw.csv = {g0:.2e} (tol {GATE0_TOL})")
    if g0 > GATE0_TOL:
        print("Gate 0 FAILED, nothing written")
        return 1
    if args.gate0_only:
        return 0

    # Primary: same 12,040 genes, GSE6269 removed from train AND test.
    X4, y4, c4 = drop_cohort(X5, y5, c5, PBMC_COHORT)
    assert_cohort_absent(c4, PBMC_COHORT)
    assert list(X4.columns) == list(X5.columns)
    results = train.run_both_schemes(X4, y4, c4)
    detail = pd.DataFrame([r.__dict__ for r in results])
    detail["gene_set"] = "fixed_5cohort"
    detail["model_seed"] = 0
    extra = [loco_rows(X4, y4, c4, model_seed=s) for s in LGBM_SEEDS[1:]]
    extra = pd.concat(extra, ignore_index=True)
    extra = extra[extra["model"] == "lightgbm"].rename(columns={"held_out": "fold"})
    extra["scheme"] = "leave_one_cohort_out"
    extra["gene_set"] = "fixed_5cohort"
    detail = pd.concat([detail, extra], ignore_index=True)

    # Secondary, descriptive: natural 4-cohort gene intersection, logistic LOCO only.
    pc4 = {k: v for k, v in per_cohort.items() if k != PBMC_COHORT}
    lm4 = {k: v for k, v in label_map.items() if k != PBMC_COHORT}
    Xn, yn, cn, notes_n = harmonize.build_matrix(pc4, lm4)
    nat = loco_rows(Xn, yn, cn)
    nat = nat[nat["model"] == "logistic_l2"].rename(columns={"held_out": "fold"})
    nat["scheme"] = "leave_one_cohort_out"
    nat["gene_set"] = f"natural_4cohort_{Xn.shape[1]}"
    detail = pd.concat([detail, nat], ignore_index=True)

    summary = train.summarize(results)
    wb = [c for c in c5.unique() if c != PBMC_COHORT]
    ref = on_file[(on_file["scheme"] == "leave_one_cohort_out") & on_file["fold"].isin(wb)]
    comp = {}
    for m in ("logistic_l2", "lightgbm"):
        five = ref[ref["model"] == m]["auc"].mean()
        four = detail[(detail["scheme"] == "leave_one_cohort_out") & (detail["model"] == m)
                      & (detail["model_seed"] == 0) & (detail["gene_set"] == "fixed_5cohort")]["auc"].mean()
        comp[m] = {"loco_5cohort_run_same4_heldout": float(five),
                   "loco_4cohort_run": float(four), "delta_train": float(four - five)}
    lg = detail[(detail["scheme"] == "leave_one_cohort_out") & (detail["model"] == "lightgbm")
                & (detail["gene_set"] == "fixed_5cohort")].groupby("model_seed")["auc"].mean()
    comp["lightgbm_loco_by_seed"] = {int(k): float(v) for k, v in lg.items()}
    comp["natural_intersection_genes"] = int(Xn.shape[1])
    comp["natural_logistic_loco_mean"] = float(nat["auc"].mean())

    RESULTS.mkdir(exist_ok=True)
    detail.to_csv(RESULTS / "whole_blood_check.csv", index=False)
    meta = {"date": date.today().isoformat(), "pbmc_cohort_removed": PBMC_COHORT,
            "gate0_max_abs_diff": g0, "n_samples": int(X4.shape[0]), "n_genes_fixed": int(X4.shape[1]),
            "summary_4cohort_fixed_genes": summary.to_dict(orient="records"),
            "comparison": comp, "notes_natural": notes_n}
    (RESULTS / "run_whole_blood_check.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(summary.to_string(index=False))
    print(json.dumps(comp, indent=2))
    print(f"wrote {RESULTS / 'whole_blood_check.csv'} and run_whole_blood_check.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
