"""v1.1: frozen published signatures, random CV vs leave-one-cohort-out (LABBOOK PLAN 2026-10-21).

    python run_frozen_signatures.py

No training. Variants: main (log2, highest-mean probe), linear (2^x), firstprobe
(lowest probe id, reads no expression value), ratio_one_term (Sweeney's gene-count
ratio on the bacterial term only), metaintegrator (calculateScore: shift if min < 0,
z-score per cohort; reads the whole cohort). Fold seeds 0 and 200-209.
Comparison: our trained models' LOCO AUC on the same cohorts (results/folds_raw.csv).
"""
from __future__ import annotations

import json
import sys
from datetime import date

import numpy as np
import pandas as pd

import run_pipeline as rp
from src import frozen_signatures as fs, harmonize, labels as lab, parse_geo

SEEDS = list(range(200, 210))
# variant -> (probe collapse, linear 2^x, gene-count ratio on bacterial term, MetaIntegrator steps)
VARIANTS = {"main": ("mean", False, False, False), "linear": ("mean", True, False, False),
            "firstprobe": ("first", False, False, False),
            "ratio_one_term": ("mean", False, True, False),  # PLAN 661ea85
            "metaintegrator": ("mean", False, False, True)}  # PLAN 661ea85
KILL_AUC = 0.75  # PLAN: below this on a cohort the paper itself reports working on = bug, stop
KILL_COHORTS = {"herberg2016_drs": ["GSE40396", "GSE60244"],
                "sweeney2016_bv7": ["GSE6269", "GSE40396", "GSE42026"]}
OUT_COHORTS = rp.RESULTS / "frozen_signatures_cohorts.csv"
OUT_SUMMARY = rp.RESULTS / "frozen_signatures_summary.csv"
META = rp.RESULTS / "run_frozen_signatures.json"


def load_cohorts(collapse: str) -> dict[str, tuple[pd.DataFrame, dict[str, str]]]:
    """Per cohort: genes x LABELLED samples (log2, same steps as run_pipeline.cmd_train)."""
    loaded = rp._load_all(rp.load_config())
    lab.refuse_mismatches([result for _, result in loaded.values()])
    out = {}
    for acc, (series, result) in loaded.items():
        annot = rp.DATA / "platforms" / f"{series.platform}.annot.gz"
        expr, _ = harmonize.maybe_log2(series.expression)
        genes = harmonize.to_symbols(expr, parse_geo.parse_platform_annot(annot), collapse)
        samples = [s for s in genes.columns if s in result.usable]
        out[acc] = (genes[samples], {s: result.usable[s] for s in samples})
    return out


def pooled_scores(cohorts, sig: fs.Signature, linear: bool, ratio: bool = False, mi: bool = False):
    """Score every cohort that carries all genes; return s, y, cohort, skipped list."""
    w = fs.ratio_weight(sig) if ratio else 1.0
    parts, ys, groups, skipped = [], [], [], []
    for acc, (genes, labels) in cohorts.items():
        if any(g not in genes.index for g in sig.genes):
            skipped.append(acc)
            continue
        m = 2.0 ** genes if linear else genes
        if mi:
            m = fs.metaintegrator_shift(m)
        fs.assert_per_sample(lambda x, g: fs.score(x, g, w), m, sig)
        s = fs.score(m, sig, w)
        s.index = [f"{acc}:{c}" for c in s.index]
        parts.append(s)
        ys.extend(1 if labels[c] == lab.BACTERIAL else 0 for c in genes.columns)
        groups.extend([acc] * len(s))
    s = pd.concat(parts)
    cohort = pd.Series(groups, index=s.index)
    if mi:  # after per-patient scoring: the one cohort-wide step
        s = fs.metaintegrator_zscore(s, cohort)
    return s, pd.Series(ys, index=s.index), cohort, skipped


def trained_loco() -> pd.DataFrame:
    f = pd.read_csv(rp.RESULTS / "folds_raw.csv")
    f = f[f["scheme"] == "leave_one_cohort_out"]
    return f.pivot(index="fold", columns="model", values="auc")


def main() -> int:
    sigs = fs.load_signatures(rp.ROOT / "config" / "signatures.json")
    ours = trained_loco()
    cohort_rows, summary_rows, kill = [], [], []
    n_by_cohort = {}
    for variant, (collapse, linear, ratio, mi) in VARIANTS.items():
        cohorts = load_cohorts(collapse)
        n_by_cohort = {a: len(l) for a, (_, l) in cohorts.items()}
        for sig in sigs:
            s, y, cohort, skipped = pooled_scores(cohorts, sig, linear, ratio, mi)
            seed0 = fs.both_schemes(s, y, cohort, seed=0)
            sweep = [fs.both_schemes(s, y, cohort, seed=k) for k in SEEDS]
            gaps = np.array([r["gap"] for r in sweep])
            rcvs = np.array([r["random_cv"] for r in sweep])
            per = seed0["per_cohort"]
            for acc, auc in per.items():
                role = "independent" if acc in sig.independent else "discovery" if acc in sig.discovery else "?"
                cohort_rows.append(dict(signature=sig.name, variant=variant, cohort=acc, role=role,
                                        n=int((cohort == acc).sum()), n_bacterial=int(y[cohort == acc].sum()),
                                        auc=auc, logistic_loco=ours.loc[acc, "logistic_l2"],
                                        lightgbm_loco=ours.loc[acc, "lightgbm"]))
                if variant == "main" and acc in KILL_COHORTS[sig.name] and auc < KILL_AUC:
                    kill.append(f"{sig.name} {acc} AUC {auc:.3f} < {KILL_AUC}")
            ind = [per[a] for a in per if a in sig.independent]
            dis = [per[a] for a in per if a in sig.discovery]
            covered = list(per)
            summary_rows.append(dict(
                signature=sig.name, variant=variant, cohorts=len(covered), n=len(s),
                skipped_missing_genes=";".join(skipped),
                random_cv_seed0=seed0["random_cv"], loco=seed0["loco"], gap_seed0=seed0["gap"],
                random_cv_sweep_mean=rcvs.mean(), gap_sweep_mean=gaps.mean(), gap_sweep_sd=gaps.std(ddof=1),
                gap_sweep_min=gaps.min(), gap_sweep_max=gaps.max(),
                loco_independent=np.mean(ind) if ind else np.nan,
                loco_discovery=np.mean(dis) if dis else np.nan,
                logistic_loco_same_cohorts=ours.loc[covered, "logistic_l2"].mean(),
                lightgbm_loco_same_cohorts=ours.loc[covered, "lightgbm"].mean()))

    cdf, sdf = pd.DataFrame(cohort_rows), pd.DataFrame(summary_rows)
    pd.set_option("display.width", 220)
    print(cdf.round(3).to_string(index=False))
    print("\n" + sdf.drop(columns=["skipped_missing_genes"]).round(3).to_string(index=False))
    print("\nKILL RULE:", "FIRED -> " + "; ".join(kill) if kill else "not fired")
    cdf.to_csv(OUT_COHORTS, index=False)
    sdf.to_csv(OUT_SUMMARY, index=False)
    META.write_text(json.dumps({
        "date": date.today().isoformat(), "plan": "LABBOOK 2026-10-21 PLAN v1.1 (c193407) + supplements PLAN (661ea85)",
        "seeds": [0] + SEEDS, "variants": {k: {"collapse": c, "linear": l, "ratio_one_term": r, "metaintegrator": m}
                     for k, (c, l, r, m) in VARIANTS.items()},
        "labelled_per_cohort": n_by_cohort, "kill_auc": KILL_AUC, "kill_fired": kill,
        "comparison": "results/folds_raw.csv LOCO, fold seed 0, model seed 0",
    }, indent=2), encoding="utf-8")
    print(f"wrote {OUT_COHORTS.name}, {OUT_SUMMARY.name}, {META.name}")
    return 1 if kill else 0


if __name__ == "__main__":
    sys.exit(main())
