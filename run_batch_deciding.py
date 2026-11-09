"""Deciding test of the batch cause (LABBOOK 2026-11-08 PLAN, the deciding batch test).

    python run_batch_deciding.py                # Gates -> Arm A (1000-draw null) -> Arm B power gate -> verdict
    python run_batch_deciding.py --gates-only   # Gates 0a/0b/0c + Gate B, writes nothing
    python run_batch_deciding.py --cache FILE   # reuse a pickled matrix

Everything is held in memory and written only at the end. Any gate miss = nothing written.
The rule in src/batch_deciding_verdict.py was committed in the PLAN entry (0e529bd) before this ran.
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
from run_batch_split import GATE0_TOL, LR, gate0
from run_gap_diagnostics import load_matrix
from src import batch_deciding as bd
from src import batch_deciding_verdict as v
from src import batch_split as bs
from src import labels as lab
from src import parse_geo
from src import transfer as tf
from src.splits import rank_normalize_within_sample

PLAN_COMMIT = "0e529bd"
DATES_ACC = "GSE40396"
OUT_NULL = rp.RESULTS / "batch_null1000.csv"
OUT_DATES = rp.RESULTS / "batch_dates_gse40396.csv"
META = rp.RESULTS / "run_batch_deciding.json"


def arm_b_counts() -> tuple[dict, dict, int]:
    """Labelled bacterial/viral arrays per hybridization date in GSE40396 (counts only, no AUC)."""
    cfg = rp.load_config()
    entry = next(e for e in cfg["cohorts"] if e["accession"] == DATES_ACC)
    series = parse_geo.parse_series_matrix(rp.DATA / "series" / f"{DATES_ACC}_series_matrix.txt.gz", DATES_ACC)
    res = lab.assign_labels(series, cfg["label_rules"], entry.get("expected_counts"),
                            entry.get("exclude_samples", {}).get("ids"))
    dates = bd.parse_dates(series.characteristics, series.samples)
    all_per_date = pd.Series(list(dates.values())).value_counts().to_dict()
    return bd.date_counts(res.usable, dates), all_per_date, len(series.samples)


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

    g0 = gate0(X, y, cohort, batches)                                   # Gate 0a (the batch-split run's)
    print(f"Gate 0a max abs diff {g0['max_abs_diff']:.6f} {g0['diffs']}", flush=True)
    if not g0["ok"]:
        print("GATE 0a MISSED: nothing written")
        return 1
    rows, yv, bvec = g0["rows"], y.to_numpy()[g0["rows"]], g0["batch"]
    Xraw = X.to_numpy()[rows]
    Xws = rank_normalize_within_sample(X.iloc[rows]).to_numpy()
    lr = tf.make_model(LR, tf.BALANCED)
    gb = bs.gate_b(lr, Xraw, bvec)
    print(f"Gate B batch-AUC {gb['mean_auc']:.4f} (pooled {gb['pooled_auc']:.4f})", flush=True)

    on_split = pd.read_csv(rp.RESULTS / "batch_split.csv")
    on_split = on_split[(on_split.part == "arm_a") & (on_split.model == LR)]
    obs, g0c = {}, {}
    for tname, Xt in (("raw", Xraw), ("within_sample", Xws)):           # Gate 0c: the four batch-out AUCs
        arm = bs.run_arm_a(lr, Xt, yv, bvec)
        for d, r in arm.items():
            old = float(on_split[(on_split["transform"] == tname) & (on_split.direction == d)].auc.iloc[0])
            obs[(tname, d)] = r["auc"]
            g0c[f"{tname}/{d}"] = abs(r["auc"] - old)
    ok_c = max(g0c.values()) <= 1e-6
    print(f"Gate 0c max abs diff {max(g0c.values()):.2e} ok={ok_c}", flush=True)
    if not ok_c or not v.gate_b_passes(gb["mean_auc"]):
        print("GATE 0c / B MISSED: nothing written")
        return 1
    if args.gates_only:
        return 0

    n_pos, n_neg = int(yv[bvec == "1"].sum()), int((1 - yv[bvec == "1"]).sum())
    frames = []
    for tname, Xt in (("raw", Xraw), ("within_sample", Xws)):
        f = bd.run_null(lr, Xt, yv, bvec, n_pos, n_neg, v.PRIMARY_N)
        f.insert(0, "transform", tname)
        frames.append(f)
        print(f"  null {tname} done {time.time() - t0:.0f}s", flush=True)
    null = pd.concat(frames, ignore_index=True)
    if null.auc.isna().any():
        raise AssertionError("NaN AUC in the null: nothing written")
    old_ctrl = pd.read_csv(rp.RESULTS / "batch_control.csv")
    g0b = bd.gate_0b(null, old_ctrl)                                   # Gate 0b (card): seeds 0-19 reproduce
    print(f"Gate 0b {g0b}", flush=True)
    if not g0b["ok"]:
        print("GATE 0b MISSED: nothing written")
        return 1

    def draws(t, d):
        return null[(null["transform"] == t) & (null.direction == d)].sort_values("draw").auc.to_numpy()

    tails = {d: v.tail_counts(obs[("raw", d)], draws("raw", d)) for d in v.DIRECTIONS}
    prefix = {d: v.tail_counts(obs[("raw", d)], draws("raw", d)[:v.CARD_N]) for d in v.DIRECTIONS}
    shares = {d: dict(share_le=tails[d]["k_le"] / tails[d]["n"], share_lt=tails[d]["k_lt"] / tails[d]["n"],
                      wilson_le_95=v.wilson(tails[d]["k_le"], tails[d]["n"])) for d in v.DIRECTIONS}
    delta = {t: {d: float(draws(t, d).mean()) - obs[(t, d)] for d in v.DIRECTIONS} for t in ("raw", "within_sample")}
    verdicts = v.decide(gb["mean_auc"], tails, delta["raw"], delta["within_sample"])
    verdict_200 = v.r_batch_n(prefix)

    counts, all_per_date, n_total = arm_b_counts()
    gate_b_arm = v.power_gate(counts)
    dates_df = pd.DataFrame([dict(date=d, n_all=all_per_date.get(d, 0), n_bacterial=b, n_viral=vi,
                                  n_labelled=b + vi, fails_gate=(b < v.MIN_BACTERIAL or b + vi < v.MIN_ARRAYS))
                             for d, (b, vi) in counts.items()])
    arm_b = "run" if gate_b_arm["passes"] else "GSE40396 dates: underpowered, not run"
    if gate_b_arm["passes"]:
        raise SystemExit("power gate passed - not expected for GSE40396; Arm B was not planned to run, stop")

    meta = dict(date=date.today().isoformat(), plan_commit=PLAN_COMMIT, card="the deciding batch test", n_draws=v.PRIMARY_N,
                gate0a=dict(max_abs_diff=g0["max_abs_diff"], tolerance=GATE0_TOL),
                gate0b=g0b, gate0c=dict(max_abs_diff=max(g0c.values()), per_row=g0c),
                gate_b=dict(mean_auc=gb["mean_auc"], pooled_auc=gb["pooled_auc"], bar=v.GATE_B_MIN),
                batch_out_auc={f"{t}/{d}": a for (t, d), a in obs.items()},
                tails_raw=tails, shares_raw=shares, tails_raw_first200=prefix, verdict_first200_descriptive=verdict_200,
                null_summary={f"{t}/{d}": dict(mean=float(draws(t, d).mean()), sd=float(draws(t, d).std(ddof=1)),
                                               min=float(draws(t, d).min()), max=float(draws(t, d).max()))
                              for t in ("raw", "within_sample") for d in v.DIRECTIONS},
                delta_raw=delta["raw"], delta_within_sample=delta["within_sample"], verdicts=verdicts,
                arm_b=dict(status=arm_b, power_gate=gate_b_arm, counts={k: list(c) for k, c in counts.items()},
                           n_arrays_total=n_total, min_bacterial=v.MIN_BACTERIAL, min_arrays=v.MIN_ARRAYS),
                runtime_s=round(time.time() - t0, 1))
    null.to_csv(OUT_NULL, index=False)
    dates_df.to_csv(OUT_DATES, index=False)
    META.write_text(json.dumps(meta, indent=2, default=lambda o: o.item() if isinstance(o, np.generic) else float(o)),
                    encoding="utf-8")
    print(json.dumps({k: meta[k] for k in ("tails_raw", "shares_raw", "tails_raw_first200", "delta_raw",
                                           "verdicts", "verdict_first200_descriptive", "arm_b")}, indent=2, default=float))
    return 0


if __name__ == "__main__":
    sys.exit(main())
