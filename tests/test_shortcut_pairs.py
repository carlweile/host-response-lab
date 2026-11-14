"""Guards for the within-cohort-pair AUC check (run_shortcut_pairs.py).

Each test builds the bug it guards against: within and cross pairs mixed, a plain fold
mean instead of pair weights, verdict edges off by one comparison, an instrument window
that lets a bad control through, a Gate 0 that can be skipped, a planted offset that
leaks into within-cohort pairs.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from sklearn.linear_model import LogisticRegression

import run_shortcut_pairs as rs
from src import shortcut_extras as ex
from src import shortcut_pairs as sp
from src.gap_diagnostics import PairAUC


def _toy(n_per=40, genes=5, seed=0):
    rng = np.random.default_rng(seed)
    cohort = np.repeat(["A", "B", "C"], n_per)
    share = {"A": 0.25, "B": 0.5, "C": 0.75}
    y = np.concatenate([(rng.random(n_per) < share[c]).astype(int) for c in "ABC"])
    X = rng.normal(size=(len(y), genes)) + y[:, None] * 0.7
    idx = [f"s{i}" for i in range(len(y))]
    return (pd.DataFrame(X, index=idx, columns=[f"g{j}" for j in range(genes)]),
            pd.Series(y, index=idx), pd.Series(cohort, index=idx))


def _folds(y, cohort):
    X = _toy()[0].to_numpy(dtype=float)
    return sp.fold_scores(X, y, LogisticRegression(max_iter=500), 0)


def test_constant_added_to_one_cohort_leaves_within_unchanged():
    _, y, cohort = _toy()
    folds = _folds(y, cohort)
    base = sp.seed_row(y, cohort, folds, "logit")
    off = np.where(cohort.to_numpy() == "A", 25.0, 0.0)
    moved = sp.seed_row(y, cohort, folds, "logit", off)
    assert moved["within"] == pytest.approx(base["within"], abs=1e-12)
    assert abs(moved["pooled"] - base["pooled"]) > 0.05     # pooled DOES see the offset


def test_aggregate_is_pair_weighted_not_a_plain_fold_mean():
    pairs = [PairAUC(0.9, 0.9, 0.9, 10, 10), PairAUC(0.5, 0.5, 0.5, 1000, 1000)]
    a = sp.aggregate(pairs)
    assert a["within"] == pytest.approx((0.9 * 10 + 0.5 * 1000) / 1010)
    assert a["within"] != pytest.approx(0.7)


def test_pooled_is_pair_weighted_mix_of_within_and_cross():
    _, y, cohort = _toy()
    r = sp.seed_row(y, cohort, _folds(y, cohort))
    mix = (r["within"] * r["n_within"] + r["cross"] * r["n_cross"]) / (r["n_within"] + r["n_cross"])
    assert r["pooled"] == pytest.approx(mix, abs=1e-12)


@pytest.mark.parametrize("d,word", [(0.01, "NOT-USED"), (0.0100001, "UNCLEAR"),
                                    (0.0199999, "UNCLEAR"), (0.02, "USED"), (-0.05, "NOT-USED")])
def test_verdict_edges(d, word):
    assert sp.verdict_word(d) == word


@pytest.mark.parametrize("m,ok", [(0.45, True), (0.55, True), (0.5, True),
                                  (0.4499, False), (0.5501, False)])
def test_instrument_window_edges(m, ok):
    assert sp.instrument_ok(m) is ok


def _ref_from(folds, y, cohort, model="m"):
    from src.gap_diagnostics import pair_auc
    yv, cv = y.to_numpy(), cohort.to_numpy()
    rows = []
    for f in folds:
        pa = pair_auc(yv[f.te], f.p, cv[f.te])
        rows.append(dict(transform="raw", model=model, scheme="random_cv", fold=f.name,
                         auc_pooled=pa.pooled, auc_within=pa.within, auc_cross=pa.cross,
                         n_pairs_within=pa.n_within, n_pairs_cross=pa.n_cross))
    return pd.DataFrame(rows)


def test_gate0_passes_on_identical_and_fails_on_nudge():
    X, y, cohort = _toy()
    folds = sp.fold_scores(X.to_numpy(dtype=float), y, LogisticRegression(max_iter=500), 0)
    ref = _ref_from(folds, y, cohort)
    assert sp.gate0_rows({"m": folds}, y, cohort, ref)["ok"]
    ref["auc_within"] += 1e-6
    assert not sp.gate0_rows({"m": folds}, y, cohort, ref)["ok"]


def test_run_cannot_skip_gate0_and_writes_nothing(tmp_path):
    X, y, cohort = _toy()
    models = {"logistic_l2": LogisticRegression(max_iter=500)}
    Xv = X.to_numpy(dtype=float)
    ref = _ref_from(sp.fold_scores(Xv, y, models["logistic_l2"], 0), y, cohort, "logistic_l2")
    ref["auc_pooled"] += 1e-6
    ref_path = tmp_path / "ref.csv"
    ref.to_csv(ref_path, index=False)
    outs = dict(out_csv=tmp_path / "a.csv", out_planted=tmp_path / "b.csv", out_json=tmp_path / "c.json")
    code = rs.run(X, y, cohort, models=models, ref_folds=ref_path, ref_null=None, seeds=[0, 1],
                  boot_b=5, **outs)
    assert code == 1
    assert not any(p.exists() for p in outs.values())


def test_run_end_to_end_on_toy_writes_all_files(tmp_path):
    X, y, cohort = _toy()
    models = {"logistic_l2": LogisticRegression(max_iter=500)}
    ref = _ref_from(sp.fold_scores(X.to_numpy(dtype=float), y, models["logistic_l2"], 0), y,
                    cohort, "logistic_l2")
    ref_path = tmp_path / "ref.csv"
    ref.to_csv(ref_path, index=False)
    outs = dict(out_csv=tmp_path / "a.csv", out_planted=tmp_path / "b.csv", out_json=tmp_path / "c.json")
    # toy data carries real signal, so the shuffled-label instrument may land outside 0.45-0.55;
    # either branch must be internally consistent
    code = rs.run(X, y, cohort, models=models, ref_folds=ref_path, ref_null=None, seeds=[0, 1],
                  boot_b=20, **outs)
    assert code == 0 and outs["out_json"].exists()
    import json
    out = json.loads(outs["out_json"].read_text())
    if out["verdict"] == "NO VERDICT (instrument)":
        assert not outs["out_csv"].exists()
    else:
        assert len(pd.read_csv(outs["out_csv"])) == 2 and outs["out_planted"].exists()


def test_planted_offset_is_constant_inside_a_cohort_and_moves_pooled_only():
    _, y, cohort = _toy()
    folds = _folds(y, cohort)
    off = ex.planted_offset(y, cohort, folds, 1.0)
    for c in "ABC":
        assert np.ptp(off[cohort.to_numpy() == c]) == 0.0
    df = ex.planted_rows(y, cohort, {0: folds})
    s = ex.planted_summary(df)
    assert s["within_shift"] < 1e-9
    assert s["sensitivity_move_at_k"] > 0.0 and s["passed"] == (s["sensitivity_move_at_k"] >= 0.02)


def test_oracle_offset_never_below_no_offset_and_leaves_within_alone():
    _, y, cohort = _toy()
    folds = _folds(y, cohort)
    base = sp.seed_row(y, cohort, folds, "logit")
    o = ex.oracle_offset_D(y, cohort, folds)
    assert o["pooled"] >= base["pooled"] and o["within"] == pytest.approx(base["within"], abs=1e-12)


def test_oracle_finds_a_planted_shortcut_when_scores_are_bad():
    # scores carry NO biology; class mix differs by hospital: the best offset must find it
    _, y, cohort = _toy()
    rng = np.random.default_rng(1)
    X = rng.normal(size=(len(y), 3))
    folds = sp.fold_scores(X, y, LogisticRegression(max_iter=500), 0)
    base = sp.seed_row(y, cohort, folds, "logit")
    assert ex.oracle_offset_D(y, cohort, folds)["D"] > base["D"] + 0.02


def test_bootstrap_interval_contains_point_estimate_and_is_ordered():
    _, y, cohort = _toy()
    folds = {0: _folds(y, cohort)}
    b = ex.bootstrap_mean_d(y, cohort, folds, B=200, seed=1)
    point = sp.seed_row(y, cohort, folds[0])["D"]
    assert b["lo"] < b["hi"] and b["lo"] - 0.05 <= point <= b["hi"] + 0.05


def test_bootstrap_with_unit_counts_reproduces_the_point_estimate():
    _, y, cohort = _toy()
    folds = _folds(y, cohort)
    yv, cv = y.to_numpy(), cohort.to_numpy()
    pairs = []
    for f in folds:
        pos, neg = yv[f.te] == 1, yv[f.te] == 0
        a, b = f.p[pos][:, None], f.p[neg][None, :]
        same = cv[f.te][pos][:, None] == cv[f.te][neg][None, :]
        pairs.append(ex._weighted_pair((a > b) + 0.5 * (a == b), same,
                                       np.ones(pos.sum()), np.ones(neg.sum())))
    assert sp.aggregate(pairs)["D"] == pytest.approx(sp.seed_row(y, cohort, folds)["D"], abs=1e-12)
