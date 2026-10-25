"""Guards for the seesaw test (LABBOOK PLAN 9daaf08, tests a-g).

Each guard was checked to go red on a deliberately broken copy of the code it
protects (see the LABBOOK RESULT entry for the list of mutations).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import seesaw as ss  # noqa: E402
from src import seesaw_verdict as sv  # noqa: E402
from src.gap_diagnostics import pair_auc  # noqa: E402
from tests.test_gap_diagnostics import _synthetic  # noqa: E402


# (a) k folds of size about n/k, and C_k follows the formula -------------------------

@pytest.mark.parametrize("k", [2, 10])
def test_runner_uses_k_folds_of_size_n_over_k(k):
    X, y, c = _synthetic([0.3, 0.7], n=40, genes=20)
    folds, _ = ss.run_arm(X, y, c, "W", ss.logistic(ss.c_for_k(k)), k, 100, {})
    assert len(folds) == k
    assert all(abs(f["n_test"] - len(y) / k) <= 1 for f in folds)


def test_c_k_formula_and_the_model_carries_it():
    assert ss.c_for_k(2) == pytest.approx(0.08)
    assert ss.c_for_k(5) == pytest.approx(0.05)
    assert ss.c_for_k(10) == pytest.approx(0.04 * 10 / 9)
    assert ss.logistic(ss.c_for_k(2)).named_steps["clf"].C == pytest.approx(0.08)


# (b) cross pairs are counted only inside one fold ------------------------------------

def test_pairs_never_cross_a_fold():
    X, y, c = _synthetic([0.3, 0.7], n=40, genes=20)
    folds, _ = ss.run_arm(X, y, c, "W", ss.logistic(0.05), 5, 100, {})
    for f in folds:
        n_pos = f["n_pos_test"]
        assert f["n_pairs_within"] + f["n_pairs_cross"] == n_pos * (f["n_test"] - n_pos)
    pooled_pairs = int(y.sum()) * int((1 - y).sum())
    assert sum(f["n_pairs_within"] + f["n_pairs_cross"] for f in folds) < pooled_pairs


def test_pair_weighting_over_folds():
    rows = pd.DataFrame([dict(set="s", part="p1", model="m", arm="W", k=2, C=0.08, seed=1,
                              auc_within=1.0, auc_cross=a, n_pairs_within=10,
                              n_pairs_cross=n, algebra=0.0, pin=0.0, n_iter=1)
                         for a, n in ((1.0, 10), (0.0, 30))])
    out = sv.pair_weighted(rows).iloc[0]
    assert out.cross == pytest.approx(0.25)
    assert out.D == pytest.approx(0.75)


# (c) center_whole zeroes cohort means; *_train twins read no test row -----------------

def test_center_whole_gives_zero_cohort_means():
    X, y, c = _synthetic([0.3, 0.7], n=30, genes=10)
    M = ss.whole_matrix("center_whole", X + 5.0 * (c == "C1").to_numpy()[:, None], c)
    means = pd.DataFrame(M).groupby(c.to_numpy()).mean()
    assert np.abs(means.to_numpy()).max() < 1e-12


@pytest.mark.parametrize("arm", ["center_train", "scale_train", "zscore_train"])
def test_train_twins_ignore_test_rows(arm):
    X, y, c = _synthetic([0.3, 0.7], n=30, genes=10)
    tr, te = np.arange(0, 60, 2), np.arange(1, 60, 2)
    A_tr, A_te = ss.train_matrices(arm, X, c, tr, te)
    X2 = X.copy()
    X2.iloc[te[0]] += 100.0
    B_tr, B_te = ss.train_matrices(arm, X2, c, tr, te)
    np.testing.assert_array_equal(A_tr, B_tr)
    np.testing.assert_array_equal(A_te[1:], B_te[1:])


# (d) scale_whole: cohort SD 1, cohort means left free ----------------------------------

def test_scale_whole_fixes_spread_not_mean():
    X, y, c = _synthetic([0.3, 0.7], n=30, genes=10)
    X = X + 3.0 * (c == "C1").to_numpy()[:, None] + 1.0
    M = pd.DataFrame(ss.whole_matrix("scale_whole", X, c))
    sds = M.groupby(c.to_numpy()).std(ddof=0).to_numpy()
    np.testing.assert_allclose(sds, 1.0, atol=1e-12)
    means = M.groupby(c.to_numpy()).mean()
    assert (means.loc["C0"] - means.loc["C1"]).abs().min() > 0.5


# (e) stumps are gene-additive: equal cohort totals on whole-cohort ranks ---------------

def _cohort_totals(model, X, y, c):
    W = ss.whole_matrix("W", X, c)
    est = model.fit(W, y.to_numpy())
    s = est.predict(W, raw_score=True)
    return pd.Series(s).groupby(c.to_numpy()).mean().to_numpy()


def test_stump_totals_match_and_deep_trees_do_not():
    X, y, c = _synthetic([0.2, 0.5, 0.8], n=60, genes=20)
    stumps = _cohort_totals(ss.lgbm(stumps=True), X, y, c)
    assert np.ptp(stumps) < 1e-9
    deep = _cohort_totals(ss.lgbm(), X, y, c)
    assert np.ptp(deep) > 1e-6


# (f) a constant on one cohort moves o_g by that constant, not the same-cohort AUC ------

def test_offset_moves_with_a_cohort_shift_and_within_auc_does_not():
    rng = np.random.default_rng(3)
    s = rng.normal(size=40)
    y = np.tile([0, 1], 20)
    c = np.repeat(["A", "B"], 20)
    shifted = s + 2.5 * (c == "A")
    assert ss.class_offset(shifted, y, c, "A") - ss.class_offset(s, y, c, "A") == \
        pytest.approx(2.5 * (1 - 0.5))   # the pooled class mean moves by half the shift
    gap = ss.class_offset(shifted, y, c, "A") - ss.class_offset(shifted, y, c, "B")
    assert gap - (ss.class_offset(s, y, c, "A") - ss.class_offset(s, y, c, "B")) ==         pytest.approx(2.5)               # the A-minus-B offset moves by the full shift
    assert pair_auc(y, shifted, c).within == pair_auc(y, s, c).within


# (g) the verdict rules on a hand-made table ---------------------------------------------

BASE = dict(X2=0.02, X5=0.05, X10=0.06, S2=0.1, S5=0.2, S10=0.3, Yc=0.04, Ys=0.005,
            X_logistic=0.05, X_lightgbm=0.01, X_lgbm_stumps=0.04,
            algebra_lightgbm=1.0, algebra_lgbm_stumps=0.1, algebra_logistic=0.1,
            pin_lightgbm=1.0, pin_lgbm_stumps=0.8, max_n_iter=100.0,
            gate_center_ok=True, gate_scale_ok=True)


def test_hand_table_passes_everything():
    assert sv.judge(BASE) == {"P1": "PASS", "P2": "PASS", "P3": "PASS"}
    assert sv.overall({"all": sv.judge(BASE), "no40396": sv.judge(BASE)}) == "SUPPORTED"


@pytest.mark.parametrize("change, part, word", [
    (dict(X5=0.0), "P1", "KILL"),
    (dict(X2=0.045), "P1", "KILL"),           # X2 >= 0.8 X5
    (dict(X2=0.03), "P1", "PARTIAL"),         # between 0.5 and 0.8 of X5
    (dict(X10=0.015), "P1", "KILL"),          # X10 < X2
    (dict(S10=0.15), "P1", "PARTIAL"),        # slope order broken
    (dict(Yc=0.0), "P2", "KILL"),
    (dict(Yc=0.01), "P2", "KILL"),            # Yc < 0.25 X5
    (dict(Ys=0.02), "P2", "KILL"),            # Ys >= 0.5 Yc
    (dict(Ys=0.015), "P2", "PARTIAL"),
    (dict(gate_center_ok=False), "P2", "VOID"),
    (dict(X_lgbm_stumps=0.01), "P3", "KILL"),  # stumps <= default
    (dict(X_lgbm_stumps=0.015), "P3", "PARTIAL"),
    (dict(pin_lgbm_stumps=0.3), "P3", "NOT INTERPRETABLE"),
    (dict(algebra_lgbm_stumps=0.9), "P3", "STOP"),
    (dict(max_n_iter=3000.0), "P1", "VOID"),
])
def test_each_rule_fires(change, part, word):
    v = sv.judge(dict(BASE, **change))
    assert v[part].startswith(word)
    assert sv.overall({"all": sv.judge(BASE), "no40396": v}).startswith("NOT SUPPORTED")
