"""Guards for the transfer matrix (LABBOOK PLAN 2026-11-05, the transfer-matrix run).

Each guard was checked to go red on a deliberately broken copy of the code it
protects (mutation list in the LABBOOK RESULT entry).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src import transfer as tf  # noqa: E402
from src import transfer_verdict as tv  # noqa: E402
from tests.test_gap_diagnostics import _synthetic  # noqa: E402

CFG = json.loads((ROOT / "config" / "cohorts.json").read_text(encoding="utf-8"))


def _data(shares=(0.3, 0.7, 0.5), n=60):
    X, y, c = _synthetic(list(shares), n=n, genes=30)
    return X.to_numpy(float), y.to_numpy().astype(int), c.to_numpy(), c


# (a) the test cohort never sits in training --------------------------------------------

def test_leaking_test_cohort_into_training_raises():
    Xv, yv, cv, _ = _data()
    tr, te = tf.pair_rows(cv, "C0", "C1")
    est = tf.make_model("logistic_l2", tf.BALANCED)
    tf.fit_auc(est, Xv, yv, cv, tr, te)                      # clean split works
    with pytest.raises(AssertionError, match="leakage"):
        tf.fit_auc(est, Xv, yv, cv, np.concatenate([tr, te[:3]]), te)


def test_pair_rows_are_disjoint_cohorts():
    _, _, cv, _ = _data()
    for a, b in tf.ordered_pairs(cv):
        tr, te = tf.pair_rows(cv, a, b)
        assert set(cv[tr]) == {a} and set(cv[te]) == {b}


# (b) the size cap is applied ------------------------------------------------------------

def test_size_cap_applied_in_every_row_and_stratified():
    Xv, yv, cv, _ = _data(n=80)
    size = tf.run_sizematched(Xv, yv, cv)
    assert (size.n_train == tf.SIZE_CAP).all()
    assert len(size) == 6 * tf.N_DRAWS
    rows = np.flatnonzero(cv == "C0")                        # share 0.3 -> 9 bacterial of 30
    d = tf.stratified_draw(yv, rows, 0)
    assert len(d) == 30 and yv[d].sum() == 9 and set(cv[d]) == {"C0"}
    assert not np.array_equal(d, tf.stratified_draw(yv, rows, 1))     # draws differ


# (c) platform / age tags ------------------------------------------------------------------

def test_real_config_tags_are_pinned():
    fam, age = tf.cohort_tags(CFG)
    assert fam == {"GSE63990": "affymetrix", "GSE6269": "affymetrix", "GSE42026": "illumina",
                   "GSE60244": "illumina", "GSE40396": "illumina"}
    assert age["GSE63990"] == "mixed" and not tf.age_testable(age)
    assert age["GSE42026"] == age["GSE40396"] == age["GSE6269"] == "children"
    assert age["GSE60244"] == "adults"
    pairs = tf.ordered_pairs(fam)
    assert len(pairs) == 20 and tf.cell_counts(pairs, fam) == {"same": 8, "cross": 12}


def test_wrong_tag_changes_the_cell_counts():
    fam, _ = tf.cohort_tags(CFG)
    pairs = tf.ordered_pairs(fam)
    bad = dict(fam, GSE63990="illumina")     # 4+1 split: 12 same pairs, not 8
    assert tf.cell_counts(pairs, bad) != tf.cell_counts(pairs, fam)


# (d) the class-weight switch is honoured ---------------------------------------------------

@pytest.mark.parametrize("model", ["logistic_l2", "lightgbm"])
def test_class_weight_switch_reaches_the_estimator(model):
    assert tf.class_weight_of(tf.make_model(model, None)) is None
    assert tf.class_weight_of(tf.make_model(model, tf.BALANCED)) == "balanced"


def test_class_weight_off_changes_the_auc_on_imbalanced_data():
    Xv, yv, cv, cohort = _data(shares=(0.1, 0.1, 0.1, 0.5), n=60)
    a = tf.run_loco(Xv, yv, cv, cohort, "logistic_l2", tf.BALANCED)
    b = tf.run_loco(Xv, yv, cv, cohort, "logistic_l2", None)
    assert (b.class_weight == "None").all() and (a.class_weight == "balanced").all()
    assert np.abs(a.auc.to_numpy() - b.auc.to_numpy()).max() > 1e-6


# (e) verdict boundaries (inclusive edges, and just beside them) -------------------------------

@pytest.mark.parametrize("main,ctrl,word", [
    (0.03, 0.02, "tracks"), (0.0299, 0.02, "undecided"), (0.03, 0.0199, "undecided"),
    (0.01, 0.0, "does not"), (0.0101, 0.0, "undecided"), (-0.05, -0.05, "does not"),
    (0.05, 0.0, "undecided")])
def test_tracks_boundaries(main, ctrl, word):
    assert tv.r_platform(main, ctrl) == word


def test_age_not_testable_and_both_and_scale():
    assert tv.r_age(0.5, 0.5, testable=False) == "not testable"
    assert tv.r_both("tracks", "tracks").startswith("confounded")
    assert tv.r_both("tracks", "does not") == "not fired"
    assert tv.r_scale("does not", 0.05, 0.0) == "not evaluated"
    assert tv.r_scale("tracks", 0.05, 0.025) == "mostly a scale offset"
    assert tv.r_scale("tracks", 0.05, 0.0251) == "not mostly a scale offset"


@pytest.mark.parametrize("deltas,word", [
    ({"a": 0.02, "b": 0.0}, "class weighting matters here"),
    ({"a": -0.02, "b": 0.0}, "class weighting matters here"),
    ({"a": 0.0099, "b": -0.0099}, "not where the drop lives"),
    ({"a": 0.01, "b": 0.0}, "undecided"), ({"a": 0.0199, "b": 0.0}, "undecided")])
def test_balance_boundaries(deltas, word):
    assert tv.r_balance(deltas) == word


def test_cell_mean_averages_draws_first_then_pairs():
    df = pd.DataFrame(dict(train=["A", "A", "A", "B"], test=["B", "B", "C", "A"],
                           auc=[0.6, 1.0, 0.5, 0.9]))
    tag = {"A": "x", "B": "x", "C": "y"}
    assert tf.cell_mean(df, tag, "same") == pytest.approx((0.8 + 0.9) / 2)
    assert tf.cell_mean(df, tag, "cross") == pytest.approx(0.5)
