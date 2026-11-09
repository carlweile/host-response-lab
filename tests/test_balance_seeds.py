"""Guards for Arm B with N LightGBM seeds (LABBOOK PLAN 2026-11-05, the 10-seed Arm B run)."""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src import transfer as tf  # noqa: E402
from src import transfer_verdict as tv  # noqa: E402
from tests.test_transfer import _data  # noqa: E402

RES = ROOT / "results"


def test_seed_count_flows_through_to_the_rows():
    Xv, yv, cv, c = _data()
    bal = tf.run_balance(Xv, yv, cv, c, range(4))
    g = bal[bal.model == "lightgbm"]
    assert sorted(g.model_seed.unique()) == [0, 1, 2, 3]
    assert len(g) == 4 * 2 * 3                       # seeds x weighting x held-out cohorts
    assert set(bal[bal.model == "logistic_l2"].model_seed) == {0}


def test_summary_is_cohort_mean_of_seed_means_and_spread_matches():
    rows = []
    for s, none_auc in enumerate((0.9, 0.7)):        # seed 0: +0.2, seed 1: 0.0 vs balanced 0.7
        for cw, a in (("balanced", 0.7), ("None", none_auc)):
            for ho in ("A", "B"):
                rows.append(dict(model="lightgbm", class_weight=cw, model_seed=s, held_out=ho, auc=a))
    out = tf.balance_summary(pd.DataFrame(rows))
    assert out["lightgbm"] == pytest.approx(0.1)
    assert out["lightgbm_per_seed"] == {0: pytest.approx(0.2), 1: pytest.approx(0.0)}
    sp = out["lightgbm_seed_spread"]
    assert (sp["min"], sp["max"]) == (pytest.approx(0.0), pytest.approx(0.2))
    assert sp["sd"] == pytest.approx(0.2 / 2 ** 0.5)


@pytest.mark.parametrize("lgbm,word", [(0.0195, "undecided"), (0.02, "class weighting matters here"),
                                       (-0.03, "class weighting matters here"), (0.009, "not where the drop lives"),
                                       (0.010, "undecided")])
def test_verdict_words_on_all_three_bands(lgbm, word):
    assert tv.r_balance({"logistic_l2": 0.0005, "lightgbm": lgbm}) == word


def test_seeds_0_2_of_the_10_seed_run_reproduce_the_committed_deltas():
    new_f = RES / "balance_control_10seeds.csv"
    if not new_f.exists():
        pytest.skip("10-seed run not written yet")
    old = pd.read_csv(RES / "balance_control.csv", keep_default_na=False)   # "None" is a string
    new = pd.read_csv(new_f, keep_default_na=False)
    assert set(old.class_weight) == {"balanced", "None"}
    for m in ("logistic_l2", "lightgbm"):
        a = tf.balance_summary(old)[m]
        b = tf.balance_summary(new[(new.model != "lightgbm") | (new.model_seed < 3)])[m]
        assert a == pytest.approx(b, abs=1e-6)
    assert sorted(new[new.model == "lightgbm"].model_seed.unique()) == list(range(10))
