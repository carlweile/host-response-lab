"""Guards for the deciding batch test (the deciding batch test). Each test builds the bug and asserts the guard fires.

Synthetic data only, plus the real count table of GSE40396 typed in as a fixture for the gate.
"""
import numpy as np
import pandas as pd
import pytest

from src import batch_deciding as bd
from src import batch_deciding_verdict as v
from src import batch_split as bs
from src import transfer as tf

SMALL = {("0", 1): 20, ("0", 0): 30, ("1", 1): 6, ("1", 0): 12}
GSE40396 = {"2010-08-03": (0, 2), "2011-07-22": (0, 8), "2012-04-24": (8, 12)}   # planning counts


def _cohort(seed: int = 0, genes: int = 20):
    rng = np.random.default_rng(seed)
    y, b = [], []
    for (batch, cls), n in SMALL.items():
        y += [cls] * n
        b += [batch] * n
    y, b = np.array(y), np.array(b)
    x = rng.normal(size=(len(y), genes))
    x[:, :5] += 1.0 * y[:, None]
    return x, y, b


def _both(k_le, k_lt, n=1000):
    return {d: dict(k_le=k_le, k_lt=k_lt, n=n) for d in v.DIRECTIONS}


# --- null: class-matched, seeds, prefix ------------------------------------------------

def test_null_draws_are_class_matched_every_time():
    x, y, b = _cohort()
    ctrl = bd.run_null(tf.make_model("logistic_l2", tf.BALANCED), x, y, b, 6, 12, n_draws=6)
    assert len(ctrl) == 12 and set(ctrl[ctrl.direction == "b0_to_b1"].n_test) == {18}


def test_null_not_class_matched_is_caught():
    _, y, b = _cohort()
    pos, neg = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    small = np.sort(np.concatenate([pos[:2], neg[:16]]))
    with pytest.raises(AssertionError, match="class-matched"):
        bs.check_partition(y, b, np.setdiff1d(np.arange(len(y)), small), small, 6, 12)


def test_longer_null_has_the_short_null_as_prefix_and_a_shifted_seed_fails_gate_0b():
    x, y, b = _cohort()
    est = tf.make_model("logistic_l2", tf.BALANCED)
    long = bd.run_null(est, x, y, b, 6, 12, n_draws=8)
    short = bs.run_control(est, x, y, b, 6, 12, n_draws=4)
    short.insert(0, "transform", "raw")
    long.insert(0, "transform", "raw")
    assert bd.gate_0b(long, short)["ok"]
    shifted = long.copy()
    shifted["draw"] = shifted["draw"] + 1                      # seed shifted by one
    assert not bd.gate_0b(shifted, short)["ok"]
    assert not bd.gate_0b(long.iloc[:6], short)["ok"]          # missing rows never pass


# --- tails and rule edges --------------------------------------------------------------

def test_tail_counts_tie_counts_as_at_or_below_but_not_below():
    t = v.tail_counts(0.9, [0.8, 0.9, 0.95])
    assert (t["k_le"], t["k_lt"], t["n"]) == (2, 1, 3)


def test_inverted_tail_lowest_auc_gets_smallest_k_le():
    draws = np.linspace(0.7, 0.99, 100)
    assert v.tail_counts(0.5, draws)["k_le"] < v.tail_counts(0.9, draws)["k_le"]
    assert v.tail_counts(0.5, draws)["k_le"] == 0 and v.tail_counts(1.0, draws)["k_le"] == 100


@pytest.mark.parametrize("n,k_ok,k_no", [(1000, 50, 51), (200, 10, 11)])
def test_costs_edge_is_inclusive_at_5_percent(n, k_ok, k_no):
    assert v.r_batch_n(_both(k_ok, k_ok, n)) == "costs AUC"
    assert v.r_batch_n(_both(k_no, k_no, n)) == "undecided"


@pytest.mark.parametrize("n,k_ok,k_no", [(1000, 250, 249), (200, 50, 49)])
def test_no_cost_edge_is_inclusive_at_25_percent(n, k_ok, k_no):
    assert v.r_batch_n(_both(k_ok + 5, k_ok, n)) == "no batch cost"
    assert v.r_batch_n(_both(k_no + 5, k_no, n)) == "undecided"


def test_both_directions_are_needed_and_the_outcomes_never_both_fire():
    t = {"b0_to_b1": dict(k_le=0, k_lt=0, n=1000), "b1_to_b0": dict(k_le=60, k_lt=59, n=1000)}
    assert v.r_batch_n(t) == "undecided"
    for le in range(0, 1001, 50):
        for lt in range(0, le + 1, 50):
            assert v.r_batch_n(_both(le, lt)) in ("costs AUC", "no batch cost", "undecided")
    with pytest.raises(ValueError):
        v.r_batch_n({"b0_to_b1": dict(k_le=0, k_lt=0, n=1000)})


def test_wilson_interval_brackets_the_share_and_is_wider_for_a_smaller_n():
    lo, hi = v.wilson(50, 1000)
    assert lo < 0.05 < hi
    lo2, hi2 = v.wilson(10, 200)
    assert hi2 - lo2 > hi - lo


def test_gate_b_first_missing_is_error_and_low_is_stop():
    with pytest.raises(ValueError):
        v.decide(None)
    assert v.decide(0.79)["r_batch"] == v.GATE_B_STOP
    with pytest.raises(ValueError):
        v.decide(0.95)
    out = v.decide(0.95, _both(0, 0), {d: 0.1 for d in v.DIRECTIONS},
                   {d: 0.02 for d in v.DIRECTIONS})
    assert (out["r_batch"], out["r_scale"]) == ("costs AUC", "mostly a scale offset")


# --- Arm B: dates, power gate -----------------------------------------------------------

def test_parse_dates_exactly_one_string():
    assert bd.parse_dates({"a": ['"hybridization date: 2012-04-24"', "x: y"]}, ["a"]) == {"a": "2012-04-24"}
    with pytest.raises(ValueError):
        bd.parse_dates({"a": ["x: y"]}, ["a"])
    with pytest.raises(ValueError):
        bd.parse_dates({"a": ["hybridization date: 1", "hybridization date: 2"]}, ["a"])


def test_date_counts_ignore_excluded_arrays():
    labels = {"a": "bacterial", "b": "viral", "c": "excluded"}
    dates = {"a": "d1", "b": "d1", "c": "d2"}
    assert bd.date_counts(labels, dates) == {"d1": (1, 1), "d2": (0, 0)}


def test_power_gate_real_gse40396_table_is_not_run():
    g = v.power_gate(GSE40396)
    assert not g["passes"]
    assert set(g["failing"]) == {"2010-08-03", "2011-07-22"}      # 2012-04-24 alone would pass (8 / 20)


def test_power_gate_edges():
    ok = {"a": (5, 5), "b": (5, 5)}
    assert v.power_gate(ok)["passes"]
    assert not v.power_gate({"a": (4, 20), "b": (5, 5)})["passes"]       # 4 bacterial
    assert not v.power_gate({"a": (5, 4), "b": (5, 5)})["passes"]        # 9 arrays
    assert not v.power_gate({"a": (50, 50)})["passes"]                   # one date is no split


class _Spy:
    calls = 0

    def fit(self, *a, **k):
        _Spy.calls += 1
        raise AssertionError("fit must not be reached")


def test_arm_b_never_fits_when_the_power_gate_fails():
    x, y, _ = _cohort()
    dv = np.array(["a"] * 40 + ["b"] * (len(y) - 40))
    _Spy.calls = 0
    with pytest.raises(RuntimeError, match="power gate"):
        bd.run_arm_b(_Spy(), x, y, dv, {"a": (4, 20), "b": (5, 5)}, n_draws=2)
    assert _Spy.calls == 0


def test_arm_b_date_leaking_into_training_is_caught():
    x, y, b = _cohort()
    dv = np.array(["a"] * 30 + ["b"] * 30 + ["c"] * (len(y) - 60))
    est = tf.make_model("logistic_l2", tf.BALANCED)
    te = np.flatnonzero(dv == "a")
    leaky = np.concatenate([np.flatnonzero(dv != "a"), te[:3]])
    with pytest.raises(AssertionError, match="leakage"):
        bs.fit_scores(est, x, y, dv, leaky, te)


def test_arm_b_runs_when_the_gate_passes_and_is_pure_dates():
    x, y, _ = _cohort()
    rng = np.random.default_rng(3)
    dv = np.array(rng.choice(["a", "b"], size=len(y)))
    counts = bd.date_counts({str(i): ("bacterial" if y[i] else "viral") for i in range(len(y))},
                            {str(i): dv[i] for i in range(len(y))})
    out = bd.run_arm_b(tf.make_model("logistic_l2", tf.BALANCED), x, y, dv, counts, n_draws=3)
    assert set(out["aucs"]) == {"a", "b"} and all(t["n"] == 3 for t in out["tails"].values())
    assert out["verdict"] in ("costs AUC", "no batch cost", "undecided")


def test_gate_0b_nan_row_missing_never_passes():
    old = pd.DataFrame(dict(transform=["raw"], direction=["b0_to_b1"], draw=[0], auc=[0.9]))
    new = pd.DataFrame(dict(transform=["raw"], direction=["b1_to_b0"], draw=[0], auc=[0.9]))
    assert not bd.gate_0b(new, old)["ok"]
