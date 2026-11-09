"""Guards for the batch-split run (the batch-split run). Each test builds the bug and asserts the guard fires.

Synthetic data only: a small two-batch cohort, no files from data/ or results/.
"""
import numpy as np
import pytest

from src import batch_split as bs
from src import batch_split_verdict as bv
from src import transfer as tf

SMALL = {("0", 1): 20, ("0", 0): 30, ("1", 1): 6, ("1", 0): 12}


def _cohort(batch_shift: float = 0.0, seed: int = 0, genes: int = 30):
    rng = np.random.default_rng(seed)
    y, b = [], []
    for (batch, cls), n in SMALL.items():
        y += [cls] * n
        b += [batch] * n
    y, b = np.array(y), np.array(b)
    x = rng.normal(size=(len(y), genes))
    x[:, :5] += 1.0 * y[:, None]
    x[:, 5:10] += batch_shift * (b == "1")[:, None]
    return x, y, b


# --- batch parser + count gate ---------------------------------------------------------

def test_parse_batches_reads_one_string_per_array():
    chars = {"a": ['"batch: 0"', "age: 3"], "b": ["batch: 1"]}
    assert bs.parse_batches(chars, ["a", "b"]) == {"a": "0", "b": "1"}


def test_parse_batches_rejects_zero_or_two_strings():
    with pytest.raises(ValueError):
        bs.parse_batches({"a": ["age: 3"]}, ["a"])
    with pytest.raises(ValueError):
        bs.parse_batches({"a": ["batch: 0", "batch: 1"]}, ["a"])


def test_count_gate_rejects_a_wrong_table():
    _, y, b = _cohort()
    bs.assert_expected_counts(y, b, SMALL)
    y2 = y.copy()
    y2[0] = 1 - y2[0]
    with pytest.raises(AssertionError):
        bs.assert_expected_counts(y2, b, SMALL)


# --- guard (i): batch must not leak into training -------------------------------------

def test_batch_out_fit_raises_when_a_batch_sits_on_both_sides():
    x, y, b = _cohort()
    est = tf.make_model("logistic_l2", tf.BALANCED)
    b1 = np.flatnonzero(b == "1")
    leaky_train = np.concatenate([np.flatnonzero(b == "0"), b1[:3]])   # 3 batch-1 arrays in training
    with pytest.raises(AssertionError, match="leakage"):
        bs.fit_scores(est, x, y, b, leaky_train, b1)


def test_arm_a_sides_are_pure_batches_and_run():
    x, y, b = _cohort()
    arm = bs.run_arm_a(tf.make_model("logistic_l2", tf.BALANCED), x, y, b)
    assert set(arm) == set(bv.DIRECTIONS)
    assert set(b[arm["b0_to_b1"]["tr"]]) == {"0"} and set(b[arm["b0_to_b1"]["te"]]) == {"1"}
    assert set(b[arm["b1_to_b0"]["tr"]]) == {"1"} and set(b[arm["b1_to_b0"]["te"]]) == {"0"}


# --- guard (ii): the control must be class-matched ------------------------------------

def test_partition_matches_the_held_out_batch_counts():
    _, y, b = _cohort()
    for d in range(5):
        big, small = bs.partition(y, d, 6, 12)
        bs.check_partition(y, b, big, small, 6, 12)
        assert (len(small), int(y[small].sum())) == (18, 6)
        assert (len(big), int(y[big].sum())) == (len(y) - 18, int(y.sum()) - 6)


def test_control_not_class_matched_fails():
    _, y, b = _cohort()
    pos, neg = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    small = np.sort(np.concatenate([pos[:3], neg[:15]]))          # 18 arrays but 3/15, not 6/12
    big = np.setdiff1d(np.arange(len(y)), small)
    with pytest.raises(AssertionError, match="class-matched"):
        bs.check_partition(y, b, big, small, 6, 12)


def test_control_overlap_or_unmixed_side_fails():
    _, y, b = _cohort()
    big, small = bs.partition(y, 0, 6, 12)
    with pytest.raises(AssertionError):
        bs.check_partition(y, b, np.concatenate([big, small[:1]]), small, 6, 12)
    only_b0 = np.flatnonzero((b == "0") & (y == 1))[:6]
    only_b0 = np.sort(np.concatenate([only_b0, np.flatnonzero((b == "0") & (y == 0))[:12]]))
    with pytest.raises(AssertionError, match="mix"):
        bs.check_partition(y, b, np.setdiff1d(np.arange(len(y)), only_b0), only_b0, 6, 12)


def test_run_control_sizes_follow_the_direction():
    x, y, b = _cohort()
    ctrl = bs.run_control(tf.make_model("logistic_l2", tf.BALANCED), x, y, b, 6, 12, n_draws=4)
    assert len(ctrl) == 8
    d1, d2 = ctrl[ctrl.direction == "b0_to_b1"], ctrl[ctrl.direction == "b1_to_b0"]
    assert set(d1.n_train) == {len(y) - 18} and set(d1.n_test) == {18}
    assert set(d2.n_train) == {18} and set(d2.n_test) == {len(y) - 18}


# --- guard (iii): Gate B ---------------------------------------------------------------

def test_gate_b_sees_a_planted_batch_shift_and_not_a_clean_cohort():
    est = tf.make_model("logistic_l2", tf.BALANCED)
    x, y, b = _cohort(batch_shift=3.0)
    assert bs.gate_b(est, x, b)["mean_auc"] >= bv.GATE_B_MIN
    x0, _, b0 = _cohort(batch_shift=0.0)
    assert bs.gate_b(est, x0, b0)["mean_auc"] < bv.GATE_B_MIN


def test_no_verdict_when_gate_b_is_missing_or_below_the_bar():
    for missing in (None, float("nan")):
        with pytest.raises(ValueError):
            bv.decide(missing)
    stop = bv.decide(0.79)
    assert stop["r_batch"] == bv.GATE_B_STOP and stop["r_scale"] == "not evaluated"
    with pytest.raises(ValueError):           # passed, but no arm results to judge
        bv.decide(0.95)


def test_gate_b_edge_is_inclusive():
    assert bv.gate_b_passes(0.80) and not bv.gate_b_passes(0.7999)


# --- guard (iv): rule edges ------------------------------------------------------------

def _both(v):
    return {d: v for d in bv.DIRECTIONS}


def test_costs_auc_edge_is_inclusive_and_needs_both_directions_and_below_all():
    assert bv.r_batch(_both(0.03), _both(True)) == "costs AUC"
    assert bv.r_batch(_both(0.0299), _both(True)) == "undecided"
    assert bv.r_batch({"b0_to_b1": 0.2, "b1_to_b0": 0.0299}, _both(True)) == "undecided"
    assert bv.r_batch(_both(0.2), {"b0_to_b1": True, "b1_to_b0": False}) == "undecided"


def test_delta_from_floats_hits_the_edge():
    assert bv.r_batch(_both(bv.delta([0.90] * 20, 0.87)), _both(True)) == "costs AUC"


def test_no_batch_cost_edge_and_negative_delta():
    assert bv.r_batch(_both(0.01), _both(False)) == "no batch cost"
    assert bv.r_batch(_both(-0.05), _both(False)) == "no batch cost"
    assert bv.r_batch(_both(0.0101), _both(False)) == "undecided"
    assert bv.r_batch({"b0_to_b1": 0.0, "b1_to_b0": 0.02}, _both(False)) == "undecided"


def test_tie_with_the_control_minimum_is_not_below_all():
    draws = [0.9, 0.95, 0.97]
    assert not bv.below_all(0.9, draws) and bv.below_all(0.8999, draws)


def test_r_scale_only_when_r_batch_fired_and_edge_is_half():
    raw = _both(0.10)
    assert bv.r_scale("undecided", raw, _both(0.0)) == "not evaluated"
    assert bv.r_scale("no batch cost", raw, _both(0.0)) == "not evaluated"
    assert bv.r_scale("costs AUC", raw, _both(0.05)) == "mostly a scale offset"
    assert bv.r_scale("costs AUC", raw, _both(0.0501)) == "not mostly a scale offset"
    assert bv.r_scale("costs AUC", raw, {"b0_to_b1": 0.0, "b1_to_b0": 0.06}) == "not mostly a scale offset"


def test_decide_reads_r_scale_only_after_a_fired_r_batch():
    out = bv.decide(0.95, _both(0.01), _both(False), _both(0.0))
    assert (out["r_batch"], out["r_scale"]) == ("no batch cost", "not evaluated")
    out = bv.decide(0.95, _both(0.10), _both(True), _both(0.02))
    assert (out["r_batch"], out["r_scale"]) == ("costs AUC", "mostly a scale offset")


# --- bootstrap -------------------------------------------------------------------------

def test_bootstrap_interval_brackets_the_point_estimate_and_counts_skips():
    rng = np.random.default_rng(1)
    y = np.array([0] * 30 + [1] * 10)
    s = rng.normal(size=40) + y
    ci = bs.bootstrap_ci(y, s, n=300)
    assert ci["lo"] <= bs.auc_of(y, s) <= ci["hi"] and ci["skipped"] == 0
    y_rare = np.array([0] * 39 + [1])
    assert bs.bootstrap_ci(y_rare, rng.normal(size=40), n=300)["skipped"] > 0
