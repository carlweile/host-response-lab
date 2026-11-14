"""Guards for the label-shuffle null (run_label_null.py).

Each test builds the bug it guards against: a shuffle that does nothing, a
within-cohort shuffle that crosses cohorts, a Gate 0 that is bypassed, a null
that does not use the headline code path, and reading-rule edges.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import run_label_null as rl
from src import label_null as ln
from src import train


def _toy(n_per=30, genes=6, seed=0):
    rng = np.random.default_rng(seed)
    cohorts = np.repeat(["A", "B", "C"], n_per)
    y = np.tile([0, 1, 1], len(cohorts) // 3)
    X = rng.normal(size=(len(y), genes)) + y[:, None] * 0.8
    idx = [f"s{i}" for i in range(len(y))]
    return (pd.DataFrame(X, index=idx, columns=[f"g{j}" for j in range(genes)]),
            pd.Series(y, index=idx), pd.Series(cohorts, index=idx))


def test_noop_shuffle_trips_guard():
    _, y, cohort = _toy()
    with pytest.raises(AssertionError, match="changed no label"):
        ln.check_shuffle(y, y.copy(), cohort, "global")


def test_global_shuffle_changes_labels_keeps_totals():
    _, y, cohort = _toy()
    y2 = ln.shuffle_global(y, 3)
    assert ln.check_shuffle(y, y2, cohort, "global") > 0
    assert y2.sum() == y.sum()


def test_global_shuffle_posing_as_within_cohort_trips_guard():
    _, y, cohort = _toy()
    # cohort A has a different class mix than the pooled labels, so a global shuffle moves counts
    y = y.copy()
    y[cohort == "A"] = 0
    y[cohort == "B"] = 1
    y2 = ln.shuffle_global(y, 1)
    with pytest.raises(AssertionError, match="crossed cohorts"):
        ln.check_shuffle(y, y2, cohort, "within_cohort")


def test_real_within_cohort_shuffle_keeps_per_cohort_counts():
    _, y, cohort = _toy()
    y2 = ln.shuffle_within_cohort(y, cohort, 5)
    ln.check_shuffle(y, y2, cohort, "within_cohort")
    assert y.groupby(cohort).sum().equals(y2.groupby(cohort).sum())


def test_shuffle_is_seeded():
    _, y, cohort = _toy()
    assert ln.shuffle_global(y, 4).equals(ln.shuffle_global(y, 4))
    assert not ln.shuffle_global(y, 4).equals(ln.shuffle_global(y, 5))


def test_random_cv_auc_equals_headline_code_path():
    X, y, cohort = _toy()
    res = train.run_both_schemes(X, y, cohort, seed=0)
    ref = [r.auc for r in res if r.model == "logistic_l2" and r.scheme == "random_cv"]
    assert ln.random_cv_auc(X, y)["mean_auc"] == pytest.approx(np.mean(ref), abs=1e-12)


def _ref_csv(tmp_path, X, y, cohort, offset=0.0):
    rcv = ln.random_cv_auc(X, y)["mean_auc"]
    loco = ln.loco_auc(X, y, cohort)
    p = tmp_path / "summary.csv"
    pd.DataFrame([dict(model="logistic_l2", scheme="random_cv", mean_auc=rcv + offset),
                  dict(model="logistic_l2", scheme="leave_one_cohort_out", mean_auc=loco)]).to_csv(p, index=False)
    return p


def test_gate0_passes_and_writes_files(tmp_path):
    X, y, cohort = _toy()
    ref = _ref_csv(tmp_path, X, y, cohort)
    out_csv, out_json = tmp_path / "o.csv", tmp_path / "o.json"
    assert rl.run(X, y, cohort, out_csv, out_json, ref, seeds=[0, 1]) == 0
    assert len(pd.read_csv(out_csv)) == 4  # 2 nulls x 2 seeds
    assert out_json.exists()


def test_gate0_bypass_is_impossible(tmp_path):
    X, y, cohort = _toy()
    ref = _ref_csv(tmp_path, X, y, cohort, offset=1e-6)
    out_csv, out_json = tmp_path / "o.csv", tmp_path / "o.json"
    assert rl.run(X, y, cohort, out_csv, out_json, ref, seeds=[0]) == 1
    assert not out_csv.exists() and not out_json.exists()


@pytest.mark.parametrize("mean,single,word", [
    (0.45, 0.50, "CLEAN"), (0.55, 0.50, "CLEAN"), (0.5501, 0.50, "LEAK"),
    (0.50, 0.65, "CLEAN"), (0.50, 0.6501, "LEAK"), (0.4499, 0.50, "ODD-LOW"),
])
def test_verdict_edges(mean, single, word):
    assert ln.verdict(mean, single) == word


@pytest.mark.parametrize("mean,word", [
    (0.5999, "SMALL"), (0.60, "MEDIUM"), (0.6999, "MEDIUM"), (0.70, "LARGE"),
])
def test_shortcut_edges(mean, word):
    assert ln.shortcut_word(mean) == word


def test_leaky_pipeline_would_be_caught():
    """If labels leaked into features, the global shuffle would still score high: verdict LEAK."""
    X, y, cohort = _toy()
    aucs = []
    for s in range(5):
        y2 = ln.shuffle_global(y, s)
        X_leak = X.copy()
        X_leak["g0"] = y2.to_numpy() * 3.0  # the bug: a feature built from the shuffled label
        aucs.append(ln.random_cv_auc(X_leak, y2)["mean_auc"])
    assert ln.verdict(float(np.mean(aucs)), max(aucs)) == "LEAK"
