"""Guards for the seed sweep (run_seed_sweep.py).

Each test builds the bug it guards against: a seed that is silently dropped, a
LOCO side that moves with the fold seed, and a kill rule that fires too easily
or not at all.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import run_seed_sweep as rs
from src import splits, train


def _toy(n_per=24, genes=8, seed=0):
    rng = np.random.default_rng(seed)
    cohorts = np.repeat(["A", "B", "C"], n_per)
    y = np.tile([0, 1], len(cohorts) // 2)
    X = rng.normal(size=(len(y), genes)) + y[:, None] * 0.8
    idx = [f"s{i}" for i in range(len(y))]
    return (pd.DataFrame(X, index=idx, columns=[f"g{j}" for j in range(genes)]),
            pd.Series(y, index=idx), pd.Series(cohorts, index=idx))


def test_fold_seed_reaches_random_splits(monkeypatch):
    X, y, cohort = _toy()
    seen = []
    real = splits.random_splits

    def spy(y_, n_splits=5, seed=0):
        seen.append(seed)
        return real(y_, n_splits=n_splits, seed=seed)

    monkeypatch.setattr(train, "random_splits", spy)
    train.run_both_schemes(X, y, cohort, seed=207)
    assert seen and set(seen) == {207}


def test_different_fold_seeds_give_different_folds():
    _, y, _ = _toy()
    a = [te.tolist() for _, (_, te) in splits.random_splits(y, seed=200)]
    b = [te.tolist() for _, (_, te) in splits.random_splits(y, seed=201)]
    assert a != b


def test_model_seed_reaches_lightgbm():
    pytest.importorskip("lightgbm")
    m = train.make_models(model_seed=203)
    assert m["lightgbm"].get_params()["random_state"] == 203


def test_loco_ignores_fold_seed():
    X, y, cohort = _toy()
    a = pd.DataFrame(rs.sweep_rows(X, y, cohort, seed=200)).set_index("model")
    b = pd.DataFrame(rs.sweep_rows(X, y, cohort, seed=209)).set_index("model")
    assert np.allclose(a["loco_auc"], b["loco_auc"])


def _sweep(gaps_logistic, gaps_lgbm):
    rows = [dict(model="logistic_l2", gap=g) for g in gaps_logistic]
    rows += [dict(model="lightgbm", gap=g) for g in gaps_lgbm]
    return pd.DataFrame(rows)


def test_kill_fires_only_when_both_models_mostly_small():
    small, big = [0.005] * 6 + [0.05] * 4, [0.07] * 10
    assert rs.kill_fired(_sweep(small, small))
    assert not rs.kill_fired(_sweep(small, big))
    assert not rs.kill_fired(_sweep(big, big))
    # exactly half small is not a majority
    half = [0.005] * 5 + [0.05] * 5
    assert not rs.kill_fired(_sweep(half, half))


def test_summary_counts_seeds_at_or_below_kill_line():
    s = rs.summarize_gaps(_sweep([0.01, 0.02, 0.0], [0.1, 0.1, 0.1])).set_index("model")
    assert s.loc["logistic_l2", "seeds_gap_le_kill"] == 2
    assert s.loc["lightgbm", "seeds_gap_le_kill"] == 0


def test_within_sample_flag_changes_numbers_and_paths(monkeypatch, tmp_path):
    """A per-patient scale shift is invisible to ranking: the flag must change the AUCs."""
    X, y, cohort = _toy()
    X = X.mul(1 + 3 * (cohort == "A").to_numpy()[:, None] * np.random.default_rng(1).uniform(size=(len(y), 1)))
    X = X + 5 * (cohort == "B").to_numpy()[:, None]
    raw = pd.DataFrame(rs.sweep_rows(X, y, cohort, seed=200)).set_index("model")
    ws = pd.DataFrame(rs.sweep_rows(splits.rank_normalize_within_sample(X), y, cohort, seed=200)).set_index("model")
    assert not np.allclose(raw["loco_auc"], ws["loco_auc"])
    assert rs.OUT_SWEEP_WS != rs.OUT_SWEEP and rs.META_WS != rs.META and rs.OUT_SUMMARY_WS != rs.OUT_SUMMARY

    # the CLI path really transforms the matrix, and gates against the within-sample reference
    seen = {}
    monkeypatch.setattr(rs, "load_matrix", lambda c: (X, y, cohort))
    monkeypatch.setattr(rs, "repro_gate", lambda X_, y_, c_, ref="summary_raw.csv": seen.update(ref=ref, X=X_) or 1.0)
    assert rs.main(["--within-sample"]) == 1  # gate "fails" -> stops before writing
    assert seen["ref"] == "summary_withinsample.csv"
    assert np.allclose(seen["X"].to_numpy(), splits.rank_normalize_within_sample(X).to_numpy())
