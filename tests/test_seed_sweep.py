"""Guards for the seed sweep (run_seed_sweep.py, the 10-seed raw sweep).

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


# --- the LightGBM model-seed check: LightGBM model-seed side check on the within-sample scale ---------------

def _ws_toy_fixture(monkeypatch, tmp_path):
    """Toy data with a per-patient scale shift, results/ redirected, gate reference = the toy's own seed-0 value."""
    pytest.importorskip("lightgbm")
    X, y, cohort = _toy()
    X = X.mul(1 + 3 * (cohort == "A").to_numpy()[:, None] * np.random.default_rng(1).uniform(size=(len(y), 1)))
    Xws = splits.rank_normalize_within_sample(X)
    ref = rs.lgbm_loco(Xws, y, cohort, 0)
    pd.DataFrame([dict(model="lightgbm", scheme="leave_one_cohort_out", mean_auc=ref)]).to_csv(
        tmp_path / "summary_withinsample.csv", index=False)
    pd.DataFrame(dict(model_seed=rs.SEEDS, model="lightgbm", loco_auc=0.80)).to_csv(tmp_path / "raw.csv", index=False)
    monkeypatch.setattr(rs.rp, "RESULTS", tmp_path)
    monkeypatch.setattr(rs, "OUT_MODEL", tmp_path / "raw.csv")
    monkeypatch.setattr(rs, "OUT_MODEL_WS", tmp_path / "new.csv")
    monkeypatch.setattr(rs, "META_MODEL_WS", tmp_path / "new.json")
    monkeypatch.setattr(rs, "load_matrix", lambda c: (X, y, cohort))
    return X, y, cohort, Xws


def test_modelseed_cli_uses_ranked_matrix_all_seeds_and_only_its_own_files(monkeypatch, tmp_path):
    X, y, cohort, Xws = _ws_toy_fixture(monkeypatch, tmp_path)
    raw_before = (tmp_path / "raw.csv").read_bytes()
    seen = []
    real = train.run_both_schemes

    def spy(X_, y_, c_, seed=0, model_seed=0):
        seen.append((model_seed, np.allclose(X_.to_numpy(), Xws.to_numpy())))
        return real(X_, y_, c_, seed=seed, model_seed=model_seed)

    monkeypatch.setattr(train, "run_both_schemes", spy)
    assert rs.main(["--within-sample", "--lgbm-modelseed"]) == 0
    assert {m for m, _ in seen} == {0, *rs.SEEDS}  # a dropped or pinned-to-0 seed goes red
    assert all(ranked for _, ranked in seen)  # falling back to the raw matrix goes red
    out = pd.read_csv(tmp_path / "new.csv")
    assert list(out.columns) == ["model_seed", "model", "loco_auc"] and len(out) == 10
    assert out["loco_auc"].nunique() > 1  # a model seed that is ignored gives ten identical rows
    assert (tmp_path / "raw.csv").read_bytes() == raw_before  # the raw side check is never rewritten
    assert sorted(p.name for p in tmp_path.iterdir()) == ["new.csv", "new.json", "raw.csv", "summary_withinsample.csv"]


def test_modelseed_gate_failure_writes_nothing(monkeypatch, tmp_path):
    _ws_toy_fixture(monkeypatch, tmp_path)
    pd.DataFrame([dict(model="lightgbm", scheme="leave_one_cohort_out", mean_auc=0.5)]).to_csv(
        tmp_path / "summary_withinsample.csv", index=False)
    assert rs.main(["--within-sample", "--lgbm-modelseed"]) == 1
    assert not (tmp_path / "new.csv").exists() and not (tmp_path / "new.json").exists()


def test_modelseed_option_refused_without_within_sample():
    with pytest.raises(SystemExit):
        rs.main(["--lgbm-modelseed"])


def test_lgbm_loco_equals_inline_mean_and_moves_with_model_seed():
    pytest.importorskip("lightgbm")
    X, y, cohort = _toy()
    res = train.run_both_schemes(X, y, cohort, seed=0, model_seed=203)
    aucs = [r.auc for r in res if r.model == "lightgbm" and r.scheme == "leave_one_cohort_out"]
    assert rs.lgbm_loco(X, y, cohort, 203) == pytest.approx(sum(aucs) / len(aucs), abs=1e-15)


def _pair(tmp_path, raw, ws):
    pd.DataFrame(dict(model_seed=range(10), model="lightgbm", loco_auc=raw)).to_csv(tmp_path / "r.csv", index=False)
    pd.DataFrame(dict(model_seed=range(10), model="lightgbm", loco_auc=ws)).to_csv(tmp_path / "w.csv", index=False)
    return rs.compare_modelseed(tmp_path / "r.csv", tmp_path / "w.csv", 0.888)


def test_reading_b_fires_only_when_paired_difference_is_inside_two_sds(tmp_path):
    noise = np.array([-.03, .03, -.02, .02, -.01, .01, -.04, .04, 0, 0])
    raw = 0.80 + noise
    clear = _pair(tmp_path, raw, raw + 0.06 + noise * 0.05)  # same shift every seed, tiny wobble
    assert not clear["reading_b"] and clear["paired_diff"]["mean"] == pytest.approx(0.06, abs=1e-3)
    muddy = _pair(tmp_path, raw, raw + 0.06 + noise * 3)  # shift buried in the pair spread
    assert muddy["reading_b"]
    # exactly on the line: mean == 2 x SD is NOT "smaller than", so (b) does not fire
    d = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    d[:] = np.linspace(-1, 1, 10)
    d = d / d.std(ddof=1) * 0.04 + 0.08  # SD 0.04, mean 0.08
    edge = _pair(tmp_path, raw, raw + d)
    assert edge["paired_diff"]["mean"] == pytest.approx(2 * edge["paired_diff"]["sd"])
    assert edge["reading_b"] == bool(edge["paired_diff"]["mean"] < 2 * edge["paired_diff"]["sd"])
    below = _pair(tmp_path, raw, raw + d - 0.001)
    assert below["reading_b"]


def test_reading_a_line_is_a_strict_0_01_shift_of_the_mean(tmp_path):
    raw = np.full(10, 0.80)
    assert not _pair(tmp_path, raw, np.full(10, 0.888 + 0.0099))["reading_a"]
    assert _pair(tmp_path, raw, np.full(10, 0.888 + 0.0101))["reading_a"]
    assert _pair(tmp_path, raw, np.full(10, 0.888 - 0.0101))["reading_a"]
