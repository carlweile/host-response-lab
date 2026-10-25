"""Exercise the harmonize -> train path on synthetic cohorts.

The synthetic data is built with a deliberate cohort offset on top of a weak
real signal, so the expected behaviour is a visible drop from random CV to
leave-one-cohort-out. If that drop ever disappears here, the measurement
machinery is broken, not the biology.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.harmonize import build_matrix, maybe_log2, to_symbols  # noqa: E402
from src.train import run_both_schemes, summarize  # noqa: E402

GENES = [f"GENE{i}" for i in range(40)]


def _fake_cohort(name: str, n: int, offset: float, rng) -> tuple[pd.DataFrame, dict[str, str]]:
    labels, cols, data = {}, [], []
    for i in range(n):
        sample = f"{name}_S{i}"
        is_bact = i % 2 == 0
        vec = rng.normal(offset, 1.0, size=len(GENES))
        vec[:4] += 1.1 if is_bact else -1.1          # the real, shared signal
        vec[4:10] += offset * 3.0                    # cohort-specific quirk
        cols.append(sample)
        data.append(vec)
        labels[sample] = "bacterial" if is_bact else "viral"
    df = pd.DataFrame(np.array(data).T, index=pd.Index(GENES, name="symbol"), columns=cols)
    return df, labels


def _make(n_cohorts: int = 3, n: int = 40):
    rng = np.random.default_rng(0)
    per_cohort, labels = {}, {}
    for k in range(n_cohorts):
        name = f"GSEfake{k}"
        df, lab = _fake_cohort(name, n, offset=float(k), rng=rng)
        per_cohort[name], labels[name] = df, lab
    return per_cohort, labels


def test_build_matrix_shapes_and_grouping():
    per_cohort, labels = _make()
    X, y, cohort, notes = build_matrix(per_cohort, labels)
    assert X.shape == (120, 40)
    assert set(cohort.unique()) == {"GSEfake0", "GSEfake1", "GSEfake2"}
    assert y.sum() == 60
    assert X.isna().sum().sum() == 0
    assert any("shared gene symbols" in n for n in notes)


def test_both_schemes_run_and_random_cv_looks_better():
    per_cohort, labels = _make()
    X, y, cohort, _ = build_matrix(per_cohort, labels)
    results = run_both_schemes(X, y, cohort, n_splits=4)
    summary = summarize(results)

    assert not summary.empty
    assert set(summary["scheme"]) == {"random_cv", "leave_one_cohort_out"}
    for row in summary.itertuples():
        assert 0.0 <= row.mean_auc <= 1.0

    rnd = summary[summary["scheme"] == "random_cv"]["mean_auc"].mean()
    loco = summary[summary["scheme"] == "leave_one_cohort_out"]["mean_auc"].mean()
    assert rnd >= loco, f"synthetic cohort offset should hurt LOCO: {rnd:.3f} vs {loco:.3f}"
    # a scorer stuck at chance (e.g. a constant 0.5) must fail here (review 2026-10-21)
    assert rnd > 0.8, f"random CV should find the shared signal: {rnd:.3f}"


def test_every_fold_reports_its_sizes():
    per_cohort, labels = _make(n_cohorts=2, n=30)
    X, y, cohort, _ = build_matrix(per_cohort, labels)
    for r in run_both_schemes(X, y, cohort, n_splits=3):
        assert r.n_train > 0 and r.n_test > 0
        assert r.auc is None or 0.0 <= r.auc <= 1.0


def test_log2_only_applied_to_linear_scale_data():
    linear = pd.DataFrame(np.full((3, 3), 5000.0))
    logged, did = maybe_log2(linear)
    assert did and logged.to_numpy().max() < 20

    already = pd.DataFrame(np.full((3, 3), 8.0))
    out, did2 = maybe_log2(already)
    assert not did2 and out.equals(already)


def test_to_symbols_keeps_highest_mean_probe_per_gene():
    expr = pd.DataFrame(
        {"S1": [1.0, 9.0, 4.0], "S2": [1.0, 9.0, 4.0]},
        index=["p_low", "p_high", "p_other"],
    )
    mapping = pd.Series({"p_low": "A", "p_high": "A", "p_other": "B"})
    out = to_symbols(expr, mapping)
    assert sorted(out.index) == ["A", "B"]
    assert out.loc["A", "S1"] == 9.0, "should keep the higher-mean probe for gene A"
