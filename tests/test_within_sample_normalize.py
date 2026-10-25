"""The within-sample transform must read one patient at a time.

The guard (`_reads_other_patients`) swaps every other row of the test cohort for
noise and checks a given row's output is unchanged. Test (b) builds the leaky
cohort-wise transform and demands the same guard fires on it.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.splits import rank_normalize_within_cohort, rank_normalize_within_sample  # noqa: E402


def _data() -> tuple[pd.DataFrame, pd.Series]:
    rng = np.random.default_rng(0)
    X = pd.DataFrame(rng.normal(size=(12, 30)), columns=[f"g{i}" for i in range(30)])
    cohort = pd.Series(["A"] * 6 + ["B"] * 6)
    return X, cohort


def _reads_other_patients(transform, row: int = 8) -> bool:
    """True if row `row`'s output changes when all other rows of its cohort change."""
    X, cohort = _data()
    before = transform(X, cohort).iloc[row].to_numpy()
    X2 = X.copy()
    rng = np.random.default_rng(99)
    others = [i for i in range(6, 12) if i != row]
    X2.iloc[others] = rng.normal(loc=5.0, scale=3.0, size=(len(others), X.shape[1]))
    after = transform(X2, cohort).iloc[row].to_numpy()
    return not np.allclose(before, after)


def test_within_sample_reads_no_other_patient():
    assert not _reads_other_patients(lambda X, c: rank_normalize_within_sample(X))


def test_guard_fires_on_the_leaky_cohort_transform():
    """Break-it case: the cohort-wise transform reads the rest of the test cohort."""
    assert _reads_other_patients(rank_normalize_within_cohort)


def test_values_in_unit_interval_and_order_preserved():
    X, _ = _data()
    out = rank_normalize_within_sample(X)
    assert out.to_numpy().min() > 0 and out.to_numpy().max() <= 1.0
    for i in range(len(X)):
        assert (np.argsort(X.iloc[i].to_numpy()) == np.argsort(out.iloc[i].to_numpy())).all()


def test_nan_filled_with_half():
    X, _ = _data()
    X.iloc[0, 0] = np.nan
    out = rank_normalize_within_sample(X)
    assert not out.isna().any().any()
    assert out.iloc[0, 0] == 0.5  # the documented fill, not just "no NaN" (review 2026-10-21)
