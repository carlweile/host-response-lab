"""Guards for run_whole_blood_check.py. Each test builds the bug it guards against."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import run_whole_blood_check as wb  # noqa: E402


def _toy(n_per=12, genes=6, seed=0):
    rng = np.random.default_rng(seed)
    cohorts = ["A", "B", "C", wb.PBMC_COHORT]
    idx = [f"{c}:{i}" for c in cohorts for i in range(n_per)]
    X = pd.DataFrame(rng.normal(size=(len(idx), genes)), index=idx,
                     columns=[f"g{j}" for j in range(genes)])
    y = pd.Series(rng.integers(0, 2, len(idx)), index=idx)
    cohort = pd.Series([c for c in cohorts for _ in range(n_per)], index=idx)
    return X, y, cohort


def test_pbmc_cohort_absent_from_train_and_test():
    X, y, c = _toy()
    X4, y4, c4 = wb.drop_cohort(X, y, c, wb.PBMC_COHORT)
    assert (c4 != wb.PBMC_COHORT).all() and len(X4) == len(y4) == len(c4) == 36
    for _, (tr, te) in wb.splits.cohort_splits(c4):
        assert wb.PBMC_COHORT not in set(c4.iloc[list(tr)]) | set(c4.iloc[list(te)])
    # the bug: forgetting to drop it -> the guard fires
    with pytest.raises(AssertionError):
        wb.assert_cohort_absent(c, wb.PBMC_COHORT)


def test_gene_set_stays_fixed_when_a_cohort_is_dropped():
    X, y, c = _toy()
    X4, _, _ = wb.drop_cohort(X, y, c, wb.PBMC_COHORT)
    assert list(X4.columns) == list(X.columns)
    # the bug: re-intersecting genes after the drop would change the columns
    assert list(X4.drop(columns=["g0"]).columns) != list(X.columns)


def test_gate0_detects_a_moved_number_and_missing_rows():
    X, y, c = _toy(n_per=20, genes=5, seed=1)
    fresh = wb.loco_rows(X, y, c)
    ref = fresh.rename(columns={"held_out": "fold"}).assign(scheme="leave_one_cohort_out")
    assert wb.gate0(X, y, c, ref) < 1e-9
    moved = ref.copy()
    moved.loc[moved.index[0], "auc"] += 0.05
    assert wb.gate0(X, y, c, moved) > wb.GATE0_TOL
    with pytest.raises(AssertionError):
        wb.gate0(X, y, c, ref.iloc[:-1])


def test_committed_result_files_agree_with_each_other():
    res = Path(__file__).resolve().parents[1] / "results"
    csv, js = res / "whole_blood_check.csv", res / "run_whole_blood_check.json"
    if not (csv.exists() and js.exists()):
        pytest.skip("whole-blood check not run")
    import json
    d = pd.read_csv(csv)
    meta = json.loads(js.read_text(encoding="utf-8"))
    assert wb.PBMC_COHORT not in set(d["fold"])
    four = d[(d.scheme == "leave_one_cohort_out") & (d.model == "logistic_l2")
             & (d.model_seed == 0) & (d.gene_set == "fixed_5cohort")]["auc"].mean()
    assert abs(four - meta["comparison"]["logistic_l2"]["loco_4cohort_run"]) < 1e-9
