"""Guards added after the three outside code reviews of 2026-10-21
(reports/review/2026-10-21/). Each test builds the bug a reviewer described and
asserts that the guard fires.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import fetch_geo, train  # noqa: E402
from src import labels as lab  # noqa: E402
from src.harmonize import build_matrix, to_symbols  # noqa: E402
from src.parse_geo import Series  # noqa: E402

RULES = json.loads((Path(__file__).resolve().parents[1] / "config" / "cohorts.json")
                   .read_text(encoding="utf-8"))["label_rules"]


def _series(chars: dict[str, str]) -> Series:
    expr = pd.DataFrame(np.zeros((2, len(chars))), columns=list(chars))
    return Series("GSEx", "GPLx", expr, {s: [c] for s, c in chars.items()})


# --- labels: a count mismatch must stop training (all three reviewers) ----------

def test_mismatch_is_refused():
    res = lab.assign_labels(_series({"a": "infection: bacterial", "b": "infection: bacterial",
                                     "c": "infection: viral"}), RULES, {"bacterial": 1, "viral": 1})
    assert res.verdict.startswith("MISMATCH")
    with pytest.raises(lab.LabelMismatch):
        lab.refuse_mismatches([res])


def test_verified_and_unverified_pass_the_refusal():
    ok = lab.assign_labels(_series({"a": "infection: bacterial", "b": "infection: viral"}),
                           RULES, {"bacterial": 1, "viral": 1})
    unverified = lab.assign_labels(_series({"a": "infection: viral"}), RULES, None)
    lab.refuse_mismatches([ok, unverified])


def test_haemophilus_influenzae_is_bacterial_not_influenza():
    res = lab.assign_labels(_series({"a": "pathogen: haemophilus influenzae"}), RULES, None)
    assert res.labels["a"] == lab.BACTERIAL


# --- platform: GPL96 is a requirement, not a preference (all three) -------------

def _fake_fetch(monkeypatch, variants):
    calls = []
    monkeypatch.setattr(fetch_geo, "list_matrix_files", lambda acc: variants)
    monkeypatch.setattr(fetch_geo, "download",
                        lambda url, dest, force=False: (calls.append(url), (True, "ok"))[1])
    return calls


def test_preferred_platform_never_takes_the_combined_matrix(monkeypatch, tmp_path):
    calls = _fake_fetch(monkeypatch, ["GSE6269-GPL2507_series_matrix.txt.gz",
                                      "GSE6269-GPL96_series_matrix.txt.gz"])
    path, _ = fetch_geo.fetch_series("GSE6269", tmp_path, prefer_platform="GPL96")
    assert path is not None and len(calls) == 1 and "GPL96_" in calls[0]


def test_missing_preferred_platform_is_refused_not_replaced(monkeypatch, tmp_path):
    calls = _fake_fetch(monkeypatch, ["GSE6269-GPL2507_series_matrix.txt.gz",
                                      "GSE6269-GPL961_series_matrix.txt.gz"])
    path, msg = fetch_geo.fetch_series("GSE6269", tmp_path, prefer_platform="GPL96")
    assert path is None and calls == [] and "no GPL96" in msg


# --- training path: rows must line up (astra 8) ---------------------------------

def test_misaligned_rows_are_refused():
    idx = [f"s{i}" for i in range(20)]
    X = pd.DataFrame(np.random.default_rng(0).normal(size=(20, 3)), index=idx)
    y = pd.Series([0, 1] * 10, index=idx)
    cohort = pd.Series(["A"] * 10 + ["B"] * 10, index=idx)
    with pytest.raises(ValueError, match="one index"):
        train.run_both_schemes(X, y.iloc[::-1], cohort)
    with pytest.raises(ValueError, match="one index"):
        train.run_both_schemes(X, y, cohort.iloc[::-1])


# --- harmonisation: imputation is on record; probe choice can avoid the data ----

def test_imputed_cell_count_is_reported():
    a = pd.DataFrame({"s1": [1.0, 2.0], "s2": [3.0, np.nan]}, index=["G1", "G2"])
    b = pd.DataFrame({"t1": [1.0, 2.0], "t2": [3.0, 4.0]}, index=["G1", "G2"])
    labels = {"A": {"s1": "viral", "s2": "bacterial"}, "B": {"t1": "viral", "t2": "bacterial"}}
    _, _, _, notes = build_matrix({"A": a, "B": b}, labels, min_genes=1)
    assert "imputed cells (cohort gene mean, else 0.0): 1" in notes


def test_first_collapse_reads_no_expression_value():
    mapping = pd.Series({"p2": "A", "p1": "A", "p3": "B"})
    expr = pd.DataFrame({"S1": [9.0, 1.0, 4.0]}, index=["p2", "p1", "p3"])
    flipped = pd.DataFrame({"S1": [1.0, 9.0, 4.0]}, index=["p2", "p1", "p3"])
    assert to_symbols(expr, mapping, "first").loc["A", "S1"] == 1.0   # p1 = lowest id
    assert to_symbols(flipped, mapping, "first").loc["A", "S1"] == 9.0
    with pytest.raises(ValueError):
        to_symbols(expr, mapping, "median")
