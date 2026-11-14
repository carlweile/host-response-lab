"""Guards the copied Bodkin 2022 numbers. Must go red on drift."""
from __future__ import annotations

import csv
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
CSV = ROOT / "results" / "bodkin_table2_comparison.csv"
GAP = ROOT / "results" / "gapdiag_cohorts.csv"


def _rows(path: Path = CSV) -> list[dict]:
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def check_rows(rows: list[dict], gap: list[dict]) -> list[str]:
    """Return a list of problems; empty = fine."""
    probs = []
    if len(rows) != 5:
        probs.append(f"expected 5 rows, got {len(rows)}")
    ref = {r["cohort"]: r for r in gap if r["transform"] == "raw" and r["model"] == "logistic_l2"}
    for r in rows:
        s = r["series"]
        if r["in_bodkin"] == "yes" and not r["source"].startswith("Table S9, row "):
            probs.append(f"{s}: missing Table S9 row source")
        if not float(r["bodkin_q1"]) <= float(r["bodkin_median"]) <= float(r["bodkin_q3"]):
            probs.append(f"{s}: bacterial Q1<=median<=Q3 broken")
        if not float(r["viral_q1"]) <= float(r["viral_median"]) <= float(r["viral_q3"]):
            probs.append(f"{s}: viral Q1<=median<=Q3 broken")
        if abs(float(r["our_random_cv"]) - float(ref[s]["rcv_auc_in_cohort"])) > 1e-4:
            probs.append(f"{s}: our_random_cv drifted from gapdiag_cohorts.csv")
        if abs(float(r["our_loco"]) - float(ref[s]["loco_auc"])) > 1e-4:
            probs.append(f"{s}: our_loco drifted from gapdiag_cohorts.csv")
    return probs


def test_comparison_csv_is_clean():
    assert check_rows(_rows(), _rows(GAP)) == []


def test_guard_fires_on_corrupted_copy():
    rows = _rows()
    rows[0]["our_random_cv"] = "0.5"
    rows[1]["source"] = ""
    rows[2]["bodkin_q1"] = "0.99"
    probs = check_rows(rows, _rows(GAP))
    assert len(probs) == 3
    with pytest.raises(AssertionError):
        assert check_rows(rows[:4], _rows(GAP)) == []
