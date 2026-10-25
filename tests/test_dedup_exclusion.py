"""GSE63990 re-run arrays (LABBOOK 2026-10-21). Each test builds the bug it guards.

The bug: GEO holds 280 arrays of 273 subjects. Without the exclusion the 7 re-runs
count as extra patients (73/117 instead of the published 70/115), and a patient
can sit on both sides of a random-CV split.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import labels as lab  # noqa: E402
from src.parse_geo import Series  # noqa: E402

CONFIG = Path(__file__).resolve().parents[1] / "config" / "cohorts.json"
RERUNS = [f"GSM{n}" for n in range(1562133, 1562140)]  # sample_274..sample_280


def _gse63990_entry() -> dict:
    cfg = json.loads(CONFIG.read_text(encoding="utf-8"))
    return next(c for c in cfg["cohorts"] if c["accession"] == "GSE63990")


def _rules() -> dict:
    return json.loads(CONFIG.read_text(encoding="utf-8"))["label_rules"]


def _series_like_gse63990() -> Series:
    """Published subjects 70 bacterial + 115 viral + 88 non-infectious, plus the 7
    re-runs carrying the surplus (3 bacterial, 2 viral, 2 non-infectious)."""
    chars = (["infection_status: bacterial"] * 70 + ["infection_status: viral"] * 115
             + ["infection_status: non-infectious illness"] * 88)
    rerun_chars = (["infection_status: bacterial"] * 3 + ["infection_status: viral"] * 2
                   + ["infection_status: non-infectious illness"] * 2)
    samples = [f"GSM{1000 + i}" for i in range(len(chars))] + RERUNS
    chars += rerun_chars
    expr = pd.DataFrame(np.ones((2, len(samples))), columns=samples)
    return Series("GSE63990", "GPL571", expr, {s: [c] for s, c in zip(samples, chars)})


def test_config_excludes_exactly_the_seven_reruns_and_expects_the_paper():
    entry = _gse63990_entry()
    assert sorted(entry["exclude_samples"]["ids"]) == RERUNS
    assert entry["expected_counts"] == {"bacterial": 70, "viral": 115}


def test_without_the_exclusion_the_reruns_come_back_as_a_mismatch():
    entry = _gse63990_entry()
    res = lab.assign_labels(_series_like_gse63990(), _rules(), entry["expected_counts"])
    assert res.counts[lab.BACTERIAL] == 73 and res.counts[lab.VIRAL] == 117
    assert res.verdict.startswith("MISMATCH"), "demonstrates why the exclusion is needed"
    with pytest.raises(lab.LabelMismatch):
        lab.refuse_mismatches([res])


def test_with_the_config_exclusion_the_parse_matches_the_paper():
    entry = _gse63990_entry()
    res = lab.assign_labels(_series_like_gse63990(), _rules(), entry["expected_counts"],
                            entry["exclude_samples"]["ids"])
    assert res.verified, res.verdict
    assert res.counts[lab.EXCLUDED] == 88 + 7
    assert not set(RERUNS) & set(res.usable)


def test_an_exclusion_id_not_in_the_series_is_refused():
    with pytest.raises(ValueError, match="not in series"):
        lab.assign_labels(_series_like_gse63990(), _rules(), None, ["GSM9999999"])
