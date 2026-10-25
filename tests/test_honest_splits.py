"""These tests must FAIL on broken code, not just pass on correct code.

Each one deliberately constructs the bug it is guarding against and asserts that
the guard fires. A test that only ever sees correct input proves nothing.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import labels as lab  # noqa: E402
from src.parse_geo import Series  # noqa: E402
from src.splits import assert_no_cohort_overlap, cohort_splits, random_splits  # noqa: E402


def _cohorts() -> pd.Series:
    return pd.Series(["A"] * 10 + ["B"] * 10 + ["C"] * 10)


def test_cohort_splits_never_share_a_cohort():
    cohort = _cohorts()
    seen = []
    for held_out, split in cohort_splits(cohort):
        assert_no_cohort_overlap(cohort, split)
        assert set(cohort.iloc[split[1]]) == {held_out}
        seen.append(held_out)
    assert sorted(seen) == ["A", "B", "C"]


def test_guard_actually_catches_a_leaky_split():
    """The break-it case: hand the guard a split that leaks and demand it raises."""
    cohort = _cohorts()
    leaky = (np.arange(0, 25), np.arange(20, 30))  # cohort C on both sides
    with pytest.raises(AssertionError, match="cohort leakage"):
        assert_no_cohort_overlap(cohort, leaky)


def test_random_splits_do_leak_cohorts_which_is_the_whole_point():
    """Confirms the two schemes really differ, so the reported gap is meaningful."""
    cohort = _cohorts()
    y = pd.Series([0, 1] * 15)
    leaked = 0
    for _fold, split in random_splits(y, n_splits=5):
        if set(cohort.iloc[split[0]]) & set(cohort.iloc[split[1]]):
            leaked += 1
    assert leaked == 5, "random CV should mix cohorts; if not, the comparison is void"


def _series_with(chars: list[str]) -> Series:
    samples = [f"GSM{i}" for i in range(len(chars))]
    expr = pd.DataFrame(np.ones((3, len(chars))), columns=samples)
    return Series("GSEtest", "GPL1", expr, {s: [c] for s, c in zip(samples, chars)})


RULES = {
    "order": ["exclude", "bacterial", "viral"],
    "exclude": ["co-infection", "healthy"],
    "bacterial": ["bacterial"],
    "viral": ["viral", "influenza"],
}


def test_coinfection_is_excluded_not_counted_as_bacterial():
    """The bug this guards: 'bacterial and viral co-infection' contains 'bacterial'."""
    series = _series_with([
        "infection: bacterial",
        "infection: viral",
        "infection: bacterial and viral co-infection",
        "status: healthy control",
    ])
    result = lab.assign_labels(series, RULES, expected=None)
    counts = result.counts
    assert counts[lab.BACTERIAL] == 1
    assert counts[lab.VIRAL] == 1
    assert counts[lab.EXCLUDED] == 2
    assert len(result.usable) == 2


def test_wrong_rule_order_produces_the_bug_so_order_is_load_bearing():
    bad_order = dict(RULES, order=["bacterial", "viral", "exclude"])
    series = _series_with(["infection: bacterial and viral co-infection"])
    result = lab.assign_labels(series, bad_order, expected=None)
    assert result.counts[lab.BACTERIAL] == 1, "demonstrates why exclusions run first"


def test_count_mismatch_is_reported_as_unverified():
    series = _series_with(["infection: bacterial", "infection: viral"])
    ok = lab.assign_labels(series, RULES, expected={"bacterial": 1, "viral": 1})
    assert ok.verified and "verified" in ok.verdict

    bad = lab.assign_labels(series, RULES, expected={"bacterial": 73, "viral": 117})
    assert not bad.verified
    assert "MISMATCH" in bad.verdict

    none = lab.assign_labels(series, RULES, expected=None)
    assert not none.verified and "UNVERIFIED" in none.verdict


def _real_rules() -> dict:
    cfg = Path(__file__).resolve().parents[1] / "config" / "cohorts.json"
    import json
    return json.loads(cfg.read_text(encoding="utf-8"))["label_rules"]


def test_afebrile_virus_positive_child_is_excluded_by_the_real_rules():
    """The bug this guards: GSE40396 lists afebrile virus carriers with a pathogen
    name, so a name-only rule counts them as viral infections (35 instead of 22)."""
    series = _series_with([
        "fever: no | pathogen: rhinovirus | group: rhinovirus",
        "fever: yes | pathogen: adenovirus | group: adenovirus",
    ])
    result = lab.assign_labels(series, _real_rules(), expected=None)
    assert result.labels["GSM0"] == lab.EXCLUDED
    assert result.labels["GSM1"] == lab.VIRAL


def test_without_the_fever_rule_the_afebrile_child_becomes_viral():
    rules = _real_rules()
    rules = dict(rules, exclude=[p for p in rules["exclude"] if p != "fever: no"])
    series = _series_with(["fever: no | pathogen: rhinovirus | group: rhinovirus"])
    result = lab.assign_labels(series, rules, expected=None)
    assert result.labels["GSM0"] == lab.VIRAL, "demonstrates why the fever rule is needed"
