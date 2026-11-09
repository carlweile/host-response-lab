import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from geo_age_check import parse_age  # noqa: E402


def test_plain_age():
    assert parse_age(["gender: F", "age: 45"]) == ("age", 45.0, False)


def test_age_with_unit_in_key():
    assert parse_age(["age (years): 7"]) == ("age (years)", 7.0, False)


def test_months_become_years():
    assert parse_age(["age (month): 6"])[1] == 0.5
    assert parse_age(["Age:4m, when sample taken"])[1] == 4 / 12
    assert parse_age(["Age: 10 y, when sample taken"])[1] == 10.0


def test_na_is_unparseable_and_missing_is_none():
    assert parse_age(["age: N/A"]) == ("age", None, True)
    assert parse_age(["stage: 3", "batch: 0"]) == (None, None, False)  # "stage" is not "age"
