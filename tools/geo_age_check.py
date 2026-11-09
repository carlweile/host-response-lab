"""Read-only check: does GEO's series matrix carry a per-sample age?

Reads the cached series matrices only (no download, no pipeline run) and writes
results/geo_age_check.csv, one row per cohort. Run from the repo root:
    python tools/geo_age_check.py
"""
from __future__ import annotations

import csv
import json
import re
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src import parse_geo  # noqa: E402

# Key contains "age" as its own word ("age", "age_years", "patient age"), not "stage"/"passage".
KEY_RE = re.compile(r"(?<![a-z])age(?![a-z])", re.I)
NUM_RE = re.compile(r"(\d+(?:\.\d+)?)\s*([a-z]*)", re.I)
MONTH_KEY_RE = re.compile(r"month", re.I)
COHORTS = ["GSE63990", "GSE42026", "GSE40396", "GSE6269", "GSE60244"]


def parse_age(strings: list[str]) -> tuple[str | None, float | None, bool]:
    """Return (key_seen, age_in_years, unparseable); months become years from one sample's strings.

    key_seen is None if no string has an age-like key. unparseable is True when
    a key exists but holds no number ("age: N/A").
    """
    for s in strings:
        key, sep, val = s.partition(":")
        if not sep or not KEY_RE.search(key):
            continue
        m = NUM_RE.search(val)
        if not m:
            return key.strip(), None, True
        age = float(m.group(1))
        unit = m.group(2).lower()
        if MONTH_KEY_RE.search(key) or unit in ("m", "mo", "month", "months"):
            age /= 12  # "age (month): 4", "Age:4m, when sample taken"
        return key.strip(), age, False
    return None, None, False


def check(acc: str, cfg_group: str) -> tuple[dict, set[str]]:
    path = ROOT / "data" / "series" / f"{acc}_series_matrix.txt.gz"
    series = parse_geo.parse_series_matrix(path, acc)
    keys: set[str] = set()
    ages: list[float] = []
    unparse = 0
    for sample in series.samples:
        key, age, bad = parse_age(series.characteristics[sample])
        if key:
            keys.add(key)
        if age is not None:
            ages.append(age)
        unparse += bad
    row = {"cohort": acc, "pinned_age_group": cfg_group, "n_samples": len(series.samples),
           "n_with_age": len(ages), "n_unparseable": unparse}
    if ages:
        row.update(min=min(ages), median=statistics.median(ages), max=max(ages),
                   n_lt16=sum(a < 16 for a in ages), n_lt18=sum(a < 18 for a in ages),
                   n_ge18=sum(a >= 18 for a in ages), n_lt21=sum(a < 21 for a in ages))
    return row, keys


def main() -> None:
    cfg = json.loads((ROOT / "config" / "cohorts.json").read_text(encoding="utf-8"))
    groups = {c["accession"]: c.get("age_group", "") for c in cfg["cohorts"]}
    fields = ["cohort", "pinned_age_group", "n_samples", "n_with_age", "n_unparseable",
              "min", "median", "max", "n_lt16", "n_lt18", "n_ge18", "n_lt21", "keys_seen"]
    rows = []
    for acc in COHORTS:
        row, keys = check(acc, groups.get(acc, ""))
        row["keys_seen"] = "|".join(sorted(keys))
        rows.append(row)
        print({k: row.get(k, "") for k in fields})
    out = ROOT / "results" / "geo_age_check.csv"
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print("wrote", out)


if __name__ == "__main__":
    main()
