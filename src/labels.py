"""Turn messy GEO characteristic strings into bacterial/viral labels, and check
the result against counts published for the same cohort.

Rule order matters and is not an implementation detail: "bacterial and viral
co-infection" contains the word "bacterial", so exclusions must be tested first
or every co-infected patient silently becomes a bacterial training example.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from .parse_geo import Series

BACTERIAL, VIRAL, EXCLUDED = "bacterial", "viral", "excluded"


@dataclass
class LabelResult:
    accession: str
    labels: dict[str, str]              # sample -> bacterial | viral | excluded
    counts: Counter
    verified: bool
    verdict: str
    distinct_strings: list[tuple[str, int]]

    @property
    def usable(self) -> dict[str, str]:
        return {s: l for s, l in self.labels.items() if l in (BACTERIAL, VIRAL)}


def _classify(blob: str, rules: dict) -> str:
    for cls in rules["order"]:
        for pattern in rules[cls]:
            if pattern in blob:
                return EXCLUDED if cls == "exclude" else cls
    return EXCLUDED


def assign_labels(series: Series, rules: dict, expected: dict | None,
                  exclude_samples: list[str] | None = None) -> LabelResult:
    """`exclude_samples`: sample ids excluded before any rule runs, e.g. re-run
    arrays of a patient already in the cohort (GSE63990, LABBOOK 2026-10-21)."""
    labels: dict[str, str] = {}
    seen: Counter = Counter()
    drop = set(exclude_samples or ())
    missing = drop - set(series.characteristics)
    if missing:  # a typo would otherwise exclude nothing, silently
        raise ValueError(f"{series.accession}: exclude_samples not in series: {sorted(missing)}")

    for sample, strings in series.characteristics.items():
        blob = " | ".join(strings).lower()
        seen[blob] += 1
        labels[sample] = EXCLUDED if sample in drop else _classify(blob, rules)

    counts = Counter(labels.values())
    verified, verdict = _verify(counts, expected)
    distinct = sorted(seen.items(), key=lambda kv: -kv[1])
    return LabelResult(series.accession, labels, counts, verified, verdict, distinct)


def _verify(counts: Counter, expected: dict | None) -> tuple[bool, str]:
    if not expected:
        return False, "UNVERIFIED: no published counts on file to check against"
    problems = []
    for cls, want in expected.items():
        got = counts.get(cls, 0)
        if got != want:
            problems.append(f"{cls}: parsed {got}, published {want}")
    if problems:
        return False, "MISMATCH: " + "; ".join(problems)
    return True, "verified against published counts"


class LabelMismatch(RuntimeError):
    """A cohort's parse disagrees with its published counts: the parse is wrong."""


def refuse_mismatches(results: list[LabelResult]) -> None:
    """Hard stop before training (review finding 2026-10-21: a MISMATCH used to be
    only a note). Unverified cohorts (no published counts) are still allowed."""
    bad = [f"{r.accession}: {r.verdict}" for r in results if r.verdict.startswith("MISMATCH")]
    if bad:
        raise LabelMismatch("refusing to train on labels that contradict published "
                            "counts: " + " | ".join(bad))


def review_text(result: LabelResult, max_rows: int = 14) -> str:
    lines = [
        f"--- {result.accession} ---",
        f"  parsed: " + ", ".join(f"{k}={v}" for k, v in sorted(result.counts.items())),
        f"  check : {result.verdict}",
        "  distinct characteristic strings (count x string):",
    ]
    for blob, n in result.distinct_strings[:max_rows]:
        short = blob if len(blob) <= 110 else blob[:107] + "..."
        lines.append(f"    {n:4d} x {short}")
    hidden = len(result.distinct_strings) - max_rows
    if hidden > 0:
        lines.append(f"    ... {hidden} more distinct strings")
    return "\n".join(lines)
