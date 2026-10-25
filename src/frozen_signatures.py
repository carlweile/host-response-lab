"""Score FROZEN published signatures on our cohorts, both validation schemes.

Nothing is fitted. A patient's score is computed from that patient's own genes
only (no labels, no other patient), with the sign the paper published: higher =
more bacterial. The AUC is never flipped to make a weak score look better.

"Random CV" here uses the pipeline's own stratified folds over all scorable
patients pooled across cohorts; with no training, each test fold is just a
pooled sample. "Leave-one-cohort-out" is the AUC inside each cohort. So the gap
measures what pooling cohorts does to a fixed score, not overfitting.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from .splits import random_splits

# Copied from the papers (see config/signatures.json for the quotes). The config is
# checked against these, so an edited gene list cannot slip through unnoticed.
PUBLISHED = {
    "herberg2016_drs": ({"FAM89A"}, {"IFI44L"}, "sum"),
    "sweeney2016_bv7": ({"HK3", "TNIP1", "GPAA1", "CTSB"}, {"IFI27", "JUP", "LAX1"}, "geomean"),
}


@dataclass(frozen=True)
class Signature:
    name: str
    bacterial_up: tuple[str, ...]
    viral_up: tuple[str, ...]
    rule: str
    independent: tuple[str, ...]  # our cohorts NOT used to choose its genes
    discovery: tuple[str, ...]    # our cohorts its genes were chosen on

    @property
    def genes(self) -> tuple[str, ...]:
        return self.bacterial_up + self.viral_up


def load_signatures(path: Path) -> list[Signature]:
    raw = json.loads(Path(path).read_text(encoding="utf-8"))["signatures"]
    sigs = []
    for s in raw:
        sig = Signature(s["name"], tuple(s["genes_bacterial_up"]), tuple(s["genes_viral_up"]),
                        s["rule"], tuple(s["validation_cohorts_ours"]),
                        tuple(s["discovery_cohorts_ours"]))
        check_published(sig)
        sigs.append(sig)
    return sigs


def check_published(sig: Signature) -> None:
    """Refuse a gene list, direction or rule that differs from the paper."""
    if sig.name not in PUBLISHED:
        raise ValueError(f"{sig.name}: not a pinned published signature")
    bact, vir, rule = PUBLISHED[sig.name]
    if set(sig.bacterial_up) != bact or set(sig.viral_up) != vir or sig.rule != rule:
        raise ValueError(f"{sig.name}: config differs from the published gene list / rule")


def _geomean(block: pd.DataFrame) -> pd.Series:
    """Row-wise geometric mean. A 0 makes it 0 (literal rule); negatives are refused."""
    if (block < 0).any().any():
        raise ValueError("geometric mean of negative expression values is undefined")
    with np.errstate(divide="ignore"):
        return np.exp(np.log(block).mean(axis=1))


def ratio_weight(sig: Signature) -> float:
    """Sweeney's 'times the ratio of the number of genes in each set', read as a factor
    on the bacterial term only, in the order the text names the sets: n_viral / n_bacterial."""
    return len(sig.viral_up) / len(sig.bacterial_up)


def score(genes_x_samples: pd.DataFrame, sig: Signature, bacterial_weight: float = 1.0) -> pd.Series:
    """One cohort's genes x samples matrix -> one score per sample (higher = bacterial).

    Raises if any signature gene is missing: a gene is never silently dropped.
    `bacterial_weight` != 1 is only the pre-registered `ratio_one_term` variant.
    """
    missing = [g for g in sig.genes if g not in genes_x_samples.index]
    if missing:
        raise KeyError(f"{sig.name}: genes missing on this platform: {missing}")
    X = genes_x_samples.loc[list(sig.genes)].T.astype("float64")
    if X.isna().any().any():
        raise ValueError(f"{sig.name}: NaN in signature genes")
    b, v = X[list(sig.bacterial_up)], X[list(sig.viral_up)]
    if sig.rule == "sum":
        out = bacterial_weight * b.sum(axis=1) - v.sum(axis=1)
    elif sig.rule == "geomean":
        out = bacterial_weight * _geomean(b) - _geomean(v)
    else:
        raise ValueError(f"unknown rule {sig.rule!r}")
    return out.rename(sig.name)


def assert_per_sample(scorer, genes_x_samples: pd.DataFrame, sig: Signature) -> None:
    """Guard: a sample's score must not change when the OTHER samples change.

    Scores the full matrix, then each sample alone; any difference means the
    scorer reads other patients (a cohort-wide fit or normalisation).
    """
    full = scorer(genes_x_samples, sig)
    for col in genes_x_samples.columns[: min(10, genes_x_samples.shape[1])]:
        alone = scorer(genes_x_samples[[col]], sig)
        if not np.isclose(alone.iloc[0], full.loc[col], rtol=0, atol=1e-12):
            raise AssertionError(f"{sig.name}: score of {col} depends on other samples")


def metaintegrator_shift(genes_x_samples: pd.DataFrame) -> pd.DataFrame:
    """MetaIntegrator calculateScore: a dataset with min < 0 is shifted by |min| + 1."""
    lo = float(np.nanmin(genes_x_samples.to_numpy()))
    return genes_x_samples + (abs(lo) + 1.0) if lo < 0 else genes_x_samples


def metaintegrator_zscore(s: pd.Series, cohort: pd.Series) -> pd.Series:
    """MetaIntegrator calculateScore's scale(): z-score inside each cohort (SD with n-1).

    Reads the whole cohort, so this is NOT a per-patient frozen score; it is the
    authors' software reading, reported as its own variant.
    """
    if not s.index.equals(cohort.index):
        raise ValueError("score and cohort must share one index in one order")
    g = s.groupby(cohort)
    return (s - g.transform("mean")) / g.transform("std")


def _auc(y: pd.Series, s: pd.Series) -> float:
    """Plain AUC with bacterial = positive class. Never flipped."""
    return float(roc_auc_score(y.to_numpy(), s.to_numpy()))


def both_schemes(s: pd.Series, y: pd.Series, cohort: pd.Series, seed: int = 0,
                 n_splits: int = 5) -> dict:
    """Random-CV mean fold AUC (pooled) and per-cohort AUC for one fixed score."""
    if not (s.index.equals(y.index) and s.index.equals(cohort.index)):
        raise ValueError("score, y and cohort must share one index in one order")
    folds = [_auc(y.iloc[te], s.iloc[te]) for _, (_, te) in random_splits(y, n_splits, seed)]
    per = {c: _auc(y[cohort == c], s[cohort == c]) for c in pd.unique(cohort)}
    rcv, loco = float(np.mean(folds)), float(np.mean(list(per.values())))
    return {"random_cv": rcv, "loco": loco, "gap": rcv - loco, "per_cohort": per}
