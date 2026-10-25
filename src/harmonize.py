"""Put cohorts on one shared feature axis: probe ids -> gene symbols -> the
symbols every cohort has.

Each cohort uses a different microarray platform, so raw probe ids are not
comparable. Collapsing many probes per gene by picking the highest-mean probe is
the common convention and is applied identically to every cohort.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def to_symbols(expr: pd.DataFrame, probe_to_symbol: pd.Series,
               collapse: str = "mean") -> pd.DataFrame:
    """probes x samples -> genes x samples, keeping the highest-mean probe per gene.

    collapse="first" keeps the lowest probe id instead and reads no expression value:
    the control for review finding 2026-10-21 (the highest-mean choice reads every
    patient of the series, test patients included).
    """
    shared = expr.index.intersection(probe_to_symbol.index)
    if len(shared) == 0:
        return pd.DataFrame()
    sub = expr.loc[shared]
    symbols = probe_to_symbol.loc[shared].astype(str).values
    if collapse == "mean":
        order = np.argsort(-sub.mean(axis=1, skipna=True).values, kind="stable")
    elif collapse == "first":
        order = np.argsort(sub.index.astype(str).values, kind="stable")
    else:
        raise ValueError(f"unknown collapse rule {collapse!r}")
    sub = sub.iloc[order]
    symbols = symbols[order]
    keep = ~pd.Index(symbols).duplicated(keep="first")
    out = sub.loc[keep]
    out.index = pd.Index(symbols[keep], name="symbol")
    return out


def maybe_log2(expr: pd.DataFrame) -> tuple[pd.DataFrame, bool]:
    """Microarray series are archived on linear or log scale. Put them all on log2."""
    finite = expr.to_numpy(dtype="float64", na_value=np.nan)
    hi = np.nanmax(finite) if np.isfinite(finite).any() else 0.0
    if hi > 100:
        return np.log2(expr.clip(lower=0) + 1.0), True
    return expr, False


def build_matrix(
    per_cohort: dict[str, pd.DataFrame],
    labels: dict[str, dict[str, str]],
    min_genes: int = 500,
) -> tuple[pd.DataFrame, pd.Series, pd.Series, list[str]]:
    """Return X (samples x shared genes), y (0 viral / 1 bacterial), cohort, gene list."""
    notes: list[str] = []
    gene_sets = [set(df.index) for df in per_cohort.values() if not df.empty]
    if not gene_sets:
        raise ValueError("no cohort produced a gene-level matrix")
    shared = sorted(set.intersection(*gene_sets))
    notes.append(f"shared gene symbols across {len(per_cohort)} cohorts: {len(shared)}")
    if len(shared) < min_genes:
        notes.append(f"WARNING: only {len(shared)} shared genes (< {min_genes})")

    blocks, ys, groups = [], [], []
    for accession, df in per_cohort.items():
        cohort_labels = labels[accession]
        samples = [s for s in df.columns if s in cohort_labels]
        if not samples:
            notes.append(f"{accession}: no labelled samples, skipped")
            continue
        block = df.loc[shared, samples].T
        block.index = [f"{accession}:{s}" for s in samples]
        blocks.append(block)
        ys.extend(1 if cohort_labels[s] == "bacterial" else 0 for s in samples)
        groups.extend([accession] * len(samples))

    X = pd.concat(blocks, axis=0)
    notes.append(f"imputed cells (cohort gene mean, else 0.0): {int(X.isna().sum().sum())}")
    # Genes absent from a given array land as NaN; impute per cohort to that
    # cohort's own gene mean so no information crosses the cohort boundary.
    X = X.groupby(pd.Series(groups, index=X.index), observed=True).transform(
        lambda col: col.fillna(col.mean())
    )
    X = X.fillna(0.0)
    y = pd.Series(ys, index=X.index, name="bacterial")
    cohort = pd.Series(groups, index=X.index, name="cohort")
    return X, y, cohort, notes
