"""Guards for the negative-gap diagnostics (LABBOOK 2026-10-19).

Each guard is shown to fire on the failure it exists for: a cohort-wide score
shift that fools a pooled AUC, an unbalanced set the balance check must reject,
a transform that reads other test rows, and labels out of row order.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from sklearn.base import clone
from sklearn.metrics import roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.gap_diagnostics import (balanced_subsample, cohort_rank_inductive,  # noqa: E402
                                 pair_auc, rank_against_reference)
from src.gap_experiments import (balance_control, both_schemes_with_scores,  # noqa: E402
                                 inductive_random_cv, score_mirror, split_half_heldout,
                                 summarize_balance)
from src.splits import cohort_splits, rank_normalize_within_cohort  # noqa: E402
from src.train import make_models  # noqa: E402

LR_ONLY = {"logistic_l2": make_models()["logistic_l2"]}


def _synthetic(shares, n=60, genes=60, seed=0):
    """Cohorts with ONE shared signal and no batch effect; only the bacterial share differs."""
    rng = np.random.default_rng(seed)
    blocks, ys, cs = [], [], []
    for k, share in enumerate(shares):
        y = (np.arange(n) < round(share * n)).astype(int)
        x = rng.normal(size=(n, genes))
        x[:, :10] += 0.8 * y[:, None]
        x[:, 10:20] -= 0.8 * y[:, None]
        blocks.append(x)
        ys.append(y)
        cs += [f"C{k}"] * n
    idx = [f"s{i}" for i in range(n * len(shares))]
    X = pd.DataFrame(np.vstack(blocks), index=idx, columns=[f"g{j}" for j in range(genes)])
    return X, pd.Series(np.concatenate(ys), index=idx), pd.Series(cs, index=idx)


# --- pair_auc --------------------------------------------------------------------------

def test_pair_auc_hand_example():
    y, s, g = [1, 0, 1, 0], [0.9, 0.1, 0.2, 0.8], ["A", "A", "B", "B"]
    pa = pair_auc(y, s, g)
    assert (pa.within, pa.n_within) == (0.5, 2)
    assert (pa.cross, pa.n_cross) == (1.0, 2)
    assert pa.pooled == 0.75 == roc_auc_score(y, s)


def test_pooled_is_sklearn_auc_and_the_pair_weighted_mix():
    rng = np.random.default_rng(1)
    y = rng.integers(0, 2, 120)
    s = rng.normal(size=120) + y
    s[:10] = s[10]                                   # ties count half
    g = rng.choice(list("ABC"), 120)
    pa = pair_auc(y, s, g)
    assert pa.pooled == pytest.approx(roc_auc_score(y, s), abs=1e-12)
    mix = (pa.within * pa.n_within + pa.cross * pa.n_cross) / (pa.n_within + pa.n_cross)
    assert pa.pooled == pytest.approx(mix, abs=1e-12)


def test_within_ignores_a_cohort_wide_shift_that_fools_pooled():
    """The suspected bug: one cohort's scores all shifted. Pooled must move (so the
    shift is a real perturbation); the same-cohort-pair AUC must not."""
    rng = np.random.default_rng(2)
    y = np.tile([1, 0], 40)
    s = rng.normal(size=80) + 1.5 * y
    g = np.repeat(["A", "B"], 40)
    before = pair_auc(y, s, g)
    after = pair_auc(y, s + np.where(g == "A", -3.0, 0.0), g)
    assert after.within == pytest.approx(before.within)
    assert after.pooled < before.pooled - 0.05


# --- balanced_subsample ----------------------------------------------------------------

def _is_balanced(y: pd.Series, cohort: pd.Series) -> bool:
    t = pd.crosstab(cohort, y)
    return bool((t[0] == t[1]).all())


def _unbalanced():
    cohort = pd.Series(["A"] * 30 + ["B"] * 20 + ["C"] * 10)
    y = pd.Series([1] * 6 + [0] * 24 + [1] * 16 + [0] * 4 + [1] * 5 + [0] * 5)
    return y, cohort


def test_balanced_subsample_balances_every_cohort():
    y, cohort = _unbalanced()
    assert not _is_balanced(y, cohort)               # break-it: the input fails the check
    idx = balanced_subsample(y, cohort, seed=0)
    assert _is_balanced(y.iloc[idx], cohort.iloc[idx])
    assert cohort.iloc[idx].value_counts().to_dict() == {"A": 12, "B": 8, "C": 10}
    assert np.array_equal(idx, balanced_subsample(y, cohort, seed=0))
    assert not np.array_equal(idx, balanced_subsample(y, cohort, seed=1))


def test_size_control_keeps_sizes_and_original_share():
    y, cohort = _unbalanced()
    idx = balanced_subsample(y, cohort, seed=0, keep_share=True)
    ys, cs = y.iloc[idx], cohort.iloc[idx]
    assert cs.value_counts().to_dict() == {"A": 12, "B": 8, "C": 10}
    assert ys[cs == "A"].sum() == 2                  # 12 x 6/30 = 2.4
    assert ys[cs == "B"].sum() == 6                  # 8 x 16/20 = 6.4
    assert not _is_balanced(ys, cs)


# --- one-row-at-a-time ranking ---------------------------------------------------------

def test_rank_against_reference_matches_pandas_rank_after_joining():
    rng = np.random.default_rng(3)
    ref = rng.integers(0, 5, size=(9, 6)).astype(float)    # many ties
    row = rng.integers(0, 5, size=(1, 6)).astype(float)
    joined = pd.DataFrame(np.vstack([ref, row])).rank(pct=True).to_numpy()[-1]
    assert np.allclose(rank_against_reference(ref, row)[0], joined)


def _reads_other_new_rows(transform) -> bool:
    rng = np.random.default_rng(4)
    ref, new = rng.normal(size=(15, 8)), rng.normal(size=(6, 8))
    before = transform(ref, new)[2]
    new2 = new.copy()
    new2[[0, 1, 3, 4, 5]] = rng.normal(5.0, 3.0, size=(5, 8))
    return not np.allclose(before, transform(ref, new2)[2])


def test_rank_against_reference_reads_no_other_new_row():
    assert not _reads_other_new_rows(rank_against_reference)


def test_guard_fires_on_ranking_new_rows_together():
    """Break-it: ranking new rows with the reference AND each other is transduction."""
    def transductive(ref, new):
        return pd.DataFrame(np.vstack([ref, new])).rank(pct=True).to_numpy()[len(ref):]
    assert _reads_other_new_rows(transductive)


def test_cohort_rank_inductive_fits_on_training_rows_only():
    rng = np.random.default_rng(5)
    X = pd.DataFrame(rng.normal(size=(20, 5)))
    cohort = pd.Series(["A"] * 12 + ["B"] * 8)
    train = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 15, 16])
    test = np.array([9, 10, 11, 17, 18, 19])
    X_tr, X_te = cohort_rank_inductive(X, cohort, train, test)
    expected = rank_normalize_within_cohort(X.iloc[train].reset_index(drop=True),
                                            cohort.iloc[train].reset_index(drop=True))
    assert np.allclose(X_tr, expected.to_numpy())
    X2 = X.copy()
    X2.iloc[9] = 100.0                                # change one TEST row only
    X_tr2, X_te2 = cohort_rank_inductive(X2, cohort, train, test)
    assert np.allclose(X_tr, X_tr2) and np.allclose(X_te[1:], X_te2[1:])
    assert not np.allclose(X_te[0], X_te2[0])


def test_cohort_rank_inductive_refuses_a_test_cohort_without_training_rows():
    X, cohort = pd.DataFrame(np.zeros((4, 2))), pd.Series(["A", "A", "B", "B"])
    with pytest.raises(ValueError, match="no training rows"):
        cohort_rank_inductive(X, cohort, np.array([0, 1]), np.array([2, 3]))


# --- experiment plumbing ---------------------------------------------------------------

def test_misaligned_labels_are_refused():
    """Break-it: labels in a different row order than X must stop the run."""
    X, y, cohort = _synthetic([0.5, 0.5])
    both_schemes_with_scores(X, y, cohort, "raw", LR_ONLY)
    with pytest.raises(AssertionError, match="not aligned"):
        both_schemes_with_scores(X, y.iloc[::-1], cohort, "raw", LR_ONLY)


def test_positive_control_unequal_balance_plus_cohort_ranking_hurts_cross_pairs():
    """The mechanism on synthetic data with no batch effect and one shared signal:
    unequal bacterial shares + cohort-wise ranking drag the pooled random-CV AUC
    below the same-cohort-pair AUC; equal shares do not."""
    def within_minus_pooled(shares):
        X, y, cohort = _synthetic(shares)
        folds, _ = both_schemes_with_scores(X, y, cohort, "cohort_rank", LR_ONLY)
        r = folds[folds["scheme"] == "random_cv"]
        return r["auc_within"].mean() - r["auc_pooled"].mean()
    assert within_minus_pooled([0.2, 0.2, 0.8, 0.8]) > 0.05
    assert abs(within_minus_pooled([0.5, 0.5, 0.5, 0.5])) < 0.03


def test_split_half_whole_cohort_auc_is_the_heldout_auc():
    X, y, cohort = _synthetic([0.3, 0.5, 0.7])
    sh = split_half_heldout(X, y, cohort, LR_ONLY, n_halves=3)
    folds, _ = both_schemes_with_scores(X, y, cohort, "cohort_rank", LR_ONLY)
    loco = folds[folds["scheme"] == "leave_one_cohort_out"].set_index("fold")["auc_pooled"]
    full = sh.groupby("cohort")["auc_full_cohort"].first()
    assert np.allclose(full.sort_index(), loco.sort_index())
    assert sh["auc_b_half_a"].notna().all()


def test_split_half_whole_column_is_the_whole_cohort_auc_on_half_b():
    """Guard for E4's key column (outside check 19.10.): auc_b_whole must be the AUC on
    half B of the scores from WHOLE-cohort ranks. If it silently became the half-A
    AUC, P5's difference would be 0 by construction and 'supported' for free."""
    X, y, cohort = _synthetic([0.3, 0.5, 0.7], n=40)
    n_halves = 6
    sh = split_half_heldout(X, y, cohort, LR_ONLY, n_halves=n_halves)
    X_cr = rank_normalize_within_cohort(X, cohort).to_numpy(dtype=float)
    yv = y.to_numpy().astype(int)
    for held_out, (tr, te) in cohort_splits(cohort):
        est = clone(LR_ONLY["logistic_l2"]).fit(X_cr[tr], yv[tr])
        p_whole = est.predict_proba(X_cr[te])[:, 1]
        for k in range(n_halves):
            b = np.random.default_rng(k).permutation(len(te))[len(te) // 2:]
            row = sh[(sh["cohort"] == held_out) & (sh["split"] == k)]
            assert row["auc_b_whole"].item() == pytest.approx(
                roc_auc_score(yv[te][b], p_whole[b]), abs=1e-12)
    # The two rankings must give different AUCs here, or this data could not tell
    # the right column from the half-A one.
    assert (sh["auc_b_whole"] - sh["auc_b_half_a"]).abs().max() > 1e-3


def test_inductive_random_cv_scores_every_fold():
    X, y, cohort = _synthetic([0.3, 0.7])
    df = inductive_random_cv(X, y, cohort, LR_ONLY)
    assert len(df) == 5 and df["auc_pooled"].between(0, 1).all()


def test_mirror_whole_cohort_ranks_fix_each_cohorts_total_score():
    """With whole-cohort ranks a linear model's train + test scores of one cohort add
    up to the same total in every equal-sized cohort (the mirror effect). With
    training-only ranks (E3) that tie must be gone, or the check proves nothing."""
    X, y, cohort = _synthetic([0.2, 0.5, 0.8])                # 3 cohorts of 60
    df = score_mirror(X, y, cohort, LR_ONLY)
    df["total"] = df["n"] * df["mean_logit"]
    per = df.groupby(["ranking", "fold", "cohort"])["total"].sum() / 60
    spread = per.groupby(level=["ranking", "fold"]).agg(lambda s: s.max() - s.min())
    assert (spread.loc["whole_cohort"] < 1e-6).all()
    assert (spread.loc["train_only"] > 1e-3).all()


def test_balance_summary_gap_is_random_cv_minus_heldout():
    X, y, cohort = _synthetic([0.3, 0.7, 0.5])
    s = summarize_balance(balance_control(X, y, cohort, "balanced", range(2), LR_ONLY))
    assert set(s["transform"]) == {"raw", "cohort_rank", "within_sample"}
    assert np.allclose(s["gap_mean"], s["random_cv"] - s["loco"])
    assert (s["subsets"] == 2).all()
