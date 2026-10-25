"""Guards for v1.1 frozen signatures. Each test builds the bug it guards against."""
from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src import frozen_signatures as fs

ROOT = Path(__file__).resolve().parent.parent
SIGS = {s.name: s for s in fs.load_signatures(ROOT / "config" / "signatures.json")}


def _cohort(n=40, seed=0, offset=0.0):
    """Synthetic genes x samples: bacterial patients have high FAM89A/HK3..., low IFI44L/IFI27..."""
    rng = np.random.default_rng(seed)
    y = np.array([1] * (n // 2) + [0] * (n - n // 2))
    genes = ["FAM89A", "IFI44L", "HK3", "TNIP1", "GPAA1", "CTSB", "IFI27", "JUP", "LAX1", "ACTB"]
    X = pd.DataFrame(rng.normal(8, 0.5, (len(genes), n)) + offset, index=genes,
                     columns=[f"s{i}" for i in range(n)])
    for g in ["FAM89A", "HK3", "TNIP1", "GPAA1", "CTSB"]:
        X.loc[g] += 1.5 * y
    for g in ["IFI44L", "IFI27", "JUP", "LAX1"]:
        X.loc[g] += 1.5 * (1 - y)
    return X, pd.Series(y, index=X.columns)


def test_real_config_pins_the_published_lists():
    assert set(SIGS["herberg2016_drs"].bacterial_up) == {"FAM89A"}
    assert set(SIGS["herberg2016_drs"].viral_up) == {"IFI44L"}
    assert set(SIGS["sweeney2016_bv7"].viral_up) == {"IFI27", "JUP", "LAX1"}
    assert set(SIGS["sweeney2016_bv7"].discovery) == {"GSE6269", "GSE40396", "GSE42026"}


def test_wrong_gene_list_in_config_is_refused(tmp_path):
    cfg = json.loads((ROOT / "config" / "signatures.json").read_text(encoding="utf-8"))
    cfg["signatures"][1]["genes_viral_up"] = ["IFI27", "JUP", "IFI44L"]  # LAX1 swapped out
    p = tmp_path / "sig.json"
    p.write_text(json.dumps(cfg), encoding="utf-8")
    with pytest.raises(ValueError, match="differs from the published"):
        fs.load_signatures(p)


def test_swapped_direction_is_refused():
    sig = SIGS["herberg2016_drs"]
    flipped = replace(sig, bacterial_up=sig.viral_up, viral_up=sig.bacterial_up)
    with pytest.raises(ValueError):
        fs.check_published(flipped)


def test_published_direction_scores_bacterial_high():
    X, y = _cohort()
    for sig in SIGS.values():
        assert fs._auc(y, fs.score(X, sig)) > 0.9


def test_reversed_score_is_reported_below_half_not_flipped():
    X, y = _cohort()
    s = -fs.score(X, SIGS["sweeney2016_bv7"])
    assert fs._auc(y, s) < 0.1


def test_missing_gene_raises_instead_of_being_dropped():
    X, _ = _cohort()
    with pytest.raises(KeyError, match="FAM89A"):
        fs.score(X.drop(index="FAM89A"), SIGS["herberg2016_drs"])


def test_per_sample_guard_fires_on_a_scorer_that_reads_other_patients():
    X, _ = _cohort()
    sig = SIGS["herberg2016_drs"]
    fs.assert_per_sample(fs.score, X, sig)  # the real scorer passes

    def leaky(m, s):  # centres on the cohort mean: a cohort-wide fit
        z = m.sub(m.mean(axis=1), axis=0)
        return fs.score(z, s)

    with pytest.raises(AssertionError, match="depends on other samples"):
        fs.assert_per_sample(leaky, X, sig)


def test_score_does_not_read_labels_or_other_rows():
    X, _ = _cohort()
    sig = SIGS["sweeney2016_bv7"]
    a = fs.score(X, sig)
    b = fs.score(X[X.columns[::-1]], sig)
    assert np.allclose(a.loc[X.columns], b.loc[X.columns])


def test_geomean_zero_is_zero_and_negative_refused():
    block = pd.DataFrame({"a": [0.0, 4.0], "b": [9.0, 1.0]})
    assert np.allclose(fs._geomean(block).to_numpy(), [0.0, 2.0])
    with pytest.raises(ValueError):
        fs._geomean(block - 1.0)


def test_pooling_offset_lowers_random_cv_but_not_loco():
    """The mechanism the PLAN predicts: a frozen score with a cohort offset."""
    parts = [_cohort(60, seed=i, offset=3.0 * i) for i in range(3)]
    sig = SIGS["herberg2016_drs"]
    s, y, c = [], [], []
    for i, (X, yy) in enumerate(parts):
        sc = fs.score(X, sig)
        # shift only the bacterial-up gene so the score itself carries the offset
        sc = sc + 2.0 * i
        sc.index = [f"c{i}:{k}" for k in sc.index]
        s.append(sc); y.extend(yy.to_list()); c.extend([f"c{i}"] * len(sc))
    s = pd.concat(s)
    r = fs.both_schemes(s, pd.Series(y, index=s.index), pd.Series(c, index=s.index))
    assert r["loco"] > 0.9
    assert r["gap"] < 0  # pooling hurts random CV; nothing was trained to exploit


# Transcribed from the supplements, independently of fs.PUBLISHED (LABBOOK PLAN 661ea85).
# Herberg eTable 3: log2 FC definite bacterial vs definite viral (+ = bacterial-up).
# Sweeney Table S1 / fig. S3: summary effect size viral vs bacterial (+ = viral-up).
HERBERG_ETABLE3_BACT_VS_VIRAL = {"FAM89A": 1.21, "IFI44L": -3.79}
SWEENEY_TABLE_S1_VIRAL_VS_BACT = {"IFI27": 2.299, "JUP": 1.008, "LAX1": 1.114, "HK3": -1.109,
                                  "TNIP1": -1.023, "GPAA1": -1.275, "CTSB": -1.080}


def _lists_from_supplement():
    h, s = HERBERG_ETABLE3_BACT_VS_VIRAL, SWEENEY_TABLE_S1_VIRAL_VS_BACT
    return {"herberg2016_drs": ({g for g, e in h.items() if e > 0}, {g for g, e in h.items() if e < 0}),
            "sweeney2016_bv7": ({g for g, e in s.items() if e < 0}, {g for g, e in s.items() if e > 0})}


def _supplement_mismatches(published):
    return [n for n, (b, v) in _lists_from_supplement().items()
            if (set(published[n][0]), set(published[n][1])) != (b, v)]


def test_published_lists_match_the_supplements():
    assert _supplement_mismatches(fs.PUBLISHED) == []
    for name, (b, v) in _lists_from_supplement().items():
        assert (set(SIGS[name].bacterial_up), set(SIGS[name].viral_up)) == (b, v)


def test_supplement_check_fires_on_a_wrong_list():
    wrong = dict(fs.PUBLISHED)
    wrong["sweeney2016_bv7"] = ({"HK3", "TNIP1", "GPAA1", "IFI44L"}, {"IFI27", "JUP", "LAX1"}, "geomean")
    assert _supplement_mismatches(wrong) == ["sweeney2016_bv7"]


def test_ratio_one_term_weights_the_bacterial_term_by_three_quarters():
    X, _ = _cohort()
    sig = SIGS["sweeney2016_bv7"]
    assert fs.ratio_weight(sig) == 3 / 4
    b = fs._geomean(X.loc[list(sig.bacterial_up)].T)
    v = fs._geomean(X.loc[list(sig.viral_up)].T)
    assert np.allclose(fs.score(X, sig, fs.ratio_weight(sig)), 0.75 * b - v)
    assert not np.allclose(fs.score(X, sig, 0.75), fs.score(X, sig))
    assert fs.ratio_weight(SIGS["herberg2016_drs"]) == 1.0  # 1 : 1 genes -> identical to main


def test_metaintegrator_zscores_per_cohort_and_keeps_within_cohort_auc():
    parts = [_cohort(40, seed=i, offset=2.0 * i) for i in range(2)]
    sig = SIGS["sweeney2016_bv7"]
    s = pd.concat([fs.score(X, sig).rename(lambda k, i=i: f"c{i}:{k}") for i, (X, _) in enumerate(parts)])
    y = pd.Series([v for _, yy in parts for v in yy], index=s.index)
    c = pd.Series([f"c{i}" for i, (X, _) in enumerate(parts) for _ in X.columns], index=s.index)
    z = fs.metaintegrator_zscore(s, c)
    assert np.allclose(z.groupby(c).mean(), 0) and np.allclose(z.groupby(c).std(), 1)
    for k in ("c0", "c1"):
        assert fs._auc(y[c == k], z[c == k]) == fs._auc(y[c == k], s[c == k])


def test_metaintegrator_shift_only_when_min_is_negative():
    X, _ = _cohort()
    assert fs.metaintegrator_shift(X) is X
    neg = X - 10.0
    shifted = fs.metaintegrator_shift(neg)
    assert np.isclose(shifted.to_numpy().min(), 1.0)
