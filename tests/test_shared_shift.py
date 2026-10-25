"""Guards for the shared-shift test (LABBOOK PLAN 2fd8865, gate 3).

Each guard was checked to go red on a deliberately broken copy of the code it
protects (mutations listed in the LABBOOK RESULT entry).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import shared_shift as sh  # noqa: E402
from src import shared_shift_verdict as hv  # noqa: E402
from tests.test_gap_diagnostics import _synthetic  # noqa: E402


def test_balanced_weights_give_each_class_half_the_mass():
    w = sh.balanced_weights(np.array([1, 0, 0, 0]))
    assert w[0] == pytest.approx(0.5) and w[1:].sum() == pytest.approx(0.5)


def test_weighted_rank_uses_class_mass_not_head_count():
    ref = np.array([[1.0], [2.0], [3.0], [4.0]])     # 1 bacterial (value 1), 3 viral
    w = sh.balanced_weights(np.array([1, 0, 0, 0]))
    # above the one bacterial, below all viral: half the mass lies below it
    assert sh.rank_weighted(ref, w, np.array([[1.5]]))[0, 0] == pytest.approx(0.5)


def test_balanced_arm_reads_no_test_row():
    X, y, c = _synthetic([0.3, 0.7], n=30, genes=10)
    tr, te = np.arange(0, 60, 2), np.arange(1, 60, 2)
    A_tr, A_te = sh.rank_balanced_inductive(X, y, c, tr, te)
    X2 = X.copy()
    X2.iloc[te[0]] += 100.0
    B_tr, B_te = sh.rank_balanced_inductive(X2, y, c, tr, te)
    np.testing.assert_array_equal(A_tr, B_tr)
    np.testing.assert_array_equal(A_te[1:], B_te[1:])


def test_recentring_uses_training_offsets_only():
    """Training scores carry no cohort shift, test scores do: recentring must leave
    the test scores alone (it would remove the shift if it read the test side)."""
    y_tr, c_tr = np.tile([0, 1], 10), np.repeat(["A", "B"], 10)
    s_tr = np.where(y_tr == 1, 1.0, -1.0)
    c_te = np.repeat(["A", "B"], 4)
    s_te = np.array([0.0, 1, 0, 1, 5, 6, 5, 6])
    np.testing.assert_array_equal(sh.recentred(s_tr, y_tr, c_tr, s_te, c_te), s_te)


def test_recentring_removes_a_shift_seen_in_training():
    y_tr, c_tr = np.tile([0, 1], 10), np.repeat(["A", "B"], 10)
    s_tr = np.where(y_tr == 1, 1.0, -1.0) + 2.0 * (c_tr == "B")
    out = sh.recentred(s_tr, y_tr, c_tr, np.array([3.0]), np.array(["B"]))
    assert out[0] == pytest.approx(2.0)          # B's training offset is +1 (half of 2)


BASE = dict(D_W=0.3, D_I=0.2, D_I_bal=0.05, D_raw=0.0, cross_W=0.70, cross_W_recentre=0.85,
            cross_raw=0.90, X2=0.10, X5=0.10, X10=0.11, R=0.75, signs_agree=True,
            slope_I=1.0, slope_I_bal=0.3, max_n_iter=50.0)


def test_hand_table_passes():
    assert hv.judge(BASE) == {"H1": "PASS", "H2": "PASS", "H3": "PASS"}
    assert hv.overall({"all": hv.judge(BASE), "no40396": hv.judge(BASE)}) == "SUPPORTED"


@pytest.mark.parametrize("change, h, word", [
    (dict(R=-0.1), "H1", "KILL"),
    (dict(R=0.3), "H1", "PARTIAL"),
    (dict(signs_agree=False), "H1", "PARTIAL"),
    (dict(cross_raw=0.6), "H1", "VOID"),
    (dict(D_I_bal=0.19), "H2", "KILL"),          # >= 0.9 D_I
    (dict(D_I_bal=0.15), "H2", "PARTIAL"),
    (dict(slope_I_bal=0.8), "H2", "PARTIAL"),
    (dict(D_I=0.0), "H2", "VOID"),
    (dict(X10=0.16), "H3", "KILL"),              # >= 1.5 X2
    (dict(X10=0.14), "H3", "PARTIAL"),
    (dict(X5=-0.01), "H3", "PARTIAL"),
    (dict(max_n_iter=3000.0), "H2", "VOID"),
])
def test_each_rule_fires(change, h, word):
    v = hv.judge(dict(BASE, **change))
    assert v[h].startswith(word)
    assert hv.overall({"all": hv.judge(BASE), "no40396": v}).startswith("NOT SUPPORTED")


def test_quantities_on_a_tiny_frame():
    rows = []
    for arm, (d_in, d_cr) in {"W": (1.0, 0.6), "W_recentre": (1.0, 0.8), "I": (1.0, 0.8),
                              "I_bal": (1.0, 0.95), "raw": (1.0, 1.0)}.items():
        rows.append(dict(set="s", part="k5", model="m", arm=arm, k=5, C=0.05, same=d_in,
                         cross=d_cr, D=d_in - d_cr, n_iter=10.0))
    summary = pd.DataFrame(rows)
    off = pd.DataFrame([dict(set="s", arm=a, cohort=g, share=sh_, pooled_share=0.5,
                             o_test=o, o_train=o)
                        for a in ("W", "I", "I_bal") for g, sh_, o in
                        (("A", 0.2, 1.0 if a != "I_bal" else 0.2),
                         ("B", 0.8, -1.0 if a != "I_bal" else -0.2))])
    q = hv.quantities(summary, off)["s"]
    assert q["R"] == pytest.approx(0.5)
    assert q["signs_agree"]
    assert q["slope_I_bal"] == pytest.approx(0.2 * q["slope_I"])
