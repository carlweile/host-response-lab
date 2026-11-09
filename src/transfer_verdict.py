"""The pre-registered rules of the 2026-11-05 PLAN (the transfer-matrix run), as pure functions.

Thresholds are copied from the PLAN entry + the pre-run clarifications entry
(commit fb2b99b), both committed before any run. Do not tune after seeing data.
"""
from __future__ import annotations

MAIN_TRACK = 0.03       # same-family minus cross-family, main arm
CONTROL_TRACK = 0.02    # same, size-matched control
NOT_TRACK = 0.01        # at or below this = "does not"
SCALE_HALF = 0.5
BAL_MATTERS = 0.02
BAL_NOT = 0.01
EPS = 1e-12             # edges are inclusive


def tracks(main: float, control: float) -> str:
    """R-platform / R-age rule."""
    if main >= MAIN_TRACK - EPS and control >= CONTROL_TRACK - EPS:
        return "tracks"
    if main <= NOT_TRACK + EPS:
        return "does not"
    return "undecided"


def r_platform(main: float, control: float) -> str:
    return tracks(main, control)


def r_age(main: float, control: float, testable: bool) -> str:
    return tracks(main, control) if testable else "not testable"


def r_both(platform: str, age: str) -> str:
    return "confounded, no cause named" if platform == age == "tracks" else "not fired"


def r_scale(platform: str, raw_diff: float, within_sample_diff: float) -> str:
    if platform != "tracks":
        return "not evaluated"
    return ("mostly a scale offset" if within_sample_diff <= SCALE_HALF * raw_diff + EPS
            else "not mostly a scale offset")


def r_balance(deltas: dict) -> str:
    """deltas: model -> mean LOCO AUC without weighting minus with weighting."""
    if any(abs(d) >= BAL_MATTERS - EPS for d in deltas.values()):
        return "class weighting matters here"
    if all(abs(d) < BAL_NOT for d in deltas.values()):
        return "not where the drop lives"
    return "undecided"
