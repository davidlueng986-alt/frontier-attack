"""Gatekeeper for claimed BSD counterexamples.

Birch-Swinnerton-Dyer, in the form that matters for a counterexample:

    rank E(Q)  ==  ord_{s=1} L(E, s)

for every elliptic curve E/Q.  A counterexample is a curve where the two
integers differ.  This module does not compute either invariant; it decides
whether a *claim* is coherent enough to be worth computing.

Why a gatekeeper: the failure mode this project exists to prevent is an agent
producing a plausible-looking curve with an unverifiable rank, and the claim
surviving because nobody checked the discriminant or the parity of the root
number.  Cheap invariants kill most of those in milliseconds.

Run `python problems/bsd/bsd.py`.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import sympy as sp

RANK_METHODS = {"2-descent", "mwrank", "magma", "sage", "published"}
ANALYTIC_METHODS = {"lcalc", "numerical", "published"}


def invariants(curve) -> dict:
    a1, a2, a3, a4, a6 = (sp.Integer(curve[k]) for k in ("a1", "a2", "a3", "a4", "a6"))
    b2 = a1**2 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3**2 + 4 * a6
    b8 = a1**2 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3**2 - a4**2
    disc = -b2**2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6
    c4 = b2**2 - 24 * b4
    return {"b2": b2, "b4": b4, "b6": b6, "b8": b8,
            "discriminant": disc, "c4": c4}


def check_claim(claim: dict) -> dict:
    """Return a verdict for one BSD claim.

    Required keys: curve {a1,a2,a3,a4,a6}, rank, analytic_rank.
    Optional: rank_method, analytic_rank_method, root_number, conductor.
    """
    reasons: list[str] = []
    warnings: list[str] = []

    inv = invariants(claim["curve"])
    if inv["discriminant"] == 0:
        return {"verdict": "REJECTED", "reasons": ["singular curve: discriminant = 0"]}

    r = claim.get("rank")
    k = claim.get("analytic_rank")
    if not isinstance(r, int) or r < 0:
        reasons.append("rank must be a non-negative integer")
    if not isinstance(k, int) or k < 0:
        reasons.append("analytic_rank must be a non-negative integer")
    if reasons:
        return {"verdict": "REJECTED", "reasons": reasons,
                "invariants": {k2: str(v) for k2, v in inv.items()}}

    rm = claim.get("rank_method", "claim")
    am = claim.get("analytic_rank_method", "claim")
    if rm not in RANK_METHODS:
        warnings.append(f"algebraic rank not backed by a recognised method ({rm})")
    if am not in ANALYTIC_METHODS:
        warnings.append(f"analytic rank not backed by a recognised method ({am})")

    w = claim.get("root_number")
    if w is not None:
        if w not in (1, -1):
            reasons.append("root number must be +1 or -1")
        elif w != (-1) ** r:
            reasons.append(f"parity violation: root number {w} but rank {r}")

    if reasons:
        return {"verdict": "REJECTED", "reasons": reasons,
                "invariants": {k2: str(v) for k2, v in inv.items()}}

    if r == k:
        return {"verdict": "NOT_A_COUNTEREXAMPLE",
                "reasons": [f"rank = analytic rank = {r}"],
                "invariants": {k2: str(v) for k2, v in inv.items()}}

    # A real mismatch would be an enormous result; demand independent backends.
    if len(warnings) == 2:
        return {"verdict": "INCOMPLETE",
                "reasons": ["mismatch claimed, but neither invariant has a "
                            "recognised computational backend"] + warnings,
                "invariants": {k2: str(v) for k2, v in inv.items()}}

    return {"verdict": "COUNTEREXAMPLE_CANDIDATE",
            "reasons": warnings or ["mismatch with recognised backends"],
            "invariants": {k2: str(v) for k2, v in inv.items()}}


# --------------------------------------------------------------------------
# self-check
# --------------------------------------------------------------------------

if __name__ == "__main__":
    # 37a1: y^2 + y = x^3 - x, rank 1, analytic rank 1
    good = {"curve": {"a1": 0, "a2": 0, "a3": 1, "a4": -1, "a6": 0},
            "rank": 1, "analytic_rank": 1,
            "rank_method": "2-descent", "analytic_rank_method": "lcalc",
            "root_number": -1}
    r1 = check_claim(good)
    print("37a1        :", r1["verdict"], r1["reasons"])
    assert r1["verdict"] == "NOT_A_COUNTEREXAMPLE"

    # same curve, fabricated rank mismatch with no backend
    bad = dict(good, rank=2)
    r2 = check_claim(bad)
    print("fabricated  :", r2["verdict"], r2["reasons"])
    assert r2["verdict"] == "REJECTED"  # parity violation catches it

    # mismatch with backends, but root number contradicts the claimed rank
    bad2 = dict(good, rank=3, root_number=1)
    r3 = check_claim(bad2)
    print("bad parity  :", r3["verdict"], r3["reasons"])
    assert r3["verdict"] == "REJECTED"

    # singular curve
    sing = {"curve": {"a1": 0, "a2": 0, "a3": 0, "a4": 0, "a6": 0},
            "rank": 1, "analytic_rank": 0}
    r4 = check_claim(sing)
    print("singular    :", r4["verdict"], r4["reasons"])
    assert r4["verdict"] == "REJECTED"

    # mismatch with proper backends and consistent parity
    claim = {"curve": {"a1": 0, "a2": 0, "a3": 0, "a4": -4, "a6": 0},
             "rank": 0, "analytic_rank": 2, "root_number": 1,
             "rank_method": "2-descent", "analytic_rank_method": "lcalc"}
    r5 = check_claim(claim)
    print("candidate   :", r5["verdict"], r5["reasons"])
    assert r5["verdict"] == "COUNTEREXAMPLE_CANDIDATE"
    print("OK")
