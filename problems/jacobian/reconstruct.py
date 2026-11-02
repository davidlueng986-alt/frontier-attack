"""Reconstruction audit for the GGHV normalisation (arXiv:2204.14178, Prop. 4.3).

The question
------------
Prop. 4.3 of Guccione--Guccione--Horruitiner--Valqui reduces a hypothetical
planar Jacobian counterexample in the paper's case (8,28) -- the remaining
(72,108) degree case of Theorem 2.1 -- to a pair

    P, Q in L^(1) = K[x, x^-1, y],   [P,Q] = x^2,

with one of two prescribed Newton polygons.  The reduction is one-directional:
it starts from a polynomial pair and applies a chain whose factors are
automorphisms of K[x,y] (the swap) or of L^(1) (the shears and phi).  The
composite is a morphism K[x,y] -> L^(1), not an automorphism of either ring.
This script asks whether the chain can be reversed to give a polynomial pair
with constant Jacobian determinant.

The published chain (Prop. 4.3 proof; PDF p.10, src.tex lines 1000-1311)
-------------------------------------------------------------------------
  T1  swap:              (x,y) -> (y,x)                    [line 1012]
  T2  shear:             y -> y + lambda_1 x^-2            [line 1074]
  T3  shear:             y -> y + lambda x^-3               [lines 1079,1082]
  T4  shear:             y -> y + alpha x^-4                [line 1132]
  T5  monomial/Laurent:  (x,y) -> (x^-1, x^4 y)             [line 1229]

The composite is A = T5 o T4 o T3 o T2 o T1 (T1 is applied first):

    A(x,y) = (x^4 y + s(x^-1), x^-1),
    s(x)   = c2 x^-2 + c3 x^-3 + c4 x^-4.

Its rational inverse is

    A^-1(u,v) = (v^-1, v^4 (u - s(v))).

So, for a normalized pair (P,Q), the only candidate reconstruction is

    P0(x,y) = P(1/y, y^4 (x - s(y))),
    Q0(x,y) = Q(1/y, y^4 (x - s(y))).

Before clearing denominators, det J(P0,Q0) = 1 exactly: the chain rule gives
det J(A) = x^2 and [P,Q] = x^2, so the x^2 cancels.  The obstruction is
polynomiality: A^-1 has v^-1 in its first component.  In particular the inverse
swap (the last inverse step) turns the x^-1 introduced by T5 into y^-1, which
is not in L^(1) and is not a polynomial.

J3 test
-------
P = a0 + a1 x,  Q = b0 + b1 x^2 y,  a1 b1 = 1,  so [P,Q] = x^2.

The inverse gives

    P0 = a0 + a1/y,                    (independent of s)
    Q0 = b0 + b1 y^2 (x - s(y)).

P0 has a pole at y=0 for every choice of c2,c3,c4 because a1 != 0.  Clearing
denominators (multiply P by y, Q by y^2) gives

    P' = a0 y + a1,
    Q' = b0 y^2 - b1 c2 y^2 - b1 c3 y - b1 c4 + b1 x y^4,
    det J(P',Q') = -a0 b1 y^4,

which is not a nonzero constant (if a0 = 0 it is 0; otherwise it is a nonzero
multiple of y^4).  So J3 is NOT reconstructible, and no choice of the shear
constants repairs it.

Verdict
-------
PARTIAL.  The inverse exists as a rational map and restores det J = 1 before
clearing, but it is not a map L^(1) -> K[x,y].  The paper supplies neither the
inverse of the composite nor a description of its image.  The shear constants
c2,c3,c4 are read off from the original polynomial P's leading forms and are
not determined by the normalized pair.  A config-2 solution therefore needs,
in addition, the polynomiality (image) conditions - the vanishing of all
negative powers of y in P0 and Q0 - which Prop. 4.3 does not state.

Run:
    python problems/jacobian/reconstruct.py
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

import sympy as sp

PAPER = (
    "arXiv:2204.14178v1, Prop. 4.3 "
    "(PDF p.10; LaTeX src.tex lines 1000-1311)"
)


@dataclass
class Stage:
    """One step of the inverse normalisation."""

    index: int
    name: str
    P: sp.Expr
    Q: sp.Expr


# ---------------------------------------------------------------------------
# basic Laurent-polynomial utilities
# ---------------------------------------------------------------------------

def neg_powers(expr: sp.Expr, y: sp.Symbol) -> Dict[int, sp.Expr]:
    """Return {k: coefficient} for the negative powers y^k (k < 0) in expr.

    The expressions produced here are rational in y with denominator a power
    of y.  Multiplying by y**N for N larger than the denominator degree turns
    the check into a finite polynomial coefficient extraction.
    """
    expr = sp.cancel(sp.together(expr))
    _num, den = sp.fraction(expr)
    dp = sp.Poly(den, y)
    d = max((m[0] for m in dp.monoms()), default=0)
    N = d + 1
    scaled = sp.expand(expr * y**N)
    p = sp.Poly(scaled, y)
    out: Dict[int, sp.Expr] = {}
    for (e,), coeff in p.terms():
        if e < N:
            out[e - N] = sp.factor(coeff)
    return out


def is_in_L1(expr: sp.Expr, y: sp.Symbol) -> bool:
    """True iff expr is a polynomial in y (x may have negative powers)."""
    try:
        sp.Poly(sp.expand(expr), y)
        return True
    except Exception:
        return False


def is_polynomial(expr: sp.Expr, x: sp.Symbol, y: sp.Symbol) -> bool:
    """True iff expr is a polynomial in x and y."""
    try:
        sp.Poly(sp.expand(expr), x, y)
        return True
    except Exception:
        return False


def det_jac(P: sp.Expr, Q: sp.Expr, x: sp.Symbol, y: sp.Symbol) -> sp.Expr:
    """Exact Jacobian determinant [P,Q] = P_x Q_y - P_y Q_x."""
    return sp.expand(sp.Matrix([P, Q]).jacobian([x, y]).det())


# ---------------------------------------------------------------------------
# the inverse chain
# ---------------------------------------------------------------------------

def forward_stages(
    P: sp.Expr,
    Q: sp.Expr,
    x: sp.Symbol,
    y: sp.Symbol,
    shears: Sequence[sp.Expr] = (0, 0, 0),
) -> Tuple[sp.Expr, sp.Expr]:
    """The published forward chain, used only for the round-trip self-test.

    Forward order: swap, then the shears, then phi.  In sympy, substitution
    composes on the left, so this order is S first and phi last.
    """
    c2, c3, c4 = shears
    s = c2 / x**2 + c3 / x**3 + c4 / x**4

    P1 = P.subs({x: y, y: x}, simultaneous=True)
    Q1 = Q.subs({x: y, y: x}, simultaneous=True)
    P2 = P1.subs(y, y + s)
    Q2 = Q1.subs(y, y + s)
    P3 = P2.subs({x: 1 / x, y: x**4 * y}, simultaneous=True)
    Q3 = Q2.subs({x: 1 / x, y: x**4 * y}, simultaneous=True)
    return P3, Q3


def inverse_stages(
    P: sp.Expr,
    Q: sp.Expr,
    x: sp.Symbol,
    y: sp.Symbol,
    shears: Sequence[sp.Expr] = (0, 0, 0),
) -> List[Stage]:
    """Undo the published chain, one transformation at a time.

    ``shears`` are (c2, c3, c4), the coefficients of x^-2, x^-3, x^-4 in the
    forward shear y -> y + c2 x^-2 + c3 x^-3 + c4 x^-4.

    The inverse order is the reverse of the forward order:
      1. undo T5 (phi is an involution),
      2. undo the shears (they commute and add),
      3. undo T1 (the swap).
    """
    c2, c3, c4 = shears
    s = c2 / x**2 + c3 / x**3 + c4 / x**4

    # Step 1: undo phi(x) = x^-1, phi(y) = x^4 y.  phi is its own inverse.
    P1 = P.subs({x: 1 / x, y: x**4 * y}, simultaneous=True)
    Q1 = Q.subs({x: 1 / x, y: x**4 * y}, simultaneous=True)

    # Step 2: undo the shears y -> y + s(x).
    P2 = P1.subs(y, y - s)
    Q2 = Q1.subs(y, y - s)

    # Step 3: undo the swap x <-> y.
    P3 = P2.subs({x: y, y: x}, simultaneous=True)
    Q3 = Q2.subs({x: y, y: x}, simultaneous=True)

    return [
        Stage(1, "undo phi: (x,y) -> (x^-1, x^4 y)", P1, Q1),
        Stage(2, "undo shears: y -> y - (c2/x^2 + c3/x^3 + c4/x^4)", P2, Q2),
        Stage(3, "undo swap: (x,y) -> (y,x)", P3, Q3),
    ]


def reconstruct(
    P: sp.Expr,
    Q: sp.Expr,
    x: sp.Symbol,
    y: sp.Symbol,
    shears: Sequence[sp.Expr] = (0, 0, 0),
) -> dict:
    """Apply the inverse chain and clear denominators.

    Returns a dict with the stages, the negative-power coefficients, the
    cleared pair, and both determinants (before and after clearing).

    The clearing powers are the largest negative y-powers actually present
    in P0 and Q0, so the cleared pair uses the minimal common denominator:
    the inverse of T5 puts x -> 1/y, and every negative power is a power of y.
    """
    stages = inverse_stages(P, Q, x, y, shears)
    P3, Q3 = stages[-1].P, stages[-1].Q

    negP = neg_powers(P3, y)
    negQ = neg_powers(Q3, y)

    DP = max([0] + [-k for k in negP if k < 0])
    DQ = max([0] + [-k for k in negQ if k < 0])

    Pc = sp.expand(y**DP * P3)
    Qc = sp.expand(y**DQ * Q3)

    det_before = sp.simplify(det_jac(P3, Q3, x, y))
    det_after = sp.expand(det_jac(Pc, Qc, x, y))

    return {
        "stages": stages,
        "negP": negP,
        "negQ": negQ,
        "DP": DP,
        "DQ": DQ,
        "P0": P3,
        "Q0": Q3,
        "Pc": Pc,
        "Qc": Qc,
        "det_before": det_before,
        "det_after": det_after,
        "P0_polynomial": is_polynomial(P3, x, y),
        "Q0_polynomial": is_polynomial(Q3, x, y),
        "Pc_polynomial": is_polynomial(Pc, x, y),
        "Qc_polynomial": is_polynomial(Qc, x, y),
        "det_constant": bool(det_after.is_number and det_after != 0),
    }


def self_test() -> None:
    """Check the forward/inverse pair on a known polynomial map.

    P0 = x, Q0 = y has det J = 1.  The forward chain must multiply the
    determinant by x^2 (the chain rule in Prop. 4.3), and the inverse chain
    must recover (x,y).  This guards against a sign or exponent error in the
    implementation; J3's failure below is then a fact about the normalisation,
    not about the code.
    """
    x, y = sp.symbols("x y")
    c2, c3, c4 = sp.symbols("c2 c3 c4")
    shears = (c2, c3, c4)

    Pn, Qn = forward_stages(x, y, x, y, shears)
    detn = sp.simplify(det_jac(Pn, Qn, x, y))
    assert sp.simplify(detn - x**2) == 0, f"forward det is {detn}, expected x^2"

    rec = inverse_stages(Pn, Qn, x, y, shears)
    assert sp.simplify(rec[-1].P - x) == 0, rec[-1].P
    assert sp.simplify(rec[-1].Q - y) == 0, rec[-1].Q

    # Positive control: a non-trivial polynomial pair in the image must be
    # reconstructed exactly, with constant determinant 1.
    P0, Q0 = x + y**2, y
    Pn, Qn = forward_stages(P0, Q0, x, y, shears)
    res = reconstruct(Pn, Qn, x, y, shears)
    assert sp.simplify(res["Pc"] - P0) == 0, res["Pc"]
    assert sp.simplify(res["Qc"] - Q0) == 0, res["Qc"]
    assert sp.simplify(res["det_after"] - 1) == 0, res["det_after"]

    print("self-test: forward det = x^2; inverse(forward(x,y)) = (x,y);")
    print("           positive control (x+y^2, y) reconstructs with det = 1  [ok]")


# ---------------------------------------------------------------------------
# the J3 family
# ---------------------------------------------------------------------------

def j3_audit() -> dict:
    """P = a0 + a1 x, Q = b0 + b1 x^2 y, a1 b1 = 1, [P,Q] = x^2."""
    x, y = sp.symbols("x y")
    a0, a1, b0, b1, c2, c3, c4 = sp.symbols("a0 a1 b0 b1 c2 c3 c4")

    P = a0 + a1 * x
    Q = b0 + b1 * x**2 * y

    result = reconstruct(P, Q, x, y, (c2, c3, c4))
    result.update({"x": x, "y": y, "P": P, "Q": Q,
                   "a0": a0, "a1": a1, "b0": b0, "b1": b1,
                   "c2": c2, "c3": c3, "c4": c4})
    return result


def representative_audit() -> dict:
    """A concrete member of J3: P = x, Q = x^2 y (a0=b0=0, a1=b1=1)."""
    x, y = sp.symbols("x y")
    return reconstruct(x, x**2 * y, x, y, (0, 0, 0))


# ---------------------------------------------------------------------------
# reporting
# ---------------------------------------------------------------------------

def print_j3_report(r: dict) -> None:
    x, y = r["x"], r["y"]
    a0, a1, b1 = r["a0"], r["a1"], r["b1"]

    print("J3 family:  P = a0 + a1*x,  Q = b0 + b1*x^2*y,  a1*b1 = 1")
    print("            [P,Q] = x^2")
    print()

    for st in r["stages"]:
        print(f"step {st.index}: {st.name}")
        print(f"  P = {sp.expand(st.P)}")
        print(f"  Q = {sp.expand(st.Q)}")
        print(f"  polynomial in x,y: {is_polynomial(st.P, x, y)} / "
              f"{is_polynomial(st.Q, x, y)}")
        print(f"  lies in L^(1)=K[x,x^-1,y]: {is_in_L1(st.P, y)} / "
              f"{is_in_L1(st.Q, y)}")
        print()

    print("negative powers after the full inverse (step 3):")
    print(f"  P: {r['negP']}")
    print(f"  Q: {r['negQ']}")
    print()

    print("exact inverse, before clearing:")
    print(f"  P0 = {r['P0']}")
    print(f"  Q0 = {r['Q0']}")
    print(f"  det J(P0,Q0) = {r['det_before']}   (= 1 under a1*b1 = 1)")
    print()

    print("cleared pair (multiply P by y^%d, Q by y^%d):" % (r["DP"], r["DQ"]))
    print(f"  P' = {r['Pc']}")
    print(f"  Q' = {r['Qc']}")
    print(f"  det J(P',Q') = {sp.factor(r['det_after'])}")
    print(f"  constant? {r['det_constant']}")
    print()

    # The exact statement for the family.
    det_rel = sp.factor(r["det_after"].subs(b1, 1 / a1))
    print("substituting b1 = 1/a1:")
    print(f"  det J(P',Q') = {det_rel}")
    print("  (a0 = 0 gives det = 0; a0 != 0 gives a nonzero multiple of y^4)")
    print()

    print("NOT RECONSTRUCTIBLE: step 3, the inverse swap turns the x^-1 from")
    print("  step 1 into y^-1.  P0 = a0 + a1/y is not a polynomial for any")
    print("  choice of c2,c3,c4, because a1 != 0.  Clearing denominators gives")
    print("  det J' = -a0*b1*y^4, which is never a nonzero constant.")
    print()

    print("Why no shear can repair it:")
    print("  P = a0 + a1*x has the Newton vertex (1,0) with coefficient a1.")
    print("  The inverse map sends x^a y^b to y^(4b-a)(x - s(y))^b.  The")
    print("  monomial (1,0) produces a1/y.  All other terms of P have b = 0")
    print("  only at (0,0), so no term can cancel the y^-1.  The shears act")
    print("  on the second coordinate and cannot touch the first one.")
    print()
    print("No other polynomial clearing helps either:")
    print("  a polynomial representative is (f*P0, g*Q0) with y | f and y^2 | g.")
    print("  Since (P0)_x = 0 and a1*b1 = 1, det(f*P0, g*Q0) = g*(f - b1*y^2*f'*P0).")
    print("  Writing f = y*h gives f - b1*y^2*f'*P0 = y*(h - b1*(a0*y+a1)*(h+y*h')),")
    print("  so the determinant is divisible by y unless it is zero.  Hence no")
    print("  polynomial clearing yields a nonzero constant determinant.")
    print()


def main() -> int:
    print("=" * 78)
    print("GGHV reconstruction audit")
    print(PAPER)
    print("=" * 78)
    print()
    self_test()
    print()

    r = j3_audit()
    print_j3_report(r)

    print("-" * 78)
    print("Concrete J3 representative: P = x, Q = x^2*y")
    rep = representative_audit()
    print(f"  P0 = {rep['P0']},  Q0 = {rep['Q0']}")
    print(f"  cleared: P' = {rep['Pc']},  Q' = {rep['Qc']}")
    print(f"  det J(P',Q') = {sp.factor(rep['det_after'])}")
    print()

    print("-" * 78)
    print("Config-2 note")
    print("  The same y^-1 mechanism is present at the vertex (1,0) of the")
    print("  config-2 polygon {(0,0),(1,0),(8,14),(8,16)}.  Unlike J3, the")
    print("  config-2 lattice also contains (1,1) and (1,2), so the coefficient")
    print("  of y^-1 can in principle be cancelled by the x^-4 shear.  Deeper")
    print("  negative powers (down to y^-8 from the (8,14),(8,16) vertices)")
    print("  give further conditions.  Prop. 4.3 states none of them, and the")
    print("  paper gives no inverse map and no characterization of the image")
    print("  A(K[x,y]) inside L^(1).")
    print()

    print("=" * 78)
    print("VERDICT: PARTIAL")
    print("  The inverse composite exists as a rational map and restores")
    print("  det J = 1 before clearing.  It is not a map to K[x,y]: the")
    print("  inverse swap introduces y^-1.  J3 fails at step 3 with the exact")
    print("  cleared determinant -a0*b1*y^4.  The published data do not supply")
    print("  the shear constants (c2,c3,c4), the inverse of the composite, or")
    print("  the polynomiality/image conditions on the normalized coefficients.")
    print("  Additional data needed: (i) the exact shear constants / branch of")
    print("  the proof; (ii) the vanishing conditions for all negative powers")
    print("  of y in P0 and Q0; (iii) a proof that the config-2 solution lies")
    print("  in A(K[x,y]).  Until then the config-2 system is necessary, not")
    print("  sufficient, and cannot reach the terminal object.")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
