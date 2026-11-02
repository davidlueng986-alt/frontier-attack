"""Exact verifier and collision hunter for the Jacobian conjecture.

The conjecture (Keller 1939): if F : C^n -> C^n is polynomial and
det J_F is a non-zero constant, then F has a polynomial inverse.

What this module does
---------------------
* `verify_map` decides, symbolically, whether a proposed counterexample is
  actually a counterexample: constant non-zero determinant, distinct domain
  points, identical image.
* `find_collision` searches one fibre for a second preimage.  This is the only
  operation that can *disprove* the conjecture, and it is cheap to check.

Nothing here proves the planar case.  It is a gate that rejects bad claims
before a human ever reads them.

Run `python problems/jacobian/jac.py` for the self-check.
"""

from __future__ import annotations

import random

import sympy as sp


def _finite(expr) -> bool:
    e = sp.simplify(expr)
    return bool(e.is_number)


def verify_map(F, vars, collision_points=None) -> dict:
    """Return a verdict dict for a candidate counterexample."""
    F = list(F)
    vars = list(vars)
    reasons: list[str] = []

    J = sp.Matrix(F).jacobian(vars)
    det = sp.simplify(J.det())

    if not _finite(det):
        return {"verdict": "REJECTED",
                "det": str(det),
                "reasons": ["Jacobian determinant is not constant"]}
    if sp.simplify(det) == 0:
        return {"verdict": "REJECTED", "det": str(det),
                "reasons": ["Jacobian determinant is zero"]}

    if collision_points is None:
        return {"verdict": "INCOMPLETE", "det": str(det),
                "reasons": ["constant non-zero determinant, but no collision supplied"]}

    pts = [tuple(sp.sympify(c) for c in p) for p in collision_points]
    if len(pts) != 2:
        reasons.append("need exactly two domain points")
    if len(set(pts)) != len(pts):
        reasons.append("domain points are not distinct")
    if len(pts[0]) != len(vars):
        reasons.append("point dimension does not match the number of variables")

    if not reasons:
        imgs = []
        for p in pts:
            sub = dict(zip(vars, p))
            imgs.append(tuple(sp.simplify(f.subs(sub)) for f in F))
        same = all(sp.simplify(a - b) == 0 for a, b in zip(imgs[0], imgs[1]))
        if not same:
            reasons.append(f"images differ: {imgs[0]} != {imgs[1]}")
        elif any(sp.simplify(a) == sp.zoo or a.has(sp.nan) for a in imgs[0]):
            reasons.append("image is not finite")

    if reasons:
        return {"verdict": "REJECTED", "det": str(det), "reasons": reasons}
    return {"verdict": "COUNTEREXAMPLE", "det": str(det), "reasons": []}


def find_collision(F, vars, samples=12, starts=6, seed=0, tol=1e-8):
    """Search one fibre at a time for a second preimage.

    Purely numerical (mpmath 30 digits via sympy.nsolve).  A hit is a *lead*,
    not a proof -- feed it to `verify_map` for the exact check.
    """
    rng = random.Random(seed)
    F = list(F)
    vars = list(vars)
    hits = []
    for _ in range(samples):
        p = tuple(sp.Float(rng.uniform(-1.5, 1.5), 30) for _ in vars)
        target = [sp.simplify(f.subs(dict(zip(vars, p)))) for f in F]
        eqs = [sp.expand(f - t) for f, t in zip(F, target)]
        for _ in range(starts):
            guess = [sp.Float(rng.uniform(-2, 2), 30) for _ in vars]
            try:
                sol = sp.nsolve(eqs, vars, guess, prec=30, maxsteps=60)
            except Exception:
                continue
            sol = tuple(sp.Float(s, 30) for s in sol)
            far = max(abs(float(a - b)) for a, b in zip(sol, p))
            if far > tol:
                hits.append({"p": [str(sp.N(c, 20)) for c in p],
                             "q": [str(sp.N(c, 20)) for c in sol],
                             "separation": far})
                break
    return hits


def alpoge_counterexample():
    """The July 2026 C^3 counterexample (Alpoge, credited to Claude Fable 5).

    Component degrees (7, 6, 4), total degree 7, det J == -2, three distinct
    points with the same image.  This is a real object, not a synthetic test:
    it is the regression that proves the verifier accepts a genuine
    counterexample.
    """
    x, y, z = sp.symbols("x y z")
    F = [
        (1 + x * y) ** 3 * z + y ** 2 * (1 + x * y) * (4 + 3 * x * y),
        y + 3 * x * (1 + x * y) ** 2 * z + 3 * x * y ** 2 * (4 + 3 * x * y),
        2 * x - 3 * x ** 2 * y - x ** 3 * z,
    ]
    pts = [
        (0, 0, sp.Rational(-1, 4)),
        (1, sp.Rational(-3, 2), sp.Rational(13, 2)),
        (-1, sp.Rational(3, 2), sp.Rational(13, 2)),
    ]
    return F, [x, y, z], pts


# --------------------------------------------------------------------------
# self-check
# --------------------------------------------------------------------------

def run_claim(path: str) -> int:
    """CLI gate: `python jac.py --claim claim.json` -> exit 0 iff COUNTEREXAMPLE."""
    import json
    from pathlib import Path

    spec = json.loads(Path(path).read_text(encoding="utf-8"))
    syms = sp.symbols(" ".join(spec["vars"]))
    vars_ = list(syms) if isinstance(syms, tuple) else [syms]
    F = [sp.sympify(s) for s in spec["F"]]
    pts = [tuple(sp.sympify(c) for c in p) for p in spec["points"]]
    result = verify_map(F, vars_, pts)
    print(json.dumps(result, indent=2))
    return 0 if result["verdict"] == "COUNTEREXAMPLE" else 1


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 2 and sys.argv[1] == "--claim":
        raise SystemExit(run_claim(sys.argv[2]))

    x, y = sp.symbols("x y")

    # 1. triangular map: det = 1, invertible, must be INCOMPLETE not a counterexample
    F1 = [x + y ** 2, y]
    r1 = verify_map(F1, [x, y])
    print("triangular :", r1)
    assert r1["verdict"] == "INCOMPLETE"

    # 2. non-constant determinant: must be REJECTED even if a collision is supplied
    F2 = [x ** 2, y]
    r2 = verify_map(F2, [x, y], collision_points=[(1, 0), (-1, 0)])
    print("x^2, y     :", r2)
    assert r2["verdict"] == "REJECTED"

    # 3. collision hunter finds the +-1 pair for x^2
    hits = find_collision(F2, [x, y], samples=4, starts=6, seed=3)
    print("collision  :", hits[:1])
    assert hits

    # 4. triangular map should yield no collision
    assert find_collision(F1, [x, y], samples=4, starts=4, seed=5) == []

    # 5. a fabricated "counterexample" with constant det but wrong points must fail
    r5 = verify_map(F1, [x, y], collision_points=[(0, 0), (1, 1)])
    print("bad claim  :", r5)
    assert r5["verdict"] == "REJECTED"

    # 6. regression: the real 2026 C^3 counterexample must be accepted
    F6, v6, pts6 = alpoge_counterexample()
    r6 = verify_map(F6, v6, pts6[:2])
    print("Alpoge 2026:", r6["verdict"], "det =", r6["det"])
    assert r6["verdict"] == "COUNTEREXAMPLE" and r6["det"] == "-2"
    r6b = verify_map(F6, v6, pts6[1:])
    assert r6b["verdict"] == "COUNTEREXAMPLE"

    print("OK")
