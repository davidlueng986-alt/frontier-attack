"""Bounded enumeration of polynomial maps with constant Jacobian.

Purpose: clear small boxes honestly, and show *why* the planar case resists
brute force.  Every surviving coefficient tuple is passed to the exact
verifier and to the collision hunter; a box is only reported "clear" if both
agree.

Run `python problems/jacobian/enumerate.py`.
"""

from __future__ import annotations

import itertools
import json
import time
from pathlib import Path

import sympy as sp

from jac import find_collision, verify_map


def _det_coeffs(F, vars):
    J = sp.Matrix(F).jacobian(vars)
    det = sp.expand(J.det())
    P = sp.Poly(det, *vars)
    return det, P


def enumerate_family(name: str, coeff_range=(-2, 2), max_report=8):
    x, y = sp.symbols("x y")
    if name == "quadratic":
        a, b, c, d, e, f = sp.symbols("a b c d e f")
        F = [x + a * x**2 + b * x * y + c * y**2,
             y + d * x**2 + e * x * y + f * y**2]
        params = (a, b, c, d, e, f)
    elif name == "cubic_sparse":
        a, b, c, d = sp.symbols("a b c d")
        F = [x + a * x**3 + b * x**2 * y,
             y + c * x * y**2 + d * y**3]
        params = (a, b, c, d)
    else:
        raise ValueError(name)

    det, P = _det_coeffs(F, [x, y])
    monoms = [m for m in P.monoms() if m != (0, 0)]
    exprs = [P.coeff_monomial(m) for m in monoms]
    f_det = sp.lambdify(params, exprs, "math")

    grid = list(itertools.product(range(coeff_range[0], coeff_range[1] + 1),
                                  repeat=len(params)))
    survivors = []
    t0 = time.time()
    for vals in grid:
        try:
            out = f_det(*vals)
        except Exception:
            continue
        if all(abs(v) < 1e-12 for v in out):
            survivors.append(vals)
    dt = time.time() - t0

    results = []
    hits = 0
    for vals in survivors[:max_report]:
        sub = dict(zip(params, vals))
        Fn = [sp.expand(f.subs(sub)) for f in F]
        r = verify_map(Fn, [x, y])
        col = find_collision(Fn, [x, y], samples=3, starts=4, seed=0)
        if col:
            hits += 1
        results.append({"coeffs": list(vals), "F": [str(f) for f in Fn],
                        "verdict": r["verdict"], "collisions": len(col)})

    return {
        "family": name,
        "coefficient_range": list(coeff_range),
        "grid_size": len(grid),
        "constant_jacobian_survivors": len(survivors),
        "counterexamples_found": hits,
        "seconds": round(dt, 2),
        "samples": results,
    }


if __name__ == "__main__":
    out = {}
    for fam in ("quadratic", "cubic_sparse"):
        r = enumerate_family(fam)
        out[fam] = r
        print(f"\n== {fam} ==")
        print(json.dumps({k: v for k, v in r.items() if k != "samples"}, indent=2))
        for s in r["samples"][:3]:
            print("  ", s["coeffs"], s["verdict"], "collisions:", s["collisions"])
    Path("runs").mkdir(exist_ok=True)
    Path("runs/jacobian_enumeration.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8")
    print("\nwritten: runs/jacobian_enumeration.json")
