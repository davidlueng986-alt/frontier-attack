"""Random search for planar cubic maps with constant Jacobian.

Samples support patterns (three monomials per component) and small integer
coefficients.  A numeric determinant prescreen keeps candidates whose Jacobian
determinant is numerically constant; survivors are handed to the exact
verifier and to the collision hunter.

Run `python problems/jacobian/cubic_search.py --samples 2000`.
"""

from __future__ import annotations

import argparse
import itertools
import json
import random
import sys
import time
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))

from jac import find_collision, verify_map  # noqa: E402


MONOMIALS = (
    lambda x, y: x**3,
    lambda x, y: x**2 * y,
    lambda x, y: x * y**2,
    lambda x, y: y**3,
    lambda x, y: x**2,
    lambda x, y: x * y,
    lambda x, y: y**2,
)


def _build_pattern(x, y, s1, s2):
    """Return (f, g, symbols, numeric determinant) for one support pattern."""
    cs = sp.symbols(f"c0:{2 * len(s1)}")
    terms1 = [c * MONOMIALS[i](x, y) for c, i in zip(cs[:len(s1)], s1)]
    terms2 = [c * MONOMIALS[i](x, y) for c, i in zip(cs[len(s1):], s2)]
    f = x + sum(terms1)
    g = y + sum(terms2)
    det = sp.expand(sp.Matrix([f, g]).jacobian([x, y]).det())
    return f, g, cs, sp.lambdify(cs, det, "math")


def search(samples: int, seed: int):
    x, y = sp.symbols("x y")
    rng = random.Random(seed)
    prescreen = random.Random(seed ^ 0x5EED)

    supports = list(itertools.combinations(range(len(MONOMIALS)), 3))
    pattern_cache: dict[tuple, tuple] = {}
    points = [(prescreen.uniform(-1.5, 1.5), prescreen.uniform(-1.5, 1.5))
              for _ in range(6)]

    survivors = []
    seen = set()
    t0 = time.time()
    for _ in range(samples):
        s1 = rng.choice(supports)
        s2 = rng.choice(supports)
        coeffs = tuple(rng.randint(-3, 3) for _ in range(6))
        key = (s1, s2, coeffs)
        if key in seen:
            continue
        seen.add(key)

        pattern = pattern_cache.get((s1, s2))
        if pattern is None:
            pattern = _build_pattern(x, y, s1, s2)
            pattern_cache[(s1, s2)] = pattern
        _, _, cs, f_det = pattern

        try:
            vals = [float(f_det(*coeffs)) for _ in points]
        except Exception:
            continue
        if max(vals) - min(vals) > 1e-9:
            continue

        survivors.append({"s1": list(s1), "s2": list(s2),
                          "coeffs": list(coeffs)})

    prescreen_seconds = time.time() - t0

    constant_det = 0
    counterexamples = 0
    for survivor in survivors:
        f, g, cs, _ = pattern_cache[(tuple(survivor["s1"]),
                                     tuple(survivor["s2"]))]
        sub = dict(zip(cs, survivor["coeffs"]))
        F = [sp.expand(f.subs(sub)), sp.expand(g.subs(sub))]
        verdict = verify_map(F, [x, y])
        hits = find_collision(F, [x, y], samples=3, starts=4, seed=seed)
        survivor["F"] = [str(F[0]), str(F[1])]
        survivor["verdict"] = verdict["verdict"]
        survivor["det"] = verdict["det"]
        survivor["collisions"] = len(hits)
        if verdict["verdict"] == "INCOMPLETE":
            constant_det += 1
        if verdict["verdict"] == "INCOMPLETE" and hits:
            counterexamples += 1

    return {
        "samples": samples,
        "constant_det": constant_det,
        "counterexamples": counterexamples,
        "survivors": survivors,
        "prescreen_seconds": round(prescreen_seconds, 2),
        "seconds": round(time.time() - t0, 2),
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args(argv)

    result = search(args.samples, args.seed)
    Path("runs").mkdir(exist_ok=True)
    Path("runs/jacobian_cubic.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8")
    print(f"cubic search: samples={result['samples']} "
          f"constant_det={result['constant_det']} "
          f"counterexamples={result['counterexamples']} "
          f"survivors={len(result['survivors'])} "
          f"seconds={result['seconds']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
