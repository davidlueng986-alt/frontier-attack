"""Acceptance gate: does a claim meet the contract in docs/DEFINITION_OF_SOLVED.md?

Exit 0 only if every clause passes. This is deliberately stricter than the
individual checkers: a true object that is off-target is rejected here.

Usage:
    python tools/accept.py --problem jac --claim claim.json
    python tools/accept.py --problem hn  --claim graph.json
    python tools/accept.py --problem bsd --claim curve.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "problems" / "jacobian"))
sys.path.insert(0, str(ROOT / "problems" / "bsd"))

import sympy as sp  # noqa: E402


def fail(reasons):
    print(json.dumps({"accepted": False, "reasons": reasons}, indent=2))
    return 1


def ok(extra=None):
    print(json.dumps({"accepted": True, **(extra or {})}, indent=2))
    return 0


# --------------------------------------------------------------------------
# Hadwiger-Nelson
# --------------------------------------------------------------------------

def accept_hn(path: Path):
    from pysat.formula import CNF
    from pysat.solvers import Solver

    spec = json.loads(path.read_text(encoding="utf-8"))
    pts = [tuple(sp.sympify(c) for c in p) for p in spec["points"]]
    n = len(pts)
    if n == 0:
        return fail(["empty graph"])

    edges = []
    bad = []
    for i in range(n):
        for j in range(i + 1, n):
            d = sp.simplify(
                sum((a - b) ** 2 for a, b in zip(pts[i], pts[j])) - 1
            )
            if d.is_number and d == 0:
                edges.append((i, j))
            elif d.is_number:
                continue
            else:
                bad.append([i, j])
    if bad:
        return fail([f"{len(bad)} pairs were not exact numbers; usage of floats is rejected"])

    # A3/A4: not 6-colourable, with a certificate.
    k = 6

    def var(v, c):
        return v * k + c + 1

    cnf = CNF()
    for v in range(n):
        cnf.append([var(v, c) for c in range(k)])
        for c1 in range(k):
            for c2 in range(c1 + 1, k):
                cnf.append([-var(v, c1), -var(v, c2)])
    for u, v in edges:
        for c in range(k):
            cnf.append([-var(u, c), -var(v, c)])

    with Solver(name="cadical153", bootstrap_with=cnf) as s:
        sat = s.solve()

    proof = spec.get("proof", {})
    cert = proof.get("path")
    cert_ok = bool(cert) and (ROOT / cert).exists() if cert else False

    if sat:
        return fail([f"graph is {k}-colourable; chi <= {k}, so this is not a solution"])
    if not cert_ok:
        return fail([
            "UNSAT confirmed in-process, but no DRAT/LRAT certificate was supplied",
            "DEFINITION_OF_SOLVED.md A4 requires a checkable trace",
        ])
    return ok({"n": n, "m": len(edges), "chi_at_least": 7, "certificate": cert})


# --------------------------------------------------------------------------
# Jacobian
# --------------------------------------------------------------------------

def _total_degree(expr, vars):
    if expr == 0:
        return -1
    return max(sum(m) for m in sp.Poly(expr, *vars).monoms())


def _is_obvious_triangular(F, vars):
    """Reject the two shear forms. Not a general composition test."""
    x, y = vars
    f, g = F
    if sp.simplify(f - x) == 0:
        return True                      # (x, y + h(x))
    if sp.simplify(g - y) == 0:
        d = sp.Poly(sp.expand(f - x), x)
        if d.degree() <= 0:
            return True                  # (x + h(y), y)
    return False


def accept_jac(path: Path):
    from jac import verify_map

    spec = json.loads(path.read_text(encoding="utf-8"))
    syms = sp.symbols(" ".join(spec["vars"]))
    vars_ = list(syms) if isinstance(syms, tuple) else [syms]
    F = [sp.sympify(s) for s in spec["F"]]
    pts = [tuple(sp.sympify(c) for c in p) for p in spec["points"]]

    r = verify_map(F, vars_, pts)
    if r["verdict"] != "COUNTEREXAMPLE":
        return fail([f"verify_map: {r['verdict']}", *r["reasons"]])

    degs = [_total_degree(f, vars_) for f in F]
    if max(degs) < 125:
        return fail([f"degree pair {degs} is below the proved bound 125; the box is empty"])

    if _is_obvious_triangular(F, vars_):
        return fail(["triangular map: invertible by construction"])

    return ok({"det": r["det"], "degrees": degs})


# --------------------------------------------------------------------------
# BSD
# --------------------------------------------------------------------------

def accept_bsd(path: Path):
    from bsd import check_claim

    claim = json.loads(path.read_text(encoding="utf-8"))
    r = check_claim(claim)
    if r["verdict"] != "COUNTEREXAMPLE_CANDIDATE":
        return fail([f"check_claim: {r['verdict']}", *r["reasons"]])

    k = claim.get("analytic_rank")
    if not isinstance(k, int) or k > 3:
        return fail(["analytic rank > 3 is not proved for any curve; k <= 3 required"])
    if claim.get("rank_method") not in {"2-descent", "mwrank", "magma", "sage", "published"}:
        return fail(["algebraic rank has no recognised backend"])
    if claim.get("analytic_rank_method") not in {"lcalc", "numerical", "published"}:
        return fail(["analytic rank has no recognised backend"])

    return ok({"rank": claim["rank"], "analytic_rank": k})


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--problem", choices=["hn", "jac", "bsd"], required=True)
    ap.add_argument("--claim", required=True)
    a = ap.parse_args()
    fn = {"hn": accept_hn, "jac": accept_jac, "bsd": accept_bsd}[a.problem]
    raise SystemExit(fn(Path(a.claim)))
