"""Independent verifier for the 509-vertex unit-distance graph (Parts 2020).

The published file is a list of 509 points in Mathematica syntax, one per line:

    {0, 0}
    {(-5 + Sqrt[33])/12, (3*Sqrt[3] - Sqrt[11])/12}

There is no edge list: edges are exactly the pairs at distance 1.  This script
parses the coordinates independently, rebuilds the edge set at 50 digits, and
decides k-colourability with a SAT solver.

Source of the file (not vendored):
    https://github.com/vasnesterov/HadwigerNelson  ->  vtx/509_parts.vtx

Usage:
    python tools/verify_509.py research/scratch/509_parts.vtx --k 4
"""

from __future__ import annotations

import argparse
import re
import time

import mpmath as mp
import sympy as sp
from pysat.formula import CNF
from pysat.solvers import Solver

mp.mp.dps = 50


def split_top(body: str):
    depth = 0
    for i, ch in enumerate(body):
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
        elif ch == "," and depth == 0:
            return body[:i], body[i + 1:]
    raise ValueError("no top-level comma")


def parse_point(line: str):
    s = line.strip()
    if not s:
        return None
    if not (s.startswith("{") and s.endswith("}")):
        raise ValueError(f"bad line: {s[:40]}")
    a, b = split_top(s[1:-1])
    out = []
    for expr in (a, b):
        # Mathematica applies functions with square brackets: Sqrt[33]
        e = expr
        while "Sqrt[" in e:
            nxt = re.sub(r"Sqrt\[([^\[\]]*)\]", r"sqrt(\1)", e)
            if nxt == e:
                break
            e = nxt
        out.append(mp.mpf(str(sp.N(sp.sympify(e), 50))))
    return tuple(out)


def build(points, tol=mp.mpf("1e-40")):
    n = len(points)
    edges = []
    worst = mp.mpf(0)
    for i in range(n):
        xi, yi = points[i]
        for j in range(i + 1, n):
            xj, yj = points[j]
            d2 = (xi - xj) ** 2 + (yi - yj) ** 2
            dev = abs(d2 - 1)
            if dev <= tol:
                edges.append((i, j))
                if dev > worst:
                    worst = dev
    return edges, worst


def colouring_cnf(n, edges, k):
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
    return cnf


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("path")
    ap.add_argument("--k", type=int, default=4)
    ap.add_argument("--tol", default="1e-40")
    a = ap.parse_args()

    pts = []
    with open(a.path, encoding="utf-8") as f:
        for line in f:
            p = parse_point(line)
            if p is not None:
                pts.append(p)

    edges, worst = build(pts, mp.mpf(a.tol))
    print(f"vertices = {len(pts)}")
    print(f"edges    = {len(edges)}")
    print(f"worst |d^2 - 1| on an edge = {mp.nstr(worst, 6)}")

    cnf = colouring_cnf(len(pts), edges, a.k)
    print(f"CNF      = {len(cnf.clauses)} clauses, k = {a.k}")
    t = time.time()
    with Solver(name="cadical153", bootstrap_with=cnf) as s:
        sat = s.solve()
    print(f"{'SATISFIABLE' if sat else 'UNSATISFIABLE'} in {time.time() - t:.1f}s")
    return 0 if not sat else 1


if __name__ == "__main__":
    raise SystemExit(main())
