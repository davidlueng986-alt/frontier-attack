#!/usr/bin/env python3
"""A verified 6-chromatic two-distance graph in the plane.

Source
------
J. De Neve, F. Vanden Kerchove, D. Colle, W. Tavernier, M. Pickavet,
"On the chromatic number of the plane with two forbidden distances",
American Mathematical Monthly 132(10) (2025) 1007-1022,
DOI 10.1080/00029890.2025.2559554.

The coordinates are published in the Appendix, Table 1 (p. 15): all 33
vertices of the smallest graph for d = sqrt(3), in the condensed notation
[a0 a1 a2 a3] meaning

    (x, y) = ( (a0 + a1*sqrt(5))/8 , (a2*sqrt(3) + a3*sqrt(15))/8 ).

The paper does not tabulate the edges.  Its Definition 1 makes the graph a
distance graph with distance set D: an edge exists *iff* the distance between
the two vertices lies in D.  So the edge set is recovered from the coordinates,
and with D = {1, sqrt(3)} it has 113 edges (81 unit + 32 of length sqrt(3)),
matching the paper's "33 vertices and 113 edges" (Figure 11, p. 12).

Terminology note
----------------
This is a two-distance *graph* (two possible edge lengths), which is what the
paper and the rest of the literature mean by the term.  It is not a strict
two-distance *set*: the 33 points realise 39 distinct pairwise distances.  No
6-chromatic strict two-distance set can exist in the plane, since a strict
two-distance set in R^2 has at most 5 points (Einhorn-Schoenberg 1966) while a
6-chromatic graph needs at least 6 vertices.  The verification below therefore
checks edge lengths, which is exactly the invariant a virtual-edge
substitution consumes.

What is verified here
---------------------
1. 33 distinct points, 113 edges.
2. Every edge has squared length exactly 1 or 3 (exact arithmetic in
   Q(sqrt(3), sqrt(5), sqrt(15)), cross-checked with sympy).
3. Every pair at squared distance 1 or 3 is present as an edge.
4. chi = 6: no 5-colouring (UNSAT from two independent SAT solvers) and an
   explicit, re-checked 6-colouring.
5. The 32 non-unit (length sqrt(3)) edges are listed.

Run:  python problems/hadwiger_nelson/two_distance.py
Writes runs/two_distance.json.  Exits 0 either way; on failure it prints
"NOT OBTAINED: <reason>" instead of writing a graph.
"""

from __future__ import annotations

import itertools
import json
import sys
import time
from pathlib import Path

import sympy as sp
from pysat.solvers import Cadical195, Glucose4

ROOT = Path(__file__).resolve().parents[2]
OUT_PATH = ROOT / "runs" / "two_distance.json"

SOURCE = (
    "J. De Neve, F. Vanden Kerchove, D. Colle, W. Tavernier, M. Pickavet, "
    "'On the chromatic number of the plane with two forbidden distances', "
    "Amer. Math. Monthly 132(10) (2025) 1007-1022, "
    "DOI 10.1080/00029890.2025.2559554; Appendix Table 1 (p. 15) and Figure 11 (p. 12)"
)

# --------------------------------------------------------------------------
# Table 1 of the appendix, transcribed verbatim (33 rows).
# --------------------------------------------------------------------------

TABLE_1 = (
    (0, 0, 0, 0),
    (0, 0, 0, 4),
    (1, 3, -1, 1),
    (1, 3, 1, 3),
    (-1, 3, 1, 1),
    (-1, 3, -1, 3),
    (1, -3, 1, 1),
    (1, -3, -1, 3),
    (2, 0, 0, 2),
    (0, 0, 2, 2),
    (0, 0, -2, 2),
    (-1, -3, 1, 3),
    (-2, 0, 0, 2),
    (-2, 0, 2, 4),
    (-1, 3, 3, 3),
    (-3, 3, -1, 1),
    (-3, 3, 1, 3),
    (3, 3, -1, 3),
    (1, 3, 3, 1),
    (3, -3, 1, 3),
    (3, -3, -1, 1),
    (4, 0, 0, 4),
    (4, 0, 0, 0),
    (-3, -3, -1, 3),
    (-3, -3, 1, 1),
    (-2, 0, -2, 4),
    (-2, 0, 2, 0),
    (-4, 0, 0, 4),
    (-4, 0, 0, 0),
    (1, -3, 3, 3),
    (1, -3, -3, 1),
    (-1, -3, -3, 3),
    (-1, -3, 3, 1),
)

DEN = 8  # coordinates have denominator 8


# --------------------------------------------------------------------------
# exact arithmetic in Q(sqrt(3), sqrt(5), sqrt(15))
# --------------------------------------------------------------------------
# Elements are 4-tuples (a, b, c, d) for a + b*sqrt(3) + c*sqrt(5) + d*sqrt(15),
# with integer numerators over a common denominator.  This makes every
# distance comparison exact and avoids floating point entirely.

def _mul(u, v):
    a1, b1, c1, d1 = u
    a2, b2, c2, d2 = v
    return (
        a1 * a2 + 3 * b1 * b2 + 5 * c1 * c2 + 15 * d1 * d2,
        a1 * b2 + b1 * a2 + 5 * c1 * d2 + 5 * d1 * c2,
        a1 * c2 + c1 * a2 + 3 * b1 * d2 + 3 * d1 * b2,
        a1 * d2 + d1 * a2 + b1 * c2 + c1 * b2,
    )


def _coords_num(a0, a1, a2, a3):
    """Numerator of (x, y) over the common denominator DEN, in the field basis."""
    return ((a0, 0, a1, 0), (0, a2, 0, a3))


def _sqdist_num(p, q):
    """DEN^2 * squared distance, as an exact field element."""
    dx = tuple(a - b for a, b in zip(p[0], q[0]))
    dy = tuple(a - b for a, b in zip(p[1], q[1]))
    s = _mul(dx, dx)
    t = _mul(dy, dy)
    return tuple(a + b for a, b in zip(s, t))


def _is_scalar(u, value):
    return u == (value, 0, 0, 0)


def sympy_points():
    """The 33 points as exact sympy expressions."""
    s3, s5, s15 = sp.sqrt(3), sp.sqrt(5), sp.sqrt(15)
    return [
        ((sp.Integer(a0) + a1 * s5) / DEN, (a2 * s3 + a3 * s15) / DEN)
        for a0, a1, a2, a3 in TABLE_1
    ]


def sympy_sqdist(p, q):
    d = sp.expand((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2)
    return sp.nsimplify(sp.simplify(d))


# --------------------------------------------------------------------------
# graph construction and verification
# --------------------------------------------------------------------------

def build_graph():
    """Recover the edge set exactly as the paper defines it: edge iff d in D."""
    raw = [_coords_num(*row) for row in TABLE_1]
    n = len(raw)
    unit, d3 = [], []
    spectrum = {}
    for i, j in itertools.combinations(range(n), 2):
        s = _sqdist_num(raw[i], raw[j])
        spectrum[s] = spectrum.get(s, 0) + 1
        if _is_scalar(s, DEN * DEN):
            unit.append((i, j))
        elif _is_scalar(s, DEN * DEN * 3):
            d3.append((i, j))
    edges = sorted(unit + d3)
    return raw, edges, unit, d3, spectrum


def verify(raw, edges, unit, d3):
    """Exact re-check, independent of the construction pass."""
    errors = []

    # 1. distinct points
    if len(set(raw)) != len(raw):
        errors.append("duplicate vertex positions")

    n = len(raw)
    edge_set = {tuple(sorted(e)) for e in edges}
    if len(edge_set) != len(edges):
        errors.append("duplicate edges")

    # 2. every edge is exactly 1 or sqrt(3), and every such pair is an edge
    for i, j in itertools.combinations(range(n), 2):
        s = _sqdist_num(raw[i], raw[j])
        is_unit = _is_scalar(s, DEN * DEN)
        is_d3 = _is_scalar(s, DEN * DEN * 3)
        present = (i, j) in edge_set
        if (is_unit or is_d3) != present:
            errors.append(f"edge/pair mismatch at ({i},{j}): sqdist={s}")
        if present and not (is_unit or is_d3):
            errors.append(f"third edge length at ({i},{j}): sqdist={s}")

    # 3. the counts the paper reports
    if (n, len(edges), len(unit), len(d3)) != (33, 113, 81, 32):
        errors.append(
            f"count mismatch: n={n} m={len(edges)} unit={len(unit)} d3={len(d3)}"
        )

    # 4. independent sympy cross-check over all pairs
    pts = sympy_points()
    for i, j in itertools.combinations(range(n), 2):
        d = sympy_sqdist(pts[i], pts[j])
        if d == 1:
            if (i, j) not in edge_set:
                errors.append(f"sympy: missing unit edge ({i},{j})")
        elif d == 3:
            if (i, j) not in edge_set:
                errors.append(f"sympy: missing sqrt(3) edge ({i},{j})")
        elif (i, j) in edge_set:
            errors.append(f"sympy: edge ({i},{j}) has sqdist {d}")

    return errors


# --------------------------------------------------------------------------
# SAT: chi = 6
# --------------------------------------------------------------------------

def _make_solver(name):
    return {"cadical195": Cadical195, "glucose4": Glucose4}[name]()


def colour_sat(k, n, edges, solver_name):
    """Return a colour list if k-colourable, else None."""
    solver = _make_solver(solver_name)
    var = lambda v, c: v * k + c + 1

    for v in range(n):
        solver.add_clause([var(v, c) for c in range(k)])
        for c1 in range(k):
            for c2 in range(c1 + 1, k):
                solver.add_clause([-var(v, c1), -var(v, c2)])

    for i, j in edges:
        for c in range(k):
            solver.add_clause([-var(i, c), -var(j, c)])

    ok = solver.solve()
    model = solver.get_model() if ok else None
    solver.delete()
    if not ok:
        return None

    positive = {x for x in model if x > 0}
    colours = [next(c for c in range(k) if var(v, c) in positive) for v in range(n)]
    return colours


def verify_colouring(colours, n, edges, k):
    if colours is None or len(colours) != n:
        return False
    if any(not (0 <= c < k) for c in colours):
        return False
    return all(colours[i] != colours[j] for i, j in edges)


# --------------------------------------------------------------------------

def main() -> int:
    t0 = time.time()
    print("6-chromatic two-distance graph -- De Neve et al. (2025), d = sqrt(3)")

    raw, edges, unit, d3, spectrum = build_graph()
    n, m = len(raw), len(edges)

    errors = verify(raw, edges, unit, d3)
    if errors:
        for e in errors[:10]:
            print("  FAIL:", e)
        print(f"NOT OBTAINED: {len(errors)} verification error(s) in the transcribed graph")
        return 0

    # point-set distance spectrum (diagnostic; not an edge-length check)
    pair_count = n * (n - 1) // 2
    other_pairs = pair_count - m
    print(f"n={n}  m={m}  (unit={len(unit)}, sqrt(3)={len(d3)})")
    print(f"edge lengths: exactly {{1, sqrt(3)}} -- no third edge length")
    print(
        f"all {pair_count} pairs: {len(spectrum)} distinct squared distances, "
        f"{other_pairs} non-edge pairs at other lengths (allowed: edge iff d in D)"
    )

    # --- chi = 6 ---------------------------------------------------------
    colours6 = colour_sat(6, n, edges, "cadical195")
    if not verify_colouring(colours6, n, edges, 6):
        print("NOT OBTAINED: no valid 6-colouring produced")
        return 0

    sat5 = {}
    for name in ("cadical195", "glucose4"):
        ts = time.time()
        sat5[name] = colour_sat(5, n, edges, name) is not None
        print(f"  {name}: 5-colouring {'SAT' if sat5[name] else 'UNSAT'} "
              f"({time.time() - ts:.2f}s)")

    if any(sat5.values()):
        print("NOT OBTAINED: graph is 5-colourable, so chi < 6")
        return 0
    if not verify_colouring(colour_sat(6, n, edges, "glucose4"), n, edges, 6):
        print("NOT OBTAINED: 6-colouring not reproducible with the second solver")
        return 0

    print("chi = 6 (5-colouring UNSAT in two independent solvers; 6-colouring SAT)")

    # --- serialise -------------------------------------------------------
    pts = sympy_points()
    payload = {
        "source": SOURCE,
        "n": n,
        "m": m,
        "d2_values": ["sqrt(3)"],
        "d2_squared": [3],
        "non_unit_edges": [[i, j] for i, j in d3],
        "chi": 6,
        "points": [[str(x), str(y)] for x, y in pts],
        "points_float": [
            [round(float(sp.N(x, 15)), 12), round(float(sp.N(y, 15)), 12)]
            for x, y in pts
        ],
        "points_table1": [list(row) for row in TABLE_1],
        "edges": [[i, j] for i, j in edges],
        "unit_edges": [[i, j] for i, j in unit],
        "edge_length_counts": {"1": len(unit), "sqrt(3)": len(d3)},
        "colouring": colours6,
        "verification": {
            "distinct_points": n,
            "pair_count": pair_count,
            "distinct_squared_distances_over_all_pairs": len(spectrum),
            "non_edge_pairs_at_other_lengths": other_pairs,
            "strict_two_distance_set": False,
            "note": (
                "Two-distance GRAPH in the sense of the source's Definition 1 "
                "(edge iff distance is in D={1,sqrt(3)}); not a strict "
                "two-distance point set. A strict two-distance set in R^2 has "
                "at most 5 points (Einhorn-Schoenberg 1966), so no 6-chromatic "
                "strict two-distance set exists."
            ),
            "edge_lengths_exact": True,
            "sat_5colouring": "UNSAT (cadical195, glucose4)",
            "sat_6colouring": "SAT (validated)",
        },
    }

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(json.dumps(payload, indent=1) + "\n", encoding="utf-8")

    print(f"non-unit edges: {len(d3)} (all length sqrt(3))")
    print(f"wrote {OUT_PATH}")
    print(f"done in {time.time() - t0:.2f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
