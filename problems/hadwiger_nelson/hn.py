"""Exact unit-distance graph toolkit -- Hadwiger-Nelson track.

A graph is accepted only if every declared edge has squared length exactly 1
under symbolic arithmetic.  The chromatic number is computed exactly by
DSATUR branch-and-bound (fine for graphs up to ~40 vertices; swap in a SAT
solver with DRAT logging for larger instances).

Run `python problems/hadwiger_nelson/hn.py` for the self-check.
"""

from __future__ import annotations

import itertools
import json
from pathlib import Path

import sympy as sp


# --------------------------------------------------------------------------
# exact arithmetic
# --------------------------------------------------------------------------

def dist2(p, q):
    """Exact squared Euclidean distance."""
    return sp.expand(
        sum((sp.sympify(a) - sp.sympify(b)) ** 2 for a, b in zip(p, q))
    )


def is_unit(p, q) -> bool:
    d = sp.simplify(dist2(p, q) - 1)
    if d.is_number:
        return bool(d == 0)
    return False


# --------------------------------------------------------------------------
# graph
# --------------------------------------------------------------------------

class UnitDistanceGraph:
    def __init__(self, points, adj, edges):
        self.points = [tuple(sp.sympify(c) for c in p) for p in points]
        self.adj = list(adj)
        self.edges = list(edges)

    # -- construction ------------------------------------------------------
    @classmethod
    def from_points(cls, points):
        pts = [tuple(sp.sympify(c) for c in p) for p in points]
        n = len(pts)
        adj = [0] * n
        edges = []
        for i, j in itertools.combinations(range(n), 2):
            if is_unit(pts[i], pts[j]):
                adj[i] |= 1 << j
                adj[j] |= 1 << i
                edges.append((i, j))
        return cls(pts, adj, edges)

    # -- invariants --------------------------------------------------------
    @property
    def n(self) -> int:
        return len(self.points)

    @property
    def m(self) -> int:
        return len(self.edges)

    def degrees(self):
        return [bin(a).count("1") for a in self.adj]

    def verify(self) -> tuple[bool, list[str]]:
        """Re-check every declared edge and every non-edge at exact precision."""
        errs = []
        for i, j in self.edges:
            if not is_unit(self.points[i], self.points[j]):
                errs.append(f"edge ({i},{j}) is not unit length")
        edge_set = {tuple(sorted(e)) for e in self.edges}
        for i, j in itertools.combinations(range(self.n), 2):
            if (i, j) not in edge_set and is_unit(self.points[i], self.points[j]):
                errs.append(f"missing edge ({i},{j})")
        return (not errs), errs

    # -- colouring ---------------------------------------------------------
    def _greedy(self) -> int:
        n, adj = self.n, self.adj
        order = sorted(range(n), key=lambda v: -bin(adj[v]).count("1"))
        col = {}
        for v in order:
            used = {col[u] for u in range(n) if ((adj[v] >> u) & 1) and u in col}
            c = 0
            while c in used:
                c += 1
            col[v] = c
        return max(col.values()) + 1 if col else 0

    def chromatic_number(self) -> int:
        """Exact chromatic number (DSATUR + branch and bound)."""
        n, adj = self.n, self.adj
        if n == 0:
            return 0
        best = [self._greedy()]
        color = [-1] * n
        sat = [0] * n

        def rec(colored: int, used: int) -> None:
            if used >= best[0]:
                return
            if colored == n:
                best[0] = used
                return
            v = max(
                (i for i in range(n) if color[i] < 0),
                key=lambda i: (bin(sat[i]).count("1"), bin(adj[i]).count("1")),
            )
            for c in range(used + 1):
                if (sat[v] >> c) & 1:
                    continue
                color[v] = c
                touched = []
                for u in range(n):
                    if ((adj[v] >> u) & 1) and color[u] < 0 and not ((sat[u] >> c) & 1):
                        sat[u] |= 1 << c
                        touched.append(u)
                rec(colored + 1, max(used, c + 1))
                for u in touched:
                    sat[u] &= ~(1 << c)
                color[v] = -1

        rec(0, 0)
        return best[0]

    def coloring(self) -> list[int]:
        """Return one optimal colouring."""
        chi = self.chromatic_number()
        n, adj = self.n, self.adj
        col = [-1] * n

        def ok(v, c):
            return all(not (((adj[v] >> u) & 1) and col[u] == c) for u in range(n))

        def rec(v):
            if v == n:
                return True
            for c in range(chi):
                if ok(v, c):
                    col[v] = c
                    if rec(v + 1):
                        return True
                    col[v] = -1
            return False

        rec(0)
        return col

    def check_coloring(self, colors, k: int | None = None) -> bool:
        """Independent verifier for a claimed colouring."""
        if len(colors) != self.n:
            return False
        for i, j in self.edges:
            if colors[i] == colors[j]:
                return False
        if k is not None and (max(colors) + 1 > k or min(colors) < 0):
            return False
        return True

    # -- serialisation -----------------------------------------------------
    def to_dict(self):
        return {
            "n": self.n,
            "m": self.m,
            "chi": self.chromatic_number(),
            "points": [str(p) for p in self.points],
            "edges": self.edges,
        }

    def save(self, path):
        Path(path).write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")


# --------------------------------------------------------------------------
# canonical small graphs
# --------------------------------------------------------------------------

def moser_spindle() -> UnitDistanceGraph:
    """The 7-vertex, 11-edge, 4-chromatic unit-distance graph.

    Two 'diamonds' (K4 minus an edge, drawn as two equilateral triangles
    sharing an edge) with a common apex O.  Their far tips R1, R2 sit at
    distance sqrt(3) from O and are exactly 1 apart, which forces
    cos(angle R1-O-R2) = 5/6.
    """
    s3, s11, s33 = sp.sqrt(3), sp.sqrt(11), sp.sqrt(33)
    half = sp.Rational(1, 2)
    O = (sp.Integer(0), sp.Integer(0))
    P1 = (s3 / 2, half)
    Q1 = (s3 / 2, -half)
    R1 = (s3, sp.Integer(0))
    P2 = (5 * s3 / 12 - s11 / 12, sp.Rational(5, 12) + s33 / 12)
    Q2 = (5 * s3 / 12 + s11 / 12, -sp.Rational(5, 12) + s33 / 12)
    R2 = (5 * s3 / 6, s33 / 6)
    return UnitDistanceGraph.from_points([O, P1, Q1, R1, P2, Q2, R2])


def triangle() -> UnitDistanceGraph:
    s3 = sp.sqrt(3)
    return UnitDistanceGraph.from_points(
        [(0, 0), (1, 0), (sp.Rational(1, 2), s3 / 2)]
    )


# --------------------------------------------------------------------------
# self-check
# --------------------------------------------------------------------------

if __name__ == "__main__":
    g = moser_spindle()
    ok, errs = g.verify()
    print(f"moser spindle: n={g.n} m={g.m} chi={g.chromatic_number()} verified={ok}")
    for e in errs:
        print("  ERR", e)
    assert ok and g.n == 7 and g.m == 11 and g.chromatic_number() == 4

    t = triangle()
    print(f"triangle:      n={t.n} m={t.m} chi={t.chromatic_number()}")
    assert t.n == 3 and t.m == 3 and t.chromatic_number() == 3
    print("OK")
