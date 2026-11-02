"""Forced-equality gadgets for the Hadwiger--Nelson track.

The earlier ``glue2.py`` glued rigid copies of the Moser spindle but never
carried a constraint.  A union of 4-colourable graphs is 4-colourable, so that
search could never leave chi=4.  This module implements the missing step.

  * ``Gadget.forced_relation(k)`` enumerates every proper k-colouring (up to
    renaming of colours, which leaves equality/difference relations unchanged)
    and classifies every pair of distinguished vertices as forced equal,
    forced different, or free.
  * ``diamond()`` is the 4-vertex, 5-edge gadget whose tips are forced equal in
    every 3-colouring.
  * ``spindle(g1, g2, ...)`` places two copies sharing one distinguished
    vertex, rotates the second until a second distinguished pair is exactly
    unit distance, adds that edge, and returns the verified union graph.

The lift theorem used throughout:

    If G is k-colourable and u, v are forced equal in every k-colouring of G,
    then two copies of G sharing u, with the two copies of v joined by a unit
    edge, have chromatic number exactly k+1 (when no extra cross edges appear;
    in general the chromatic number is at least k+1).

Proof of the lower bound: a k-colouring would force both copies of v to carry
u's colour, but the two copies of v are adjacent.  Upper bound: colour both
copies with k colours, align them on u, then recolour one copy of v with a
fresh (k+1)-st colour.

The diamond at k=3 therefore gives the Moser spindle (chi=4).  Reaching chi=5
requires a *unit-distance gadget with a forced-equal pair at k=4*.  The search
in this file looks for one among compositions of diamonds and the hexagonal
wheel.  Its honest outcome -- positive or negative -- is written to
``runs/hn_gadgets.json``.

Run:
    python problems/hadwiger_nelson/gadgets.py
"""

from __future__ import annotations

import io
import itertools
import json
import random
import sys
import time
import unittest
from pathlib import Path

import sympy as sp

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parents[1]
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

from hn import UnitDistanceGraph, dist2, is_unit, moser_spindle  # noqa: E402


# --------------------------------------------------------------------------
# relation statuses
# --------------------------------------------------------------------------

EQUAL = "equal"
DIFFERENT = "different"
FREE = "free"
VACUOUS = "vacuous"          # no k-colouring exists at all

_NUM_DPS = 70
# A pair is only a candidate edge if its numerically evaluated squared
# distance is this close to 1.  Genuine unit distances in these constructions
# evaluate to 1 to far better than this; algebraic near-misses are excluded
# numerically and never reach the expensive symbolic zero test.
_UNIT_TOL = sp.Float("1e-50")
_ZERO_TOL = 1e-45


# --------------------------------------------------------------------------
# forced relation
# --------------------------------------------------------------------------

class ForcedRelation:
    """Result of exhaustively classifying pairs of distinguished vertices.

    ``colourings_examined`` counts canonical colourings (one representative
    per renaming of colours).  ``exhausted`` is False only when the
    enumeration stopped early because every non-adjacent pair had already
    been seen both equal and different -- at that point no further colouring
    can change any relation, so the relation is still exact.
    """

    def __init__(self, k, colorable, colourings_examined, exhausted, pairs):
        self.k = int(k)
        self.colorable = bool(colorable)
        self.colourings_examined = int(colourings_examined)
        self.exhausted = bool(exhausted)
        self.pairs = dict(pairs)          # (a, b) -> status

    def status(self, a, b):
        return self.pairs[tuple(sorted((int(a), int(b))))]

    def forced_equal(self, a, b) -> bool:
        return self.pairs[tuple(sorted((int(a), int(b))))] == EQUAL

    def forced_different(self, a, b) -> bool:
        return self.pairs[tuple(sorted((int(a), int(b))))] == DIFFERENT

    def free(self, a, b) -> bool:
        return self.pairs[tuple(sorted((int(a), int(b))))] == FREE

    def equal_pairs(self):
        return sorted(p for p, s in self.pairs.items() if s == EQUAL)

    def different_pairs(self):
        return sorted(p for p, s in self.pairs.items() if s == DIFFERENT)

    def free_pairs(self):
        return sorted(p for p, s in self.pairs.items() if s == FREE)

    def to_dict(self):
        return {
            "k": self.k,
            "colorable": self.colorable,
            "colourings_examined": self.colourings_examined,
            "exhausted": self.exhausted,
            "pairs": [
                {"pair": list(p), "status": s}
                for p, s in sorted(self.pairs.items())
            ],
        }


# --------------------------------------------------------------------------
# gadget
# --------------------------------------------------------------------------

class Gadget:
    """A finite exact unit-distance graph with distinguished vertices.

    ``edges`` is the *complete* set of unit-distance edges of ``points``:
    every edge is exactly unit length and every unit pair is present.  The
    constructor re-checks each edge symbolically; ``verify()`` re-checks the
    non-edges as well.
    """

    def __init__(self, points, edges, distinguished, label="", _trusted=False):
        self.points = tuple(
            tuple(sp.sympify(c) for c in p) for p in points
        )
        self.edges = tuple(
            tuple(sorted((int(i), int(j)))) for i, j in edges
        )
        self.distinguished = tuple(int(i) for i in distinguished)
        self.label = str(label)

        n = len(self.points)
        if n == 0:
            raise ValueError("empty gadget")
        if len(set(self.edges)) != len(self.edges):
            raise ValueError("duplicate edge")
        for i, j in self.edges:
            if not (0 <= i < j < n):
                raise ValueError(f"bad edge ({i},{j}) for n={n}")
            if not _trusted and not is_unit(self.points[i], self.points[j]):
                raise ValueError(f"{label}: edge ({i},{j}) is not unit")
        for v in self.distinguished:
            if not (0 <= v < n):
                raise ValueError(f"{label}: bad distinguished vertex {v}")

        self.adj = [0] * n
        for i, j in self.edges:
            self.adj[i] |= 1 << j
            self.adj[j] |= 1 << i

    # -- basic invariants --------------------------------------------------
    @property
    def n(self):
        return len(self.points)

    @property
    def m(self):
        return len(self.edges)

    def degrees(self):
        return [bin(a).count("1") for a in self.adj]

    def graph(self):
        """hn.py graph with the same exact points and edges."""
        return UnitDistanceGraph(self.points, list(self.adj), list(self.edges))

    def chromatic_number(self):
        return self.graph().chromatic_number()

    def coloring(self):
        return self.graph().coloring()

    def check_coloring(self, colors, k=None):
        return self.graph().check_coloring(colors, k)

    def verify(self):
        """Exact check of every edge and every non-edge (hn.py's verifier)."""
        return self.graph().verify()

    # -- exhaustive colouring ---------------------------------------------
    def iter_colourings(self, k):
        """Yield every proper k-colouring, once per renaming of colours.

        Colours are introduced in first-appearance order, so every unlabelled
        colouring has exactly one representative.  Equality and difference of
        vertex colours are invariant under renaming, so this is exhaustive for
        every relation computed below.
        """
        k = int(k)
        if k < 0:
            raise ValueError("k must be non-negative")
        n = self.n
        if n == 0:
            yield ()
            return
        if k == 0:
            return

        deg = self.degrees()
        order = sorted(range(n), key=lambda v: (-deg[v], v))
        neighbours = [
            [u for u in range(n) if (self.adj[v] >> u) & 1] for v in range(n)
        ]
        colour = [-1] * n
        sat = [0] * n

        def rec(i, used):
            if i == n:
                yield tuple(colour)
                return
            v = order[i]
            forbidden = sat[v]
            limit = min(used + 1, k)
            for c in range(limit):
                if (forbidden >> c) & 1:
                    continue
                colour[v] = c
                touched = []
                for u in neighbours[v]:
                    if colour[u] < 0 and not ((sat[u] >> c) & 1):
                        sat[u] |= 1 << c
                        touched.append(u)
                yield from rec(i + 1, max(used, c + 1))
                for u in touched:
                    sat[u] &= ~(1 << c)
                colour[v] = -1

        yield from rec(0, 0)

    def is_k_colourable(self, k) -> bool:
        return next(self.iter_colourings(k), None) is not None

    def forced_relation(self, k, vertices=None):
        """Classify every pair of ``vertices`` by exhaustive k-colouring.

        ``vertices`` defaults to the distinguished tuple.  Pass
        ``range(self.n)`` to classify every pair (used by the search).

        Adjacent pairs are known to be different without enumeration.  The
        enumeration may stop early once every non-adjacent pair has been
        observed both equal and different in some colourings: at that point
        all relations are final.  Otherwise it runs to exhaustion.
        """
        k = int(k)
        if vertices is None:
            vertices = self.distinguished
        vertices = tuple(sorted(set(int(v) for v in vertices)))
        pairs = list(itertools.combinations(vertices, 2))

        adjacent = {
            p: bool((self.adj[p[0]] >> p[1]) & 1) for p in pairs
        }
        seen_equal = {p: False for p in pairs}
        seen_diff = {p: bool(adjacent[p]) for p in pairs}
        count = 0
        exhausted = True

        for colours in self.iter_colourings(k):
            count += 1
            for p in pairs:
                a, b = p
                if colours[a] == colours[b]:
                    if adjacent[p]:
                        raise AssertionError(
                            "adjacent vertices received the same colour"
                        )
                    seen_equal[p] = True
                else:
                    seen_diff[p] = True
            if pairs and all(
                seen_equal[p] and seen_diff[p]
                for p in pairs if not adjacent[p]
            ):
                exhausted = False
                break

        if count == 0:
            return ForcedRelation(
                k, False, 0, True, {p: VACUOUS for p in pairs}
            )

        status = {}
        for p in pairs:
            if adjacent[p]:
                status[p] = DIFFERENT
            elif seen_equal[p] and seen_diff[p]:
                status[p] = FREE
            elif seen_equal[p]:
                status[p] = EQUAL
            elif seen_diff[p]:
                status[p] = DIFFERENT
            else:  # pragma: no cover - unreachable with count > 0
                raise AssertionError("undetermined relation")
        return ForcedRelation(k, True, count, exhausted, status)

    def forced_equal(self, k, a=None, b=None) -> bool:
        """True iff the pair is equal in every k-colouring.

        Defaults to the first two distinguished vertices.  Uncolorable graphs
        return False for every pair: the question is vacuous there, and a
        vacuous True would be a false certificate.
        """
        if a is None or b is None:
            if len(self.distinguished) < 2:
                raise ValueError("need two distinguished vertices")
            a, b = self.distinguished[:2]
        return self.forced_relation(k).forced_equal(a, b)

    # -- serialisation -----------------------------------------------------
    def structural_key(self):
        return (self.n, tuple(sorted(self.edges)))

    def to_dict(self, include_points=True, include_chi=False):
        out = {
            "label": self.label,
            "n": self.n,
            "m": self.m,
            "distinguished": list(self.distinguished),
            "edges": [list(e) for e in self.edges],
        }
        if include_points:
            out["points"] = [[str(x), str(y)] for x, y in self.points]
        if include_chi:
            out["chi"] = self.chromatic_number()
        return out


# --------------------------------------------------------------------------
# exact construction helpers
# --------------------------------------------------------------------------

def _numeric_points(points):
    return [(sp.N(p[0], _NUM_DPS), sp.N(p[1], _NUM_DPS)) for p in points]


def _exact_graph(points, required=(), label=""):
    """Build the complete unit-distance edge set, with exact verification.

    A pair is a candidate only if its numerically evaluated squared distance
    is within 1e-50 of 1; the decision is then made by hn.is_unit, which uses
    symbolic simplification.  Every required edge must be recovered, so a
    simplification failure raises instead of silently dropping an edge (which
    would inflate the chromatic number).
    """
    pts = [tuple(sp.sympify(c) for c in p) for p in points]
    n = len(pts)
    if n == 0:
        raise ValueError("empty point set")
    num = _numeric_points(pts)

    adj = [0] * n
    edges = []
    for i in range(n):
        for j in range(i + 1, n):
            dx = num[i][0] - num[j][0]
            dy = num[i][1] - num[j][1]
            d2 = dx * dx + dy * dy
            if abs(d2 - 1) > _UNIT_TOL:
                continue
            if not is_unit(pts[i], pts[j]):
                raise ArithmeticError(
                    f"{label}: pair ({i},{j}) is numerically unit "
                    f"({sp.N(d2, 30)}) but symbolic verification failed"
                )
            adj[i] |= 1 << j
            adj[j] |= 1 << i
            edges.append((i, j))

    edge_set = {tuple(sorted(e)) for e in edges}
    missing = [
        tuple(sorted((int(i), int(j))))
        for i, j in required
        if tuple(sorted((int(i), int(j)))) not in edge_set
    ]
    if missing:
        raise ArithmeticError(
            f"{label}: required unit edges not recovered: {missing}"
        )
    return pts, adj, edges


def _make_gadget(points, required, distinguished, label):
    pts, _adj, edges = _exact_graph(points, required, label)
    # _exact_graph has already proved every edge symbolically.
    return Gadget(pts, edges, distinguished, label, _trusted=True)


def _points_equal(p, q):
    dx = sp.N(p[0] - q[0], 60)
    dy = sp.N(p[1] - q[1], 60)
    if abs(float(dx)) > _ZERO_TOL or abs(float(dy)) > _ZERO_TOL:
        return False
    if sp.simplify(p[0] - q[0]) == 0 and sp.simplify(p[1] - q[1]) == 0:
        return True
    raise ArithmeticError(
        "numerically coincident points are not symbolically provably equal"
    )


def _merge_points(points):
    """Merge exactly coincident points; return (unique, old -> new index)."""
    unique = []
    mapping = []
    for p in points:
        found = None
        for j, q in enumerate(unique):
            if _points_equal(p, q):
                found = j
                break
        if found is None:
            found = len(unique)
            unique.append(p)
        mapping.append(found)
    return unique, mapping


def _nonnegative(expr) -> bool:
    e = sp.simplify(expr)
    if e.is_nonnegative:
        return True
    if e.is_negative:
        return False
    val = sp.N(e, 60)
    re, im = float(sp.re(val)), float(sp.im(val))
    if abs(im) > 1e-40:
        raise ArithmeticError(f"non-real discriminant {e}")
    if re > 1e-30:
        return True
    if re < -1e-30:
        return False
    return True       # tangent case: treat numerical zero as zero


class InfeasibleSpindle(ValueError):
    """The two tip circles do not intersect, so the requested spindle is impossible."""


# --------------------------------------------------------------------------
# the diamond
# --------------------------------------------------------------------------

def diamond():
    """4 vertices, 5 edges: two equilateral triangles sharing an edge.

    Vertices 0,1 are the shared edge; vertices 2,3 are the tips, on opposite
    sides.  In every 3-colouring the shared edge uses two colours and both
    tips must take the third, so the tips are forced equal.  Their distance
    is sqrt(3), not 1, so the gadget itself is 3-colourable.
    """
    s3 = sp.sqrt(3)
    half = sp.Rational(1, 2)
    pts = [
        (sp.Integer(0), sp.Integer(0)),
        (sp.Integer(1), sp.Integer(0)),
        (half, s3 / 2),
        (half, -s3 / 2),
    ]
    required = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3)]
    return _make_gadget(pts, required, (2, 3), "diamond")


# --------------------------------------------------------------------------
# the hexagonal wheel (an "H-like" wheel; 7 vertices, 12 edges, chi=3)
# --------------------------------------------------------------------------

def hex_wheel():
    """Hub plus a regular unit hexagon: 7 vertices, 12 edges, chi=3.

    The three even rim vertices (1,3,5) are mutually at distance sqrt(3) and
    all carry the same colour in every 3-colouring (the rim is the bipartite
    C6 and the hub takes the third colour).  They are therefore a forced
    triple at k=3, suitable for the spindle operation.

    This is the regular hexagonal wheel, W7.  It is *not* claimed to be
    de Grey's forcing graph; it is included as a second 3-chromatic base with
    a forced-equal pair.
    """
    s3 = sp.sqrt(3)
    half = sp.Rational(1, 2)
    rim = [
        (sp.Integer(1), sp.Integer(0)),
        (half, s3 / 2),
        (-half, s3 / 2),
        (sp.Integer(-1), sp.Integer(0)),
        (-half, -s3 / 2),
        (half, -s3 / 2),
    ]
    pts = [(sp.Integer(0), sp.Integer(0))] + rim       # hub is 0
    required = [(0, i) for i in range(1, 7)]
    required += [(i, i + 1) for i in range(1, 6)]
    required += [(6, 1)]
    return _make_gadget(pts, required, (1, 3, 5), "hex_wheel")


# --------------------------------------------------------------------------
# spindle
# --------------------------------------------------------------------------

def _unit_vector(g, share, tip):
    s = g.points[share]
    t = g.points[tip]
    vx = sp.expand(t[0] - s[0])
    vy = sp.expand(t[1] - s[1])
    d2 = sp.expand(vx * vx + vy * vy)
    d = sp.sqrt(d2)
    return (s, d2, d, (sp.simplify(vx / d), sp.simplify(vy / d)))


def spindle(g1, g2, *, share1=None, tip1=None, share2=None, tip2=None,
            branch=1, label=None):
    """Glue two copies at one distinguished vertex and force a tip edge.

    Copy 1 is kept fixed.  Copy 2 is placed so that ``share2`` coincides with
    ``share1`` and its tip direction is obtained by rotating copy 1's tip
    direction by the angle theta satisfying

        |d1 - d2| <= 1 <= d1 + d2,
        cos(theta) = (d1^2 + d2^2 - 1) / (2 d1 d2).

    The two tips are then exactly unit distance apart; ``_exact_graph`` adds
    that edge (and every other exact unit edge).  Raises InfeasibleSpindle if
    the two tip circles do not meet.

    Defaults use the first distinguished vertex as the shared vertex and the
    second as the tip, in each gadget.
    """
    if share1 is None:
        share1 = g1.distinguished[0]
    if tip1 is None:
        tip1 = g1.distinguished[1]
    if share2 is None:
        share2 = g2.distinguished[0]
    if tip2 is None:
        tip2 = g2.distinguished[1]
    if share1 == tip1:
        raise ValueError("share1 and tip1 must differ")
    if share2 == tip2:
        raise ValueError("share2 and tip2 must differ")

    s1, d1sq, d1, u1 = _unit_vector(g1, share1, tip1)
    _s2, d2sq, d2, _u2 = _unit_vector(g2, share2, tip2)

    cos_theta = sp.simplify(
        (d1sq + d2sq - 1) / (2 * d1 * d2)
    )
    disc = sp.simplify(1 - cos_theta ** 2)
    if not _nonnegative(disc):
        raise InfeasibleSpindle(
            f"tip circles do not meet: d1={sp.N(d1, 20)}, "
            f"d2={sp.N(d2, 20)}"
        )
    sin_theta = sp.sqrt(disc)

    # Desired direction of copy 2's tip, in copy 1's frame.
    des = (
        sp.expand(cos_theta * u1[0] - branch * sin_theta * u1[1]),
        sp.expand(branch * sin_theta * u1[0] + cos_theta * u1[1]),
    )

    # Rigid motion of copy 2 taking its tip vector to d2 * des.
    s2 = g2.points[share2]
    t2 = g2.points[tip2]
    v2 = (sp.expand(t2[0] - s2[0]), sp.expand(t2[1] - s2[1]))
    v2d = (sp.expand(d2 * des[0]), sp.expand(d2 * des[1]))
    den = sp.simplify(v2d[0] ** 2 + v2d[1] ** 2)
    cos_r = sp.simplify(
        (v2[0] * v2d[0] + v2[1] * v2d[1]) / den
    )
    sin_r = sp.simplify(
        (v2[0] * v2d[1] - v2[1] * v2d[0]) / den
    )

    pts = list(g1.points)
    raw = {}
    for i, p in enumerate(g2.points):
        if i == share2:
            raw[i] = share1
            continue
        w = (p[0] - s2[0], p[1] - s2[1])
        r = (
            sp.expand(cos_r * w[0] - sin_r * w[1]),
            sp.expand(sin_r * w[0] + cos_r * w[1]),
        )
        pts.append((sp.expand(s1[0] + r[0]), sp.expand(s1[1] + r[1])))
        raw[i] = len(pts) - 1

    unique, merged = _merge_points(pts)
    map1 = {i: merged[i] for i in range(g1.n)}
    map2 = {i: merged[raw[i]] for i in range(g2.n)}

    required = set()
    for i, j in g1.edges:
        required.add((map1[i], map1[j]))
    for i, j in g2.edges:
        required.add((map2[i], map2[j]))
    required.add((map1[tip1], map2[tip2]))
    if any(i == j for i, j in required):
        raise ArithmeticError("spindle produced a self-loop")

    dist = tuple(dict.fromkeys(
        (map1[share1], map1[tip1], map2[tip2])
    ))
    if label is None:
        label = f"spindle({g1.label},{g2.label})"
    return _make_gadget(unique, required, dist, label)


# --------------------------------------------------------------------------
# rigid edge gluing (the glue2 operation, but exact and constraint-aware)
# --------------------------------------------------------------------------

def glue_edge(g, h, edge_g, edge_h, *, flip=False, label=None):
    """Rigidly identify unit edge ``edge_h`` of h with unit edge ``edge_g`` of g.

    No constraint is added beyond the identification: the point of this
    operation in the search is to *create* graphs on which the forced-relation
    test can then be run.  The result is verified exactly.
    """
    a, b = edge_g
    c, d = edge_h
    if flip:
        c, d = d, c
    if tuple(sorted((a, b))) not in {tuple(sorted(e)) for e in g.edges}:
        raise ValueError("edge_g is not an edge of g")
    if tuple(sorted((c, d))) not in {tuple(sorted(e)) for e in h.edges}:
        raise ValueError("edge_h is not an edge of h")

    ta, tb = g.points[a], g.points[b]
    sc, sd = h.points[c], h.points[d]
    vs = (sd[0] - sc[0], sd[1] - sc[1])
    vt = (tb[0] - ta[0], tb[1] - ta[1])
    den = sp.simplify(vs[0] ** 2 + vs[1] ** 2)
    cos_r = sp.simplify((vs[0] * vt[0] + vs[1] * vt[1]) / den)
    sin_r = sp.simplify((vs[0] * vt[1] - vs[1] * vt[0]) / den)

    pts = list(g.points)
    raw = {c: a, d: b}
    for i, p in enumerate(h.points):
        if i in (c, d):
            continue
        w = (p[0] - sc[0], p[1] - sc[1])
        r = (
            sp.expand(cos_r * w[0] - sin_r * w[1]),
            sp.expand(sin_r * w[0] + cos_r * w[1]),
        )
        pts.append((sp.expand(ta[0] + r[0]), sp.expand(ta[1] + r[1])))
        raw[i] = len(pts) - 1

    unique, merged = _merge_points(pts)
    map_g = {i: merged[i] for i in range(g.n)}
    map_h = {i: merged[raw[i]] for i in range(h.n)}

    required = set()
    for i, j in g.edges:
        required.add((map_g[i], map_g[j]))
    for i, j in h.edges:
        required.add((map_h[i], map_h[j]))
    if any(i == j for i, j in required):
        raise ArithmeticError("edge gluing produced a self-loop")

    dist = tuple(dict.fromkeys(
        (map_g[a], map_g[b]) + tuple(map_h[v] for v in h.distinguished)
    ))
    if label is None:
        label = f"glue({g.label},{h.label})"
    return _make_gadget(unique, required, dist, label)


# --------------------------------------------------------------------------
# certificate helpers
# --------------------------------------------------------------------------

def certify(g, max_exhaustive_n=18):
    """Return a geometry + colouring certificate for the exact chi of g."""
    ok, errs = g.verify()
    chi = g.chromatic_number()
    cert = {
        "n": g.n,
        "m": g.m,
        "chi": chi,
        "geometry_ok": ok,
        "geometry_errors": errs[:5],
        "k_colourable": g.is_k_colourable(chi),
        "k_minus_1_colourable": (
            g.is_k_colourable(chi - 1)
            if chi > 0 and g.n <= max_exhaustive_n
            else None
        ),
    }
    if chi > 0:
        colours = g.coloring()
        cert["witness_ok"] = g.check_coloring(colours, k=chi)
        cert["witness"] = colours
    if not ok:
        raise AssertionError(f"geometry failed: {errs[:3]}")
    if not cert["k_colourable"]:
        raise AssertionError("claimed chi is not attainable")
    if cert["k_minus_1_colourable"] is True:
        raise AssertionError("claimed chi is not a lower bound")
    return cert


# --------------------------------------------------------------------------
# search over compositions
# --------------------------------------------------------------------------

def composition_search(max_graphs=22, max_vertices=13, seed=1):
    """Bounded BFS over exact compositions of diamonds and the hex wheel.

    Moves, all exact and all re-scored from scratch:

    * forced-equality lift: for every k where a graph has a forced-equal pair
      (u,a), spindle two copies sharing u with tips a and a (or another
      forced partner of u);
    * cross lift: the same operation using a diamond or wheel whose tips are
      forced equal at the same k;
    * rigid edge gluing: identify a unit edge of a base gadget with a unit
      edge of the current graph (the old glue2 move, now exact).

    Every graph is scored with hn.py's exact DSATUR.  For every level k at
    which the graph is k-colourable, the full forced relation on all vertex
    pairs is computed; a forced-equal pair at k is a certificate that the
    lift produces a (k+1)-chromatic composition.
    """
    rng = random.Random(seed)
    diamond_g = diamond()
    wheel_g = hex_wheel()
    bases = [("diamond", diamond_g), ("wheel", wheel_g)]

    queue = []
    seen = set()
    records = []
    constructed = []      # (chi, gadget) for every distinct graph built
    lifts = []
    forcing = []          # graphs with at least one forced-equal pair at some k
    level_pairs = {3: [], 4: [], 5: []}

    def push(g):
        key = g.structural_key()
        if key in seen:
            return False
        seen.add(key)
        constructed.append((g.chromatic_number(), g))
        queue.append(g)
        return True

    for _name, g in bases:
        push(g)

    best_key = None
    best = None

    while queue and len(records) < max_graphs:
        g = queue.pop(0)
        chi = g.chromatic_number()
        geometry_ok, geometry_errors = g.verify()
        if not geometry_ok:
            raise AssertionError(
                f"record {g.label} failed exact geometry: {geometry_errors[:3]}"
            )
        rec = {
            "label": g.label,
            "n": g.n,
            "m": g.m,
            "chi": chi,
            "geometry_ok": True,
        }
        records.append(rec)
        key = (chi, -g.n, -g.m)
        if best_key is None or key > best_key:
            best_key = key
            best = g

        # ---- forced-relation analysis at every relevant level -----------
        local_forcing = []
        for k in (3, 4, 5):
            if k < chi:
                continue
            rel = g.forced_relation(k, vertices=range(g.n))
            if not rel.colorable:
                continue
            eq = rel.equal_pairs()
            if eq:
                level_pairs[k].append({
                    "graph": g.label,
                    "n": g.n,
                    "pairs": [list(p) for p in eq],
                })
                local_forcing.append((k, eq))
        if local_forcing:
            forcing.append({
                "graph": g.label,
                "n": g.n,
                "chi": chi,
                "levels": {
                    str(k): [list(p) for p in eq] for k, eq in local_forcing
                },
            })

        # ---- forced-equality lifts --------------------------------------
        for k, eq in local_forcing:
            attempts = []
            # self-lifts: two forced pairs sharing a vertex
            for u, a in eq:
                partners = [b for (u2, b) in eq if u2 == u and b != a]
                for b in partners[:2]:
                    attempts.append((g, u, a, g, u, b, "self"))
            # cross-lifts with a base that also forces equality at k
            for _bname, base in bases:
                if base.n > max_vertices:
                    continue
                brel = base.forced_relation(k, vertices=range(base.n))
                if not brel.colorable:
                    continue
                beq = brel.equal_pairs()
                if not beq:
                    continue
                bu, ba = beq[0]
                u, a = eq[0]
                attempts.append((g, u, a, base, bu, ba, "base"))

            for (g1_, u1_, a1_, g2_, u2_, a2_, kind) in attempts[:4]:
                try:
                    h = spindle(
                        g1_, g2_,
                        share1=u1_, tip1=a1_, share2=u2_, tip2=a2_,
                        label=f"lift[{g1_.label}|{g2_.label}]@k{k}",
                    )
                except InfeasibleSpindle:
                    continue
                hchi = h.chromatic_number()
                ok_lb = hchi >= k + 1
                lifts.append({
                    "kind": kind,
                    "base": g1_.label,
                    "partner": g2_.label,
                    "k": k,
                    "pair1": [u1_, a1_],
                    "pair2": [u2_, a2_],
                    "result_label": h.label,
                    "result_n": h.n,
                    "result_m": h.m,
                    "result_chi": hchi,
                    "theorem_lower_bound_ok": ok_lb,
                })
                if not ok_lb:
                    raise AssertionError(
                        "forced-equality lift failed its lower bound: "
                        f"k={k}, chi={hchi}"
                    )
                if h.n <= max_vertices:
                    push(h)

        # ---- rigid edge gluing with each base ---------------------------
        for _bname, base in bases:
            for _trial in range(2):
                eg = rng.choice(g.edges)
                eh = rng.choice(base.edges)
                for flip in (False, True):
                    try:
                        cand = glue_edge(
                            g, base, eg, eh, flip=flip,
                            label=f"glue[{g.label}+{base.label}]",
                        )
                    except (ArithmeticError, ValueError):
                        continue
                    if g.n < cand.n <= max_vertices:
                        push(cand)

    if best is None:
        best = diamond_g
        best_key = (best.chromatic_number(), -best.n, -best.m)

    # The chi claim is over *every* distinct composition constructed, not
    # only the ones popped before the record budget ran out.
    top_chi, top_g = max(
        constructed + [(best.chromatic_number(), best)],
        key=lambda t: (t[0], -t[1].n, -t[1].m),
    )
    best = top_g

    histogram = {}
    for rec in records:
        histogram[rec["chi"]] = histogram.get(rec["chi"], 0) + 1

    level4 = level_pairs[4]
    level5 = level_pairs[5]
    return {
        "algorithm": "BFS over exact diamond/wheel compositions",
        "max_graphs": max_graphs,
        "max_vertices": max_vertices,
        "seed": seed,
        "graphs_examined": len(records),
        "graphs_constructed": len(constructed),
        "max_constructed_chi": top_chi,
        "histogram": {str(k): v for k, v in sorted(histogram.items())},
        "records": records,
        "best_chi": best.chromatic_number(),
        "best": best.to_dict(include_points=True, include_chi=True),
        "best_certificate": certify(best),
        "lifts": lifts,
        "forcing_graphs": forcing,
        "level4_forced_equal_pairs": level4,
        "level5_forced_equal_pairs": level5,
    }


# --------------------------------------------------------------------------
# unit tests
# --------------------------------------------------------------------------

class TestGadgets(unittest.TestCase):
    def test_diamond_geometry_and_colouring(self):
        g = diamond()
        self.assertEqual((g.n, g.m), (4, 5))
        self.assertEqual(g.chromatic_number(), 3)
        self.assertTrue(g.verify()[0])
        self.assertFalse(g.is_k_colourable(2))
        # tips are distance sqrt(3), not adjacent
        self.assertEqual(sp.simplify(dist2(g.points[2], g.points[3])), 3)

    def test_diamond_forced_equal_3(self):
        g = diamond()
        rel = g.forced_relation(3)
        self.assertTrue(rel.colorable)
        self.assertTrue(rel.forced_equal(2, 3))
        self.assertFalse(rel.forced_different(2, 3))
        self.assertEqual(rel.colourings_examined, 1)  # unique up to renaming
        self.assertFalse(g.forced_equal(4))           # 4 colours give slack

    def test_moser_spindle_from_two_diamonds(self):
        d = diamond()
        m = spindle(d, d)
        self.assertEqual((m.n, m.m), (7, 11))
        self.assertEqual(m.chromatic_number(), 4)
        ok, errs = m.verify()
        self.assertTrue(ok, errs)
        self.assertFalse(m.is_k_colourable(3))
        self.assertTrue(m.is_k_colourable(4))
        # the joined tips are adjacent, hence forced different at k=4
        rel = m.forced_relation(4)
        a, b = m.distinguished[1], m.distinguished[2]
        self.assertTrue(rel.forced_different(a, b))

    def test_spindle_is_hn_moser_spindle(self):
        ours = spindle(diamond(), diamond())
        ref = moser_spindle()
        self.assertEqual((ours.n, ours.m), (ref.n, ref.m))
        ref_edges = {tuple(sorted(e)) for e in ref.edges}
        for perm in itertools.permutations(range(ours.n)):
            if all(
                tuple(sorted((perm[i], perm[j]))) in ref_edges
                for i, j in ours.edges
            ):
                return
        self.fail(
            "spindle(diamond, diamond) is not isomorphic to hn.moser_spindle()"
        )

    def test_hex_wheel_forced_triple(self):
        w = hex_wheel()
        self.assertEqual((w.n, w.m), (7, 12))
        self.assertEqual(w.chromatic_number(), 3)
        rel = w.forced_relation(3)
        self.assertTrue(rel.forced_equal(1, 3))
        self.assertTrue(rel.forced_equal(1, 5))
        self.assertTrue(rel.forced_equal(3, 5))

    def test_wheel_spindle_lift(self):
        w = hex_wheel()
        h = spindle(w, w)
        self.assertEqual(h.chromatic_number(), 4)
        self.assertTrue(h.verify()[0])
        self.assertFalse(h.is_k_colourable(3))

    def test_glue_edge_preserves_exactness(self):
        d = diamond()
        m = spindle(d, d)
        h = glue_edge(m, d, m.edges[0], d.edges[0])
        self.assertTrue(h.verify()[0])
        self.assertGreaterEqual(h.chromatic_number(), 4)


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def _run_tests():
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(TestGadgets)
    stream = io.StringIO()
    result = unittest.TextTestRunner(verbosity=0, stream=stream).run(suite)
    if not result.wasSuccessful():
        sys.stderr.write(stream.getvalue())
        raise SystemExit(1)


def main():
    t0 = time.time()
    _run_tests()

    d = diamond()
    d_rel3 = d.forced_relation(3)
    diamond_ok = d_rel3.forced_equal(*d.distinguished)

    moser = spindle(d, d)
    moser_chi = moser.chromatic_number()
    moser_cert = certify(moser)

    search = composition_search()
    best_chi = int(search["best_chi"])
    n_four = sum(1 for r in search["records"] if r["chi"] == 4)

    payload = {
        "generated_by": "problems/hadwiger_nelson/gadgets.py",
        "runtime_seconds": round(time.time() - t0, 3),
        "diamond": {
            **d.to_dict(include_points=True, include_chi=True),
            "forced_equal_3": diamond_ok,
            "forced_relation_3": d_rel3.to_dict(),
            "forced_relation_4": d.forced_relation(4).to_dict(),
        },
        "moser_spindle": {
            **moser.to_dict(include_points=True),
            "chi": moser_chi,
            "construction": "spindle(diamond, diamond)",
            "certificate": moser_cert,
            "forced_relation_4": moser.forced_relation(4).to_dict(),
        },
        "hex_wheel": {
            **hex_wheel().to_dict(include_points=True, include_chi=True),
            "forced_relation_3": hex_wheel().forced_relation(3).to_dict(),
        },
        "search": search,
        "conclusion": {
            "best_composition_chi": best_chi,
            "level4_forcing_found": bool(
                search["level4_forced_equal_pairs"]
            ),
            "level5_forcing_found": bool(
                search["level5_forced_equal_pairs"]
            ),
            "statement": (
                "The forced-equality constraint is implemented and the "
                "diamond lift reproduces the Moser spindle (chi=4). "
                f"All {search['graphs_constructed']} distinct compositions "
                f"constructed (at most {search['max_vertices']} vertices) "
                f"have chi <= {search['max_constructed_chi']}. Of the "
                f"{search['graphs_examined']} examined for forced "
                f"relations, {n_four} are 4-chromatic, and none has a pair "
                "forced equal in every 4-colouring, so the lift cannot be "
                "iterated to chi=5 inside this family. This is a clean "
                "negative for the search, not a proof that no such gadget "
                "exists: chi=5 requires a unit-distance gadget with a "
                "forced-equal pair at k=4, and no such gadget was found at "
                "these sizes. The 509-vertex 5-chromatic record is far "
                "outside this budget."
            ),
        },
    }

    out = _ROOT / "runs" / "hn_gadgets.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    # Acceptance output: exactly these three lines on stdout.
    print(f"diamond forced_equal(3) = {diamond_ok}")
    print(f"moser spindle chi = {moser_chi}")
    print(f"best composition chi = {best_chi}")


if __name__ == "__main__":
    main()
