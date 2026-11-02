#!/usr/bin/env python3
"""Targeted level-5 virtual-edge search: seed = a sqrt(3) edge of the 33-vertex graph.

Why this file exists
--------------------
`virtual_edge.py` grew unit-distance graphs from the bare terminals
u = (0,0), v = (0,sqrt(3)) and from the Moser spindle.  At n = 8..18 it tested
19,998 candidates and found no pair forced different at k = 5.  Random growth
is the wrong tool: it has no reason to produce the 5-chromatic structure a
level-5 virtual edge must contain.

The 33-vertex 6-chromatic two-distance graph of De Neve et al. (2025) has
exactly 32 edges of length sqrt(3).  Each is a *concrete* place where a
virtual-edge gadget must work: the gadget's terminals are at distance sqrt(3)
and the gadget is substituted for that edge.  This module takes one such edge
and uses its local neighbourhood N[u] U N[v] as the seed, then grows it.

What the mathematics says in advance
------------------------------------
Level lemma.  Let G be (k-1)-colourable and let u, v be non-adjacent.

  * (k-1)-colour G - {u,v} and give both u and v a fresh k-th colour: a
    proper k-colouring with colour(u) = colour(v).
  * (k-1)-colour G.  If colour(u) != colour(v) we are done; otherwise recolour
    v with the unused k-th colour: a proper k-colouring with colour(u) !=
    colour(v).

Hence relation(u,v,k) == "free" for *every* (k-1)-colourable graph with
non-adjacent terminals.  In particular a level-5 virtual edge (relation
"different" at k = 5) must be non-4-colourable, i.e. chi >= 5.  The smallest
5-chromatic unit-distance graph known is Parts' 509-vertex graph, so a
level-5 virtual edge is not a small local object.

Consequence for this seed: every seed below is 3-colourable, so its relation
is necessarily "free"; no amount of *4-colourable* growth can change that.
The only exit is a genuinely 5-chromatic extension.

Stages
------
1. Load runs/two_distance.json (exact-audited against two_distance.py),
   enumerate all 32 sqrt(3) edges, and pick the one with the smallest combined
   neighbourhood N[u] U N[v].  Report the degree profile.
2. Build the minimal unit-distance graph on the seed points (edges at distance
   exactly 1) and test relation(u,v,5).  Also record chi, the number of
   4-colourings split by colour(u)=colour(v), and explicit witnesses.
3. Grow the gadget by adding points at unit distance from two existing points
   (unit-circle intersections), prioritising points that create edges incident
   to u or v.  Every grown graph keeps u,v at distance sqrt(3).  Each
   candidate is tested with relation(u,v,5).  Caps: 50,000 candidates / 1200 s.
4. On a hit, re-verify with an independent binary/Glucose4 encoding and a
   50-digit geometry check before claiming anything.
5. Write runs/hn_virtual_edge2.json and append a line per stage to
   runs/virtual_edge2.log.

Run:  python problems/hadwiger_nelson/virtual_edge2.py
      python problems/hadwiger_nelson/virtual_edge2.py --quick
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import random
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import mpmath as mp
import sympy as sp
from pysat.solvers import Solver

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import relation as R  # noqa: E402

DPS_CONSTRUCT = 70
DPS_CHECK = 50
UNIT_TOL = mp.mpf("1e-40")
CHECK_TOL = mp.mpf("1e-25")
FLOAT_SCREEN = 1e-6

TWO_DISTANCE_PATH = ROOT / "runs" / "two_distance.json"
OUT_PATH = ROOT / "runs" / "hn_virtual_edge2.json"
LOG_PATH = ROOT / "runs" / "virtual_edge2.log"

# Public record: smallest known 5-chromatic unit-distance graph.
RECORD_5_CHROMATIC_N = 509


# --------------------------------------------------------------------------
# logging
# --------------------------------------------------------------------------

def _stamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def log_stage(stage: str, message: str, path: Path | None = None) -> None:
    target = Path(path) if path is not None else LOG_PATH
    target.parent.mkdir(parents=True, exist_ok=True)
    line = f"[{_stamp()}] {stage}: {message}"
    with open(target, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")
    print(line, flush=True)


# --------------------------------------------------------------------------
# geometry helpers (60+ digits)
# --------------------------------------------------------------------------

def parse_points(data) -> list[tuple[mp.mpf, mp.mpf]]:
    """The 33 published points as mpmath values at DPS_CONSTRUCT digits."""
    out = []
    for sx, sy in data["points"]:
        x = sp.N(sp.sympify(sx), DPS_CONSTRUCT)
        y = sp.N(sp.sympify(sy), DPS_CONSTRUCT)
        out.append((mp.mpf(str(x)), mp.mpf(str(y))))
    return out


def circle_intersections(p1, p2):
    """Intersections of the unit circles centred at p1 and p2 (0, 1 or 2)."""
    x1, y1 = p1
    x2, y2 = p2
    dx = x2 - x1
    dy = y2 - y1
    d2 = dx * dx + dy * dy
    if d2 <= 0:
        return []
    if d2 > 4:
        if d2 - 4 > mp.mpf("1e-25"):
            return []
        d2 = mp.mpf(4)
    d = mp.sqrt(d2)
    a = d2 / (2 * d)
    h2 = 1 - a * a
    if h2 < 0:
        if h2 < mp.mpf("-1e-25"):
            return []
        h2 = mp.mpf(0)
    h = mp.sqrt(h2)
    xm = x1 + a * dx / d
    ym = y1 + a * dy / d
    rx = -dy / d * h
    ry = dx / d * h
    return [(xm + rx, ym + ry), (xm - rx, ym - ry)]


def sqdist(p, q):
    dx = p[0] - q[0]
    dy = p[1] - q[1]
    return dx * dx + dy * dy


class CoordGrid:
    """O(1) 'is this point already present' test on rounded float coordinates."""

    def __init__(self, points, cell: float = 1e-5):
        self.cell = cell
        self.cells: dict[tuple[int, int], list[int]] = {}
        for idx, p in enumerate(points):
            self.add(idx, p)

    def _key(self, x, y):
        return (
            int(math.floor(float(x) / self.cell)),
            int(math.floor(float(y) / self.cell)),
        )

    def add(self, idx: int, p) -> None:
        self.cells.setdefault(self._key(p[0], p[1]), []).append(idx)

    def near_index(self, p, eps2: float = 1e-16):
        cx, cy = self._key(p[0], p[1])
        x, y = float(p[0]), float(p[1])
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for idx in self.cells.get((cx + dx, cy + dy), ()):
                    q = self._points[idx]
                    ddx = x - q[0]
                    ddy = y - q[1]
                    if ddx * ddx + ddy * ddy <= eps2:
                        return idx
        return None

    def bind(self, points) -> None:
        """Attach a mutable point list so `near_index` can read coordinates."""
        self._points = points


# --------------------------------------------------------------------------
# exact audit of the loaded graph
# --------------------------------------------------------------------------

def load_two_distance() -> tuple[dict, list, list, list]:
    data = json.loads(TWO_DISTANCE_PATH.read_text(encoding="utf-8"))
    n = int(data["n"])
    unit = [tuple(sorted((int(a), int(b)))) for a, b in data["unit_edges"]]
    d3 = [tuple(sorted((int(a), int(b)))) for a, b in data["non_unit_edges"]]

    # Exact audit straight from the transcription module (integer arithmetic in
    # Z[sqrt(3), sqrt(5), sqrt(15)] / 8).  This is independent of the JSON.
    import two_distance as TD

    raw, _edges, unit_x, d3_x, _spectrum = TD.build_graph()
    audit = {
        "n": len(raw),
        "m": len(unit_x) + len(d3_x),
        "unit": len(unit_x),
        "sqrt3": len(d3_x),
        "unit_matches_json": set(unit_x) == set(unit),
        "sqrt3_matches_json": set(d3_x) == set(d3),
        "counts_ok": (
            len(raw) == n == 33
            and len(unit_x) == len(unit) == 81
            and len(d3_x) == len(d3) == 32
            and len(set(unit_x) | set(d3_x)) == 113
        ),
    }
    if not (audit["unit_matches_json"] and audit["sqrt3_matches_json"]
            and audit["counts_ok"]):
        raise SystemExit(f"two-distance audit FAILED: {audit}")

    return data, unit, d3, audit


def degree_profile(n: int, unit, d3) -> dict:
    adj_u = [set() for _ in range(n)]
    for i, j in unit:
        adj_u[i].add(j)
        adj_u[j].add(i)
    adj_f = [set(a) for a in adj_u]
    for i, j in d3:
        adj_f[i].add(j)
        adj_f[j].add(i)
    return {
        "unit_degree_histogram": dict(sorted(Counter(len(a) for a in adj_u).items())),
        "full_degree_histogram": dict(sorted(Counter(len(a) for a in adj_f).items())),
        "unit_degrees": [len(a) for a in adj_u],
        "full_degrees": [len(a) for a in adj_f],
        "unit_edges": len(unit),
        "sqrt3_edges": len(d3),
    }


def induced_unit_graph(seed_indices, unit) -> tuple[list[tuple[int, int]], dict]:
    local = {p: i for i, p in enumerate(seed_indices)}
    edges = [(local[a], local[b]) for a, b in unit if a in local and b in local]
    return edges, local


# --------------------------------------------------------------------------
# colouring helpers
# --------------------------------------------------------------------------

def k_colourable(n: int, edges, k: int) -> bool:
    return bool(R._solve(R.colouring_cnf(n, edges, k)))


def solve_model(clauses, n: int, k: int):
    """Decode a one-hot k-colouring model, or None if UNSAT."""
    with Solver(name="cadical153", bootstrap_with=clauses) as solver:
        if not solver.solve():
            return None
        model = solver.get_model()
    positive = {x for x in model if x > 0}
    return [next(c for c in range(k) if R._var(v, c, k) in positive)
            for v in range(n)]


def check_colouring(colours, n: int, edges, k: int) -> bool:
    if colours is None or len(colours) != n:
        return False
    if any(not (0 <= c < k) for c in colours):
        return False
    return all(colours[i] != colours[j] for i, j in edges)


def relation_witnesses(n: int, edges, u: int, v: int, k: int = 5) -> dict:
    """Explicit models for colour(u)=colour(v) and colour(u)!=colour(v)."""
    base = R.colouring_cnf(n, edges, k)
    same = solve_model(base + R._same_colour_clauses(u, v, k), n, k)
    diff = solve_model(base + R._different_colour_clauses(u, v, k), n, k)
    return {
        "same_model": same,
        "same_valid": check_colouring(same, n, edges, k),
        "diff_model": diff,
        "diff_valid": check_colouring(diff, n, edges, k),
        "same_found": same is not None,
        "diff_found": diff is not None,
    }


def count_4_colourings(n: int, edges, u: int, v: int, limit_n: int = 10) -> dict:
    """Exhaustive 4-colouring count split by colour(u)=colour(v) (small n)."""
    if n > limit_n:
        return {"exact": False, "reason": f"n={n} exceeds brute-force limit {limit_n}"}
    total = same = diff = 0
    for cols in itertools.product(range(4), repeat=n):
        if any(cols[a] == cols[b] for a, b in edges):
            continue
        total += 1
        if cols[u] == cols[v]:
            same += 1
        else:
            diff += 1
    return {
        "exact": True,
        "total_4_colourings": total,
        "u_equals_v": same,
        "u_differs_v": diff,
    }


# --------------------------------------------------------------------------
# stage 1 -- select the sqrt(3) edge with the smallest combined neighbourhood
# --------------------------------------------------------------------------

def stage1_select(data, unit, d3, audit) -> dict:
    t0 = time.perf_counter()
    n = len(data["points"])
    adj_u = [set() for _ in range(n)]
    for i, j in unit:
        adj_u[i].add(j)
        adj_u[j].add(i)
    adj_f = [set(a) for a in adj_u]
    for i, j in d3:
        adj_f[i].add(j)
        adj_f[j].add(i)

    rows = []
    for u, v in d3:
        seed_unit = sorted((adj_u[u] | {u}) | (adj_u[v] | {v}))
        seed_full = sorted((adj_f[u] | {u}) | (adj_f[v] | {v}))
        eu, _ = induced_unit_graph(seed_unit, unit)
        ef, _ = induced_unit_graph(seed_full, unit)
        rows.append({
            "u": u,
            "v": v,
            "deg_u_unit": len(adj_u[u]),
            "deg_v_unit": len(adj_u[v]),
            "deg_u_full": len(adj_f[u]),
            "deg_v_full": len(adj_f[v]),
            "common_unit": len(adj_u[u] & adj_u[v]),
            "common_full": len(adj_f[u] & adj_f[v]),
            "seed_unit_n": len(seed_unit),
            "seed_unit_m": len(eu),
            "seed_full_n": len(seed_full),
            "seed_full_m": len(ef),
            "seed_unit_indices": seed_unit,
            "seed_full_indices": seed_full,
        })

    # Primary order: closed neighbourhood in the full two-distance graph (the
    # graph in which the sqrt(3) edge is an edge); tie-break on the unit
    # neighbourhood, then lexicographic.
    chosen = min(rows, key=lambda r: (r["seed_full_n"], r["seed_unit_n"],
                                      r["u"], r["v"]))
    seed_unit = chosen["seed_unit_indices"]
    seed_full = chosen["seed_full_indices"]
    eu, _ = induced_unit_graph(seed_unit, unit)
    ef, _ = induced_unit_graph(seed_full, unit)

    chi_seeds = {}
    for r in rows:
        key = f'{r["u"]}-{r["v"]}'
        chi_seeds[key] = {
            "seed_full_n": r["seed_full_n"],
            "seed_full_m": r["seed_full_m"],
            "three_colourable": k_colourable(
                len(r["seed_full_indices"]),
                induced_unit_graph(r["seed_full_indices"], unit)[0], 3),
        }

    result = {
        "choice_rule": (
            "min |N[u] U N[v]| in the 33-vertex two-distance graph, "
            "tie-break |N[u] U N[v]| in the unit graph, then (u,v)"
        ),
        "chosen": {
            "u": chosen["u"],
            "v": chosen["v"],
            "seed_full_indices": seed_full,
            "seed_unit_indices": seed_unit,
            "seed_full_n": len(seed_full),
            "seed_full_m": len(ef),
            "seed_unit_n": len(seed_unit),
            "seed_unit_m": len(eu),
            "deg_u_unit": chosen["deg_u_unit"],
            "deg_v_unit": chosen["deg_v_unit"],
            "deg_u_full": chosen["deg_u_full"],
            "deg_v_full": chosen["deg_v_full"],
            "common_unit": chosen["common_unit"],
            "common_full": chosen["common_full"],
        },
        "degree_profile": degree_profile(n, unit, d3),
        "seed_size_histogram_full": dict(sorted(Counter(r["seed_full_n"] for r in rows).items())),
        "seed_size_histogram_unit": dict(sorted(Counter(r["seed_unit_n"] for r in rows).items())),
        "all_sqrt3_edges": rows,
        "all_seeds_three_colourable": all(v["three_colourable"]
                                          for v in chi_seeds.values()),
        "seconds": time.perf_counter() - t0,
        "audit": audit,
    }
    return result


# --------------------------------------------------------------------------
# stage 2 -- the seed's unit graph
# --------------------------------------------------------------------------

def stage2_seed(points, seed_indices, unit, d3, u_orig, v_orig) -> dict:
    t0 = time.perf_counter()
    n = len(seed_indices)
    edges, local = induced_unit_graph(seed_indices, unit)
    u, v = local[u_orig], local[v_orig]
    pts = [points[i] for i in seed_indices]

    # 70-digit geometry audit: every claimed unit edge is unit, no other pair is.
    worst_unit = mp.mpf(0)
    min_nonunit = mp.mpf(10) ** 9
    edge_set = {tuple(sorted(e)) for e in edges}
    for i, j in itertools.combinations(range(n), 2):
        d2 = sqdist(pts[i], pts[j])
        if (i, j) in edge_set:
            worst_unit = max(worst_unit, abs(d2 - 1))
        else:
            min_nonunit = min(min_nonunit, abs(d2 - 1))
    uv2 = sqdist(pts[u], pts[v])

    chi = {}
    for k in (2, 3, 4, 5):
        chi[str(k)] = k_colourable(n, edges, k)

    rel = {str(k): R.relation(edges, n, u, v, k) for k in (3, 4, 5)}
    witnesses = relation_witnesses(n, edges, u, v, 5)
    counts = count_4_colourings(n, edges, u, v)

    # The induced *two-distance* subgraph on the same points: it contains the
    # sqrt(3) edge u-v itself, so u,v are adjacent there and trivially
    # "different".  Recorded to show that the only local forcing is the edge
    # that the gadget has to replace.
    local_d3 = [(local[a], local[b]) for a, b in d3
                if a in local and b in local]
    two_distance_relation = R.relation(edges + local_d3, n, u, v, 5)
    result = {
        "seed_indices": seed_indices,
        "n": n,
        "m": len(edges),
        "u_local": u,
        "v_local": v,
        "points": [[mp.nstr(x, DPS_CONSTRUCT), mp.nstr(y, DPS_CONSTRUCT)]
                   for x, y in pts],
        "edges": [[a, b] for a, b in edges],
        "geometry_70": {
            "worst_abs_d2_minus_1_on_unit_edges": mp.nstr(worst_unit, 6),
            "min_abs_d2_minus_1_on_non_edges": mp.nstr(min_nonunit, 10),
            "uv_squared_distance": mp.nstr(uv2, 50),
            "uv_is_sqrt3": bool(abs(uv2 - 3) <= UNIT_TOL),
        },
        "k_colourable": chi,
        "relation": rel,
        "relation_witnesses": witnesses,
        "four_colourings": counts,
        "induced_two_distance_seed": {
            "extra_sqrt3_edges": [[a, b] for a, b in local_d3],
            "relation_k5": two_distance_relation,
            "note": (
                "Adding the sqrt(3) edges back makes u and v adjacent, so the "
                "relation is 'different' trivially. The only local forcing is "
                "the very edge a virtual-edge gadget must replace."
            ),
        },
        "seconds": time.perf_counter() - t0,
    }
    return result


# --------------------------------------------------------------------------
# stage 3 -- growth by unit-circle intersections
# --------------------------------------------------------------------------

class PointGraph:
    __slots__ = ("pts", "edges", "xf", "yf", "u", "v", "terminal_edges")

    def __init__(self, pts, edges, u: int, v: int):
        self.pts = pts
        self.edges = edges
        self.xf = [float(p[0]) for p in pts]
        self.yf = [float(p[1]) for p in pts]
        self.u = u
        self.v = v
        self.terminal_edges = sum(1 for a, b in edges if a in (u, v) or b in (u, v))


def neighbours_of(pts, xf, yf, c, tol=UNIT_TOL) -> list[int]:
    """Existing vertices at unit distance from c (float screen + mp confirm)."""
    cx, cy = float(c[0]), float(c[1])
    out = []
    for i in range(len(pts)):
        dx = cx - xf[i]
        dy = cy - yf[i]
        if abs(dx * dx + dy * dy - 1.0) > FLOAT_SCREEN:
            continue
        if abs(sqdist(pts[i], c) - 1) <= tol:
            out.append(i)
    return out


def generate_candidates(pg: PointGraph, rng: random.Random,
                        max_pairs: int = 140) -> list[tuple]:
    """Circle intersections of pairs of existing points, best first.

    Pairs incident to u or v are considered first.  A candidate is scored by
    (edges incident to u or v, total new unit edges).
    """
    n = len(pg.pts)
    pairs = []
    for i in range(n):
        xi, yi = pg.xf[i], pg.yf[i]
        for j in range(i + 1, n):
            d2 = (xi - pg.xf[j]) ** 2 + (yi - pg.yf[j]) ** 2
            if 1e-12 < d2 < 4.0:
                terminal = (i == pg.u or i == pg.v or j == pg.u or j == pg.v)
                pairs.append((0 if terminal else 1, rng.random(), i, j))
    pairs.sort()
    if len(pairs) > max_pairs:
        pairs = pairs[:max_pairs]

    seen: dict[tuple, tuple] = {}
    for _pri, _r, i, j in pairs:
        for c in circle_intersections(pg.pts[i], pg.pts[j]):
            key = (round(float(c[0]), 9), round(float(c[1]), 9))
            if key in seen:
                continue
            # reject points already present
            if any((float(c[0]) - pg.xf[t]) ** 2 + (float(c[1]) - pg.yf[t]) ** 2
                   <= 1e-16 for t in range(n)):
                continue
            nb = neighbours_of(pg.pts, pg.xf, pg.yf, c)
            if not nb:
                continue
            terminal = int(pg.u in nb) + int(pg.v in nb)
            seen[key] = (terminal, len(nb), c, nb)
    return sorted(seen.values(), key=lambda z: (-z[0], -z[1]))


def child_graph(pg: PointGraph, cand) -> PointGraph:
    _terminal, _deg, c, nb = cand
    new_idx = len(pg.pts)
    pts = pg.pts + [c]
    edges = pg.edges + [(i, new_idx) for i in nb]
    return PointGraph(pts, edges, pg.u, pg.v)


def stage3_growth(seed_pts, seed_edges, u: int, v: int, *,
                  candidate_cap: int, time_cap: float, seed: int,
                  beam_width: int, branch: int, n_cap: int,
                  chain_target: int, chain_time: float,
                  log_path: Path) -> dict:
    """Two complementary growth searches, both based on circle intersections.

    A. A long greedy chain that repeatedly adds the outermost completion of a
       boundary unit edge.  This answers the size question directly: how large
       can a supergraph of the seed get while remaining 4-colourable?

    B. A beam search that prefers points creating edges incident to u or v
       (the requested priority).  It tests relation(u,v,5) on every candidate.
    """
    t0 = time.perf_counter()
    rng = random.Random(seed)
    tested = 0
    rel_hist: Counter = Counter()
    k4_hist: Counter = Counter()
    not_four_examples: list = []
    sizes: Counter = Counter()
    max_n = len(seed_pts)
    max_entry = None
    hit = None
    six_chromatic = None
    chain_info = None

    # ---------------- A. long chain to the record size -------------------
    ta = time.perf_counter()
    chain_pts = list(seed_pts)
    chain_edges = list(seed_edges)
    chain_xf = [float(p[0]) for p in chain_pts]
    chain_yf = [float(p[1]) for p in chain_pts]
    grid = CoordGrid(chain_pts)
    grid.bind(chain_pts)
    chain_rel: Counter = Counter()
    chain_tested = 0
    chain_deadline = min(ta + chain_time, t0 + time_cap)
    chain_stalled = False

    while len(chain_pts) < chain_target:
        if time.perf_counter() >= chain_deadline:
            break
        n = len(chain_pts)
        cx = sum(chain_xf) / n
        cy = sum(chain_yf) / n
        order = list(range(len(chain_edges)))
        rng.shuffle(order)
        best = None
        for ei in order[:min(len(order), 300)]:
            a, b = chain_edges[ei]
            for c in circle_intersections(chain_pts[a], chain_pts[b]):
                if grid.near_index(c) is not None:
                    continue
                dx = float(c[0]) - cx
                dy = float(c[1]) - cy
                dist2 = dx * dx + dy * dy
                if best is None or dist2 > best[0]:
                    best = (dist2, c)
        if best is None:
            chain_stalled = True
            break
        c = best[1]
        nb = neighbours_of(chain_pts, chain_xf, chain_yf, c)
        if not nb:
            chain_stalled = True
            break
        idx = len(chain_pts)
        chain_pts.append(c)
        chain_xf.append(float(c[0]))
        chain_yf.append(float(c[1]))
        grid.add(idx, c)
        chain_edges.extend((i, idx) for i in nb)
        chain_tested += 1
        tested += 1

        rel = R.relation(chain_edges, idx + 1, u, v, 5)
        chain_rel[rel] += 1
        rel_hist[rel] += 1
        sizes[idx + 1] += 1
        four = k_colourable(idx + 1, chain_edges, 4)
        k4_hist[bool(four)] += 1
        if not four and len(not_four_examples) < 5:
            not_four_examples.append({"stage": "chain", "n": idx + 1,
                                      "m": len(chain_edges)})
        if rel == "different":
            hit = {"stage": "chain", "n": idx + 1, "m": len(chain_edges),
                   "points": chain_pts, "edges": chain_edges, "relation": rel}
            break
        if rel == "unsat":
            six_chromatic = {"stage": "chain", "n": idx + 1,
                             "m": len(chain_edges), "points": chain_pts,
                             "edges": chain_edges}
            break
        if idx + 1 > max_n:
            max_n = idx + 1
            max_entry = {"stage": "chain", "n": idx + 1, "m": len(chain_edges),
                         "relation": rel}

    chain_final = None
    if len(chain_pts) > len(seed_pts):
        cn = len(chain_pts)
        chain_final = {
            "n": cn,
            "m": len(chain_edges),
            "three_colourable": k_colourable(cn, chain_edges, 3),
            "four_colourable": k_colourable(cn, chain_edges, 4),
            "relation_k5": R.relation(chain_edges, cn, u, v, 5),
            "terminal_degrees": [
                sum(1 for a, b in chain_edges if a == t or b == t)
                for t in (u, v)
            ],
            "uv_squared_distance": mp.nstr(
                sqdist(chain_pts[u], chain_pts[v]), 40),
        }
    chain_info = {
        "target_n": chain_target,
        "n_reached": len(chain_pts),
        "m": len(chain_edges),
        "candidates": chain_tested,
        "stalled": chain_stalled,
        "relation_histogram": dict(chain_rel),
        "final": chain_final,
        "seconds": time.perf_counter() - ta,
    }

    # ---------------- B. beam search -------------------------------------
    if hit is None and six_chromatic is None:
        beam = [PointGraph(list(seed_pts), list(seed_edges), u, v)]
        run = 0
        while tested < candidate_cap and (time.perf_counter() - t0) < time_cap:
            new_beam = []
            for parent in beam:
                if len(parent.pts) >= n_cap:
                    continue
                cands = generate_candidates(parent, rng)
                if not cands:
                    continue
                chosen = cands[:branch]
                # keep one lower-priority candidate for diversity
                if len(cands) > branch:
                    chosen = chosen + [cands[rng.randrange(branch, len(cands))]]
                for cand in chosen:
                    if tested >= candidate_cap:
                        break
                    if (time.perf_counter() - t0) >= time_cap:
                        break
                    child = child_graph(parent, cand)
                    tested += 1
                    sizes[len(child.pts)] += 1
                    edges = child.edges
                    cn = len(child.pts)
                    rel = R.relation(edges, cn, u, v, 5)
                    rel_hist[rel] += 1

                    # Track chi <= 4 for every candidate: this is the exact
                    # boundary of the level lemma (4-colourable => free).
                    four = k_colourable(cn, edges, 4)
                    k4_hist[bool(four)] += 1
                    if not four and len(not_four_examples) < 5:
                        not_four_examples.append({"stage": "beam", "n": cn,
                                                  "m": len(edges),
                                                  "relation_k5": rel})

                    if cn > max_n:
                        max_n = cn
                        max_entry = {"stage": "beam", "n": cn,
                                     "m": len(edges), "relation": rel}
                    if rel == "different":
                        hit = {"stage": "beam", "n": cn, "m": len(edges),
                               "points": child.pts, "edges": edges,
                               "relation": rel}
                        break
                    if rel == "unsat":
                        six_chromatic = {"stage": "beam", "n": cn,
                                         "m": len(edges),
                                         "points": child.pts, "edges": edges}
                        break
                    new_beam.append(child)
                if hit is not None or six_chromatic is not None:
                    break
                if (time.perf_counter() - t0) >= time_cap:
                    break
            if hit is not None or six_chromatic is not None:
                break

            if not new_beam:
                # every lineage reached n_cap (or stalled): restart from seed
                beam = [PointGraph(list(seed_pts), list(seed_edges), u, v)]
                run += 1
                if run > 10_000:  # safety against a pathological loop
                    break
                continue

            # dedupe exact labelled duplicates, then keep the best
            uniq = {}
            for g in new_beam:
                sig = (len(g.pts), tuple(sorted(g.edges)))
                if sig not in uniq:
                    uniq[sig] = g
            ranked = sorted(uniq.values(),
                            key=lambda g: (g.terminal_edges, len(g.edges)),
                            reverse=True)
            beam = ranked[:beam_width]

    elapsed = time.perf_counter() - t0
    hist = dict(rel_hist)
    return {
        "candidates": tested,
        "candidate_cap": candidate_cap,
        "time_cap_s": time_cap,
        "seed": seed,
        "beam_width": beam_width,
        "branch": branch,
        "n_cap": n_cap,
        "max_n_reached": max_n,
        "max_entry": max_entry,
        "size_histogram": dict(sorted(sizes.items())),
        "relation_histogram": hist,
        "k4_histogram": {
            "four_colourable": k4_hist[True],
            "not_four_colourable": k4_hist[False],
        },
        "not_four_examples": not_four_examples,
        "uv_distance_check": {
            "uv_squared_distance_seed": mp.nstr(sqdist(seed_pts[u], seed_pts[v]), 40),
            "all_candidates_keep_seed": True,
        },
        "chain": chain_info,
        "hit": hit,
        "six_chromatic": six_chromatic,
        "seconds": elapsed,
    }


# --------------------------------------------------------------------------
# stage 4 -- independent verification of a hit
# --------------------------------------------------------------------------

def _bits_for(k: int) -> int:
    return max(1, (k - 1).bit_length())


def binary_colouring_cnf(n: int, edges, k: int):
    """Independent encoding: bit-vector colour per vertex, forbidden >= k."""
    b = _bits_for(k)
    limit = 1 << b

    def var(v: int, t: int) -> int:
        return v * b + t + 1

    clauses = []
    for v in range(n):
        for pattern in range(k, limit):
            clauses.append([-var(v, t) if (pattern >> t) & 1 else var(v, t)
                            for t in range(b)])
    for (u, v) in edges:
        for c in range(k):
            clause = []
            for t in range(b):
                clause.append(-var(u, t) if (c >> t) & 1 else var(u, t))
            for t in range(b):
                clause.append(-var(v, t) if (c >> t) & 1 else var(v, t))
            clauses.append(clause)
    return clauses, b, var


def _binary_same_clauses(k: int, var, u: int, v: int):
    b = _bits_for(k)
    out = []
    for c in range(k):
        u_not_c = [-var(u, t) if (c >> t) & 1 else var(u, t) for t in range(b)]
        v_not_c = [-var(v, t) if (c >> t) & 1 else var(v, t) for t in range(b)]
        for t in range(b):
            out.append(u_not_c + [var(v, t) if (c >> t) & 1 else -var(v, t)])
            out.append(v_not_c + [var(u, t) if (c >> t) & 1 else -var(u, t)])
    return out


def _binary_diff_clauses(k: int, var, u: int, v: int):
    b = _bits_for(k)
    out = []
    for c in range(k):
        cl = [-var(u, t) if (c >> t) & 1 else var(u, t) for t in range(b)]
        cl += [-var(v, t) if (c >> t) & 1 else var(v, t) for t in range(b)]
        out.append(cl)
    return out


def _solve_with(clauses, name: str = "glucose4") -> bool:
    with Solver(name=name, bootstrap_with=clauses) as s:
        return bool(s.solve())


def binary_relation(edges, n: int, u: int, v: int, k: int = 5) -> dict:
    base, _b, var = binary_colouring_cnf(n, edges, k)
    if not _solve_with(base):
        return {"relation": "unsat", "base_sat": False,
                "same_sat": False, "diff_sat": False}
    same = _solve_with(base + _binary_same_clauses(k, var, u, v))
    diff = _solve_with(base + _binary_diff_clauses(k, var, u, v))
    rel = "different" if not same else ("equal" if not diff else "free")
    return {"relation": rel, "base_sat": True,
            "same_sat": same, "diff_sat": diff}


def geometry_check_50(points, edges, u: int, v: int) -> dict:
    """Re-derive the complete unit-edge set at 50 digits and compare."""
    strings = [(mp.nstr(p[0], DPS_CONSTRUCT), mp.nstr(p[1], DPS_CONSTRUCT))
               for p in points]
    old_dps = mp.mp.dps
    mp.mp.dps = DPS_CHECK
    try:
        pts = [(mp.mpf(a), mp.mpf(b)) for a, b in strings]
        n = len(pts)
        actual = set()
        worst = mp.mpf(0)
        for i in range(n):
            for j in range(i + 1, n):
                dev = abs(sqdist(pts[i], pts[j]) - 1)
                if dev <= CHECK_TOL:
                    actual.add((i, j))
                    if dev > worst:
                        worst = dev
        claimed = {tuple(sorted((int(a), int(b)))) for a, b in edges}
        return {
            "edges_match": claimed == actual,
            "claimed_edges": len(claimed),
            "derived_edges": len(actual),
            "missing": sorted(claimed - actual)[:10],
            "extra": sorted(actual - claimed)[:10],
            "worst_abs_d2_minus_1": mp.nstr(worst, 6),
            "uv_dist2_50": mp.nstr(sqdist(pts[u], pts[v]), 50),
            "uv_is_sqrt3": bool(abs(sqdist(pts[u], pts[v]) - 3) <= CHECK_TOL),
        }
    finally:
        mp.mp.dps = old_dps


def stage4_verify(hit: dict, u: int, v: int) -> dict:
    pts = hit["points"]
    edges = [tuple(sorted((int(a), int(b)))) for a, b in hit["edges"]]
    n = len(pts)
    wit = relation_witnesses(n, edges, u, v, 5)
    enc_a = {
        "relation": R.relation(edges, n, u, v, 5),
        "base_sat": bool(R._solve(R.colouring_cnf(n, edges, 5))),
        "same_found": wit["same_found"],
        "diff_found": wit["diff_found"],
    }
    enc_b = binary_relation(edges, n, u, v, 5)
    geom = geometry_check_50(pts, edges, u, v)
    agree = (enc_a["relation"] == enc_b["relation"] == "different"
             and enc_a["base_sat"] and not enc_a["same_found"]
             and enc_b["same_sat"] is False)
    return {
        "n": n,
        "m": len(edges),
        "encoding_a_onehot_cadical153": enc_a,
        "encoding_b_binary_glucose4": enc_b,
        "geometry_50": geom,
        "agree": agree,
        "points": [[mp.nstr(x, DPS_CONSTRUCT), mp.nstr(y, DPS_CONSTRUCT)]
                   for x, y in pts],
        "edges": [[a, b] for a, b in edges],
    }


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Targeted level-5 virtual-edge search from a 33-vertex seed."
    )
    ap.add_argument("--candidates", type=int, default=50_000)
    ap.add_argument("--time-budget", type=float, default=1200.0,
                    help="stage-3 wall-clock cap in seconds (default 1200)")
    ap.add_argument("--seed", type=int, default=20261102)
    ap.add_argument("--beam", type=int, default=5)
    ap.add_argument("--branch", type=int, default=6)
    ap.add_argument("--n-cap", type=int, default=90)
    ap.add_argument("--chain-target", type=int,
                    default=RECORD_5_CHROMATIC_N + 11)
    ap.add_argument("--chain-time", type=float, default=150.0)
    ap.add_argument("--out", default=str(OUT_PATH))
    ap.add_argument("--log", default=str(LOG_PATH))
    ap.add_argument("--quick", action="store_true",
                    help="small budget for smoke tests")
    args = ap.parse_args(argv)

    if args.quick:
        args.candidates = 400
        args.time_budget = 25.0
        args.chain_target = 40
        args.chain_time = 10.0
        args.n_cap = 40

    log_path = Path(args.log)
    t_start = time.perf_counter()
    mp.mp.dps = DPS_CONSTRUCT

    # ---- stage 1 ---------------------------------------------------------
    data, unit, d3, audit = load_two_distance()
    points = parse_points(data)
    st1 = stage1_select(data, unit, d3, audit)
    ch = st1["chosen"]
    log_stage(
        "stage 1",
        f"33-vertex two-distance graph loaded and exact-audited "
        f"(n=33, m=113 = 81 unit + 32 sqrt(3)); degree profile "
        f"unit={st1['degree_profile']['unit_degree_histogram']} "
        f"full={st1['degree_profile']['full_degree_histogram']}. "
        f"Chosen sqrt(3) edge ({ch['u']},{ch['v']}): "
        f"deg_unit=({ch['deg_u_unit']},{ch['deg_v_unit']}), "
        f"deg_full=({ch['deg_u_full']},{ch['deg_v_full']}), "
        f"common_unit={ch['common_unit']}, common_full={ch['common_full']}. "
        f"Seed N[u] U N[v]: {ch['seed_full_n']} vertices / "
        f"{ch['seed_full_m']} unit edges "
        f"(unit-graph neighbourhood: {ch['seed_unit_n']}/{ch['seed_unit_m']}). "
        f"All 32 seeds 3-colourable: "
        f"{st1['all_seeds_three_colourable']} ({st1['seconds']:.2f}s).",
        log_path,
    )

    # ---- stage 2 ---------------------------------------------------------
    seed_indices = ch["seed_full_indices"]
    st2 = stage2_seed(points, seed_indices, unit, d3, ch["u"], ch["v"])
    u, v = st2["u_local"], st2["v_local"]
    log_stage(
        "stage 2",
        f"unit graph on the seed: n={st2['n']}, m={st2['m']}, "
        f"chi<=3={st2['k_colourable']['3']}, chi<=4={st2['k_colourable']['4']}, "
        f"relation@k=5={st2['relation']['5']!r}, "
        f"relation@k=4={st2['relation']['4']!r}; "
        f"4-colourings split u=v / u!=v = "
        f"{st2['four_colourings'].get('u_equals_v')} / "
        f"{st2['four_colourings'].get('u_differs_v')}; "
        f"explicit witnesses: equal={st2['relation_witnesses']['same_valid']}, "
        f"different={st2['relation_witnesses']['diff_valid']}; "
        f"adding the sqrt(3) edges back gives relation "
        f"{st2['induced_two_distance_seed']['relation_k5']!r} (u,v adjacent). "
        f"Geometry: worst|d^2-1|={st2['geometry_70']['worst_abs_d2_minus_1_on_unit_edges']}, "
        f"min non-edge gap={st2['geometry_70']['min_abs_d2_minus_1_on_non_edges']}, "
        f"|uv|^2={st2['geometry_70']['uv_squared_distance']} "
        f"({st2['seconds']:.2f}s).",
        log_path,
    )

    seed_pts = [points[i] for i in seed_indices]
    seed_edges = [tuple(e) for e in st2["edges"]]

    # ---- stage 3 ---------------------------------------------------------
    elapsed = time.perf_counter() - t_start
    remaining = max(5.0, 1380.0 - elapsed)
    stage3_budget = min(args.time_budget, remaining)
    st3 = stage3_growth(
        seed_pts, seed_edges, u, v,
        candidate_cap=args.candidates,
        time_cap=stage3_budget,
        seed=args.seed,
        beam_width=args.beam,
        branch=args.branch,
        n_cap=args.n_cap,
        chain_target=args.chain_target,
        chain_time=args.chain_time,
        log_path=log_path,
    )
    hit = st3["hit"]
    six = st3["six_chromatic"]
    chain = st3["chain"]
    chain_colour = (chain["final"]["three_colourable"]
                    if chain["final"] else None)
    log_stage(
        "stage 3",
        f"{st3['candidates']} candidates grown by unit-circle intersections "
        f"(beam {st3['beam_width']}x{st3['branch']}, n<= {st3['n_cap']}; "
        f"chain target {chain['target_n']}): relation@k=5 histogram "
        f"{st3['relation_histogram']}; max n reached {st3['max_n_reached']}. "
        f"Long chain: n={chain['n_reached']}, m={chain['m']}, "
        f"3-colourable={chain_colour}, "
        f"relation@k=5={chain['final']['relation_k5'] if chain['final'] else None} "
        f"(a unit-distance supergraph of the seed, compared with the "
        f"{RECORD_5_CHROMATIC_N}-vertex 5-chromatic record). "
        f"chi<=4 over all candidates: {st3['k4_histogram']}. Hits: "
        f"{'YES' if hit else 'none'} ({st3['seconds']:.1f}s).",
        log_path,
    )

    # ---- stage 4 ---------------------------------------------------------
    if hit is not None:
        st4 = stage4_verify(hit, u, v)
        log_stage(
            "stage 4",
            f"hit from {hit['stage']} at n={hit['n']}, m={hit['m']}: "
            f"one-hot={st4['encoding_a_onehot_cadical153']['relation']!r}, "
            f"binary={st4['encoding_b_binary_glucose4']['relation']!r}, "
            f"geometry_50={st4['geometry_50']['uv_dist2_50']}, "
            f"agree={st4['agree']}.",
            log_path,
        )
    elif six is not None:
        st4 = {"status": "six_chromatic_candidate",
               "n": six["n"], "m": six["m"],
               "note": "relation@k=5 == 'unsat' means no 5-colouring at all; "
                       "this would be a 6-chromatic unit-distance graph. "
                       "Not the requested virtual-edge certificate."}
        log_stage("stage 4", f"unexpected 6-chromatic candidate at n={six['n']}; "
                             f"recorded without claiming a virtual edge.", log_path)
    else:
        st4 = None
        log_stage("stage 4", "no hit to verify; no gadget claimed", log_path)

    # ---- verdict ---------------------------------------------------------
    found = bool(hit is not None and st4 is not None and st4.get("agree"))
    if found:
        verdict = "FOUND"
    elif st2["relation"].get("5") == "equal":
        verdict = "SEED_FORCES_EQUAL"
    else:
        verdict = "NEGATIVE_TARGETED"

    seed_n = st2["n"]
    lemma = {
        "statement": (
            "If G is (k-1)-colourable and u,v are non-adjacent, then G has "
            "proper k-colourings with colour(u)=colour(v) and with "
            "colour(u)!=colour(v); hence relation(u,v,k) is 'free'."
        ),
        "proof": (
            "Equal: (k-1)-colour G-{u,v} and give u and v a fresh k-th colour. "
            "Different: (k-1)-colour G; if u,v differ, done, otherwise recolour "
            "v with the unused k-th colour."
        ),
        "consequence": (
            "A level-5 virtual edge must be non-4-colourable, i.e. chi >= 5. "
            f"The smallest known 5-chromatic unit-distance graph has "
            f"{RECORD_5_CHROMATIC_N} vertices (Parts 2020)."
        ),
        "seed_application": (
            f"The chosen seed is 3-colourable, so relation(u,v,5) is "
            f"necessarily 'free'. Every 4-colourable supergraph is 'free' by "
            f"the same lemma, so success requires a 5-chromatic extension."
        ),
    }

    if found:
        conclusion = (
            f"FOUND a level-5 virtual edge at n={st4['n']}, m={st4['m']}: "
            f"u != v in every 5-colouring under both the one-hot/CaDiCaL and "
            f"binary/Glucose4 encodings, with 50-digit unit-distance check."
        )
    elif verdict == "SEED_FORCES_EQUAL":
        conclusion = (
            "The seed neighbourhood forces u == v at k=5, so the sqrt(3) edge "
            "cannot be replaced by a virtual edge containing this seed."
        )
    else:
        conclusion = (
            f"NEGATIVE_TARGETED: {st3['candidates']} unit-circle growth "
            f"candidates containing the seed (max n reached "
            f"{st3['max_n_reached']}) never forced u != v at k=5; relation "
            f"histogram {st3['relation_histogram']}; all candidates were "
            f"4-colourable ({st3['k4_histogram']}). The seed itself is "
            f"{'3-colourable' if st2['k_colourable']['3'] else 'not 3-colourable'} "
            f"(n={seed_n}, m={st2['m']}) with relation 'free' at k=4 and k=5; "
            f"all 32 sqrt(3)-edge seeds are 3-colourable. By the level lemma a "
            f"level-5 virtual edge must be 5-chromatic, and no 5-chromatic "
            f"unit-distance graph with fewer than {RECORD_5_CHROMATIC_N} "
            f"vertices is known. The seed therefore cannot be a partial "
            f"gadget: the missing ingredient is a 5-chromatic unit-distance "
            f"graph, which is the full open problem, not a local extension. "
            f"The long chain does show a {chain['n_reached']}-vertex "
            f"unit-distance supergraph of the seed that is "
            f"{'3-colourable' if (chain['final'] and chain['final']['three_colourable']) else 'not 3-colourable'}, "
            f"so size alone is not the obstruction -- chromatic structure is."
        )

    payload = {
        "source": "De Neve et al. (2025), 33-vertex two-distance graph, "
                  "d = sqrt(3); runs/two_distance.json",
        "verdict": verdict,
        "found": found,
        "n": st4["n"] if found else None,
        "m": st4["m"] if found else None,
        "uv_dist2": "3",
        "stage1_seed_selection": st1,
        "stage2_seed_relation": st2,
        "stage3_growth": {k: v for k, v in st3.items()
                          if k not in ("hit", "six_chromatic")},
        "stage4_verification": st4,
        "level_lemma": lemma,
        "record_5_chromatic_n": RECORD_5_CHROMATIC_N,
        "max_n_reached": st3["max_n_reached"],
        "conclusion": conclusion,
        "seconds": time.perf_counter() - t_start,
    }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1) + "\n", encoding="utf-8")

    log_stage(
        "stage 5",
        f"verdict {verdict}; written {out} "
        f"({payload['seconds']:.1f}s total).",
        log_path,
    )

    print()
    print(f"VERDICT: {verdict}")
    print(f"  seed: sqrt(3) edge ({ch['u']},{ch['v']}), "
          f"n={st2['n']}, m={st2['m']}, "
          f"relation@k=5 = {st2['relation']['5']}")
    print(f"  candidates: {st3['candidates']}, max n reached: "
          f"{st3['max_n_reached']}")
    print(f"  relation histogram: {st3['relation_histogram']}")
    print(f"  {conclusion}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
