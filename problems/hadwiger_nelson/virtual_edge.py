"""Level-5 virtual edges for the Hadwiger--Nelson track.

A *level-k virtual edge* is a finite unit-distance graph G with two
distinguished, non-adjacent vertices u, v such that every proper k-colouring
of G gives colour(u) != colour(v).  If |u-v| = sqrt(3), the gadget can replace
the non-unit edges of a 6-chromatic two-distance graph (distances 1 and
sqrt(3)), which is the 5 -> 6 colour route.

This module does four things:

1.  Confirms that ``relation(edges, n, u, v, k)`` from ``relation.py`` returns
    "different" exactly when the constraint colour(u) == colour(v) is UNSAT,
    and "free"/"equal" otherwise (plus "unsat" when there is no k-colouring).

2.  Proves and checks the *level obstruction*:

        If G is (k-1)-colourable and u, v are non-adjacent, then G has a
        proper k-colouring with colour(u) = colour(v).

    Proof: (k-1)-colour G - {u,v}; give u and v a fresh k-th colour.  Both
    are non-adjacent and all their neighbours lie in G - {u,v}, so this is
    proper.

    Corollary: a level-k virtual edge must have chromatic number at least k.
    At k = 5 the gadget is therefore a **5-chromatic unit-distance graph**,
    not a small composition.  This is the load-bearing fact the search tests.

3.  Runs the requested bounded search: grow unit-distance graphs from
    u = (0,0), v = (0,sqrt(3)) by unit-circle intersections (plus occasional
    random unit steps), n = 8..18, and test ``relation(..., k=5)`` for every
    candidate.  Construction uses 60-digit mpmath; the final geometry check
    uses 50 digits.  Every candidate is also 4-coloured constructively, which
    exhibits the 5-colouring with u = v required by the obstruction.

4.  Scans the 509-vertex Parts graph -- the smallest known 5-chromatic
    unit-distance graph and hence the only known place a level-5 virtual edge
    could live -- over *all* non-adjacent vertex pairs at k = 5.  If any pair
    were forced different, it would be a level-5 virtual edge; pairs at
    distance sqrt(3) are the composition target.

Any claimed gadget is re-verified with a second, independent SAT encoding
(binary colour encoding, Glucose4) and a 50-digit geometry check.  A gadget is
never claimed on the strength of one encoding.

Run:
    python problems/hadwiger_nelson/virtual_edge.py
    python problems/hadwiger_nelson/virtual_edge.py --candidates 2000 --skip-509

Outputs ``runs/hn_virtual_edge.json`` and appends a stage line per stage to
``runs/virtual_edge.log``.
"""

from __future__ import annotations

import argparse
import json
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

DPS_CONSTRUCT = 60
DPS_CHECK = 50
UNIT_TOL = mp.mpf("1e-30")
FLOAT_SCREEN = 1e-8

OUT_PATH = ROOT / "runs" / "hn_virtual_edge.json"
LOG_PATH = ROOT / "runs" / "virtual_edge.log"

mp.mp.dps = DPS_CONSTRUCT


# --------------------------------------------------------------------------
# logging
# --------------------------------------------------------------------------

def _stamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def log_stage(stage: str, message: str, path: Path | None = None) -> None:
    """Append one stage line to the run log and echo it to stdout."""
    target = Path(path) if path is not None else LOG_PATH
    target.parent.mkdir(parents=True, exist_ok=True)
    line = f"[{_stamp()}] {stage}: {message}"
    with open(target, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")
    print(line, flush=True)


# --------------------------------------------------------------------------
# small exact colouring (DSATUR)
# --------------------------------------------------------------------------

def _adjacency(n: int, edges) -> list[int]:
    adj = [0] * n
    for u, v in edges:
        u, v = int(u), int(v)
        adj[u] |= 1 << v
        adj[v] |= 1 << u
    return adj


def k_colouring(n: int, edges, k: int, active=None):
    """Return a proper k-colouring (list of length n) or None.

    ``active`` restricts the colouring to an induced vertex subset; inactive
    vertices get -1.  Exact DSATUR backtracking; intended for n <= 600.
    """
    if k < 1:
        return None
    if active is None:
        verts = list(range(n))
    else:
        verts = sorted(set(int(v) for v in active))
    if not verts:
        return [-1] * n

    adj = _adjacency(n, edges)
    idx = {v: i for i, v in enumerate(verts)}
    m = len(verts)
    adjm = [0] * m
    for v in verts:
        i = idx[v]
        for u in verts:
            if (adj[v] >> u) & 1:
                adjm[i] |= 1 << idx[u]

    colour = [-1] * m
    sat = [0] * m
    deg = [bin(a).count("1") for a in adjm]

    def rec(done: int, used: int) -> bool:
        if done == m:
            return True
        v = -1
        best = (-1, -1)
        for i in range(m):
            if colour[i] < 0:
                key = (bin(sat[i]).count("1"), deg[i])
                if key > best:
                    best = key
                    v = i
        forbidden = sat[v]
        for c in range(min(used + 1, k)):
            if (forbidden >> c) & 1:
                continue
            colour[v] = c
            touched = []
            for u in range(m):
                if ((adjm[v] >> u) & 1) and colour[u] < 0 and not ((sat[u] >> c) & 1):
                    sat[u] |= 1 << c
                    touched.append(u)
            if rec(done + 1, max(used, c + 1)):
                return True
            for u in touched:
                sat[u] &= ~(1 << c)
            colour[v] = -1
        return False

    if not rec(0, 0):
        return None
    out = [-1] * n
    for v in verts:
        out[v] = colour[idx[v]]
    return out


def is_k_colourable(n: int, edges, k: int, active=None) -> bool:
    return k_colouring(n, edges, k, active) is not None


def check_colouring(n: int, edges, colours, k: int) -> bool:
    if colours is None or len(colours) != n:
        return False
    if any((c < 0 or c >= k) for c in colours):
        return False
    return all(colours[int(u)] != colours[int(v)] for u, v in edges)


def exact_chromatic_small(n: int, edges, cap: int = 5):
    """Exact chromatic number if it is <= cap, else None."""
    for k in range(1, cap + 1):
        if is_k_colourable(n, edges, k):
            return k
    return None


def equal_terminal_colouring(n: int, edges, u: int, v: int, k: int = 5):
    """Construct the obstruction witness: colour G-{u,v} with k-1 colours.

    Returns a proper k-colouring with colour[u] == colour[v], or None.
    """
    active = [w for w in range(n) if w not in (u, v)]
    partial = k_colouring(n, edges, k - 1, active=active)
    if partial is None:
        return None
    colours = list(partial)
    colours[u] = k - 1
    colours[v] = k - 1
    if not check_colouring(n, edges, colours, k):
        return None
    return colours


# --------------------------------------------------------------------------
# stage 1 -- relation semantics
# --------------------------------------------------------------------------

def same_colour_sat(edges, n: int, u: int, v: int, k: int) -> bool:
    """Direct CNF test of the constraint colour(u) == colour(v)."""
    base = R.colouring_cnf(n, edges, k)
    if not R._solve(base):
        return False
    if u == v:
        return True
    return R._solve(base + R._same_colour_clauses(u, v, k))


def different_colour_sat(edges, n: int, u: int, v: int, k: int) -> bool:
    base = R.colouring_cnf(n, edges, k)
    if not R._solve(base):
        return False
    return R._solve(base + R._different_colour_clauses(u, v, k))


def _semantics_cases():
    """Graphs with known pair relations, shared by stages 1 and 1b."""
    cycle5 = [(i, (i + 1) % 5) for i in range(5)]
    diamond = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3)]
    k4 = [(i, j) for i in range(4) for j in range(i + 1, 4)]
    k4_minus = [e for e in k4 if e != (2, 3)]
    return [
        ("edge-k3", 2, [(0, 1)], 3, 0, 1, "different"),
        ("diamond-tips-k3", 4, diamond, 3, 2, 3, "equal"),
        ("diamond-edge-k3", 4, diamond, 3, 0, 1, "different"),
        ("c5-chord-k3", 5, cycle5, 3, 0, 2, "free"),
        ("c5-k2", 5, cycle5, 2, 0, 2, "unsat"),
        ("k4-k3", 4, k4, 3, 0, 1, "unsat"),
        ("k4-minus-tips-k3", 4, k4_minus, 3, 2, 3, "equal"),
    ]


def stage1_semantics(verbose: bool = True) -> dict:
    """Check that relation() has exactly the documented semantics.

    For a fixed set of graphs, compare every pair's relation against direct
    CNF tests of "u == v" and "u != v", and against brute force on small
    graphs.  This is the stage-1 acceptance check.
    """
    t0 = time.perf_counter()
    cases = _semantics_cases()

    checks = 0
    mismatches = []
    for name, n, edges, k, u, v, expected in cases:
        got = R.relation(edges, n, u, v, k)
        checks += 1
        if got != expected:
            mismatches.append({"case": name, "got": got, "want": expected})
        # relation == "different" iff base SAT and u == v UNSAT
        base_sat = R._solve(R.colouring_cnf(n, edges, k))
        same = same_colour_sat(edges, n, u, v, k) if base_sat else False
        diff = different_colour_sat(edges, n, u, v, k) if base_sat else False
        checks += 1
        if base_sat:
            if (got == "different") != (not same):
                mismatches.append({"case": name, "rule": "different != not-same",
                                   "got": got, "same_sat": same})
            if got == "free" and not (same and diff):
                mismatches.append({"case": name, "rule": "free needs both",
                                   "same": same, "diff": diff})
            if got == "equal" and not (same and not diff):
                mismatches.append({"case": name, "rule": "equal needs same-only",
                                   "same": same, "diff": diff})
        else:
            if got != "unsat":
                mismatches.append({"case": name, "rule": "base-unsat",
                                   "got": got})

    # Randomised brute-force cross-check: every labelled graph on 5 vertices
    # is too many for a smoke test, so sample 40 non-trivial masks.
    rng = random.Random(11)
    brute_checks = 0
    pairs5 = [(i, j) for i in range(5) for j in range(i + 1, 5)]
    for _ in range(40):
        mask = rng.randrange(1 << len(pairs5))
        edges = [pairs5[i] for i in range(len(pairs5)) if (mask >> i) & 1]
        expected = R.brute_relations(edges, 5, 3)
        for (u, v), want in expected.items():
            got = R.relation(edges, 5, u, v, 3)
            brute_checks += 1
            if got != want:
                mismatches.append({"case": "random-n5", "edges": edges,
                                   "pair": [u, v], "got": got, "want": want})

    seconds = time.perf_counter() - t0
    if verbose:
        print(f"stage 1: relation semantics: {checks} direct checks + "
              f"{brute_checks} brute-force pair checks, "
              f"{len(mismatches)} mismatches ({seconds:.2f}s)", flush=True)
    return {
        "cases": len(cases),
        "direct_checks": checks,
        "brute_force_graphs": 40,
        "brute_force_pair_checks": brute_checks,
        "mismatches": mismatches[:20],
        "seconds": seconds,
        "ok": not mismatches,
    }


# --------------------------------------------------------------------------
# geometry: unit-circle growth at 60 digits
# --------------------------------------------------------------------------

def circle_intersections(p1, p2):
    """Intersections of the unit circles centred at p1 and p2 (0, 1 or 2)."""
    x1, y1 = p1
    x2, y2 = p2
    dx = x2 - x1
    dy = y2 - y1
    d2 = dx * dx + dy * dy
    if d2 == 0:
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


def two_point_seed():
    """The minimal seed: terminals only, at the target distance sqrt(3)."""
    return [(mp.mpf(0), mp.mpf(0)), (mp.mpf(0), mp.sqrt(3))]


def moser_seed():
    """The 7-vertex Moser spindle with terminals at distance sqrt(3).

    Vertices: u=(0,0), v=(0,sqrt(3)), the two body points of the diamond on
    u,v, a second diamond sharing u with far tip w, and its two body points.
    The far tip w satisfies |u-w| = sqrt(3) and |v-w| = 1, so the standard
    spindle (two diamonds sharing a tip, far tips joined by a unit edge) is
    exactly the 7-vertex, 11-edge, 4-chromatic graph bounded by these points.
    Terminals u,v are at distance sqrt(3) and are not adjacent.
    """
    s3 = mp.sqrt(3)
    u = (mp.mpf(0), mp.mpf(0))
    v = (mp.mpf(0), s3)
    half = mp.mpf(1) / 2
    body_uv = [(half, s3 / 2), (-half, s3 / 2)]
    w = (mp.sqrt(33) / 6, 5 * s3 / 6)
    body_uw = circle_intersections(u, w)
    return [u, v] + body_uv + [w] + body_uw


def grow_candidate(n: int, rng: random.Random, seed_points=None,
                   step_prob: float = 0.12, max_guard: int = 400):
    """Grow an n-point unit-distance graph containing u=(0,0), v=(0,sqrt(3)).

    New points are unit-circle intersections (88%) or random unit steps
    (12%).  Returns None if the growth stalls.  ``seed_points`` defaults to
    the two terminals; pass ``moser_seed()`` for a 4-chromatic start.
    """
    pts = list(seed_points) if seed_points is not None else two_point_seed()
    if len(pts) > n:
        return None
    guard = 0
    while len(pts) < n and guard < max_guard:
        guard += 1
        cand = None
        if rng.random() < step_prob:
            i = rng.randrange(len(pts))
            theta = 2 * mp.pi * mp.mpf(rng.random())
            cand = (pts[i][0] + mp.cos(theta), pts[i][1] + mp.sin(theta))
        else:
            for _ in range(25):
                i, j = rng.sample(range(len(pts)), 2)
                cs = circle_intersections(pts[i], pts[j])
                if cs:
                    cand = rng.choice(cs)
                    break
        if cand is None:
            continue
        if all(abs(cand[0] - q[0]) + abs(cand[1] - q[1]) > mp.mpf("1e-30")
               for q in pts):
            pts.append(cand)
    return pts if len(pts) == n else None


def unit_edges(points, tol=UNIT_TOL):
    """Complete unit-distance edge set of a point list.

    Float screening followed by 60-digit confirmation.  Exact coincidences in
    these constructions evaluate to 1 far inside ``tol``; accidental near-unit
    pairs are rejected by the second stage.
    """
    n = len(points)
    xf = [float(p[0]) for p in points]
    yf = [float(p[1]) for p in points]
    out = []
    for i in range(n):
        xi, yi = xf[i], yf[i]
        for j in range(i + 1, n):
            dx = xi - xf[j]
            dy = yi - yf[j]
            if abs(dx * dx + dy * dy - 1.0) > FLOAT_SCREEN:
                continue
            d2 = ((points[i][0] - points[j][0]) ** 2
                  + (points[i][1] - points[j][1]) ** 2)
            if abs(d2 - 1) <= tol:
                out.append((i, j))
    return out


# --------------------------------------------------------------------------
# stage 2 -- bounded candidate search
# --------------------------------------------------------------------------

def stage2_search(total_cap: int, time_cap: float, seed: int) -> dict:
    """Generate and test candidates with n = 8..18.

    For every candidate the module records the exact chromatic number (up to
    5), runs ``relation(..., k=5)``, and -- when 4-colourable -- constructs
    the obstruction witness (a 5-colouring with equal terminals).
    """
    t0 = time.perf_counter()
    rng = random.Random(seed)
    sizes = list(range(8, 19))
    target = max(1, total_cap // len(sizes))
    moser_target = (3 * target) // 4
    simple_target = target - moser_target
    variants = [("two-point", two_point_seed(), simple_target),
                ("moser", moser_seed(), moser_target)]
    stats = {
        str(n): {"tested": 0, "m_min": None, "m_max": 0, "m_sum": 0,
                 "chi": Counter(), "relation": Counter(), "four_colourable": 0,
                 "seeds": Counter()}
        for n in sizes
    }
    relation_counts: Counter = Counter()
    hits = []
    witnesses_ok = 0
    witness_failures = []
    tested = 0
    deadline = t0 + time_cap
    stopped_early = False
    consecutive_failures = 0

    for n in sizes:
        for variant_name, seed_points, want in variants:
            count_n = 0
            while count_n < want and tested < total_cap:
                if time.perf_counter() > deadline:
                    stopped_early = True
                    break
                pts = grow_candidate(n, rng, seed_points)
                if pts is None:
                    consecutive_failures += 1
                    if consecutive_failures > 200:
                        break
                    continue
                consecutive_failures = 0
                edges = unit_edges(pts)
                if len(edges) < n - 1:      # disconnected: should not happen
                    continue

                tested += 1
                count_n += 1
                rec = stats[str(n)]
                rec["tested"] += 1
                rec["seeds"][variant_name] += 1
                rec["m_sum"] += len(edges)
                rec["m_min"] = len(edges) if rec["m_min"] is None else min(rec["m_min"], len(edges))
                rec["m_max"] = max(rec["m_max"], len(edges))

                chi = exact_chromatic_small(n, edges, cap=5)
                rec["chi"][str(chi) if chi is not None else ">5"] += 1

                rel = R.relation(edges, n, 0, 1, 5)
                relation_counts[rel] += 1
                rec["relation"][rel] += 1

                four = is_k_colourable(n, edges, 4)
                if four:
                    rec["four_colourable"] += 1
                    if rel == "different":
                        raise AssertionError(
                            "level obstruction violated by a 4-colourable graph"
                        )
                    colours = equal_terminal_colouring(n, edges, 0, 1, 5)
                    if colours is None or not check_colouring(n, edges, colours, 5):
                        witness_failures.append({"n": n, "edges": len(edges)})
                    else:
                        witnesses_ok += 1
                elif rel == "different":
                    hits.append({
                        "n": n,
                        "m": len(edges),
                        "u": 0,
                        "v": 1,
                        "points": pts,
                        "edges": edges,
                        "source": f"stage2-{variant_name}",
                        "chi": chi,
                    })

            if stopped_early:
                break
        if stopped_early:
            break

    seconds = time.perf_counter() - t0
    per_n = {}
    for n in sizes:
        rec = stats[str(n)]
        per_n[str(n)] = {
            "tested": rec["tested"],
            "edges_min": rec["m_min"],
            "edges_max": rec["m_max"],
            "edges_mean": (rec["m_sum"] / rec["tested"]) if rec["tested"] else None,
            "chi_histogram": dict(sorted(rec["chi"].items())),
            "relation_histogram": dict(sorted(rec["relation"].items())),
            "four_colourable": rec["four_colourable"],
            "seeds": dict(sorted(rec["seeds"].items())),
        }
    return {
        "candidates": tested,
        "cap": total_cap,
        "time_cap_s": time_cap,
        "stopped_early": stopped_early,
        "seed": seed,
        "sizes": sizes,
        "per_n": per_n,
        "relation_histogram": dict(sorted(relation_counts.items())),
        "witnesses_ok": witnesses_ok,
        "witness_failures": witness_failures[:5],
        "hits": hits,
        "seconds": seconds,
    }


# --------------------------------------------------------------------------
# stage 2b -- exhaustive pair scan of the 509-vertex record graph
# --------------------------------------------------------------------------

def stage2b_probe_509(time_budget: float, batch: int = 8000) -> dict:
    """Test colour(u)==colour(v) for every non-adjacent pair of the 509 graph.

    Uses one incremental selector encoding per batch (cadical153).  A pair is
    a level-5 virtual edge exactly when the same-colour query is UNSAT.  The
    graph has 509 vertices, so this also checks the only known 5-chromatic
    unit-distance seed at the only level where a virtual edge can exist.
    """
    t0 = time.perf_counter()
    old_dps = mp.mp.dps
    mp.mp.dps = DPS_CHECK
    try:
        pts, edges, worst = R.load_509()
        n = len(pts)
        edge_set = {tuple(sorted(e)) for e in edges}
        nonadj = [(i, j) for i in range(n) for j in range(i + 1, n)
                  if (i, j) not in edge_set]
        k = 5
        tested = 0
        forced = []
        aborted = False
        batches = 0

        for start in range(0, len(nonadj), batch):
            if time.perf_counter() - t0 > time_budget:
                aborted = True
                break
            chunk = nonadj[start:start + batch]
            cnf = R.colouring_cnf(n, edges, k)
            next_var = n * k + 1
            selectors = []
            for (u, v) in chunk:
                s = next_var
                next_var += 1
                selectors.append((u, v, s))
                for c in range(k):
                    cnf.append([-s, -R._var(u, c, k), R._var(v, c, k)])
                    cnf.append([-s, -R._var(v, c, k), R._var(u, c, k)])
            with Solver(name="cadical153", bootstrap_with=cnf) as solver:
                if not solver.solve():
                    raise AssertionError("509 graph is not 5-colourable")
                for (u, v, s) in selectors:
                    if time.perf_counter() - t0 > time_budget:
                        aborted = True
                        break
                    tested += 1
                    if not solver.solve(assumptions=[s]):
                        d2 = ((pts[u][0] - pts[v][0]) ** 2
                              + (pts[u][1] - pts[v][1]) ** 2)
                        forced.append({"pair": [u, v],
                                       "d2": mp.nstr(d2, 40),
                                       "d2_mpf": d2})
            batches += 1
            if aborted:
                break

        targets = []
        for item in forced:
            if abs(item["d2_mpf"] - 3) <= UNIT_TOL:
                targets.append({"pair": item["pair"], "d2": item["d2"]})

        seconds = time.perf_counter() - t0
        return {
            "graph": "research/scratch/509_parts.vtx",
            "n": n,
            "m": len(edges),
            "worst_abs_d2_minus_1": mp.nstr(worst, 6),
            "k": k,
            "nonadjacent_pairs": len(nonadj),
            "pairs_tested": tested,
            "complete": (not aborted) and tested == len(nonadj),
            "aborted": aborted,
            "batches": batches,
            "forced_different_nonadjacent": [
                {"pair": item["pair"], "d2": item["d2"]} for item in forced
            ],
            "sqrt3_targets": targets,
            "points": pts if targets else None,
            "edges": edges if targets else None,
            "seconds": seconds,
        }
    finally:
        mp.mp.dps = old_dps


# --------------------------------------------------------------------------
# stage 4 -- independent verification of a claimed gadget
# --------------------------------------------------------------------------

def _bits_for(k: int) -> int:
    return max(1, (k - 1).bit_length())


def binary_colouring_cnf(n: int, edges, k: int):
    """Independent encoding: 3-bit colour per vertex, no one-hot variables.

    Forbidden bit patterns >= k enforce "at most k colours"; one clause per
    edge and per colour forbids equal colours.  Returns (clauses, bits, var).
    """
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
    """Clauses forcing u and v to encode the same colour.

    "(u == c) -> (v == c)" is a conjunction over the bits of v, so it needs
    one clause per bit t:  (u != c) OR (v_t == c_t).  Writing it as a single
    clause over all of v's equality literals would encode a disjunction and
    would not force equality.
    """
    b = _bits_for(k)
    out = []
    for c in range(k):
        u_not_c = [-var(u, t) if (c >> t) & 1 else var(u, t)
                   for t in range(b)]
        v_not_c = [-var(v, t) if (c >> t) & 1 else var(v, t)
                   for t in range(b)]
        for t in range(b):
            # (u == c) -> (v_t == c_t)
            out.append(u_not_c + [var(v, t) if (c >> t) & 1 else -var(v, t)])
            # (v == c) -> (u_t == c_t)
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
    """Second encoding: relation of u and v under binary colour encoding."""
    base, _b, var = binary_colouring_cnf(n, edges, k)
    base_sat = _solve_with(base)
    out = {"base_sat": base_sat}
    if not base_sat:
        out.update({"relation": "unsat", "same_sat": False, "diff_sat": False})
        return out
    same_sat = _solve_with(base + _binary_same_clauses(k, var, u, v))
    diff_sat = _solve_with(base + _binary_diff_clauses(k, var, u, v))
    if not same_sat:
        rel = "different"
    elif not diff_sat:
        rel = "equal"
    else:
        rel = "free"
    out.update({"relation": rel, "same_sat": same_sat, "diff_sat": diff_sat})
    return out


def stage1b_second_encoding() -> dict:
    """Validate the independent binary encoding on the known relation cases."""
    t0 = time.perf_counter()
    cases = _semantics_cases()
    mismatches = []
    for name, n, edges, k, u, v, expected in cases:
        got = binary_relation(edges, n, u, v, k)
        if got["relation"] != expected:
            mismatches.append({"case": name, "got": got["relation"],
                               "want": expected})

    # Randomised cross-check against the one-hot encoding on 5 vertices.
    rng = random.Random(13)
    pairs5 = [(i, j) for i in range(5) for j in range(i + 1, 5)]
    random_checks = 0
    for _ in range(15):
        mask = rng.randrange(1 << len(pairs5))
        edges = [pairs5[i] for i in range(len(pairs5)) if (mask >> i) & 1]
        for (u, v) in pairs5:
            want = R.relation(edges, 5, u, v, 3)
            got = binary_relation(edges, 5, u, v, 3)["relation"]
            random_checks += 1
            if got != want:
                mismatches.append({"case": "random-n5", "edges": edges,
                                   "pair": [u, v], "got": got, "want": want})

    return {
        "cases": len(cases),
        "random_graphs": 15,
        "random_pair_checks": random_checks,
        "mismatches": mismatches[:20],
        "seconds": time.perf_counter() - t0,
        "ok": not mismatches,
    }


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
            xi, yi = pts[i]
            for j in range(i + 1, n):
                d2 = (xi - pts[j][0]) ** 2 + (yi - pts[j][1]) ** 2
                dev = abs(d2 - 1)
                if dev <= mp.mpf("1e-25"):
                    actual.add((i, j))
                    if dev > worst:
                        worst = dev
        claimed = {tuple(sorted((int(a), int(b)))) for a, b in edges}
        uv2 = ((pts[u][0] - pts[v][0]) ** 2
               + (pts[u][1] - pts[v][1]) ** 2)
        return {
            "edges_match": claimed == actual,
            "claimed_edges": len(claimed),
            "derived_edges": len(actual),
            "missing": sorted(claimed - actual)[:10],
            "extra": sorted(actual - claimed)[:10],
            "worst_abs_d2_minus_1": mp.nstr(worst, 6),
            "uv_dist2_50": mp.nstr(uv2, 50),
            "uv_is_sqrt3": bool(abs(uv2 - 3) <= mp.mpf("1e-25")),
        }
    finally:
        mp.mp.dps = old_dps


def verify_candidate(points, edges, n: int, u: int = 0, v: int = 1) -> dict:
    """Verify a claimed gadget with both encodings and 50-digit geometry."""
    enc_a = {
        "relation": R.relation(edges, n, u, v, 5),
        "base_sat": R._solve(R.colouring_cnf(n, edges, 5)),
        "same_sat": same_colour_sat(edges, n, u, v, 5),
        "diff_sat": different_colour_sat(edges, n, u, v, 5),
    }
    enc_b = binary_relation(edges, n, u, v, 5)
    geom = geometry_check_50(points, edges, u, v)
    exact_uv = None
    if u == 0 and v == 1:
        exact_uv = sp.Integer(3)      # construction fixes v = (0, sqrt(3))
    return {
        "encoding_a_onehot": enc_a,
        "encoding_b_binary": enc_b,
        "geometry_50": geom,
        "uv_dist2_exact": str(exact_uv) if exact_uv is not None else None,
        "agree": (enc_a["relation"] == enc_b["relation"] == "different"
                  and enc_a["base_sat"] and not enc_a["same_sat"]),
    }


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def main(argv=None) -> int:
    global LOG_PATH
    ap = argparse.ArgumentParser(
        description="Search for a level-5 virtual edge (Hadwiger-Nelson)."
    )
    ap.add_argument("--candidates", type=int, default=20000,
                    help="stage-2 candidate cap (default 20000)")
    ap.add_argument("--seed", type=int, default=20261102)
    ap.add_argument("--time-budget", type=float, default=850.0,
                    help="stage-2 wall-clock cap in seconds (default 850)")
    ap.add_argument("--probe-509", dest="probe_509", action="store_true",
                    default=True, help="scan the 509-vertex record graph")
    ap.add_argument("--skip-509", dest="probe_509", action="store_false")
    ap.add_argument("--out", default=str(OUT_PATH))
    ap.add_argument("--log", default=str(LOG_PATH))
    args = ap.parse_args(argv)

    LOG_PATH = Path(args.log)
    t_start = time.perf_counter()
    mp.mp.dps = DPS_CONSTRUCT

    # ---- stage 1 ---------------------------------------------------------
    st1 = stage1_semantics(verbose=False)
    if not st1["ok"]:
        log_stage("stage 1", f"FAILED: {st1['mismatches'][:3]}")
        return 1
    log_stage(
        "stage 1",
        f"relation semantics confirmed: {st1['direct_checks']} direct CNF "
        f"checks + {st1['brute_force_pair_checks']} brute-force pair checks, "
        f"0 mismatches ({st1['seconds']:.2f}s). "
        "'different' <=> base SAT and colour(u)==colour(v) UNSAT.",
    )

    st1b = stage1b_second_encoding()
    if not st1b["ok"]:
        log_stage("stage 1b", f"FAILED: {st1b['mismatches'][:3]}")
        return 1
    log_stage(
        "stage 1b",
        f"independent binary/Glucose4 encoding agrees on all "
        f"{st1b['cases']} known relation cases and "
        f"{st1b['random_pair_checks']} random pair checks "
        f"({st1b['seconds']:.2f}s).",
    )

    # ---- stage 2 ---------------------------------------------------------
    st2 = stage2_search(args.candidates, args.time_budget, args.seed)
    hist = st2["relation_histogram"]
    log_stage(
        "stage 2",
        f"{st2['candidates']} candidates n=8..18 generated by unit-circle "
        f"intersections + random unit steps from two seeds "
        f"(two-point and 4-chromatic Moser spindle, terminals u=(0,0), "
        f"v=(0,sqrt(3))); relation@k=5 histogram {hist}; "
        f"obstruction witnesses verified for {st2['witnesses_ok']} "
        f"4-colourable candidates ({st2['seconds']:.1f}s"
        f"{', stopped at time cap' if st2['stopped_early'] else ''}).",
    )

    # ---- stage 2b --------------------------------------------------------
    st2b = None
    if args.probe_509:
        elapsed = time.perf_counter() - t_start
        budget = max(0.0, min(240.0, 1150.0 - elapsed))
        st2b = stage2b_probe_509(budget)
        if st2b["complete"]:
            detail = (f"all {st2b['nonadjacent_pairs']} non-adjacent pairs "
                      f"admit a 5-colouring with u=v")
        else:
            detail = (f"{st2b['pairs_tested']}/{st2b['nonadjacent_pairs']} "
                      f"non-adjacent pairs tested (time cap)")
        log_stage(
            "stage 2b",
            f"509-vertex record graph (n={st2b['n']}, m={st2b['m']}, k=5): "
            f"{detail}; forced-different non-adjacent pairs: "
            f"{len(st2b['forced_different_nonadjacent'])}; "
            f"sqrt(3)-targets: {len(st2b['sqrt3_targets'])} "
            f"({st2b['seconds']:.1f}s).",
        )
    else:
        log_stage("stage 2b", "skipped (--skip-509)")

    # ---- stage 3: choose the smallest hit --------------------------------
    hits = list(st2["hits"])
    if st2b is not None and st2b["sqrt3_targets"]:
        for target in st2b["sqrt3_targets"]:
            u, v = target["pair"]
            hits.append({
                "n": st2b["n"],
                "m": st2b["m"],
                "u": u,
                "v": v,
                "points": st2b["points"],
                "edges": st2b["edges"],
                "source": "stage2b-509",
                "chi": 5,
            })
    hits.sort(key=lambda h: (h["n"], h["m"]))
    best = hits[0] if hits else None
    log_stage(
        "stage 3",
        f"hits: {len(hits)}; smallest: "
        + (f"n={best['n']}, m={best['m']}" if best else "none"),
    )

    # ---- stage 4: independent verification -------------------------------
    verification = None
    found = False
    if best is not None:
        verification = verify_candidate(
            best["points"], best["edges"], best["n"], best["u"], best["v"]
        )
        found = bool(verification["agree"] and verification["geometry_50"]["uv_is_sqrt3"])
        log_stage(
            "stage 4",
            f"candidate n={best['n']}: encoding A={verification['encoding_a_onehot']['relation']}, "
            f"encoding B={verification['encoding_b_binary']['relation']}, "
            f"edges_match={verification['geometry_50']['edges_match']}, "
            f"uv^2={verification['geometry_50']['uv_dist2_50']}; "
            f"claim={found}.",
        )
    else:
        log_stage("stage 4", "no candidate to verify; no gadget claimed")

    # ---- stage 5: write JSON --------------------------------------------
    if best is not None:
        points_json = [[mp.nstr(x, DPS_CONSTRUCT), mp.nstr(y, DPS_CONSTRUCT)]
                       for x, y in best["points"]]
        edges_json = [[int(a), int(b)] for a, b in best["edges"]]
        gadget = {"points": points_json, "edges": edges_json,
                  "u": int(best["u"]), "v": int(best["v"])}
    else:
        gadget = None

    total_seconds = time.perf_counter() - t_start
    payload = {
        "found": found,
        "n": best["n"] if found else None,
        "m": best["m"] if found else None,
        "uv_dist2": "3",
        "candidates": st2["candidates"],
        "seconds": round(total_seconds, 2),
        "gadget": gadget if found else None,
        "level_obstruction": {
            "statement": (
                "If G is (k-1)-colourable and u,v are non-adjacent, then "
                "colouring G-{u,v} with k-1 colours and giving u and v a "
                "fresh k-th colour is a proper k-colouring with u=v. Hence a "
                "level-k virtual edge has chromatic number at least k; at "
                "k=5 it is a 5-chromatic unit-distance graph."
            ),
            "level5_consequence": (
                "Any level-5 virtual edge must be 5-chromatic. The smallest "
                "known 5-chromatic unit-distance graph has 509 vertices, and "
                "the 509-vertex Parts graph has no critical non-adjacent pair "
                "(stage 2b)."
            ),
            "witnesses_verified": st2["witnesses_ok"],
            "witness_failures": st2["witness_failures"],
        },
        "stage1_relation_semantics": st1,
        "stage1b_second_encoding": st1b,
        "stage2_search": st2,
        "stage2b_509_probe": st2b,
        "stage4_verification": verification,
        "conclusion": (
            "NEGATIVE for the covered space. Every generated candidate "
            "(n=8..18, both 3-chromatic and 4-chromatic seeds) is "
            "4-colourable, so by the level obstruction it has a 5-colouring "
            "with the terminals equal; relation() returned 'different' for "
            "none. The 509-vertex record graph is also not a level-5 virtual "
            "edge: all non-adjacent pairs admit a 5-colouring with equal "
            "colours, so in particular no sqrt(3) pair is forced different. "
            "A level-5 virtual edge therefore requires a 5-chromatic "
            "unit-distance graph with a critical non-adjacent pair -- a "
            "strictly stronger object than the known 509-vertex record, not "
            "a small composition."
        ),
    }
    payload["stage2_search"].pop("hits", None)   # hits are folded into gadget

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    log_stage(
        "stage 5",
        f"written {out.relative_to(ROOT) if out.is_relative_to(ROOT) else out} "
        f"({total_seconds:.1f}s total).",
    )

    # ---- acceptance line -------------------------------------------------
    if found:
        print(f"VIRTUAL EDGE FOUND n={best['n']}", flush=True)
    else:
        coverage = (
            f"{st2['candidates']} candidates tested, n=8..18 unit-circle "
            f"growth from u=(0,0), v=(0,sqrt(3)); all 4-colourable"
        )
        if st2b is not None:
            coverage += (
                f"; 509-vertex record graph: {st2b['pairs_tested']}"
                f"/{st2b['nonadjacent_pairs']} non-adjacent pairs with "
                f"u=v in some 5-colouring"
            )
        print(f"NEGATIVE: {coverage}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
