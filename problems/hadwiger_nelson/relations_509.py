"""Census of forced colour relations in the 509-vertex Parts graph at k = 5.

Why k = 5.  The graph is 5-chromatic: its 4-colouring instance is UNSAT and its
5-colouring instance is SAT (research/hn-509-verify.md).  At k = 3 and k = 4
every pair therefore answers "unsat", which says nothing about the colouring.
The informative question is what a proper *5*-colouring is forced to do, and
that is what this script measures.

What is measured
----------------
1.  The graph is loaded with the existing `relation.load_509` and independently
    rebuilt with the exact distance code in `tools/verify_509.py`; the two edge
    sets must agree.
2.  Pairs probed at k = 5, in three strata:
        all 2442 edges
        all 14720 pairs at graph distance 2
        2000 distinct uniformly random pairs (fixed seed)
    Every pair gets `relation(..., k=5)` semantics: "equal", "different",
    "free" (or "unsat", which cannot happen here and is checked).
3.  The forced-equality graph on the 509 vertices has an edge for every pair
    with relation "equal"; its connected components are reported.
    Equality is transitive in every colouring, so a component of size >= 3 is
    a chain of forced equalities, i.e. a genuinely stronger constraint.
4.  Diamonds: every 4-vertex induced K4 minus an edge whose two apices are at
    squared distance 3.  A diamond forces its apices to share a colour at
    k = 3 (r, s are adjacent and use two colours; p, q must use the third).
    Each diamond is checked as a 4-vertex gadget at k = 3, k = 4 and then as a
    pair inside the full 509-vertex graph at k = 5.

Method and cost
---------------
`relation()` rebuilds and re-solves the whole CNF per call (~25-30 ms here).
The bulk loop therefore uses one incremental solver with selector literals
(same CNF, same variable numbering, assumption-driven), which is semantically
identical and several times faster (~8 ms per pair at this size, including the
301k-clause selector build).  The shortcut is not trusted: a stratified
sample of pairs is re-probed with `relation()` itself and the two answers must
agree, and all diamond apex pairs are re-probed with `relation()` as well.
The total wall clock is capped at 900 s (default); the cap, the elapsed time
and whether it was hit are recorded in the JSON.

Usage
-----
    python problems/hadwiger_nelson/relations_509.py
    python problems/hadwiger_nelson/relations_509.py --cap 300 --random-pairs 500

Writes runs/hn_509_relations.json.  Exit code 0 means the census ran to
completion and every internal check passed.
"""

from __future__ import annotations

import argparse
import collections
import itertools
import json
import random
import sys
import time
from pathlib import Path

import mpmath as mp

mp.mp.dps = 50
TOL = mp.mpf("1e-40")

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for _p in (str(HERE), str(ROOT)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from relation import (  # noqa: E402
    SOLVER_NAME,
    VTX_PATH,
    _var,
    colouring_cnf,
    load_509,
    relation,
)
from pysat.solvers import Solver  # noqa: E402
from tools.verify_509 import build as build_unit_edges  # noqa: E402
from tools.verify_509 import parse_point as parse_point_509  # noqa: E402

OUT_PATH = ROOT / "runs" / "hn_509_relations.json"

CAP_SECONDS = 900.0
RANDOM_PAIRS = 2000
SEED = 20261102
CROSSCHECK_PAIRS = 200
PROGRESS_EVERY = 4000
K = 5


# --------------------------------------------------------------------------
# time budget
# --------------------------------------------------------------------------

class Budget:
    """Wall-clock cap for the whole run, checked between phases and queries."""

    def __init__(self, cap_seconds: float):
        if cap_seconds <= 0:
            raise ValueError("cap must be positive")
        self.cap = float(cap_seconds)
        self.t0 = time.perf_counter()

    def elapsed(self) -> float:
        return time.perf_counter() - self.t0

    def remaining(self) -> float:
        return max(0.0, self.cap - self.elapsed())

    def expired(self) -> bool:
        return self.remaining() <= 0.0


# --------------------------------------------------------------------------
# graph loading: existing loader + independent rebuild with tools/verify_509.py
# --------------------------------------------------------------------------

def load_graph() -> dict:
    t0 = time.perf_counter()
    points, edges, worst = load_509()
    load_seconds = time.perf_counter() - t0

    t0 = time.perf_counter()
    pts2 = []
    with open(VTX_PATH, encoding="utf-8") as f:
        for line in f:
            p = parse_point_509(line)
            if p is not None:
                pts2.append(p)
    edges2, worst2 = build_unit_edges(pts2, TOL)
    verify_seconds = time.perf_counter() - t0
    assert worst2 <= TOL

    e1 = {(min(u, v), max(u, v)) for u, v in edges}
    e2 = {(min(u, v), max(u, v)) for u, v in edges2}
    identical = e1 == e2 and len(points) == len(pts2)

    adj = [set() for _ in range(len(points))]
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)
    deg = [len(a) for a in adj]

    print(
        f"graph: n={len(points)} m={len(edges)} "
        f"worst|d^2-1|={mp.nstr(worst, 6)} "
        f"edge-set identical to tools/verify_509.py: {identical} "
        f"({load_seconds + verify_seconds:.1f}s)",
        flush=True,
    )
    return {
        "points": points,
        "edges": edges,
        "adj": adj,
        "degree": deg,
        "worst": worst,
        "identical_to_verify_tool": identical,
        "m_verify_tool": len(edges2),
        "load_seconds": load_seconds,
        "verify_seconds": verify_seconds,
    }


def squared_distance_pairs(points, target_sq, tol: "mp.mpf" = TOL):
    """Pairs with |d^2 - target_sq| <= tol.

    Same arithmetic, precision (50 digits) and tolerance as
    tools/verify_509.build(), with the target distance parameterised.
    """
    n = len(points)
    hits: list[tuple[int, int]] = []
    min_gap = None
    for i in range(n):
        xi, yi = points[i]
        for j in range(i + 1, n):
            xj, yj = points[j]
            dev = abs((xi - xj) ** 2 + (yi - yj) ** 2 - target_sq)
            if dev <= tol:
                hits.append((i, j))
            elif min_gap is None or dev < min_gap:
                min_gap = dev
    return hits, min_gap


def distance_two_pairs(adj) -> list[tuple[int, int]]:
    out: set[tuple[int, int]] = set()
    for u in range(len(adj)):
        au = adj[u]
        for v in au:
            for w in adj[v]:
                if w != u and w not in au:
                    out.add((min(u, w), max(u, w)))
    return sorted(out)


def random_pairs(n: int, count: int, seed: int) -> list[tuple[int, int]]:
    """`count` distinct pairs drawn uniformly without replacement."""
    total = n * (n - 1) // 2
    count = min(count, total)
    rng = random.Random(seed)
    out: list[tuple[int, int]] = []
    seen: set[tuple[int, int]] = set()
    while len(out) < count:
        u, v = rng.sample(range(n), 2)
        key = (u, v) if u < v else (v, u)
        if key not in seen:
            seen.add(key)
            out.append(key)
    return out


# --------------------------------------------------------------------------
# incremental prober: same CNF and same semantics as relation(), one solver
# --------------------------------------------------------------------------

class IncrementalProber:
    """Selector-literal prober, semantically identical to `relation()`.

    Base formula: colouring_cnf(n, edges, k), which for k >= 2 has no unit
    clauses.  For each pair (u, v) two selector variables s_same and s_diff are
    introduced; the clauses guarded by s_same force x[u][c] <-> x[v][c], the
    clauses guarded by s_diff forbid a common colour.  Without the assumption
    the selectors are unconstrained and every guarded clause is satisfied.

    Query semantics, identical to relation():
        s_same UNSAT                    -> "different"  (u != v is forced)
        s_same SAT and s_diff UNSAT     -> "equal"      (u == v is forced)
        both SAT                        -> "free"
    An UNSAT base means "unsat" for every pair.
    """

    def __init__(self, edges, n: int, k: int, pairs):
        self.n, self.k = n, k
        self.pairs = list(pairs)
        t0 = time.perf_counter()
        clauses = list(colouring_cnf(n, edges, k))
        next_var = n * k + 1
        self._index: dict[tuple[int, int], tuple[int, int]] = {}
        for u, v in self.pairs:
            s_same, s_diff = next_var, next_var + 1
            next_var += 2
            self._index[(u, v)] = (s_same, s_diff)
            for c in range(k):
                clauses.append([-s_same, -_var(u, c, k), _var(v, c, k)])
                clauses.append([-s_same, -_var(v, c, k), _var(u, c, k)])
                clauses.append([-s_diff, -_var(u, c, k), -_var(v, c, k)])
        self.clauses = clauses
        self.n_vars = next_var - 1
        self.build_seconds = time.perf_counter() - t0

        t0 = time.perf_counter()
        self.solver = Solver(name=SOLVER_NAME, bootstrap_with=clauses)
        self.base_sat = bool(self.solver.solve())
        self.base_seconds = time.perf_counter() - t0
        self.model = self.solver.get_model() if self.base_sat else None

        self.queries = 0
        self.solve_calls = 0
        self.seconds = 0.0
        self.max_query_seconds = 0.0
        self.slowest: list[tuple[float, int, int]] = []

    def close(self) -> None:
        self.solver.delete()

    def ask(self, u: int, v: int) -> str:
        if not self.base_sat:
            return "unsat"
        key = (u, v) if u < v else (v, u)
        s_same, s_diff = self._index[key]
        t0 = time.perf_counter()
        self.solve_calls += 1
        if not self.solver.solve(assumptions=[s_same]):
            out = "different"
        else:
            self.solve_calls += 1
            out = "equal" if not self.solver.solve(assumptions=[s_diff]) else "free"
        dt = time.perf_counter() - t0
        self.queries += 1
        self.seconds += dt
        if dt > self.max_query_seconds:
            self.max_query_seconds = dt
            self.slowest.append((dt, u, v))
            self.slowest = sorted(self.slowest, reverse=True)[:5]
        return out


def colouring_from_model(model, n: int, k: int):
    """Decode a one-hot model into a colour per vertex, or None if malformed."""
    if model is None:
        return None
    colours = []
    for v in range(n):
        chosen = [c for c in range(k) if model[v * k + c] > 0]
        if len(chosen) != 1:
            return None
        colours.append(chosen[0])
    return colours


def colouring_valid(colours, edges) -> bool:
    return colours is not None and all(colours[u] != colours[v] for u, v in edges)


# --------------------------------------------------------------------------
# forced-equality graph
# --------------------------------------------------------------------------

def equality_components(n: int, eq_pairs) -> dict:
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for u, v in eq_pairs:
        ru, rv = find(u), find(v)
        if ru != rv:
            parent[ru] = rv

    groups: dict[int, list[int]] = collections.defaultdict(list)
    for v in range(n):
        groups[find(v)].append(v)
    comps = sorted((sorted(g) for g in groups.values()), key=lambda g: (-len(g), g[0]))
    sizes = [len(g) for g in comps]
    return {
        "vertices": n,
        "edges": len(eq_pairs),
        "components": len(comps),
        "component_size_histogram": {str(s): sizes.count(s) for s in sorted(set(sizes))},
        "component_sizes_desc": sizes[:25],
        "max_component_size": max(sizes) if sizes else 0,
        "components_with_at_least_2": sum(1 for s in sizes if s >= 2),
        "components_larger_than_2": sum(1 for s in sizes if s > 2),
        "largest_nontrivial_component": max((s for s in sizes if s >= 2), default=0),
        "non_singleton_components": [g for g in comps if len(g) >= 2][:25],
        "components_full": comps,
    }


# --------------------------------------------------------------------------
# diamonds
# --------------------------------------------------------------------------

def find_diamonds(adj, sq3_pairs) -> dict:
    """4-vertex K4 minus an edge whose apices are at squared distance 3.

    If p, q are at distance sqrt(3), the points at distance 1 from both are the
    two circle intersections and are themselves at distance 1.  So a diamond is
    exactly a squared-distance-3 pair with two common neighbours.
    """
    diamonds: list[tuple[int, int, int, int]] = []
    per_apex: collections.Counter = collections.Counter()
    two_common = 0
    common_not_adjacent = 0
    max_common = 0
    for p, q in sq3_pairs:
        common = sorted(adj[p] & adj[q])
        max_common = max(max_common, len(common))
        if len(common) != 2:
            continue
        two_common += 1
        a, b = common
        if b not in adj[a]:
            common_not_adjacent += 1
            continue
        diamonds.append((p, q, a, b))
        per_apex[(p, q)] += 1
    return {
        "diamonds": diamonds,
        "per_apex": per_apex,
        "apex_pairs_with_two_common": two_common,
        "apex_pairs_two_common_not_adjacent": common_not_adjacent,
        "max_common_neighbours_of_sq3_pair": max_common,
    }


def verify_diamond(edge_set, p: int, q: int, a: int, b: int) -> bool:
    """True iff {p,q,a,b} spans exactly K4 minus the edge pq."""
    vs = sorted({p, q, a, b})
    if len(vs) != 4:
        return False
    present = 0
    missing = []
    for x, y in itertools.combinations(vs, 2):
        if (x, y) in edge_set:
            present += 1
        else:
            missing.append((x, y))
    return present == 5 and missing == [(p, q)]


def gadget_relation(p: int, q: int, a: int, b: int, k: int) -> str:
    """relation() on the isolated 4-vertex diamond, apices p and q."""
    vs = sorted({p, q, a, b})
    local = {v: i for i, v in enumerate(vs)}
    es = [
        (local[p], local[a]),
        (local[p], local[b]),
        (local[q], local[a]),
        (local[q], local[b]),
        (local[a], local[b]),
    ]
    return relation(es, 4, local[p], local[q], k)


# --------------------------------------------------------------------------
# main census
# --------------------------------------------------------------------------

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Forced colour relations in the 509-vertex Parts graph at k=5."
    )
    ap.add_argument("--cap", type=float, default=CAP_SECONDS,
                    help="wall-clock cap in seconds for the whole run")
    ap.add_argument("--random-pairs", type=int, default=RANDOM_PAIRS)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--crosscheck", type=int, default=CROSSCHECK_PAIRS)
    ap.add_argument("--no-direct-diamonds", action="store_true",
                    help="skip the per-diamond relation() recheck at k=5")
    ap.add_argument("--out", default=str(OUT_PATH))
    args = ap.parse_args(argv)

    budget = Budget(args.cap)
    timings: dict[str, float] = {}
    problems: list[str] = []

    print("== 509-vertex Parts graph: forced colour relations at k = 5 ==", flush=True)

    # ---- load ------------------------------------------------------------
    g = load_graph()
    points, edges, adj = g["points"], g["edges"], g["adj"]
    n = len(points)
    edge_set = {(min(u, v), max(u, v)) for u, v in edges}
    timings["load"] = g["load_seconds"] + g["verify_seconds"]
    if not g["identical_to_verify_tool"]:
        problems.append("edge set differs from tools/verify_509.py")

    if budget.expired():
        problems.append("cap hit during graph loading")

    # ---- strata ----------------------------------------------------------
    d2 = distance_two_pairs(adj)
    d2_set = set(d2)
    rnd = random_pairs(n, args.random_pairs, args.seed)
    rnd_new = [p for p in rnd if p not in edge_set and p not in d2_set]

    strata = {
        "edges": sorted(edge_set),
        "distance2": d2,
        "random": rnd,
    }
    pair_tags: dict[tuple[int, int], list[str]] = collections.defaultdict(list)
    for name, ps in strata.items():
        for p in ps:
            pair_tags[p].append(name)
    all_pairs = sorted(pair_tags)

    print(
        f"strata: edges={len(strata['edges'])} distance2={len(d2)} "
        f"random={len(rnd)} (random pairs not in the other strata: {len(rnd_new)}) "
        f"-> {len(all_pairs)} unique pairs",
        flush=True,
    )

    # ---- squared distance 3 pairs (diamond apices) ------------------------
    t0 = time.perf_counter()
    sq3, sq3_gap = squared_distance_pairs(points, mp.mpf(3))
    timings["sq3_search"] = time.perf_counter() - t0
    print(
        f"squared-distance-3 pairs: {len(sq3)} "
        f"(nearest non-match |d^2-3| = "
        f"{mp.nstr(sq3_gap, 6) if sq3_gap is not None else 'n/a'}, "
        f"{timings['sq3_search']:.1f}s)",
        flush=True,
    )

    # ---- k = 5 base + bulk scan -----------------------------------------
    prober = IncrementalProber(edges, n, K, all_pairs)
    timings["cnf_build"] = prober.build_seconds
    timings["base_solve"] = prober.base_seconds
    base = "sat" if prober.base_sat else "unsat"
    print(
        f"k=5: base={base.upper()} in {prober.base_seconds:.2f}s "
        f"({prober.n_vars} vars, {len(prober.clauses)} clauses)",
        flush=True,
    )
    if not prober.base_sat:
        problems.append("k=5 base is UNSAT; the graph is not 5-colourable")

    colours = colouring_from_model(prober.model, n, K)
    valid_colouring = colouring_valid(colours, edges)
    if not valid_colouring:
        problems.append("witness 5-colouring from the solver is not proper")
    colour_class_sizes = (
        [sum(1 for c in colours if c == i) for i in range(K)] if colours else []
    )

    relations: dict[tuple[int, int], str] = {}
    scan_capped = False
    t0 = time.perf_counter()
    for i, p in enumerate(all_pairs):
        if budget.expired():
            scan_capped = True
            break
        relations[p] = prober.ask(*p)
        if (i + 1) % PROGRESS_EVERY == 0:
            eq = sum(1 for x in relations.values() if x == "equal")
            print(
                f"  scan {i + 1}/{len(all_pairs)} "
                f"({time.perf_counter() - t0:.0f}s, equal so far: {eq})",
                flush=True,
            )
    timings["scan"] = time.perf_counter() - t0
    prober.close()

    if scan_capped:
        problems.append("cap hit during the pair scan; tallies are partial")

    # ---- tallies ---------------------------------------------------------
    def tally(pairs):
        c = collections.Counter(relations[p] for p in pairs if p in relations)
        return {
            "pairs": len(pairs),
            "pairs_examined": sum(c.values()),
            "equal": c.get("equal", 0),
            "different": c.get("different", 0),
            "free": c.get("free", 0),
            "unsat": c.get("unsat", 0),
        }

    tallies = {name: tally(ps) for name, ps in strata.items()}
    total_tally = tally(all_pairs)
    equal_pairs = [p for p in all_pairs if relations.get(p) == "equal"]
    nonedge_different = [
        p for p in all_pairs if relations.get(p) == "different" and p not in edge_set
    ]

    print(f"k=5 tally: {total_tally['equal']} equal, "
          f"{total_tally['different']} different, {total_tally['free']} free "
          f"({total_tally['pairs_examined']}/{len(all_pairs)} pairs examined)", flush=True)
    for name in ("edges", "distance2", "random"):
        t = tallies[name]
        print(
            f"  {name:<10} equal={t['equal']:<5} different={t['different']:<5} "
            f"free={t['free']:<6} unsat={t['unsat']}",
            flush=True,
        )
    print(f"  forced-different non-edges: {len(nonedge_different)}", flush=True)

    # ---- forced-equality graph ------------------------------------------
    t0 = time.perf_counter()
    eq_graph = equality_components(n, equal_pairs)
    timings["equality_graph"] = time.perf_counter() - t0
    print(
        f"forced-equality graph: {eq_graph['edges']} edges -> "
        f"{eq_graph['components']} components, sizes "
        f"{eq_graph['component_size_histogram']}, max={eq_graph['max_component_size']}, "
        f"components with >2 vertices: {eq_graph['components_larger_than_2']}",
        flush=True,
    )

    # ---- diamonds --------------------------------------------------------
    t0 = time.perf_counter()
    d = find_diamonds(adj, sq3)
    diamonds = d["diamonds"]
    print(
        f"diamonds (K4-e, apices at squared distance 3): {len(diamonds)} "
        f"[apex pairs with two common neighbours: {d['apex_pairs_with_two_common']}, "
        f"of which the common neighbours are non-adjacent: "
        f"{d['apex_pairs_two_common_not_adjacent']}]",
        flush=True,
    )

    verified = sum(1 for p, q, a, b in diamonds if verify_diamond(edge_set, p, q, a, b))
    if verified != len(diamonds):
        problems.append(f"{len(diamonds) - verified} diamonds are not K4 minus an edge")
    per_apex_max = max(d["per_apex"].values(), default=0)
    if per_apex_max > 1:
        problems.append("an apex pair supports more than one diamond")

    # local gadget relation at k=3 (and k=4 for contrast): 4 vertices, trivial
    t0 = time.perf_counter()
    local_k3: collections.Counter = collections.Counter()
    local_k4: collections.Counter = collections.Counter()
    for p, q, a, b in diamonds:
        local_k3[gadget_relation(p, q, a, b, 3)] += 1
        local_k4[gadget_relation(p, q, a, b, 4)] += 1
    timings["diamond_local_k3_k4"] = time.perf_counter() - t0
    if local_k3.get("equal", 0) != len(diamonds):
        problems.append("a diamond does not force its apices equal at k=3")

    # apices inside the full 509-vertex graph at k=5 (from the bulk scan)
    apex_pairs = [(p, q) for p, q, _a, _b in diamonds]
    apex_in_scan = [x for x in apex_pairs if x in relations]
    if len(apex_in_scan) != len(apex_pairs):
        problems.append("an apex pair was not covered by the distance-2 stratum")
    full_k5: collections.Counter = collections.Counter(
        relations[x] for x in apex_in_scan
    )

    # independent re-probe of every apex pair with relation() itself
    direct_k5: collections.Counter = collections.Counter()
    direct_mismatches: list[dict] = []
    direct_seconds = 0.0
    direct_skipped = None
    if args.no_direct_diamonds:
        direct_skipped = "--no-direct-diamonds"
    elif budget.remaining() < 60.0:
        direct_skipped = "cap"
    else:
        t0 = time.perf_counter()
        for p, q in apex_pairs:
            got = relation(edges, n, p, q, K)
            direct_k5[got] += 1
            if relations.get((p, q)) != got:
                direct_mismatches.append(
                    {"pair": [p, q], "bulk": relations.get((p, q)), "direct": got}
                )
        direct_seconds = time.perf_counter() - t0
        timings["diamond_direct_k5"] = direct_seconds
        if direct_mismatches:
            problems.append("relation() and the bulk prober disagree on a diamond apex")

    print(
        f"  diamond apices as a 4-vertex gadget: k=3 {dict(local_k3)}, "
        f"k=4 {dict(local_k4)}",
        flush=True,
    )
    print(
        f"  diamond apices in the full graph at k=5: {dict(full_k5)} "
        f"(direct relation() recheck: {dict(direct_k5) if not direct_skipped else 'skipped: ' + direct_skipped}, "
        f"{direct_seconds:.1f}s)",
        flush=True,
    )

    # ---- validation: direct relation() against the bulk prober ------------
    t0 = time.perf_counter()
    rng = random.Random(args.seed + 1)
    per_stratum = max(1, args.crosscheck // 3)
    check: list[tuple[int, int]] = []
    for name in ("edges", "distance2", "random"):
        ps = strata[name]
        take = min(per_stratum, len(ps))
        check.extend(rng.sample(ps, take) if ps else [])
    check = sorted(set(check))
    cross_mismatches: list[dict] = []
    for p in check:
        if budget.expired():
            problems.append("cap hit during the crosscheck")
            break
        got = relation(edges, n, p[0], p[1], K)
        if relations.get(p) != got:
            cross_mismatches.append(
                {"pair": list(p), "bulk": relations.get(p), "direct": got}
            )
    timings["crosscheck"] = time.perf_counter() - t0
    if cross_mismatches:
        problems.append(f"{len(cross_mismatches)} crosscheck mismatches")
    edges_all_different = all(relations.get(p) == "different" for p in edge_set)

    print(
        f"validation: {len(check)} pairs re-probed with relation(...,k=5), "
        f"mismatches={len(cross_mismatches)}; all edges different: {edges_all_different}",
        flush=True,
    )

    # ---- timing / cap ----------------------------------------------------
    timings["total"] = budget.elapsed()
    capped = budget.expired() or scan_capped
    print(
        f"timings: " + ", ".join(f"{k}={v:.1f}s" for k, v in timings.items())
        + f" | cap={budget.cap:.0f}s hit={capped} "
        f"max single query={1000 * prober.max_query_seconds:.1f} ms",
        flush=True,
    )

    # ---- JSON ------------------------------------------------------------
    payload = {
        "solver": SOLVER_NAME,
        "k": K,
        "cap_seconds": budget.cap,
        "capped": capped,
        "cap_hit": capped,
        "problems": problems,
        "timings_s": timings,
        "graph": {
            "path": str(VTX_PATH.relative_to(ROOT)).replace("\\", "/"),
            "n": n,
            "m": len(edges),
            "worst_abs_d2_minus_1": float(g["worst"]),
            "edge_set_matches_tools_verify_509": g["identical_to_verify_tool"],
            "m_from_tools_verify_509": g["m_verify_tool"],
            "degree_min": min(g["degree"]),
            "degree_max": max(g["degree"]),
            "degree_mean": sum(g["degree"]) / n,
        },
        "k5": {
            "base": base,
            "witness_colouring_proper": valid_colouring,
            "colours_used": len(set(colours)) if colours else 0,
            "colour_class_sizes": colour_class_sizes,
            "pairs_total_c509_2": n * (n - 1) // 2,
            "unique_pairs_examined": total_tally["pairs_examined"],
            "random_pairs_drawn": len(rnd),
            "random_pairs_not_in_other_strata": len(rnd_new),
            "random_seed": args.seed,
            "strata": tallies,
            "tally": total_tally,
            "coverage": {
                "edges": "exhaustive (all 2442 edges)",
                "distance2": "exhaustive (all pairs at graph distance 2)",
                "random": (
                    f"{len(rnd)} distinct pairs, uniform without replacement, "
                    f"seed {args.seed}"
                ),
                "pairs_not_examined": n * (n - 1) // 2 - total_tally["pairs_examined"],
            },
            "equal_pairs": [list(p) for p in equal_pairs],
            "nonedge_different_pairs": [list(p) for p in nonedge_different],
            "edges_all_different": edges_all_different,
        },
        "forced_equality_graph": eq_graph,
        "diamonds": {
            "definition": "4 vertices spanning K4 minus one edge, the non-adjacent "
                          "pair (apices) at squared distance exactly 3",
            "count": len(diamonds),
            "verified_as_k4_minus_edge": verified,
            "distinct_apex_pairs": len(d["per_apex"]),
            "max_diamonds_per_apex_pair": per_apex_max,
            "squared_distance_3_pairs": len(sq3),
            "nearest_nonmatch_abs_d2_minus_3": (
                float(sq3_gap) if sq3_gap is not None else None
            ),
            "apex_pairs_with_two_common_neighbours": d["apex_pairs_with_two_common"],
            "apex_pairs_two_common_not_adjacent": d[
                "apex_pairs_two_common_not_adjacent"
            ],
            "max_common_neighbours_of_squared_distance_3_pair": d[
                "max_common_neighbours_of_sq3_pair"
            ],
            "local_gadget_k3": dict(local_k3),
            "local_gadget_k4": dict(local_k4),
            "local_gadget_k3_all_equal": local_k3.get("equal", 0) == len(diamonds),
            "full_graph_k5": dict(full_k5),
            "full_graph_k5_direct_relation": (
                None if direct_skipped else dict(direct_k5)
            ),
            "full_graph_k5_direct_seconds": direct_seconds,
            "full_graph_k5_direct_skipped": direct_skipped,
            "full_graph_k5_direct_mismatches": direct_mismatches,
        },
        "validation": {
            "edge_set_matches_tools_verify_509": g["identical_to_verify_tool"],
            "witness_colouring_proper": valid_colouring,
            "edges_all_different": edges_all_different,
            "diamonds_all_k4_minus_edge": verified == len(diamonds),
            "diamonds_all_force_equal_at_k3": local_k3.get("equal", 0) == len(diamonds),
            "crosscheck_with_relation": {
                "pairs": len(check),
                "mismatches": cross_mismatches,
                "seconds": timings["crosscheck"],
            },
            "bulk_prober": {
                "vars": prober.n_vars,
                "clauses": len(prober.clauses),
                "base_seconds": prober.base_seconds,
                "queries": prober.queries,
                "solve_calls": prober.solve_calls,
                "scan_seconds": timings["scan"],
                "max_query_seconds": prober.max_query_seconds,
                "slowest_queries": [
                    {"seconds": s, "pair": [u, v]} for s, u, v in prober.slowest
                ],
            },
        },
        "notes": [
            "k=3 and k=4 are UNSAT for the full 509 graph, so a relation() call "
            "there returns 'unsat' for every pair; the diamond k=3/k=4 checks are "
            "therefore made on the isolated 4-vertex gadget.",
            "relation 'equal' is transitive in every colouring, so a connected "
            "component of the forced-equality graph is a clique of forced "
            "equalities; a component of size >= 3 is a forced-equality chain.",
            "For a non-adjacent pair, 'equal' means adding the edge uv makes the "
            "graph 6-chromatic, and 'different' means identifying u and v makes "
            "it 6-chromatic.",
            "Coverage is exhaustive for edges and for graph distance 2 and random "
            "for the rest; forced relations between pairs at graph distance >= 3 "
            "are not excluded by this census, they are simply not sampled "
            "exhaustively.",
        ],
    }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"written: {out}", flush=True)

    if problems:
        for msg in problems:
            print(f"FAIL: {msg}", flush=True)
        return 1
    print("all internal checks passed", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
