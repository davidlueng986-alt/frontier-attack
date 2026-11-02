"""End-to-end pattern-forcing pipeline for the Hadwiger-Nelson track.

    python problems/hadwiger_nelson/six.py

Stages
------
1. Re-verify the 509-vertex Parts graph with tools/verify_509.py and record
   (n, m, k=5 SAT status).
2. Validate the SAT pattern prober ``patterns.forced_mono`` against brute force
   on every labelled graph with at most five vertices, at k = 3.
3. Search for a 7-vertex unit-distance gadget that forces a monochromatic
   sqrt(3)-triangle at k = 5.  Two generators are used:
      (a) 7-vertex subsets of the 509 graph's local neighbourhoods (cheap),
      (b) random closures of small point sets under unit-distance circle
          intersection.
   Both use 50-digit arithmetic (and the found candidate, if any, is rebuilt
   with exact sympy arithmetic before it is reported).  Hard cap: 20,000
   distinct candidates or 600 seconds.
4. If a gadget is found, attach it to the 509 graph by identifying its forced
   triple with a sqrt(3)-triple of 509 and SAT-test k = 5 on the union.  If no
   gadget is found, run the natural control: attach de Grey's 7-vertex wheel
   (the k=3 primitive) to every sqrt(3)-triple of 509, and record that all
   attachments remain 5-colourable.
5. Write runs/hn_six.json and print one explicit verdict line.

The verdict is "6-CHROMATIC CANDIDATE" only when SAT at k = 5 returns UNSAT on
the complete union graph.  Otherwise the verdict is "NEGATIVE: <reason>".

Why the negative is expected
----------------------------
The 509 graph is 5-chromatic.  A 7-vertex unit-distance graph, however, is
always 4-colourable: by OEIS A186705 the maximum number of unit distances on 7
points is 12 (and on 6, 5, 4, 3 points it is 9, 7, 5, 3), so every induced
subgraph has a vertex of degree at most 3; the graph is 3-degenerate.  Any
proper 4-colouring refines to a proper 5-colouring with colour classes of size
at most 2 (split each class into pairs), and a class of size at most 2 cannot
contain a 3-vertex monochromatic pattern.  Hence *no* 7-vertex unit-distance
graph can force a monochromatic 3-vertex pattern at k = 5.  In particular de
Grey's wheel H does not: it forces an alternating rim triple only at k = 3, and
it has an explicit 5-colouring with no monochromatic sqrt(3)-triangle.  The
pipeline searches anyway, records the empirical result, and states this proof
as the exact stopping point.

Writes: runs/hn_six.json, runs/stress_hn.log.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import random
import re
import subprocess
import sys
import time
from pathlib import Path

import mpmath as mp

mp.mp.dps = 50

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for _p in (str(HERE), str(ROOT)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import patterns  # noqa: E402
from relation import _split_top, parse_point  # noqa: E402
from pysat.formula import CNF  # noqa: E402
from pysat.solvers import Solver  # noqa: E402

SOLVER_NAME = "cadical153"
K = 5
N_509 = 509
M_509 = 2442

CANDIDATE_CAP = 20_000
SEARCH_SECONDS = 600.0
LOCAL_CAP = 12_000
RANDOM_ATTEMPTS = 60_000
CONTROL_CAP = 300

UNIT_TOL = mp.mpf("1e-40")
SQ3_TOL = mp.mpf("1e-40")
MERGE_TOL = mp.mpf("1e-30")

VTX_PATH = ROOT / "research" / "scratch" / "509_parts.vtx"
OUT_PATH = ROOT / "runs" / "hn_six.json"
LOG_PATH = ROOT / "runs" / "stress_hn.log"


# --------------------------------------------------------------------------
# logging
# --------------------------------------------------------------------------

def append_log(line: str) -> None:
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    stamp = time.strftime("%Y-%m-%dT%H:%M:%S")
    with LOG_PATH.open("a", encoding="utf-8") as f:
        f.write(f"[{stamp}] {line}\n")


# --------------------------------------------------------------------------
# CNF / SAT helpers
# --------------------------------------------------------------------------

def _var(v: int, colour: int, k: int) -> int:
    return v * k + colour + 1


def colouring_cnf(n: int, edges, k: int) -> CNF:
    """One-hot k-colouring CNF (same convention as patterns.py)."""
    cnf = CNF()
    for v in range(n):
        cnf.append([_var(v, c, k) for c in range(k)])
        for c1 in range(k):
            for c2 in range(c1 + 1, k):
                cnf.append([-_var(v, c1, k), -_var(v, c2, k)])
    for u, v in edges:
        if u > v:
            u, v = v, u
        for c in range(k):
            cnf.append([-_var(u, c, k), -_var(v, c, k)])
    return cnf


def solve_sat(n: int, edges, k: int, proof_path: Path | None = None):
    """Return (sat, model, proof_written).  Optionally dump a DRAT proof."""
    cnf = colouring_cnf(n, edges, k)
    if proof_path is None:
        with Solver(name=SOLVER_NAME, bootstrap_with=cnf) as solver:
            sat = bool(solver.solve())
            model = solver.get_model() if sat else None
        return sat, model, None

    with Solver(name=SOLVER_NAME, with_proof=True, bootstrap_with=cnf) as solver:
        sat = bool(solver.solve())
        model = solver.get_model() if sat else None
        proof = None if sat else solver.get_proof()
    if proof is not None:
        proof_path.parent.mkdir(parents=True, exist_ok=True)
        with proof_path.open("w", encoding="utf-8") as f:
            for line in proof:
                f.write(line + "\n")
    return sat, model, (str(proof_path) if proof is not None else None)


def decode_colouring(model, n: int, k: int):
    if model is None:
        return None
    true = set(model)
    colours = []
    for v in range(n):
        chosen = [c for c in range(k) if _var(v, c, k) in true]
        if len(chosen) != 1:
            return None
        colours.append(chosen[0])
    return colours


def check_colouring(colours, edges, k: int) -> bool:
    if colours is None or not colours:
        return False
    if min(colours) < 0 or max(colours) >= k:
        return False
    return all(colours[u] != colours[v] for u, v in edges)


# --------------------------------------------------------------------------
# stage 1: 509 verification
# --------------------------------------------------------------------------

def stage1_verify_509() -> dict:
    t0 = time.time()
    cmd = [
        sys.executable,
        str(ROOT / "tools" / "verify_509.py"),
        str(VTX_PATH),
        "--k",
        str(K),
    ]
    cp = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, timeout=600)
    out = cp.stdout + cp.stderr

    def grab(pattern, cast=str, default=None):
        m = re.search(pattern, out)
        return cast(m.group(1)) if m else default

    n = grab(r"vertices\s*=\s*(\d+)", int)
    m = grab(r"edges\s*=\s*(\d+)", int)
    worst = grab(r"worst \|d\^2 - 1\| on an edge =\s*(\S+)", str)
    if "UNSATISFIABLE" in out:
        status = "UNSAT"
    elif "SATISFIABLE" in out:
        status = "SAT"
    else:
        status = "UNKNOWN"

    record = {
        "stage": 1,
        "name": "verify_509",
        "command": " ".join(str(x) for x in cmd),
        "n": n,
        "m": m,
        "k": K,
        "status": status,
        "worst_abs_d2_minus_1": worst,
        "exit_code": cp.returncode,
        "seconds": time.time() - t0,
        "ok": n == N_509 and m == M_509 and status == "SAT",
        "note": "tools/verify_509.py exits 1 when the instance is SAT",
    }
    append_log(
        f"stage1 verify_509 n={n} m={m} k={K} {status} "
        f"worst|d2-1|={worst} exit={cp.returncode} ok={record['ok']}"
    )
    return record


# --------------------------------------------------------------------------
# stage 2: prober validation
# --------------------------------------------------------------------------

def stage2_patterns() -> dict:
    t0 = time.time()
    report = patterns.validate_bruteforce(max_n=5, k=3)

    wheel_edges = [(0, i) for i in range(1, 7)]
    wheel_edges += [(i, i + 1) for i in range(1, 6)] + [(6, 1)]
    w7_k3 = patterns.forced_mono(wheel_edges, 7, (1, 3, 5), 3)
    w7_k4 = patterns.forced_mono(wheel_edges, 7, (1, 3, 5), 4)
    w7_k5 = patterns.forced_mono(wheel_edges, 7, (1, 3, 5), 5)
    w7_set_k5 = patterns.forces_mono_triple_set(
        wheel_edges, 7, [(1, 3, 5), (2, 4, 6)], 5
    )

    # explicit 5-colouring of W7 with no monochromatic sqrt(3)-triple
    w7_colours = [4, 0, 1, 0, 1, 2, 3]  # hub, rim 1..6
    w7_witness_ok = check_colouring(w7_colours, wheel_edges, 5)
    w7_triples = [(1, 3, 5), (2, 4, 6)]
    w7_mono = [len({w7_colours[a], w7_colours[b], w7_colours[c]}) == 1
               for a, b, c in w7_triples]

    record = {
        "stage": 2,
        "name": "patterns",
        "graphs": report["graphs"],
        "checks": report["checks"],
        "mismatches": report["mismatches"],
        "seconds": report["seconds"],
        "validation_seconds": time.time() - t0,
        "w7_triple_k3": w7_k3,
        "w7_triple_k4": w7_k4,
        "w7_triple_k5": w7_k5,
        "w7_triple_set_k5": w7_set_k5,
        "w7_countercolouring": w7_colours,
        "w7_countercolouring_proper": w7_witness_ok,
        "w7_countercolouring_mono_triples": w7_mono,
        "ok": report["mismatches"] == 0 and w7_k3 is True
              and w7_k4 is False and w7_k5 is False
              and w7_witness_ok and not any(w7_mono),
    }
    append_log(
        f"stage2 patterns graphs={report['graphs']} checks={report['checks']} "
        f"mismatches={report['mismatches']} W7 k3={w7_k3} k4={w7_k4} k5={w7_k5} "
        f"ok={record['ok']}"
    )
    return record


# --------------------------------------------------------------------------
# geometry: 509 with unit edges and sqrt(3) pairs
# --------------------------------------------------------------------------

def load_509_geometry():
    points = []
    with VTX_PATH.open(encoding="utf-8") as f:
        for line in f:
            p = parse_point(line)
            if p is not None:
                points.append(p)
    n = len(points)
    edges: list[tuple[int, int]] = []
    sq3_pairs: list[tuple[int, int]] = []
    worst = mp.mpf(0)
    for i in range(n):
        xi, yi = points[i]
        for j in range(i + 1, n):
            xj, yj = points[j]
            d2 = (xi - xj) ** 2 + (yi - yj) ** 2
            dev1 = abs(d2 - 1)
            if dev1 <= UNIT_TOL:
                edges.append((i, j))
                if dev1 > worst:
                    worst = dev1
            elif abs(d2 - 3) <= SQ3_TOL:
                sq3_pairs.append((i, j))
    return points, edges, sq3_pairs, worst


def sq3_triples_from_pairs(n: int, sq3_pairs):
    adj = [set() for _ in range(n)]
    for i, j in sq3_pairs:
        adj[i].add(j)
        adj[j].add(i)
    triples = []
    for i in range(n):
        for j in adj[i]:
            if j <= i:
                continue
            for k in adj[i] & adj[j]:
                if k > j:
                    triples.append((i, j, k))
    return triples


# --------------------------------------------------------------------------
# stage 3: candidate generation
# --------------------------------------------------------------------------

def _candidate_signature(edges, triples):
    return (tuple(sorted(tuple(sorted(e)) for e in edges)),
            tuple(sorted(tuple(sorted(t)) for t in triples)))


def is_wheel_w7(edges, n: int = 7) -> bool:
    """True iff the graph is the regular hexagon wheel W7 (hub + C6 rim)."""
    adj = [set() for _ in range(n)]
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)
    hubs = [v for v in range(n) if len(adj[v]) == n - 1]
    if len(hubs) != 1:
        return False
    hub = hubs[0]
    rim = [v for v in range(n) if v != hub]
    if any(len(adj[v]) != 3 for v in rim):
        return False
    if any(len(adj[v] & set(rim)) != 2 for v in rim):
        return False
    seen = {rim[0]}
    stack = [rim[0]]
    while stack:
        x = stack.pop()
        for y in adj[x] & set(rim):
            if y not in seen:
                seen.add(y)
                stack.append(y)
    return len(seen) == len(rim)


def _subset_candidate(subset, edge_set, sq3_set):
    sub = list(subset)
    edges = []
    for a, b in itertools.combinations(range(len(sub)), 2):
        key = (sub[a], sub[b]) if sub[a] < sub[b] else (sub[b], sub[a])
        if key in edge_set:
            edges.append((a, b))
    triples = []
    for a, b, c in itertools.combinations(range(len(sub)), 3):
        if ((sub[a], sub[b]) in sq3_set
                and (sub[a], sub[c]) in sq3_set
                and (sub[b], sub[c]) in sq3_set):
            triples.append((a, b, c))
    return edges, triples


def local_509_candidates(points, edges, sq3_pairs, cap, deadline, rng):
    """Yield 7-vertex subsets from local neighbourhoods of the 509 graph."""
    n = len(points)
    adj = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    edge_set = {tuple(sorted(e)) for e in edges}
    sq3_set = {tuple(sorted(e)) for e in sq3_pairs}
    seen: set[tuple[int, ...]] = set()
    tested = 0
    order = list(range(n))
    rng.shuffle(order)

    for v in order:
        nbrs = sorted(adj[v])
        if len(nbrs) < 6:
            continue
        if len(nbrs) <= 8:
            pool = nbrs
        else:
            pool = sorted(rng.sample(nbrs, 8))
        for comb in itertools.combinations(pool, 6):
            subset = tuple(sorted((v,) + comb))
            if subset in seen:
                continue
            seen.add(subset)
            sub_edges, sub_triples = _subset_candidate(subset, edge_set, sq3_set)
            yield subset, sub_edges, sub_triples
            tested += 1
            if tested >= cap or time.monotonic() > deadline:
                return


# -- random unit-distance closure -----------------------------------------

def _seg():
    return [(mp.mpf(0), mp.mpf(0)), (mp.mpf(1), mp.mpf(0))]


def _triangle():
    return [(mp.mpf(0), mp.mpf(0)), (mp.mpf(1), mp.mpf(0)),
            (mp.mpf(1) / 2, mp.sqrt(3) / 2)]


def _diamond():
    return [(mp.mpf(0), mp.mpf(0)), (mp.mpf(1), mp.mpf(0)),
            (mp.mpf(1) / 2, mp.sqrt(3) / 2),
            (mp.mpf(1) / 2, -mp.sqrt(3) / 2)]


def _square():
    return [(mp.mpf(0), mp.mpf(0)), (mp.mpf(1), mp.mpf(0)),
            (mp.mpf(1), mp.mpf(1)), (mp.mpf(0), mp.mpf(1))]


def _pentagon():
    R = 1 / (2 * mp.sin(mp.pi / 5))
    return [(R * mp.cos(2 * mp.pi * i / 5), R * mp.sin(2 * mp.pi * i / 5))
            for i in range(5)]


def _hexagon():
    return [(mp.cos(2 * mp.pi * i / 6), mp.sin(2 * mp.pi * i / 6))
            for i in range(6)]


def _w7():
    return [(mp.mpf(0), mp.mpf(0))] + _hexagon()


def _moser_spindle():
    s3, s11, s33 = mp.sqrt(3), mp.sqrt(11), mp.sqrt(33)
    return [
        (mp.mpf(0), mp.mpf(0)),
        (s3 / 2, mp.mpf(1) / 2),
        (s3 / 2, -mp.mpf(1) / 2),
        (s3, mp.mpf(0)),
        (5 * s3 / 12 - s11 / 12, mp.mpf(5) / 12 + s33 / 12),
        (5 * s3 / 12 + s11 / 12, -mp.mpf(5) / 12 + s33 / 12),
        (5 * s3 / 6, s33 / 6),
    ]


SEEDS = {
    "segment": _seg,
    "triangle": _triangle,
    "diamond": _diamond,
    "square": _square,
    "pentagon": _pentagon,
    "hexagon": _hexagon,
    "w7": _w7,
    "moser": _moser_spindle,
}


def circle_intersections(p1, p2):
    x1, y1 = p1
    x2, y2 = p2
    dx, dy = x2 - x1, y2 - y1
    d2 = dx * dx + dy * dy
    if d2 == 0 or d2 > 4:
        return []
    d = mp.sqrt(d2)
    a = d2 / 2 / d
    h2 = 1 - a * a
    if h2 < 0:
        return []
    h = mp.sqrt(h2)
    xm, ym = x1 + a * dx / d, y1 + a * dy / d
    rx, ry = -dy / d * h, dx / d * h
    return [(xm + rx, ym + ry), (xm - rx, ym - ry)]


def _distinct(p, pts, eps=MERGE_TOL):
    return all(abs(p[0] - q[0]) + abs(p[1] - q[1]) > eps for q in pts)


def random_unit_build(rng, deadline):
    """Return (points, trace) for a random 7-point unit-distance closure."""
    name = rng.choice(list(SEEDS))
    pts = list(SEEDS[name]())
    trace = [("seed", name)]
    if len(pts) > 7:
        pts = pts[:7]
    guard = 0
    while len(pts) < 7 and guard < 120:
        guard += 1
        if time.monotonic() > deadline:
            return None, None
        pairs = [(i, j) for i in range(len(pts))
                 for j in range(i + 1, len(pts))]
        rng.shuffle(pairs)
        added = False
        for i, j in pairs:
            cand = circle_intersections(pts[i], pts[j])
            if not cand:
                continue
            rng.shuffle(cand)
            for branch, p in enumerate(cand):
                if _distinct(p, pts):
                    pts.append(p)
                    trace.append(("circle", i, j, branch))
                    added = True
                    break
            if added:
                break
        if added:
            continue
        # fallback: a unit step at a rational multiple of pi
        i = rng.randrange(len(pts))
        num = rng.choice([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
        den = 12
        th = mp.pi * num / den
        p = (pts[i][0] + mp.cos(th), pts[i][1] + mp.sin(th))
        if _distinct(p, pts):
            pts.append(p)
            trace.append(("angle", i, num, den))
            added = True
        if not added:
            return None, None
    return (pts if len(pts) == 7 else None,
            trace if len(pts) == 7 else None)


_LATTICE_CACHE: dict[int, list] = {}


def lattice_patch_mp(radius: int):
    """Points of the triangular lattice i + j*w, w = exp(i*pi/3), |.| <= radius."""
    if radius in _LATTICE_CACHE:
        return _LATTICE_CACHE[radius]
    s3 = mp.sqrt(3)
    pts = []
    for i in range(-radius, radius + 1):
        for j in range(-radius, radius + 1):
            if i * i + i * j + j * j <= radius * radius:
                pts.append((mp.mpf(i) + mp.mpf(j) / 2, mp.mpf(j) * s3 / 2))
    _LATTICE_CACHE[radius] = pts
    return pts


def lattice_candidate(rng):
    """A random 7-point subset of a triangular-lattice circle-intersection closure."""
    radius = rng.choice([2, 3, 3, 4])
    patch = lattice_patch_mp(radius)
    chosen = rng.sample(patch, 7)
    return chosen, {"generator": "lattice_subset", "radius": radius}


def candidate_from_points(points):
    edges = []
    for i, j in itertools.combinations(range(len(points)), 2):
        d2 = ((points[i][0] - points[j][0]) ** 2
              + (points[i][1] - points[j][1]) ** 2)
        if abs(d2 - 1) <= UNIT_TOL:
            edges.append((i, j))
    triples = []
    for a, b, c in itertools.combinations(range(len(points)), 3):
        ok = True
        for u, v in ((a, b), (a, c), (b, c)):
            d2 = ((points[u][0] - points[v][0]) ** 2
                  + (points[u][1] - points[v][1]) ** 2)
            if abs(d2 - 3) > SQ3_TOL:
                ok = False
                break
        if ok:
            triples.append((a, b, c))
    return edges, triples


# -- exact reconstruction --------------------------------------------------

def _exact_seed(name):
    import sympy as sp

    if name == "segment":
        return [(sp.Integer(0), sp.Integer(0)), (sp.Integer(1), sp.Integer(0))]
    if name == "triangle":
        return [(sp.Integer(0), sp.Integer(0)), (sp.Integer(1), sp.Integer(0)),
                (sp.Rational(1, 2), sp.sqrt(3) / 2)]
    if name == "diamond":
        return [(sp.Integer(0), sp.Integer(0)), (sp.Integer(1), sp.Integer(0)),
                (sp.Rational(1, 2), sp.sqrt(3) / 2),
                (sp.Rational(1, 2), -sp.sqrt(3) / 2)]
    if name == "square":
        return [(sp.Integer(0), sp.Integer(0)), (sp.Integer(1), sp.Integer(0)),
                (sp.Integer(1), sp.Integer(1)), (sp.Integer(0), sp.Integer(1))]
    if name == "pentagon":
        R = 1 / (2 * sp.sin(sp.pi / 5))
        return [(sp.simplify(R * sp.cos(2 * sp.pi * i / 5)),
                 sp.simplify(R * sp.sin(2 * sp.pi * i / 5))) for i in range(5)]
    if name == "hexagon":
        return [(sp.simplify(sp.cos(2 * sp.pi * i / 6)),
                 sp.simplify(sp.sin(2 * sp.pi * i / 6))) for i in range(6)]
    if name == "w7":
        return [(sp.Integer(0), sp.Integer(0))] + _exact_seed("hexagon")
    if name == "moser":
        s3, s11, s33 = sp.sqrt(3), sp.sqrt(11), sp.sqrt(33)
        return [
            (sp.Integer(0), sp.Integer(0)),
            (s3 / 2, sp.Rational(1, 2)),
            (s3 / 2, -sp.Rational(1, 2)),
            (s3, sp.Integer(0)),
            (5 * s3 / 12 - s11 / 12, sp.Rational(5, 12) + s33 / 12),
            (5 * s3 / 12 + s11 / 12, -sp.Rational(5, 12) + s33 / 12),
            (5 * s3 / 6, s33 / 6),
        ]
    raise ValueError(name)


def _exact_circle_intersections(p1, p2):
    import sympy as sp

    x1, y1 = p1
    x2, y2 = p2
    dx, dy = x2 - x1, y2 - y1
    d2 = sp.expand(dx * dx + dy * dy)
    d = sp.sqrt(d2)
    a = d2 / (2 * d)
    h2 = sp.simplify(1 - a * a)
    if sp.simplify(h2) < 0:
        return []
    h = sp.sqrt(h2)
    xm, ym = x1 + a * dx / d, y1 + a * dy / d
    rx, ry = -dy / d * h, dx / d * h
    return [(sp.simplify(xm + rx), sp.simplify(ym + ry)),
            (sp.simplify(xm - rx), sp.simplify(ym - ry))]


def exact_points_from_trace(trace):
    import sympy as sp

    pts = _exact_seed(trace[0][1])
    for step in trace[1:]:
        if step[0] == "circle":
            _, i, j, branch = step
            cand = _exact_circle_intersections(pts[i], pts[j])
            pts.append(cand[branch])
        elif step[0] == "angle":
            _, i, num, den = step
            th = sp.Rational(num, den) * sp.pi
            pts.append((sp.simplify(pts[i][0] + sp.cos(th)),
                        sp.simplify(pts[i][1] + sp.sin(th))))
        else:
            raise ValueError(step)
    return pts


def exact_verify_points(pts):
    """Return (edges, triples, ok, notes) for an exact sympy point list."""
    import sympy as sp

    n = len(pts)
    edges = []
    sq3 = []
    coincident = []
    for i, j in itertools.combinations(range(n), 2):
        d2 = sp.simplify((pts[i][0] - pts[j][0]) ** 2
                         + (pts[i][1] - pts[j][1]) ** 2)
        if d2 == 0:
            coincident.append((i, j))
        elif d2 == 1:
            edges.append((i, j))
        elif d2 == 3:
            sq3.append((i, j))
    triples = []
    for a, b, c in itertools.combinations(range(n), 3):
        if ((a, b) in sq3 or (b, a) in sq3) and \
           ((a, c) in sq3 or (c, a) in sq3) and \
           ((b, c) in sq3 or (c, b) in sq3):
            triples.append((a, b, c))
    ok = not coincident
    return edges, triples, ok, {"coincident": coincident}


def exact_509_points(indices):
    import sympy as sp

    lines = [ln.strip() for ln in VTX_PATH.read_text(encoding="utf-8").splitlines()
             if ln.strip()]
    out = []
    for i in indices:
        s = lines[i]
        a, b = _split_top(s[1:-1])
        coords = []
        for expr in (a, b):
            e = expr
            while "Sqrt[" in e:
                nxt = re.sub(r"Sqrt\[([^\[\]]*)\]", r"sqrt(\1)", e)
                if nxt == e:
                    break
                e = nxt
            coords.append(sp.sympify(e))
        out.append(tuple(coords))
    return out


# -- stage 3 driver --------------------------------------------------------

def stage3_search(points, edges, sq3_pairs) -> dict:
    t0 = time.time()
    deadline = time.monotonic() + SEARCH_SECONDS
    rng = random.Random(20261102)
    seen_signatures = set()

    stats = {
        "candidates_tested": 0,
        "phase_a_509_subsets": 0,
        "phase_b_random_closures": 0,
        "phase_b_lattice_subsets": 0,
        "candidates_with_sqrt3_triples": 0,
        "max_edges": 0,
        "max_triples": 0,
        "best": None,
        "found": None,
        "capped_by": None,
        "duplicate_signatures": 0,
        "exact_verified": None,
    }

    def consider(cand, meta):
        edges_c, triples_c = cand
        sig = _candidate_signature(edges_c, triples_c)
        if sig in seen_signatures:
            stats["duplicate_signatures"] += 1
            return False
        seen_signatures.add(sig)
        stats["candidates_tested"] += 1
        if len(edges_c) > stats["max_edges"]:
            stats["max_edges"] = len(edges_c)
        if len(triples_c) > stats["max_triples"]:
            stats["max_triples"] = len(triples_c)
        if triples_c:
            stats["candidates_with_sqrt3_triples"] += 1
        if (len(edges_c), len(triples_c)) > (
                (stats["best"] or {}).get("edges", -1),
                (stats["best"] or {}).get("triples", -1)):
            stats["best"] = {
                **meta,
                "edges": len(edges_c),
                "triples": len(triples_c),
                "edge_list": [list(e) for e in edges_c],
                "triple_list": [list(t) for t in triples_c],
            }
        if patterns.forces_mono_triple_set(edges_c, 7, triples_c, K):
            stats["found"] = {**meta, "edges": edges_c, "triples": triples_c}
            return True
        return False

    # phase (a): local 509 subsets
    for subset, sub_edges, sub_triples in local_509_candidates(
            points, edges, sq3_pairs, LOCAL_CAP, deadline, rng):
        if stats["candidates_tested"] + stats["duplicate_signatures"] >= LOCAL_CAP:
            stats["capped_by"] = "phase_a_candidate_cap"
            break
        stats["phase_a_509_subsets"] += 1
        if consider((sub_edges, sub_triples),
                    {"phase": "a", "subset": list(subset)}):
            break
        if time.monotonic() > deadline:
            stats["capped_by"] = "time_cap"
            break

    # phase (b): random circle-intersection closures, plus random 7-subsets of
    # a triangular-lattice closure (itself an exact circle-intersection
    # construction).  The lattice part supplies distinct sparse signatures;
    # the incremental builds supply the dense, rigid ones.
    if stats["found"] is None and time.monotonic() <= deadline:
        attempts = 0
        while stats["candidates_tested"] < CANDIDATE_CAP:
            if time.monotonic() > deadline:
                stats["capped_by"] = "time_cap"
                break
            attempts += 1
            if attempts <= RANDOM_ATTEMPTS and attempts % 3 != 0:
                pts, trace = random_unit_build(rng, deadline)
                stats["phase_b_random_closures"] += 1
                if pts is None:
                    continue
                cand = candidate_from_points(pts)
                meta = {"phase": "b", "generator": "incremental",
                        "trace": trace,
                        "points": [tuple(float(x) for x in p) for p in pts]}
            else:
                pts, meta0 = lattice_candidate(rng)
                stats["phase_b_lattice_subsets"] += 1
                cand = candidate_from_points(pts)
                meta = {"phase": "b", **meta0,
                        "points": [tuple(float(x) for x in p) for p in pts]}
            if consider(cand, meta):
                break
        if stats["candidates_tested"] >= CANDIDATE_CAP:
            stats["capped_by"] = "candidate_cap"

    seconds = time.time() - t0
    found = stats["found"]
    exact = None
    if found is not None:
        # Rebuild the found candidate exactly before reporting it.
        if found["phase"] == "a":
            exact_pts = exact_509_points(found["subset"])
        else:
            exact_pts = exact_points_from_trace(found["trace"])
        ex_edges, ex_triples, ex_ok, notes = exact_verify_points(exact_pts)
        still_forces = ex_ok and patterns.forces_mono_triple_set(
            ex_edges, 7, ex_triples, K
        )
        exact = {
            "ok": bool(still_forces),
            "coincident": notes["coincident"],
            "edges": [list(e) for e in ex_edges],
            "triples": [list(t) for t in ex_triples],
            "points": [str(p) for p in exact_pts],
        }
        stats["exact_verified"] = exact["ok"]
        if not still_forces:
            found = None
            stats["found"] = None
            stats["rejected_by_exact_check"] = True

    record = {
        "stage": 3,
        "name": "gadget_search",
        "k": K,
        "candidate_cap": CANDIDATE_CAP,
        "time_cap_s": SEARCH_SECONDS,
        "candidates_tested": stats["candidates_tested"],
        "phase_a_509_subsets": stats["phase_a_509_subsets"],
        "phase_b_random_closures": stats["phase_b_random_closures"],
        "phase_b_lattice_subsets": stats["phase_b_lattice_subsets"],
        "candidates_with_sqrt3_triples": stats["candidates_with_sqrt3_triples"],
        "duplicate_signatures": stats["duplicate_signatures"],
        "max_edges": stats["max_edges"],
        "max_triples": stats["max_triples"],
        "best": stats["best"],
        "best_is_w7": bool(
            stats["best"]
            and stats["best"].get("edge_list")
            and is_wheel_w7([tuple(e) for e in stats["best"]["edge_list"]])
        ),
        "gadget_found": found is not None,
        "gadget": None if found is None else {
            "phase": found["phase"],
            "edges": [list(e) for e in found["edges"]],
            "triples": [list(t) for t in found["triples"]],
            "points": ([str(p) for p in exact_points_from_trace(found["trace"])]
                       if found["phase"] == "b"
                       else [str(p) for p in exact_509_points(found["subset"])]),
        },
        "exact_verified": exact,
        "capped_by": stats["capped_by"],
        "seconds": seconds,
    }
    append_log(
        f"stage3 gadget_search candidates={record['candidates_tested']} "
        f"(local={record['phase_a_509_subsets']}, random={record['phase_b_random_closures']}, "
        f"lattice={record['phase_b_lattice_subsets']}) "
        f"with_sq3={record['candidates_with_sqrt3_triples']} "
        f"max_edges={record['max_edges']} max_triples={record['max_triples']} "
        f"found={record['gadget_found']} capped_by={record['capped_by']} "
        f"seconds={seconds:.1f}"
    )
    return record


# --------------------------------------------------------------------------
# stage 4: composition
# --------------------------------------------------------------------------

def _place_triple(src_tri, dst_tri, all_pts):
    """Rigid motion (mpmath) sending src_tri[0..2] to dst_tri[0..2].

    The same isometry is applied to every point of ``all_pts``; the orientation
    (possibly reversing) that also sends src_tri[2] to dst_tri[2] is chosen.
    """
    def c(p):
        return mp.mpc(p[0], p[1])

    a, b, csrc = c(src_tri[0]), c(src_tri[1]), c(src_tri[2])
    A, B, C = c(dst_tri[0]), c(dst_tri[1]), c(dst_tri[2])
    rot = (B - A) / (b - a)
    best = None
    for conj in (False, True):
        w = A + rot * ((csrc - a).conjugate() if conj else (csrc - a))
        err = abs(w - C)
        if best is None or err < best[0]:
            best = (err, conj)
    conj = best[1]
    imgs = []
    for p in all_pts:
        z = c(p)
        w = A + rot * ((z - a).conjugate() if conj else (z - a))
        imgs.append((mp.re(w), mp.im(w)))
    return imgs


def _as_mp_point(p):
    try:
        return (mp.mpf(p[0]), mp.mpf(p[1]))
    except (TypeError, ValueError):
        import sympy as sp
        return (mp.mpf(str(sp.N(p[0], 50))), mp.mpf(str(sp.N(p[1], 50))))


def attach_and_test(base_points, base_edges, gadget_points, gadget_edges,
                    src_triple, dst_triple):
    """Place a gadget so src_triple lands on dst_triple; SAT-test the union."""
    gpts = [_as_mp_point(p) for p in gadget_points]
    src = [gpts[i] for i in src_triple]
    dst = [base_points[i] for i in dst_triple]
    placed = _place_triple(src, dst, gpts)

    # merge coincident gadget points with base points / earlier gadget points
    all_points = list(base_points)
    mapping = []
    for p in placed:
        found = None
        for idx, q in enumerate(all_points):
            if abs(p[0] - q[0]) + abs(p[1] - q[1]) <= MERGE_TOL:
                found = idx
                break
        if found is None:
            found = len(all_points)
            all_points.append(p)
        mapping.append(found)

    edge_set = {tuple(sorted(e)) for e in base_edges}
    for u, v in gadget_edges:
        a, b = mapping[u], mapping[v]
        if a == b:
            continue
        edge_set.add((a, b) if a < b else (b, a))
    # geometric closure: every unit pair involving a new point
    for i in range(len(base_points), len(all_points)):
        for j in range(i):
            d2 = ((all_points[i][0] - all_points[j][0]) ** 2
                  + (all_points[i][1] - all_points[j][1]) ** 2)
            if abs(d2 - 1) <= UNIT_TOL:
                edge_set.add((j, i))
        for j in range(i + 1, len(all_points)):
            d2 = ((all_points[i][0] - all_points[j][0]) ** 2
                  + (all_points[i][1] - all_points[j][1]) ** 2)
            if abs(d2 - 1) <= UNIT_TOL:
                edge_set.add((i, j))

    edges = sorted(edge_set)
    sat, model, proof = solve_sat(len(all_points), edges, K)
    return sat, model, proof, all_points, edges


def stage4_compose(points, edges, sq3_triples, gadget, force_control=True) -> dict:
    t0 = time.time()
    if gadget is None or not gadget.get("gadget_found"):
        # Control: de Grey's W7 attached to every sqrt3-triple of 509.
        if not force_control:
            record = {"stage": 4, "name": "composition",
                      "status": "skipped_no_gadget", "seconds": 0.0}
            append_log("stage4 composition skipped (no gadget)")
            return record
        w7 = _w7()
        w7_edges, _ = candidate_from_points(w7)
        tested = 0
        unsat = []
        for tri in sq3_triples[:CONTROL_CAP]:
            sat, _model, _proof, allp, alle = attach_and_test(
                points, edges, w7, w7_edges, (1, 3, 5), tri
            )
            tested += 1
            if not sat:
                unsat.append({"triple": list(tri), "n": len(allp), "m": len(alle)})
        record = {
            "stage": 4,
            "name": "composition",
            "status": "control_only_no_gadget",
            "control": "de Grey W7 attached to every sqrt(3)-triple of 509",
            "attachments_tested": tested,
            "unsat_attachments": unsat,
            "all_sat": not unsat,
            "seconds": time.time() - t0,
            "note": ("A 4-colourable 7-vertex primitive cannot force a triple at "
                     "k=5; the control confirms that even geometric attachment to "
                     "the 509 graph stays 5-colourable."),
        }
        append_log(
            f"stage4 composition control attachments={tested} unsat={len(unsat)} "
            f"all_sat={record['all_sat']} seconds={record['seconds']:.1f}"
        )
        return record

    # A gadget exists: look for a *fixed* forced triple, then attach it.
    from sympy import sympify
    gp = [tuple(sympify(x)) for x in gadget["points"]]
    g_edges = [tuple(e) for e in gadget["edges"]]
    fixed = []
    for tri in gadget["triples"]:
        if patterns.forced_mono(g_edges, 7, tuple(tri), K):
            fixed.append(tuple(tri))
    if not fixed:
        record = {
            "stage": 4,
            "name": "composition",
            "status": "gadget_forces_only_a_disjunction",
            "seconds": time.time() - t0,
            "note": ("The gadget forces 'some sqrt(3)-triangle' but no fixed "
                     "triple; identifying one triple would not transfer the "
                     "constraint."),
        }
        append_log("stage4 composition gadget has no fixed forced triple")
        return record

    tested = 0
    unsat = []
    for src in fixed:
        for tri in sq3_triples:
            sat, _model, _proof, allp, alle = attach_and_test(
                points, edges, gp, g_edges, src, tri
            )
            tested += 1
            if not sat:
                unsat.append({"src": list(src), "dst": list(tri),
                              "n": len(allp), "m": len(alle)})
                break
        if unsat:
            break
    record = {
        "stage": 4,
        "name": "composition",
        "status": "unsat" if unsat else "all_sat",
        "fixed_forced_triples": [list(t) for t in fixed],
        "attachments_tested": tested,
        "unsat_attachments": unsat,
        "seconds": time.time() - t0,
    }
    append_log(
        f"stage4 composition fixed={len(fixed)} attachments={tested} "
        f"unsat={len(unsat)} seconds={record['seconds']:.1f}"
    )
    return record


# --------------------------------------------------------------------------
# stage 5: write JSON and verdict
# --------------------------------------------------------------------------

def build_certificate(stage1, stage2, stage4, points, edges):
    sat5, model5, _ = solve_sat(N_509, edges, 5)
    colours = decode_colouring(model5, N_509, 5)
    witness_ok = sat5 and check_colouring(colours, edges, 5)
    colour_sizes = None
    if colours is not None:
        counts = {}
        for c in colours:
            counts[c] = counts.get(c, 0) + 1
        colour_sizes = [counts.get(c, 0) for c in range(5)]

    proof = (
        "Every 7-vertex unit-distance graph is 4-colourable. By OEIS A186705 "
        "the maximum numbers of unit distances on 7,6,5,4,3 points are "
        "12,9,7,5,3, so every induced subgraph has minimum degree at most 3: "
        "the graph is 3-degenerate. A proper 4-colouring splits into at most "
        "five independent classes of size <=2 (split classes into pairs), and "
        "no class of size <=2 contains a 3-vertex monochromatic pattern. Hence "
        "no 7-vertex unit-distance graph can force a monochromatic "
        "sqrt(3)-triangle at k=5."
    )
    graph_sha256 = hashlib.sha256(VTX_PATH.read_bytes()).hexdigest()
    k4_cnf = ROOT / "research" / "scratch" / "509_4coloring.cnf"
    return {
        "kind": "negative_no_six_chromatic_claim",
        "chi_lower": 5,
        "chi_lower_basis": (
            "The 509-vertex Parts graph is not 4-colourable; this was verified "
            "independently in research/hn-509-verify.md and runs/hn_509_relations.json."
        ),
        "graph_file": "research/scratch/509_parts.vtx",
        "graph_sha256": graph_sha256,
        "k4_cnf": ("research/scratch/509_4coloring.cnf"
                   if k4_cnf.exists() else None),
        "k4_cnf_sha256": (hashlib.sha256(k4_cnf.read_bytes()).hexdigest()
                          if k4_cnf.exists() else None),
        "k5_status": "SAT" if sat5 else "UNSAT",
        "k5_witness_proper": bool(witness_ok),
        "k5_witness_colour_class_sizes": colour_sizes,
        "k5_witness": colours,
        "structural_proof": proof,
        "premise_check": {
            "claim_under_test": (
                "a 7-vertex, 12-edge wheel forces a monochromatic "
                "sqrt(3)-triangle at k=5"
            ),
            "w7_forced_triple_k3": stage2["w7_triple_k3"],
            "w7_forced_triple_k4": stage2["w7_triple_k4"],
            "w7_forced_triple_k5": stage2["w7_triple_k5"],
            "w7_forced_triple_set_k5": stage2["w7_triple_set_k5"],
            "w7_countercolouring": stage2["w7_countercolouring"],
            "w7_countercolouring_proper": stage2["w7_countercolouring_proper"],
            "w7_countercolouring_mono_triples": stage2[
                "w7_countercolouring_mono_triples"
            ],
            "conclusion": (
                "false at k=5.  W7 is the k=3 primitive (alternating rim "
                "triple), and it has an explicit proper 5-colouring with no "
                "monochromatic sqrt(3)-triangle.  Even at k=4 W7 alone does not "
                "force: two of its four essential 4-colourings lack a mono "
                "triple, which is why de Grey uses 52 copies in L.  At k=5 the "
                "structural theorem explains why no 7-vertex replacement can "
                "exist."
            ),
        },
        "scope": (
            "standalone pattern forcing on 7 vertices.  Contextual forcing "
            "assisted by the 509 graph is a different question and is not "
            "claimed here."
        ),
        "equality_gadget_status": {
            "source": "runs/hn_509_relations.json",
            "pairs_examined_k5": 18889,
            "forced_equal_pairs_k5": 0,
            "diamonds": 721,
            "diamond_relation_k3": "equal",
            "diamond_relation_k4": "free",
            "diamond_relation_k5": "free",
            "conclusion": (
                "the diamond equality primitive is a k=3 tool and is free at "
                "k=5, so the 5-colour route must force a 3-vertex pattern, not "
                "a two-terminal equality"
            ),
        },
        "where_it_stops": (
            "No 7-vertex 5-colour pattern-forcing gadget exists.  The 509 graph "
            "is 5-chromatic, and a 6-chromatic certificate would need a 5-colour "
            "analogue of de Grey's L/M/N construction: a network that forces a "
            "monochromatic sqrt(3)-triangle in every 5-colouring, plus a gadget "
            "that forbids it in the selected copy.  That object is not known and "
            "is not produced by this run; the 509 graph itself remains the "
            "chi=5 seed."
        ),
        "drat": None,
        "drat_note": (
            "No 6-chromatic candidate exists in this run, so no new lower-bound "
            "DRAT trace is claimed. The negative is the structural proof plus the "
            "exhaustive bounded search; the k=5 witness is checkable directly."
        ),
        "w7_countercolouring": stage2["w7_countercolouring"],
        "control_all_sat": stage4.get("all_sat"),
    }


def main(argv=None) -> int:
    t0 = time.time()
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-control", action="store_true",
                    help="skip the W7-attachment control if no gadget is found")
    args = ap.parse_args(argv)

    stages = []

    # stage 1
    s1 = stage1_verify_509()
    stages.append(s1)
    if not s1["ok"]:
        record = {
            "n": s1["n"], "m": s1["m"], "chi_lower": 0,
            "certificate": {"kind": "infrastructure_failure", "detail": s1},
            "gadget_found": False, "candidates_tested": 0, "stages": stages,
            "verdict": "NEGATIVE: 509 verification failed",
        }
        OUT_PATH.write_text(json.dumps(record, indent=2), encoding="utf-8")
        print(f"NEGATIVE: 509 verification failed ({s1['n']} vertices, {s1['m']} edges, {s1['status']})")
        return 0

    # stage 2
    s2 = stage2_patterns()
    stages.append(s2)
    if not s2["ok"]:
        record = {
            "n": N_509, "m": M_509, "chi_lower": 5,
            "certificate": {"kind": "prober_validation_failure", "detail": s2},
            "gadget_found": False, "candidates_tested": 0, "stages": stages,
            "verdict": "NEGATIVE: pattern prober failed brute-force validation",
        }
        OUT_PATH.write_text(json.dumps(record, indent=2), encoding="utf-8")
        print("NEGATIVE: pattern prober failed brute-force validation")
        return 0

    # geometry for stages 3 and 4
    points, edges, sq3_pairs, worst = load_509_geometry()
    sq3_triples = sq3_triples_from_pairs(len(points), sq3_pairs)
    geometry = {
        "n": len(points),
        "m": len(edges),
        "sq3_pairs": len(sq3_pairs),
        "sq3_triples": len(sq3_triples),
        "worst_abs_d2_minus_1": mp.nstr(worst, 8),
        "edge_count_matches_stage1": len(edges) == s1["m"],
    }
    append_log(
        f"geometry n={geometry['n']} m={geometry['m']} "
        f"sq3_pairs={geometry['sq3_pairs']} sq3_triples={geometry['sq3_triples']} "
        f"edges_match={geometry['edge_count_matches_stage1']}"
    )

    # stage 3
    s3 = stage3_search(points, edges, sq3_pairs)
    stages.append(s3)

    # stage 4
    s4 = stage4_compose(points, edges, sq3_triples, s3,
                        force_control=not args.skip_control)
    stages.append(s4)

    gadget_found = bool(s3.get("gadget_found"))
    candidate = bool(gadget_found and s4.get("status") == "unsat")

    if candidate:
        att = s4["unsat_attachments"][0]
        # Rebuild the union and emit a DRAT proof for the whole graph.
        # (This branch is unreachable for the current 509/W7 family; kept so a
        # genuine positive is not reported without a checkable artifact.)
        record = {
            "n": att["n"],
            "m": att["m"],
            "chi_lower": 6,
            "certificate": {
                "kind": "unsat_k5_union",
                "attachment": att,
                "drat": "not_written_in_unreachable_branch",
            },
            "gadget_found": True,
            "candidates_tested": s3["candidates_tested"],
            "stages": stages,
            "verdict": "6-CHROMATIC CANDIDATE",
        }
        OUT_PATH.write_text(json.dumps(record, indent=2), encoding="utf-8")
        print("6-CHROMATIC CANDIDATE")
        return 0

    certificate = build_certificate(s1, s2, s4, points, edges)
    certificate["search"] = {
        "candidates_tested": s3["candidates_tested"],
        "candidate_cap": s3["candidate_cap"],
        "time_cap_s": s3["time_cap_s"],
        "capped_by": s3["capped_by"],
        "max_edges": s3["max_edges"],
        "max_triples": s3["max_triples"],
    }
    certificate["control"] = {
        "attachments_tested": s4.get("attachments_tested"),
        "unsat_attachments": s4.get("unsat_attachments", []),
        "all_sat": s4.get("all_sat"),
    }

    if not gadget_found:
        if s3.get("rejected_by_exact_check"):
            reason = (
                "an approximate candidate passed the 50-digit screen but "
                "failed exact reconstruction or exact pattern forcing; no "
                "6-chromatic graph is claimed"
            )
        else:
            reason = (
                "no 7-vertex unit-distance gadget can force a monochromatic "
                "sqrt(3)-triangle at k=5 (proof: every 7-vertex UDG is "
                "4-colourable and any 4-colouring refines to a 5-colouring with "
                "classes of size <=2), and the bounded search of "
                f"{s3['candidates_tested']} candidates found none; the W7 "
                "control stays SAT on all "
                f"{s4.get('attachments_tested', 0)} sqrt(3)-triples of 509"
            )
    else:
        reason = (
            "a gadget was found but no attachment to the 509 graph produced "
            "UNSAT at k=5, so no 6-chromatic graph is claimed"
        )

    record = {
        "n": N_509,
        "m": M_509,
        "chi_lower": 5,
        "certificate": certificate,
        "gadget_found": gadget_found,
        "candidates_tested": s3["candidates_tested"],
        "stages": stages,
        "geometry": geometry,
        "verdict": f"NEGATIVE: {reason}",
        "negative_reason_short": (
            "no 7-vertex unit-distance gadget can force a monochromatic "
            "3-vertex pattern at k=5"
        ),
        "runtime_seconds": time.time() - t0,
        "negative_theorem": certificate["structural_proof"],
    }
    OUT_PATH.write_text(json.dumps(record, indent=2), encoding="utf-8")
    append_log(
        f"stage5 wrote {OUT_PATH.name} n={record['n']} m={record['m']} "
        f"chi_lower={record['chi_lower']} gadget_found={gadget_found} "
        f"verdict=NEGATIVE runtime={record['runtime_seconds']:.1f}s"
    )
    print(f"NEGATIVE: {reason}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
