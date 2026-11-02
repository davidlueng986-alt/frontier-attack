"""Newton-polygon edge cascade at the (75,125) Jacobian corner, F_2(3,5), j=1.

CORNER DATA (unchanged from ``config75.py``)
--------------------------------------------
    Delta' = {(0,0),(3,0),(4,1),(0,5)}
    N(P)   = 3*Delta' = {(0,0),(9,0),(12,3),(0,15)}     124 lattice points
    N(Q)   = 5*Delta' = {(0,0),(15,0),(20,5),(0,25)}    321 lattice points
    [P,Q]  = x^2, reduced degrees (15,25)

The full-lattice coefficient system has 445 unknowns and 723 equations, which
no Groebner basis will touch.  This script implements the *edge cascade*: the
Newton polygons are decomposed into edges, each edge carries a weight vector,
and the leading part of ``det J`` with respect to that weight must vanish.

WHAT IS COMPUTED
----------------
Stage 1  Edge decomposition.  For each of N(P), N(Q): the four hull edges, the
         lattice points on each edge, the monomials, the primitive outward
         normal ``w``, the edge value ``d = w . alpha``, the primitive edge
         direction ``t = x^u y^v`` and the univariate edge form
         ``P_e = x^p0 y^q0 f(t)``.  The two polygons have the same normal fan,
         so the edges pair up canonically into four "edge slots".

Stage 2  The top edge.  "Top" is the slot of maximal filtration degree
         ``d_P + d_Q - w.(1,1)``, i.e. the face carrying the top-degree part of
         det J; here it is ``w = (1,1)``, ``d_P = 15``, ``d_Q = 25``,
         filtration degree 38 (against ``w.(2,0) = 2`` for the target x^2).

         Leading part: for *any* weight w,
             LP_w(det J) = det J(P_w, Q_w)
         where P_w, Q_w are the maximal-w-weight parts of P, Q.  On the top
         edge P_w = y^15 p(u), Q_w = y^25 q(u) with u = x/y, and
             det J(P_w, Q_w) = y^38 * (25 p'(u) q(u) - 15 p(u) q'(u)).
         The top-edge condition is therefore the polynomial identity
             25 p' q - 15 p q' = 0 .
         It is 31 homogeneous quadratic equations in the 13 + 21 = 34
         top-edge coefficients.

         A *direct* Groebner basis of that 34-variable system is attempted in
         a killable child process and the outcome (including a timeout) is
         reported verbatim.  The solution dimension is then obtained exactly
         from the structure of the equation: the Wronskian identity
             W(p^5, q^3) = p^4 q^2 (5 p' q - 3 p q')
         shows that 25 p' q = 15 p q' forces p^5 = C q^3; unique
         factorisation then gives p = c1 h^3, q = c2 h^5 with deg h = 4.
         So the strict-chart solution set is exactly
             {(c1 h^3, c2 h^5) : h in Sym^4, c1, c2 != 0},
         of dimension 4 + 1 + 1 = 6; fixing the two vertex coefficients
         a_(12,3) = c1 and b_(20,5) = c2 to 1 (both are vertices of the Newton
         polygons, hence nonzero, and both equations are homogeneous) gives
         dimension 4.  The family is verified by exact symbolic substitution.

Stage 3  Descent.  The edge order is the boundary cycle starting at the top
         edge and running clockwise, which is also decreasing filtration
         degree: (1,1), (1,-1), (0,-1), (-1,0).  Each further edge is solved in
         its own edge coefficients, with the shared vertex coefficients
         substituted from the previous solution (the top and the second edge
         share exactly the two vertices (12,3) and (20,5) of P and Q).  The
         descent stops at the first empty edge or after three edges.

         The two diagonal edges again obey the same Wronskian theorem: with
         g = gcd(d_P, d_Q), det J(P_e,Q_e) = 0 iff P_e = c1 G^(d_P/g) and
         Q_e = c2 G^(d_Q/g) for one edge form G of weighted degree g.  The
         horizontal and vertical edges give identically zero leading forms
         (both edge forms are univariate in the *same* variable), so they are
         vacuous: they impose no condition on their own coefficients.  That is
         reported as such, not as a pass.

WHAT IS NOT CLAIMED
-------------------
* An edge condition is *necessary* for det J = x^2, never sufficient.  A
  non-empty cascade says only that no contradiction was found in the first
  three leading-form conditions; it does not exhibit a Jacobian pair.
* The dimension reported for the top edge is the dimension of the strict
  Newton-polygon solution set (all vertex coefficients nonzero, two scalings
  fixed).  The literal variety of the 34-variable system also contains the
  degenerate component {p = 0} (all 13 top-edge P-coefficients zero, q free,
  dimension 21); both numbers are reported.
* No Groebner basis is reported where none was computed.  A timed-out call is
  recorded as ``timeout`` with the wall-clock cap.

Run:
    python problems/jacobian/cascade75.py
"""

from __future__ import annotations

import argparse
import json
import math
import random
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from config2_system import lattice_points  # noqa: E402
from config75 import (  # noqa: E402
    DELTA_PRIME,
    P_CORNERS,
    Q_CORNERS,
    det_coefficients,
)

DEFAULT_OUT = ROOT / "runs" / "jacobian_75_125_cascade.json"
LOG_PATH = ROOT / "runs" / "cascade75.log"

TARGET = (2, 0)
MAX_CASCADE_EDGES = 3

TOP_GB_Q_CAP_S = 90.0       # direct Groebner over Q, killed at this cap
TOP_GB_MODP_CAP_S = 60.0    # direct Groebner over GF(p), killed at this cap
EDGE_GB_CAP_S = 120.0       # per-descent-edge Groebner cap
MODULUS = 32003

X, Y, U, V = sp.symbols("x y u v")


# --------------------------------------------------------------------------
# logging
# --------------------------------------------------------------------------

def log_stage(stage, payload):
    """Append exactly one line per stage to runs/cascade75.log."""
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    line = "%s stage=%s %s\n" % (
        stamp, stage, json.dumps(payload, sort_keys=True, default=str))
    with LOG_PATH.open("a", encoding="utf-8") as fh:
        fh.write(line)


# --------------------------------------------------------------------------
# stage 1: edge decomposition
# --------------------------------------------------------------------------

def _primitive(vector):
    g = math.gcd(abs(int(vector[0])), abs(int(vector[1])))
    if g == 0:
        raise ValueError("zero direction")
    return (int(vector[0]) // g, int(vector[1]) // g)


def _signed_area2(points):
    n = len(points)
    return sum(
        points[i][0] * points[(i + 1) % n][1]
        - points[(i + 1) % n][0] * points[i][1]
        for i in range(n)
    )


def polygon_edges(corners):
    """Hull edges of a counter-clockwise lattice polygon.

    Returns one record per edge: start, end, primitive direction, primitive
    *outward* normal ``w`` and the edge value ``d = w . alpha`` (constant on
    the edge, and the maximum of ``w . alpha`` over the closed polygon).
    """
    pts = [tuple(int(c) for c in p) for p in corners]
    if _signed_area2(pts) <= 0:
        raise ValueError("corners must be counter-clockwise")
    edges = []
    n = len(pts)
    for i in range(n):
        a, b = pts[i], pts[(i + 1) % n]
        d = _primitive((b[0] - a[0], b[1] - a[1]))
        w = _primitive((d[1], -d[0]))          # outward normal for CCW
        value = w[0] * a[0] + w[1] * a[1]
        edges.append({
            "index": i,
            "start": a,
            "end": b,
            "direction": d,
            "normal": w,
            "value": value,
        })
    return edges


def edge_lattice_points(points, edge):
    """Lattice points of the polygon lying on the edge, ordered along it."""
    w, value = edge["normal"], edge["value"]
    on = [p for p in points if w[0] * p[0] + w[1] * p[1] == value]
    start, d = edge["start"], edge["direction"]
    on.sort(key=lambda p: (p[0] - start[0]) * d[0] + (p[1] - start[1]) * d[1])

    # the points must form a contiguous arithmetic progression along d
    for k, p in enumerate(on):
        expected = (start[0] + k * d[0], start[1] + k * d[1])
        if p != expected:
            raise AssertionError(
                "edge points are not contiguous: %s != %s" % (p, expected))
    return on


def edge_record(edge, points):
    """One edge, with its lattice points, monomials and univariate form."""
    w = edge["normal"]
    value = edge["value"]
    base = points[0]
    d = edge["direction"]
    monomials = [[int(p[0]), int(p[1])] for p in points]
    return {
        "index": int(edge["index"]),
        "start": list(edge["start"]),
        "end": list(edge["end"]),
        "direction": list(d),
        "normal": list(w),
        "value": int(value),
        "weight_grading": "w.(i,j) with w=(%d,%d); maximal value %d on this edge"
                          % (w[0], w[1], value),
        "n_points": len(points),
        "n_monomials": len(points),
        "points": monomials,
        "monomials": ["x^%d*y^%d" % (p[0], p[1]) for p in points],
        "edge_polynomial": " + ".join(
            "a_%d_%d*x^%d*y^%d" % (p[0], p[1], p[0], p[1]) for p in points),
        "univariate": {
            "base_point": list(base),
            "base_monomial": "x^%d*y^%d" % (base[0], base[1]),
            "t_monomial": "x^%d*y^%d" % (d[0], d[1]),
            "form": "P_e = x^%d*y^%d * f(t),  t = x^%d*y^%d,  "
                    "f = sum_k c_k t^k"
                    % (base[0], base[1], d[0], d[1]),
            "coefficient_lattice_points": monomials,
        },
        "is_horizontal": d[1] == 0,
        "is_vertical": d[0] == 0,
    }


def build_decomposition():
    """Edge decomposition of N(P) and N(Q) plus the paired edge slots."""
    lattice_p = sorted(lattice_points(P_CORNERS))
    lattice_q = sorted(lattice_points(Q_CORNERS))
    edges_p = polygon_edges(P_CORNERS)
    edges_q = polygon_edges(Q_CORNERS)

    records_p, records_q = [], []
    for edge in edges_p:
        pts = edge_lattice_points(lattice_p, edge)
        records_p.append(edge_record(edge, pts))
    for edge in edges_q:
        pts = edge_lattice_points(lattice_q, edge)
        records_q.append(edge_record(edge, pts))

    by_normal_q = {tuple(r["normal"]): r for r in records_q}
    normals_p = {tuple(r["normal"]) for r in records_p}
    if normals_p != set(by_normal_q):
        raise AssertionError(
            "normal fans differ: %s vs %s" % (normals_p, set(by_normal_q)))

    slots = []
    for rp in records_p:
        rq = by_normal_q[tuple(rp["normal"])]
        w = tuple(rp["normal"])
        d_p, d_q = rp["value"], rq["value"]
        filtration = d_p + d_q - w[0] - w[1]   # weight of the top face of det J
        target_weight = w[0] * TARGET[0] + w[1] * TARGET[1]
        slots.append({
            "normal": list(w),
            "d_P": int(d_p),
            "d_Q": int(d_q),
            "points_P": rp["points"],
            "points_Q": rq["points"],
            "n_points_P": rp["n_points"],
            "n_points_Q": rq["n_points"],
            "edge_P": rp,
            "edge_Q": rq,
            "filtration_degree": int(filtration),
            "target_weight": int(target_weight),
            "condition_required": bool(filtration > target_weight),
            "vertices_P": [rp["start"], rp["end"]],
            "vertices_Q": [rq["start"], rq["end"]],
            "leading_form_is_identically_zero": bool(
                rp["is_horizontal"] or rp["is_vertical"]),
        })

    # cascade order: start at the top (max filtration degree) and walk the
    # hull cycle clockwise, i.e. backwards through the CCW edge list
    n = len(slots)
    top = max(range(n), key=lambda k: (slots[k]["filtration_degree"],
                                       slots[k]["d_P"] + slots[k]["d_Q"]))
    order = [(top - k) % n for k in range(n)]
    for k in range(1, n):
        a = slots[order[k - 1]]["filtration_degree"]
        b = slots[order[k]]["filtration_degree"]
        if b > a:
            raise AssertionError("cascade order is not top-down: %s > %s" % (b, a))

    return {
        "lattice_points_P": lattice_p,
        "lattice_points_Q": lattice_q,
        "edges_P": records_p,
        "edges_Q": records_q,
        "slots": slots,
        "order": order,
        "top_slot": top,
    }


def _boundary_check(records, lattice_points_all):
    """Every boundary lattice point must lie on exactly one edge; count it."""
    counts = {}
    for r in records:
        for p in map(tuple, r["points"]):
            counts[p] = counts.get(p, 0) + 1
    boundary = {p for p in counts}
    interior = [p for p in lattice_points_all if tuple(p) not in boundary]
    return {
        "edge_point_incidences": int(sum(counts.values())),
        "distinct_boundary_points": len(boundary),
        "vertices_double_counted": int(sum(1 for v in counts.values() if v == 2)),
        "interior_points": len(interior),
        "min_incidence": min(counts.values()) if counts else None,
        "max_incidence": max(counts.values()) if counts else None,
    }


# --------------------------------------------------------------------------
# exact edge systems
# --------------------------------------------------------------------------

def sym_a(i, j):
    return sp.Symbol("a_%d_%d" % (int(i), int(j)))


def sym_b(i, j):
    return sp.Symbol("b_%d_%d" % (int(i), int(j)))


def build_edge_system(points_p, points_q):
    """Exact coefficient equations of det J restricted to two edge forms.

    ``det_coefficients`` (config75) accumulates det J as an integer bilinear
    form; here each monomial of det J becomes one equation.
    """
    det, pair_products = det_coefficients(points_p, points_q)
    a = [sym_a(*p) for p in points_p]
    b = [sym_b(*p) for p in points_q]
    equations, monomials, term_counts = [], [], []
    for mono in sorted(det):
        poly = sp.Integer(0)
        for (i, j), coeff in sorted(det[mono].items()):
            poly += int(coeff) * a[i] * b[j]
        equations.append(sp.expand(poly))
        monomials.append([int(mono[0]), int(mono[1])])
        term_counts.append(len(det[mono]))
    return {
        "a": a,
        "b": b,
        "equations": equations,
        "det_map": det,
        "monomials": monomials,
        "term_counts": term_counts,
        "pair_products": int(pair_products),
    }


def apply_fix(system, fix):
    """Substitute ``{symbol_name: value}`` and drop the fixed variables."""
    if not fix:
        return system["equations"], system["a"] + system["b"]
    subs = {sp.Symbol(k): sp.Integer(v) for k, v in fix.items()}
    eqs = [sp.expand(e.subs(subs)) for e in system["equations"]]
    gens = [g for g in system["a"] + system["b"] if g not in subs]
    return eqs, gens


# --------------------------------------------------------------------------
# Groebner driver: killable child process
# --------------------------------------------------------------------------

def gb_child(spec_path):
    """Child entry point: rebuild the requested system and compute a basis."""
    spec = json.loads(Path(spec_path).read_text(encoding="utf-8"))
    decomp = build_decomposition()
    slot = None
    for cand in decomp["slots"]:
        if list(cand["normal"]) == list(spec["normal"]):
            slot = cand
            break
    if slot is None:
        raise SystemExit("unknown edge normal %s" % (spec["normal"],))

    system = build_edge_system(slot["points_P"], slot["points_Q"])
    eqs, gens = apply_fix(system, spec.get("fix") or {})

    kwargs = {"order": spec.get("order", "grevlex")}
    if spec.get("modulus"):
        kwargs["modulus"] = int(spec["modulus"])
    if spec.get("method"):
        kwargs["method"] = spec["method"]

    t0 = time.perf_counter()
    basis = sp.groebner(eqs, *gens, **kwargs)
    seconds = time.perf_counter() - t0

    # Leading monomials must be taken with the order that was actually used.
    # SymPy rebuilds Poly objects in the default (lex) order, so taking
    # ``basis.polys[i].monoms()[0]`` would silently give the lex leader.
    from sympy.polys.orderings import grevlex as _grevlex_key

    order_name = kwargs["order"]
    supports = []
    for expr in basis.exprs:
        monoms = sp.Poly(expr, *gens).monoms()
        if order_name in ("grevlex", "degrevlex"):
            lead = max(monoms, key=_grevlex_key)
        elif order_name == "grlex":
            lead = max(monoms, key=lambda m: (sum(m), m))
        else:
            lead = monoms[0]
        supports.append([k for k, e in enumerate(lead) if e])

    sample = []
    for expr in basis.exprs[:8]:
        text = str(expr)
        sample.append(text if len(text) <= 400 else text[:400] + "...")

    sys.stdout.write(json.dumps({
        "status": "ok",
        "basis_size": len(basis.exprs),
        "lm_supports": supports,
        "gens": [str(g) for g in gens],
        "seconds": seconds,
        "sample": sample,
        "n_vars": len(gens),
        "n_equations": len(eqs),
        "modulus": spec.get("modulus"),
        "method": spec.get("method") or "buchberger",
        "order": kwargs["order"],
    }))
    sys.stdout.flush()
    return 0


def groebner_capped(normal, fix, cap_seconds, modulus=None, method=None):
    """Run the child with a hard wall-clock cap; never raise on timeout."""
    spec = {
        "normal": list(normal),
        "fix": fix or {},
        "modulus": modulus,
        "method": method,
        "order": "grevlex",
    }
    import os
    import tempfile
    spec_path = Path(tempfile.gettempdir()) / (
        "_cascade75_spec_%d.json" % os.getpid())
    spec_path.write_text(json.dumps(spec), encoding="utf-8")
    t0 = time.perf_counter()
    try:
        proc = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()),
             "--gb-child", str(spec_path)],
            capture_output=True, text=True, encoding="utf-8",
            timeout=cap_seconds, cwd=str(HERE),
        )
    except subprocess.TimeoutExpired:
        return {
            "status": "timeout",
            "cap_seconds": cap_seconds,
            "seconds": round(time.perf_counter() - t0, 3),
            "modulus": modulus,
            "method": method or "buchberger",
            "basis_size": None,
            "dimension": None,
            "empty": None,
        }
    finally:
        try:
            spec_path.unlink()
        except OSError:
            pass

    if proc.returncode != 0:
        return {
            "status": "error",
            "cap_seconds": cap_seconds,
            "seconds": round(time.perf_counter() - t0, 3),
            "error": (proc.stderr or "")[-2000:],
            "modulus": modulus,
            "method": method or "buchberger",
            "basis_size": None,
            "dimension": None,
            "empty": None,
        }
    text = (proc.stdout or "").strip()
    try:
        payload = json.loads(text.splitlines()[-1])
    except (ValueError, IndexError):
        return {
            "status": "error",
            "cap_seconds": cap_seconds,
            "seconds": round(time.perf_counter() - t0, 3),
            "error": "unparseable child output: " + text[-1000:],
            "modulus": modulus,
            "method": method or "buchberger",
            "basis_size": None,
            "dimension": None,
            "empty": None,
        }
    payload["cap_seconds"] = cap_seconds
    payload["seconds"] = round(time.perf_counter() - t0, 3)
    return payload


# --------------------------------------------------------------------------
# Krull dimension from a leading-monomial ideal
# --------------------------------------------------------------------------

class _Budget(Exception):
    pass


def _greedy_hitting(sets):
    used, hit = 0, [False] * len(sets)
    remaining = list(range(len(sets)))
    while remaining:
        counts = {}
        for idx in remaining:
            for v in sets[idx]:
                counts[v] = counts.get(v, 0) + 1
        v = max(counts, key=counts.get)
        used += 1
        remaining = [idx for idx in remaining if v not in sets[idx]]
    return used


def _greedy_disjoint(sets):
    taken, used = [], 0
    for s in sorted(sets, key=len):
        if all(s.isdisjoint(t) for t in taken):
            taken.append(s)
            used += 1
    return used


def krull_dimension_lm(supports, nvars, node_cap=400000, time_cap=60.0):
    """Krull dimension of a monomial ideal from its minimal generators.

    ``supports`` are the variable sets of the leading monomials.  For a
    monomial ideal M, dim(k[x]/M) is the largest size of a set S of variables
    such that no minimal generator of M is supported inside S, i.e.
    ``nvars - (minimum hitting set of the supports)``.  A constant leading
    monomial (empty support) means the unit ideal, dimension -1.
    """
    sets = []
    for support in supports:
        fs = frozenset(int(i) for i in support)
        if not fs:
            return {"dimension": -1, "unit_ideal": True, "empty": True,
                    "method": "constant leading monomial"}
        sets.append(fs)
    uniq = sorted(set(sets), key=len)
    minimal = []
    for s in uniq:
        if not any(t <= s for t in minimal):
            minimal.append(s)

    t0 = time.perf_counter()
    if nvars <= 18:
        masks = [sum(1 << i for i in s) for s in minimal]
        best = 0
        for mask in range(1 << nvars):
            k = bin(mask).count("1")
            if k <= best:
                continue
            if all((mask & sm) != sm for sm in masks):
                best = k
        return {"dimension": best, "unit_ideal": False, "empty": False,
                "method": "exhaustive subset scan over %d variables" % nvars,
                "min_supports": len(minimal), "seconds": round(
                    time.perf_counter() - t0, 3)}

    best0 = _greedy_hitting(minimal)
    best = [best0]
    nodes = [0]
    deadline = time.perf_counter() + time_cap
    sys.setrecursionlimit(20000)

    def rec(remaining, used):
        nodes[0] += 1
        if nodes[0] > node_cap or time.perf_counter() > deadline:
            raise _Budget
        if not remaining:
            if used < best[0]:
                best[0] = used
            return
        if used + 1 >= best[0]:
            return
        if used + _greedy_disjoint(remaining) >= best[0]:
            return
        pivot = min(remaining, key=len)
        for v in pivot:
            rest = [s for s in remaining if v not in s]
            rec(rest, used + 1)

    status = "ok"
    try:
        rec(minimal, 0)
    except _Budget:
        status = "budget"
    except RecursionError:
        status = "recursion"
    return {
        "dimension": nvars - best[0] if status == "ok" else None,
        "unit_ideal": False,
        "empty": False,
        "method": "branch-and-bound minimum hitting set over %d variables"
                  % nvars,
        "min_supports": len(minimal),
        "hitting_set": best[0],
        "nodes": nodes[0],
        "search_status": status,
        "seconds": round(time.perf_counter() - t0, 3),
    }


# --------------------------------------------------------------------------
# exact certificates
# --------------------------------------------------------------------------

def certificate_leading_form(slot, trials=2, seed=20260101):
    """Verify det J(P_e,Q_e) = M * [(p0 s0 - q0 r0) f g
                                    + (p0 v - q0 u) t f g'
                                    + (u s0 - v r0) t f' g].

    ``M = x^(p0+r0-1) y^(q0+s0-1)``, ``t = x^u y^v``, edge forms
    ``P_e = x^p0 y^q0 f(t)``, ``Q_e = x^r0 y^s0 g(t)``.  Checked on exact
    random integer specialisations of the edge coefficients.
    """
    pts_p, pts_q = slot["points_P"], slot["points_Q"]
    system = build_edge_system(pts_p, pts_q)
    a, b = system["a"], system["b"]
    u_, v_ = slot["edge_P"]["direction"]
    p0, q0 = pts_p[0]
    r0, s0 = pts_q[0]
    ts = sp.Symbol("t")
    rng = random.Random(seed)
    checks = []
    for _ in range(trials):
        while True:
            vals_a = [rng.randint(-9, 9) for _ in a]
            vals_b = [rng.randint(-9, 9) for _ in b]
            if all(vals_a) and all(vals_b):
                break
        sub = dict(zip(a, vals_a))
        sub.update(dict(zip(b, vals_b)))
        det_expr = sp.expand(sum(
            eq.subs(sub) * X ** mono[0] * Y ** mono[1]
            for eq, mono in zip(system["equations"], system["monomials"])))
        fa = sum(vals_a[k] * ts ** k for k in range(len(a)))
        gb = sum(vals_b[k] * ts ** k for k in range(len(b)))
        univariate = sp.expand(
            (p0 * s0 - q0 * r0) * fa * gb
            + (p0 * v_ - q0 * u_) * ts * fa * sp.diff(gb, ts)
            + (u_ * s0 - v_ * r0) * ts * sp.diff(fa, ts) * gb)
        t = X ** u_ * Y ** v_
        expected = sp.expand(
            X ** (p0 + r0 - 1) * Y ** (q0 + s0 - 1) * univariate.subs(ts, t))
        checks.append(bool(sp.expand(det_expr - expected) == 0))
    return {
        "identity": "det J(P_e,Q_e) = x^(p0+r0-1) y^(q0+s0-1) * "
                    "[(p0*s0-q0*r0)*f*g + (p0*v-q0*u)*t*f*g' "
                    "+ (u*s0-v*r0)*t*f'*g]",
        "trials": trials,
        "all_pass": all(checks),
        "per_trial": checks,
    }


def certificate_wronskian(seed=20260101):
    """Verify W(p^5, q^3) = p^4 q^2 (5 p' q - 3 p q') on exact instances."""
    rng = random.Random(seed)
    ok = True
    for _ in range(2):
        p = sum(sp.Integer(rng.randint(-9, 9)) * U ** k for k in range(13))
        q = sum(sp.Integer(rng.randint(-9, 9)) * U ** k for k in range(21))
        if p == 0 or q == 0:
            continue
        left = sp.expand(sp.diff(p ** 5, U) * q ** 3 - p ** 5 * sp.diff(q ** 3, U))
        right = sp.expand(p ** 4 * q ** 2 * (5 * sp.diff(p, U) * q
                                             - 3 * p * sp.diff(q, U)))
        ok = ok and (left == right)
    identity = "W(p^5,q^3) = p^4 q^2 (5 p' q - 3 p q')"
    return {
        "identity": identity,
        "identity_expanded": "W(p^5,q^3) = p^4 q^2 (25 p' q - 15 p q') / 5",
        "checked_on_random_degrees_12_20": bool(ok),
        "consequence": (
            "25 p' q - 15 p q' = 0 with p q != 0 forces W(p^5,q^3) = 0, "
            "hence p^5/q^3 is constant, hence p^5 = C q^3; unique "
            "factorisation gives p = c1 h^3, q = c2 h^5."),
    }


def certificate_family(slot, seed=20260101):
    """Verify that p = h^3, q = h^5 solves the top-edge equations exactly.

    Two independent checks:
      * symbolic: 25 p' q - 15 p q' vanishes identically for h generic;
      * coefficient-level: substitute p = h^3, q = h^5 into the actual 31
        coefficient equations at an exact random quartic h and check that
        every equation is zero.
    """
    c0, c1, c2, c3 = sp.symbols("c0 c1 c2 c3")
    h = U ** 4 + c3 * U ** 3 + c2 * U ** 2 + c1 * U + c0
    p = sp.expand(h ** 3)
    q = sp.expand(h ** 5)
    wronskian_symbolic = sp.expand(25 * sp.diff(p, U) * q - 15 * p * sp.diff(q, U))

    pts_p, pts_q = slot["points_P"], slot["points_Q"]
    system = build_edge_system(pts_p, pts_q)
    rng = random.Random(seed)
    trial_values = [rng.randint(-6, 6) for _ in range(4)]
    if trial_values[0] == 0:
        trial_values[0] = 1
    point = {c0: trial_values[0], c1: trial_values[1],
             c2: trial_values[2], c3: trial_values[3]}
    coeffs_p = sp.Poly(p.subs(point), U).all_coeffs()[::-1]
    coeffs_q = sp.Poly(q.subs(point), U).all_coeffs()[::-1]
    coeffs_p += [sp.Integer(0)] * (len(system["a"]) - len(coeffs_p))
    coeffs_q += [sp.Integer(0)] * (len(system["b"]) - len(coeffs_q))
    sub = dict(zip(system["a"], coeffs_p))
    sub.update(dict(zip(system["b"], coeffs_q)))
    residuals = [sp.expand(eq.subs(sub)) for eq in system["equations"]]
    return {
        "family": "p = h^3, q = h^5, h = u^4 + c3 u^3 + c2 u^2 + c1 u + c0",
        "wronskian_symbolically_zero": bool(wronskian_symbolic == 0),
        "equations_checked": len(residuals),
        "all_equations_vanish": bool(all(r == 0 for r in residuals)),
        "trial_point": dict(zip(["c0", "c1", "c2", "c3"], trial_values)),
        "nonzero_residuals": int(sum(1 for r in residuals if r != 0)),
    }


def certificate_degenerate_component(slot):
    """The {p = 0} component: all top-edge P-coefficients zero, q arbitrary."""
    pts_p, pts_q = slot["points_P"], slot["points_Q"]
    system = build_edge_system(pts_p, pts_q)
    sub = {sym: sp.Integer(0) for sym in system["a"]}
    residuals = [sp.expand(eq.subs(sub)) for eq in system["equations"]]
    return {
        "component": "{all top-edge P coefficients = 0, Q edge form free}",
        "equations_vanish": bool(all(r == 0 for r in residuals)),
        "dimension": len(system["b"]),
    }


# --------------------------------------------------------------------------
# stage 2: the top edge
# --------------------------------------------------------------------------

def analyse_top_edge(slot, next_slot, cap_q, cap_modp):
    pts_p, pts_q = slot["points_P"], slot["points_Q"]
    system = build_edge_system(pts_p, pts_q)
    a, b = system["a"], system["b"]

    n_unknowns = len(a) + len(b)
    payload = {
        "slot_normal": slot["normal"],
        "d_P": slot["d_P"],
        "d_Q": slot["d_Q"],
        "filtration_degree": slot["filtration_degree"],
        "target_weight": slot["target_weight"],
        "unknowns": n_unknowns,
        "equations": len(system["equations"]),
        "pair_products": system["pair_products"],
        "det_monomials": system["monomials"],
        "max_terms_per_equation": max(system["term_counts"]),
        "all_equations_homogeneous_degree_2": bool(all(
            sp.Poly(eq, *(a + b)).total_degree() == 2
            for eq in system["equations"])),
        "edge_points_P": pts_p,
        "edge_points_Q": pts_q,
        "base_point_P": pts_p[0],
        "base_point_Q": pts_q[0],
        "t_monomial_P": "x^%d*y^%d" % slot["edge_P"]["direction"],
        "leading_form": (
            "det J(P_e,Q_e) = y^38 * (25 p'(u) q(u) - 15 p(u) q'(u)), "
            "u = x/y, p = sum a_(i,15-i) u^i, q = sum b_(j,25-j) u^j; "
            "equivalently x^12*y^3*f(y/x) and x^20*y^5*g(y/x)"),
        "leading_form_identity": certificate_leading_form(slot),
        "wronskian": certificate_wronskian(),
        "condition": "25 p' q - 15 p q' = 0 identically  (31 equations)",
    }

    # strict chart: normalise the two vertex coefficients that the top edge
    # shares with the next edge of the cascade (they are nonzero because they
    # are vertices of N(P), N(Q)); this is exactly the substitution used
    # again when descending
    shared_p, shared_q = shared_vertices(slot, next_slot)
    fix = {}
    for p in shared_p:
        fix["a_%d_%d" % tuple(p)] = 1
    for p in shared_q:
        fix["b_%d_%d" % tuple(p)] = 1
    if not fix:
        raise AssertionError("the top edge has no shared vertex with the "
                             "next edge; cannot fix the scaling consistently")
    payload["strict_chart"] = {
        "normalization": fix,
        "normalized_vertices": {
            "P": [list(p) for p in shared_p],
            "Q": [list(p) for p in shared_q],
        },
        "why": (
            "the coefficients of lattice vertices are nonzero because those "
            "points are vertices of N(P), N(Q); the equations are homogeneous "
            "of degree 1 in each block, so setting the two shared vertex "
            "coefficients to 1 removes exactly the two scalings, and the "
            "normalisation is the same one used on the descending edge"),
        "unknowns": n_unknowns - len(fix),
        "equations": len(system["equations"]),
    }
    gb_q = groebner_capped(slot["normal"], fix, cap_q)
    if gb_q["status"] == "ok":
        dim = krull_dimension_lm(gb_q["lm_supports"], gb_q["n_vars"])
        gb_q["dimension"] = dim["dimension"]
        gb_q["empty"] = dim["empty"]
        gb_q["dimension_detail"] = dim
    payload["groebner_over_Q"] = gb_q

    # (2) direct Groebner over GF(p) with the F5B engine, same cap policy
    gb_p = groebner_capped(slot["normal"], fix, cap_modp,
                           modulus=MODULUS, method="f5b")
    if gb_p["status"] == "ok":
        dim = krull_dimension_lm(gb_p["lm_supports"], gb_p["n_vars"])
        gb_p["dimension"] = dim["dimension"]
        gb_p["empty"] = dim["empty"]
        gb_p["dimension_detail"] = dim
    payload["groebner_mod_p"] = gb_p

    # (3) exact solution dimension
    g = math.gcd(int(slot["d_P"]), int(slot["d_Q"]))
    n_p = len(pts_p) - 1
    n_q = len(pts_q) - 1
    if n_p * g % slot["d_P"] or n_q * g % slot["d_Q"]:
        raise AssertionError("edge forms are not compatible powers")
    dim_normalized = n_p * g // slot["d_P"]
    if dim_normalized != n_q * g // slot["d_Q"]:
        raise AssertionError("power exponents disagree between P and Q")

    family = certificate_family(slot)
    degenerate = certificate_degenerate_component(slot)

    groebner_computed = [
        r for r in (gb_q, gb_p) if r.get("status") == "ok"
    ]
    if groebner_computed:
        dims = [r["dimension"] for r in groebner_computed]
        if any(d is not None for d in dims):
            exercised = [d for d in dims if d is not None]
            if any(d != dim_normalized for d in exercised):
                raise AssertionError(
                    "computed dimension %s disagrees with the structure "
                    "theorem %d" % (exercised, dim_normalized))
            dimension_method = (
                "Groebner basis (strict chart) confirms the structure "
                "theorem: %s" % dims)
        else:
            dimension_method = "structure theorem (Groebner dimension search capped)"
    else:
        dimension_method = (
            "structure theorem (Wronskian): direct Groebner basis timed out "
            "at the stated caps; see groebner_over_Q / groebner_mod_p")

    payload["dimension"] = int(dim_normalized)
    payload["dimension_scaled"] = int(dim_normalized + 2)
    payload["dimension_literal_full_system"] = int(
        max(dim_normalized + 2, degenerate["dimension"]))
    payload["dimension_method"] = dimension_method
    payload["power_structure"] = {
        "g": int(g),
        "claim": "det J(P_e,Q_e) = 0  =>  P_e = c1 G^(d_P/g), "
                 "Q_e = c2 G^(d_Q/g), deg G = g",
        "parameter_count": (
            "deg G = %d gives %d coefficients; the shape P_e = x^p0 y^q0 f(t) "
            "forces G = x^(p0/(d_P/g)) y^(q0/(d_P/g)) k(t) with k monic of "
            "degree %d, so %d free parameters for G; plus c1, c2 gives %d, "
            "and fixing the two vertex coefficients leaves %d"
            % (g, g + 1, dim_normalized, dim_normalized,
               dim_normalized + 2, dim_normalized)),
        "wronskian_reduction": (
            "in the edge parameter t: d_Q f g' - d_P f' g = 0, i.e. "
            "(f^(d_Q/g))' (g^(d_P/g)) - ... = 0, whose solutions are "
            "f = c k^(d_P/g), g = c' k^(d_Q/g) by the same ODE argument"),
    }
    payload["family_certificate"] = family
    payload["degenerate_component"] = degenerate
    payload["empty"] = False
    return payload


# --------------------------------------------------------------------------
# stage 3: descend
# --------------------------------------------------------------------------

def shared_vertices(slot_a, slot_b):
    """Vertices common to two consecutive edges, per polygon."""
    if slot_a is None or slot_b is None:
        return [], []
    va = {tuple(p) for p in slot_a["vertices_P"]}
    vb = {tuple(p) for p in slot_b["vertices_P"]}
    qa = {tuple(p) for p in slot_a["vertices_Q"]}
    qb = {tuple(p) for p in slot_b["vertices_Q"]}
    return (sorted(va & vb, key=lambda z: (z[0], z[1])),
            sorted(qa & qb, key=lambda z: (z[0], z[1])))


def analyse_descent_edge(slot, index, previous, cap, fixed_symbols=None):
    """Solve one further edge condition in its own edge coefficients."""
    pts_p, pts_q = slot["points_P"], slot["points_Q"]
    system = build_edge_system(pts_p, pts_q)

    # lattice points shared with the previously solved edge (the common
    # vertex of two consecutive hull edges), for each polygon separately
    shared_p, shared_q = [], []
    if previous is not None:
        prev_p = set(map(tuple, previous["points_P"]))
        prev_q = set(map(tuple, previous["points_Q"]))
        shared_p = sorted(set(map(tuple, pts_p)) & prev_p)
        shared_q = sorted(set(map(tuple, pts_q)) & prev_q)

    record = {
        "edge_index": index,
        "normal": slot["normal"],
        "d_P": slot["d_P"],
        "d_Q": slot["d_Q"],
        "filtration_degree": slot["filtration_degree"],
        "target_weight": slot["target_weight"],
        "unknowns": len(system["a"]) + len(system["b"]),
        "equations": len(system["equations"]),
        "edge_points_P": pts_p,
        "edge_points_Q": pts_q,
        "n_points_P": len(pts_p),
        "n_points_Q": len(pts_q),
        "det_monomials": system["monomials"],
        "shared_points_with_previous": {
            "P": [list(p) for p in shared_p],
            "Q": [list(p) for p in shared_q],
        },
        "vertices_P": slot["vertices_P"],
        "vertices_Q": slot["vertices_Q"],
    }

    if not system["equations"]:
        record.update({
            "leading_form_identically_zero": True,
            "reason": (
                "both edge forms are univariate in the same variable "
                "(%s edge), so det J(P_e,Q_e) = 0 identically for every "
                "choice of edge coefficients"
                % ("horizontal" if slot["edge_P"]["is_horizontal"]
                   else "vertical")),
            "vacuous": True,
            "dimension": record["unknowns"],
            "dimension_own_space": record["unknowns"],
            "empty": False,
            "groebner": {"status": "not_run",
                         "reason": "edge system is empty (0 equations)"},
            "note": (
                "the edge-local condition is automatically satisfied; the "
                "first non-vacuous condition at this weight is of lower "
                "filtration order and involves interior coefficients, which "
                "is outside the 3-edge cascade"),
        })
        return record

    record["leading_form_identically_zero"] = False
    g = math.gcd(int(slot["d_P"]), int(slot["d_Q"]))
    n_p = len(pts_p) - 1
    dim_normalized = n_p * g // slot["d_P"]
    record["power_structure_g"] = int(g)
    record["dimension_structure_theorem"] = int(dim_normalized)

    # substitute the solved values at the shared vertices: the top-edge chart
    # fixes a_(12,3) = c1 = 1 and b_(20,5) = c2 = 1, and those two points are
    # exactly the vertices shared by the top and the second edge
    fixed_symbols = fixed_symbols or set()
    fix = {}
    for p in shared_p:
        name = "a_%d_%d" % p
        if name in fixed_symbols:
            fix[name] = 1
    for p in shared_q:
        name = "b_%d_%d" % p
        if name in fixed_symbols:
            fix[name] = 1
    record["substitution"] = {
        "fix": fix,
        "values": "1 (the top-edge family gives a_(12,3) = c1 and "
                  "b_(20,5) = c2; both were normalised to 1)",
        "consistent": True,
        "shared_but_not_fixed": {
            "P": [list(p) for p in shared_p
                  if "a_%d_%d" % tuple(p) not in fixed_symbols],
            "Q": [list(p) for p in shared_q
                  if "b_%d_%d" % tuple(p) not in fixed_symbols],
        },
    }

    gb = groebner_capped(slot["normal"], fix, cap)
    if gb.get("status") == "ok":
        dim = krull_dimension_lm(gb["lm_supports"], gb["n_vars"])
        gb["dimension"] = dim["dimension"]
        gb["empty"] = dim["empty"]
        gb["dimension_detail"] = dim
    record["groebner"] = gb
    if gb.get("status") == "ok" and gb.get("dimension") is not None:
        record["dimension"] = gb["dimension"]
        record["dimension_method"] = "Groebner basis (grevlex) + LM ideal"
        if gb["dimension"] != dim_normalized:
            record["dimension_disagreement"] = {
                "groebner": gb["dimension"],
                "structure_theorem": dim_normalized,
            }
    else:
        record["dimension"] = int(dim_normalized)
        record["dimension_method"] = (
            "structure theorem (Groebner %s)" % gb.get("status"))
    record["empty"] = bool(gb.get("empty") is True)
    record["family_certificate"] = _edge_family_check(slot)
    return record


def _edge_family_check(slot):
    """Exact check that the power family solves the edge equations."""
    pts_p, pts_q = slot["points_P"], slot["points_Q"]
    system = build_edge_system(pts_p, pts_q)
    g = math.gcd(int(slot["d_P"]), int(slot["d_Q"]))
    exp_p = int(slot["d_P"]) // g
    exp_q = int(slot["d_Q"]) // g
    t = sp.Symbol("t")
    # G has the shape x^a y^b k(t) with k monic of degree dim
    dim = (len(pts_p) - 1) * g // int(slot["d_P"])
    r = sp.Symbol("r")
    if dim == 1:
        k = t - r
    else:
        coeffs = sp.symbols("k0:%d" % dim)
        k = t ** dim + sum(coeffs[i] * t ** i for i in range(dim))
    # P_e / Q_e are determined by the monomial factor of the edge supports
    a0, b0 = pts_p[0], pts_q[0]
    gp = X ** (a0[0] // exp_p) * Y ** (a0[1] // exp_p)
    gq = X ** (b0[0] // exp_q) * Y ** (b0[1] // exp_q)
    up, vp = slot["edge_P"]["direction"]
    tp = X ** up * Y ** vp
    P_edge = sp.expand(gp ** exp_p * k.subs(t, tp) ** exp_p)
    Q_edge = sp.expand(gq ** exp_q * k.subs(t, tp) ** exp_q)
    det = sp.expand(sp.diff(P_edge, X) * sp.diff(Q_edge, Y)
                    - sp.diff(P_edge, Y) * sp.diff(Q_edge, X))
    ok = bool(det == 0)
    return {
        "family": "P_e = G^%d, Q_e = G^%d, G = x^%d*y^%d*k(t), t = x^%d*y^%d, "
                  "deg k = %d" % (exp_p, exp_q, a0[0] // exp_p,
                                  a0[1] // exp_p, up, vp, dim),
        "det_J_identically_zero": ok,
        "free_parameters_in_G": int(dim),
    }


def run_descent(decomp, top_payload, cap, fixed_symbols=None):
    order = decomp["order"]
    slots = decomp["slots"]
    records = []
    first_empty = None
    previous_slot = slots[order[0]]
    processed = 1

    for step in range(1, MAX_CASCADE_EDGES):
        if step >= len(order):
            break
        slot = slots[order[step]]
        rec = analyse_descent_edge(slot, step + 1, previous_slot, cap,
                                   fixed_symbols)
        records.append(rec)
        processed += 1
        previous_slot = slot
        if rec.get("empty"):
            first_empty = "edge %d (normal %s)" % (step + 1, slot["normal"])
            break
    return records, first_empty, processed


# --------------------------------------------------------------------------
# report
# --------------------------------------------------------------------------

def build_report(cap_top_q=TOP_GB_Q_CAP_S, cap_top_p=TOP_GB_MODP_CAP_S,
                 cap_edge=EDGE_GB_CAP_S):
    t0 = time.perf_counter()
    decomp = build_decomposition()
    slots = decomp["slots"]
    order = decomp["order"]
    top = slots[order[0]]

    edges_p = decomp["edges_P"]
    edges_q = decomp["edges_Q"]

    stage1 = {
        "num_edges_NP": len(edges_p),
        "num_edges_NQ": len(edges_q),
        "num_edges_total": len(edges_p) + len(edges_q),
        "points_per_edge_NP": [e["n_points"] for e in edges_p],
        "points_per_edge_NQ": [e["n_points"] for e in edges_q],
        "monomials_per_edge_NP": [e["n_monomials"] for e in edges_p],
        "monomials_per_edge_NQ": [e["n_monomials"] for e in edges_q],
        "boundary_check_NP": _boundary_check(edges_p, decomp["lattice_points_P"]),
        "boundary_check_NQ": _boundary_check(edges_q, decomp["lattice_points_Q"]),
        "lattice_points_NP": len(decomp["lattice_points_P"]),
        "lattice_points_NQ": len(decomp["lattice_points_Q"]),
        "edges_NP": edges_p,
        "edges_NQ": edges_q,
        "edge_slots": [{
            "normal": s["normal"],
            "d_P": s["d_P"],
            "d_Q": s["d_Q"],
            "n_points_P": s["n_points_P"],
            "n_points_Q": s["n_points_Q"],
            "monomials_P": s["n_points_P"],
            "monomials_Q": s["n_points_Q"],
            "filtration_degree": s["filtration_degree"],
            "target_weight": s["target_weight"],
            "condition_required": s["condition_required"],
            "leading_form_is_identically_zero":
                s["leading_form_is_identically_zero"],
        } for s in slots],
        "cascade_order": [
            {"edge": k + 1, "hull_index": order[k],
             "normal": slots[order[k]]["normal"],
             "filtration_degree": slots[order[k]]["filtration_degree"]}
            for k in range(len(order))
        ],
        "grading": {
            "definition": "deg_w(x^i y^j) = w.(i,j) for the outward normal w "
                          "of the edge; the edge is the maximal face",
            "leading_form": "LP_w(det J) = det J(P_w,Q_w) for any weight w",
            "filtration_degree_of_det_J_face":
                "d_P + d_Q - w.(1,1)",
            "target": "x^2 has w-weight 2*w_1",
        },
    }
    log_stage("1_edges", {
        "num_edges_NP": stage1["num_edges_NP"],
        "num_edges_NQ": stage1["num_edges_NQ"],
        "points_per_edge_NP": stage1["points_per_edge_NP"],
        "points_per_edge_NQ": stage1["points_per_edge_NQ"],
        "cascade_order": [e["normal"] for e in stage1["cascade_order"]],
    })

    top_payload = analyse_top_edge(top, slots[order[1]], cap_top_q, cap_top_p)
    log_stage("2_top_edge", {
        "normal": top_payload["slot_normal"],
        "unknowns": top_payload["unknowns"],
        "equations": top_payload["equations"],
        "groebner_over_Q": top_payload["groebner_over_Q"]["status"],
        "groebner_mod_p": top_payload["groebner_mod_p"]["status"],
        "top_dim": top_payload["dimension"],
        "empty": top_payload["empty"],
    })

    descent, first_empty, processed = run_descent(
        decomp, top_payload, cap_edge,
        set(top_payload["strict_chart"]["normalization"]))

    if top_payload["empty"]:
        verdict = "TOP_EDGE_EMPTY"
    elif first_empty is not None:
        verdict = "INCONSISTENT_AT_EDGE %d" % processed
    else:
        verdict = "DESCENDED %d EDGES" % processed

    log_stage("3_descent", {
        "edges_processed": processed,
        "edge_dims": [top_payload["dimension"]] + [
            r.get("dimension") for r in descent],
        "first_empty_edge": first_empty,
        "verdict": verdict,
    })

    if verdict == "TOP_EDGE_EMPTY":
        verdict_detail = "the top leading-form condition has no solution"
    elif verdict.startswith("INCONSISTENT_AT_EDGE"):
        verdict_detail = (
            "edge %d of the cascade is empty: its leading-form condition "
            "has no solution compatible with the coefficients fixed by the "
            "previous edges" % processed)
    else:
        verdict_detail = (
            "no contradiction was found in the first %d edge conditions; "
            "each edge condition is necessary for det J = x^2, so this is a "
            "filter, not a solution" % processed)

    report = {
        "problem": "jacobian-75-125-cascade",
        "case": "F_2(3,5), j=1",
        "method": "Newton-polygon edge filtration / leading-form cascade",
        "corner": {
            "Delta_prime": [list(p) for p in DELTA_PRIME],
            "N_P_corners": [list(p) for p in P_CORNERS],
            "N_Q_corners": [list(p) for p in Q_CORNERS],
            "normalisation": "[P,Q] = x^2, reduced degrees (15,25)",
            "full_lattice_unknowns": 445,
            "full_lattice_equations": 723,
        },
        "edges": stage1,
        "points_per_edge": {
            "N(P)": stage1["points_per_edge_NP"],
            "N(Q)": stage1["points_per_edge_NQ"],
        },
        "top_edge": top_payload,
        "top_dim": top_payload["dimension"],
        "descent": descent,
        "first_empty_edge": first_empty,
        "verdict": verdict,
        "verdict_detail": verdict_detail,
        "conventions": {
            "top_edge": "the edge of maximal filtration degree "
                        "d_P + d_Q - w.(1,1); here w = (1,1)",
            "leading_part": "LP_w(det J) = det J(P_w, Q_w)",
            "strict_chart": "the vertex coefficients (12,3) and (20,5) are "
                            "nonzero (they are vertices of N(P), N(Q)) and "
                            "are normalised to 1; every equation is "
                            "homogeneous of degree 1 in each block",
            "dimension": "affine dimension of the solution set in the strict "
                         "chart",
            "edge_order": "boundary cycle from the top edge, clockwise; "
                          "equals decreasing filtration degree",
        },
        "timings": {
            "stage_seconds": {},
            "total_seconds": None,
        },
        "elapsed_seconds": None,
    }
    report["timings"]["total_seconds"] = round(time.perf_counter() - t0, 3)
    report["elapsed_seconds"] = report["timings"]["total_seconds"]
    return report


def print_report(rep):
    e = rep["edges"]
    print("== (75,125) Newton-polygon edge cascade ==")
    print("N(P): %d edges, points/edge %s, %d lattice points"
          % (e["num_edges_NP"], e["points_per_edge_NP"],
             e["lattice_points_NP"]))
    print("N(Q): %d edges, points/edge %s, %d lattice points"
          % (e["num_edges_NQ"], e["points_per_edge_NQ"],
             e["lattice_points_NQ"]))
    for s in e["edge_slots"]:
        print("  slot w=%-8s d_P=%-3d d_Q=%-3d points %dx%d "
              "filtration degree %-2d target weight %d%s"
              % (s["normal"], s["d_P"], s["d_Q"], s["n_points_P"],
                 s["n_points_Q"], s["filtration_degree"], s["target_weight"],
                 "  (leading form identically zero)"
                 if s["leading_form_is_identically_zero"] else ""))
    top = rep["top_edge"]
    print()
    print("top edge  w=%s: %d unknowns, %d equations"
          % (top["slot_normal"], top["unknowns"], top["equations"]))
    print("  leading form: %s" % top["leading_form"])
    print("  Groebner over Q     : %s (cap %.0f s)"
          % (top["groebner_over_Q"]["status"],
             top["groebner_over_Q"]["cap_seconds"]))
    print("  Groebner over GF(%d): %s (cap %.0f s)"
          % (MODULUS, top["groebner_mod_p"]["status"],
             top["groebner_mod_p"]["cap_seconds"]))
    print("  dimension (strict)  : %d" % top["dimension"])
    print("  dimension (scaled)  : %d" % top["dimension_scaled"])
    print("  dimension (literal, incl. degenerate {p=0}): %d"
          % top["dimension_literal_full_system"])
    print("  family              : %s" % top["family_certificate"]["family"])
    print("  family verified     : %s"
          % top["family_certificate"]["all_equations_vanish"])
    print()
    for rec in rep["descent"]:
        print("edge %d  w=%s: %s" % (rec["edge_index"], rec["normal"],
                                     "VACUOUS" if rec.get("vacuous")
                                     else "equations %d" % rec["equations"]))
        if rec.get("vacuous"):
            print("   %s" % rec["reason"])
        else:
            print("   dimension %s (%s)"
                  % (rec.get("dimension"), rec.get("dimension_method")))
    print()
    print("first empty edge: %s" % rep["first_empty_edge"])
    print("VERDICT: %s" % rep["verdict"])
    print("elapsed: %.3f s" % rep["elapsed_seconds"])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--no-write", action="store_true")
    ap.add_argument("--cap-top-q", type=float, default=TOP_GB_Q_CAP_S)
    ap.add_argument("--cap-top-modp", type=float, default=TOP_GB_MODP_CAP_S)
    ap.add_argument("--cap-edge", type=float, default=EDGE_GB_CAP_S)
    ap.add_argument("--gb-child", type=str, default=None,
                    help=argparse.SUPPRESS)
    args = ap.parse_args(argv)

    if args.gb_child:
        return gb_child(args.gb_child)

    t0 = time.perf_counter()
    rep = build_report(cap_top_q=args.cap_top_q,
                       cap_top_p=args.cap_top_modp,
                       cap_edge=args.cap_edge)
    print_report(rep)

    if not args.no_write:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(rep, indent=2), encoding="utf-8")
        try:
            shown = args.out.relative_to(ROOT)
        except ValueError:
            shown = args.out
        print("\nwritten: %s" % shown)
    print("total wall clock: %.3f s" % (time.perf_counter() - t0))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
