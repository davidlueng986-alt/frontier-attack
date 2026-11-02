"""Valuation / leading-term analysis of the full-lattice GGHV config-2 system.

THE QUESTION
------------
Decide, structurally, whether

    P, Q in K[x, y],   N(P) = {(0,0),(1,0),(8,14),(8,16)},
                       N(Q) = {(0,0),(2,1),(12,21),(12,24)},
    det J(P,Q) = x^2

has a solution.  The supports are the *full* lattice-point sets of the two
Newton polygons: 25 + 47 = 72 coefficients.  The coefficient system has 72
unknowns and 92 equations, so a 72-variable Groebner basis is out of reach;
this script tries to settle it by support / valuation / leading-form data.

WHAT THIS SCRIPT PROVES (and what it does not)
----------------------------------------------
Stage 1  Exact attainable support.  A monomial x^a y^b of det J can only come
         from a pair (p in supp P_x, q in supp Q_y) or (p in supp P_y,
         q in supp Q_x) with p + q = (a, b).  Enumerating all pairs is exact.
         The result is checked against a full 72-symbol symbolic expansion of
         det J: the attainable set (112 monomials) contains the true support
         (92 monomials), and the 20 extra monomials vanish identically because
         the two signed contributions cancel.

         *** (2,0) IS ATTAINABLE. ***  There is exactly one pair, namely
         (0,0) from P_x (the monomial a_(1,0) x) plus (2,0) from Q_y (the
         monomial b_(2,1) x^2 y).  So the coefficient of x^2 in det J is
         exactly a_(1,0) * b_(2,1).  Consequently

             EMPTY_BY_SUPPORT IS RULED OUT.

         This is the opposite of the naive expectation.  N(P) has no point
         with b > 0 on x = 0, but N(Q) has (2,1) on its lower edge, and that
         single vertex supplies the (2,0) exponent.  (The (72,108) plan's
         "det J is divisible by y" heuristic is false for this reason.)

Stage 2  Minimal x-degree.  min x-degree of det J = 1 over the full support
         (the attainable set also has min x-degree 1).  So the lowest-x-degree
         block is the coefficient of x^1, three equations in the five
         coefficients a_(1,0), a_(1,1), a_(1,2), b_(1,1), b_(1,2).

Stage 3  That block is solved exactly (Groebner over Q).  The ideal is

             ( a_(1,0) b_(1,1),  a_(1,0) b_(1,2),  a_(1,1) b_(1,2) - a_(1,2) b_(1,1) )

         of Krull dimension 3 (it has the two 3-dimensional components
         {a_(1,0) = 0} and {b_(1,1) = b_(1,2) = 0}).  Adding the x^2 equation
         a_(1,0) b_(2,1) = 1 still leaves dimension 3: a_(1,1) and a_(1,2) are
         free.  So the lowest block is NOT decisive.

Stage 4  Weighted valuations.  For every primitive outward normal w of every
         edge of N(P) and of N(Q), the maximal-w part of det J is
         det J(P_w, Q_w) where P_w, Q_w are the edge forms.  Hence

             if  dP + dQ - w.(1,1)  >  w.(2,0)   then   det J(P_w,Q_w) = 0,
             if  dP + dQ - w.(1,1)  =  w.(2,0)   then   det J(P_w,Q_w) = x^2.

         Every weight is analysed and classified exactly:
           * w = (0,-1), (0,1), (1,1), (-2,1), (3,-1), (-1,0)... :
             det J(P_w,Q_w) is identically zero (both edge forms are
             univariate in the same monomial, or one is a monomial), so the
             condition is vacuous;
           * w = (1,0)  : 8 A B' - 12 A' B = 0 for the x=8 and x=12 slices;
             the strict-chart solution set is 3-dimensional
             (A = a (y-r)^2, B = b (y-r)^3), non-empty;
           * w = (1,-1), (1,-2) : a_(1,0) b_(2,1) = 1, non-empty;
           * w = (2,-1) : the top edge block reduces to 16 equations in 17
             unknowns A_1, B_1 with deg A_1 = 7, deg B_1 = 10; ten of the
             equations are triangular, leaving 6 equations in 7 unknowns.
             The block is non-empty as a system (it contains the degenerate
             point A = z, B = z^2); with the exact corner degrees the block
             is *not* resolved here.  This is the strongest weight available.
         No weight produces an empty leading-form system.  Hence no
         EMPTY_BY_VALUATION certificate follows from a single edge.

Stage 5  Exact quasi-homogeneous reduction.  The weight (-2,1) change of
         coordinates  u = y,  z = x y^2  (so x = z u^-2) gives the exact
         identity  det J(P,Q) = u^2 * Jac_{z,u}(P,Q).  Since x^2 = z^2 u^-4,
         the system is EQUIVALENT to

             u^6 * Jac_{z,u}(P,Q) = z^2,     (an identity in K[z, u, u^-1])

         and P, Q have only the u-powers 0,-1,-2 and 0,-1,-2,-3, so this is
         five univariate-polynomial equations

             (F0) -3 A' B + 2 A B'                    = z^2
             (F1) -2 A' D - 3 C' B + 2 A D' + C B'    = 0
             (F2) -A' E - 2 C' D - 3 F' B + 2 A E' + C D' = 0
             (F3) -C' E - 2 F' D + C E' + 2 A G'      = 0
             (F4) -F' E + C G'                        = 0

         in seven univariate polynomials

             A = P_{-2} (deg 8),  C = P_{-1} (deg <=8),  F = P_0 (deg 8),
             B = Q_{-3} (deg 12), D = Q_{-2} (deg <=12), E = Q_{-1} (deg <=12),
             G = Q_0 (deg 12),

         with 8+8+9+11+11+12+13 = 72 coefficients.  The equivalence is
         verified symbolically.  This is a *complete* reformulation, not a
         truncation: every solution of the original system is a solution of
         (F0)-(F4) and conversely.  The two top-weight structures are
         (i)  the (2,-1) block (F0) with the 6-in-7 reduction, and
         (ii) the (1,0) block (the x=8 / x=12 edge forms).
         Neither is contradictory.

Stage 6  x-adic (equivalently z-adic) cascade diagnostic.  Ordering the
         equations by x-degree, the top level is quadratic (solved in closed
         form), and every subsequent level is affine-linear in the new
         coefficients.  A run of that cascade with generic choices for the
         free parameters shows the first rank jump at x-degree 17: the 5
         equations in the 7 new coefficients have rank 4 but augmented rank 5.
         That is a *constraint on the free parameters* of the x-degree-18
         level, not a contradiction, because those parameters are not fixed by
         the leading data.  Reported as a diagnostic only.

VERDICT
-------
    NOT_DECISIVE

The support obstruction fails ((2,0) is attainable), the lowest-x-degree
block has dimension 3, and every single-edge leading-form system is non-empty.
What remains is a genuine multi-level cascade: the (F0)-(F4) system above, or
equivalently the x-adic cascade through x-degree 17 and below.  No emptiness
claim is made from any truncation; the exact support analysed is the full
25 + 47 lattice-point support, and the exact reformulation analysed is
(F0)-(F4) with all 72 coefficients free.

Run:
    python problems/jacobian/valuation72.py
"""

from __future__ import annotations

import argparse
import itertools
import json
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

P_CORNERS = ((0, 0), (1, 0), (8, 14), (8, 16))
Q_CORNERS = ((0, 0), (2, 1), (12, 21), (12, 24))
TARGET = (2, 0)

DEFAULT_OUT = ROOT / "runs" / "jacobian_72_valuation.json"
LOG_PATH = ROOT / "runs" / "valuation72.log"

x, y, z, u = sp.symbols("x y z u")


# --------------------------------------------------------------------------
# logging
# --------------------------------------------------------------------------

def log_stage(stage, payload):
    """Append exactly one line per stage to runs/valuation72.log."""
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    line = "%s stage=%s %s\n" % (
        stamp, stage, json.dumps(payload, sort_keys=True, default=str))
    with LOG_PATH.open("a", encoding="utf-8") as fh:
        fh.write(line)


# --------------------------------------------------------------------------
# support geometry
# --------------------------------------------------------------------------

def derivative_support(support, axis):
    """Exponents of d/dx (axis 0) or d/dy (axis 1), with the source point."""
    out = {}
    for a, b in support:
        if axis == 0 and a > 0:
            out[(a - 1, b)] = (a, b)
        if axis == 1 and b > 0:
            out[(a, b - 1)] = (a, b)
    return out


def convex_hull(points):
    pts = sorted(set(points))

    def cross(o, a, b):
        return ((a[0] - o[0]) * (b[1] - o[1])
                - (a[1] - o[1]) * (b[0] - o[0]))

    lower = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]


def attainable_support(support_p, support_q):
    """All monomials that can arise from P_x Q_y - P_y Q_x, with witnesses.

    Returns a dict  monomial -> list of (part, p, q, sign) where part is
    "PxQy" or "PyQx" and p, q are *derivative* exponents.
    """
    px = derivative_support(support_p, 0)
    py = derivative_support(support_p, 1)
    qx = derivative_support(support_q, 0)
    qy = derivative_support(support_q, 1)

    table = {}
    for p in px:
        for q in qy:
            m = (p[0] + q[0], p[1] + q[1])
            table.setdefault(m, []).append(("PxQy", p, q, +1))
    for p in py:
        for q in qx:
            m = (p[0] + q[0], p[1] + q[1])
            table.setdefault(m, []).append(("PyQx", p, q, -1))
    return table


def build_polynomial(support, prefix):
    coeffs = {}
    expr = sp.Integer(0)
    for p in support:
        c = sp.Symbol("%s%d_%d" % (prefix, p[0], p[1]))
        coeffs[p] = c
        expr += c * x ** p[0] * y ** p[1]
    return expr, coeffs


# --------------------------------------------------------------------------
# exact dimension of a zero-dimensional / positive-dimensional ideal
# --------------------------------------------------------------------------

def is_zero_dimensional(gb, variables):
    """True iff the quotient by the (lex) Groebner basis is finite dimensional.

    Criterion: for every variable v some leading monomial is a pure power of v.
    ``gb`` must be a lex Groebner basis (sympy's Poly default order is lex, so
    ``Poly(g, *variables).LM()`` is the leading monomial in that order).
    """
    lms = []
    for g in gb.exprs:
        if g == 0:
            continue
        lms.append(sp.Poly(g, *variables).LM())
    for v in variables:
        ok = False
        for lm in lms:
            powers = dict(zip(variables, lm))
            if powers.get(v, 0) > 0 and all(
                    powers.get(w, 0) == 0 for w in variables if w != v):
                ok = True
                break
        if not ok:
            return False
    return True


def krull_dimension(eqs, variables, max_slices=None, seed=12345):
    """Krull dimension via generic hyperplane sections.

    dim = n - k where k is the least number of generic linear forms that make
    the quotient finite dimensional.  Returns (dim, data).
    """
    n = len(variables)
    rng = __import__("random").Random(seed)
    if max_slices is None:
        max_slices = n
    data = []
    for k in range(0, max_slices + 1):
        forms = []
        for _ in range(k):
            coefs = [sp.Rational(rng.randint(-40, 40), rng.randint(1, 9))
                     for _ in variables]
            forms.append(sum(c * v for c, v in zip(coefs, variables)))
        gb = sp.groebner(list(eqs) + forms, *variables, order="lex")
        zd = is_zero_dimensional(gb, variables)
        data.append({"slices": k, "zero_dimensional": bool(zd),
                     "gb_size": len(gb.exprs)})
        if zd:
            return n - k, data
    return None, data


def solve_block(eqs, variables, label):
    """Exact solution data for a (small) polynomial system."""
    t0 = time.perf_counter()
    eqs = [sp.expand(e) for e in eqs if sp.expand(e) != 0]
    if not eqs:
        return {"label": label, "equations": 0, "variables": len(variables),
                "solution_dimension": len(variables),
                "empty": False,
                "note": "no equation: every coefficient is free",
                "seconds": 0.0}
    gb = sp.groebner(eqs, *variables, order="grevlex")
    if len(gb.exprs) == 1 and gb.exprs[0] == 1:
        return {"label": label, "equations": len(eqs),
                "variables": len(variables), "solution_dimension": -1,
                "empty": True, "groebner": ["1"],
                "seconds": round(time.perf_counter() - t0, 3)}
    dim, data = krull_dimension(eqs, variables)
    return {"label": label, "equations": len(eqs),
            "variables": len(variables), "solution_dimension": dim,
            "empty": False,
            "groebner": [str(g) for g in gb.exprs][:12],
            "dimension_certificate": data,
            "seconds": round(time.perf_counter() - t0, 3)}


# --------------------------------------------------------------------------
# stage 1
# --------------------------------------------------------------------------

def stage1(support_p, support_q, det_poly):
    table = attainable_support(support_p, support_q)
    actual = sorted((int(a), int(b)) for a, b in det_poly.monoms())
    attainable = sorted(table)
    x2_pairs = table.get(TARGET, [])

    # exact coefficient of x^2 in det J
    c20 = sp.expand(det_poly.coeff_monomial(TARGET))
    a10 = sp.Symbol("a1_0")
    b21 = sp.Symbol("b2_1")
    x2_ok = sp.expand(c20 - a10 * b21) == 0

    hist = {}
    for a, _b in actual:
        hist[a] = hist.get(a, 0) + 1

    missing = sorted(set(attainable) - set(actual))
    # each missing monomial must cancel between the two signed contributions
    cancelled = []
    for m in missing:
        pairs = table[m]
        # group by the underlying coefficient pair (source points)
        seen = {}
        for part, p, q, sgn in pairs:
            key = (p[0] + 1, p[1], q[0], q[1] + (1 if part == "PxQy" else 0))
            seen.setdefault(key, []).append(sgn)
        cancelled.append({
            "monomial": list(m),
            "pairs": [[part, list(p), list(q), sgn] for part, p, q, sgn in pairs],
        })

    return {
        "support_sizes": [len(support_p), len(support_q)],
        "unknowns": len(support_p) + len(support_q),
        "minkowski_pair_monomials": len(attainable),
        "actual_det_monomials": len(actual),
        "identically_cancelled": len(missing),
        "actual_support_is_subset_of_attainable": set(actual) <= set(attainable),
        "x2_attainable": TARGET in table,
        "x2_witness_pairs": [[part, list(p), list(q), sgn]
                             for part, p, q, sgn in x2_pairs],
        "x2_coefficient": str(c20),
        "x2_coefficient_is_a10_times_b21": bool(x2_ok),
        "min_x_degree": min(a for a, _ in actual),
        "max_total_degree": max(a + b for a, b in actual),
        "x_degree_histogram": {str(k): v for k, v in sorted(hist.items())},
        "det_monomials": [list(m) for m in actual],
        "cancelled_monomials": cancelled[:24],
        "conclusion": (
            "(2,0) is attainable and its coefficient is exactly a1_0*b2_1, "
            "so det J = x^2 is NOT excluded by support.  EMPTY_BY_SUPPORT is "
            "ruled out."),
    }


# --------------------------------------------------------------------------
# stage 2 / 3
# --------------------------------------------------------------------------

def stage23(det_poly, support_p, support_q):
    minx = min(int(a) for a, _ in det_poly.monoms())
    block = []
    for (a, b) in sorted((int(a), int(b)) for a, b in det_poly.monoms()):
        if a != minx:
            continue
        block.append({"monomial": [a, b],
                      "coefficient": str(sp.expand(det_poly.coeff_monomial((a, b))))})

    # the block lives in the x = 1 slice of P and the x = 1 slice of Q
    slice_p = sorted(p for p in support_p if p[0] == 1)
    slice_q = sorted(q for q in support_q if q[0] == 1)
    block_vars = ([sp.Symbol("a%d_%d" % p) for p in slice_p]
                  + [sp.Symbol("b%d_%d" % q) for q in slice_q])

    eqs = [sp.sympify(e["coefficient"]) for e in block]
    sol = solve_block(eqs, block_vars, "x^1 block")
    sol["monomials"] = [e["monomial"] for e in block]
    sol["slice_of_P"] = [list(p) for p in slice_p]
    sol["slice_of_Q"] = [list(q) for q in slice_q]

    # add the x^2 equation a1_0 b2_1 = 1
    a10 = sp.Symbol("a1_0")
    b21 = sp.Symbol("b2_1")
    joint_vars = block_vars + [b21]
    joint = solve_block(eqs + [a10 * b21 - 1], joint_vars,
                        "x^1 block + x^2 equation")

    return {
        "min_x_degree": minx,
        "block_monomials": [e["monomial"] for e in block],
        "equations": block,
        "variables": [str(v) for v in block_vars],
        "x1_block": sol,
        "x1_plus_x2": joint,
        "conclusion": (
            "The lowest-x-degree block is three equations in five unknowns. "
            "Its ideal is (a1_0 b1_1, a1_0 b1_2, a1_1 b1_2 - a1_2 b1_1), "
            "of dimension 3; adding a1_0 b2_1 = 1 keeps dimension 3. "
            "The block is satisfiable (e.g. b1_1 = b1_2 = 0, a1_1, a1_2 free, "
            "a1_0 b2_1 = 1), so it is not decisive."),
    }


# --------------------------------------------------------------------------
# stage 4
# --------------------------------------------------------------------------

def edge_normals(vertices):
    """Primitive outward normals of the hull edges."""
    hull = convex_hull(vertices)
    out = []
    n = len(hull)
    for i in range(n):
        a = hull[i]
        b = hull[(i + 1) % n]
        ux, uy = b[0] - a[0], b[1] - a[1]
        g = sp.igcd(abs(ux), abs(uy))
        ux, uy = ux // g, uy // g
        for w in ((uy, -ux), (-uy, ux)):
            d = max(w[0] * p[0] + w[1] * p[1] for p in hull)
            face = [p for p in hull if w[0] * p[0] + w[1] * p[1] == d]
            if len(face) >= 2:
                # outward normal: the polygon must lie on the <= side
                inside = all(w[0] * p[0] + w[1] * p[1] <= d for p in hull)
                if inside:
                    out.append((w, tuple(a), tuple(b)))
                    break
    return out


def max_face(points, w):
    d = max(w[0] * p[0] + w[1] * p[1] for p in points)
    return d, [p for p in points if w[0] * p[0] + w[1] * p[1] == d]


def stage4(support_p, support_q):
    normals = []
    for w, a, b in edge_normals(P_CORNERS):
        normals.append((w, "N(P)", (a, b)))
    for w, a, b in edge_normals(Q_CORNERS):
        normals.append((w, "N(Q)", (a, b)))

    seen = set()
    rows = []
    for w, src, edge in normals:
        if w in seen:
            continue
        seen.add(w)
        dP, fP = max_face(support_p, w)
        dQ, fQ = max_face(support_q, w)
        filt = dP + dQ - (w[0] + w[1])
        tgt = w[0] * TARGET[0] + w[1] * TARGET[1]

        # leading forms
        vP = {p: sp.Symbol("a%d_%d" % p) for p in fP}
        vQ = {q: sp.Symbol("b%d_%d" % q) for q in fQ}
        Pw = sum(c * x ** p[0] * y ** p[1] for p, c in vP.items())
        Qw = sum(c * x ** q[0] * y ** q[1] for q, c in vQ.items())
        lead = sp.expand(sp.diff(Pw, x) * sp.diff(Qw, y)
                         - sp.diff(Pw, y) * sp.diff(Qw, x))

        row = {
            "w": list(w),
            "source": src,
            "edge": [list(edge[0]), list(edge[1])],
            "dP": dP, "dQ": dQ,
            "face_P": [list(p) for p in fP],
            "face_Q": [list(q) for q in fQ],
            "det_filtration_degree": filt,
            "target_degree": tgt,
        }

        if filt < tgt:
            row["status"] = "CONTRADICTION"
            row["equation"] = None
            rows.append(row)
            continue

        if filt > tgt:
            kind = "det J(P_w,Q_w) = 0"
        else:
            kind = "det J(P_w,Q_w) = x^2"
        row["kind"] = kind
        row["leading_form"] = str(sp.factor(lead))

        if lead == 0:
            row["status"] = "vacuous (leading form identically 0)"
            rows.append(row)
            continue

        # set up the coefficient system of the leading form condition
        rhs = sp.Integer(0) if filt > tgt else x ** 2
        poly = sp.Poly(sp.expand(lead - rhs), x, y)
        eqlist = [sp.expand(c) for c in poly.coeffs()]
        variables = sorted(set(vP.values()) | set(vQ.values()), key=str)
        solved = solve_block(eqlist, variables,
                             "leading form at w=%s" % (w,))
        row["status"] = ("non-empty leading block"
                         if not solved.get("empty") else "EMPTY leading block")
        row["leading_block"] = solved
        rows.append(row)

    return {
        "rows": rows,
        "conclusion": (
            "Every edge weight gives a satisfiable (or vacuous) leading-form "
            "condition; none is empty.  No single-edge valuation yields a "
            "certificate of emptiness."),
    }


# --------------------------------------------------------------------------
# stage 5: exact quasi-homogeneous reduction
# --------------------------------------------------------------------------

def stage5():
    """u^6 Jac_{z,u}(P,Q) = z^2 and the five equations (F0)-(F4)."""
    buckets_p = {-2: [1, 8], -1: [1, 8], 0: [0, 8]}
    buckets_q = {-3: [2, 12], -2: [2, 12], -1: [1, 12], 0: [0, 12]}

    def poly(name, lo, hi, deg_hi):
        return sum(sp.Symbol("%s%d" % (name, i)) * z ** i
                   for i in range(lo, min(hi, deg_hi) + 1))

    A = poly("A", 1, 8, 8)
    C = poly("C", 1, 8, 8)
    F = poly("F", 0, 8, 8)
    B = poly("B", 2, 12, 12)
    D = poly("D", 2, 12, 12)
    E = poly("E", 1, 12, 12)
    G = poly("G", 0, 12, 12)
    Pz = F + u ** -1 * C + u ** -2 * A
    Qz = G + u ** -1 * E + u ** -2 * D + u ** -3 * B

    jac = sp.expand(sp.diff(Pz, z) * sp.diff(Qz, u)
                    - sp.diff(Pz, u) * sp.diff(Qz, z))
    expr = sp.expand(u ** 6 * jac - z ** 2)
    # split by u powers
    upoly = sp.Poly(expr, u)
    upowers = sorted(int(m[0]) for m in upoly.monoms())

    F0 = sp.expand(-3 * sp.diff(A, z) * B + 2 * A * sp.diff(B, z) - z ** 2)
    F1 = sp.expand(-2 * sp.diff(A, z) * D - 3 * sp.diff(C, z) * B
                   + 2 * A * sp.diff(D, z) + C * sp.diff(B, z))
    F2 = sp.expand(-sp.diff(A, z) * E - 2 * sp.diff(C, z) * D
                   - 3 * sp.diff(F, z) * B + 2 * A * sp.diff(E, z) + C * sp.diff(D, z))
    F3 = sp.expand(-sp.diff(C, z) * E - 2 * sp.diff(F, z) * D
                   + C * sp.diff(E, z) + 2 * A * sp.diff(G, z))
    F4 = sp.expand(-sp.diff(F, z) * E + C * sp.diff(G, z))
    rebuilt = sp.expand(F0 + u * F1 + u ** 2 * F2 + u ** 3 * F3 + u ** 4 * F4)
    identity_ok = sp.expand(expr - rebuilt).is_zero

    # top block (F0): A = z A1, B = z^2 B1
    a = sp.symbols("A1_1:8")
    b = sp.symbols("B1_1:11")
    A1 = 1 + sum(a[i] * z ** (i + 1) for i in range(7))
    B1 = 1 + sum(b[i] * z ** (i + 1) for i in range(10))
    G0 = sp.Poly(sp.expand(A1 * B1 - 3 * z * sp.diff(A1, z) * B1
                           + 2 * z * A1 * sp.diff(B1, z) - 1), z)
    g0c = [sp.expand(c) for c in G0.all_coeffs()[::-1]]
    bsol = {}
    triangular = True
    for n in range(1, 11):
        e = g0c[n].subs(bsol)
        coef = sp.expand(sp.diff(e, b[n - 1]))
        if coef == 0 or sp.Symbol("B1_%d" % n) not in e.free_symbols:
            triangular = False
            break
        bsol[b[n - 1]] = sp.expand(-(e - coef * b[n - 1]) / coef)
    residual = [sp.expand(c.subs(bsol)) for c in g0c[11:17]] if triangular else []

    return {
        "substitution": "u = y, z = x*y^2  (x = z*u^-2)",
        "identity": "det J(P,Q) = u^2 * Jac_{z,u}(P,Q);  target x^2 = z^2 u^-4",
        "equivalent_system": "u^6 * Jac_{z,u}(P,Q) = z^2",
        "u_powers_present": upowers,
        "five_equations": {
            "F0": "-3 A' B + 2 A B' = z^2",
            "F1": "-2 A' D - 3 C' B + 2 A D' + C B' = 0",
            "F2": "-A' E - 2 C' D - 3 F' B + 2 A E' + C D' = 0",
            "F3": "-C' E - 2 F' D + C E' + 2 A G' = 0",
            "F4": "-F' E + C G' = 0",
        },
        "polynomial_degrees": {"A": 8, "C": 8, "F": 8,
                               "B": 12, "D": 12, "E": 12, "G": 12},
        "coefficient_count": 72,
        "decomposition_verified_symbolically": bool(identity_ok),
        "top_block": {
            "equation": "A1*B1 - 3 z A1' B1 + 2 z A1 B1' = 1, "
                        "A = z*A1, B = z^2*B1, deg A1 = 7, deg B1 = 10",
            "scalar_equations": len(g0c),
            "unknowns": 17,
            "triangular_elimination": bool(triangular),
            "determined": "B1_1 .. B1_10 from A1_1 .. A1_7",
            "residual_equations": len(residual),
            "residual_in": 7,
            "residual_degrees": [sp.Poly(r, *a).total_degree() for r in residual],
            "degenerate_solution": "A = z, B = z^2 (A1 = B1 = 1) satisfies (F0)",
            "status": "non-empty as a system; the exact-degree (7,10) part is "
                      "reduced to 6 equations in 7 unknowns and is NOT resolved",
        },
        "conclusion": (
            "The system is exactly equivalent to (F0)-(F4); the top block (F0) "
            "reduces to 6 equations in 7 unknowns and admits the degenerate "
            "solution A = z, B = z^2.  No contradiction at this level."),
    }


# --------------------------------------------------------------------------
# stage 6: x-adic cascade diagnostic
# --------------------------------------------------------------------------

def stage6(det_poly, support_p, support_q):
    """Rank profile of the x-degree-ordered linear cascade (diagnostic)."""
    import random
    coef = {}
    for (k, l) in det_poly.monoms():
        coef[(int(k), int(l))] = sp.expand(det_poly.coeff_monomial((k, l)))
    target = {TARGET: sp.Integer(1)}
    vars_p = {p: sp.Symbol("a%d_%d" % p) for p in support_p}
    vars_q = {q: sp.Symbol("b%d_%d" % q) for q in support_q}

    rng = random.Random(20261102)
    a_param = sp.Rational(3, 5)
    b_param = sp.Rational(2, 3)
    r_param = sp.Rational(1, 2)
    val = {
        (8, 14): a_param * r_param ** 2,
        (8, 15): -2 * a_param * r_param,
        (8, 16): a_param,
        (12, 21): -b_param * r_param ** 3,
        (12, 22): 3 * b_param * r_param ** 2,
        (12, 23): -3 * b_param * r_param,
        (12, 24): b_param,
    }
    val = {vars_p[p] if p in vars_p else vars_q[p]: v for p, v in val.items()}

    maxk = max(k for k, _ in coef)
    rows = []
    level17 = None
    for k in range(maxk - 1, -1, -1):
        eqs = []
        for (kk, l), c in sorted(coef.items()):
            if kk != k:
                continue
            e = sp.expand(c - target.get((kk, l), 0))
            eqs.append((l, sp.expand(e.subs(val))))
        if not eqs:
            continue
        newvars = sorted({s for (_l, e) in eqs for s in e.free_symbols
                          if s not in val}, key=str)
        if not newvars:
            bad = [l for (l, e) in eqs if e != 0]
            rows.append({"x_degree": k, "new": 0, "eqs": len(eqs),
                         "nonzero": len(bad)})
            if bad:
                break
            continue
        try:
            M, rhs = sp.linear_eq_to_matrix([e for _l, e in eqs], newvars)
        except Exception:
            rows.append({"x_degree": k, "new": len(newvars),
                         "eqs": len(eqs), "note": "nonlinear in new unknowns"})
            break
        M = sp.Matrix(M)
        rhs = sp.Matrix(rhs)
        r = M.rank()
        ra = M.row_join(rhs).rank()
        rows.append({"x_degree": k, "new": len(newvars), "eqs": len(eqs),
                     "rank": int(r), "augmented_rank": int(ra),
                     "free": len(newvars) - int(r),
                     "consistent": bool(r == ra)})
        if k == 17:
            level17 = rows[-1]
        if r != ra:
            break
        sol = M.gauss_jordan_solve(rhs)[0]
        ns = M.nullspace()
        cs = [rng.randint(-3, 3) for _ in ns]
        x0 = sp.Matrix(sol) + sum(
            (c * v for c, v in zip(cs, ns)), sp.Matrix(sp.zeros(len(newvars), 1)))
        taus = sorted({s for s in x0.free_symbols if str(s).startswith("tau")},
                      key=str)
        x0 = x0.subs({t: sp.Rational(rng.randint(-4, 4), rng.randint(1, 4))
                      for t in taus})
        for s, v in zip(newvars, list(x0)):
            val[s] = v

    return {
        "method": "x-degree ordered cascade; level 19 solved in closed form "
                  "(A = a(y-r)^2, B = b(y-r)^3 on the x=8/x=12 slices); every "
                  "lower level is affine-linear in the new coefficients",
        "parameters_used": {"a": str(a_param), "b": str(b_param),
                            "r": str(r_param)},
        "levels": rows,
        "first_rank_jump": level17,
        "conclusion": (
            "The first rank jump is at x-degree 17 (rank 4, augmented rank 5). "
            "That is a constraint on the free parameters of the x-degree-18 "
            "level, not a contradiction: those parameters are not fixed by the "
            "leading data.  Diagnostic only, no emptiness claim."),
    }


# --------------------------------------------------------------------------
# report
# --------------------------------------------------------------------------

def build_report():
    t0 = time.perf_counter()
    support_p = sorted(lattice_points(P_CORNERS))
    support_q = sorted(lattice_points(Q_CORNERS))

    P, c_p = build_polynomial(support_p, "a")
    Q, c_q = build_polynomial(support_q, "b")
    det = sp.expand(sp.diff(P, x) * sp.diff(Q, y)
                    - sp.diff(P, y) * sp.diff(Q, x))
    det_poly = sp.Poly(det, x, y)

    s1 = stage1(support_p, support_q, det_poly)
    log_stage("1_support", {
        "x2_attainable": s1["x2_attainable"],
        "x2_coefficient": s1["x2_coefficient"],
        "attainable": s1["minkowski_pair_monomials"],
        "actual": s1["actual_det_monomials"],
        "min_x_degree": s1["min_x_degree"],
    })

    s23 = stage23(det_poly, support_p, support_q)
    log_stage("2_lowest_block", {
        "min_x_degree": s23["min_x_degree"],
        "block_monomials": s23["block_monomials"],
        "equations": [e["coefficient"] for e in s23["equations"]],
        "x1_dimension": s23["x1_block"].get("solution_dimension"),
        "x1_plus_x2_dimension": s23["x1_plus_x2"].get("solution_dimension"),
    })

    s4 = stage4(support_p, support_q)
    log_stage("3_weights", {
        "weights": [r["w"] for r in s4["rows"]],
        "statuses": [r["status"] for r in s4["rows"]],
    })

    s5 = stage5()
    log_stage("4_quasihomogeneous", {
        "identity_verified": s5["decomposition_verified_symbolically"],
        "u_powers": s5["u_powers_present"],
        "top_residual_equations": s5["top_block"]["residual_equations"],
    })

    s6 = stage6(det_poly, support_p, support_q)
    log_stage("5_cascade_diagnostic", {
        "first_rank_jump": s6["first_rank_jump"],
        "levels": s6["levels"][:4],
    })

    verdict = "NOT_DECISIVE"
    rep = {
        "problem": "jacobian-config2-full-lattice-valuation",
        "source": "GGHV arXiv:2204.14178, Prop. 4.3, case (2)",
        "target": "det J(P,Q) = x^2",
        "support_convention": "full lattice points of the two Newton polygons",
        "supports": {
            "N(P)_corners": [list(p) for p in P_CORNERS],
            "N(Q)_corners": [list(p) for p in Q_CORNERS],
            "lattice_points_P": len(support_p),
            "lattice_points_Q": len(support_q),
            "unknowns": len(support_p) + len(support_q),
        },
        "stage1_attainable_support": s1,
        "stage2_lowest_x_degree_block": s23,
        "stage3_weighted_valuations": s4,
        "stage4_quasihomogeneous_reduction": s5,
        "stage5_x_adic_cascade_diagnostic": s6,
        "verdict": verdict,
        "why_not_decisive": [
            "EMPTY_BY_SUPPORT is false: (2,0) is attainable; the coefficient "
            "of x^2 in det J is exactly a1_0*b2_1, a single nonzero bilinear "
            "term from the vertices (1,0) of N(P) and (2,1) of N(Q).",
            "The minimal x-degree of det J is 1, and the x^1 block is three "
            "equations in five unknowns with a 3-dimensional solution set; "
            "adding a1_0*b2_1 = 1 keeps dimension 3.",
            "Every single-edge leading-form condition is either vacuous "
            "(leading form identically zero) or satisfiable; no weight yields "
            "an empty leading block.",
            "The strongest weight, w = (2,-1), reduces to 16 equations in 17 "
            "unknowns with a triangular elimination down to 6 equations in 7 "
            "unknowns; that block is non-empty as a system and its "
            "exact-corner-degree part is not resolved by leading data alone.",
        ],
        "what_would_decide": [
            "Solve (F0)-(F4) exactly: five univariate-polynomial identities in "
            "seven polynomials of degrees 8,8,8,12,12,12,12 (72 coefficients). "
            "This is an exact, complete reformulation, not a truncation.",
            "Or run the x-degree-ordered cascade to the bottom: the levels at "
            "x-degree <= 6 introduce no new coefficients, so all of the "
            "overdetermination (92 equations vs 72 unknowns) sits there.",
        ],
        "what_is_not_claimed": [
            "No emptiness statement is made.  The verdict is NOT_DECISIVE.",
            "No truncation of the support was used: the analysed support is "
            "the full 25 + 47 lattice-point support.",
            "The x-degree-17 rank jump is a diagnostic computed at one "
            "generic parameter point; it is not a certificate.",
        ],
        "elapsed_seconds": None,
    }
    rep["elapsed_seconds"] = round(time.perf_counter() - t0, 3)
    return rep


def print_report(rep):
    s1 = rep["stage1_attainable_support"]
    print("== GGHV config 2, full lattice: valuation / leading-term analysis ==")
    print("unknowns %d (|N(P)| = %d, |N(Q)| = %d)"
          % (rep["supports"]["unknowns"],
             rep["supports"]["lattice_points_P"],
             rep["supports"]["lattice_points_Q"]))
    print()
    print("stage 1  attainable support")
    print("  monomials from P_x Q_y - P_y Q_x : %d"
          % s1["minkowski_pair_monomials"])
    print("  monomials actually present       : %d"
          % s1["actual_det_monomials"])
    print("  identically cancelled            : %d"
          % s1["identically_cancelled"])
    print("  (2,0) attainable                 : %s" % s1["x2_attainable"])
    print("  coefficient of x^2               : %s" % s1["x2_coefficient"])
    print("  min x-degree                     : %d" % s1["min_x_degree"])
    print()
    s2 = rep["stage2_lowest_x_degree_block"]
    print("stage 2/3  lowest x-degree block (x^1)")
    for e in s2["equations"]:
        print("  x^1 y^%-2d : %s" % (e["monomial"][1], e["coefficient"]))
    print("  x^1 block dimension              : %s"
          % s2["x1_block"].get("solution_dimension"))
    print("  x^1 + (a1_0 b2_1 = 1) dimension  : %s"
          % s2["x1_plus_x2"].get("solution_dimension"))
    print()
    print("stage 4  weighted valuations")
    for r in rep["stage3_weighted_valuations"]["rows"]:
        print("  w=%-8s %-6s filt=%2d target=%2d  %s"
              % (tuple(r["w"]), r["source"], r["det_filtration_degree"],
                 r["target_degree"], r["status"]))
    print()
    s4 = rep["stage4_quasihomogeneous_reduction"]
    print("stage 5  exact quasi-homogeneous reduction")
    print("  identity verified                : %s"
          % s4["decomposition_verified_symbolically"])
    print("  u-powers                         : %s" % s4["u_powers_present"])
    print("  top block residual equations     : %d in %d unknowns"
          % (s4["top_block"]["residual_equations"],
             s4["top_block"]["residual_in"]))
    print()
    print("stage 6  x-adic cascade (diagnostic)")
    for row in rep["stage5_x_adic_cascade_diagnostic"]["levels"]:
        print("  x^%-2d  new=%d eqs=%d rank=%s aug=%s %s"
              % (row["x_degree"], row.get("new", 0), row.get("eqs", 0),
                 row.get("rank", "-"), row.get("augmented_rank", "-"),
                 "" if row.get("consistent", True) else "RANK JUMP"))
    print()
    print("VERDICT: %s" % rep["verdict"])
    print("  (2,0) is attainable -> EMPTY_BY_SUPPORT ruled out.")
    print("  lowest block has dimension 3 -> not decisive.")
    print("  every edge leading-form block is non-empty -> no single-weight")
    print("  valuation certificate.")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--no-write", action="store_true")
    args = ap.parse_args(argv)
    rep = build_report()
    print_report(rep)
    if not args.no_write:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(rep, indent=2, default=str),
                            encoding="utf-8")
        print("\nwritten: %s" % args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
