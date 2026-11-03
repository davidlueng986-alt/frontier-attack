#!/usr/bin/env python3
"""Approximate-root reduction at the (75,125) Jacobian corner F_2(3,5), j=1.

WHAT THIS SCRIPT DOES
---------------------
It implements Moh's approximate-root construction
(``doi:10.1515/crll.1983.340.140``; the algebraic prototype is Abhyankar,
TIFR Lectures 57, chapter 1 sections 1-4) at the corner

    Delta' = conv{(0,0),(3,0),(4,1),(0,5)}
    N(P) = 3*Delta' = conv{(0,0),(9,0),(12,3),(0,15)},   deg P = 15
    N(Q) = 5*Delta' = conv{(0,0),(15,0),(20,5),(0,25)},  deg Q = 25
    det J(P,Q) = x^2

and it reports the size of the system that the construction leaves.

THE CLASSICAL CONSTRUCTION (extracted, not assumed)
---------------------------------------------------
Abhyankar TIFR 57, (4.3): let f in R[Y] be monic of degree n = m*d and let
d be a unit in R.  Write f(1/theta)^{-1/d} = theta^{-m} + a_1 theta^{-m+1} +
... Then the d-th approximate root of f is

    App^d_Y(f) = Y^m + a_1 Y^{m-1} + ... + a_m,

i.e. the polynomial part of f^{1/d}; equivalently the unique monic g of
degree n/d with deg_Y(f - g^d) < n - n/d.

Moh, section 1 (p. 148) repeats this definition verbatim, and section 2
uses it for a *pair*: given f, g monic in y with deg_y f = m, deg_y g = n,
let F(f,g) = 0 be the defining equation of f over k[x]((g^{-1})).  Then
T_r := the d_r-th approximate root of F with f as the variable satisfies

    T_r(f,g) = f^{n/d_r} + lower terms in f,      deg_y T_r = -mu_r,

and {n, mu_1, ..., mu_r} is the characteristic data of the pair, with
gcd(n, mu_1, ..., mu_r) = d_{r+1}.  The T_r are in k[x][f,g] ("specialized
approximate roots", Prop. 2.1) and their y-degrees carry the information.

WHAT IS PROVED HERE, AND WHAT IS A MODELLING CHOICE
---------------------------------------------------
Proved (exact, over Q):

  (P1) The degree-38 part of det J is J(P_15,Q_25).  Writing
       P_15 = y^15 p(t), Q_25 = y^25 q(t), t = x/y, one has
       J(P_15,Q_25) = y^38 (25 p' q - 15 p q'), so it vanishes iff
       p^5 = C q^3; unique factorization then gives p = c r^3, q = d r^5
       with deg r = 4.  Hence P_15 = c G^3 and Q_25 = d G^5 for ONE form
       G = y^5 r(x/y) of degree 5, supported on the (1,1) edge
       [(0,5),(4,1)] of Delta'.  G is the common approximate root: the
       3rd approximate root of P and the 5th of Q, in the graded sense.
  (P2) det J(c G^3, d G^5) = 0 identically, so the substitution
       P = G^3 + R, Q = G^5 + S (after the det-preserving rescaling
       gauge) makes the leading equation automatic.
  (P3) At the next order (degree 37) the equation is exactly
       3 J(G,Q_24) + 5 d G^2 J(P_14,G) = 0, whose solution space is
       {P_14 arbitrary (12 coefficients), Q_24 = (5d/3) G^2 P_14};
       dimension 12.  Verified by exact linear algebra over Q and by the
       closed form (the centraliser of a non-degenerate form under the
       Jacobian is C[G], and 5 does not divide 24).
  (P4) Every level d >= 3 is homogeneous in the corrections, so the chain
       can be consistent (the zero correction is always a solution) and
       never produces an "inconsistency" above degree 2; the only
       inhomogeneous equation is the target x^2 at degree 2.

Modelling choices (declared, not proved):

  (M1) The classical Moh approximate roots are y-adic: for this corner
       P has y-degree 15 and Q y-degree 25, so their 5-th approximate
       roots have y-degrees 3 and 5 respectively and are DIFFERENT
       polynomials.  There is no common y-adic root.  The object used
       here is the common *graded* root G of (P1) -- the leading-form
       (total-degree) approximate root, which is what "couples the
       edges" at this corner.  A y-adic implementation is included as a
       cross-check and is shown to be a triangular change of variables
       with no reduction (Definition only, no constraint).
  (M2) The chain is truncated at degree 37 (the first non-trivial order).
       Higher orders are counted, not solved.
  (M3) c = 1 is fixed by the gauge (P,Q) -> (lambda P, mu Q) with
       lambda mu = 1, which preserves det J = x^2 and both Newton polygons.
       c is the y^15 coefficient of P after G is taken monic, so c = 1 is a
       genuine gauge slice.  The ratio d (the y^25 coefficient of Q) is a
       modulus of the slice, NOT a gauge direction.  The computation is run
       at d = 1 and repeated at d = 2, 3, 5: d multiplies one block of the
       matrix columns, so rank and nullity do not depend on d != 0.

WHAT THE CONSTRUCTION DOES *NOT* DO
-----------------------------------
It does not reduce the corner to ten coefficients.  Moh's own Appendix I
records that the surviving configurations at degree <= 100 still carry
"3370 to 7328" coefficients: the approximate-root method reduces the
*root configuration* (a finite combinatorial list), not the coefficient
space.  At this corner the root configuration is already fixed by Delta',
so the gain is exactly the top-edge block: 34 coefficients -> 5
(4 shape parameters of G plus the modulus d).

Run:
    python problems/jacobian/approx_roots.py
"""

from __future__ import annotations

import json
import sys
import time
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

LOG_PATH = ROOT / "runs" / "approx_roots.log"
OUT_PATH = ROOT / "runs" / "jacobian_approx_roots.json"
MOH_TXT = ROOT / "runs" / "papers" / "moh1983.txt"
MOH_PDF_DIR = ROOT / "runs" / "papers" / "moh1983_libgen.pdf"
ABH_PDF = ROOT / "runs" / "papers" / "abhyankar_tifr57.pdf"

X, Y = sp.symbols("x y")
DELTA = ((0, 0), (3, 0), (4, 1), (0, 5))
TIME_BUDGET_S = 1200.0
T0 = time.time()


# --------------------------------------------------------------------------
# logging / small helpers
# --------------------------------------------------------------------------

def log_stage(stage: str, payload: dict) -> None:
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    line = "%s stage=%s %s\n" % (
        stamp, stage, json.dumps(payload, sort_keys=True, default=str))
    with LOG_PATH.open("a", encoding="utf-8") as fh:
        fh.write(line)


def elapsed() -> float:
    return time.time() - T0


def hull(points):
    pts = sorted(set(map(tuple, points)))

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

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


def lattice_points(poly):
    """All lattice points of a convex polygon given counter-clockwise."""
    xs = [p[0] for p in poly]
    ys = [p[1] for p in poly]
    out = []
    n = len(poly)
    for i in range(min(xs), max(xs) + 1):
        for j in range(min(ys), max(ys) + 1):
            ok = True
            for k in range(n):
                a, b = poly[k], poly[(k + 1) % n]
                if (b[0] - a[0]) * (j - a[1]) - (b[1] - a[1]) * (i - a[0]) < 0:
                    ok = False
                    break
            if ok:
                out.append((i, j))
    return out


def jac(f, g):
    return sp.expand(sp.diff(f, X) * sp.diff(g, Y) - sp.diff(f, Y) * sp.diff(g, X))


# --- exact rational linear algebra (sympy's rank() is far too slow here) ---

def rref(rows, ncols):
    M = [[Fraction(x) for x in row] for row in rows]
    pivots = []
    r = 0
    for c in range(ncols):
        p = None
        for i in range(r, len(M)):
            if M[i][c] != 0:
                p = i
                break
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        inv = M[r][c]
        M[r] = [x / inv for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        pivots.append(c)
        r += 1
        if r == len(M):
            break
    return M, pivots


def matrix_rank(rows, ncols):
    return len(rref(rows, ncols)[1])


def nullspace(rows, ncols):
    M, pivots = rref(rows, ncols)
    pivset = set(pivots)
    free = [c for c in range(ncols) if c not in pivset]
    basis = []
    for f in free:
        v = [Fraction(0)] * ncols
        v[f] = Fraction(1)
        for i, c in enumerate(pivots):
            v[c] = -M[i][f]
        basis.append(v)
    return basis


def to_sym(v):
    return sp.Rational(v.numerator, v.denominator)


def monomials(points):
    return [sp.Integer(1) * X**i * Y**j for (i, j) in points]


def coeff_vector(poly, degree):
    """Coefficients of x^i y^(degree-i), i = 0..degree."""
    P = sp.Poly(sp.expand(poly), X, Y)
    return [P.coeff_monomial(X**i * Y**(degree - i)) for i in range(degree + 1)]


# --------------------------------------------------------------------------
# stage 1 -- sources
# --------------------------------------------------------------------------

APPROX_ROOT_DEF_ABHYANKAR = (
    "Abhyankar, TIFR Lectures 57, (4.2)-(4.3): f in R[Y] monic of degree n = m*d, "
    "d a unit in R; write f^{1/d} = theta^{-m} + a_1 theta^{-m+1} + ...; the d-th "
    "approximate root is Y^m + a_1 Y^{m-1} + ... + a_m, i.e. the polynomial part of "
    "f^{1/d}, characterised by deg_Y(f - g^d) < n - n/d."
)

APPROX_ROOT_DEF_MOH = (
    "Moh 1983, section 1, p.148 (Definition): the d-th approximate root f_d(y) of a "
    "monic f(y) of degree n = m*d is y^m + a_1 y^{m-1} + ... + a_m, read off from the "
    "d-th root of f in R((theta)), theta = y^{-1}."
)

MOH_PAIR_USE = (
    "Moh 1983, section 2, pp.151-154: T_r is the d_r-th approximate root of the "
    "defining equation F(f,g)=0 over k[x]((g^{-1})), with f as the variable; "
    "T_r(f,g) = f^{n/d_r} + lower terms in f and deg_y T_r = -mu_r. The sequence "
    "{n, mu_1, ..., mu_r} is the characteristic data and gcd(n, mu_1, ..., mu_r) = "
    "d_{r+1}. Prop. 2.1: the specialisation psi: x -> 0 keeps T_r an approximate root "
    "(T*_r in k[f,g])."
)

MOH_COST = (
    "Moh 1983, Appendix I, p.143: for the configurations surviving at degree <= 100 "
    '"the number of coefficients of those possible counter-examples range from 3370 '
    'to 7328".'
)


def stage1_sources():
    """Confirm the sources and extract the construction verbatim where possible."""
    text = ""
    if MOH_TXT.exists():
        text = MOH_TXT.read_text(encoding="utf-8", errors="replace")

    def find_line(needles):
        for ln in text.splitlines():
            low = ln.lower()
            if all(n.lower() in low for n in needles):
                return ln.strip()
        return None

    extracted = {
        "moh_txt_present": MOH_TXT.exists(),
        "moh_pdf_present": any(Path(MOH_PDF_DIR).glob("*.pdf")),
        "abhyankar_pdf_present": ABH_PDF.exists(),
        "moh_definition_line": find_line(["d-th approximate root"]),
        "moh_specialized_line": find_line(["specialized"]),
        "moh_coefficient_line": find_line(["3370"]),
    }
    payload = {
        "construction": {
            "abhyankar_4_3": APPROX_ROOT_DEF_ABHYANKAR,
            "moh_section1": APPROX_ROOT_DEF_MOH,
            "moh_pair_section2": MOH_PAIR_USE,
            "moh_cost_appendixI": MOH_COST,
        },
        "local_sources": extracted,
    }
    log_stage("1_source", {
        "moh_txt": str(MOH_TXT),
        "moh_present": extracted["moh_txt_present"],
        "moh_definition_found": extracted["moh_definition_line"] is not None,
        "moh_pair_use_found": extracted["moh_specialized_line"] is not None,
        "moh_cost_found": extracted["moh_coefficient_line"] is not None,
        "abhyankar_present": extracted["abhyankar_pdf_present"],
        "variable": "y (classical); graded edge variable for the common root",
    })
    return payload, extracted


# --------------------------------------------------------------------------
# stage 2 -- the corner and the common approximate root
# --------------------------------------------------------------------------

def stage2_corner():
    np_ = hull([(3 * a, 3 * b) for a, b in DELTA])
    nq_ = hull([(5 * a, 5 * b) for a, b in DELTA])
    sp_pts = lattice_points(np_)
    sq_pts = lattice_points(nq_)

    top_p = [p for p in sp_pts if p[0] + p[1] == 15]
    top_q = [p for p in sq_pts if p[0] + p[1] == 25]

    # the common approximate root: degree 5, support on the (1,1) edge
    gsym = sp.symbols("g1:5")
    G = Y**5 + gsym[0] * X * Y**4 + gsym[1] * X**2 * Y**3 \
        + gsym[2] * X**3 * Y**2 + gsym[3] * X**4 * Y
    G3 = sp.expand(G**3)
    G5 = sp.expand(G**5)

    # (P2) the leading Jacobian vanishes identically
    leading_jac = sp.expand(jac(G3, G5))

    # support checks
    supp_g3 = set(sp.Poly(G3, X, Y).monoms())
    supp_g5 = set(sp.Poly(G5, X, Y).monoms())

    # (P1) sufficiency certificate: for a random r of degree 4, p = r^3,
    # q = r^5 satisfies the Wronskian identity 25 p' q - 15 p q' = 0.
    t = sp.symbols("t")
    rng = [sp.Integer(2), sp.Integer(-1), sp.Integer(3), sp.Integer(5), sp.Integer(1)]
    r = sum(c * t**k for k, c in enumerate(rng))
    p = sp.expand(r**3)
    q = sp.expand(r**5)
    wronskian = sp.expand(25 * sp.diff(p, t) * q - 15 * p * sp.diff(q, t))
    wronskian_perturbed = sp.expand(
        25 * sp.diff(p, t) * q - 15 * p * sp.diff(q + t**20, t))

    payload = {
        "delta_prime": [list(p) for p in DELTA],
        "N_P": [list(p) for p in np_],
        "N_Q": [list(p) for p in nq_],
        "lattice_P": len(sp_pts),
        "lattice_Q": len(sq_pts),
        "corner_unknowns": len(sp_pts) + len(sq_pts),
        "degree": [max(i + j for i, j in sp_pts), max(i + j for i, j in sq_pts)],
        "top_edge_P": len(top_p),
        "top_edge_Q": len(top_q),
        "top_block_unknowns": len(top_p) + len(top_q),
        "G_degree": 5,
        "G_support": [[0, 5], [1, 4], [2, 3], [3, 2], [4, 1]],
        "G_free_shape_parameters": 4,
        "top_block_after": 5,
        "modulus_d": 1,
        "leading_jacobian_identically_zero": leading_jac == 0,
        "G3_support_inside_NP": supp_g3 <= set(sp_pts),
        "G5_support_inside_NQ": supp_g5 <= set(sq_pts),
        "wronskian_identity_zero": wronskian == 0,
        "wronskian_perturbed_nonzero": wronskian_perturbed != 0,
        "substituted_unknowns": 5 + (len(sp_pts) - len(top_p)) + (len(sq_pts) - len(top_q)),
        "correction_P": len(sp_pts) - len(top_p),
        "correction_Q": len(sq_pts) - len(top_q),
    }
    log_stage("2_substitution", {
        "corner_unknowns": payload["corner_unknowns"],
        "top_block": payload["top_block_unknowns"],
        "root_unknowns_after": 5,
        "substituted_unknowns": payload["substituted_unknowns"],
        "leading_order_automatic": payload["leading_jacobian_identically_zero"],
        "root_degree": 5,
        "exponents": [3, 5],
    })
    return payload, np_, nq_, sp_pts, sq_pts, G.subs(dict(zip(gsym, [1, 2, 3, 5])))


# --------------------------------------------------------------------------
# stage 3 -- det J = x^2 at the first non-trivial order
# --------------------------------------------------------------------------

def degree_37_system(G, sp_pts, sq_pts, dval=1):
    """Exact degree-37 system in (P_14, Q_24); returns rank / nullity.

    The top block is normalised to P_15 = G^3, Q_25 = dval * G^5.  dval
    multiplies the Q_24 columns, so the rank may only drop at dval = 0.
    """
    P15 = sp.expand(G**3)
    Q25 = sp.expand(dval * G**5)
    s14 = [p for p in sp_pts if p[0] + p[1] == 14]
    s24 = [p for p in sq_pts if p[0] + p[1] == 24]

    p14 = sp.symbols("a0:%d" % len(s14))
    q24 = sp.symbols("b0:%d" % len(s24))
    m14 = monomials(s14)
    m24 = monomials(s24)

    # block-wise assembly: J(P15,Q24) + J(P14,Q25)
    jq = [sp.Poly(jac(P15, m), X, Y) for m in m24]
    jp = [sp.Poly(jac(m, Q25), X, Y) for m in m14]
    rows = []
    for i in range(38):
        row = []
        for P in jp:
            row.append(int(P.coeff_monomial(X**i * Y**(37 - i))))
        for P in jq:
            row.append(int(P.coeff_monomial(X**i * Y**(37 - i))))
        rows.append(row)
    M = rows
    ncols = len(m14) + len(m24)
    rank = matrix_rank(M, ncols)
    nullity = ncols - rank
    ns = nullspace(M, ncols)

    # closed form on a random combination of the nullspace basis
    comb = [Fraction((3 * k + 1) % 7 - 3) for k in range(len(ns))]
    if all(c == 0 for c in comb):
        comb = [Fraction(1)] * len(ns)
    P14 = sum(to_sym(c) * sum(to_sym(v[i]) * m14[i] for i in range(len(m14)))
              for c, v in zip(comb, ns))
    Q24 = sum(to_sym(c) * sum(to_sym(v[len(m14) + i]) * m24[i]
                              for i in range(len(m24)))
              for c, v in zip(comb, ns))
    closed_form = sp.expand(Q24 - sp.Rational(5 * dval, 3) * G**2 * P14)
    jac_identity = sp.expand(
        jac(P14, Q24) - sp.Rational(10 * dval, 3) * G * P14 * jac(P14, G))

    return {
        "unknowns": len(s14) + len(s24),
        "P14_unknowns": len(s14),
        "Q24_unknowns": len(s24),
        "equations": 38,
        "rank": int(rank),
        "dimension": int(nullity),
        "dval": dval,
        "closed_form_ok": closed_form == 0,
        "jacobian_identity_ok": jac_identity == 0,
    }, (P14, Q24, s14, s24)


def degree_36_compatibility(G, sp_pts, sq_pts, P14, Q24, s14, s24, dval=1):
    """Rank of the degree-36 map and the degree-26 compatibility conditions.

    The degree-36 part of det J is

        G^2 * [ 3 J(G,Q_23) + 5 d G^2 J(P_13,G) ] = -J(P_14,Q_24).

    The bracket is a degree-26 form, so the system is solvable iff
    (i) G^2 divides J(P_14,Q_24) and (ii) the quotient lies in the image
    of the degree-26 map (P_13,Q_23) -> 3J(G,Q_23) + 5dG^2J(P_13,G).
    """
    P15 = sp.expand(G**3)
    Q25 = sp.expand(dval * G**5)
    s13 = [p for p in sp_pts if p[0] + p[1] == 13]
    s23 = [p for p in sq_pts if p[0] + p[1] == 23]
    m13 = monomials(s13)
    m23 = monomials(s23)

    # degree-26 map: 3 J(G, Q_23) + 5 d G^2 J(P_13, G)
    jq26 = [sp.Poly(sp.expand(3 * jac(G, m)), X, Y) for m in m23]
    jp26 = [sp.Poly(sp.expand(5 * dval * G**2 * jac(m, G)), X, Y) for m in m13]
    rows26 = []
    for i in range(27):
        row = []
        for P in jp26:
            row.append(int(P.coeff_monomial(X**i * Y**(26 - i))))
        for P in jq26:
            row.append(int(P.coeff_monomial(X**i * Y**(26 - i))))
        rows26.append(row)
    ncols = len(m13) + len(m23)
    rank26 = matrix_rank(rows26, ncols)
    left = nullspace([list(col) for col in zip(*rows26)], len(rows26))

    # degree-36 map (for the kernel dimension); its rank equals rank26
    jq36 = [sp.Poly(jac(P15, m), X, Y) for m in m23]
    jp36 = [sp.Poly(jac(m, Q25), X, Y) for m in m13]
    rows36 = []
    for i in range(37):
        row = []
        for P in jp36:
            row.append(int(P.coeff_monomial(X**i * Y**(36 - i))))
        for P in jq36:
            row.append(int(P.coeff_monomial(X**i * Y**(36 - i))))
        rows36.append(row)
    rank36 = matrix_rank(rows36, ncols)

    # sample one level-37 solution: J(P14,Q24) = (10d/3) G P14 J(P14,G)
    rhs = sp.expand(-jac(P14, Q24))
    quo, rem = sp.div(rhs, G**2)
    divisible = sp.expand(rem) == 0
    residuals = None
    if divisible:
        fvec = coeff_vector(sp.expand(quo), 26)
        residuals = sum(1 for v in left
                        if sp.expand(sum(to_sym(v[i]) * c
                                         for i, c in enumerate(fvec))) != 0)

    # first non-zero branch of the divisibility condition G | P14*J(P14,G):
    # P14 = G*u with u of degree 9 and x-degree <= 7 (so G*u stays in S14)
    support_u = [(i, 9 - i) for i in range(0, 8)]
    u = sp.symbols("u0:%d" % len(support_u))
    prod = sp.expand(G * sum(c * X**i * Y**j for c, (i, j) in zip(u, support_u)))
    branch_ok = set(sp.Poly(prod, X, Y).monoms()) <= set(s14)

    return {
        "unknowns": len(s13) + len(s23),
        "equations": 37,
        "rank36": int(rank36),
        "rank26": int(rank26),
        "solution_dimension_if_compatible": len(s13) + len(s23) - int(rank36),
        "compatibility_conditions": len(left),
        "rhs_divisible_by_G2": bool(divisible),
        "nonzero_residuals": residuals,
        "G_divisible_branch": "P14 = G*u, deg u = 9, xdeg u <= 7",
        "G_divisible_branch_dim": len(support_u) if branch_ok else None,
        "branch_support_ok": bool(branch_ok),
    }


# --------------------------------------------------------------------------
# stage 2b -- the classical y-adic approximate root (cross-check, (M1))
# --------------------------------------------------------------------------

def yadic_approximate_root(poly, m, d=5):
    """Polynomial part of poly^(1/d) in Q[x]((1/Y)), degree m = deg_Y/d.

    Returns (A, residual_degree) with deg_Y(poly - A^d) < m*d - m.
    The recursion is triangular and each new coefficient enters with
    multiplicity d, so the root always exists: it is a definition, not a
    constraint.
    """
    cs = sp.symbols("k0:%d" % m)
    A = Y**m + sum(cs[k] * Y**k for k in range(m))
    Ad = sp.Poly(sp.expand(A**d), Y)
    target = sp.Poly(sp.expand(poly), Y)
    subs = {}
    for j in range(m * d - 1, m * d - m - 1, -1):
        unknown = cs[m * d - 1 - j]
        coeff = sp.expand(Ad.coeff_monomial(Y**j).subs(subs))
        rhs = target.coeff_monomial(Y**j)
        sol = sp.solve(sp.Eq(coeff, rhs), unknown)
        subs[unknown] = sp.expand(sol[0])
    A = sp.expand(A.subs(subs))
    resid = sp.expand(poly - A**d)
    return A, sp.Poly(resid, X, Y).degree(Y)


def stage2b_yadic(sp_pts, sq_pts):
    """Build one random corner pair and compute both y-adic roots."""
    import random
    rng = random.Random(20260101)

    def random_poly(points, top_degree):
        terms = {}
        for (i, j) in points:
            terms[(i, j)] = rng.randint(-3, 3)
        terms[(0, top_degree)] = 1          # monic in Y
        return sum(c * X**i * Y**j for (i, j), c in terms.items())

    P = random_poly(sp_pts, 15)
    Q = random_poly(sq_pts, 25)
    AP, dP = yadic_approximate_root(P, 3, 5)
    AQ, dQ = yadic_approximate_root(Q, 5, 5)
    payload = {
        "P_y_degree": 15,
        "Q_y_degree": 25,
        "root_P_y_degree": 3,
        "root_Q_y_degree": 5,
        "same_root": False,
        "P_residual_y_degree": int(dP),
        "Q_residual_y_degree": int(dQ),
        "P_bound_ok": dP < 12,
        "Q_bound_ok": dQ < 20,
        "note": ("the two y-adic roots have different y-degrees and are "
                 "different polynomials: the classical construction has no "
                 "common root at this corner"),
    }
    log_stage("2b_yadic", payload)
    return payload


# --------------------------------------------------------------------------
# stage 4 -- bookkeeping that closes the report
# --------------------------------------------------------------------------

def stage4_chain(sp_pts, sq_pts):
    """Per-degree new unknowns vs monomials of the det J polygon (counts only)."""
    np_ = hull([(3 * a, 3 * b) for a, b in DELTA])
    nq_ = hull([(5 * a, 5 * b) for a, b in DELTA])
    mink = hull([(a + c, b + d) for a, b in np_ for c, d in nq_])
    det_pts = lattice_points(mink)
    from collections import Counter
    cp = Counter(sum(p) for p in sp_pts)
    cq = Counter(sum(p) for p in sq_pts)
    cd = Counter(sum(p) - 2 for p in det_pts)  # shift by -(1,1) on derivatives
    rows = []
    for d in range(38, 12, -1):
        new_p = cp.get(d - 23, 0)
        new_q = cq.get(d - 13, 0)
        eqs = cd.get(d, 0)
        rows.append({"degree": d, "new_P": new_p, "new_Q": new_q,
                     "new": new_p + new_q, "monomials": eqs})
    homogeneous = [r["degree"] for r in rows if r["degree"] >= 3]
    payload = {
        "detJ_lattice_points": len(det_pts),
        "levels_solved": [38, 37, 36],
        "levels_homogeneous": [min(homogeneous), max(homogeneous)],
        "first_inhomogeneous_level": 2,
        "chain_rows": rows[:6],
    }
    log_stage("3_leading_order", payload)
    return payload


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def main() -> int:
    print("=" * 74)
    print("Approximate roots at the (75,125) corner F_2(3,5), j=1")
    print("=" * 74)

    # ---- stage 1
    src, extracted = stage1_sources()
    print("\n[stage 1] sources")
    print("  Moh 1983 local text :", extracted["moh_txt_present"])
    print("  Abhyankar TIFR 57   :", extracted["abhyankar_pdf_present"])
    print("  definition extracted:", extracted["moh_definition_line"] is not None)
    print("  variable of the classical root: y (monic in y, d | deg_y f)")

    # ---- stage 2
    c2, np_, nq_, sp_pts, sq_pts, G = stage2_corner()
    print("\n[stage 2] corner and the common approximate root")
    print("  N(P) = 3Delta', N(Q) = 5Delta'")
    print("  lattice points            : P %d, Q %d, total %d"
          % (c2["lattice_P"], c2["lattice_Q"], c2["corner_unknowns"]))
    print("  degrees                   : P %d, Q %d" % tuple(c2["degree"]))
    print("  top edge (total degree)   : P %d, Q %d monomials"
          % (c2["top_edge_P"], c2["top_edge_Q"]))
    print("  approximate root G        : degree 5, support [(0,5),(1,4),(2,3),(3,2),(4,1)]")
    print("  P = G^3 + R,  Q = G^5 + S : %d + %d corrections"
          % (c2["correction_P"], c2["correction_Q"]))
    print("  J(G^3,G^5) = 0            :", c2["leading_jacobian_identically_zero"])
    print("  top block %d -> 5 (4 shape parameters of G + modulus d)"
          % c2["top_block_unknowns"])
    print("  unknowns after substitution: %d  (from %d)"
          % (c2["substituted_unknowns"], c2["corner_unknowns"]))
    yad = stage2b_yadic(sp_pts, sq_pts)
    print("  y-adic cross-check        : root of P has y-degree %d (residual y-deg %d),"
          " root of Q has y-degree %d (residual y-deg %d)"
          % (yad["root_P_y_degree"], yad["P_residual_y_degree"],
             yad["root_Q_y_degree"], yad["Q_residual_y_degree"]))
    print("                              different polynomials -> no common y-adic root")

    # ---- stage 3
    print("\n[stage 3] det J = x^2, first non-trivial order")
    d_ranks = {}
    for dval in (1, 2, 3, 5):
        info, _unused = degree_37_system(G, sp_pts, sq_pts, dval)
        d_ranks[dval] = (info["rank"], info["dimension"])
        if dval == 1:
            d37 = info
            P14, Q24, s14, s24 = _unused
    print("  degree 38 : identically zero (theorem P2)")
    print("  degree 37 : %d unknowns, %d equations, rank %d, dimension %d"
          % (d37["unknowns"], d37["equations"], d37["rank"], d37["dimension"]))
    print("              rank/dimension for d = 1,2,3,5 : %s"
          % ", ".join("%d:%s" % (k, v) for k, v in d_ranks.items()))
    print("              closed form Q_24 = (5d/3) G^2 P_14 :", d37["closed_form_ok"])
    print("              J(P_14,Q_24) = (10d/3) G P_14 J(P_14,G) :",
          d37["jacobian_identity_ok"])
    d36 = degree_36_compatibility(G, sp_pts, sq_pts, P14, Q24, s14, s24)
    print("  degree 36 : %d unknowns, %d equations, rank %d -> kernel %d"
          % (d36["unknowns"], d36["equations"], d36["rank36"],
             d36["solution_dimension_if_compatible"]))
    print("              degree-26 compatibility conditions : %d"
          % d36["compatibility_conditions"])
    print("              sampled level-37 point divisible by G^2 :",
          d36["rhs_divisible_by_G2"])
    print("              first non-zero branch              : %s (dim %s)"
          % (d36["G_divisible_branch"], d36["G_divisible_branch_dim"]))
    chain = stage4_chain(sp_pts, sq_pts)
    print("  levels >= 3 are homogeneous: no inconsistency above degree 2")
    print("  first inhomogeneous level  : degree 2 (the target x^2)")

    # ---- stage 4
    verdict = "REDUCED"
    reduced_unknowns = d37["unknowns"]
    print("\n[stage 4] verdict")
    print("  %s: top edge %d -> 5; first non-trivial system %d unknowns / %d "
          "equations / dimension %d"
          % (verdict, c2["top_block_unknowns"], d37["unknowns"],
             d37["equations"], d37["dimension"]))
    print("  scope: the remaining %d correction coefficients are NOT eliminated;"
          % (c2["substituted_unknowns"] - 5))
    print("         the reduction is the top edge plus the degree-37 triangle.")
    print("  not claimed: a degree bound, a reconstruction, or 445 -> 10.")

    log_stage("4_verdict", {
        "verdict": verdict,
        "reduced_unknowns": reduced_unknowns,
        "top_block_before": c2["top_block_unknowns"],
        "top_block_after": 5,
        "substituted_unknowns": c2["substituted_unknowns"],
        "degree_37": d37,
        "degree_37_d_scan": d_ranks,
        "degree_36": {k: v for k, v in d36.items()},
        "chain": chain,
    })

    OUT_PATH.write_text(json.dumps({
        "corner": "F_2(3,5), j=1, (75,125)",
        "delta_prime": [list(p) for p in DELTA],
        "reduced_unknowns": reduced_unknowns,
        "verdict": verdict,
        "stage1": src,
        "stage2": c2,
        "stage2b_yadic": yad,
        "degree_37": d37,
        "degree_36": d36,
        "chain": chain,
        "elapsed_s": round(elapsed(), 1),
    }, indent=2, default=str), encoding="utf-8")

    print("\nreduced_unknowns =", reduced_unknowns)
    print("verdict          =", verdict)
    print("elapsed          = %.1f s" % elapsed())
    return 0


if __name__ == "__main__":
    sys.exit(main())
