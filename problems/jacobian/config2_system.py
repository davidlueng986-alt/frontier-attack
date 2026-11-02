"""Size the GGHV config-2 coefficient system for the (72,108) Jacobian pair.

Source: J. A. Guccione, J. J. Guccione, R. Horruitiner, C. Valqui,
"Increasing the degree of a possible counterexample to the Jacobian Conjecture
from 100 to 108", arXiv:2204.14178, Proposition 4.3 (Case (8,28)), case (2):

    N(P) = {(0,0), (1,0), (8,14), (8,16)}
    N(Q) = {(0,0), (2,1), (12,21), (12,24)}

SUPPORT CONVENTION (the decision this script makes)
---------------------------------------------------
The displayed sets are **corner lists**, not supports.  Two independent reasons:

1. Prop. 4.3's proof says "the corners of the polygons of P and Q in this case
   are given by ...", and the proposition is stating a normal form for the
   *Newton polygons*, not a monomial ansatz.
2. GGHV use the standard definition elsewhere (arXiv:1401.1784, Sec. 1):
   "we let H(P) denote the convex hull of the support of P.  As it is well
   known, H(P) is a polygon, called the Newton polygon of P."  The Newton
   polygon is the hull; the support is a subset of it.

A *generic* polynomial with Newton polygon N therefore has every lattice point
of N in its support.  That is the primary convention here:

    |N(P) ∩ Z^2| = 25,   |N(Q) ∩ Z^2| = 47,   72 unknown coefficients.

The vertices-only reading (4 + 4 = 8 unknowns) is kept as a cross-check under
`variants.vertices_only`.  It is the minimal-support specialization -- the
reading the plan's "< 30 unknowns" estimate assumed -- but it is not generic:
its solutions are a strict subset of the full-lattice solutions, and it imposes
extra monomial equations (e.g. -20 a2 b1 = 0) that the generic ansatz does not.

NORMALIZATION (a correction to the task text)
---------------------------------------------
Prop. 4.3 normalizes [P,Q] = x^2, not [P,Q] = 1.  This is not cosmetic.
Neither polygon contains (0,1) as a point of N(P), nor (1,0)/(0,1) as points of
N(Q); hence P has no linear y-term and Q has no linear part at all, so

    det J(P,Q)|_(0,0) = P_x(0,0) Q_y(0,0) - P_y(0,0) Q_x(0,0) = 0

identically on these supports.  Consequently det J = 1 is impossible on either
support convention: the requested system contains the coefficient equation
"-1 = 0".  The script still builds det J - 1 exactly as requested, and also
builds det J - x^2, which is the system the paper actually needs.

DEGREES
-------
These are the *reduced* supports: deg P = 24 and deg Q = 36 (total degree), not
72/108.  The (72,108) pair of the theorem is obtained from this reduced pair by
further transformations; no degree-108 polynomial is expanded here.  Expansion
is cheap (measured: ~2 s end to end, 92 monomials in det J), so no truncation is
needed.  A Minkowski-sum support bound is used automatically only if a future
support is large enough to make the dense expansion expensive.

Run:
    python problems/jacobian/config2_system.py --report
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / "runs" / "jacobian_config2_size.json"

ARXIV_ID = "2204.14178"
ARXIV_URL = f"https://arxiv.org/abs/{ARXIV_ID}"
AR5IV_URL = f"https://ar5iv.labs.arxiv.org/html/{ARXIV_ID}"

P_CORNERS = ((0, 0), (1, 0), (8, 14), (8, 16))
Q_CORNERS = ((0, 0), (2, 1), (12, 21), (12, 24))

# Dense expansion is skipped above this many derivative-pair products.
MAX_PAIR_PRODUCTS = 2_000_000


# --------------------------------------------------------------------------
# lattice geometry
# --------------------------------------------------------------------------

def convex_hull(points):
    """Andrew monotone chain; returns vertices counter-clockwise, no collinear."""
    pts = sorted(set(points))

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


def lattice_points(vertices):
    """All lattice points of the closed convex polygon, boundary included."""
    hull = convex_hull(vertices)
    xs = [p[0] for p in hull]
    ys = [p[1] for p in hull]
    out = []
    for i in range(min(xs), max(xs) + 1):
        for j in range(min(ys), max(ys) + 1):
            ok = True
            for k in range(len(hull)):
                a = hull[k]
                b = hull[(k + 1) % len(hull)]
                if (b[0] - a[0]) * (j - a[1]) - (b[1] - a[1]) * (i - a[0]) < 0:
                    ok = False
                    break
            if ok:
                out.append((i, j))
    return out


def pick_check(vertices):
    """Independent lattice-point count via Pick's theorem."""
    hull = convex_hull(vertices)
    twice_area = abs(sum(
        hull[i][0] * hull[(i + 1) % len(hull)][1]
        - hull[(i + 1) % len(hull)][0] * hull[i][1]
        for i in range(len(hull))
    ))
    boundary = sum(
        sp.igcd(abs(hull[(i + 1) % len(hull)][0] - hull[i][0]),
                abs(hull[(i + 1) % len(hull)][1] - hull[i][1]))
        for i in range(len(hull))
    )
    interior = (twice_area - boundary + 2) // 2
    return {
        "twice_area": int(twice_area),
        "boundary_points": int(boundary),
        "interior_points": int(interior),
        "total_points": int(interior + boundary),
    }


def derivative_support(support, axis):
    out = set()
    for a, b in support:
        if axis == 0 and a > 0:
            out.add((a - 1, b))
        if axis == 1 and b > 0:
            out.add((a, b - 1))
    return out


def minkowski_sum(left, right):
    return {(a + c, b + d) for a, b in left for c, d in right}


def det_support_bound(support_p, support_q):
    """Exact monomial superset of det J (before cancellation).

    Every monomial of P_x Q_y - P_y Q_x lies in one of the two Minkowski sums;
    cancellations can only remove monomials, so len(bound) >= number of terms.
    """
    px = derivative_support(support_p, 0)
    py = derivative_support(support_p, 1)
    qx = derivative_support(support_q, 0)
    qy = derivative_support(support_q, 1)
    return minkowski_sum(px, qy) | minkowski_sum(py, qx)


# --------------------------------------------------------------------------
# polynomial system
# --------------------------------------------------------------------------

def build_polynomial(support, prefix, x, y):
    coeffs = [sp.Symbol(f"{prefix}{i}") for i in range(len(support))]
    expr = sp.Add(*[c * x ** a * y ** b for c, (a, b) in zip(coeffs, support)])
    return expr, coeffs


def _degree_in(c, unknown_set):
    """Total degree of c in the coefficient symbols.

    det J is bilinear, so every term is a product of at most two symbols with
    exponent one; counting powers avoids building a Poly over 72 symbols for
    every one of the ~90 equations (which cost ~50 s per system).
    """
    best = 0
    for term in sp.Add.make_args(c):
        powers = term.as_powers_dict()
        best = max(best, sum(int(powers.get(u, 0)) for u in unknown_set))
    return best


def size_system(support_p, support_q, target, label, x, y, cross_check=False):
    """Build P, Q on the supports, expand det J - target, count equations."""
    t0 = time.perf_counter()
    unknown_count = len(support_p) + len(support_q)
    pair_products = (
        len(derivative_support(support_p, 0)) * len(derivative_support(support_q, 1))
        + len(derivative_support(support_p, 1)) * len(derivative_support(support_q, 0))
    )

    P, cp = build_polynomial(support_p, "a", x, y)
    Q, cq = build_polynomial(support_q, "b", x, y)
    unknowns = cp + cq
    unknown_set = set(unknowns)
    t_build = time.perf_counter() - t0

    record = {
        "label": label,
        "support_sizes": [len(support_p), len(support_q)],
        "unknowns": unknown_count,
        "pair_products": pair_products,
        "degs_PQ": [max(a + b for a, b in support_p),
                    max(a + b for a, b in support_q)],
        "target": str(target),
        "mode": "dense",
    }

    if pair_products > MAX_PAIR_PRODUCTS:
        # Honest fallback: no expansion, only exact support bound.
        bound = det_support_bound(support_p, support_q)
        record.update({
            "mode": "support_bound_only",
            "reason": f"{pair_products} derivative pairs exceeds guard "
                      f"{MAX_PAIR_PRODUCTS}",
            "det_monomials": len(bound),
            "equations": len(bound),
            "equations_are_upper_bound": True,
            "max_total_degree": max(a + b for a, b in bound),
            "max_total_degree_in_coefficients": 2,
            "seconds": round(time.perf_counter() - t0, 3),
            "determinedness": ("overdetermined" if len(bound) > unknown_count
                               else "underdetermined"),
        })
        return record

    t1 = time.perf_counter()
    det = sp.expand(P.diff(x) * Q.diff(y) - P.diff(y) * Q.diff(x))
    t_expand = time.perf_counter() - t1

    if cross_check:
        det_matrix = sp.expand(sp.Matrix([P, Q]).jacobian([x, y]).det())
        if sp.expand(det_matrix - det) != 0:
            raise RuntimeError("manual Jacobian != sympy Matrix.jacobian determinant")

    t2 = time.perf_counter()
    det_poly = sp.Poly(det, x, y)
    det_monomials = len(det_poly.monoms())
    det_has_constant = (0, 0) in det_poly.monoms()
    poly = sp.Poly(det - target, x, y)
    monoms = [(int(m[0]), int(m[1])) for m in poly.monoms()]
    equations = []
    max_coeff_degree = 0
    max_terms = 0
    reconstruction = sp.Integer(0)
    for index, m in enumerate(monoms):
        c = poly.coeff_monomial(m)
        reconstruction += c * x ** m[0] * y ** m[1]
        n_terms = len(sp.Add.make_args(c))
        max_terms = max(max_terms, n_terms)
        deg_c = _degree_in(c, unknown_set)
        max_coeff_degree = max(max_coeff_degree, deg_c)
        equations.append({
            "monomial": [m[0], m[1]],
            "total_degree": m[0] + m[1],
            "n_terms": n_terms,
            "degree_in_coefficients": int(deg_c),
            # factor only a few: 72-symbol factorisation of all ~90
            # coefficients is wasted work when the JSON keeps samples.
            "coefficient": (str(sp.factor(c))[:400]
                            if index < 6 or m == (0, 0) else None),
        })
    t_extract = time.perf_counter() - t2

    reconstruction_ok = sp.expand(reconstruction - (det - target)) == 0
    if not reconstruction_ok:
        raise RuntimeError(
            "coefficient extraction failed: reconstructed polynomial != det - target")

    constant = next((e for e in equations if e["monomial"] == [0, 0]), None)
    record.update({
        "det_monomials": det_monomials,
        "equations": len(equations),
        "equations_are_upper_bound": False,
        "max_total_degree": max(m[0] + m[1] for m in monoms) if monoms else 0,
        "min_total_degree": min(m[0] + m[1] for m in monoms) if monoms else 0,
        "max_total_degree_in_coefficients": int(max_coeff_degree),
        "max_terms_per_equation": max_terms,
        "constant_equation": constant,
        "reconstruction_ok": bool(reconstruction_ok),
        "det_has_constant_term": bool(det_has_constant),
        "determinedness": ("overdetermined" if len(equations) > unknown_count
                           else "underdetermined" if len(equations) < unknown_count
                           else "square"),
        "seconds_build": round(t_build, 3),
        "seconds_expand": round(t_expand, 3),
        "seconds_extract": round(t_extract, 3),
        "seconds": round(time.perf_counter() - t0, 3),
        "equations_sample": equations[:6],
        "equation_monomials": [e["monomial"] for e in equations],
    })
    if cross_check:
        record["jacobian_matrix_ok"] = True
    return record


def source_evidence(check_live=False):
    evidence = {
        "arxiv": ARXIV_ID,
        "url": ARXIV_URL,
        "html": AR5IV_URL,
        "proposition": "4.3 (Case (8,28)), case (2)",
        "verified": "2026-11-02",
        "quote": ("(2) N(P)={(0,0),(1,0),(8,14),(8,16)}, "
                  "N(Q)={(0,0),(2,1),(12,21),(12,24)}"),
        "corner_language": ("Proof: 'the corners of the polygons of P and Q in "
                            "this case are given by {(0,0),(1,0),(8,28),(0,4)}, "
                            "multiplied by (m,n)=(3,2)'"),
        "newton_polygon_definition": (
            "arXiv:1401.1784, Sec. 1: H(P) is 'the convex hull of the support "
            "of P', 'called the Newton polygon of P'"),
        "normalization": ("Prop. 4.3: 'there exist P,Q in L^(1) with "
                          "[P,Q] = x^2'"),
    }
    if check_live:
        import urllib.request

        info = {"attempted": True, "reachable": False}
        try:
            req = urllib.request.Request(AR5IV_URL, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=40) as resp:
                body = resp.read().decode("utf-8", "replace")
            flat = " ".join(body.split())
            info["reachable"] = True
            info["status"] = 200
            info["bytes"] = len(body)
            info["prop_4_3_present"] = "Proposition 4.3" in flat
            info["quote_present"] = ("(0,0),(1,0),(8,14),(8,16)" in flat
                                     and "(0,0),(2,1),(12,21),(12,24)" in flat)
        except Exception as exc:  # network failure must not fail --report
            info["error"] = f"{type(exc).__name__}: {exc}"
        evidence["live_check"] = info
    return evidence


# --------------------------------------------------------------------------
# report
# --------------------------------------------------------------------------

def build_report(check_live=True):
    x, y = sp.symbols("x y")
    full_p = lattice_points(P_CORNERS)
    full_q = lattice_points(Q_CORNERS)

    full = size_system(full_p, full_q, sp.Integer(1), "full_lattice/detJ-1", x, y,
                       cross_check=True)
    full_x2 = size_system(full_p, full_q, x ** 2, "full_lattice/detJ-x^2", x, y)
    vert = size_system(list(P_CORNERS), list(Q_CORNERS), sp.Integer(1),
                       "vertices_only/detJ-1", x, y)
    vert_x2 = size_system(list(P_CORNERS), list(Q_CORNERS), x ** 2,
                          "vertices_only/detJ-x^2", x, y)

    infeasible = (not full["det_has_constant_term"]
                  and full["constant_equation"] is not None)
    effective_unknowns = full["unknowns"] - 2  # two additive constants

    return {
        "problem": "jacobian-config2",
        "arxiv": ARXIV_ID,
        "proposition": "4.3, case (2)",
        "support_convention": {
            "choice": "full lattice points of the Newton polygon (generic)",
            "why": (
                "Prop. 4.3 lists corners ('the corners of the polygons ... are "
                "given by'), and GGHV define the Newton polygon as the convex "
                "hull of the support (arXiv:1401.1784 Sec. 1). A generic "
                "polynomial with that Newton polygon has every lattice point of "
                "the polygon in its support. Vertices-only is a minimal-support "
                "specialization, reported as a variant, not the generic model."
            ),
            "vertices_are_supports": False,
        },
        "supports": {
            "N(P)_corners": [list(p) for p in P_CORNERS],
            "N(Q)_corners": [list(p) for p in Q_CORNERS],
            "N(P)_hull": [list(p) for p in convex_hull(P_CORNERS)],
            "N(Q)_hull": [list(p) for p in convex_hull(Q_CORNERS)],
            "lattice_points_P": len(full_p),
            "lattice_points_Q": len(full_q),
            "lattice_points_total": len(full_p) + len(full_q),
        },
        "unknowns": full["unknowns"],
        "equations": full["equations"],
        "max_total_degree": full["max_total_degree"],
        "max_total_degree_in_coefficients": full["max_total_degree_in_coefficients"],
        "determinedness": full["determinedness"],
        "effective": {
            "unknowns": effective_unknowns,
            "note": "the additive constants of P and Q never enter det J",
            "nontrivial_equations": full["equations"] - 1,
            "equations_minus_effective_unknowns": (full["equations"] - 1) - effective_unknowns,
        },
        "feasibility": {
            "det_J_equals_1": "infeasible",
            "reason": (
                "det J has no constant term on these supports: P has no linear "
                "y-term and Q has no linear part, so the (0,0) coefficient "
                "equation is -1 = 0."
            ),
            "det_J_equals_x2": "not solved (per paper normalization [P,Q]=x^2)",
            "infeasible_by_constant_term": bool(infeasible),
        },
        "requested_system": {
            "target": "det J(P,Q) - 1",
            "unknowns": full["unknowns"],
            "equations": full["equations"],
            "max_total_degree": full["max_total_degree"],
            "constant_equation": full["constant_equation"],
            "seconds": full["seconds"],
        },
        "normalization_corrected": {
            "target": "det J(P,Q) - x^2",
            "source": "Prop. 4.3 states [P,Q] = x^2",
            "unknowns": full_x2["unknowns"],
            "equations": full_x2["equations"],
            "max_total_degree": full_x2["max_total_degree"],
            "max_total_degree_in_coefficients": full_x2[
                "max_total_degree_in_coefficients"],
            "seconds": full_x2["seconds"],
            "sample_equations": full_x2["equations_sample"],
        },
        "variants": {
            "full_lattice": {
                "detJ_minus_1": _strip(full),
                "detJ_minus_x2": _strip(full_x2),
            },
            "vertices_only": {
                "detJ_minus_1": _strip(vert),
                "detJ_minus_x2": _strip(vert_x2),
            },
        },
        "timings": {
            "full_lattice_detJ_minus_1_seconds": full["seconds"],
            "full_lattice_detJ_minus_x2_seconds": full_x2["seconds"],
            "vertices_only_seconds": round(vert["seconds"] + vert_x2["seconds"], 3),
            "expanded_dense": True,
            "blow_up": False,
            "note": (
                "The reduced supports have deg P = 24, deg Q = 36, not 72/108; "
                "the degree-108 pair is obtained from this reduced pair by "
                "further transformations. Dense expansion stayed small, so no "
                "truncated bound was needed."
            ),
        },
        "degrees": {
            "deg_P_total": max(a + b for a, b in P_CORNERS),
            "deg_Q_total": max(a + b for a, b in Q_CORNERS),
            "detJ_max_total_degree": full["max_total_degree"],
            "detJ_min_total_degree": full["min_total_degree"],
        },
        "source_evidence": source_evidence(check_live=check_live),
        "verification": {
            "pick_P": pick_check(P_CORNERS),
            "pick_Q": pick_check(Q_CORNERS),
            "lattice_enumeration_matches_pick": (
                pick_check(P_CORNERS)["total_points"] == len(full_p)
                and pick_check(Q_CORNERS)["total_points"] == len(full_q)
            ),
            "jacobian_matrix_matches_manual": bool(
                full.get("jacobian_matrix_ok", False)),
            "reconstruction_ok_all": all(
                r.get("reconstruction_ok")
                for r in (full, full_x2, vert, vert_x2)),
        },
        "notes": [
            "No solving or elimination was performed.",
            "Every coefficient polynomial was reconstructed into det J - target "
            "and the identity checked exactly (reconstruction_ok).",
            "Constants of P and Q are free (translation invariance of det J).",
            "The plan's '< 30 unknowns' holds only for the vertices-only "
            "variant (8 unknowns); the generic full-lattice system has 72.",
            "Vertices-only is a specialization: its solutions form a subset of "
            "the full-lattice solutions, but a generic polynomial with this "
            "Newton polygon can use every lattice point of the polygon.",
        ],
        "elapsed_seconds": None,
    }


def _strip(record):
    return {k: v for k, v in record.items()
            if k not in ("equations_sample", "equation_monomials")}


def print_report(rep):
    print("== GGHV 2022 arXiv:%s, Prop. 4.3, case (2) ==" % ARXIV_ID)
    print("support convention : %s" % rep["support_convention"]["choice"])
    print("  N(P): %d corners, %d lattice points"
          % (len(P_CORNERS), rep["supports"]["lattice_points_P"]))
    print("  N(Q): %d corners, %d lattice points"
          % (len(Q_CORNERS), rep["supports"]["lattice_points_Q"]))
    print("  deg P = %d, deg Q = %d (total degree)"
          % (rep["degrees"]["deg_P_total"], rep["degrees"]["deg_Q_total"]))
    print()
    print("requested system   : det J(P,Q) - 1")
    print("unknowns           : %d" % rep["requested_system"]["unknowns"])
    print("equations          : %d" % rep["requested_system"]["equations"])
    print("max total degree   : %d (in x,y); %d (in coefficients)"
          % (rep["requested_system"]["max_total_degree"],
             rep["max_total_degree_in_coefficients"]))
    print("determinedness     : %s (%d equations / %d unknowns; "
          "%d effective unknowns)"
          % (rep["determinedness"], rep["equations"], rep["unknowns"],
             rep["effective"]["unknowns"]))
    print("timing             : %.3f s (dense build+expand+extract, no blow-up)"
          % rep["requested_system"]["seconds"])
    print()
    print("WARNING: %s" % rep["feasibility"]["reason"])
    print("         det J = 1 is infeasible on this support; the paper's")
    print("         normalization is [P,Q] = x^2, which gives %d equations."
          % rep["normalization_corrected"]["equations"])
    print()
    v = rep["variants"]["vertices_only"]["detJ_minus_1"]
    print("cross-check (vertices-only reading): %d unknowns, %d equations, "
          "max degree %d, %.3f s"
          % (v["unknowns"], v["equations"], v["max_total_degree"], v["seconds"]))
    print("  -> the plan's '< 30 unknowns' is true only for this specialization.")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--report", action="store_true",
                    help="print the size report (default action)")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT,
                    help="JSON output path")
    ap.add_argument("--no-source-check", action="store_true",
                    help="skip the arXiv:2204.14178 reachability check")
    ap.add_argument("--no-write", action="store_true", help="do not write JSON")
    args = ap.parse_args(argv)

    t0 = time.perf_counter()
    rep = build_report(check_live=not args.no_source_check)
    rep["elapsed_seconds"] = round(time.perf_counter() - t0, 3)

    print_report(rep)

    if not args.no_write:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(rep, indent=2), encoding="utf-8")
        print("\nwritten: %s" % args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
