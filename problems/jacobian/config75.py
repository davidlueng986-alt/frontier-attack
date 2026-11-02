"""Size the coefficient system at the (75,125) Jacobian corner, F_2(3,5), j=1.

CORNER DATA
-----------
Verbatim from the wstrinz repo, cited in
``research/jacobian-degree-pairs-ge125.md``:

    Delta' = {(0,0),(3,0),(4,1),(0,5)}
    N(P) = 3*Delta' = {(0,0),(9,0),(12,3),(0,15)}
    N(Q) = 5*Delta' = {(0,0),(15,0),(20,5),(0,25)}
    [P,Q] = x^2      reduced degrees (15,25)

Source: ``research/scratch/wstrinz/d2/PASSPORT_75_125_REPAIR.md:29-33``.
The case dossier is ``research/scratch/wstrinz/d2/case_dossier_F2_j1_75_125.json``
(``case``: A0=[5,20], A0'=[1,0], A1="(7\\5,2)", a=3, b=5, m=3, n=5, k=1).

SUPPORT CONVENTION
------------------
Two readings are sized:

* ``full_lattice``: every lattice point of the closed Newton polygon.  This is
  the generic polynomial with that Newton polygon, and it is the primary
  convention here (the same convention as ``config2_system.py``).
* ``vertices_only``: only the four corners of each polygon.  This is a
  minimal-support specialization, not the generic support.

The polygons are closed, boundary included.  Pick's theorem gives
124 points for N(P) (85 interior + 39 boundary) and 321 for N(Q)
(256 interior + 65 boundary), hence 445 unknowns in the full-lattice variant.

METHOD
------
``det J`` is bilinear in the coefficients.  Instead of constructing a SymPy
polynomial in 445 symbols and expanding it, the script accumulates each
coefficient equation as an integer combination of products ``a_i*b_j``.
That is exact (no truncation), much smaller than the dense symbolic object,
and gives the monomial support, the coefficient degrees, and the constant
term directly.  The vertices-only result is cross-checked against
``config2_system.size_system`` (SymPy) in the report.  A full SymPy
cross-check is available behind ``--slow-cross-check`` (about 95-105 s measured).

The target is ``det J - x^2 = 0``.  The script also sizes ``det J - 1 = 0``
for comparison; at this corner the two systems have the same number of
equations, but the constant equations differ.

Run:
    python problems/jacobian/config75.py --report
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import sympy as sp

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from config2_system import (  # noqa: E402
    MAX_PAIR_PRODUCTS,
    convex_hull,
    derivative_support,
    det_support_bound,
    lattice_points,
    pick_check,
    size_system,
)

ROOT = HERE.parents[1]
DEFAULT_OUT = ROOT / "runs" / "jacobian_75_125_size.json"

DELTA_PRIME = ((0, 0), (3, 0), (4, 1), (0, 5))
P_CORNERS = ((0, 0), (9, 0), (12, 3), (0, 15))
Q_CORNERS = ((0, 0), (15, 0), (20, 5), (0, 25))

TARGET_X2 = (2, 0)
TARGET_ONE = (0, 0)


# --------------------------------------------------------------------------
# exact bilinear determinant
# --------------------------------------------------------------------------

def _points_hash(points):
    payload = json.dumps(sorted(points), separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _indexed_support(support, prefix):
    """``[(a, b, prefix, i)]`` in the same order as ``build_polynomial``."""
    return [(int(a), int(b), prefix, i) for i, (a, b) in enumerate(support)]


def _derivative_terms(indexed, axis):
    """Terms of d/dx (axis=0) or d/dy (axis=1), with integer exponents."""
    out = []
    for a, b, prefix, i in indexed:
        if axis == 0 and a > 0:
            out.append((a - 1, b, a, prefix, i))
        elif axis == 1 and b > 0:
            out.append((a, b - 1, b, prefix, i))
    return out


def _pair_products(support_p, support_q):
    """Number of derivative-pair products in the dense expansion."""
    return (
        len(derivative_support(support_p, 0))
        * len(derivative_support(support_q, 1))
        + len(derivative_support(support_p, 1))
        * len(derivative_support(support_q, 0))
    )


def det_coefficients(support_p, support_q):
    """Exact coefficient map of det J.

    Returns ``(det, pair_products)``, where ``det[(p, q)]`` is a dict
    ``{(i, j): integer}``.  The coefficient of ``x^p y^q`` is then

        sum_{i,j} integer * a_i * b_j .

    ``det J`` is bilinear, so every monomial has coefficient degree exactly
    two (or is absent); no higher products can occur.  Zero coefficients are
    dropped, so the monomial count is exact, not a Minkowski upper bound.
    """
    p_terms = _indexed_support(support_p, "a")
    q_terms = _indexed_support(support_q, "b")
    px = _derivative_terms(p_terms, 0)
    py = _derivative_terms(p_terms, 1)
    qx = _derivative_terms(q_terms, 0)
    qy = _derivative_terms(q_terms, 1)

    pair_products = len(px) * len(qy) + len(py) * len(qx)
    det = defaultdict(lambda: defaultdict(int))
    for ax, ay, ca, _, ia in px:
        for bx, by, cb, _, ib in qy:
            det[(ax + bx, ay + by)][(ia, ib)] += ca * cb
    for ax, ay, ca, _, ia in py:
        for bx, by, cb, _, ib in qx:
            det[(ax + bx, ay + by)][(ia, ib)] -= ca * cb

    out = {}
    for monomial, terms in det.items():
        cleaned = {pair: coeff for pair, coeff in terms.items() if coeff}
        if cleaned:
            out[monomial] = cleaned
    return out, pair_products


def _format_coeff(eq):
    """Readable coefficient polynomial from the ``{(i, j): integer}`` map."""
    if not eq:
        return "0"
    parts = []
    for (i, j), coeff in sorted(eq.items()):
        mono = f"a{i}*b{j}"
        if coeff == 1:
            parts.append(mono)
        elif coeff == -1:
            parts.append("-" + mono)
        else:
            parts.append(f"{coeff}*{mono}")
    return " + ".join(parts).replace("+ -", "- ")


def _format_equation(eq, subtract_one):
    """Coefficient side of a target equation; the RHS is zero."""
    base = _format_coeff(eq)
    if not subtract_one:
        return base
    return "-1" if base == "0" else base + " - 1"


def _equation_record(monomial, eq, subtract_one):
    p, q = monomial
    n_terms = len(eq) + (1 if subtract_one else 0)
    return {
        "monomial": [int(p), int(q)],
        "total_degree": int(p + q),
        "n_terms": int(n_terms),
        "degree_in_coefficients": 2 if eq else 0,
        "coefficient": _format_equation(eq, subtract_one)[:500],
    }


# --------------------------------------------------------------------------
# sizing
# --------------------------------------------------------------------------

def _support_bound_det(support_p, support_q):
    bound = det_support_bound(support_p, support_q)
    return {
        "mode": "support_bound_only",
        "monomials": len(bound),
        "monomials_are_upper_bound": True,
        "support_bound_monomials": len(bound),
        "cancellations_vs_support_bound": None,
        "cancelled_monomials": [],
        "max_total_degree": max(p + q for p, q in bound),
        "min_total_degree": min(p + q for p, q in bound),
        "max_total_degree_in_coefficients": 2,
        "max_terms_per_equation": None,
        "constant_term_present": None,
        "constant_term_zero": None,
        "constant_coefficient": None,
        "pair_products": None,
        "seconds": None,
    }


def _bound_target_record(target, unknowns, bound, label):
    """Honest upper-bound record used when the dense expansion is skipped."""
    target_monomial = (int(target[0]), int(target[1]))
    monomials = set(bound) | {target_monomial}
    equations = len(monomials)
    effective = max(0, unknowns - 2)
    return {
        "label": label,
        "target": "x^2" if target_monomial == TARGET_X2 else "1",
        "target_support": [[target_monomial[0], target_monomial[1]]],
        "unknowns": unknowns,
        "equations": equations,
        "equations_are_upper_bound": True,
        "det_monomials": None,
        "max_total_degree": max(p + q for p, q in monomials),
        "min_total_degree": min(p + q for p, q in monomials),
        "max_total_degree_in_coefficients": 2,
        "max_terms_per_equation": None,
        "coefficient_degree_histogram": {},
        "constant_equation": None,
        "target_monomial_equation": None,
        "constant_term_of_det_zero": None,
        "target_monomial_present_in_det": None,
        "constant_term_test": "not evaluated (support-bound mode)",
        "constant_term_test_infeasible": None,
        "target_monomial_test": "not evaluated (support-bound mode)",
        "infeasible_by_target_monomial": None,
        "effective_unknowns": effective,
        "nontrivial_equations": None,
        "equations_minus_effective_unknowns": equations - effective,
        "determinedness": (
            "overdetermined" if equations > unknowns
            else "underdetermined" if equations < unknowns
            else "square"
        ),
        "determinedness_effective": (
            "overdetermined" if equations > effective
            else "underdetermined" if equations < effective
            else "square"
        ),
        "sample_equations": [],
        "equation_monomials": [
            [m[0], m[1]]
            for m in sorted(monomials, key=lambda z: (z[0] + z[1], z[0], z[1]))
        ],
        "equation_summaries": [],
    }


def _target_record(target, det, det_record, unknowns, label):
    """Size one target system, det J - target = 0."""
    target_monomial = (int(target[0]), int(target[1]))
    monomials = set(det) | {target_monomial}

    summaries = []
    degree_hist = Counter()
    max_terms = 0
    max_coeff_degree = 0
    nontrivial = 0
    for m in sorted(monomials, key=lambda z: (z[0] + z[1], z[0], z[1])):
        eq = det.get(m, {})
        subtract = (m == target_monomial)
        n_terms = len(eq) + (1 if subtract else 0)
        coeff_degree = 2 if eq else 0
        max_terms = max(max_terms, n_terms)
        max_coeff_degree = max(max_coeff_degree, coeff_degree)
        if n_terms:
            nontrivial += 1
        degree_hist[coeff_degree] += 1
        summaries.append({
            "monomial": [m[0], m[1]],
            "total_degree": m[0] + m[1],
            "n_terms": n_terms,
            "degree_in_coefficients": coeff_degree,
        })

    def sample_for(m):
        eq = det.get(m, {})
        return _equation_record(m, eq, m == target_monomial)

    # A small, useful sample: for small systems keep every equation; for the
    # 445-unknown system keep the target monomial, the constant equation, the
    # densest equation, and the highest-degree equations.
    if len(monomials) <= 12:
        sample_monomials = sorted(
            monomials, key=lambda z: (-(z[0] + z[1]), z[0], z[1]))
    else:
        sample_monomials = [target_monomial]
        if (0, 0) in monomials:
            sample_monomials.append((0, 0))
        densest = max(
            monomials,
            key=lambda m: (
                len(det.get(m, {})) + (1 if m == target_monomial else 0)),
        )
        sample_monomials.append(densest)
        sample_monomials.extend(sorted(
            monomials, key=lambda z: (-(z[0] + z[1]), z[0], z[1])
        )[:3])
    samples = []
    seen = set()
    for m in sample_monomials:
        if m in seen:
            continue
        seen.add(m)
        samples.append(sample_for(m))

    constant_eq = None
    if (0, 0) in monomials:
        constant_eq = _equation_record(
            (0, 0), det.get((0, 0), {}), target_monomial == (0, 0))

    target_eq = sample_for(target_monomial)
    target_present = target_monomial in det

    # The constant-term test.  For target x^2 the target has no constant
    # term, so a missing det constant term is 0 = 0, not a contradiction.
    # For target 1 the constant equation is -1 = 0 when det has none.
    det_constant_zero = not det_record["constant_term_present"]
    if target_monomial == TARGET_ONE:
        infeasible_by_constant_term = det_constant_zero
        constant_test = (
            "-1 = 0 (det J has no constant term, target is 1)"
            if det_constant_zero else
            "c00 - 1 = 0 with c00 nonzero (no contradiction)"
        )
    else:
        infeasible_by_constant_term = False
        constant_test = (
            "0 = 0 (det J and x^2 both have zero constant term; vacuous)"
            if det_constant_zero else
            "c00 = 0 with c00 nonzero (a genuine equation, no contradiction)"
        )

    infeasible_by_target_monomial = not target_present
    if infeasible_by_target_monomial:
        target_test = (
            f"x^{target_monomial[0]}*y^{target_monomial[1]}: absent from det J, "
            "so the equation is -1 = 0"
        )
    else:
        target_test = (
            f"x^{target_monomial[0]}*y^{target_monomial[1]}: present in det J "
            f"({len(det[target_monomial])} coefficient terms), so the equation "
            "is nontrivial"
        )

    effective_unknowns = max(0, unknowns - 2)
    equations = len(monomials)
    return {
        "label": label,
        "target": "x^2" if target_monomial == TARGET_X2 else "1",
        "target_support": [[target_monomial[0], target_monomial[1]]],
        "unknowns": unknowns,
        "equations": equations,
        "equations_are_upper_bound": False,
        "det_monomials": det_record["monomials"],
        "max_total_degree": max(m[0] + m[1] for m in monomials),
        "min_total_degree": min(m[0] + m[1] for m in monomials),
        "max_total_degree_in_coefficients": int(max_coeff_degree),
        "max_terms_per_equation": int(max_terms),
        "coefficient_degree_histogram": {
            str(k): v for k, v in sorted(degree_hist.items())
        },
        "constant_equation": constant_eq,
        "target_monomial_equation": target_eq,
        "constant_term_of_det_zero": bool(det_constant_zero),
        "target_monomial_present_in_det": bool(target_present),
        "constant_term_test": constant_test,
        "constant_term_test_infeasible": bool(infeasible_by_constant_term),
        "target_monomial_test": target_test,
        "infeasible_by_target_monomial": bool(infeasible_by_target_monomial),
        "effective_unknowns": effective_unknowns,
        "nontrivial_equations": int(nontrivial),
        "equations_minus_effective_unknowns": int(equations - effective_unknowns),
        "determinedness": (
            "overdetermined" if equations > unknowns
            else "underdetermined" if equations < unknowns
            else "square"
        ),
        "determinedness_effective": (
            "overdetermined" if equations > effective_unknowns
            else "underdetermined" if equations < effective_unknowns
            else "square"
        ),
        "sample_equations": samples,
        "equation_monomials": [
            [m[0], m[1]]
            for m in sorted(monomials, key=lambda z: (z[0] + z[1], z[0], z[1]))
        ],
        "equation_summaries": summaries,
    }


def build_variant(label, support_p, support_q):
    """Build the exact coefficient map and both target systems for a support."""
    t0 = time.perf_counter()
    unknowns = len(support_p) + len(support_q)
    pair_products = _pair_products(support_p, support_q)
    record = {
        "label": label,
        "mode": "dense_bilinear",
        "support_sizes": [len(support_p), len(support_q)],
        "unknowns": unknowns,
        "pair_products": pair_products,
        "degs_PQ": [
            max(a + b for a, b in support_p),
            max(a + b for a, b in support_q),
        ],
    }

    if pair_products > MAX_PAIR_PRODUCTS:
        bound = det_support_bound(support_p, support_q)
        record.update({
            "mode": "support_bound_only",
            "reason": (
                f"{pair_products} derivative pairs exceeds guard "
                f"{MAX_PAIR_PRODUCTS}; not expanded"
            ),
            "det_J": _support_bound_det(support_p, support_q),
            "target_x2": _bound_target_record(
                TARGET_X2, unknowns, bound, f"{label}/detJ-x^2"),
            "target_1": _bound_target_record(
                TARGET_ONE, unknowns, bound, f"{label}/detJ-1"),
            "seconds": round(time.perf_counter() - t0, 3),
        })
        return record

    det, pair_products = det_coefficients(support_p, support_q)
    bound = det_support_bound(support_p, support_q)
    constant_present = (0, 0) in det
    det_record = {
        "mode": "dense_bilinear",
        "monomials": len(det),
        "monomials_are_upper_bound": False,
        "support_bound_monomials": len(bound),
        "cancellations_vs_support_bound": len(bound) - len(det),
        "cancelled_monomials": [
            [int(p), int(q)] for p, q in sorted(bound - set(det))
        ],
        "max_total_degree": max(p + q for p, q in det),
        "min_total_degree": min(p + q for p, q in det),
        "max_total_degree_in_coefficients": 2,
        "max_terms_per_equation": max(len(eq) for eq in det.values()),
        "constant_term_present": bool(constant_present),
        "constant_term_zero": not constant_present,
        "constant_coefficient": (
            _format_coeff(det[(0, 0)]) if constant_present else "0"
        ),
        "pair_products": pair_products,
        "seconds": round(time.perf_counter() - t0, 3),
    }
    record["det_J"] = det_record
    record["target_x2"] = _target_record(
        TARGET_X2, det, det_record, unknowns, f"{label}/detJ-x^2")
    record["target_1"] = _target_record(
        TARGET_ONE, det, det_record, unknowns, f"{label}/detJ-1")
    record["seconds"] = round(time.perf_counter() - t0, 3)
    return record


# --------------------------------------------------------------------------
# checks and report
# --------------------------------------------------------------------------

def _vertices_sympy_cross_check(vtx_record):
    """Independent SymPy check of the 8-unknown vertices-only system."""
    x, y = sp.symbols("x y")
    slow_x2 = size_system(
        list(P_CORNERS), list(Q_CORNERS), x ** 2,
        "vertices_only/detJ-x^2", x, y)
    slow_1 = size_system(
        list(P_CORNERS), list(Q_CORNERS), sp.Integer(1),
        "vertices_only/detJ-1", x, y)

    fast_x2 = vtx_record["target_x2"]
    fast_1 = vtx_record["target_1"]
    checks = {
        "unknowns": slow_x2["unknowns"] == vtx_record["unknowns"] == 8,
        "equations_x2": slow_x2["equations"] == fast_x2["equations"],
        "monomials_x2": (
            {tuple(m) for m in slow_x2["equation_monomials"]}
            == {tuple(m) for m in fast_x2["equation_monomials"]}
        ),
        "equations_1": slow_1["equations"] == fast_1["equations"],
        "monomials_1": (
            {tuple(m) for m in slow_1["equation_monomials"]}
            == {tuple(m) for m in fast_1["equation_monomials"]}
        ),
        "det_constant_term": (
            slow_x2["det_has_constant_term"]
            == vtx_record["det_J"]["constant_term_present"]
        ),
        "max_total_degree": (
            slow_x2["max_total_degree"] == fast_x2["max_total_degree"]
        ),
        "max_degree_in_coefficients": (
            slow_x2["max_total_degree_in_coefficients"]
            == fast_x2["max_total_degree_in_coefficients"]
        ),
        "constant_equation_1": (
            slow_1["constant_equation"]["coefficient"]
            == fast_1["constant_equation"]["coefficient"]
            == "-1"
        ),
        "sample_28_4": (
            next(e["coefficient"] for e in fast_x2["sample_equations"]
                 if e["monomial"] == [28, 4])
            == slow_x2["equations_sample"][0]["coefficient"]
            == "45*a1*b2"
        ),
    }
    return {
        "ok": all(checks.values()),
        "checks": {k: bool(v) for k, v in checks.items()},
        "sympy": {
            "target_x2": {
                "unknowns": slow_x2["unknowns"],
                "equations": slow_x2["equations"],
                "det_monomials": slow_x2["det_monomials"],
                "max_total_degree": slow_x2["max_total_degree"],
                "max_total_degree_in_coefficients": (
                    slow_x2["max_total_degree_in_coefficients"]),
                "constant_equation": slow_x2["constant_equation"],
            },
            "target_1": {
                "unknowns": slow_1["unknowns"],
                "equations": slow_1["equations"],
                "det_monomials": slow_1["det_monomials"],
                "constant_equation": slow_1["constant_equation"],
            },
        },
    }


def _full_sympy_cross_check(full_record, full_p, full_q):
    """Optional slow SymPy check (about 95-105 s measured)."""
    x, y = sp.symbols("x y")
    t0 = time.perf_counter()
    slow = size_system(
        full_p, full_q, x ** 2, "full_lattice/detJ-x^2", x, y)
    seconds = round(time.perf_counter() - t0, 3)
    fast = full_record["target_x2"]
    checks = {
        "unknowns": slow["unknowns"] == full_record["unknowns"] == 445,
        "equations": slow["equations"] == fast["equations"],
        "det_monomials": slow["det_monomials"] == full_record["det_J"]["monomials"],
        "monomials": (
            {tuple(m) for m in slow["equation_monomials"]}
            == {tuple(m) for m in fast["equation_monomials"]}
        ),
        "constant_n_terms": (
            slow["constant_equation"]["n_terms"]
            == fast["constant_equation"]["n_terms"]
        ),
        "max_total_degree": slow["max_total_degree"] == fast["max_total_degree"],
        "max_degree_in_coefficients": (
            slow["max_total_degree_in_coefficients"]
            == fast["max_total_degree_in_coefficients"]
        ),
    }
    return {
        "ran": True,
        "ok": all(checks.values()),
        "checks": {k: bool(v) for k, v in checks.items()},
        "seconds": seconds,
        "sympy": {
            "unknowns": slow["unknowns"],
            "equations": slow["equations"],
            "det_monomials": slow["det_monomials"],
            "max_total_degree": slow["max_total_degree"],
            "max_total_degree_in_coefficients": (
                slow["max_total_degree_in_coefficients"]),
        },
    }


def build_report(slow_cross_check=False):
    t0 = time.perf_counter()
    full_p = lattice_points(P_CORNERS)
    full_q = lattice_points(Q_CORNERS)

    full = build_variant("full_lattice", full_p, full_q)
    vtx = build_variant("vertices_only", list(P_CORNERS), list(Q_CORNERS))

    pick_p = pick_check(P_CORNERS)
    pick_q = pick_check(Q_CORNERS)
    bound_p = det_support_bound(full_p, full_q)

    verification = {
        "pick_P": pick_p,
        "pick_Q": pick_q,
        "lattice_enumeration_matches_pick": (
            pick_p["total_points"] == len(full_p)
            and pick_q["total_points"] == len(full_q)
        ),
        "hull_matches_corners": (
            convex_hull(P_CORNERS) == list(P_CORNERS)
            and convex_hull(Q_CORNERS) == list(Q_CORNERS)
        ),
        "full_det_within_support_bound": (
            full["det_J"]["monomials"] <= len(bound_p)
        ),
        "full_support_bound_monomials": len(bound_p),
        "full_det_cancellations_vs_support_bound": (
            full["det_J"]["cancellations_vs_support_bound"]
        ),
        "vertices_only_sympy_cross_check": _vertices_sympy_cross_check(vtx),
    }
    if slow_cross_check and full["mode"] != "support_bound_only":
        verification["full_lattice_sympy_cross_check"] = _full_sympy_cross_check(
            full, full_p, full_q)
    elif slow_cross_check:
        verification["full_lattice_sympy_cross_check"] = {
            "ran": False,
            "note": "skipped: the full system was not expanded",
        }
    else:
        verification["full_lattice_sympy_cross_check"] = {
            "ran": False,
            "note": (
                "not run by default; --slow-cross-check runs "
                "config2_system.size_system on the full system "
                "(about 95-105 s measured)"
            ),
        }

    theorem = {
        "file": "research/scratch/wstrinz/d2/PROOF_72_108.md",
        "lines": [858, 902, 906, 908, 922],
        "quotes": [
            {
                "line": 858,
                "text": (
                    "Section 11 is a by-product section: 'None is needed for "
                    "Theorem A. All three are stated because they are genuine "
                    "theorems.'"
                ),
            },
            {
                "line": 902,
                "text": (
                    "At $(75,125)$: $80/29$, so $\\alpha = 80$, "
                    "$q_{\\mathrm{window}} = 29$, and the $y$-order cap "
                    "$L(w) = \\lceil \\alpha w/q\\rceil$ is superadditive off "
                    "the period: $L(5) + L(24) = 14 + 67 = 81 > 80 = L(29) = "
                    "\\mathrm{ord}_y\\Phi$, so no relation $c\\Phi = eB$ can "
                    "exist there."
                ),
            },
            {
                "line": 906,
                "text": (
                    "Because $q_{\\mathrm{window}} = 29 = M$ exactly, no split "
                    "with $1 \\le w < M$ is a multiple of $q_{\\mathrm{window}}$, "
                    "so the carry is $1$ on every admissible split and there is "
                    "no weight at which the obstruction could vanish. ... "
                    "So $(75,125)$ is more firmly obstructed by this mechanism "
                    "than the earlier numbers claimed, not less."
                ),
            },
            {
                "line": 908,
                "text": (
                    "Theorem 11.4 (criterion). Let $\\alpha/q$ be "
                    "$\\mathrm{ord}_y(\\Phi)/M$ in lowest terms, so "
                    "$\\gcd(\\alpha,q) = 1$ and $q \\mid M$; let "
                    "$L(w) = \\lceil\\alpha w/q\\rceil$ and "
                    "$\\mathrm{carry}(w_1,w_2) := L(w_1)+L(w_2)-L(w_1+w_2)$. "
                    "For any candidate relation $c\\Phi = eB$ with $c$ a nonzero "
                    "constant and $w(e)+w(B) = M$, $1 \\le w(e) < M$: "
                    "$\\mathrm{carry}(w(e), M-w(e)) = 0 \\iff "
                    "$q_{\\mathrm{window}} \\mid w(e)$, and the carry is exactly "
                    "$0$ or $1$, never more."
                ),
            },
            {
                "line": 922,
                "text": (
                    "Theorem 11.4 is a statement about the obstruction, not "
                    "about existence: it shows that if the relation existed at "
                    "$(75,125)$ then $(75,125)$ would have no lifts at all; the "
                    "load-bearing negative there is a separate exhaustive search "
                    "over exact identities, which is complete for ideal "
                    "identities but not for membership in a radical or a "
                    "saturation."
                ),
            },
        ],
        "classification": "tool",
        "verdict_on_exclusion": "does NOT exclude (75,125)",
        "is_obstruction_to_case": False,
        "excludes_75_125": False,
        "headline": (
            "Theorem 11.4 does NOT exclude (75,125). It is a criterion for the "
            "existence of the K-syzygy relation c*Phi = e*B. At (75,125), "
            "q_window = 29 = M and w(e) = 5, so no admissible split has "
            "carry 0; the relation cannot exist. The theorem therefore removes "
            "the (72,108) syzygy route at this corner; it says nothing about "
            "whether a Jacobian pair with these Newton polygons exists."
        ),
        "role": (
            "negative/diagnostic tool for the K-syzygy mechanism, not an "
            "obstruction to the (75,125) coefficient system"
        ),
        "what_it_does_not_say": (
            "It does not say det J = x^2 is infeasible, and it does not "
            "exclude any (P,Q) with N(P)=3*Delta', N(Q)=5*Delta'."
        ),
        "research_note_claim": (
            "research/jacobian-degree-pairs-ge125.md:81 calls it 'a proved "
            "obstruction'. That is accurate only for the K-syzygy mechanism; "
            "it is not an exclusion of (75,125)."
        ),
        "repo_status": {
            "file": "research/scratch/wstrinz/README.md",
            "lines": "201-203",
            "quote": (
                "Theorem 2.1 is strict (max >= 125), so F_2(3,5)/125 sits *at* "
                "the bound and correctly stays open -- it is the unique row "
                "there, and the next increment of the bound depends on it "
                "alone."
            ),
            "meaning": (
                "The wstrinz repo itself lists (75,125) as open, not excluded."
            ),
        },
    }

    full_x2 = full["target_x2"]
    full_1 = full["target_1"]
    vtx_x2 = vtx["target_x2"]
    vtx_1 = vtx["target_1"]

    report = {
        "problem": "jacobian-75-125",
        "case": "F_2(3,5), j=1",
        "source": {
            "corner_file": (
                "research/scratch/wstrinz/d2/PASSPORT_75_125_REPAIR.md"
            ),
            "corner_lines": "29-33",
            "corner_quote": (
                "Delta' = {(0,0),(3,0),(4,1),(0,5)}; "
                "N(P) = 3*Delta' = {(0,0),(9,0),(12,3),(0,15)}; "
                "N(Q) = 5*Delta' = {(0,0),(15,0),(20,5),(0,25)}; "
                "[P,Q] = x^2; reduced degrees (15,25)"
            ),
            "dossier": (
                "research/scratch/wstrinz/d2/case_dossier_F2_j1_75_125.json"
            ),
            "dossier_case": {
                "A0": [5, 20],
                "A0prime": [1, 0],
                "A1": "(7\\5,2)",
                "a": 3,
                "b": 5,
                "m": 3,
                "n": 5,
                "k": 1,
                "j": 1,
                "v11": 25,
                "t": 4,
                "kappa": 2,
                "q": 1,
            },
        },
        "support_convention": {
            "choice": "full lattice points of the Newton polygon (generic)",
            "why": (
                "The corner data gives the Newton polygons N(P)=3*Delta' and "
                "N(Q)=5*Delta'. A generic polynomial with a given Newton "
                "polygon has every lattice point of the closed polygon in its "
                "support; vertices-only is a minimal-support specialization, "
                "reported as a separate variant."
            ),
            "vertices_are_supports": False,
            "boundary_included": True,
        },
        "supports": {
            "Delta_prime": [list(p) for p in DELTA_PRIME],
            "N(P)_corners": [list(p) for p in P_CORNERS],
            "N(Q)_corners": [list(p) for p in Q_CORNERS],
            "N(P)_hull": [list(p) for p in convex_hull(P_CORNERS)],
            "N(Q)_hull": [list(p) for p in convex_hull(Q_CORNERS)],
            "lattice_points_P": len(full_p),
            "lattice_points_Q": len(full_q),
            "lattice_points_total": len(full_p) + len(full_q),
            "interior_points_P": pick_p["interior_points"],
            "boundary_points_P": pick_p["boundary_points"],
            "interior_points_Q": pick_q["interior_points"],
            "boundary_points_Q": pick_q["boundary_points"],
            "interior_points_total": (
                pick_p["interior_points"] + pick_q["interior_points"]
            ),
            "boundary_points_total": (
                pick_p["boundary_points"] + pick_q["boundary_points"]
            ),
            "support_sha256_P": _points_hash(full_p),
            "support_sha256_Q": _points_hash(full_q),
            "symbols": {
                "P": "a0..a123",
                "Q": "b0..b320",
                "note": (
                    "0-based, matching config2_system.build_polynomial"
                ),
            },
        },
        "degrees": {
            "deg_P_total": max(a + b for a, b in P_CORNERS),
            "deg_Q_total": max(a + b for a, b in Q_CORNERS),
            "reduced_degrees": [15, 25],
            "detJ_max_total_degree": full["det_J"]["max_total_degree"],
            "detJ_min_total_degree": full["det_J"]["min_total_degree"],
            "target": "x^2",
        },
        "unknowns": full["unknowns"],
        "full_unknowns": full["unknowns"],
        "vtx_unknowns": vtx["unknowns"],
        "full_lattice_unknowns": full["unknowns"],
        "vertices_only_unknowns": vtx["unknowns"],
        "equations": full_x2["equations"],
        "full_equations": full_x2["equations"],
        "vtx_equations": vtx_x2["equations"],
        "thm114": theorem["classification"],
        "thm114_excludes_75_125": theorem["excludes_75_125"],
        "max_total_degree": full_x2["max_total_degree"],
        "max_total_degree_in_coefficients": (
            full_x2["max_total_degree_in_coefficients"]),
        "determinedness": full_x2["determinedness"],
        "effective": {
            "unknowns": full_x2["effective_unknowns"],
            "note": "the additive constants of P and Q never enter det J",
            "nontrivial_equations": full_x2["nontrivial_equations"],
            "equations_minus_effective_unknowns": (
                full_x2["equations_minus_effective_unknowns"]
            ),
        },
        "feasibility": {
            "det_J_equals_1": {
                "full_lattice": {
                    "constant_equation": full_1["constant_equation"],
                    "verdict": (
                        "not evaluated (support-bound mode)"
                        if full_1["constant_term_test_infeasible"] is None else
                        "not excluded by the constant-term test: "
                        "c00 - 1 = 0 is a genuine equation, with "
                        "c00 = -a1*b26 + a16*b1"
                    ),
                    "infeasible_by_constant_term": (
                        full_1["constant_term_test_infeasible"]
                    ),
                },
                "vertices_only": {
                    "constant_equation": vtx_1["constant_equation"],
                    "verdict": (
                        "not evaluated (support-bound mode)"
                        if vtx_1["constant_term_test_infeasible"] is None else
                        "infeasible: det J has no constant term, so the "
                        "constant equation is -1 = 0"
                    ),
                    "infeasible_by_constant_term": (
                        vtx_1["constant_term_test_infeasible"]
                    ),
                },
                "note": (
                    "The paper normalization is [P,Q] = x^2, not 1; this "
                    "entry is only a support-level comparison."
                ),
            },
            "det_J_equals_x2": {
                "full_lattice": {
                    "det_constant_term_zero": full_x2[
                        "constant_term_of_det_zero"],
                    "constant_term_test": full_x2["constant_term_test"],
                    "target_monomial_test": full_x2["target_monomial_test"],
                    "verdict": (
                        "not evaluated (support-bound mode)"
                        if full_x2["infeasible_by_target_monomial"] is None else
                        "not excluded by the constant-term test: c00 = 0 is a "
                        "genuine equation and (2,0) is present"
                    ),
                    "infeasible_by_target_monomial": (
                        full_x2["infeasible_by_target_monomial"]
                    ),
                },
                "vertices_only": {
                    "det_constant_term_zero": vtx_x2[
                        "constant_term_of_det_zero"],
                    "constant_term_test": vtx_x2["constant_term_test"],
                    "target_monomial_test": vtx_x2["target_monomial_test"],
                    "verdict": (
                        "not evaluated (support-bound mode)"
                        if vtx_x2["infeasible_by_target_monomial"] is None else
                        "infeasible: the x^2 monomial is absent from det J, so "
                        "the (2,0) equation is -1 = 0"
                    ),
                    "infeasible_by_target_monomial": (
                        vtx_x2["infeasible_by_target_monomial"]
                    ),
                },
            },
            "detJ_minus_x2_feasible_by_constant_term": {
                "full_lattice": (
                    None if full_x2["constant_term_test_infeasible"] is None
                    else not full_x2["constant_term_test_infeasible"]
                ),
                "vertices_only": (
                    None if vtx_x2["constant_term_test_infeasible"] is None
                    else not vtx_x2["constant_term_test_infeasible"]
                ),
                "note": (
                    "For vertices_only the pure constant-term test is vacuous "
                    "(0 = 0); the system is excluded instead by the x^2 "
                    "monomial test."
                ),
            },
        },
        "variants": {
            "full_lattice": full,
            "vertices_only": vtx,
        },
        "summary": {
            "full_unknowns": full["unknowns"],
            "vtx_unknowns": vtx["unknowns"],
            "full_equations_x2": full_x2["equations"],
            "vtx_equations_x2": vtx_x2["equations"],
            "full_equations_1": full_1["equations"],
            "vtx_equations_1": vtx_1["equations"],
            "full_max_total_degree": full_x2["max_total_degree"],
            "vtx_max_total_degree": vtx_x2["max_total_degree"],
            "full_max_coefficient_degree": (
                full_x2["max_total_degree_in_coefficients"]),
            "vtx_max_coefficient_degree": (
                vtx_x2["max_total_degree_in_coefficients"]),
            "full_det_constant_term_zero": full["det_J"][
                "constant_term_zero"],
            "vtx_det_constant_term_zero": vtx["det_J"][
                "constant_term_zero"],
            "full_detJ_minus_x2_infeasible": (
                full_x2["infeasible_by_target_monomial"]
            ),
            "vtx_detJ_minus_x2_infeasible": (
                vtx_x2["infeasible_by_target_monomial"]
            ),
        },
        "theorem_11_4": theorem,
        "verification": verification,
        "timings": {
            "full_lattice_seconds": full["seconds"],
            "vertices_only_seconds": vtx["seconds"],
            "expanded_dense": True,
            "blow_up": False,
            "note": (
                "The exact bilinear accumulator expands the full 445-unknown "
                "system in well under a second; the measured SymPy expansion "
                "of the same determinant is about 13 s, and the full SymPy "
                "coefficient extraction is about 95-105 s (optional "
                "--slow-cross-check)."
            ),
        },
        "notes": [
            "No solving, elimination, Groebner basis, or collision search was "
            "performed.",
            "det J is bilinear in the coefficients, so every nonconstant "
            "coefficient equation has total degree 2 in the unknowns.",
            "The full-lattice system is overdetermined: 723 equations in 445 "
            "unknowns, or 723 in 443 after the two additive constants are "
            "removed.",
            "The vertices-only system is 8 unknowns / 7 equations, but it is "
            "infeasible for det J = x^2 because the x^2 monomial is absent "
            "from det J; this is a support obstruction, not a Groebner "
            "obstruction.",
            "No nonvanishing constraint is imposed on the vertex coefficients; "
            "solutions whose support drops to a smaller Newton polygon are "
            "included in the generic ansatz as boundary strata.",
            "Theorem 11.4 is not an exclusion of (75,125); see "
            "theorem_11_4 above.",
        ],
        "elapsed_seconds": None,
    }
    report["elapsed_seconds"] = round(time.perf_counter() - t0, 3)
    return report


def print_report(rep):
    full = rep["variants"]["full_lattice"]
    vtx = rep["variants"]["vertices_only"]
    full_x2 = full["target_x2"]
    vtx_x2 = vtx["target_x2"]
    print("== Jacobian corner (75,125), F_2(3,5), j=1 ==")
    print("Delta'          : {(0,0),(3,0),(4,1),(0,5)}")
    print("N(P)=3*Delta'   : {(0,0),(9,0),(12,3),(0,15)}  "
          "(%d lattice points)" % rep["supports"]["lattice_points_P"])
    print("N(Q)=5*Delta'   : {(0,0),(15,0),(20,5),(0,25)}  "
          "(%d lattice points)" % rep["supports"]["lattice_points_Q"])
    print("normalisation   : [P,Q] = x^2 ; reduced degrees (15,25)")
    print()
    print("full lattice    : %d unknowns, %d equations (det J - x^2)"
          % (full_x2["unknowns"], full_x2["equations"]))
    print("                  max total degree %s (x,y); %s (coefficients); "
          "%s terms max/equation"
          % (full_x2["max_total_degree"],
             full_x2["max_total_degree_in_coefficients"],
             full_x2["max_terms_per_equation"]))
    print("                  %s (raw); %s (effective, %d unknowns)"
          % (full_x2["determinedness"],
             full_x2["determinedness_effective"],
             full_x2["effective_unknowns"]))
    print("vertices only   : %d unknowns, %d equations (det J - x^2)"
          % (vtx_x2["unknowns"], vtx_x2["equations"]))
    print("                  max total degree %d (x,y); %d (coefficients)"
          % (vtx_x2["max_total_degree"],
             vtx_x2["max_total_degree_in_coefficients"]))
    print()
    print("constant term   : full det J = %s"
          % full["det_J"]["constant_coefficient"])
    print("                  (zero? %s); vertices det J zero? %s"
          % (full["det_J"]["constant_term_zero"],
             vtx["det_J"]["constant_term_zero"]))
    print("det J - x^2     : full = %s"
          % rep["feasibility"]["det_J_equals_x2"]["full_lattice"]["verdict"])
    print("                  vertices = %s"
          % rep["feasibility"]["det_J_equals_x2"]["vertices_only"]["verdict"])
    print()
    print("Theorem 11.4    : NOT an exclusion of (75,125). It rules out the")
    print("                  K-syzygy relation c*Phi = e*B at this corner")
    print("                  (q_window = 29 = M); the coefficient system")
    print("                  remains live.")
    print()
    print("elapsed         : %.3f s" % rep["elapsed_seconds"])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--report", action="store_true",
                    help="print the size report (default action)")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT,
                    help="JSON output path")
    ap.add_argument("--no-write", action="store_true",
                    help="do not write the JSON")
    ap.add_argument("--slow-cross-check", action="store_true",
                    help="also run the ~100 s full SymPy cross-check")
    args = ap.parse_args(argv)

    t0 = time.perf_counter()
    rep = build_report(slow_cross_check=args.slow_cross_check)
    rep["elapsed_seconds"] = round(time.perf_counter() - t0, 3)
    print_report(rep)

    if not args.no_write:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(rep, indent=2), encoding="utf-8")
        try:
            shown = args.out.relative_to(ROOT)
        except ValueError:
            shown = args.out
        print("\nwritten: %s" % shown)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
