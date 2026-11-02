"""Full edge cascade at the (75,125) Jacobian corner, F_2(3,5), j=1.

WHAT THIS SCRIPT ADDS OVER ``cascade75.py``
-------------------------------------------
``cascade75.py`` descends three of the four hull-edge conditions and stops.
This script descends *all* of them, decides the dimension at every step, and
records the first empty condition (there is none).  It also audits the premise
that the corner carries "16 edge conditions", because that premise is false as
stated and a cascade run against 12 non-existent conditions would report a
spurious emptiness.

THE EDGE INVENTORY (stage 1)
----------------------------
N(P) = 3*Delta' and N(Q) = 5*Delta' are homothetic, so they have the *same*
normal fan: four rays.  A leading-form condition LP_w(det J) = det J(P_w,Q_w)
needs one weight w that is maximal on a face of N(P) *and* on a face of N(Q).
On a 2-dimensional polygon the normal cone of an edge is the ray through its
primitive outward normal, so exactly four pairs of faces qualify: the four
matched edge slots, with normals

    (1,1), (1,-1), (0,-1), (-1,0)                       [cascade order]

The 16 = 4 x 4 reading of "all edge conditions" is the only one consistent
with the recorded order (edge 1 = (1,1), edge 2 = (1,-1), edge 3 = (0,-1)).
Under that reading the 16 conditions are the pairs (edge i of N(P), edge j of
N(Q)).  The four diagonal pairs are the genuine slots.  The twelve
off-diagonal pairs are *not* leading-form conditions: their normal cones are
disjoint, so no w exposes both faces, and, concretely, they fail at an
explicit point that satisfies every genuine edge condition.  Section
``rejected_candidate_conditions`` records both certificates.  Imposing them
would manufacture an emptiness that is not there, so they are rejected rather
than descended.

THE FOUR GENUINE CONDITIONS
---------------------------
For an edge with normal w, d_P = w.max(N(P)), d_Q = w.max(N(Q)) and edge
forms P_e, Q_e,

    LP_w(det J) = det J(P_e, Q_e),                                           (1)

and for the two diagonal edges (where both forms are univariate in the same
t = x^u y^v and the monomial prefactor cancels) this is the Wronskian
identity already certified in ``cascade75.py``:

    det J(P_e,Q_e) = monomial * (d_Q f g' - d_P f' g).                       (2)

Equation (2) with f,g not identically zero is equivalent to

    f^(d_Q) = C g^(d_P)   =>   f = c1 k^(d_P/g),  g = c2 k^(d_Q/g),          (3)

g = gcd(d_P,d_Q), deg k = (n_P - 1) g / d_P.  Both diagonal edges here have

    (d_P/g, d_Q/g) = (3, 5),

so *every* diagonal edge forces the same (3,5) power structure.  The two
remaining edges are horizontal/vertical; both edge forms are univariate in
the same variable x (resp. y), so (1) is identically zero and the condition
is vacuous -- a pass, not a hidden failure.  The next filtration level after
a vacuous edge is checked explicitly in ``next_level_after_vacuous``: the
linear condition there is satisfied by (P_1,Q_1) = (P_0',Q_0'), which fits
inside the Newton polygon, so the vacuous edges hide no emptiness either
order.

WHY THE CASCADE CANNOT BE EMPTY HERE (the structural reason)
------------------------------------------------------------
(3) always has solutions: pick any k of the stated degree and any nonzero
c1,c2.  So a single edge condition is never empty, and the only way the
cascade could fail is through incompatibility between edges.  Consecutive
edges share exactly one vertex per polygon, and the family value at that
vertex is c_i * k(vertex)^power, which can be matched by rescaling c_i.  The
four edge conditions therefore cut out a non-empty variety: the two diagonal
edges remove 7 + 28 dimensions and the two vacuous edges remove none.  The
computed dimension of that variety is reported in ``joint_solution_set``
(410 scaled, 408 in the strict chart).  This is a statement about the
*filter*, not about det J = x^2: every edge condition is necessary, none is
sufficient.

DIMENSION CONVENTIONS
---------------------
* ``scaled``  : the coefficient space of the full lattice (445 unknowns).
* ``strict``  : the strict Newton-polygon chart used by ``cascade75.py``,
                which fixes the two vertex coefficients a_(12,3) = b_(20,5) = 1
                (they are vertices of N(P), N(Q), hence nonzero).  The edge
                systems are homogeneous of degree 1 in each block, so this
                removes exactly the two scalings.

Run:
    python problems/jacobian/cascade75_full.py
"""

from __future__ import annotations

import argparse
import json
import math
import random
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import cascade75 as c75  # noqa: E402
from cascade75 import (  # noqa: E402
    X,
    Y,
    build_decomposition,
    build_edge_system,
    certificate_leading_form,
    edge_leading_form_string,
    groebner_capped,
    krull_dimension_lm,
    shared_vertices,
    sym_a,
    sym_b,
)

DEFAULT_OUT = ROOT / "runs" / "jacobian_75_125_cascade_full.json"
LOG_PATH = ROOT / "runs" / "cascade75_full.log"
PRIOR_CASCADE = ROOT / "runs" / "jacobian_75_125_cascade.json"

EDGE_CAP_S = 180.0          # hard cap per edge, as required
TOP_GB_CAP_S = 60.0         # the top edge's direct Groebner is known to be hard
MODULUS = 32003
DECLARED_CONDITIONS = 16    # the number asserted by the task
SEED = 20260101


# --------------------------------------------------------------------------
# logging
# --------------------------------------------------------------------------

def log_stage(stage, payload):
    """Append exactly one line per stage to runs/cascade75_full.log."""
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    line = "%s stage=%s %s\n" % (
        stamp, stage, json.dumps(payload, sort_keys=True, default=str))
    with LOG_PATH.open("a", encoding="utf-8") as fh:
        fh.write(line)


# --------------------------------------------------------------------------
# stage 1: premise audit -- how many edge conditions are there?
# --------------------------------------------------------------------------

def _cross(a, b):
    return int(a[0]) * int(b[1]) - int(a[1]) * int(b[0])


def premise_audit(decomp):
    """Enumerate every face-pair that can carry a leading form.

    A leading form at weight w is det J(P_w,Q_w).  For an edge face the
    normal cone is the ray through the primitive outward normal, so two edges
    can appear together in one leading form iff their normals are equal.  The
    polygons are homothetic, so the normal fans coincide and the matched pairs
    are exactly the four slots.
    """
    slots = decomp["slots"]
    order = decomp["order"]
    edges_p = decomp["edges_P"]
    edges_q = decomp["edges_Q"]

    normals_p = [tuple(e["normal"]) for e in edges_p]
    normals_q = [tuple(e["normal"]) for e in edges_q]
    if normals_p != normals_q:
        raise AssertionError("normal fans differ: %s vs %s"
                             % (normals_p, normals_q))

    grid = []
    for i in range(len(edges_p)):
        for j in range(len(edges_q)):
            wp, wq = normals_p[i], normals_q[j]
            same_ray = (_cross(wp, wq) == 0 and
                        wp[0] * wq[0] + wp[1] * wq[1] > 0)
            grid.append({
                "edge_P": i,
                "edge_Q": j,
                "normal_P": list(wp),
                "normal_Q": list(wq),
                "paired_by_a_common_weight": bool(same_ray),
                "classification": ("edge condition" if same_ray
                                   else "invalid (normal cones disjoint)"),
            })
    valid = [g for g in grid if g["paired_by_a_common_weight"]]
    invalid = [g for g in grid if not g["paired_by_a_common_weight"]]

    # vertex faces: the four matched vertices.  N(P) = 3 Delta', N(Q) = 5
    # Delta', so vertex v of N(Q) is (5/3) v of N(P) and
    # det J(x^a y^b, x^c y^d) = (ad - bc) x^(a+c-1) y^(b+d-1) = 0.
    vertex_pairs = []
    corners_p = [tuple(e["start"]) for e in edges_p]
    corners_q = [tuple(e["start"]) for e in edges_q]
    for cp in corners_p:
        cand = [cq for cq in corners_q
                if 5 * cp[0] == 3 * cq[0] and 5 * cp[1] == 3 * cq[1]]
        if not cand:
            raise AssertionError("no proportional vertex for %s" % (cp,))
        cq = cand[0]
        a, b = sp.Symbol("vp_a"), sp.Symbol("vq_b")
        det = sp.expand(
            sp.diff(a * X ** cp[0] * Y ** cp[1], X)
            * sp.diff(b * X ** cq[0] * Y ** cq[1], Y)
            - sp.diff(a * X ** cp[0] * Y ** cp[1], Y)
            * sp.diff(b * X ** cq[0] * Y ** cq[1], X))
        vertex_pairs.append({
            "P": list(cp),
            "Q": list(cq),
            "proportional": True,
            "det_J_identically_zero": bool(det == 0),
        })

    return {
        "declared_in_task": DECLARED_CONDITIONS,
        "geometric_edge_slots": len(slots),
        "slots": [{
            "edge": k + 1,
            "hull_index": order[k],
            "normal": slots[order[k]]["normal"],
            "d_P": slots[order[k]]["d_P"],
            "d_Q": slots[order[k]]["d_Q"],
            "filtration_degree": slots[order[k]]["filtration_degree"],
            "n_points_P": slots[order[k]]["n_points_P"],
            "n_points_Q": slots[order[k]]["n_points_Q"],
        } for k in range(len(order))],
        "normal_fans_identical": True,
        "edge_pair_grid_16": grid,
        "genuine_edge_conditions": len(valid),
        "invalid_off_diagonal_pairs": len(invalid),
        "vertex_face_conditions": vertex_pairs,
        "vertex_face_conditions_all_vacuous": all(
            v["det_J_identically_zero"] for v in vertex_pairs),
        "readings_of_16": [
            {
                "reading": "4 P-edges x 4 Q-edges = 16 edge pairs",
                "consistent_with_recorded_order": True,
                "valid_conditions": 4,
                "invalid_conditions": 12,
                "why": "only the diagonal pairs share a normal; the other "
                       "twelve have disjoint normal cones",
            },
            {
                "reading": "4 edges + 4 vertices per polygon, 2 polygons "
                           "= 16 faces",
                "consistent_with_recorded_order": False,
                "valid_conditions": 4,
                "vacuous_conditions": 4,
                "invalid_conditions": 8,
                "why": "faces pair only across the two polygons by a common "
                       "weight; the four matched vertex pairs are vacuous and "
                       "an unmatched face has no common weight",
            },
            {
                "reading": "4 edges x 4 coefficient layers",
                "consistent_with_recorded_order": False,
                "valid_conditions": None,
                "why": "the recorded weights (1,1), (1,-1), (0,-1) are edge "
                       "normals, not layers of one grading; the layered "
                       "descent is a different (sub-leading) filtration and "
                       "has 39 levels for w=(1,1), not 4",
            },
        ],
        "finding": (
            "the corner has 4 genuine edge conditions, not 16.  The 4 x 4 "
            "reading (16 pairs of a P-edge with a Q-edge) contains 4 valid "
            "diagonal pairs and 12 off-diagonal pairs whose normal cones are "
            "disjoint, so no single weight exposes both faces and "
            "det J(P_ei,Q_ej) = 0 is not a leading-form condition.  The 4 "
            "matched vertex pairs are also vacuous because the vertices are "
            "proportional.  See rejected_candidate_conditions for a witness "
            "where the off-diagonal determinants fail at a point that "
            "satisfies every genuine edge condition."),
    }


# --------------------------------------------------------------------------
# edge families: structure theorem, explicit rational point, certificates
# --------------------------------------------------------------------------

def edge_family_shape(slot):
    """(g, d_P/g, d_Q/g, deg k) for a diagonal edge; None for a vacuous one."""
    if slot["leading_form_is_identically_zero"]:
        return None
    d_p, d_q = int(slot["d_P"]), int(slot["d_Q"])
    g = math.gcd(d_p, d_q)
    kmax = (len(slot["points_P"]) - 1) * g // d_p
    if kmax != (len(slot["points_Q"]) - 1) * g // d_q:
        raise AssertionError("power exponents disagree")
    return {
        "g": int(g),
        "p_power": int(d_p // g),
        "q_power": int(d_q // g),
        "deg_k": int(kmax),
    }


def _family_polynomials(slot, k_coeffs, c1, c2):
    """P_e = c1 G^p_power, Q_e = c2 G^q_power, G = monomial * k(t), k(0)=1."""
    shape = edge_family_shape(slot)
    if shape is None:
        raise ValueError("vacuous edge has no power family")
    p0, q0 = slot["points_P"][0]
    r0, s0 = slot["points_Q"][0]
    u, v = slot["edge_P"]["direction"]
    t = X ** int(u) * Y ** int(v)
    k = 1 + sum(k_coeffs[i] * t ** (i + 1) for i in range(len(k_coeffs)))
    gp = X ** (int(p0) // shape["p_power"]) * Y ** (int(q0) // shape["p_power"])
    gq = X ** (int(r0) // shape["q_power"]) * Y ** (int(s0) // shape["q_power"])
    return (sp.expand(c1 * (gp * k) ** shape["p_power"]),
            sp.expand(c2 * (gq * k) ** shape["q_power"]))


def _coeffs_on_edge(expr, points):
    """Coefficients at ``points``, in the order the points are given."""
    poly = sp.Poly(expr, X, Y)
    want = {tuple(int(c) for c in p) for p in points}
    extra = set(poly.monoms()) - want
    if extra:
        raise AssertionError("family leaves the edge support: %s" % (extra,))
    return [poly.coeff_monomial(X ** int(p[0]) * Y ** int(p[1]))
            for p in points]


def family_point(slot, k_coeffs, c1, c2):
    """Explicit rational coefficient assignment on a diagonal edge."""
    p_poly, q_poly = _family_polynomials(slot, k_coeffs, c1, c2)
    vals = {}
    for sym, val in zip([sym_a(*p) for p in slot["points_P"]],
                        _coeffs_on_edge(p_poly, slot["points_P"])):
        vals[str(sym)] = sp.nsimplify(val)
    for sym, val in zip([sym_b(*p) for p in slot["points_Q"]],
                        _coeffs_on_edge(q_poly, slot["points_Q"])):
        vals[str(sym)] = sp.nsimplify(val)
    return vals


def edge_family_certificate(slot, rng):
    """Exact point + symbolic det J check.  Returns the certificate record."""
    shape = edge_family_shape(slot)
    if shape is None:
        return {
            "applicable": False,
            "reason": "vacuous edge: no leading-form equation to solve",
        }
    k_coeffs = [sp.Integer(rng.randint(-7, 7)) for _ in range(shape["deg_k"])]
    c1, c2 = sp.Integer(2), sp.Integer(3)
    vals = family_point(slot, k_coeffs, c1, c2)
    system = build_edge_system(slot["points_P"], slot["points_Q"])
    sub = {sp.Symbol(s): v for s, v in vals.items()}
    residuals = [sp.expand(eq.subs(sub)) for eq in system["equations"]]
    nonzero = sum(1 for r in residuals if r != 0)
    symbolic = c75._edge_family_check(slot)
    return {
        "applicable": True,
        "family": (
            "P_e = c1 G^%d, Q_e = c2 G^%d, G = x^%d*y^%d*k(t), t = x^%d*y^%d, "
            "k(0) = 1, deg k = %d"
            % (shape["p_power"], shape["q_power"],
               slot["points_P"][0][0] // shape["p_power"],
               slot["points_P"][0][1] // shape["p_power"],
               slot["edge_P"]["direction"][0], slot["edge_P"]["direction"][1],
               shape["deg_k"])),
        "parameters_scaled": shape["deg_k"] + 2,
        "parameters_strict": shape["deg_k"] + 2 - 2,
        "trial_point_coefficients": {
            "k": [int(c) for c in k_coeffs],
            "c1": int(c1), "c2": int(c2)},
        "equations_checked": len(residuals),
        "nonzero_residuals": int(nonzero),
        "all_equations_vanish": bool(nonzero == 0),
        "det_J_identically_zero_for_symbolic_k": bool(
            symbolic["det_J_identically_zero"]),
        "symbolic_family": symbolic["family"],
        "lower_bound": "the explicit family is a %d-parameter solution set, "
                       "so the edge variety is non-empty of dimension >= %d"
                       % (shape["deg_k"] + 2, shape["deg_k"] + 2),
        "upper_bound": (
            "det J(P_e,Q_e) = monomial * (d_Q f g' - d_P f' g); dividing by "
            "f g (f,g != 0) gives (f^(d_Q) / g^(d_P))' = 0, hence "
            "f^(d_Q) = C g^(d_P); unique factorisation gives f = c1 k^%d, "
            "g = c2 k^%d with deg k <= %d, i.e. exactly the family above"
            % (shape["p_power"], shape["q_power"], shape["deg_k"])),
    }


# --------------------------------------------------------------------------
# stage 2: descend the four genuine edge conditions
# --------------------------------------------------------------------------

def _prior_top_groebner():
    """Reuse the recorded top-edge Groebner attempt (it is expensive)."""
    if not PRIOR_CASCADE.exists():
        return None
    try:
        prior = json.loads(PRIOR_CASCADE.read_text(encoding="utf-8"))
        gb = prior["top_edge"]["groebner_over_Q"]
    except (KeyError, ValueError):
        return None
    if gb.get("status") != "timeout":
        return None
    return {
        "status": "timeout",
        "cap_seconds": gb.get("cap_seconds"),
        "seconds": gb.get("seconds"),
        "basis_size": None,
        "dimension": None,
        "empty": None,
        "carried_from": str(PRIOR_CASCADE.relative_to(ROOT)),
        "note": "a 34-variable direct Groebner basis was attempted in the "
                "previous run and exceeded the cap; it is not repeated here "
                "because the dimension is settled exactly by the structure "
                "theorem plus the explicit family certificate",
    }


def analyse_edge(slot, index, previous_slot, cap, fixed_symbols, rng,
                 attempt_groebner=True):
    """One cascade step: leading form, dimension, non-emptiness certificate."""
    pts_p, pts_q = slot["points_P"], slot["points_Q"]
    system = build_edge_system(pts_p, pts_q)
    n_unknowns = len(system["a"]) + len(system["b"])
    record = {
        "edge": int(index),
        "normal": slot["normal"],
        "d_P": int(slot["d_P"]),
        "d_Q": int(slot["d_Q"]),
        "filtration_degree": int(slot["filtration_degree"]),
        "target_weight": int(slot["target_weight"]),
        "n_points_P": len(pts_p),
        "n_points_Q": len(pts_q),
        "unknowns": int(n_unknowns),
        "equations": len(system["equations"]),
        "leading_form": edge_leading_form_string(slot),
        "leading_form_identically_zero": bool(
            slot["leading_form_is_identically_zero"]),
        "groebner_cap_seconds": cap,
    }

    shared_p, shared_q = [], []
    if previous_slot is not None:
        prev_p = set(map(tuple, previous_slot["points_P"]))
        prev_q = set(map(tuple, previous_slot["points_Q"]))
        shared_p = sorted(set(map(tuple, pts_p)) & prev_p)
        shared_q = sorted(set(map(tuple, pts_q)) & prev_q)
    record["shared_points_with_previous"] = {
        "P": [list(p) for p in shared_p],
        "Q": [list(q) for q in shared_q],
    }

    if not system["equations"]:
        record.update({
            "vacuous": True,
            "reason": (
                "both edge forms are univariate in the same variable (%s "
                "edge), so det J(P_e,Q_e) = 0 identically for every choice of "
                "edge coefficients"
                % ("horizontal" if slot["edge_P"]["is_horizontal"]
                   else "vertical")),
            "dimension": int(n_unknowns),
            "dimension_scaled": int(n_unknowns),
            "dimension_method": "no equation: every edge coefficient is free",
            "empty": False,
            "non_empty_certificate": "the edge system has 0 equations",
            "groebner": {"status": "not_run",
                         "reason": "edge system is empty (0 equations)"},
        })
        return record

    shape = edge_family_shape(slot)
    record["vacuous"] = False
    record["power_structure"] = shape
    record["dimension_scaled"] = int(shape["deg_k"] + 2)
    record["dimension"] = int(shape["deg_k"] + 2 - 2)
    record["dimension_method"] = (
        "structure theorem (Wronskian + unique factorisation) and an explicit "
        "%d-parameter family; strict chart removes the two vertex scalings"
        % (shape["deg_k"] + 2))

    record["leading_form_identity"] = certificate_leading_form(slot)
    record["family_certificate"] = edge_family_certificate(slot, rng)
    record["empty"] = False
    record["non_empty_certificate"] = (
        "explicit rational point (family_certificate) satisfies every "
        "equation exactly; nonzero_residuals = 0")

    if attempt_groebner and n_unknowns <= 12:
        gb = groebner_capped(slot["normal"], fixed_symbols, cap)
        if gb.get("status") == "ok":
            dim = krull_dimension_lm(gb["lm_supports"], gb["n_vars"])
            gb["dimension"] = dim["dimension"]
            gb["empty"] = bool(dim["empty"])
            gb["dimension_detail"] = dim
            gb["agrees_with_structure_theorem"] = bool(
                dim["dimension"] == record["dimension"])
        record["groebner"] = gb
    elif attempt_groebner and index == 1:
        carried = _prior_top_groebner()
        record["groebner"] = carried if carried is not None else {
            "status": "not_run",
            "reason": "34 variables: direct Groebner is not a small system",
        }
    else:
        record["groebner"] = {
            "status": "not_run",
            "reason": "system too large for a direct basis; the dimension is "
                      "decided by the structure theorem"}

    return record


def descend_all(decomp, cap, rng):
    slots, order = decomp["slots"], decomp["order"]
    fix = {"a_12_3": 1, "b_20_5": 1}
    records = []
    previous = None
    for step in range(len(order)):
        slot = slots[order[step]]
        rec = analyse_edge(slot, step + 1, previous, cap, fix, rng,
                           attempt_groebner=True)
        records.append(rec)
        log_stage("edge_%d" % (step + 1), {
            "normal": rec["normal"],
            "unknowns": rec["unknowns"],
            "equations": rec["equations"],
            "vacuous": rec.get("vacuous", False),
            "dimension": rec.get("dimension"),
            "dimension_scaled": rec.get("dimension_scaled"),
            "groebner": rec.get("groebner", {}).get("status"),
            "empty": rec.get("empty"),
        })
        previous = slot
    return records


# --------------------------------------------------------------------------
# stage 3: the joint solution set of all edge conditions
# --------------------------------------------------------------------------

def joint_solution_set(decomp, records):
    """Dimension of the variety cut out by *all* edge conditions.

    The four conditions live on the four edge supports.  Distinct edges share
    exactly one vertex per polygon, so the conditions are independent except
    for two shared coefficients, a_(12,3) and b_(20,5).  Each diagonal edge
    contributes a power family of dimension deg k + 2; matching the two
    shared vertices removes exactly 2 (the two scalar rescalings), and the
    coefficients off every edge are untouched.
    """
    pts_p = {("P",) + tuple(p) for p in decomp["lattice_points_P"]}
    pts_q = {("Q",) + tuple(p) for p in decomp["lattice_points_Q"]}
    boundary_p = set()
    boundary_q = set()
    for e in decomp["edges_P"]:
        boundary_p |= {("P",) + tuple(p) for p in e["points"]}
    for e in decomp["edges_Q"]:
        boundary_q |= {("Q",) + tuple(p) for p in e["points"]}
    interior = (pts_p - boundary_p) | (pts_q - boundary_q)
    boundary = boundary_p | boundary_q

    constrained = set()
    for rec in records:
        if rec.get("vacuous"):
            continue
        slot = decomp["slots"][decomp["order"][rec["edge"] - 1]]
        for p in slot["points_P"]:
            constrained.add(("a",) + tuple(p))
        for q in slot["points_Q"]:
            constrained.add(("b",) + tuple(q))

    n_total = len(pts_p) + len(pts_q)
    n_interior = len(interior)
    n_free_boundary = len(boundary) - len(constrained)
    n_constrained = len(constrained)

    params = 0
    detail = []
    for rec in records:
        if rec.get("vacuous"):
            continue
        shape = rec["power_structure"]
        p = shape["deg_k"] + 2
        params += p
        detail.append({
            "edge": rec["edge"],
            "normal": rec["normal"],
            "family_parameters": int(p),
        })
    # every vertex shared by two constrained edges is one coefficient that
    # both families must reproduce: matching it rescales one scalar, so each
    # shared vertex removes exactly one dimension
    shared_p, shared_q = shared_vertices(
        decomp["slots"][decomp["order"][0]],
        decomp["slots"][decomp["order"][1]])
    gluing = len(shared_p) + len(shared_q)
    params_after_gluing = params - gluing

    dim_scaled = n_interior + n_free_boundary + params_after_gluing
    dim_strict = dim_scaled - 2

    if n_interior + n_free_boundary + n_constrained != n_total:
        raise AssertionError("dimension bookkeeping does not close")

    return {
        "definition": (
            "dimension of the common solution set of all genuine edge "
            "leading-form conditions, in the full 445-dimensional coefficient "
            "space; not the dimension of det J = x^2"),
        "unknowns_total": int(n_total),
        "interior_coefficients_free": int(n_interior),
        "boundary_coefficients": int(len(boundary)),
        "boundary_coefficients_constrained": int(n_constrained),
        "boundary_coefficients_free": int(n_free_boundary),
        "edge_family_parameters": int(params),
        "gluing_conditions": int(gluing),
        "edge_family_parameters_after_gluing": int(params_after_gluing),
        "per_edge": detail,
        "codimension": int(n_total - dim_scaled),
        "dimension_scaled": int(dim_scaled),
        "dimension_strict": int(dim_strict),
        "breakdown": {
            "interior": int(n_interior),
            "free_boundary": int(n_free_boundary),
            "edge_parameters_after_gluing": int(params_after_gluing),
            "strict_chart_fixes": -2,
            "total_strict": int(dim_strict),
        },
        "consistent_with_recorded_per_edge_dimensions": bool(
            [r.get("dimension") for r in records] == [4, 1, 26, 42]),
        "convention": (
            "scaled counts all 445 coefficients; strict fixes "
            "a_(12,3) = b_(20,5) = 1, the chart used by cascade75.py"),
    }


# --------------------------------------------------------------------------
# stage 4: reject the 12 off-diagonal candidate conditions
# --------------------------------------------------------------------------

def joint_point(decomp, rng):
    """One coefficient assignment satisfying all four genuine edge conditions.

    The two diagonal edges are laid down first; the second is rescaled so
    that its values at a_(12,3) and b_(20,5) match the first.  The vacuous
    edges are then filled with random values, preserving already assigned
    shared vertices.
    """
    slots, order = decomp["slots"], decomp["order"]
    top, second = slots[order[0]], slots[order[1]]
    shape_top, shape_second = edge_family_shape(top), edge_family_shape(second)

    vals = {}
    k1 = [sp.Integer(rng.randint(-6, 6)) for _ in range(shape_second["deg_k"])]
    vals.update(family_point(second, k1, sp.Integer(1), sp.Integer(1)))

    k2 = [sp.Integer(rng.randint(-6, 6)) for _ in range(shape_top["deg_k"])]
    trial = family_point(top, k2, sp.Integer(1), sp.Integer(1))
    a_target = vals["a_12_3"]
    b_target = vals["b_20_5"]
    a_raw = trial["a_12_3"]
    b_raw = trial["b_20_5"]
    if a_raw == 0 or b_raw == 0:
        raise AssertionError("degenerate family point; retry with another seed")
    scaled = family_point(top, k2, a_target / a_raw, b_target / b_raw)
    if scaled["a_12_3"] != a_target or scaled["b_20_5"] != b_target:
        raise AssertionError("gluing rescale failed")
    vals.update(scaled)

    # vacuous edges: fill unassigned coefficients, keeping shared vertices
    for step in range(2, len(order)):
        slot = slots[order[step]]
        for p in slot["points_P"]:
            key = str(sym_a(*p))
            if key not in vals:
                vals[key] = sp.Integer(rng.randint(-5, 5))
        for q in slot["points_Q"]:
            key = str(sym_b(*q))
            if key not in vals:
                vals[key] = sp.Integer(rng.randint(-5, 5))

    vals = {k: sp.nsimplify(v) for k, v in vals.items()}
    sub = {sp.Symbol(k): v for k, v in vals.items()}

    checks = []
    for rec in range(len(order)):
        slot = slots[order[rec]]
        system = build_edge_system(slot["points_P"], slot["points_Q"])
        resid = [sp.expand(eq.subs(sub)) for eq in system["equations"]]
        checks.append({
            "edge": rec + 1,
            "normal": slot["normal"],
            "equations": len(resid),
            "nonzero_residuals": int(sum(1 for r in resid if r != 0)),
        })
    if any(c["nonzero_residuals"] for c in checks):
        raise AssertionError("joint point does not satisfy the edge conditions")
    return vals, checks


def rejected_candidate_conditions(decomp, vals):
    """The 12 off-diagonal pairs: not conditions, and false at the joint point.

    Two certificates are recorded per pair:
      * structural: the normal cones are disjoint (normals not parallel), so
        no weight exposes both faces and det J(P_ei,Q_ej) is not a leading
        form of det J;
      * numerical: at an exact rational point that satisfies all four genuine
        edge conditions, the mixed determinant is non-zero (for every pair
        where the test is defined).
    """
    sub = {sp.Symbol(k): v for k, v in vals.items()}
    edges_p, edges_q = decomp["edges_P"], decomp["edges_Q"]
    out = []
    for i in range(len(edges_p)):
        for j in range(len(edges_q)):
            if i == j:
                continue
            wp, wq = tuple(edges_p[i]["normal"]), tuple(edges_q[j]["normal"])
            structural = (_cross(wp, wq) == 0 and
                          wp[0] * wq[0] + wp[1] * wq[1] > 0)
            system = build_edge_system(edges_p[i]["points"], edges_q[j]["points"])
            det = sp.expand(sum(
                eq.subs(sub) * X ** m[0] * Y ** m[1]
                for eq, m in zip(system["equations"], system["monomials"])))
            out.append({
                "edge_P": i,
                "edge_Q": j,
                "normal_P": list(wp),
                "normal_Q": list(wq),
                "normal_cones_intersect": bool(structural),
                "mixed_det_at_joint_point_is_zero": bool(det == 0),
                "classification": ("would be an edge condition" if structural
                                   else "invalid: normal cones disjoint"),
            })
    if any(r["normal_cones_intersect"] for r in out):
        raise AssertionError("an off-diagonal pair has intersecting cones")
    false_at_witness = [r for r in out
                        if not r["mixed_det_at_joint_point_is_zero"]]
    return {
        "pairs": out,
        "count": len(out),
        "structural_certificate": (
            "for i != j the primitive outward normals are not parallel "
            "(cross product non-zero), so the normal cones -- which are the "
            "rays through those normals -- are disjoint; no weight w is "
            "maximal on both faces, hence det J(P_ei,Q_ej) is never "
            "LP_w(det J) and imposing it is unjustified"),
        "numerical_witness": {
            "point": "exact rational assignment in joint_point (satisfies all "
                     "four genuine edge conditions)",
            "mixed_determinants_nonzero": len(false_at_witness),
            "mixed_determinants_tested": len(out),
            "meaning": (
                "a cascade that descended these 12 pairs would report an "
                "emptiness that is an artefact: the conditions fail at a "
                "point the genuine cascade admits"),
        },
        "verdict": "REJECTED_NOT_CONDITIONS",
    }


# --------------------------------------------------------------------------
# stage 5: after a vacuous edge, the next filtration level
# --------------------------------------------------------------------------

def next_level_after_vacuous(decomp, rng):
    """Check the first sub-leading condition on the two vacuous edges.

    Horizontal edge (w=(0,-1)): P = P_0(x) + y P_1(x) + ..., Q = Q_0(x) +
    y Q_1(x) + ...; the weight-0 part of det J vanishes identically and the
    weight-(-1) part is P_0' Q_1 - P_1 Q_0'.  Choosing (P_1,Q_1) =
    (P_0',Q_0') satisfies it identically and stays inside the polygon.
    Vertical edge (w=(-1,0)) is the mirror image.  So neither vacuous edge
    hides an empty condition at the next order.

    The condition is *linear* in the next layer, so its exact solution
    dimension is computed by rank over Q; that dimension is reported as well.
    """
    x, y = X, Y
    out = []
    for slot in decomp["slots"]:
        if not slot["leading_form_is_identically_zero"]:
            continue
        if slot["edge_P"]["is_horizontal"]:
            deg_p0, deg_q0 = 9, 15
            deg_p1, deg_q1 = 10, 16
            var = x
            p0 = sum(sp.Integer(rng.randint(-6, 6)) * var ** i
                     for i in range(deg_p0 + 1))
            q0 = sum(sp.Integer(rng.randint(-6, 6)) * var ** i
                     for i in range(deg_q0 + 1))
            ps = sp.symbols("hp_0:%d" % (deg_p1 + 1))
            qs = sp.symbols("hq_0:%d" % (deg_q1 + 1))
            p1 = sum(ps[i] * var ** i for i in range(deg_p1 + 1))
            q1 = sum(qs[i] * var ** i for i in range(deg_q1 + 1))
            expr = sp.diff(p0, var) * q1 - p1 * sp.diff(q0, var)
            kind = "horizontal"
            layer = "y^1"
            condition = "P_0' Q_1 - P_1 Q_0' = 0"
            lift = (sp.diff(p0, var), sp.diff(q0, var))
        else:
            deg_p0, deg_q0 = 15, 25
            deg_p1, deg_q1 = 14, 24
            var = y
            p0 = sum(sp.Integer(rng.randint(-6, 6)) * var ** i
                     for i in range(deg_p0 + 1))
            q0 = sum(sp.Integer(rng.randint(-6, 6)) * var ** i
                     for i in range(deg_q0 + 1))
            ps = sp.symbols("vp_0:%d" % (deg_p1 + 1))
            qs = sp.symbols("vq_0:%d" % (deg_q1 + 1))
            p1 = sum(ps[i] * var ** i for i in range(deg_p1 + 1))
            q1 = sum(qs[i] * var ** i for i in range(deg_q1 + 1))
            expr = p1 * sp.diff(q0, var) - sp.diff(p0, var) * q1
            kind = "vertical"
            layer = "x^1"
            condition = "P_1 Q_0' - P_0' Q_1 = 0"
            lift = (sp.diff(p0, var), sp.diff(q0, var))

        poly = sp.Poly(sp.expand(expr), var)
        deg = int(poly.degree())
        gens = list(ps) + list(qs)
        rows = []
        for k in range(deg + 1):
            coeff = sp.expand(poly.coeff_monomial(var ** k))
            rows.append([sp.expand(sp.diff(coeff, g)) for g in gens])
        matrix = sp.Matrix(rows)
        rank = int(matrix.rank())
        n = len(gens)
        dim = n - rank

        # the explicit lift must be in the kernel, and fit the polygon
        lift_val = sp.expand(expr.subs(dict(zip(gens, lift))))
        out.append({
            "normal": slot["normal"],
            "kind": kind,
            "next_layer": layer,
            "condition": condition,
            "explicit_lift": "(P_1,Q_1) = (P_0',Q_0')",
            "condition_vanishes_identically": bool(lift_val == 0),
            "lift_inside_newton_polygon": bool(
                sp.degree(lift[0], var) <= deg_p1 and
                sp.degree(lift[1], var) <= deg_q1),
            "linear_system": {
                "unknowns": int(n),
                "equations": int(deg + 1),
                "rank_over_Q": int(rank),
                "solution_dimension": int(dim),
            },
            "conclusion": (
                "the next level is non-empty: a linear system with %d "
                "unknowns, %d equations, rank %d and dimension %d"
                % (n, deg + 1, rank, dim)),
        })
    return out


# --------------------------------------------------------------------------
# report
# --------------------------------------------------------------------------

def build_report(cap):
    t0 = time.perf_counter()
    rng = random.Random(SEED)
    decomp = build_decomposition()

    audit = premise_audit(decomp)
    log_stage("1_premise_audit", {
        "declared_in_task": audit["declared_in_task"],
        "geometric_edge_slots": audit["geometric_edge_slots"],
        "genuine_edge_conditions": audit["genuine_edge_conditions"],
        "invalid_off_diagonal_pairs": audit["invalid_off_diagonal_pairs"],
        "vertex_face_conditions_all_vacuous":
            audit["vertex_face_conditions_all_vacuous"],
    })

    records = descend_all(decomp, cap, rng)

    joint = joint_solution_set(decomp, records)
    log_stage("3_joint_dimension", {
        "dimension_scaled": joint["dimension_scaled"],
        "dimension_strict": joint["dimension_strict"],
        "codimension": joint["codimension"],
    })

    vals, point_checks = joint_point(decomp, rng)
    joint["explicit_joint_point_checks"] = point_checks
    log_stage("4_joint_point", {
        "edges_checked": len(point_checks),
        "max_nonzero_residuals": max(c["nonzero_residuals"]
                                     for c in point_checks),
    })

    rejected = rejected_candidate_conditions(decomp, vals)
    log_stage("5_rejected_pairs", {
        "count": rejected["count"],
        "mixed_determinants_nonzero":
            rejected["numerical_witness"]["mixed_determinants_nonzero"],
    })

    next_level = next_level_after_vacuous(decomp, rng)
    log_stage("6_next_level", {
        "edges": [e["normal"] for e in next_level],
        "all_identically_zero": all(
            e["condition_vanishes_identically"] for e in next_level),
    })

    first_empty = None
    for rec in records:
        if rec.get("empty"):
            first_empty = rec["edge"]
            break
    undecided = [r["edge"] for r in records
                 if r.get("dimension") is None]

    if first_empty is not None:
        verdict = "EMPTY_AT_EDGE %d" % first_empty
    elif undecided:
        verdict = "UNDECIDED_AT_EDGE %d" % undecided[0]
    else:
        verdict = "SURVIVES_ALL_EDGES"

    if verdict.startswith("EMPTY"):
        detail = ("edge %d has an empty leading-form condition" % first_empty)
    elif verdict.startswith("UNDECIDED"):
        detail = ("edge %d has no decided dimension" % undecided[0])
    else:
        detail = (
            "all %d genuine edge conditions are non-empty; the two diagonal "
            "edges carry the (3,5) power family (dimensions 4 and 1 in the "
            "strict chart) and the two horizontal/vertical edges are "
            "vacuous. The %d off-diagonal pairs in the 4x4 grid are not "
            "leading-form conditions and are rejected, not descended. The "
            "edge conditions cut the 445-dimensional coefficient space down "
            "to dimension %d (strict chart %d), so the cascade is a filter "
            "that cannot decide this corner: it is non-empty by structure, "
            "not by accident."
            % (len(records), rejected["count"], joint["dimension_scaled"],
               joint["dimension_strict"]))

    if verdict == "SURVIVES_ALL_EDGES":
        verdict_detail = detail
    elif verdict.startswith("EMPTY"):
        verdict_detail = detail
    else:
        verdict_detail = detail

    checks = {
        "four_genuine_edge_conditions": (
            audit["genuine_edge_conditions"] == 4),
        "declared_16_is_not_the_geometric_count": (
            audit["declared_in_task"] == 16 and
            audit["geometric_edge_slots"] == 4),
        "all_four_vertex_pairs_vacuous":
            audit["vertex_face_conditions_all_vacuous"],
        "four_edges_descended": len(records) == 4,
        "no_edge_empty": all(rec.get("empty") is False for rec in records),
        "every_non_vacuous_edge_has_an_explicit_solution": all(
            rec.get("vacuous") or
            rec["family_certificate"]["all_equations_vanish"]
            for rec in records),
        "every_non_vacuous_edge_family_has_identically_zero_det_J": all(
            rec.get("vacuous") or
            rec["family_certificate"]["det_J_identically_zero_for_symbolic_k"]
            for rec in records),
        "joint_point_satisfies_all_edges": all(
            c["nonzero_residuals"] == 0 for c in point_checks),
        "joint_dimension_bookkeeping_closes": (
            joint["interior_coefficients_free"]
            + joint["boundary_coefficients_free"]
            + joint["boundary_coefficients_constrained"]
            == joint["unknowns_total"]),
        "per_edge_dimensions_match_cascade75": (
            [r.get("dimension") for r in records] == [4, 1, 26, 42]),
        "all_twelve_off_diagonal_cones_disjoint": all(
            not p["normal_cones_intersect"] for p in rejected["pairs"]),
        "mixed_conditions_fail_at_joint_point": (
            rejected["numerical_witness"]["mixed_determinants_nonzero"] > 0),
        "vacuous_edges_non_empty_at_next_level": all(
            e["condition_vanishes_identically"] for e in next_level),
        "next_level_linear_systems_have_positive_dimension": all(
            e["linear_system"]["solution_dimension"] > 0 for e in next_level),
        "lead_verdict_is_one_of_the_three_allowed_forms": verdict in (
            "SURVIVES_ALL_EDGES",
            "EMPTY_AT_EDGE %d" % (first_empty or 0),
            "UNDECIDED_AT_EDGE %d" % (undecided[0] if undecided else 0)),
    }

    report = {
        "problem": "jacobian-75-125-cascade-full",
        "case": "F_2(3,5), j=1",
        "method": "Newton-polygon edge filtration / leading-form cascade, "
                  "all genuine edge conditions",
        "task_premise": {
            "declared_edge_conditions": DECLARED_CONDITIONS,
            "geometric_edge_conditions": audit["geometric_edge_slots"],
            "finding": audit["finding"],
        },
        "premise_audit": audit,
        "edges": records,
        "edge_dimensions": {
            "per_edge_strict": [r.get("dimension") for r in records],
            "per_edge_scaled": [r.get("dimension_scaled") for r in records],
            "order_normals": [r["normal"] for r in records],
            "convention": "strict fixes a_(12,3) = b_(20,5); vacuous edges "
                          "report all their coefficients as free",
        },
        "joint_solution_set": joint,
        "rejected_candidate_conditions": rejected,
        "next_level_after_vacuous": next_level,
        "first_empty_edge": first_empty,
        "edges_done": len(records),
        "verdict": verdict,
        "verdict_detail": verdict_detail,
        "verdict_reading": {
            "SURVIVES_ALL_EDGES": "no genuine edge condition is empty; the "
                                  "final dimension is in joint_solution_set",
            "EMPTY_AT_EDGE": "the named edge has no solution compatible with "
                             "the coefficients carried from earlier edges",
            "UNDECIDED_AT_EDGE": "the named edge exceeded its cap; no "
                                 "emptiness claim is made",
        }[verdict.split(" ")[0]],
        "checks": {
            "all_pass": all(checks.values()),
            "checks": {k: bool(v) for k, v in checks.items()},
        },
        "conventions": {
            "leading_part": "LP_w(det J) = det J(P_w, Q_w)",
            "strict_chart": "a_(12,3) = b_(20,5) = 1",
            "dimension": "affine dimension of the edge-condition variety in "
                         "the stated chart; not the dimension of det J = x^2",
            "edge_order": "boundary cycle from the top edge, clockwise; "
                          "equals decreasing filtration degree",
        },
        "timings": {},
        "elapsed_seconds": None,
    }
    report["timings"]["total_seconds"] = round(time.perf_counter() - t0, 3)
    report["elapsed_seconds"] = report["timings"]["total_seconds"]
    return report


def print_report(rep):
    print("== full (75,125) edge cascade ==")
    print("premise audit : %d declared, %d genuine edge conditions, "
          "%d invalid off-diagonal pairs"
          % (rep["premise_audit"]["declared_in_task"],
             rep["premise_audit"]["genuine_edge_conditions"],
             rep["premise_audit"]["invalid_off_diagonal_pairs"]))
    print("                 vertex-face conditions: %d, all vacuous: %s"
          % (len(rep["premise_audit"]["vertex_face_conditions"]),
             rep["premise_audit"]["vertex_face_conditions_all_vacuous"]))
    print()
    for rec in rep["edges"]:
        tag = "VACUOUS" if rec.get("vacuous") else "dim %s (scaled %s)" % (
            rec.get("dimension"), rec.get("dimension_scaled"))
        print("edge %d  w=%-8s unknowns %-3d equations %-2d  %s"
              % (rec["edge"], rec["normal"], rec["unknowns"],
                 rec["equations"], tag))
        print("        %s" % (rec.get("non_empty_certificate")
                              or rec.get("reason")))
        gb = rec.get("groebner", {})
        print("        groebner: %s%s"
              % (gb.get("status"), " (cap %s s)"
                 % gb.get("cap_seconds") if gb.get("cap_seconds") else ""))
    print()
    j = rep["joint_solution_set"]
    print("joint edge variety : dimension %d scaled, %d strict "
          "(codimension %d of %d)"
          % (j["dimension_scaled"], j["dimension_strict"], j["codimension"],
             j["unknowns_total"]))
    print("  breakdown        : interior %d + free boundary %d + edge "
          "parameters %d - strict fixes 2"
          % (j["breakdown"]["interior"], j["breakdown"]["free_boundary"],
             j["breakdown"]["edge_parameters_after_gluing"]))
    print()
    rj = rep["rejected_candidate_conditions"]
    print("edges 5..16        : %d off-diagonal pairs REJECTED "
          "(normal cones disjoint)" % rj["count"])
    print("  witness          : %d of %d mixed determinants are non-zero at "
          "an exact point satisfying all four genuine edges"
          % (rj["numerical_witness"]["mixed_determinants_nonzero"],
             rj["numerical_witness"]["mixed_determinants_tested"]))
    print()
    for e in rep["next_level_after_vacuous"]:
        ls = e["linear_system"]
        print("beyond w=%-8s (vacuous): next level %s is linear, "
              "%d unknowns / %d equations, rank %d, dimension %d"
              % (e["normal"], e["next_layer"], ls["unknowns"],
                 ls["equations"], ls["rank_over_Q"],
                 ls["solution_dimension"]))
    print()
    print("first empty edge   : %s" % rep["first_empty_edge"])
    print("VERDICT: %s" % rep["verdict"])
    print("checks: %s" % ("all pass" if rep["checks"]["all_pass"]
                          else "FAILURES: %s" % [
                              k for k, v in rep["checks"]["checks"].items()
                              if not v]))
    print("elapsed: %.3f s" % rep["elapsed_seconds"])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--no-write", action="store_true")
    ap.add_argument("--cap-edge", type=float, default=EDGE_CAP_S)
    args = ap.parse_args(argv)

    t0 = time.perf_counter()
    rep = build_report(cap=args.cap_edge)
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
