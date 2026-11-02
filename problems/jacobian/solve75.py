"""Settle the (75,125) vertices-only coefficient system.

CORNER DATA
-----------
Verbatim from ``research/jacobian-degree-pairs-ge125.md`` (itself from the
wstrinz repo, ``d2/PASSPORT_75_125_REPAIR.md:29-33``):

    Delta' = {(0,0),(3,0),(4,1),(0,5)}
    N(P) = 3*Delta' = {(0,0),(9,0),(12,3),(0,15)}
    N(Q) = 5*Delta' = {(0,0),(15,0),(20,5),(0,25)}
    [P,Q] = x^2      reduced degrees (15,25)

This script takes the *vertices-only* reading: the support of P is exactly the
four corners of N(P), and the support of Q is exactly the four corners of N(Q).
That is 4 + 4 = 8 unknown coefficients and, after subtracting x^2, 7 coefficient
equations (det J has 6 monomials and x^2 contributes one more).

WHAT IS DECIDED
---------------
The system is INCONSISTENT over Q, hence over C.  The reduced grevlex Groebner
basis is ``[1]``.  There is also a support-level proof that does not need a
basis: every monomial of det J has x-degree at least 8, so the coefficient of
x^2 in det J is identically zero, while det J - x^2 asks it to be 1.  The
(2,0) equation is literally ``-1 = 0``.

This is a statement about the minimal-support specialization only.  It does
NOT exclude the full-lattice system of ``config75.py``: there x^2 is present in
det J, and the coefficient system remains live.

METHOD
------
* Build P, Q with ``config2_system.build_polynomial`` (the same helper used by
  ``config2_system.py``, ``config75.py`` and ``solve_vertices.py``).
* Cross-check every coefficient of det J against the independent bilinear
  accumulator in ``config75.det_coefficients``.
* Run sympy's grevlex Groebner basis over Q in a killable child process with a
  hard 900 s cap, plus a 60 s deadline-bounded Buchberger run that records the
  intermediate basis size (needed only if the main run times out).
* If the basis is not ``[1]``: extract solutions, verify det J = x^2 exactly,
  reject the two triangular shear forms, hunt collisions with ``jac.py``, and
  submit the first exact collision to ``tools/accept.py``.

Every stage appends one line to ``runs/stress_jac.log``.  The machine-readable
record is ``runs/jacobian_75_125_solve.json``.

Run:
    python problems/jacobian/solve75.py
"""

from __future__ import annotations

import argparse
import json
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

from config2_system import build_polynomial  # noqa: E402
from config75 import (  # noqa: E402
    DELTA_PRIME,
    P_CORNERS,
    Q_CORNERS,
    build_variant,
    det_coefficients,
)
from solve_vertices import buchberger_progress  # noqa: E402

import jac  # noqa: E402

DEFAULT_OUT = ROOT / "runs" / "jacobian_75_125_solve.json"
LOG_PATH = ROOT / "runs" / "stress_jac.log"
CLAIM_PATH = HERE / "claims" / "solve75_candidate.json"

GROEBNER_TIMEOUT_S = 900.0      # per instruction: cap at 900 s
TOTAL_BUDGET_S = 1500.0         # acceptance harness allows 1500 s
PROGRESS_TIMEOUT_S = 60.0       # independent Buchberger progress run
CERTIFICATE_MAX_S = 120.0       # only recompute in-parent for a cheap basis

SYMBOL_NAMES = ["a0", "a1", "a2", "a3", "b0", "b1", "b2", "b3"]
TARGET = (2, 0)


# --------------------------------------------------------------------------
# logging
# --------------------------------------------------------------------------

def log_stage(stage, payload):
    """Append exactly one line for one stage.  Never used for bulk output."""
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    line = "%s stage=%s %s\n" % (
        stamp, stage, json.dumps(payload, sort_keys=True, default=str))
    with LOG_PATH.open("a", encoding="utf-8") as fh:
        fh.write(line)


# --------------------------------------------------------------------------
# the system
# --------------------------------------------------------------------------

def build_system():
    """det J(P,Q) - x^2 on the four corners of each Newton polygon.

    Construction is delegated to ``config2_system.build_polynomial`` so the
    supports, coefficient names and order cannot drift from the sized systems.
    """
    x, y = sp.symbols("x y")
    P, cp = build_polynomial(list(P_CORNERS), "a", x, y)
    Q, cq = build_polynomial(list(Q_CORNERS), "b", x, y)
    unknowns = list(cp) + list(cq)

    names = [s.name for s in unknowns]
    if names != SYMBOL_NAMES:
        raise RuntimeError("unexpected coefficient order: %s" % names)

    det = sp.expand(P.diff(x) * Q.diff(y) - P.diff(y) * Q.diff(x))
    det_poly = sp.Poly(det, x, y)
    det_monomials = sorted(
        (int(m[0]), int(m[1])) for m in det_poly.monoms())

    system_poly = sp.Poly(det - x ** 2, x, y)
    monomials = sorted(
        ((int(m[0]), int(m[1])) for m in system_poly.monoms()),
        key=lambda m: (m[0] + m[1], m[0], m[1]))
    equations = [sp.expand(system_poly.coeff_monomial(m)) for m in monomials]

    reconstructed = sp.Add(*[
        eq * x ** p * y ** q for eq, (p, q) in zip(equations, monomials)])
    if sp.expand(reconstructed - (det - x ** 2)) != 0:
        raise RuntimeError("coefficient reconstruction failed")

    if TARGET not in monomials:
        raise RuntimeError("x^2 is missing from det J - x^2; build is wrong")

    target_equation = equations[monomials.index(TARGET)]

    return {
        "x": x,
        "y": y,
        "P": P,
        "Q": Q,
        "det": det,
        "unknowns": unknowns,
        "equations": equations,
        "monomials": monomials,
        "det_monomials": det_monomials,
        "target_equation": target_equation,
        "supports": {
            "P": [list(p) for p in P_CORNERS],
            "Q": [list(q) for q in Q_CORNERS],
        },
        "degrees": [
            max(a + b for a, b in P_CORNERS),
            max(a + b for a, b in Q_CORNERS),
        ],
    }


def _sympy_coefficient_map(data):
    """{(p,q): {(i,j): int}} for det J, read from the expanded sympy object."""
    det_poly = sp.Poly(data["det"], data["x"], data["y"])
    unknowns = data["unknowns"]
    out = {}
    for m in data["det_monomials"]:
        coeff = sp.expand(det_poly.coeff_monomial(m))
        poly = sp.Poly(coeff, *unknowns)
        terms = {}
        for monom, value in zip(poly.monoms(), poly.coeffs()):
            nz = [i for i, e in enumerate(monom) if e]
            if len(nz) == 2 and all(monom[i] == 1 for i in nz):
                ia, ib = nz
                if ia < 4 <= ib:
                    terms[(ia, ib - 4)] = int(value)
                else:
                    terms[("unexpected", tuple(nz))] = int(value)
            elif not nz:
                terms[("constant", ())] = int(value)
            else:
                terms[("unexpected", tuple(monom))] = int(value)
        out[m] = terms
    return out


def support_obstruction(data):
    """The x-degree argument, independent of any Groebner computation."""
    det_monomials = data["det_monomials"]
    min_x = min(p for p, _ in det_monomials)
    target_present = TARGET in det_monomials
    decisive = (not target_present) and min_x > TARGET[0]
    return {
        "det_monomials": [list(m) for m in det_monomials],
        "min_x_degree_of_det_J": int(min_x),
        "target": "x^2",
        "target_x_degree": int(TARGET[0]),
        "target_monomial_present_in_det_J": bool(target_present),
        "target_coefficient_equation": sp.sstr(data["target_equation"]),
        "decisive": bool(decisive),
        "argument": (
            "Every monomial of det J has x-degree at least %d, so the "
            "coefficient of x^2 in det J is identically zero.  The identity "
            "det J = x^2 requires that coefficient to be 1, so the (2,0) "
            "coefficient equation is -1 = 0.  No Groebner basis is needed for "
            "this contradiction." % min_x
        ),
        "scope": (
            "This is a support obstruction of the vertices-only "
            "specialization.  It says nothing about the full-lattice system, "
            "where x^2 is present in det J."
        ),
    }


def stage1(data, seconds):
    """Build-time record plus two independent exact cross-checks."""
    bilinear, pair_products = det_coefficients(
        list(P_CORNERS), list(Q_CORNERS))
    sympy_map = _sympy_coefficient_map(data)
    bilinear_map = {
        m: {tuple(k): int(v) for k, v in terms.items()}
        for m, terms in bilinear.items()
    }
    variant = build_variant("vertices_only", list(P_CORNERS), list(Q_CORNERS))

    min_x = min(p for p, _ in data["det_monomials"])
    target_present = TARGET in data["det_monomials"]
    checks = {
        "unknowns_is_8": len(data["unknowns"]) == 8,
        "equations_is_7": len(data["equations"]) == 7,
        "reconstruction_ok": True,
        "det_monomial_count_matches_bilinear_accumulator": (
            len(data["det_monomials"]) == len(bilinear)),
        "det_coefficient_map_matches_bilinear_accumulator": (
            sympy_map == bilinear_map),
        "config75_agrees_on_unknowns": variant["unknowns"] == 8,
        "config75_agrees_on_equations": (
            variant["target_x2"]["equations"] == len(data["equations"])),
        "config75_agrees_target_absent": (
            variant["target_x2"]["target_monomial_present_in_det"] is False),
        "config75_agrees_target_equation": (
            variant["target_x2"]["target_monomial_equation"]["coefficient"]
            == "-1"),
        "target_equation_is_minus_one": data["target_equation"] == -1,
        "pair_products": pair_products,
    }
    critical = [
        "unknowns_is_8",
        "equations_is_7",
        "det_monomial_count_matches_bilinear_accumulator",
        "det_coefficient_map_matches_bilinear_accumulator",
        "config75_agrees_on_unknowns",
        "config75_agrees_on_equations",
        "config75_agrees_target_absent",
        "config75_agrees_target_equation",
    ]
    bad = [name for name in critical if not checks[name]]
    if bad:
        raise RuntimeError("stage-1 cross-checks failed: %s" % bad)

    return {
        "status": "ok",
        "support_convention": "vertices_only (four corners per component)",
        "Delta_prime": [list(p) for p in DELTA_PRIME],
        "support_P": data["supports"]["P"],
        "support_Q": data["supports"]["Q"],
        "unknowns": len(data["unknowns"]),
        "equations": len(data["equations"]),
        "equation_monomials": [list(m) for m in data["monomials"]],
        "det_monomials": [list(m) for m in data["det_monomials"]],
        "det_J": sp.sstr(data["det"]),
        "degrees": data["degrees"],
        "target": "det J(P,Q) - x^2 = 0",
        "min_x_degree_of_det_J": int(min_x),
        "target_monomial_present_in_det_J": bool(target_present),
        "target_coefficient_equation": {
            "monomial": list(TARGET),
            "polynomial": sp.sstr(data["target_equation"]),
        },
        "cross_checks": {k: (v if isinstance(v, (bool, int)) else str(v))
                         for k, v in checks.items()},
        "all_critical_checks_pass": not bad,
        "seconds": round(seconds, 3),
    }


# --------------------------------------------------------------------------
# Groebner basis in a killable child process
# --------------------------------------------------------------------------

def _gb_child():
    data = build_system()
    t0 = time.perf_counter()
    G = sp.groebner(
        data["equations"], *data["unknowns"], order="grevlex")
    sys.stdout.write(json.dumps({
        "basis": [sp.sstr(e) for e in G.exprs],
        "basis_size": len(G.exprs),
        "order": "grevlex",
        "seconds": time.perf_counter() - t0,
    }))
    sys.stdout.flush()
    return 0


def groebner_capped(timeout):
    """Run sympy's grevlex basis over Q, hard-killed after `timeout` seconds."""
    t0 = time.perf_counter()
    try:
        proc = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "--gb-child"],
            capture_output=True, text=True, encoding="utf-8",
            timeout=timeout, cwd=str(HERE))
    except subprocess.TimeoutExpired:
        return {"timed_out": True, "basis": None, "basis_size": None,
                "seconds": round(time.perf_counter() - t0, 3)}
    if proc.returncode != 0:
        raise RuntimeError(
            "groebner child failed:\n" + (proc.stderr or "")[-4000:])
    text = (proc.stdout or "").strip()
    payload = json.loads(text.splitlines()[-1] if text else "{}")
    payload["timed_out"] = False
    payload["seconds"] = round(time.perf_counter() - t0, 3)
    return payload


def stage2(data, gb, progress, obstruction, gb_timeout):
    """Record the decisive basis and a cheap independent certificate."""
    if gb["timed_out"]:
        return {
            "status": "timeout",
            "order": "grevlex",
            "timeout_seconds": round(gb_timeout, 3),
            "timed_out": True,
            "basis": None,
            "basis_size": None,
            "intermediate_basis_size": progress.get("basis_size"),
            "intermediate_basis_trace": progress.get("trace"),
            "intermediate_pairs_processed": progress.get("pairs_processed"),
            "intermediate_pairs_left": progress.get("pairs_left"),
            "seconds": gb["seconds"],
            "progress_cross_check": _progress_record(progress),
            "support_obstruction": obstruction,
            "inconsistent": None,
            "note": (
                "The Groebner run exceeded the cap; the intermediate basis "
                "size is from the deadline-bounded Buchberger cross-check."
            ),
        }

    basis_strings = list(gb["basis"])
    basis_exprs = [
        sp.sympify(s, locals={s.name: s for s in data["unknowns"]})
        for s in basis_strings
    ]
    inconsistent = basis_exprs == [sp.Integer(1)]

    certificate = None
    if gb["seconds"] <= CERTIFICATE_MAX_S:
        t0 = time.perf_counter()
        G = sp.groebner(
            data["equations"], *data["unknowns"], order="grevlex")
        certificate = {
            "each_equation_reduces_to_zero": all(
                G.reduce(e)[1] == 0 for e in data["equations"]),
            "recomputed_basis": [sp.sstr(e) for e in G.exprs],
            "recomputed_basis_size": len(G.exprs),
            "matches_returned_basis": (
                [sp.sstr(e) for e in G.exprs] == basis_strings),
            "seconds": round(time.perf_counter() - t0, 3),
        }
        if not certificate["each_equation_reduces_to_zero"]:
            raise RuntimeError("certificate failed: an equation does not reduce")
        if not certificate["matches_returned_basis"]:
            raise RuntimeError("certificate failed: recomputed basis differs")

    return {
        "status": "ok",
        "order": "grevlex",
        "timeout_seconds": round(gb_timeout, 3),
        "timed_out": False,
        "basis": basis_strings,
        "basis_size": len(basis_strings),
        "inconsistent": bool(inconsistent),
        "seconds": gb["seconds"],
        "progress_cross_check": _progress_record(progress),
        "certificate": certificate,
        "support_obstruction": obstruction,
    }


def _progress_record(progress):
    out = {
        "ok": progress.get("ok"),
        "timed_out": progress.get("timed_out"),
        "basis_size": progress.get("basis_size"),
        "basis": progress.get("basis"),
        "pairs_processed": progress.get("pairs_processed"),
        "pairs_left": progress.get("pairs_left"),
        "seconds": progress.get("seconds"),
    }
    if progress.get("error"):
        out["error"] = progress["error"]
    return out


# --------------------------------------------------------------------------
# solution extraction and exact verification (only reached if basis != [1])
# --------------------------------------------------------------------------

def _monomials_of(expr, x, y):
    expr = sp.expand(expr)
    if expr == 0:
        return []
    return sorted(
        (int(a), int(b)) for a, b in sp.Poly(expr, x, y).monoms())


def triangular_reason(P, Q, x, y):
    """Mirror tools/accept.py's two obvious-shear tests."""
    if sp.simplify(P - x) == 0:
        return "P == x: the (x, y + h(x)) shear form"
    if sp.simplify(Q - y) == 0:
        degree_in_x = sp.Poly(sp.expand(P - x), x).degree()
        if degree_in_x <= 0:
            return "Q == y and P has no x-dependence: the (x + h(y), y) form"
    return None


def _specialize(expressions, params):
    """Pick small positive integers for free parameters; keep denominators safe."""
    params = sorted(params, key=lambda s: s.name)
    for base in (1, 2, 3, 5, 7, 11):
        values = {s: sp.Integer(base + i) for i, s in enumerate(params)}
        out = [sp.expand(e.subs(values)) for e in expressions]
        if all(not e.has(*params) for e in out):
            return values, out
    return None, None


def _critical_collision(P, Q, x, y):
    """Exact collision on the critical line x = 0, if both components are flat."""
    P0 = sp.expand(P.subs(x, 0))
    Q0 = sp.expand(Q.subs(x, 0))
    if P0.has(y) or Q0.has(y):
        return None
    p = (sp.Integer(0), sp.Integer(0))
    q = (sp.Integer(0), sp.Integer(1))
    img_p = (sp.simplify(P.subs({x: p[0], y: p[1]})),
             sp.simplify(Q.subs({x: p[0], y: p[1]})))
    img_q = (sp.simplify(P.subs({x: q[0], y: q[1]})),
             sp.simplify(Q.subs({x: q[0], y: q[1]})))
    if all(sp.simplify(a - b) == 0 for a, b in zip(img_p, img_q)):
        return {
            "points": [[sp.sstr(c) for c in p], [sp.sstr(c) for c in q]],
            "distinct": p != q,
            "images": [[sp.sstr(c) for c in img_p],
                       [sp.sstr(c) for c in img_q]],
            "same_image": True,
            "locus": "x = 0, where det J = x^2 vanishes",
        }
    return None


def analyze_solutions(data):
    gens = list(data["unknowns"])
    eqs = data["equations"]
    x, y = data["x"], data["y"]
    try:
        branches = sp.solve(eqs, gens, dict=True)
    except Exception as exc:  # diagnostics must not hide behind an exception
        return {"status": "error",
                "error": "%s: %s" % (type(exc).__name__, exc),
                "number_of_branches": 0, "branches": []}

    records = []
    for index, sub in enumerate(branches):
        residuals = [sp.simplify(sp.expand(e.subs(sub))) for e in eqs]
        P_sub = sp.expand(data["P"].subs(sub))
        Q_sub = sp.expand(data["Q"].subs(sub))
        det_sub = sp.expand(data["det"].subs(sub))
        det_ok = sp.expand(det_sub - x ** 2) == 0
        tri = triangular_reason(P_sub, Q_sub, x, y)

        p_monomials = _monomials_of(P_sub, x, y)
        q_monomials = _monomials_of(Q_sub, x, y)
        uses_all = (
            p_monomials == sorted(tuple(p) for p in P_CORNERS)
            and q_monomials == sorted(tuple(q) for q in Q_CORNERS))

        record = {
            "index": index,
            "parametrization": {
                s.name: (sp.sstr(sub[s]) if s in sub else "free")
                for s in gens
            },
            "free_parameters": [s.name for s in gens if s not in sub],
            "residuals_zero": all(r == 0 for r in residuals),
            "P": sp.sstr(P_sub),
            "Q": sp.sstr(Q_sub),
            "det_J": sp.sstr(det_sub),
            "det_J_minus_x2": sp.sstr(sp.expand(det_sub - x ** 2)),
            "det_J_equals_x2": bool(det_ok),
            "triangular": tri is not None,
            "triangular_reason": tri,
            "support_used": {
                "P": [list(m) for m in p_monomials],
                "Q": [list(m) for m in q_monomials],
            },
            "uses_all_four_corners": bool(uses_all),
            "non_triangular": bool(
                det_ok and tri is None and all(r == 0 for r in residuals)),
            "representative": None,
        }

        if record["non_triangular"]:
            params = (P_sub.free_symbols | Q_sub.free_symbols) - {x, y}
            values, specialized = _specialize(
                [P_sub, Q_sub, det_sub], params)
            if specialized is not None:
                P_num, Q_num = specialized[0], specialized[1]
                det_num = sp.expand(
                    P_num.diff(x) * Q_num.diff(y)
                    - P_num.diff(y) * Q_num.diff(x))
                exact = _critical_collision(P_num, Q_num, x, y)
                hits = jac.find_collision(
                    [P_num, Q_num], [x, y], samples=8, starts=6, seed=0)
                record["representative"] = {
                    "parameter_values": {
                        s.name: sp.sstr(v) for s, v in values.items()},
                    "P": sp.sstr(P_num),
                    "Q": sp.sstr(Q_num),
                    "det_J": sp.sstr(det_num),
                    "det_J_equals_x2": bool(
                        sp.expand(det_num - x ** 2) == 0),
                    "exact_collision": exact,
                    "numerical_hits": hits,
                    "numerical_hit_count": len(hits),
                    "note": (
                        "numerical hits are leads, not proofs; only an exact "
                        "collision is carried to accept.py"
                    ),
                }
            else:
                record["representative"] = {
                    "error": "could not specialize free parameters"}
        records.append(record)

    return {
        "status": "ok",
        "number_of_branches": len(branches),
        "branches": records,
    }


def stage3(solutions):
    if solutions is None:
        return {
            "status": "skipped",
            "reason": "Groebner basis == [1]: no solutions over C, so there "
                      "is nothing to extract",
        }
    return solutions


def stage4(solutions):
    if not solutions or solutions.get("status") != "ok":
        return {"status": "skipped", "reason": "no solution set to search"}
    candidates = [b for b in solutions["branches"] if b.get("non_triangular")]
    if not candidates:
        return {
            "status": "skipped",
            "reason": "no non-triangular solution",
            "search_engine": "jac.find_collision (mpmath/nsolve, numerical)",
            "collisions": [],
        }
    collisions = []
    for branch in candidates:
        rep = branch.get("representative") or {}
        collisions.append({
            "solution_index": branch["index"],
            "P": rep.get("P"),
            "Q": rep.get("Q"),
            "det_J": rep.get("det_J"),
            "exact": rep.get("exact_collision"),
            "numerical_hits": rep.get("numerical_hits", []),
            "numerical_hit_count": rep.get("numerical_hit_count", 0),
        })
    return {
        "status": "ok",
        "search_engine": "jac.find_collision (mpmath/nsolve, numerical)",
        "candidates": len(candidates),
        "collisions": collisions,
        "exact_collision_count": sum(1 for c in collisions if c["exact"]),
    }


def stage5(collision_stage):
    exact = [c for c in collision_stage.get("collisions", []) if c["exact"]]
    if not exact:
        return {
            "status": "skipped",
            "reason": "no exact collision; numerical leads cannot be a claim",
            "ran": False,
            "accepted": None,
            "exit_code": None,
            "stdout": None,
        }
    collision = exact[0]
    CLAIM_PATH.parent.mkdir(parents=True, exist_ok=True)
    claim = {
        "note": (
            "Generated by problems/jacobian/solve75.py for the (75,125) "
            "vertices-only system.  det J = x^2, which is not a non-zero "
            "constant; this is expected to fail the Jacobian gate J1."
        ),
        "vars": ["x", "y"],
        "F": [collision["P"], collision["Q"]],
        "points": collision["exact"]["points"],
        "det_J": collision["det_J"],
        "target": "det J = x^2",
    }
    CLAIM_PATH.write_text(json.dumps(claim, indent=2), encoding="utf-8")
    cmd = [sys.executable, str(ROOT / "tools" / "accept.py"),
           "--problem", "jac", "--claim", str(CLAIM_PATH)]
    proc = subprocess.run(
        cmd, capture_output=True, text=True, encoding="utf-8", cwd=str(ROOT))
    return {
        "status": "ran",
        "ran": True,
        "claim_path": str(CLAIM_PATH.relative_to(ROOT)),
        "command": cmd,
        "exit_code": proc.returncode,
        "accepted": proc.returncode == 0,
        "stdout": proc.stdout,
        "stderr": proc.stderr,
        "verdict_verbatim": proc.stdout.strip(),
    }


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------

def run(out_path: Path, timeout: float, write: bool = True):
    t_start = time.perf_counter()
    stages = {}

    # ---- stage 1: build ----------------------------------------------------
    t0 = time.perf_counter()
    data = build_system()
    build_seconds = time.perf_counter() - t0
    obstruction = support_obstruction(data)
    stages["1_build"] = stage1(data, build_seconds)
    log_stage(1, {
        "status": "ok",
        "unknowns": stages["1_build"]["unknowns"],
        "equations": stages["1_build"]["equations"],
        "det_monomials": len(data["det_monomials"]),
        "min_x_degree": obstruction["min_x_degree_of_det_J"],
        "target_present": obstruction["target_monomial_present_in_det_J"],
        "seconds": stages["1_build"]["seconds"],
    })

    # ---- stage 2: Groebner basis ------------------------------------------
    elapsed = time.perf_counter() - t_start
    progress_budget = max(
        1.0, min(PROGRESS_TIMEOUT_S, TOTAL_BUDGET_S - elapsed - 30.0))
    progress = buchberger_progress(
        data["equations"], data["unknowns"],
        deadline=time.monotonic() + progress_budget)

    elapsed = time.perf_counter() - t_start
    gb_timeout = max(1.0, min(timeout, TOTAL_BUDGET_S - elapsed - 20.0))
    gb = groebner_capped(gb_timeout)
    stages["2_groebner"] = stage2(
        data, gb, progress, obstruction, gb_timeout)
    log_stage(2, {
        "status": stages["2_groebner"]["status"],
        "timed_out": stages["2_groebner"]["timed_out"],
        "basis_size": stages["2_groebner"]["basis_size"],
        "intermediate_basis_size": stages["2_groebner"].get(
            "intermediate_basis_size"),
        "inconsistent": stages["2_groebner"]["inconsistent"],
        "seconds": stages["2_groebner"]["seconds"],
    })

    # ---- stages 3-5: conditional ------------------------------------------
    gb_completed = stages["2_groebner"]["status"] == "ok"
    basis_is_one = (
        gb_completed
        and stages["2_groebner"]["basis_size"] == 1
        and stages["2_groebner"]["inconsistent"])

    if not gb_completed:
        solutions_record = None
        stages["3_solutions"] = {
            "status": "skipped",
            "reason": "Groebner basis timed out; no basis to solve",
        }
        log_stage(3, {
            "status": "skipped",
            "reason": "groebner timeout",
            "solutions": 0,
        })
    elif basis_is_one:
        solutions_record = None
        stages["3_solutions"] = stage3(None)
        log_stage(3, {
            "status": "skipped",
            "reason": "basis == [1]",
            "solutions": 0,
        })
    else:
        solutions_record = analyze_solutions(data)
        stages["3_solutions"] = stage3(solutions_record)
        log_stage(3, {
            "status": stages["3_solutions"]["status"],
            "branches": stages["3_solutions"].get("number_of_branches"),
            "non_triangular": sum(
                1 for b in stages["3_solutions"].get("branches", [])
                if b.get("non_triangular")),
        })

    if solutions_record is None:
        collision_record = stage4(None)
        stages["4_collisions"] = collision_record
        log_stage(4, {"status": "skipped",
                      "reason": stages["3_solutions"]["status"],
                      "exact_collisions": 0})
        stages["5_accept"] = stage5(collision_record)
        log_stage(5, {"status": "skipped", "reason": "no candidate",
                      "accepted": None, "exit_code": None})
    else:
        collision_record = stage4(solutions_record)
        stages["4_collisions"] = collision_record
        log_stage(4, {
            "status": collision_record["status"],
            "candidates": collision_record.get("candidates", 0),
            "exact_collisions": collision_record.get(
                "exact_collision_count", 0),
            "numerical_hits": sum(
                c.get("numerical_hit_count", 0)
                for c in collision_record.get("collisions", [])),
        })
        stages["5_accept"] = stage5(collision_record)
        log_stage(5, {
            "status": stages["5_accept"]["status"],
            "accepted": stages["5_accept"].get("accepted"),
            "exit_code": stages["5_accept"].get("exit_code"),
        })

    # ---- verdict -----------------------------------------------------------
    timed_out = stages["2_groebner"]["timed_out"]
    if basis_is_one:
        verdict = "INCONSISTENT"
        inconsistency_method = "reduced grevlex Groebner basis is [1]"
        inconsistent_flag = True
    elif timed_out and obstruction["decisive"]:
        verdict = "INCONSISTENT"
        inconsistency_method = (
            "support obstruction is decisive; the Groebner run timed out")
        inconsistent_flag = True
    elif timed_out:
        verdict = "UNKNOWN (timeout)"
        inconsistency_method = "Groebner run timed out; not decided"
        inconsistent_flag = None
    else:
        inconsistent_flag = False
        branches = stages["3_solutions"].get("branches", [])
        non_tri = [b for b in branches if b.get("non_triangular")]
        if not non_tri:
            verdict = "NO NONTRIVIAL SOLUTION"
        else:
            verdict = "COUNTEREXAMPLE CANDIDATE"
        inconsistency_method = "Groebner basis is not [1]"

    accept_stage = stages["5_accept"]
    record = {
        "problem": "jacobian-75-125-vertices",
        "case": "F_2(3,5), j=1",
        "corner": {
            "Delta_prime": [list(p) for p in DELTA_PRIME],
            "N_P": [list(p) for p in P_CORNERS],
            "N_Q": [list(q) for q in Q_CORNERS],
            "normalisation": "[P,Q] = x^2",
            "reduced_degrees": [15, 25],
        },
        "support_convention": "vertices_only (four corners per component)",
        "target": "det J(P,Q) - x^2 = 0",
        "unknowns": len(data["unknowns"]),
        "equations": len(data["equations"]),
        "equation_monomials": [list(m) for m in data["monomials"]],
        "det_monomials": [list(m) for m in data["det_monomials"]],
        "basis_size": stages["2_groebner"]["basis_size"],
        "basis": stages["2_groebner"]["basis"],
        "intermediate_basis_size": stages["2_groebner"].get(
            "intermediate_basis_size"),
        "inconsistent": inconsistent_flag,
        "inconsistency_method": inconsistency_method,
        "timed_out": bool(timed_out),
        "support_obstruction": obstruction,
        "solutions": (
            [] if solutions_record is None
            else stages["3_solutions"].get("branches", [])),
        "collisions": stages["4_collisions"].get("collisions", []),
        "accept": accept_stage,
        "accept_verdict": (
            accept_stage.get("stdout") if accept_stage.get("ran") else None),
        "accept_accepted": accept_stage.get("accepted"),
        "stages": stages,
        "verdict": verdict,
        "verdict_detail": _verdict_detail(verdict, obstruction, timed_out),
        "seconds": None,
    }
    record["seconds"] = round(time.perf_counter() - t_start, 3)

    if write:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(record, indent=2, default=str),
                            encoding="utf-8")
    log_stage(6, {
        "status": "written",
        "path": str(out_path),
        "verdict": verdict,
        "basis_size": record["basis_size"],
        "inconsistent": record["inconsistent"],
        "seconds": record["seconds"],
    })
    return record


def _verdict_detail(verdict, obstruction, timed_out):
    if verdict == "INCONSISTENT":
        return (
            "det J = x^2 has no solution on the vertices-only support.  The "
            "reduced grevlex basis is [1].  Independently, every monomial of "
            "det J has x-degree >= %d, so the x^2 coefficient equation is "
            "-1 = 0.  This is a support obstruction of the minimal-support "
            "specialization and does not exclude the full-lattice system."
            % obstruction["min_x_degree_of_det_J"]
        )
    if verdict == "UNKNOWN (timeout)":
        return (
            "The Groebner computation exceeded its cap and the support "
            "obstruction is not decisive; no verdict is claimed."
        )
    if verdict == "NO NONTRIVIAL SOLUTION":
        return (
            "The system is consistent, but every extracted solution is "
            "triangular or fails exact verification."
        )
    return (
        "A non-triangular solution with an exact collision was found and "
        "submitted to tools/accept.py; see accept.stdout for the verbatim "
        "gate verdict."
    )


def print_summary(record):
    print("== Jacobian corner (75,125), vertices only ==")
    print("Delta'            : {(0,0),(3,0),(4,1),(0,5)}")
    print("N(P)=3*Delta'     : {(0,0),(9,0),(12,3),(0,15)}")
    print("N(Q)=5*Delta'     : {(0,0),(15,0),(20,5),(0,25)}")
    print("target            : det J(P,Q) = x^2")
    print("unknowns          : %d" % record["unknowns"])
    print("equations         : %d" % record["equations"])
    print("Groebner order    : grevlex over Q")
    print("basis size        : %s" % record["basis_size"])
    print("basis             : %s" % record["basis"])
    print("inconsistent      : %s" % record["inconsistent"])
    print("min x-degree detJ : %d"
          % record["support_obstruction"]["min_x_degree_of_det_J"])
    print("x^2 in det J      : %s"
          % record["support_obstruction"]["target_monomial_present_in_det_J"])
    print("verdict           : %s" % record["verdict"])
    print("detail            : %s" % record["verdict_detail"])
    if record.get("accept", {}).get("ran"):
        print("accept exit code  : %s" % record["accept"]["exit_code"])
        print("accept stdout     : %s" % record["accept"]["stdout"].strip())
    print("seconds           : %.3f" % record["seconds"])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--timeout", type=float, default=GROEBNER_TIMEOUT_S,
                    help="seconds allowed for the Groebner basis (default 900)")
    ap.add_argument("--no-write", action="store_true")
    ap.add_argument("--gb-child", action="store_true", help=argparse.SUPPRESS)
    args = ap.parse_args(argv)

    if args.gb_child:
        return _gb_child()

    record = run(args.out, args.timeout, write=not args.no_write)
    print_summary(record)
    if not args.no_write:
        try:
            shown = args.out.relative_to(ROOT)
        except ValueError:
            shown = args.out
        print("\nwritten: %s" % shown)
        print("log    : %s" % LOG_PATH)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
