"""Solve the vertices-only GGHV config-2 coefficient system exactly.

Source: J. A. Guccione, J. J. Guccione, R. Horruitiner, C. Valqui,
"Increasing the degree of a possible counterexample to the Jacobian Conjecture
from 100 to 108", arXiv:2204.14178, Proposition 4.3 (Case (8,28)), case (2):

    N(P) = {(0,0), (1,0), (8,14), (8,16)}
    N(Q) = {(0,0), (2,1), (12,21), (12,24)}

This script uses the *minimal-support* (vertices-only) reading kept as a variant
in `config2_system.py`:

    P = a0 + a1*x + a2*x^8*y^14 + a3*x^8*y^16
    Q = b0 + b1*x^2*y + b2*x^12*y^21 + b3*x^12*y^24

and builds the coefficient system of the paper's normalisation [P,Q] = x^2,
i.e. det J(P,Q) - x^2 = 0.  That is 8 unknowns and 7 equations.

Result (proved by the Groebner basis below, not merely searched):

    det J - x^2 = (a1*b1 - 1)*x^2 - 20*a2*b1*x^9*y^14
                  - 24*a3*b1*x^9*y^16 + 21*a1*b2*x^12*y^20
                  + 24*a1*b3*x^12*y^23 - 24*a3*b2*x^19*y^36
                  + 24*a2*b3*x^19*y^37

    reduced grevlex basis over Q:  [a1*b1 - 1, a2, a3, b2, b3]

so the system is *consistent* and the solution set is a 3-dimensional family

    a0, b0 free,  b1 != 0 free,  a1 = 1/b1,  a2 = a3 = b2 = b3 = 0:

    P = a0 + x/b1,   Q = b0 + b1*x^2*y,   det J = x^2.

Every solution is therefore degenerate: the four genuine vertices (8,14),
(8,16) of N(P) and (12,21), (12,24) of N(Q) all carry coefficient zero, so no
solution has the required Newton polygons.  Saturating by a2*a3*b2*b3 gives
[1], i.e. the strict-support system (all eight coefficients nonzero) is
inconsistent.

The family is also not a Keller map: det J = x^2 is not a nonzero constant.
It is injective off the critical line x = 0 and collapses that entire line to
the single point (a0, b0), so it does carry an exact collision -- but on the
locus where the Jacobian drops rank, which is outside the Jacobian conjecture.

Reproduce:  python problems/jacobian/solve_vertices.py
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DEFAULT_OUT = ROOT / "runs" / "jacobian_config2_vertices.json"

# Reuse the exact supports and polynomial builder of config2_system.py so the
# system here cannot drift from the sized system.
sys.path.insert(0, str(HERE))
from config2_system import P_CORNERS, Q_CORNERS, build_polynomial  # noqa: E402

GROEBNER_TIMEOUT_S = 600.0          # hard cap on the decisive computation
TOTAL_BUDGET_S = 690.0             # acceptance harness allows 700 s
BUCHBERGER_PROGRESS_S = 60.0        # independent progress-tracking run

SYMBOL_NAMES = ["a0", "a1", "a2", "a3", "b0", "b1", "b2", "b3"]


# --------------------------------------------------------------------------
# the system
# --------------------------------------------------------------------------

def coefficients():
    return tuple(sp.Symbol(n) for n in SYMBOL_NAMES)


def build_system():
    """det J(P,Q) - x^2 on the four corners of each polygon.

    Same construction as config2_system.size_system(..., target=x**2):
    build P and Q from the support lists, expand the Jacobian determinant,
    then read off one equation per monomial of det J - x^2.
    """
    x, y = sp.symbols("x y")
    P, cp = build_polynomial(list(P_CORNERS), "a", x, y)
    Q, cq = build_polynomial(list(Q_CORNERS), "b", x, y)
    unknowns = cp + cq
    det = sp.expand(P.diff(x) * Q.diff(y) - P.diff(y) * Q.diff(x))
    poly = sp.Poly(det - x ** 2, x, y)
    monomials = [(int(m[0]), int(m[1])) for m in poly.monoms()]
    equations = [sp.expand(poly.coeff_monomial(m)) for m in monomials]
    return {
        "x": x, "y": y, "P": P, "Q": Q, "det": det,
        "unknowns": unknowns, "equations": equations, "monomials": monomials,
    }


# --------------------------------------------------------------------------
# independent Buchberger run (deadline-bounded, records intermediate sizes)
# --------------------------------------------------------------------------

def _mono_key(m):
    """grevlex key: bigger total degree first, then smaller last exponent."""
    return (sum(m), tuple(-e for e in reversed(m)))


def _from_expr(expr, gens):
    p = sp.Poly(expr, *gens)
    return {tuple(int(e) for e in mon): sp.Rational(c)
            for mon, c in zip(p.monoms(), p.coeffs())}


def _add(p, q, sign=1):
    out = dict(p)
    for m, c in q.items():
        v = out.get(m, sp.Integer(0)) + sign * c
        if v == 0:
            out.pop(m, None)
        else:
            out[m] = v
    return out


def _shift_scale(p, shift, factor):
    return {tuple(a + b for a, b in zip(m, shift)): c * factor
            for m, c in p.items()}


def _normal_form(p, basis):
    r, rem = dict(p), {}
    while r:
        m = max(r, key=_mono_key)
        c = r[m]
        for g in basis:
            mg = max(g, key=_mono_key)
            if all(a >= b for a, b in zip(m, mg)):
                shift = tuple(a - b for a, b in zip(m, mg))
                r = _add(r, _shift_scale(g, shift, c / g[mg]), sign=-1)
                break
        else:
            rem[m] = c
            r.pop(m)
    return rem


def _monic(p):
    lc = p[max(p, key=_mono_key)]
    return {m: c / lc for m, c in p.items()}


def _spoly(f, g):
    mf, mg = max(f, key=_mono_key), max(g, key=_mono_key)
    lcm = tuple(max(a, b) for a, b in zip(mf, mg))
    sf = tuple(a - b for a, b in zip(lcm, mf))
    sg = tuple(a - b for a, b in zip(lcm, mg))
    return _add(_shift_scale(f, sf, 1), _shift_scale(g, sg, 1), sign=-1)


def buchberger_progress(eqs, gens, deadline):
    """Run Buchberger over Q with grevlex, recording the basis size as it grows.

    This is a belt-and-braces check on sympy's answer and, more importantly,
    the only way to report an *intermediate* basis size if the main run hits
    the time cap.  Returns a dict, never raises.
    """
    try:
        t0 = time.monotonic()
        basis = []
        for e in eqs:
            r = _normal_form(_from_expr(e, gens), basis)
            if r:
                basis.append(_monic(r))
        pairs = [(i, j) for i in range(len(basis)) for j in range(i)]
        trace, processed = [len(basis)], 0
        while pairs:
            if time.monotonic() > deadline:
                return {"ok": True, "timed_out": True, "basis_size": len(basis),
                        "basis": [sp.sstr(sp.expand(sum(
                            c * sp.prod(g ** e for g, e in zip(gens, m))
                            for m, c in p.items()))) for p in basis],
                        "trace": trace, "pairs_processed": processed,
                        "pairs_left": len(pairs),
                        "seconds": round(time.monotonic() - t0, 4)}
            i, j = pairs.pop(0)
            processed += 1
            r = _normal_form(_spoly(basis[i], basis[j]), basis)
            if r:
                basis.append(_monic(r))
                pairs.extend((k, len(basis) - 1) for k in range(len(basis) - 1))
                trace.append(len(basis))
        return {"ok": True, "timed_out": False, "basis_size": len(basis),
                "basis": [sp.sstr(sp.expand(sum(
                    c * sp.prod(g ** e for g, e in zip(gens, m))
                    for m, c in p.items()))) for p in basis],
                "trace": trace, "pairs_processed": processed, "pairs_left": 0,
                "seconds": round(time.monotonic() - t0, 4)}
    except Exception as exc:  # a diagnostic must never break the run
        return {"ok": False, "error": f"{type(exc).__name__}: {exc}"}


# --------------------------------------------------------------------------
# sympy Groebner basis in a killable child process
# --------------------------------------------------------------------------

def _gb_child() -> int:
    data = build_system()
    t0 = time.perf_counter()
    G = sp.groebner(data["equations"], *data["unknowns"], order="grevlex")
    sys.stdout.write(json.dumps({
        "basis": [str(e) for e in G.exprs],
        "basis_size": len(G.exprs),
        "order": "grevlex",
        "seconds": time.perf_counter() - t0,
    }))
    sys.stdout.flush()
    return 0


def groebner_capped(timeout):
    """sympy `groebner` over Q, grevlex, killed at `timeout` seconds."""
    t0 = time.perf_counter()
    try:
        proc = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "--gb-child"],
            capture_output=True, text=True, encoding="utf-8",
            timeout=timeout, cwd=str(HERE),
        )
    except subprocess.TimeoutExpired:
        return {"timed_out": True, "basis": None, "basis_size": None,
                "seconds": round(time.perf_counter() - t0, 3)}
    if proc.returncode != 0:
        raise RuntimeError("groebner child failed:\n" + (proc.stderr or "")[-4000:])
    text = (proc.stdout or "").strip()
    payload = json.loads(text.splitlines()[-1] if text else "{}")
    payload["timed_out"] = False
    payload["seconds"] = round(time.perf_counter() - t0, 3)
    return payload


# --------------------------------------------------------------------------
# analysis of a computed basis
# --------------------------------------------------------------------------

def _as_expr(text, gens):
    return sp.sympify(text, locals={s.name: s for s in gens})


def krull_dimension(basis_exprs, gens):
    """dim Q[a]/I via the leading-term ideal: max independent set of variables.

    For a monomial ideal the dimension is the largest set S of variables such
    that no minimal generator of the ideal is supported inside S.
    """
    supports = []
    for e in basis_exprs:
        p = sp.Poly(sp.expand(e), *gens)
        for m in p.monoms():
            s = {i for i, k in enumerate(m) if k}
            if s:                      # ignore the constant monomial
                supports.append(s)
    n = len(gens)
    best = 0
    for mask in range(1 << n):
        chosen = {i for i in range(n) if mask >> i & 1}
        if all(not s <= chosen for s in supports):
            best = max(best, len(chosen))
    return best


def eliminate(basis_exprs, keep, gens):
    """Elimination ideal I cap Q[keep], read off a lex basis (elimination order)."""
    return [sp.sstr(e) for e in basis_exprs
            if e.free_symbols <= set(keep)]


# --------------------------------------------------------------------------
# solution extraction and verification
# --------------------------------------------------------------------------

def solve_family(data, basis_exprs):
    """Extract the parametrised solution set and verify it exactly."""
    gens = list(data["unknowns"])
    eqs = data["equations"]
    x, y = data["x"], data["y"]

    branches = sp.solve(eqs, gens, dict=True)
    out = []
    for index, sub in enumerate(branches):
        residuals = [sp.simplify(sp.expand(e.subs(sub))) for e in eqs]
        free = [s for s in gens if s not in sub]
        # exactness certificate: the ideal of the parametrisation must have the
        # same reduced Groebner basis as the system itself
        claim, denominators = [], []
        for v in gens:
            if v in sub:
                num, den = sp.fraction(sp.together(v - sub[v]))
                denominators.append(den)
                num = sp.expand(num)
                if num != 0 and not num.is_number:
                    claim.append(num)
        param_symbols = set()
        for d in denominators:
            param_symbols |= d.free_symbols
        claim_basis = sp.groebner(claim, *gens, order="grevlex").exprs if claim else []
        same_ideal = list(claim_basis) == list(basis_exprs)

        P_sub = sp.simplify(data["P"].subs(sub))
        Q_sub = sp.simplify(data["Q"].subs(sub))
        det_sub = sp.expand(data["det"].subs(sub))

        def used(expr):
            p = sp.Poly(sp.expand(expr), x, y)
            return [[int(a), int(b)] for a, b in p.monoms()]

        out.append({
            "index": index,
            "parametrization": {s.name: (sp.sstr(sub[s]) if s in sub else "free")
                                for s in gens},
            "free_parameters": [s.name for s in free],
            "nonzero_conditions": [f"{s.name} != 0" for s in free
                                   if s in param_symbols],
            "P": sp.sstr(P_sub),
            "Q": sp.sstr(Q_sub),
            "det_J": sp.sstr(det_sub),
            "det_J_minus_x2": sp.sstr(sp.expand(det_sub - x ** 2)),
            "det_J_equals_x2": sp.expand(det_sub - x ** 2) == 0,
            "residuals_zero": all(r == 0 for r in residuals),
            "parametrization_ideal_equals_system_ideal": bool(same_ideal),
            "support_used": {"P": used(P_sub), "Q": used(Q_sub)},
            "uses_all_four_corners": (
                set(map(tuple, used(P_sub))) == set(P_CORNERS)
                and set(map(tuple, used(Q_sub))) == set(Q_CORNERS)),
        })
    return out


def representative_check(values, data):
    """Build P, Q at one exact point of the family and test everything."""
    import jac

    x, y = data["x"], data["y"]
    gens = list(data["unknowns"])
    sub = dict(zip(gens, [sp.sympify(values[s.name]) for s in gens]))
    P, Q = sp.expand(data["P"].subs(sub)), sp.expand(data["Q"].subs(sub))

    det = sp.expand(P.diff(x) * Q.diff(y) - P.diff(y) * Q.diff(x))
    residuals = [sp.simplify(e.subs(sub)) for e in data["equations"]]

    # exact collision on the critical line x = 0: (0,0) and (0,1)
    p, q = (sp.Integer(0), sp.Integer(0)), (sp.Integer(0), sp.Integer(1))
    img_p = tuple(sp.simplify(f.subs({x: p[0], y: p[1]})) for f in (P, Q))
    img_q = tuple(sp.simplify(f.subs({x: q[0], y: q[1]})) for f in (P, Q))

    hits = jac.find_collision([P, Q], [x, y], samples=6, starts=4, seed=0)
    verdict = jac.verify_map([P, Q], [x, y])

    return {
        "values": {s.name: sp.sstr(sub[s]) for s in gens},
        "P": sp.sstr(P),
        "Q": sp.sstr(Q),
        "det_J": sp.sstr(det),
        "residuals_zero": all(r == 0 for r in residuals),
        "det_J_equals_x2": sp.expand(det - x ** 2) == 0,
        "exact_collision": {
            "points": [[sp.sstr(c) for c in p], [sp.sstr(c) for c in q]],
            "distinct": p != q,
            "images": [[sp.sstr(c) for c in img_p],
                       [sp.sstr(c) for c in img_q]],
            "same_image": img_p == img_q,
        },
        "numerical_fibre_search": {"hits": hits, "samples": 6, "starts": 4},
        "jac_verify_map": verdict,
    }


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------

def run(out_path: Path, timeout: float):
    t_start = time.perf_counter()
    data = build_system()
    gens = list(data["unknowns"])
    eqs = data["equations"]

    progress = buchberger_progress(
        eqs, gens, deadline=time.monotonic() + BUCHBERGER_PROGRESS_S)

    remaining = TOTAL_BUDGET_S - (time.perf_counter() - t_start)
    gb = groebner_capped(max(1.0, min(timeout, remaining)))

    record = {
        "problem": "jacobian-config2-vertices",
        "arxiv": "2204.14178",
        "proposition": "4.3, case (2)",
        "variant": "vertices_only",
        "support_convention": "the four corners of each Newton polygon are the supports",
        "supports": {
            "P": [list(p) for p in P_CORNERS],
            "Q": [list(q) for q in Q_CORNERS],
        },
        "target": "det J(P,Q) = x^2",
        "target_polynomial": "det J(P,Q) - x**2",
        "unknowns": len(gens),
        "equations": len(eqs),
        "equation_monomials": [list(m) for m in data["monomials"]],
        "equation_list": [
            {"monomial": list(m), "polynomial": sp.sstr(e)}
            for m, e in zip(data["monomials"], eqs)
        ],
        "groebner_order": "grevlex",
        "timeout_seconds": timeout,
        "timed_out": bool(gb["timed_out"]),
        "seconds": None,
    }

    if gb["timed_out"]:
        record.update({
            "inconsistent": "unknown",
            "basis_size": None,
            "basis": None,
            "intermediate_basis_size": progress.get("basis_size"),
            "intermediate_basis_trace": progress.get("trace"),
            "solutions": None,
            "note": ("Groebner computation exceeded the cap; the intermediate "
                     "basis size comes from the deadline-bounded Buchberger "
                     "cross-check."),
        })
        record["seconds"] = round(time.perf_counter() - t_start, 3)
        return record

    basis_exprs = [_as_expr(s, gens) for s in gb["basis"]]
    record["basis"] = [sp.sstr(e) for e in basis_exprs]
    record["basis_size"] = len(basis_exprs)

    # --- the decisive branch ------------------------------------------------
    if basis_exprs == [sp.Integer(1)]:
        record.update({
            "inconsistent": True,
            "solutions": [],
            "solution_dimension": -1,
            "note": "reduced Groebner basis is [1]: no solution over C",
        })
        record["seconds"] = round(time.perf_counter() - t_start, 3)
        return record

    record["inconsistent"] = False

    # certificate that the returned basis really is a basis of the system
    G = sp.groebner(eqs, *gens, order="grevlex")
    record["certificate"] = {
        "each_equation_reduces_to_zero": all(G.reduce(e)[1] == 0 for e in eqs),
        "basis_is_reduced_groebner_basis": list(G.exprs) == basis_exprs,
        "independently_recomputed_basis_size": len(G.exprs),
    }

    lex = sp.groebner(eqs, *gens, order="lex")
    dim = krull_dimension(basis_exprs, gens)
    record["solution_dimension"] = dim
    record["elimination_ideal"] = {
        "keep_a1_b1": eliminate(lex.exprs, [sp.Symbol("a1"), sp.Symbol("b1")], gens),
        "keep_a2_a3_b2_b3": eliminate(
            lex.exprs, [sp.Symbol("a2"), sp.Symbol("a3"),
                        sp.Symbol("b2"), sp.Symbol("b3")], gens),
        "note": "read off a lex basis (elimination order)",
    }

    branches = solve_family(data, basis_exprs)
    record["solutions"] = {
        "dimension": dim,
        "number_of_branches": len(branches),
        "branches": branches,
    }

    # --- strict support: all eight corner coefficients nonzero --------------
    t_sat = sp.Symbol("t")
    a2, a3, b2, b3 = (sp.Symbol(n) for n in ("a2", "a3", "b2", "b3"))
    sat = sp.groebner(eqs + [1 - t_sat * a2 * a3 * b2 * b3], *gens, t_sat,
                      order="grevlex")
    strict_inconsistent = list(sat.exprs) == [sp.Integer(1)]
    record["strict_support"] = {
        "condition": "a2*a3*b2*b3 != 0 (all four corners of each polygon used)",
        "saturation_basis": [sp.sstr(e) for e in sat.exprs],
        "inconsistent": bool(strict_inconsistent),
        "meaning": ("no solution of det J = x^2 has all four corners of N(P) "
                    "and all four corners of N(Q) present"
                    if strict_inconsistent else
                    "a strict-support solution exists"),
    }

    # --- representatives, det J check, collision test ----------------------
    reps = []
    for values in (
        {"a0": 0, "a1": 1, "a2": 0, "a3": 0, "b0": 0, "b1": 1, "b2": 0, "b3": 0},
        {"a0": 1, "a1": 2, "a2": 0, "a3": 0, "b0": 3,
         "b1": sp.Rational(1, 2), "b2": 0, "b3": 0},
    ):
        reps.append(representative_check(values, data))
    record["solutions"]["representatives"] = reps
    record["solutions"]["full_support_solutions"] = not strict_inconsistent

    branch = branches[0] if branches else {}
    collision_ok = all(r["exact_collision"]["same_image"] for r in reps)
    record["collision"] = {
        "critical_locus": "x = 0, where det J = x^2 vanishes",
        "exact_collision_on_critical_locus": bool(collision_ok),
        "generic_fibre_hits": [r["numerical_fibre_search"]["hits"] for r in reps],
        "generic_fibre_search_note": (
            "jac.find_collision samples random targets.  Every fibre with "
            "x != 0 is a single point, so no generic hit is expected; the "
            "collision lives on the critical value (a0, b0)."),
        "injective_off_critical_locus": True,
        "injectivity_argument": (
            "P = a0 + a1*x with a1 != 0 determines x; for x != 0, "
            "Q = b0 + b1*x^2*y determines y. Hence the only collisions are on "
            "x = 0, where the whole line maps to (a0, b0)."
        ),
        "keller_counterexample": False,
        "reason": (
            "det J = x^2 is not a non-zero constant, so the Jacobian "
            "conjecture does not apply; the collision lies on the critical "
            "locus where the Jacobian drops rank."
        ),
        "det_J_in_branch": branch.get("det_J"),
    }

    record["verdict"] = (
        "CONSISTENT_BUT_DEGENERATE: det J = x^2 has a 3-parameter solution "
        "family, every member of which sets a2 = a3 = b2 = b3 = 0, so the four "
        "genuine vertices are absent and no P, Q have the required Newton "
        "polygons.  The family is P = a0 + a1*x, Q = b0 + b1*x^2*y with "
        "a1*b1 = 1: it is not a Keller map (det J = x^2) and it collapses the "
        "line x = 0, so it is not a counterexample either."
    )
    record["seconds"] = round(time.perf_counter() - t_start, 3)
    return record


def print_summary(rec):
    print("== GGHV 2022 arXiv:2204.14178, Prop. 4.3 case (2), vertices only ==")
    print("unknowns            : %d" % rec["unknowns"])
    print("equations           : %d" % rec["equations"])
    print("target              : %s" % rec["target"])
    print("groebner order      : %s" % rec["groebner_order"])
    if rec.get("timed_out"):
        print("STATUS              : TIMEOUT after %.1f s" % rec["seconds"])
        print("intermediate basis  : %s" % rec.get("intermediate_basis_size"))
        return
    print("basis size          : %d" % rec["basis_size"])
    print("basis               : %s" % ", ".join(rec["basis"]))
    print("inconsistent        : %s" % rec["inconsistent"])
    if not rec["inconsistent"]:
        print("solution dimension  : %d" % rec["solution_dimension"])
        for br in rec["solutions"]["branches"]:
            print("  P = %s" % br["P"])
            print("  Q = %s" % br["Q"])
            print("  det J = %s   (equals x^2: %s)"
                  % (br["det_J"], br["det_J_equals_x2"]))
        print("strict support      : inconsistent = %s"
              % rec["strict_support"]["inconsistent"])
        for rep in rec["solutions"]["representatives"]:
            print("  rep %s -> det J = %s, exact collision %s"
                  % (rep["values"], rep["det_J"], rep["exact_collision"]["same_image"]))
        print("verdict             : %s" % rec["verdict"])
    print("seconds             : %.3f" % rec["seconds"])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--timeout", type=float, default=GROEBNER_TIMEOUT_S,
                    help="seconds allowed for the Groebner basis (default 600)")
    ap.add_argument("--no-write", action="store_true")
    ap.add_argument("--gb-child", action="store_true", help=argparse.SUPPRESS)
    args = ap.parse_args(argv)

    if args.gb_child:
        return _gb_child()

    rec = run(args.out, args.timeout)
    print_summary(rec)
    if not args.no_write:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(rec, indent=2), encoding="utf-8")
        print("\nwritten: %s" % args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
