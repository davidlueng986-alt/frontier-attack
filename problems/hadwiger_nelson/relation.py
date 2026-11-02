"""Forced colour relations in k-colourings -- Hadwiger-Nelson track.

Given a graph and two vertices, decide what relation their colours have in
*every* proper k-colouring:

    "unsat"      the graph has no proper k-colouring at all
    "equal"      k-colourable, but every k-colouring gives u and v one colour
    "different"  k-colourable, but every k-colouring gives u and v distinct colours
    "free"       both possibilities occur

Encoding: one boolean variable x[v][c] per vertex/colour, exactly one true per
vertex, and for every edge uv and colour c: not (x[u][c] and x[v][c]).
SAT solving is delegated to pysat's cadical153.

The pair scan uses selector literals and assumptions, so one incremental solver
answers many pair queries without rebuilding the formula.  For k >= 2 the
formula contains no unit clauses; this keeps CaDiCaL's model checker quiet on
UNSAT instances.  The one-colour case is decided directly: a graph is
1-colourable exactly when it has no edges.

Usage:
    python problems/hadwiger_nelson/relation.py --selftest
    python problems/hadwiger_nelson/relation.py --graph509 --pair-cap 2000

The self-test compares `relation` with exhaustive enumeration on all 1024
labelled graphs on five vertices for k = 3.  `--selftest` also runs the 509
graph application, as required by the acceptance check.
"""

from __future__ import annotations

import argparse
import itertools
import json
import re
import time
from pathlib import Path

import mpmath as mp
import sympy as sp
from pysat.solvers import Solver

ROOT = Path(__file__).resolve().parents[2]
VTX_PATH = ROOT / "research" / "scratch" / "509_parts.vtx"
OUT_PATH = ROOT / "runs" / "hn_relations.json"

SOLVER_NAME = "cadical153"

# Pair-loop guards for the 509-vertex application.  Both are reported in the
# JSON; the wall-clock guard is what actually keeps a run below ten minutes.
DEFAULT_PAIR_CAP = 2000
DEFAULT_TIME_BUDGET = 420.0


# --------------------------------------------------------------------------
# CNF encoding
# --------------------------------------------------------------------------

def _var(v: int, c: int, k: int) -> int:
    """1-based DIMACS variable for 'vertex v has colour c'."""
    return v * k + c + 1


def _normalise_edges(edges, n: int) -> list[tuple[int, int]]:
    if n < 0:
        raise ValueError(f"n must be non-negative, got {n}")
    out: list[tuple[int, int]] = []
    seen: set[tuple[int, int]] = set()
    for e in edges:
        if len(e) != 2:
            raise ValueError(f"edge must be a pair, got {e!r}")
        u, v = int(e[0]), int(e[1])
        if not (0 <= u < n and 0 <= v < n):
            raise ValueError(f"edge {(u, v)} out of range for n={n}")
        if u == v:
            raise ValueError(f"self-loop at vertex {u}")
        if u > v:
            u, v = v, u
        if (u, v) not in seen:
            seen.add((u, v))
            out.append((u, v))
    return out


def colouring_cnf(n: int, edges, k: int) -> list[list[int]]:
    """Plain one-hot k-colouring CNF (no unit clauses for k >= 2)."""
    if k < 1:
        raise ValueError(f"k must be at least 1, got {k}")
    es = _normalise_edges(edges, n)
    clauses: list[list[int]] = []
    for v in range(n):
        clauses.append([_var(v, c, k) for c in range(k)])
        for c1 in range(k):
            for c2 in range(c1 + 1, k):
                clauses.append([-_var(v, c1, k), -_var(v, c2, k)])
    for u, v in es:
        for c in range(k):
            clauses.append([-_var(u, c, k), -_var(v, c, k)])
    return clauses


# --------------------------------------------------------------------------
# relation prober
# --------------------------------------------------------------------------

def _same_colour_clauses(u: int, v: int, k: int) -> list[list[int]]:
    """Clauses forcing x[u][c] <-> x[v][c] for every colour c."""
    out: list[list[int]] = []
    for c in range(k):
        out.append([-_var(u, c, k), _var(v, c, k)])
        out.append([-_var(v, c, k), _var(u, c, k)])
    return out


def _different_colour_clauses(u: int, v: int, k: int) -> list[list[int]]:
    """Clauses forcing u and v to receive different colours."""
    return [[-_var(u, c, k), -_var(v, c, k)] for c in range(k)]


def _solve(clauses, solver_name: str = SOLVER_NAME) -> bool:
    with Solver(name=solver_name, bootstrap_with=clauses) as solver:
        return bool(solver.solve())


def relation(edges, n: int, u: int, v: int, k: int) -> str:
    """Forced colour relation of vertices u and v in every proper k-colouring.

    Returns one of "equal", "different", "free", "unsat".
    """
    if k < 1:
        raise ValueError(f"k must be at least 1, got {k}")
    if not (0 <= u < n and 0 <= v < n):
        raise ValueError(f"vertex ({u}, {v}) out of range for n={n}")

    es = _normalise_edges(edges, n)
    if k == 1:
        # One colour only: colourable iff there are no edges, and then every
        # pair of vertices is forced equal.
        return "unsat" if es else "equal"

    base = colouring_cnf(n, es, k)

    if u == v:
        # The constraint u != v is unsatisfiable, so the pair is forced equal.
        return "equal" if _solve(base) else "unsat"

    # If u == v is already impossible, every colouring separates the pair.
    same = base + _same_colour_clauses(u, v, k)
    if not _solve(same):
        return "different" if _solve(base) else "unsat"

    # u == v is possible, so the graph is colourable.  If u != v is impossible,
    # every colouring identifies the pair.
    diff = base + _different_colour_clauses(u, v, k)
    if not _solve(diff):
        return "equal"

    return "free"


# --------------------------------------------------------------------------
# batch scan
# --------------------------------------------------------------------------

def scan_equal_pairs(
    edges,
    n: int,
    k: int,
    pair_cap: int = DEFAULT_PAIR_CAP,
    time_budget: float = DEFAULT_TIME_BUDGET,
) -> dict:
    """Scan vertex pairs and collect those forced equal, with explicit caps.

    The pair loop stops at `pair_cap` pairs or `time_budget` seconds, whichever
    comes first.  If the graph is not k-colourable no pair loop is run: every
    pair has relation "unsat".

    All selector clauses are added before the first solve.  The solver is then
    queried only through assumptions, which keeps the model checker quiet and
    the incremental state unambiguous.
    """
    if pair_cap < 0:
        raise ValueError("pair_cap must be non-negative")
    if time_budget <= 0:
        raise ValueError("time_budget must be positive")

    t0 = time.perf_counter()
    total = n * (n - 1) // 2
    cap = min(pair_cap, total)
    pairs = list(itertools.islice(itertools.combinations(range(n), 2), cap))
    es = _normalise_edges(edges, n)

    info: dict = {
        "k": k,
        "pairs_total": total,
        "pair_cap": pair_cap,
        "pairs_in_scope": len(pairs),
        "time_budget_s": time_budget,
        "pairs_examined": 0,
        "equal_pairs": [],
        "equal_pair_count": 0,
        "capped": cap < total,
        "cap_reason": "pair-cap" if cap < total else None,
        "base": None,
        "base_seconds": 0.0,
        "seconds": 0.0,
    }

    if k == 1:
        # One colour only: colourable iff there are no edges, and then every
        # pair in scope is forced equal.
        if es:
            info["base"] = "unsat"
            info["capped"] = False
            info["cap_reason"] = "base-unsat"
        else:
            info["base"] = "sat"
            info["pairs_examined"] = len(pairs)
            info["equal_pairs"] = list(pairs)
            info["equal_pair_count"] = len(pairs)
        info["seconds"] = time.perf_counter() - t0
        return info

    base = colouring_cnf(n, es, k)
    next_var = n * k + 1
    clauses = list(base)
    selectors: list[tuple[int, int, int, int]] = []
    for u, v in pairs:
        s_same = next_var
        s_diff = next_var + 1
        next_var += 2
        selectors.append((u, v, s_same, s_diff))
        for c in range(k):
            clauses.append([-s_same, -_var(u, c, k), _var(v, c, k)])
            clauses.append([-s_same, -_var(v, c, k), _var(u, c, k)])
            clauses.append([-s_diff, -_var(u, c, k), -_var(v, c, k)])

    with Solver(name=SOLVER_NAME, bootstrap_with=clauses) as solver:
        tb = time.perf_counter()
        base_sat = bool(solver.solve())
        info["base_seconds"] = time.perf_counter() - tb
        info["base"] = "sat" if base_sat else "unsat"

        if not base_sat:
            info["capped"] = False
            info["cap_reason"] = "base-unsat"
            info["seconds"] = time.perf_counter() - t0
            return info

        for u, v, s_same, s_diff in selectors:
            if time.perf_counter() - t0 >= time_budget:
                info["capped"] = True
                info["cap_reason"] = "time-budget"
                break
            info["pairs_examined"] += 1
            if not solver.solve(assumptions=[s_same]):
                continue  # forced different
            if not solver.solve(assumptions=[s_diff]):
                info["equal_pairs"].append((u, v))

        info["equal_pair_count"] = len(info["equal_pairs"])
        info["seconds"] = time.perf_counter() - t0

    return info


def forced_equal_pairs(
    edges,
    n: int,
    k: int,
    pair_cap: int = DEFAULT_PAIR_CAP,
    time_budget: float = DEFAULT_TIME_BUDGET,
) -> list[tuple[int, int]]:
    """All pairs (u, v) with relation "equal", subject to the loop caps."""
    return scan_equal_pairs(edges, n, k, pair_cap, time_budget)["equal_pairs"]


# --------------------------------------------------------------------------
# self-test against brute force
# --------------------------------------------------------------------------

def _all_graphs(n: int):
    pairs = list(itertools.combinations(range(n), 2))
    for mask in range(1 << len(pairs)):
        yield [pairs[i] for i in range(len(pairs)) if (mask >> i) & 1]


def brute_relations(edges, n: int, k: int) -> dict[tuple[int, int], str]:
    """Exhaustive relation table over all k**n colourings."""
    pairs = list(itertools.combinations(range(n), 2))
    has_eq = {p: False for p in pairs}
    has_ne = {p: False for p in pairs}
    proper = 0
    for cols in itertools.product(range(k), repeat=n):
        if any(cols[u] == cols[v] for u, v in edges):
            continue
        proper += 1
        for u, v in pairs:
            if cols[u] == cols[v]:
                has_eq[(u, v)] = True
            else:
                has_ne[(u, v)] = True

    if proper == 0:
        return {p: "unsat" for p in pairs}

    table: dict[tuple[int, int], str] = {}
    for p in pairs:
        if has_eq[p] and has_ne[p]:
            table[p] = "free"
        elif has_eq[p]:
            table[p] = "equal"
        elif has_ne[p]:
            table[p] = "different"
        else:  # pragma: no cover - impossible when a proper colouring exists
            raise AssertionError(f"no relation for pair {p}")
    return table


def selftest(k: int = 3, n: int = 5, verbose: bool = True) -> dict:
    """Compare `relation` with brute force on every labelled n-vertex graph."""
    t0 = time.perf_counter()
    graphs = 0
    queries = 0
    mismatches: list[dict] = []

    for edges in _all_graphs(n):
        graphs += 1
        expected = brute_relations(edges, n, k)
        for (u, v), want in expected.items():
            got = relation(edges, n, u, v, k)
            queries += 1
            if got != want:
                mismatches.append(
                    {
                        "edges": list(edges),
                        "pair": [u, v],
                        "got": got,
                        "want": want,
                    }
                )

    elapsed = time.perf_counter() - t0
    if verbose:
        print(
            f"self-test: graphs={graphs} n={n} k={k} "
            f"pair_queries={queries} mismatches={len(mismatches)} "
            f"({elapsed:.2f}s)",
            flush=True,
        )
        for m in mismatches[:10]:
            print(f"  MISMATCH {m}", flush=True)

    return {
        "graphs": graphs,
        "n": n,
        "k": k,
        "pair_queries": queries,
        "mismatches": mismatches,
        "seconds": elapsed,
    }


# --------------------------------------------------------------------------
# 509-vertex Parts graph
# --------------------------------------------------------------------------

def _split_top(body: str):
    depth = 0
    for i, ch in enumerate(body):
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
        elif ch == "," and depth == 0:
            return body[:i], body[i + 1 :]
    raise ValueError("no top-level comma")


def parse_point(line: str):
    """Parse one Mathematica-syntax point; same grammar as tools/verify_509.py."""
    s = line.strip()
    if not s:
        return None
    if not (s.startswith("{") and s.endswith("}")):
        raise ValueError(f"bad line: {s[:40]}")
    a, b = _split_top(s[1:-1])
    out = []
    for expr in (a, b):
        e = expr
        while "Sqrt[" in e:
            nxt = re.sub(r"Sqrt\[([^\[\]]*)\]", r"sqrt(\1)", e)
            if nxt == e:
                break
            e = nxt
        out.append(mp.mpf(str(sp.N(sp.sympify(e), 50))))
    return tuple(out)


def build_unit_edges(points, tol=mp.mpf("1e-40")):
    n = len(points)
    edges: list[tuple[int, int]] = []
    worst = mp.mpf(0)
    for i in range(n):
        xi, yi = points[i]
        for j in range(i + 1, n):
            xj, yj = points[j]
            d2 = (xi - xj) ** 2 + (yi - yj) ** 2
            dev = abs(d2 - 1)
            if dev <= tol:
                edges.append((i, j))
                if dev > worst:
                    worst = dev
    return edges, worst


def load_509(path: Path = VTX_PATH):
    mp.mp.dps = 50
    points = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            p = parse_point(line)
            if p is not None:
                points.append(p)
    edges, worst = build_unit_edges(points)
    return points, edges, worst


def run_509(
    pair_cap: int = DEFAULT_PAIR_CAP,
    time_budget: float = DEFAULT_TIME_BUDGET,
) -> dict:
    print(f"509 graph: parsing {VTX_PATH} ...", flush=True)
    t0 = time.perf_counter()
    points, edges, worst = load_509()
    print(
        f"509 graph: n={len(points)} m={len(edges)} "
        f"worst|d^2-1|={mp.nstr(worst, 6)} "
        f"({time.perf_counter() - t0:.1f}s)",
        flush=True,
    )

    results: dict[str, dict] = {}
    unsat_at = None

    # Solve k=4 first.  If that is UNSAT, k=3 is UNSAT by monotonicity, so the
    # second (potentially expensive) base solve can be skipped rigorously.
    for k in (4, 3):
        if unsat_at is not None:
            total = len(points) * (len(points) - 1) // 2
            results[str(k)] = {
                "k": k,
                "pairs_total": total,
                "pair_cap": pair_cap,
                "pairs_in_scope": 0,
                "time_budget_s": time_budget,
                "pairs_examined": 0,
                "equal_pairs": [],
                "equal_pair_count": 0,
                "capped": False,
                "cap_reason": "base-unsat-inferred",
                "base": "unsat",
                "base_seconds": 0.0,
                "seconds": 0.0,
                "inferred_from_k": unsat_at,
            }
            print(
                f"  k={k}: base=unsat (inferred from k={unsat_at}) "
                f"equal_pairs=0 pairs_examined=0/{total} cap={pair_cap}",
                flush=True,
            )
            continue

        print(f"  k={k}: solving base CNF ...", flush=True)
        info = scan_equal_pairs(
            edges, len(points), k, pair_cap=pair_cap, time_budget=time_budget
        )
        results[str(k)] = info
        print(
            f"  k={k}: base={info['base']} "
            f"equal_pairs={info['equal_pair_count']} "
            f"pairs_examined={info['pairs_examined']}/{info['pairs_total']} "
            f"cap={pair_cap} ({info['seconds']:.1f}s)",
            flush=True,
        )
        if info["base"] == "unsat":
            unsat_at = k

    summary = "; ".join(
        f"k={k}: {results[str(k)]['equal_pair_count']} equal pairs "
        f"({results[str(k)]['base']})"
        for k in (3, 4)
    )
    print(
        f"509 result: {summary} "
        f"[pair_cap={pair_cap}, time_budget={time_budget:.0f}s]",
        flush=True,
    )

    return {
        "path": str(VTX_PATH.relative_to(ROOT)).replace("\\", "/"),
        "n": len(points),
        "m": len(edges),
        "worst_abs_d2_minus_1": float(worst),
        "pair_cap": pair_cap,
        "time_budget_s": time_budget,
        "k": results,
    }


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Forced colour relations in k-colourings (Hadwiger-Nelson)."
    )
    ap.add_argument(
        "--selftest",
        action="store_true",
        help="run the brute-force self-test and the 509 graph application",
    )
    ap.add_argument(
        "--graph509",
        action="store_true",
        help="run only the 509 graph application",
    )
    ap.add_argument(
        "--selftest-only",
        action="store_true",
        help="run only the brute-force self-test",
    )
    ap.add_argument("--pair-cap", type=int, default=DEFAULT_PAIR_CAP)
    ap.add_argument("--time-budget", type=float, default=DEFAULT_TIME_BUDGET)
    ap.add_argument("--out", default=str(OUT_PATH))
    args = ap.parse_args(argv)

    if not (args.selftest or args.graph509 or args.selftest_only):
        args.selftest = True  # default: both parts

    do_selftest = args.selftest or args.selftest_only
    do_509 = args.selftest or args.graph509

    payload: dict = {"solver": SOLVER_NAME}
    st = None

    if do_selftest:
        st = selftest()
        payload["selftest"] = {
            "graphs": st["graphs"],
            "n": st["n"],
            "k": st["k"],
            "pair_queries": st["pair_queries"],
            "mismatches": len(st["mismatches"]),
            "mismatch_details": st["mismatches"][:20],
            "seconds": st["seconds"],
        }

    if do_509:
        payload["graph_509"] = run_509(args.pair_cap, args.time_budget)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"written: {out}", flush=True)

    if st is not None and st["mismatches"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
