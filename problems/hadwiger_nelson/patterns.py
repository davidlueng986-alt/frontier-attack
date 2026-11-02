"""Pattern forcing for finite graphs, with a SAT prober and a brute-force oracle.

The central question of this module is not "is this graph k-colourable" but
"does every proper k-colouring make a prescribed vertex pattern monochromatic?".

For a triple T = (a, b, c) the encoding is direct.  A colouring makes T
monochromatic iff there is a colour ``col`` with x[a][col] = x[b][col] =
x[c][col] = true.  Since the CNF uses one-hot colour variables, the negation is
the k clauses

    (-x[a][col] OR -x[b][col] OR -x[c][col])        col = 0 .. k-1

so a single extra SAT call decides whether a non-monochromatic colouring of T
exists.  This is exact; no enumeration and no tolerance is involved.

Conventions
-----------
``forced_mono`` returns True iff

  * the graph has at least one proper k-colouring, and
  * every proper k-colouring makes the triple monochromatic.

An uncolourable graph returns False.  Logically "every colouring" would be
vacuously true, but that answer carries no information about the pattern and
would make every UNSAT instance look like a gadget.  The brute-force oracle
below uses the same convention, so the validation is a genuine cross-check.

The module is intentionally dependency-light: python-sat is the only non-standard
import.  ``validate_bruteforce`` is the acceptance test for the prober; it
enumerates *all* labelled graphs on n <= 5 vertices and all triples, at k = 3,
and compares the SAT answer with direct enumeration of the 3^n colourings.

Run ``python problems/hadwiger_nelson/patterns.py`` for the self-check.
"""

from __future__ import annotations

import itertools
import time

from pysat.formula import CNF
from pysat.solvers import Solver

SOLVER_NAME = "cadical153"


# --------------------------------------------------------------------------
# CNF helpers
# --------------------------------------------------------------------------

def _var(v: int, colour: int, k: int) -> int:
    """1-based DIMACS variable for 'vertex v has colour colour'."""
    return v * k + colour + 1


def _normalise_edges(edges, n: int) -> list[tuple[int, int]]:
    """Validate and canonicalise an edge list (deduplicated, u < v)."""
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


def colouring_cnf(n: int, edges, k: int) -> CNF:
    """One-hot CNF for proper k-colourings of (V, edges)."""
    if k < 1:
        raise ValueError(f"k must be at least 1, got {k}")
    es = _normalise_edges(edges, n)
    cnf = CNF()
    for v in range(n):
        cnf.append([_var(v, c, k) for c in range(k)])
        for c1 in range(k):
            for c2 in range(c1 + 1, k):
                cnf.append([-_var(v, c1, k), -_var(v, c2, k)])
    for u, v in es:
        for c in range(k):
            cnf.append([-_var(u, c, k), -_var(v, c, k)])
    return cnf


def _check_triple(triple, n: int) -> tuple[int, int, int]:
    if len(triple) != 3:
        raise ValueError(f"triple must have three vertices, got {triple!r}")
    a, b, c = (int(x) for x in triple)
    if not all(0 <= x < n for x in (a, b, c)):
        raise ValueError(f"triple {(a, b, c)} out of range for n={n}")
    if len({a, b, c}) != 3:
        raise ValueError(f"triple vertices must be distinct, got {(a, b, c)}")
    return a, b, c


def _not_mono_clauses(triple, k: int) -> list[list[int]]:
    a, b, c = triple
    return [[-_var(a, col, k), -_var(b, col, k), -_var(c, col, k)]
            for col in range(k)]


# --------------------------------------------------------------------------
# the prober
# --------------------------------------------------------------------------

def forced_mono(edges, n: int, triple, k: int) -> bool:
    """True iff every proper k-colouring makes ``triple`` monochromatic.

    Returns False when the graph is not k-colourable (see module docstring).
    """
    a, b, c = _check_triple(triple, n)
    cnf = colouring_cnf(n, edges, k)

    with Solver(name=SOLVER_NAME, bootstrap_with=cnf) as solver:
        if not solver.solve():
            return False
        solver.append_formula(_not_mono_clauses((a, b, c), k))
        return not solver.solve()


def forced_equal(edges, n: int, u: int, v: int, k: int) -> bool:
    """True iff every proper k-colouring gives u and v the same colour.

    This is the two-terminal relation of the diamond.  It is included because
    the 509 census is phrased in these terms; it is not used by the triple
    gadget search.  Returns False when the graph is not k-colourable.
    """
    if not (0 <= u < n and 0 <= v < n):
        raise ValueError(f"pair {(u, v)} out of range for n={n}")
    cnf = colouring_cnf(n, edges, k)
    with Solver(name=SOLVER_NAME, bootstrap_with=cnf) as solver:
        if not solver.solve():
            return False
        if u == v:
            return True
        solver.append_formula([[-_var(u, c, k), -_var(v, c, k)]
                               for c in range(k)])
        return not solver.solve()


def forces_mono_triple_set(edges, n: int, triples, k: int) -> bool:
    """True iff every proper k-colouring makes at least one triple monochromatic.

    This is the pattern-forcing primitive used by the gadget search: a gadget
    "forces a monochromatic sqrt(3)-triangle" when this returns True for the
    family of its three-point sqrt(3) triangles.  The empty family returns False.

    Returns False when the graph is not k-colourable.
    """
    clean: list[tuple[int, int, int]] = []
    seen: set[tuple[int, int, int]] = set()
    for t in triples:
        a, b, c = _check_triple(t, n)
        key = (a, b, c)
        if key not in seen:
            seen.add(key)
            clean.append(key)
    if not clean:
        return False

    cnf = colouring_cnf(n, edges, k)
    with Solver(name=SOLVER_NAME, bootstrap_with=cnf) as solver:
        if not solver.solve():
            return False
        extra: list[list[int]] = []
        for t in clean:
            extra.extend(_not_mono_clauses(t, k))
        solver.append_formula(extra)
        return not solver.solve()


# --------------------------------------------------------------------------
# brute-force oracle
# --------------------------------------------------------------------------

def brute_forced_mono(edges, n: int, triple, k: int) -> bool:
    """Direct enumeration oracle with the same convention as forced_mono."""
    a, b, c = _check_triple(triple, n)
    es = _normalise_edges(edges, n)
    found_proper = False
    for colouring in itertools.product(range(k), repeat=n):
        if any(colouring[u] == colouring[v] for u, v in es):
            continue
        found_proper = True
        if not (colouring[a] == colouring[b] == colouring[c]):
            return False
    return found_proper


def all_graphs(n: int):
    """Yield every labelled simple graph on n vertices as an edge list."""
    pairs = list(itertools.combinations(range(n), 2))
    for mask in range(1 << len(pairs)):
        yield [pairs[i] for i in range(len(pairs)) if (mask >> i) & 1]


def validate_bruteforce(max_n: int = 5, k: int = 3) -> dict:
    """Compare forced_mono with brute force on all graphs n <= max_n.

    Returns a report dict.  ``mismatches`` must be 0.
    """
    t0 = time.time()
    graphs = 0
    checks = 0
    mismatches = 0
    examples: list[dict] = []

    for n in range(3, max_n + 1):
        triples = list(itertools.combinations(range(n), 3))
        for edges in all_graphs(n):
            graphs += 1
            for triple in triples:
                checks += 1
                sat_answer = forced_mono(edges, n, triple, k)
                brute_answer = brute_forced_mono(edges, n, triple, k)
                if sat_answer != brute_answer:
                    mismatches += 1
                    if len(examples) < 10:
                        examples.append({
                            "n": n,
                            "edges": edges,
                            "triple": triple,
                            "sat": sat_answer,
                            "brute": brute_answer,
                        })

    return {
        "max_n": max_n,
        "k": k,
        "graphs": graphs,
        "checks": checks,
        "mismatches": mismatches,
        "mismatch_examples": examples,
        "seconds": time.time() - t0,
    }


# --------------------------------------------------------------------------
# self-check
# --------------------------------------------------------------------------

def _selftest() -> None:
    # diamond: the two tips are equal in every 3-colouring, free at 4.
    # (This is the two-terminal primitive, not a three-vertex pattern; the
    # triple prober is exercised on W7 below.)
    diamond_edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3)]
    assert forced_equal(diamond_edges, 4, 2, 3, 3) is True
    assert forced_equal(diamond_edges, 4, 2, 3, 4) is False

    # W7 (hub + regular hexagon): the alternating rim triple is forced at k=3,
    # but not at k=4 or k=5.
    wheel_edges = [(0, i) for i in range(1, 7)]
    wheel_edges += [(i, i + 1) for i in range(1, 6)] + [(6, 1)]
    assert forced_mono(wheel_edges, 7, (1, 3, 5), 3) is True
    assert forced_mono(wheel_edges, 7, (1, 3, 5), 4) is False
    assert forced_mono(wheel_edges, 7, (1, 3, 5), 5) is False
    assert forces_mono_triple_set(wheel_edges, 7, [(1, 3, 5), (2, 4, 6)], 4) is False
    assert forces_mono_triple_set(wheel_edges, 7, [(1, 3, 5), (2, 4, 6)], 5) is False

    report = validate_bruteforce(max_n=5, k=3)
    print("forced_mono brute-force validation")
    print(f"  graphs  = {report['graphs']}")
    print(f"  checks  = {report['checks']}")
    print(f"  mismatches = {report['mismatches']}")
    print(f"  seconds = {report['seconds']:.1f}")
    assert report["mismatches"] == 0, report["mismatch_examples"]
    print("OK")


if __name__ == "__main__":
    _selftest()
