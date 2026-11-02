"""Scratch: build the De Neve 33-vertex distance graph and SAT-check chi."""
import itertools, json, re, time
import sympy as sp
from pysat.solvers import Cadical195

TXT = r"runs/papers/deNeve2025.txt"
lines = open(TXT, encoding="utf-8").read().splitlines()
i1 = next(i for i, l in enumerate(lines) if l.startswith("Table 1."))
i2 = next(i for i, l in enumerate(lines) if l.startswith("Table 2."))
rows = [tuple(map(int, m.groups()))
        for l in lines[i1:i2]
        if (m := re.fullmatch(r"\[(-?\d+) (-?\d+) (-?\d+) (-?\d+)\]", l.strip()))]

s3, s5, s15 = sp.sqrt(3), sp.sqrt(5), sp.sqrt(15)
pts = [((a0 + a1 * s5) / 8, (a2 * s3 + a3 * s15) / 8) for a0, a1, a2, a3 in rows]

def d2(p, q):
    return sp.nsimplify(sp.expand(sum((p[k] - q[k]) ** 2 for k in range(2))))

pairs = list(itertools.combinations(range(len(pts)), 2))
unit = [(i, j) for i, j in pairs if d2(pts[i], pts[j]) == 1]
d3 = [(i, j) for i, j in pairs if d2(pts[i], pts[j]) == 3]
edges = sorted(unit + d3)
print("n", len(pts), "unit", len(unit), "sqrt3", len(d3), "edges", len(edges))

def sat_colourable(k, edges, n):
    solver = Cadical195()
    var = lambda v, c: v * k + c + 1
    for v in range(n):
        solver.add_clause([var(v, c) for c in range(k)])
    for v in range(n):
        for c1 in range(k):
            for c2 in range(c1 + 1, k):
                solver.add_clause([-var(v, c1), -var(v, c2)])
    for i, j in edges:
        for c in range(k):
            solver.add_clause([-var(i, c), -var(j, c)])
    t = time.time()
    ok = solver.solve()
    print(f"k={k}: {'SAT' if ok else 'UNSAT'} ({time.time()-t:.2f}s)")
    if ok:
        model = set(x for x in solver.get_model() if x > 0)
        cols = [next(c for c in range(k) if var(v, c) in model) for v in range(n)]
        assert all(cols[i] != cols[j] for i, j in edges)
        print("  colouring:", cols)
    solver.delete()
    return ok

sat_colourable(6, edges, len(pts))
sat_colourable(5, edges, len(pts))
