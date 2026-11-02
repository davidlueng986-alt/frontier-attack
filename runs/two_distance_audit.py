"""Independent audit of runs/two_distance.json (does not import the module)."""
import itertools, json, math
import sympy as sp
from pysat.solvers import Cadical195

D = json.load(open(r"runs/two_distance.json", encoding="utf-8"))
n = D["n"]
pts = [(sp.sympify(x), sp.sympify(y)) for x, y in D["points"]]
print("source:", D["source"][:110], "...")
print("n,m,chi:", n, D["m"], D["chi"], "d2:", D["d2_values"])

assert n == 33 and D["m"] == 113 and D["chi"] == 6

# exact pairwise distances
cls = {}
for i, j in itertools.combinations(range(n), 2):
    d2 = sp.nsimplify(sp.simplify(sp.expand((pts[i][0]-pts[j][0])**2 + (pts[i][1]-pts[j][1])**2)))
    cls[(i, j)] = d2
unit = {e for e, d in cls.items() if d == 1}
d3 = {e for e, d in cls.items() if d == 3}
E = {tuple(e) for e in D["edges"]}
NU = {tuple(e) for e in D["non_unit_edges"]}
print("unit", len(unit), "sqrt3", len(d3), "edges", len(E), "non_unit", len(NU))
assert E == unit | d3
assert NU == d3
assert len(E) == 113 and len(NU) == 32

# points match the published Table 1 encoding
T = D["points_table1"]
assert len(T) == 33
for (a0, a1, a2, a3), (x, y) in zip(T, pts):
    assert sp.simplify(x - (a0 + a1*sp.sqrt(5))/8) == 0
    assert sp.simplify(y - (a2*sp.sqrt(3) + a3*sp.sqrt(15))/8) == 0

# float coordinates: edge-length error
worst = 0.0
for i, j in E:
    px, py = D["points_float"][i]
    qx, qy = D["points_float"][j]
    d = math.hypot(px-qx, py-qy)
    target = 1.0 if (i, j) in unit else math.sqrt(3)
    worst = max(worst, abs(d-target))
print("max |float edge length - target| =", worst)
assert worst < 1e-9

# colouring
col = D["colouring"]
assert len(col) == n and set(col) <= set(range(6))
assert all(col[i] != col[j] for i, j in E)
print("6-colouring valid; colour classes:", {c: col.count(c) for c in sorted(set(col))})

# 5-colouring UNSAT, independent of the module
k = 5
var = lambda v, c: v*k + c + 1
s = Cadical195()
for v in range(n):
    s.add_clause([var(v, c) for c in range(k)])
    for c1 in range(k):
        for c2 in range(c1+1, k):
            s.add_clause([-var(v, c1), -var(v, c2)])
for i, j in E:
    for c in range(k):
        s.add_clause([-var(i, c), -var(j, c)])
print("independent 5-colouring:", "SAT" if s.solve() else "UNSAT")
s.delete()

# connectivity / degrees
deg = [0]*n
for i, j in E:
    deg[i] += 1; deg[j] += 1
seen = {0}; stack = [0]
adj = {i: set() for i in range(n)}
for i, j in E:
    adj[i].add(j); adj[j].add(i)
while stack:
    v = stack.pop()
    for w in adj[v]:
        if w not in seen:
            seen.add(w); stack.append(w)
print("connected:", len(seen) == n, "min/max degree:", min(deg), max(deg))
print("AUDIT OK")
