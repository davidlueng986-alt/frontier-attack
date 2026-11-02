"""Scratch 3: exact det-J support + G0 (top weight block) solution dimension."""
import sys, time
from pathlib import Path
import sympy as sp

sys.path.insert(0, str(Path(r"C:\Users\daive\Downloads\frontier-attack\problems\jacobian")))
from config2_system import lattice_points, P_CORNERS, Q_CORNERS

z = sp.Symbol("z")
x, y = sp.symbols("x y")

SP = lattice_points(P_CORNERS); SQ = lattice_points(Q_CORNERS)

# ---------------- exact attainable support of P_xQ_y - P_yQ_x ----------------
def dsup(supp, axis):
    if axis == 0:
        return {(a - 1, b): (a, b) for a, b in supp if a > 0}
    return {(a, b - 1): (a, b) for a, b in supp if b > 0}

Px, Py = dsup(SP, 0), dsup(SP, 1)
Qx, Qy = dsup(SQ, 0), dsup(SQ, 1)

pairs = {}   # monomial -> list of (which, p, q)
for p in Px:
    for q in Qy:
        m = (p[0] + q[0], p[1] + q[1])
        pairs.setdefault(m, []).append(("PxQy", p, q, +1))
for p in Py:
    for q in Qx:
        m = (p[0] + q[0], p[1] + q[1])
        pairs.setdefault(m, []).append(("PyQx", p, q, -1))
print("attainable monomials:", len(pairs))
print("(2,0) attainable:", (2, 0) in pairs)
print("terms at (2,0):", pairs.get((2, 0)))
print("min x-degree:", min(m[0] for m in pairs))
print("x=1 monomials:", sorted([m for m in pairs if m[0] == 1]))
mindeg = min(m[0] for m in pairs)

# symbolic check on the full 72-coefficient system
def build(supp, pre):
    cs = {}
    expr = 0
    for p in supp:
        c = sp.Symbol("%s%d_%d" % (pre, p[0], p[1]))
        cs[p] = c
        expr += c * x ** p[0] * y ** p[1]
    return expr, cs

P, cP = build(SP, "a"); Q, cQ = build(SQ, "b")
t0 = time.time()
det = sp.expand(sp.diff(P, x) * sp.diff(Q, y) - sp.diff(P, y) * sp.diff(Q, x))
print("expand seconds %.2f" % (time.time() - t0))
dp = sp.Poly(det, x, y)
mons = sorted((int(a), int(b)) for a, b in dp.monoms())
print("det monomials:", len(mons), "min xdeg:", min(m[0] for m in mons),
      "max deg:", max(a + b for a, b in mons))
print("(2,0) in det:", (2, 0) in mons)
c20 = sp.expand(dp.coeff_monomial((2, 0)))
print("coeff of x^2:", c20)
print("support == attainable set:", set(mons) == set(pairs))
