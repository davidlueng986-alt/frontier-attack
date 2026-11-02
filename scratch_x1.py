"""Scratch 4: x^1 block dimension, G0 block dimension, edge leading-form systems."""
import sys, time, itertools
from pathlib import Path
import sympy as sp

sys.path.insert(0, str(Path(r"C:\Users\daive\Downloads\frontier-attack\problems\jacobian")))
from config2_system import lattice_points, P_CORNERS, Q_CORNERS

x, y, z = sp.symbols("x y z")
SP = lattice_points(P_CORNERS); SQ = lattice_points(Q_CORNERS)

def build(supp, pre):
    cs = {}
    expr = 0
    for p in supp:
        c = sp.Symbol("%s%d_%d" % (pre, p[0], p[1]))
        cs[p] = c
        expr += c * x ** p[0] * y ** p[1]
    return expr, cs

P, cP = build(SP, "a"); Q, cQ = build(SQ, "b")
det = sp.expand(sp.diff(P, x) * sp.diff(Q, y) - sp.diff(P, y) * sp.diff(Q, x))
dp = sp.Poly(det, x, y)

print("=== x^1 block ===")
blk = []
for t in range(0, 8):
    if dp.coeff_monomial((1, t)) != 0 or (1, t) in [(m[0], m[1]) for m in dp.monoms()]:
        c = sp.expand(dp.coeff_monomial((1, t)))
        if c != 0:
            print("  y^%d: %s" % (t, c))
            blk.append(c)
x2 = sp.expand(dp.coeff_monomial((2, 0)))
print("  x^2 coeff:", x2)

# dimension of the x^1 block alone, and with x^2=1
syms = [sp.Symbol("a1_0"), sp.Symbol("a1_1"), sp.Symbol("a1_2"),
        sp.Symbol("b1_1"), sp.Symbol("b1_2")]
t0 = time.time()
G = sp.groebner(blk, *syms, order="grevlex")
print("  GB(x1 block):", G)
print("  seconds %.2f" % (time.time() - t0))
# dimension via the leading monomials
lm = [sp.Poly(g, *syms, order="grevlex").LM() for g in G.exprs]
print("  leading monomials:", lm)

t0 = time.time()
G2 = sp.groebner(blk + [x2 - 1], *syms, order="grevlex")
print("  GB(x1 block + x2=1) in 6 vars:", G2)
print("  seconds %.2f" % (time.time() - t0))
