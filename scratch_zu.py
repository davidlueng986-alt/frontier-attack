"""Scratch: verify (z,u) reduction and probe the degree-(7,10) leading block."""
import sys, time, json
from pathlib import Path
import sympy as sp

sys.path.insert(0, str(Path(r"C:\Users\daive\Downloads\frontier-attack\problems\jacobian")))
from config2_system import lattice_points, P_CORNERS, Q_CORNERS

x, y, z, u = sp.symbols("x y z u")

SP = lattice_points(P_CORNERS)
SQ = lattice_points(Q_CORNERS)
print("supports", len(SP), len(SQ))

# ---- build P,Q in (z,u) coordinates with the exact slice structure ----
def slice_of(pt):
    a, b = pt
    return a, b - 2 * a          # (z-exponent, u-exponent)

from collections import defaultdict
buckets_P = defaultdict(list)
buckets_Q = defaultdict(list)
for p in SP:
    a, c = slice_of(p)
    buckets_P[c].append(a)
for q in SQ:
    a, c = slice_of(q)
    buckets_Q[c].append(a)
print("P u-slices:", {k: (min(v), max(v), len(v)) for k, v in sorted(buckets_P.items())})
print("Q u-slices:", {k: (min(v), max(v), len(v)) for k, v in sorted(buckets_Q.items())})

# named polynomials A=P_{-2}, C=P_{-1}, F=P_0 ; B=Q_{-3}, D=Q_{-2}, E=Q_{-1}, G=Q_0
def poly(name, exps):
    c = {e: sp.Symbol("%s%d" % (name, e)) for e in sorted(exps)}
    return sum(v * z ** e for e, v in c.items()), c

A, cA = poly("A", buckets_P[-2])
C, cC = poly("C", buckets_P[-1])
F, cF = poly("F", buckets_P[0])
B, cB = poly("B", buckets_Q[-3])
D, cD = poly("D", buckets_Q[-2])
E, cE = poly("E", buckets_Q[-1])
G, cG = poly("G", buckets_Q[0])
print("A terms", sp.Poly(A, z).monoms()[:3], "deg", sp.degree(A, z))
print("B terms", sp.Poly(B, z).monoms()[:3], "deg", sp.degree(B, z))

Pz = F + u**-1 * C + u**-2 * A
Qz = G + u**-1 * E + u**-2 * D + u**-3 * B

# exact identity: u^6 * Jac_{z,u}(P,Q) = z^2  <=>  det J = x^2
Jac = sp.expand(sp.diff(Pz, z) * sp.diff(Qz, u) - sp.diff(Pz, u) * sp.diff(Qz, z))
expr = sp.expand(u ** 6 * Jac)
print("u-powers present:", sorted({sp.degree(t, u) for t in sp.Add.make_args(expr)}))
E0 = sp.expand(-3 * sp.diff(A, z) * B + 2 * A * sp.diff(B, z) - z**2)
E1 = sp.expand(-2 * sp.diff(A, z) * D - 3 * sp.diff(C, z) * B + 2 * A * sp.diff(D, z) + C * sp.diff(B, z))
E2 = sp.expand(-sp.diff(A, z) * E - 2 * sp.diff(C, z) * D - 3 * sp.diff(F, z) * B
               + 2 * A * sp.diff(E, z) + C * sp.diff(D, z))
E3 = sp.expand(-sp.diff(C, z) * E - 2 * sp.diff(F, z) * D + C * sp.diff(E, z) + 2 * A * sp.diff(G, z))
E4 = sp.expand(-sp.diff(F, z) * E + C * sp.diff(G, z))
rebuilt = sp.expand(E0 + u * E1 + u**2 * E2 + u**3 * E3 + u**4 * E4)
print("identity u^6 Jac == sum u^k E_k :", sp.simplify(expr - rebuilt) == 0)

# ---- leading block (G0): A = z*A1, B = z^2*B1 ----
a = sp.symbols("a0:8")     # A1 coefficients z^0..z^7
b = sp.symbols("b0:11")    # B1 coefficients z^0..z^10
A1 = sum(a[i] * z**i for i in range(8))
B1 = sum(b[i] * z**i for i in range(11))
G0 = sp.expand(A1 * B1 - 3 * z * sp.diff(A1, z) * B1 + 2 * z * A1 * sp.diff(B1, z) - 1)
polyG0 = sp.Poly(G0, z)
coeffsG0 = [sp.expand(c) for c in polyG0.all_coeffs()[::-1]]
print("G0 degree", polyG0.degree(), "ncoeffs", len(coeffsG0))
for i, c in enumerate(coeffsG0):
    print("  z^%d: %s" % (i, c))
