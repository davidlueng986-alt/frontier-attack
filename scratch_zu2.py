"""Scratch 2: debug the u-expansion identity."""
import sys
from pathlib import Path
import sympy as sp
from collections import defaultdict

sys.path.insert(0, str(Path(r"C:\Users\daive\Downloads\frontier-attack\problems\jacobian")))
from config2_system import lattice_points, P_CORNERS, Q_CORNERS

x, y, z, u = sp.symbols("x y z u")

SP = lattice_points(P_CORNERS); SQ = lattice_points(Q_CORNERS)

def slices(supp):
    d = defaultdict(list)
    for a, b in supp:
        d[b - 2 * a].append(a)
    return d

buP, buQ = slices(SP), slices(SQ)

def poly(name, exps):
    return sum(sp.Symbol("%s%d" % (name, e)) * z ** e for e in sorted(exps))

A = poly("A", buP[-2]); C = poly("C", buP[-1]); F = poly("F", buP[0])
B = poly("B", buQ[-3]); D = poly("D", buQ[-2]); E = poly("E", buQ[-1]); G = poly("G", buQ[0])

Pz = F + u**-1 * C + u**-2 * A
Qz = G + u**-1 * E + u**-2 * D + u**-3 * B
Jac = sp.expand(sp.diff(Pz, z) * sp.diff(Qz, u) - sp.diff(Pz, u) * sp.diff(Qz, z))
expr = sp.expand(u ** 6 * Jac)
p = sp.Poly(expr, u)
print("u-degrees:", p.degree())
for (k,), co in sorted(p.terms()):
    print("u^%d: %s" % (k, sp.expand(co)))
print()
print("expected E0..E4:")
print("E0", sp.expand(-3 * sp.diff(A, z) * B + 2 * A * sp.diff(B, z) - z**2))
print("E1", sp.expand(-2 * sp.diff(A, z) * D - 3 * sp.diff(C, z) * B + 2 * A * sp.diff(D, z) + C * sp.diff(B, z)))
print("E2", sp.expand(-sp.diff(A, z) * E - 2 * sp.diff(C, z) * D - 3 * sp.diff(F, z) * B + 2 * A * sp.diff(E, z) + C * sp.diff(D, z)))
print("E3", sp.expand(-sp.diff(C, z) * E - 2 * sp.diff(F, z) * D + C * sp.diff(E, z) + 2 * A * sp.diff(G, z)))
print("E4", sp.expand(-sp.diff(F, z) * E + C * sp.diff(G, z)))
