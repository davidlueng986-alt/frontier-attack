"""Scratch 5: G0 block (top weight (2,-1)) solution dimension; edge systems."""
import sys, time
from pathlib import Path
import sympy as sp

z = sp.Symbol("z")
a = sp.symbols("a1:8")     # A1 coeffs z^1..z^7  (a1 = coeff of z^1)
b = sp.symbols("b1:11")    # B1 coeffs
A1 = 1 + sum(a[i] * z ** (i + 1) for i in range(7))
B1 = 1 + sum(b[i] * z ** (i + 1) for i in range(10))
G0 = sp.Poly(sp.expand(A1 * B1 - 3 * z * sp.diff(A1, z) * B1
                       + 2 * z * A1 * sp.diff(B1, z) - 1), z)
cs = G0.all_coeffs()[::-1]
print("num coeffs", len(cs), "degrees of coeffs in unknowns:",
      [sp.Poly(c, *a, *b).total_degree() if c != 0 else 0 for c in cs])
for n, c in enumerate(cs):
    if c != 0:
        print(" z^%d: %s" % (n, c))
