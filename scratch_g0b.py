"""Scratch 6: G0 residual (6 eqs in a1..a7): dimension and structure."""
import sympy as sp, time
z = sp.Symbol("z")
a = sp.symbols("a1:8")
b = sp.symbols("b1:11")
A1 = 1 + sum(a[i] * z ** (i + 1) for i in range(7))
B1 = 1 + sum(b[i] * z ** (i + 1) for i in range(10))
G0 = sp.Poly(sp.expand(A1 * B1 - 3 * z * sp.diff(A1, z) * B1
                       + 2 * z * A1 * sp.diff(B1, z) - 1), z)
cs = G0.all_coeffs()[::-1]

# solve b1..b10 linearly from z^1..z^10
sol = {}
eqs = [sp.expand(c) for c in cs[1:11]]
vs = list(b)
M, rhs = sp.linear_eq_to_matrix(eqs, vs)
Ms = sp.Matrix(M)
print("matrix det (should be nonzero):", sp.factor(Ms.det()))
bsol = Ms.LUsolve(-rhs)
subs = dict(zip(vs, list(bsol)))
res = [sp.expand(c.subs(subs)) for c in cs[11:17]]
print("residual equations in a1..a7:")
for r in res:
    print("  ", sp.factor(r))

t0 = time.time()
G = sp.groebner(res, *a, order="grevlex")
print("GB seconds %.2f" % (time.time() - t0))
print("GB:", G)
