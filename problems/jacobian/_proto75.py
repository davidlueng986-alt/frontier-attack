import sys, time
from collections import defaultdict
sys.path.insert(0, 'problems/jacobian')
import sympy as sp
import config2_system as c

P = ((0, 0), (9, 0), (12, 3), (0, 15))
Q = ((0, 0), (15, 0), (20, 5), (0, 25))


def poly_terms(support, prefix):
    return [(a, b, prefix, i) for i, (a, b) in enumerate(support)]


def deriv_terms(poly, axis):
    out = []
    for a, b, sym, i in poly:
        if axis == 0 and a > 0:
            out.append((a - 1, b, a, sym, i))
        elif axis == 1 and b > 0:
            out.append((a, b - 1, b, sym, i))
    return out


def fast_det(sp, sq):
    pt = poly_terms(sp, 'a')
    qt = poly_terms(sq, 'b')
    px, py = deriv_terms(pt, 0), deriv_terms(pt, 1)
    qx, qy = deriv_terms(qt, 0), deriv_terms(qt, 1)
    det = defaultdict(lambda: defaultdict(int))
    for ax, ay, ca, _, ia in px:
        for bx, by, cb, _, ib in qy:
            det[(ax + bx, ay + by)][(ia, ib)] += ca * cb
    for ax, ay, ca, _, ia in py:
        for bx, by, cb, _, ib in qx:
            det[(ax + bx, ay + by)][(ia, ib)] -= ca * cb
    out = {}
    for k, v in det.items():
        vv = {kk: vvv for kk, vvv in v.items() if vvv}
        if vv:
            out[k] = vv
    return out, len(px) * len(qy) + len(py) * len(qx)


def fmt_eq(eq, p, q):
    return ' + '.join(
        f'{co}*a{i}*b{j}' if co != 1 else f'a{i}*b{j}'
        for (i, j), co in sorted(eq.items()))


for label, sp_, sq_ in [('vtx', list(P), list(Q)), ('full', c.lattice_points(P), c.lattice_points(Q))]:
    t = time.perf_counter()
    d, pairs = fast_det(sp_, sq_)
    print(label, 'points', len(sp_), len(sq_), 'pairs', pairs, 'monoms', len(d),
          'maxdeg', max(a + b for a, b in d), 'time', round(time.perf_counter() - t, 3))
    print('  c00', d.get((0, 0)), 'c20', d.get((2, 0)))
    print('  eq00', fmt_eq(d.get((0, 0), {}), sp_, sq_) if (0, 0) in d else 'ABSENT')
    # compare vertices to sympy
    if label == 'vtx':
        x, y = sp.symbols('x y')
        r = c.size_system(sp_, sq_, x**2, 'v', x, y)
        print('  sympy monoms', r['equation_monomials'])
        print('  sympy sample', r['equations_sample'][0])
        print('  fast sample', fmt_eq(d[(28, 4)], sp_, sq_))
