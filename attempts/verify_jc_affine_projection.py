"""Independent chain-rule reconstruction and exact Bezout/minor proof replay."""
import copy
import json
import sys
from pathlib import Path
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'problems/jacobian'))
from jac import alpoge_counterexample


def check_minor_identity(columns, spec, u, v, t):
    pairs = ((0, 1), (0, 2), (1, 2))
    polys = [s.Poly(s.expand(columns[a, 0]*columns[b, 1]-columns[a, 1]*columns[b, 0]), u, v)
             for a, b in pairs]

    def minor(monomials):
        assert len(monomials) == 3 and len(set(map(tuple, monomials))) == 3
        assert all(len(m) == 2 and all(type(i) is int and i >= 0 for i in m)
                   and tuple(m) != (0, 0) for m in monomials)
        return s.expand(s.Matrix([[p.coeff_monomial(u**a*v**b) for p in polys]
                                  for a, b in monomials]).det())

    def poly(coeffs):
        assert isinstance(coeffs, list) and coeffs
        return sum(s.Rational(a)*t**i for i, a in enumerate(coeffs))

    if isinstance(spec, list):
        total = 0
        for item in spec:
            determinant = minor(item['monomials'])
            assert s.expand(determinant-poly(item['determinant'])) == 0
            total += determinant*poly(item['multiplier'])
        assert s.expand(total-1) == 0
    else:
        determinant = minor(spec['monomials'])
        assert determinant == s.Rational(spec['determinant']) != 0


def verify(report):
    assert report['schema'] == 'jc-three-source-pairs-v1' and report['planar_JC'] == 'OPEN'
    assert sorted(tuple(c['collision_pair']) for c in report['cases']) == [(0, 1), (0, 2), (1, 2)]
    assert len(report['cases']) == 3
    F, xyz, points = alpoge_counterexample()
    J = s.Matrix(F).jacobian(xyz)
    assert s.expand(J.det()) == -2
    u, v, t = s.symbols('u v t')
    for case in report['cases']:
        p, q = [s.Matrix(points[i]) for i in case['collision_pair']]
        assert case['collision_points'] == [[str(x) for x in point] for point in (p, q)]
        assert p != q and p[0] != q[0]
        assert all(s.expand(f.subs(dict(zip(xyz, p)))-f.subs(dict(zip(xyz, q)))) == 0 for f in F)
        assert case['finite_chart_direction'] == ['0', '1', 't']
        assert case['infinity_chart_direction'] == ['0', '0', '1']
        direction = q-p
        for name, w in [('finite', s.Matrix([0, 1, t])), ('infinity', s.Matrix([0, 0, 1]))]:
            H = p + u*direction + v*w
            # Producer differentiated composed polynomials; verifier instead uses the chain rule.
            columns = J.subs(dict(zip(xyz, H)), simultaneous=True)*s.Matrix.hstack(direction, w)
            spec = case['finite_bezout'] if name == 'finite' else case['infinity']
            check_minor_identity(columns, spec, u, v, t)
    return True


def negative_controls(report):
    cases = []
    bad = copy.deepcopy(report)
    bad['cases'][0]['finite_bezout'][0]['multiplier'][0] = '0'
    cases.append(bad)
    bad = copy.deepcopy(report)
    bad['cases'][0]['infinity']['determinant'] = '1'
    cases.append(bad)
    bad = copy.deepcopy(report)
    bad['cases'][0]['finite_bezout'][0]['monomials'][0] = [0, 0]
    cases.append(bad)
    for bad in cases:
        try:
            verify(bad)
        except AssertionError:
            continue
        raise AssertionError('Corrupted certificate accepted')
    # Positive-control map (u,v,u+v) DOES admit a constant-Jacobian linear projection.
    # It must not be excluded by the nonconstant-coefficient obstruction certificate.
    u, v, t = s.symbols('u v t')
    columns = s.Matrix([[1, 0], [0, 1], [1, 1]])
    try:
        check_minor_identity(columns, report['cases'][0]['finite_bezout'], u, v, t)
    except AssertionError:
        return True
    raise AssertionError('Valid projection positive control was wrongly excluded')


if __name__ == '__main__':
    report = json.loads((ROOT/'runs/jc_affine_projection_certificate.json').read_text())
    assert verify(report) and negative_controls(report)
    print('EXACT_ALL_THREE_PAIRS_ALL_AFFINE_PLANES_ALL_LINEAR_PROJECTIONS_PASS')
    print('3 certificate mutations rejected; valid projection positive control preserved; planar JC OPEN')
