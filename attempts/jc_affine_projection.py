"""Exact closed attempt: all affine collision planes and all linear output projections.
A negative result is restricted to this construction family, never planar JC.
"""
import itertools,json,sys
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'problems/jacobian'))
from jac import alpoge_counterexample
u,v,t=s.symbols('u v t')


def matrix(chart, pair=(0,1)):
 F,xyz,points=alpoge_counterexample();p,q=[points[i] for i in pair]
 w=(0,1,t) if chart=='finite' else (0,0,1)
 H=[p[i]+u*(q[i]-p[i])+v*w[i] for i in range(3)]
 G=[s.expand(f.subs(dict(zip(xyz,H)),simultaneous=True)) for f in F]
 polys=[s.Poly(s.expand(s.diff(G[a],u)*s.diff(G[b],v)-s.diff(G[a],v)*s.diff(G[b],u)),u,v)
        for a,b in ((0,1),(0,2),(1,2))]
 terms=sorted(set().union(*(set(P.monoms()) for P in polys))-{(0,0)})
 B=s.Matrix([[P.coeff_monomial(u**a*v**b) for P in polys] for a,b in terms])
 return terms,B


def coefficients(p):
 P=s.Poly(p,t,domain=s.QQ)
 return [str(P.nth(i)) for i in range(max(0,P.degree())+1)]


def generate(pair=(0,1)):
 terms,B=matrix('finite',pair);g=s.Poly(0,t,domain=s.QQ);chosen=[];multipliers=[]
 for ids in itertools.combinations(range(B.rows),3):
  d=s.Poly(B[list(ids),:].det(),t,domain=s.QQ)
  if d.is_zero:continue
  a,b,h=s.gcdex(g,d)
  multipliers=[s.expand(a.as_expr()*x) for x in multipliers]+[b.as_expr()]
  chosen.append((ids,d));g=h
  if g.degree()==0:break
 assert g.as_expr()==1,'No full-rank-for-all-parameters certificate found; do not infer failure.'
 finite=[]
 for (ids,d),mult in zip(chosen,multipliers):
  if mult!=0:finite.append({'monomials':[list(terms[i]) for i in ids],
                          'determinant':coefficients(d),'multiplier':coefficients(mult)})
 keys,A=matrix('infinity',pair);pivots=A.T.rref()[1]
 assert len(pivots)==3
 ids=list(pivots);det=s.expand(A[ids,:].det());assert det!=0
 return {'schema':'jc-collision-plane-projection-v1','source':'problems/jacobian/jac.py:alpoge_counterexample',
         'source_dimension':3,'collision_pair':list(pair),'collision_points':[[str(x) for x in alpoge_counterexample()[2][i]] for i in pair],
         'finite_chart_direction':['0','1','t'],'infinity_chart_direction':['0','0','1'],
         'finite_bezout':finite,'infinity':{'monomials':[list(keys[i]) for i in ids],'determinant':str(det)},
         'conclusion':'NO_NONZERO_CONSTANT_JACOBIAN_LINEAR_PROJECTION_ON_ANY_AFFINE_PLANE_THROUGH_THESE_TWO_POINTS',
         'root_planar_JC':'OPEN','scope':'One supplied 3D source and one collision pair; all complex affine planes and all constant linear output projections'}


if __name__=='__main__':
 reports=[generate(pair) for pair in itertools.combinations(range(3),2)]
 out=ROOT/'runs/jc_affine_projection_certificate.json'
 out.write_text(json.dumps({'schema':'jc-three-source-pairs-v1','cases':reports,'planar_JC':'OPEN'},indent=2)+'\n')
 print(json.dumps([{'pair':r['collision_pair'],'finite_minor_count':len(r['finite_bezout']),
 'infinity_det':r['infinity']['determinant'],'root_planar_JC':r['root_planar_JC']} for r in reports]))
