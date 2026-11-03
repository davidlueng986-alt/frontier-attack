# Closed attempt: affine-plane projection of the supplied three-dimensional collision source

2026-11-03. **Planar Jacobian conjecture remains OPEN.**

This is a counterexample-construction calibration with a complete negative answer for its
specified family. It is not new main-goal progress: affine domain restriction and linear output
projection leave degree at most7, already within the established planar low-degree exclusion range.
That theorem should have screened out this family before symbolic search. The artifacts are retained
as exact verification controls, not a promising new research lead or a degree-bound improvement.

## Source and candidate family

The polynomial source is exactly `alpoge_counterexample()` in `problems/jacobian/jac.py`.
The present run independently expanded its determinant as -2 and verified the listed pairwise
collisions. This source is three-dimensional and is not itself a planar counterexample.
The historical attribution is retained in that source file; no new authorship claim is made.

For each of its three listed collision pairs p,q, let d=q-p. All have d_x≠0.
Every affine plane through that collision line has, after an invertible change of plane parameters,
one of the forms:

    H_t(u,v)=p+u*d+v*(0,1,t), t∈C,
    H_infinity(u,v)=p+u*d+v*(0,0,1).

Indeed subtracting a multiple of d from the second plane direction makes its x-coordinate zero;
its y-coordinate is either nonzero (scale to1) or zero (scale the z-coordinate to1).
Changing plane coordinates only composes with an affine automorphism and rescales the Jacobian,
so cannot hide a nonzero constant-Jacobian projection.

Consider **every** constant linear map L:C³→C² and F=L∘source∘H.
It is polynomial and automatically retains the collision H(0,0)=p, H(1,0)=q.
A nonzero constant determinant would therefore immediately be a genuine planar counterexample.

## Exact exhaustive certificate for the family

Write the three 2×2 minors of D(source∘H) as M12,M13,M23.
For any linear output projection, Cauchy–Binet gives

    det DF = a*M12+b*M13+c*M23

for three constants a,b,c. If det DF is constant, every nonconstant monomial coefficient must vanish.
Let B(t) be the matrix whose columns are those coefficients in the three minors. Thus B(t)*(a,b,c)=0.

For each finite chart, the artifact supplies 3×3 minors D_i(t) of B(t) and rational polynomials A_i(t)
with the exact identity Σ A_i(t)D_i(t)=1. Hence B(t) has column rank3 for **every complex t**,
not just every sampled rational value. Thus a=b=c=0 and det DF=0, contradicting a Keller projection.
The infinity chart has an explicitly nonzero rational 3×3 minor.

| Source pair | Finite-chart nonzero Bezout terms | Infinity-chart minor |
|---|---:|---:|
| 0,1 | 1 | 189/2 |
| 0,2 | 1 | 189/2 |
| 1,2 | 2 | -84375/2 |

For the first pair the coefficients at (u^0 v^1),(u^0 v^3),(u^1 v^1) already have determinant59808,
independent of t. The certificate covers all three pairs and both projective charts.
No output projection of this affine family can be a nonzero constant-Jacobian map.

## Verification and controls

- Producer differentiates the composed polynomials.
- Independent verifier imports no producer code; it reconstructs D(source)(H)*DH by the chain rule,
  recomputes the named coefficients/minors, and checks the rational polynomial identities exactly.
- It also checks the source determinant and the exact collisions.
- Corrupt Bezout coefficient, incorrect infinity determinant, and a constant-monomial substitution
  are all rejected.
- The positive-control triple (u,v,u+v) has a valid identity projection and is not excluded by this
  obstruction test.

```bash
python3 attempts/jc_affine_projection.py
python3 attempts/verify_jc_affine_projection.py
```

Artifact: `runs/jc_affine_projection_certificate.json`. Standard Python plus SymPy only.

## Decision

Stop affine-plane/linear-projection sampling for these three source collision pairs. Nonlinear
embeddings and other pairs are not covered by this computation, but any new family must first
survive known degree bounds. Quadratic embeddings give degree at most14 and are also excluded
as new research leads. Higher degree alone is not sufficient: require a genuine curved polynomial
image and the full determinant identity. An unrestricted proof or counterexample has not been obtained.
