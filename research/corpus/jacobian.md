# Planar Jacobian: annotated corpus

The terminal object is a polynomial pair `(P,Q)` with `J(P,Q) ≡ 1`, `(deg P, deg Q) = (72,108)`
(or a later pair), and two distinct points with the same image. The corpus has three jobs:
fix the degree ladder, fix the Newton polygons, and separate the 2026 dimension-3
counterexample from the still-open plane case.

## 1. The degree ladder

| bound | statement | source |
|---|---|---|
| ≤ 2 | every Keller map of degree ≤ 2 is invertible, all n | `doi:10.1016/0021-8693(80)90233-1` (Wang 1980) |
| ≤ 3 | complete classification of planar Keller pairs, up to linear conjugacy | `doi:10.15393/j3.art.2016.3510` (Starkov 2016) |
| ≤ 4 | complete classification; new quartic family `F_{c,s}` | `doi:10.15393/j3.art.2017.3870` (Starkov 2017) |
| ≤ 100 | no planar counterexample of total degree ≤ 100 | `doi:10.1515/crll.1983.340.140` (Moh 1983) |
| ≤ 104 | complex JC(2) proved through degree 104 | `arXiv:1902.05923` (Nguyen) |
| < 125 | only `(72,108)` and `(108,72)` survive | `arXiv:2204.14178` (GGHV 2022, unrefereed) |

**Read the last row carefully.** GGHV's abstract says the lower bound rises "from 100 to
108", not to 125. Theorem 2.1 says every counterexample has `max{deg P, deg Q} ≥ 125`, **or**
one of the two pairs `(72,108)`, `(108,72)`. The pair is not excluded; it is the single
surviving hole. The bound is 108 because `max(72,108) = 108`.

`arXiv:1708.07936` (GGHV 2017) is the enumerator behind the table: it computes all
admissible corners for `gcd(deg P, deg Q) ≤ 35` and all degree pairs with
`max(deg P, deg Q) ≤ 150`. It is the only published enumerator of this kind; the paper
describes a C++/PostgreSQL pipeline but releases no code.

## 2. Newton polygons and the two configurations

**`abhyankar:tifr57`.** Theorem (19.4) is the dictionary the whole search rests on. For a
Keller pair `(f,g)` in `A = k[x,y]`, the following are equivalent:

1. `k[f,g] = A`;
2. `f` has only one point at infinity;
3. `N(f)` is a triangle with vertices `(0,n)`, `(0,0)`, `(m,0)`;
4. `deg f | deg g` or `deg g | deg f`.

Theorem 18.13 gives the general constraint: each component of a Keller pair has at most two
points at infinity. A counterexample must therefore sit in the narrow band where both
components have exactly two, and neither Newton polygon is a triangle.

**`arXiv:2204.14178`, Proposition 4.3 (p. 23).** For `(72,108)`, after the GGHV
normalisation, exactly two corner sets survive:

```
config 1 (may be excluded by the no-bisector condition):
  N(P) = {(0,0),(1,0),(8,14),(8,16),(0,8)}
  N(Q) = {(0,0),(2,1),(12,21),(12,24),(0,12)}

config 2 (no slope-1 edge):
  N(P) = {(0,0),(1,0),(8,14),(8,16)}
  N(Q) = {(0,0),(2,1),(12,21),(12,24)}
```

These are the **reduced** polygons: `deg P = 24`, `deg Q = 36`, scaled by `(m,n) = (3,2)`
to reach `(72,108)`. The normal form is `[P,Q] = x²`, so the system to solve is
`det J(P,Q) = x²`, not `det J = 1`.

**The reconstruction gap.** The GGHV object lives in `L^(1) = K[x, x⁻¹, y]`, and the
reduction from a polynomial counterexample uses automorphisms of `L^(1)` that are not
automorphisms of `K[x,y]`. The implication is one-directional:

```
polynomial counterexample  ⇒  normalised pair with [P,Q] = x²
```

It does not reverse. A solution of the config-2 system is not a polynomial map with
`det J = 1`, and `jac.py` rejects it for exactly that reason. The corpus records this as a
negative result about the current pipeline, not as a defect of GGHV.

## 3. Methods and reductions

**`doi:10.1007/BF01459126` (Drużkowski).** Reduces JC(n) to cubic-linear maps
`F = x + (Ax)^{*3}` with `A` nilpotent. This is the standard finite-dimensional reduction:
it turns a statement about all degrees into a statement about one cubic homogeneous map.

**`doi:10.1090/s0002-9939-05-07570-2` (de Bondt–van den Essen).** Reduces JC to the case
where the Jacobian matrix is symmetric. The coefficient space of a search shrinks; the
degree does not.

**`doi:10.13140/2.1.3906.3844` (Nowicki).** The derivation-side formulation: a polynomial
map is a Keller map iff the associated derivations have the right rings of constants. The
Nowicki conjecture on Weitzenböck derivations lives here.

**`doi:10.1007/978-3-030-60535-3` (van den Essen, 2nd ed. 2021).** The standard monograph.
Use it for normal forms and the equivalence list; do not use it as a source of a degree
bound (Moh and Nguyen are the primary sources).

## 4. The July 2026 counterexample in dimension 3

**`zenodo:21461572`.** An independent verification and geometric analysis of the July 2026
map `F : ℂ³ → ℂ³`. It records the constant determinant and the three-point collision, and
contains a Python verifier.

**`tao:jacobian-2026`.** The cleanest derivation: component degrees `(7,6,4)`, total degree
7, `det JF ≡ −2`, generic fibre of 3 points. A linear change of one component normalises
the determinant to 1, so this is a genuine Keller map with a collision — in dimension 3.

**`speyer:jacobian-2026`.** The algebraic-geometry reading: the map is built from a
resultant/degree-7 construction; the collision is explicit and the monodromy is `S₃`. The
comments also record a separate 2026 cubic-homogeneous map in 24 variables with
`det(DG) = 1` and a rational collision.

**What this does and does not settle.** It settles JC(3) negatively. It does **not** settle
JC(2), and it does not produce a planar pair. The repo's `accept.py` correctly rejects the
C³ map on the degree filter: its degree is 7, and the planar bound is 125. The plane
problem is a separate open problem, not a corollary of the dimension-3 counterexample.

## 5. The (72,108) exclusion claim

**`github:wstrinz/plane-jacobian-72-108`.** Theorem A states: for a field `K` of
characteristic 0, there is no pair `(P,Q)` in `K[x,x⁻¹,y]²` with `[P,Q] = x²` and either
configuration (1) or configuration (2) of GGHV 2022, Proposition 4.3. Corollary B,
conditional on two stated inputs, gives `max{deg P, deg Q} ≥ 125` for every plane
counterexample.

The README also states that Billel Helali published an independent exclusion on 21 July
2026 (Zenodo DOI `10.5281/zenodo.21479814`) and that the repo's audit adjudicates Helali's
chain as subsuming its own. The claim is unrefereed. The repo's own audit
(`research/jacobian-wstrinz-audit.md`) found the argument substantive but noted internal
inconsistencies: `CURRENT_STATUS.md` still says `OPEN`, two cited close-out documents are
missing, and the stated proof hash is absent from the public history.

**How to use it.** Treat Theorem A as a strong prior, not a theorem. If it holds, the
`(72,108)` target is dead and the next region is the 24 degree pairs with
`max_deg ≥ 125` listed by the repo's J4 pass; the smallest is `(75,125)`. If it fails, the
counterexample is in the two configurations above and the coefficient system is the place
to look. Either way, the reconstruction gap must be closed before a config-2 solution can
become a claim.

## 6. Analogues

**`doi:10.24033/bsmf.1667` (Dixmier 1968).** The Weyl-algebra endomorphism problem.
`Dixmier(n)` is equivalent to `JC(n)` (Tsuchimoto; Belov-Kanel–Kontsevich). It is open in
every dimension, including n = 2, and is the closest non-commutative sibling.

**`arXiv:2608.02634` (Mondello 2026).** A dimension-2 counterexample to the **separable**
Jacobian conjecture in characteristic 2. Over the algebraic closure of `𝔽₂`:

```
P = x + x²y + x⁴ + x⁶y²
Q = y + x⁵ + x⁶y + x⁷y² + x⁸y³
det J(P,Q) = 1
(0,1), (1,0), (1,1) have the same image
```

The function-field extension has degree 3 and is separable. This is the cleanest "the
analogue is false" data point in the corpus: characteristic 2 breaks the theorem, and the
break is explicit at degree 3.

**`doi:10.1090/s0894-0347-03-00440-5` (Shestakov–Umirbaev 2003).** The tame generators
problem in dimension 3: `Aut(k[x,y,z])` contains wild automorphisms. A structural question
that looked like a stepping stone to JC was settled negatively, by a method (Nagy's
ridigity/degree argument) that does not touch the Jacobian condition.

**`doi:10.1007/978-3-030-42136-6_14` (van den Essen 2020).** Mathieu–Zhao spaces. JC implies
the Mathieu conjecture for the relevant class; the Mathieu–Zhao conjecture is strictly
weaker and still open. This is the widest currently open statement that JC is known to
imply.

**`doi:10.4064/ap-62-1-45-74` (Feßler 1995).** The Markus–Yamabe conjecture in dimension 2:
a local diffeomorphism of `ℝ²` with negative-trace Jacobian everywhere has 0 as a global
attractor. This is the dynamical analogue of JC(2) and it is **true**.

**`doi:10.1006/aima.1997.1673` (Cima–van den Essen–Gasull–Hubbers 1997).** The polynomial
counterexample in dimension 3. The planar theorem does not lift. This pair is the best
" solved analogue with a sharp boundary" story in the corpus: a conjecture true in
dimension 2, false in dimension 3, with an explicit polynomial map as the witness.

## 7. What the corpus implies

1. **The degree ladder is a theorem, not a heuristic, up to 104.** Below 104 the search is
   provably empty. Between 105 and 124 the only surviving pair is `(72,108)`, conditional
   on GGHV 2022.
2. **The reconstruction gap is the actual wall.** A config-2 solution has determinant `x²`
   in `L^(1)`, not a constant in `K[x,y]`. Until the inverse normalisation is written and
   checked, no amount of elimination produces the terminal object.
3. **The dimension-3 counterexample changes the context, not the target.** JC is false for
   n = 3 in this record. The plane case is still open, and the existing degree bounds are
   what keep it small.
4. **The `(72,108)` exclusion is the one claim that could close the branch.** Read
   `github:wstrinz/plane-jacobian-72-108` with `research/jacobian-wstrinz-audit.md` open;
   the executable checkers are the part worth trusting, the prose is not refereed.
