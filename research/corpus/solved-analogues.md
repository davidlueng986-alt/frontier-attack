# Solved analogues

Problems that are structurally close to the three targets and are **already settled**. Each
entry names the method that settled it and, more importantly, what part of that method does
not transfer. The pattern is deliberately repetitive: every genuine solution below ends in
a finite, independently checkable certificate.

| # | problem | settled | method | relation to our targets |
|---:|---|---|---|---|
| 1 | Falconer distance conjecture, all d | 2026 | harmonic analysis / projection theory | distance sets in ℝ^d; same objects, different question |
| 2 | Kakeya set conjecture in ℝ³ | 2025 | tube volume estimates | incidence geometry, no colouring |
| 3 | separable JC in characteristic 2, dimension 2 | 2026 | explicit map + function-field degree | the Jacobian condition, broken by the field |
| 4 | tame generators problem, dimension 3 | 2003 | rigidity / degree estimates | automorphism structure, not the Jacobian |
| 5 | Markus–Yamabe, dimension 2 | 1995 | planar dynamics / divergence | the dynamical analogue of JC(2) |
| 6 | Sato–Tate for elliptic curves over Q | 2011 | potential automorphy | local statistics of one curve |
| 7 | Fermat's Last Theorem | 1995 | modularity lifting | "big theorem, complete method" |
| 8 | Gauss class number problem | 1976–1986 | analytic rank ≥ 3 + Heegner points | BSD-type input settles a classical problem |
| 9 | congruent numbers (families) | 2014 | Heegner points on `m y² = x³ − x` | conditional on BSD in general |
| 10 | JC(2) degree ≤ 104 | 1983–2019 | Newton polygon / root configurations | the degree ladder itself |
| 11 | minimum 6-chromatic 3-distance graph | 2026 | exact arithmetic + SAT + catalog filter | finite combinatorial certificate |
| 12 | χ(ℝ²) ≥ 5 | 2018 | gadget composition + SAT | the HN method that actually worked |

---

## 1. Falconer distance conjecture (2026)

**Solved.** For every integer `d ≥ 2`, a compact `E ⊂ ℝ^d` with `dim_H E > d/2` determines a
distance set of positive Lebesgue measure (`openai-math:falconer`, Theorem 1.1). The
threshold `d/2` is sharp.

**Method.** Harmonic analysis: a profile decomposition, an angular projection gain, and
incidence estimates applied at the scale where the energy hypotheses hold. The proof is
infinitary; it does not list points.

**Why it worked.** The conjecture had a sharp conjectured threshold and a well-developed
intermediate machinery (restriction estimates, projection theorems). The 2026 proof closes
the last gap rather than inventing a new field.

**What does not transfer.** HN is a *chromatic* statement, not a measure statement. Falconer
gives "the distance set is large"; HN needs "this finite graph is not k-colourable". A
measure-theoretic existence proof supplies neither a graph nor an algorithm. The one
transferable idea is the threshold: a sharp geometric inequality can settle an entire
dimension range.

## 2. Kakeya in ℝ³ (2025)

**Solved.** Every Kakeya set in ℝ³ has Minkowski and Hausdorff dimension 3
(`arXiv:2502.17655`).

**Method.** A volume estimate for unions of `δ`-tubes with a convex-set avoidance
condition, proved by partitioning the tubes according to the convex sets that contain them.
The result is a quantitative lower bound on the union volume, then a limiting argument.

**Why it worked.** The 3D case has enough structure for the tube partition to close: the
proof reduces to a volume estimate for unions of tubes that avoid a common convex set, then
iterates that estimate over scales.

**What does not transfer.** There is no colouring content. The technique is a model for
"prove a lower bound by a quantitative incidence estimate", which is the same shape as the
tiling/density bounds behind `v₇ > 6,992`, but it does not produce a graph. Kakeya in
ℝ^n for n ≥ 4 remains open.

## 3. Separable Jacobian conjecture in characteristic 2 (2026)

**Solved (negatively).** Over the algebraic closure of `𝔽₂`, the map

```
F = (x + x²y + x⁴ + x⁶y²,  y + x⁵ + x⁶y + x⁷y² + x⁸y³)
```

has `det JF = 1` and collapses `(0,1)`, `(1,0)`, `(1,1)` to one point; the function-field
extension has degree 3 and is separable (`arXiv:2608.02634`).

**Method.** An explicit map, then a function-field degree computation: show the extension
`k(x,y)/k(P,Q)` has the stated degree and is separable, and read the collision off the
three points. The proof is a finite certificate plus an irreducibility argument.

**Why it worked.** Characteristic 2 removes the separability that the classical argument
uses, and the explicit map is small enough to verify by hand. The counterexample descends
from a three-variable map of Huq-Kuruvilla by a coordinate permutation.

**What does not transfer.** Characteristic 0 is exactly the case where the separable and
inseparable phenomena collapse; the map does not lift. It is a warning, not a construction:
the same "find a small explicit map" strategy is what the plane search is doing, and here
it succeeded because the field was wrong.

## 4. Tame generators problem in dimension 3 (2003)

**Solved (negatively).** `Aut(k[x,y,z])` contains wild automorphisms
(`doi:10.1090/s0894-0347-03-00440-5`). The tame subgroup is strictly smaller than the full
automorphism group.

**Method.** A rigidity/degree argument: assign a degree to a hypothetical counterexample,
show it can be reduced, and derive a contradiction from the leading terms. The proof is
elementary in the sense that it uses no heavy machinery, but it is a genuine structural
theorem.

**Why it worked.** The three-variable automorphism group has enough room for a
counterexample and enough structure to rule out the tame case. The two-variable case is
classical (every automorphism is tame); dimension 3 is the first place the phenomenon
appears.

**What does not transfer.** The Jacobian condition is not used. A wild automorphism is
still invertible; the tame problem is about generation, not about `det J`. The transferable
lesson is negative: a structural conjecture that looks like a stepping stone can be false
while the analytic one remains open.

## 5. Markus–Yamabe, dimension 2 (1995) — and its failure in dimension 3 (1997)

**Solved.** In dimension 2, a local diffeomorphism of `ℝ²` whose Jacobian has negative trace
everywhere has 0 as a global attractor (`doi:10.4064/ap-62-1-45-74`).

**Method.** Planar dynamics: index theory, the Poincaré–Bendixson theorem, and a careful
analysis of the divergence. The proof exploits dimension 2 in an essential way (no chaotic
invariant sets, no 3D homoclinic tangles).

**Why it worked.** Dimension 2 is rigid enough that the divergence sign controls the global
phase portrait.

**The boundary.** In dimension 3 the conjecture is **false**: Cima, van den Essen, Gasull
and Hubbers exhibit a polynomial counterexample (`doi:10.1006/aima.1997.1673`). The planar
theorem does not lift.

**Why this is the best analogue in the corpus.** It is the exact shape of the risk in the
Jacobian workstream: a conjecture true in dimension 2, with a proof that uses dimension 2,
and a counterexample in dimension 3. It is also a warning about "the plane case is the hard
one": sometimes the plane case is the *true* one and the higher-dimensional case is where
the counterexample lives. The Jacobian record has the same structure in 2026: JC(3) is
false, JC(2) is open.

## 6. Sato–Tate (2011)

**Solved.** For an elliptic curve `E/Q` without CM, the Frobenius traces equidistribute
with respect to the semicircle law (`doi:10.2977/prims/31`, building on Clozel–Harris–
Shepherd-Barron–Taylor).

**Method.** Potential automorphy: prove that a symmetric power of the Galois representation
becomes automorphic after a finite extension, then apply the Arthur–Clozel classification
and the Chebotarev density theorem.

**Why it worked.** The conjecture is a statement about *all* primes, and potential
automorphy gives exactly the finite list of symmetric powers needed. The method is
modularity lifting in a different guise.

**What does not transfer.** BSD is a statement about a single value `s = 1`; Sato–Tate is a
statement about the distribution of local data. The automorphy method does not compute
`L(E,1)` or the rank. It is a model for "prove the local statistics, then use a global
density theorem", not a route to the rank.

## 7. Fermat's Last Theorem (1995)

**Solved.** `x^n + y^n = z^n` has no positive integer solutions for `n ≥ 3`
(`doi:10.2307/2118559`).

**Method.** Modularity lifting: every semistable elliptic curve over Q is modular; a
hypothetical FLT solution produces a semistable curve whose mod-3 representation is
impossible; the contradiction is a deformation-theoretic computation (the `R = T`
theorem).

**Why it worked.** The problem was reduced to a finite, rigid deformation problem with a
computable obstruction. The proof is long but every step is a statement about a finite
Galois representation.

**What does not transfer.** There is no direct bridge from modularity to the Jacobian
conjecture or to BSD's rank. The transferable meta-lesson is the shape: reduce to a finite
deformation problem, prove the deformation ring is the Hecke ring, compute the obstruction.
The Jacobian's Newton-polygon cascade is the analogous finite reduction; the missing piece
is the analogue of `R = T`.

## 8. Gauss class number problem (1976–1986)

**Solved.** The number of imaginary quadratic fields of class number 1 is finite, and the
list is known (Heegner 1952; Baker 1966; Stark 1967 for the list; Goldfeld 1976 and
Gross–Zagier 1986 for the effective lower bound).

**Method.** Goldfeld reduced the effective class-number bound to the existence of an
elliptic curve with a zero of order ≥ 3 at `s = 1` (`goldfeld:numdam-1976`). Gross–Zagier
supplied the rank-3 curve via Heegner points (`doi:10.1007/BF01388809`). The effective
bound is then a height computation.

**Why it worked.** A problem about class numbers was converted into a problem about the
*derivative* of an L-function, and the derivative had a geometric interpretation as a
Heegner point height. The analytic input is exactly the GZK mechanism.

**What does not transfer.** The class-number problem has a one-sided inequality
(`h(d) → ∞`), so a single rank-3 curve suffices. BSD needs an equality for every curve, and
a single curve cannot supply it. The transferable part is the reduction: find an analytic
quantity that is also a height, then bound the height.

## 9. Congruent numbers (families, 2014)

**Partially solved.** Infinitely many `n` are congruent; Tian's Heegner-point construction
covers integers with many prime divisors (`arXiv:1210.8231`). The full problem is
equivalent to BSD for the family `E_n : n y² = x³ − x`.

**Method.** Heegner points on the congruent-number curve, with a height computation to show
the point is non-torsion. Elkies does the same for `|D| < 10^6` (`arXiv:math/0208056`).

**Why it worked.** The curve has CM, so Heegner points are explicit and their heights are
computable. The sign of the functional equation is determined by `n mod 8`, which gives a
free supply of odd-rank curves.

**What does not transfer.** The CM structure is special. A general elliptic curve has no
explicit Heegner point, and the construction does not produce a rank-mismatch certificate.
The transferable lesson is the calibration value: the congruent-number family is where a
rank pipeline should be tested first, because the answer is known in large ranges.

## 10. JC(2) degree ≤ 104 (1983–2019)

**Solved.** Every planar Keller map of total degree ≤ 104 is invertible
(`doi:10.1515/crll.1983.340.140` for 100; `arXiv:1902.05923` for 104).

**Method.** Newton polygon and root-configuration analysis. For each degree pair, the
leading forms of `P` and `Q` are constrained by the Jacobian identity; the constraints
reduce to a finite list of root configurations, and each is excluded by a valuation
argument or an explicit computation.

**Why it worked.** The Jacobian identity is a polynomial identity, so the top-degree block
is a finite determinantal condition. The Newton polygon is the bookkeeping device that
makes the induction terminate.

**What does not transfer.** The method gives a *bound*, not a counterexample. At degree 105
the number of configurations grows, and the (72,108) case is the first where the finite
list has more than one member. The cascade is the right algorithm, but the output is a
degree bound unless a configuration survives.

## 11. Minimum 6-chromatic 3-distance graph (2026)

**Solved.** `min{|P| : χ(G_D(P)) ≥ 6} = 14` for `D = {1, 1/√3, 2}`
(`github:HeliCorgi/fourteen-points-six-colors`).

**Method.** Three layers:

1. an explicit 14-point set `G14` with χ = 6, verified by exact arithmetic in `ℚ(√3)` and
   three independent SAT solvers;
2. a structural lemma that the distance triple admits no K₅ anywhere in the plane, so the
   6-chromaticity is not clique-driven;
3. a complete catalog argument: 1,497,264 edge-6-critical graphs on ≤ 13 vertices
   (Lalonde 2025), filtered to 85,431 candidates, all certified non-embeddable by interval
   arithmetic.

**Why it worked.** The lower bound is a finite catalog plus a checkable non-embeddability
certificate. The upper bound is a 14-point set. Nothing is asymptotic.

**What does not transfer.** The catalog is only complete for ≤ 13 vertices and this
distance triple. For a 6-chromatic *unit-distance* graph, the vertex count is unknown and
the catalog is far beyond reach. The transferable lesson is the three-layer structure:
explicit object, structural lemma, exhaustive lower bound with certificates.

## 12. χ(ℝ²) ≥ 5 (2018)

**Solved.** There is a finite unit-distance graph that is not 4-colourable
(`arXiv:1804.02385`).

**Method.** Gadget composition:

1. `H` (7 vertices, 12 edges) has four essential 4-colourings, two with a monochromatic
   √3-triangle;
2. `L` (52 copies of `H`) forces at least one copy to carry such a triple;
3. `M` contains `H` and forbids the triple in that copy;
4. assembling 52 copies of `M` against `L` gives a non-4-colourable graph.

**Why it worked.** The gadget is a *pattern-forcing* object, not an equality gadget. The
composition is finite and the final graph is explicit; the colouring obstruction is a
finite case analysis, later replaced by a SAT proof.

**What does not transfer to χ = 7.** The same composition pattern at level 5 would need a
virtual edge (a unit-distance gadget forcing two non-adjacent vertices to differ in every
5-colouring). No such gadget is known. The 2018 method produced 5 from 4 because the
pattern (a monochromatic √3-triangle) was already known to be forceable; at level 5 the
pattern is missing.

---

## The common shape

Every solved case above has three components:

1. **A finite object.** A graph, a map, a point set, a catalog, a deformation ring. No
   solution in this list is an existence proof by compactness.
2. **A finite check.** Exact arithmetic, SAT + DRAT, interval arithmetic, a Gröbner basis,
   or a finite case split. The check is reproducible by a third party.
3. **A reduction theorem.** A statement that the original problem is equivalent to the
   finite object. The reduction is where the mathematics lives; the check is where the
   trust lives.

The three targets fail at different components:

- **HN χ = 7** has a reduction (de Bruijn–Erdős) but no finite object at the required size;
  the 6,993-vertex floor is the obstruction.
- **Jacobian** has a finite object (the two Newton configurations) but the reduction is
  one-directional; the reconstruction gap means the finite object is not equivalent to the
  target.
- **BSD** has a finite object (a curve and two certificates) but no known instance; the
  reduction is a theorem (GZK) that rules out the easy region, and the hard region has no
  positive example.

That is the honest summary of the corpus: the solved analogues are solved because all three
components are present, and each of our targets is missing exactly one.
