# Open analogues

Sibling problems that are still open, with the strongest bounds the corpus can support.
"Open" here means no finite certificate is known; it does not mean nothing is known.

## 1. The chromatic number of the plane

```
6 ≤ χ(ℝ²) ≤ 7
```

- **Lower bound 6**: `openai-math:euclidean-five`, Theorem 1.1. No proper 5-colouring, even
  with arbitrary colour classes.
- **Upper bound 7**: the hexagonal tiling (Isbell 1950), also in the 2026 paper.
- **The gap**: exactly 6 versus 7. No explicit 6-chromatic unit-distance graph is known,
  and no 7-colouring better than the hexagon tiling is known.

**Why it is hard.** The 2026 lower bound is infinitary: it proves that a 5-colouring cannot
exist but does not exhibit a finite obstruction. By de Bruijn–Erdős, a finite 6-chromatic
unit-distance graph must exist, but the proof gives no size bound. The counting bound says
a 7-chromatic graph needs at least 6,993 vertices (`arXiv:2010.12668`,
`arXiv:2303.14714`), so the search space for the terminal object is enormous and the only
known 5-chromatic seed is 509 vertices.

**What would settle it.** For χ = 7: an explicit 7-chromatic unit-distance graph plus a
DRAT/LRAT proof that it is not 6-colourable. For χ = 6: an explicit 6-colouring of the
plane, which is a construction, not a finite graph.

## 2. The chromatic number of ℝ³

```
6 ≤ χ(ℝ³) ≤ 15
```

- **Lower bound 6**: 59- and 60-vertex 6-chromatic unit-distance graphs in ℝ³
  (`mathworld:deGreyGraphs`).
- **Upper bound 15**: the standard 3D grid colouring (15 colours).
- **Refinement**: a 61-vertex triangle-free 5-chromatic graph in ℝ³ (2026), which shows the
  dimension-3 constructions are still improving.

**Why it is easier than the plane.** Three dimensions give more room for gadgets and more
freedom for colourings; the lower bound is explicit, unlike the plane's 6. The gap 6–15 is
wide, and no 7-chromatic example is known.

**Relation to our target.** None directly. Do not import an ℝ³ count into a plane argument:
a graph that is 6-chromatic in ℝ³ need not embed in ℝ², and the plane's unit-distance
constraints are strictly stronger.

## 3. The measurable chromatic number of the plane

```
5 ≤ χ_m(ℝ²) ≤ 7
```

- **Lower bound 5**: Falconer 1981, surveyed in `doi:10.1007/978-0-387-74642-5_9`.
- **Upper bound 7**: the hexagon tiling is measurable.
- **The gap**: 5 versus 7. The unrestricted number is 6 or 7; the measurable number could
  be 5, 6, or 7.

**Why it matters.** The 2026 proof transfers arbitrary colourings to a weak measurable
condition (`openai-math:euclidean-five`, Theorem 1.3) and then obstructs five labels. That
does **not** prove `χ_m(ℝ²) ≥ 6`; the transfer is one-directional. A measurable 5-colouring
remains logically possible, and finding one would be a major result in the opposite
direction.

**Relation to our target.** The measurable variant is the standard test of whether a
colouring argument can be made geometric. A 7-chromatic finite graph would settle the
unrestricted number and say nothing about the measurable one.

## 4. The Erdős unit-distance problem

```
u(n) = max number of unit distances among n points in the plane
conjecture: u(n) = n^{1+o(1)}
```

- **Trivial upper bound**: `O(n^{3/2})` (incidence geometry).
- **2026 result**: `u(n) ≤ C n^β` for some `β < 4/3` (`openai-math:unit-distances`,
  `theorem main`, 137 lines, no `sorry`).
- **Lower bound**: `n^{1+c/log log n}` (the integer lattice construction).

**Why it is hard.** The gap between `n^{1+o(1)}` and `n^{4/3}` is a polynomial gap. The
2026 bound is the first sub-3/2 exponent with a formal proof, but it is not the conjecture.

**Relation to our target.** The same objects (unit-distance graphs in the plane) and the
same extremal flavour, but a counting question, not a colouring question. A counterexample
to the counting conjecture would be a point set with too many unit distances; it would not
give a 7-chromatic graph.

## 5. Kakeya in dimensions ≥ 4

```
Kakeya set conjecture in ℝ^n: every Kakeya set has Hausdorff dimension n
```

- **Solved in ℝ³**: `arXiv:2502.17655` (2025).
- **Open in ℝ^n for n ≥ 4**: the tube-volume method does not yet close; the known bounds
  are below n.

**Why it is hard.** The 3D proof uses a partition of tubes by the convex sets that contain
them, and the combinatorics of that partition degrades in higher dimensions. The polynomial
method gives weaker bounds.

**Relation to our target.** Kakeya is the model for a quantitative incidence bound settling
a conjecture in a fixed dimension. HN's counting bound (`v₇ > 6,992`) is the same style of
result, but it is a lower bound on a graph, not an incidence theorem.

## 6. The planar Jacobian conjecture above degree 125

```
JC(2): every Keller map (P,Q) with J(P,Q) ∈ C* is invertible
```

- **Proved for total degree ≤ 104**: `arXiv:1902.05923`.
- **Below 125**: only `(72,108)` and `(108,72)` survive (`arXiv:2204.14178`, conditional on
  the unrefereed GGHV preprint).
- **The (72,108) claim**: `github:wstrinz/plane-jacobian-72-108`, Theorem A, asserts both
  configurations are empty; if true, the bound rises to 125. Unrefereed.
- **Above 125**: 24 degree pairs are listed by the repo's J4 pass; only 9 have a Newton
  configuration written down. The smallest is `(75,125)`.

**Why it is hard.** The GGHV reduction is one-directional: a counterexample normalises to
`[P,Q] = x²` in `K[x,x⁻¹,y]`, but a solution of the normalised system does not reconstruct
a polynomial map with `det J = 1`. The (72,108) system has 72 unknowns / 92 equations on the
full lattice; Gröbner elimination is not feasible, and the edge cascade is the intended
method.

**Relation to our target.** This *is* our target. The open analogue is the same problem at
a different degree: the dimension-3 case is settled (negatively), the plane case is not.

## 7. The Dixmier conjecture

```
Every endomorphism of the Weyl algebra A_n(C) is an automorphism.
```

- **Open for n ≥ 2**; the first Weyl algebra `A_1` is the classical settled case.
- **Equivalent to JC(n)** (Tsuchimoto; Belov-Kanel–Kontsevich).
- **Evidence**: the conjecture holds for the first Weyl algebra and for some classes of
  endomorphisms.

**Why it is hard.** The Weyl algebra is non-commutative; the Jacobian condition becomes a
commutator condition. The equivalence with JC means a solution of one solves the other, but
the non-commutative formulation has not produced a new degree bound.

**Relation to our target.** A parallel formulation, not a shortcut. If the plane JC is
settled, Dixmier(2) follows; the converse is also true.

## 8. The Mathieu–Zhao conjecture

```
The image of the substitution map f -> f(F) is a Mathieu-Zhao subspace
for every polynomial map F with F(0) = 0.
```

- **Open** for the relevant class of maps.
- **JC implies the Mathieu conjecture** for the appropriate class
  (`doi:10.1007/978-3-030-42136-6_14`).
- **Evidence**: verified for monomial maps and several low-dimensional families.

**Why it is hard.** It is a statement about the image of a polynomial map, not about its
Jacobian; the implication from JC is one-directional. A counterexample to Mathieu–Zhao
would not refute JC. The precise formulation and the known cases are in
`doi:10.1007/978-3-030-42136-6_14`.

**Relation to our target.** The widest currently open statement that JC implies. It is a
candidate for a "weaker theorem that might be provable first", not a route to the terminal
object.

## 9. BSD for rank ≥ 2

```
ord_{s=1} L(E,s) = rank E(Q)
```

- **Proved for analytic rank ≤ 1**: Gross–Zagier + Kolyvagin (`doi:10.1007/BF01388809`,
  `doi:10.1007/978-0-8176-4575-5_11`, `arXiv:0707.0032` Theorem 1.1).
- **Open for rank ≥ 2**: no general theorem; the rank-2 case is the first open one.
- **Verified**: 16,714 of 16,725 curves of conductor < 5,000 with rank 0 or 1
  (`doi:10.1112/s1461157011000180`); LMFDB covers all conductors < 500,000.
- **No mismatch has ever been observed** in any peer-reviewed computation.

**Why it is hard.** The algebraic rank is the expensive certificate. A 2-descent gives an
interval `r1 ≤ r ≤ r2`; it is exact only when `r1 = r2`. For rank ≥ 2 the Selmer group can
have large Sha, so the interval often does not collapse. The analytic rank is certifiable by
zero-counting, but the algebraic rank is not.

**Relation to our target.** This is our target. The open analogue is the rank-2 case,
which is exactly where the first admissible mismatch lives.

## 10. The class number problem, effective form

```
h(d) → ∞ effectively as d → −∞
```

- **Solved** by Goldfeld 1976 plus Gross–Zagier 1986 (`goldfeld:numdam-1976`,
  `doi:10.1007/BF01388809`).
- **The method**: reduce to an elliptic curve with a zero of order ≥ 3 at `s = 1`, then use
  a Heegner point height.
- **What remains**: explicit numerical bounds and the full list of small class numbers
  (Heegner/Baker/Stark for class number 1).

**Why it is listed as open.** The *existence* of an effective bound is settled; the
*explicit* bound is still being improved. The analogue is the shape of the argument, not the
current status.

**Relation to our target.** The best example of a BSD-type analytic input settling a
classical Diophantine problem. It works because the problem is one-sided; BSD is not.

## 11. The (72,108) exclusion

```
No (P,Q) in K[x,x^-1,y]^2 with [P,Q] = x^2 and configuration (1) or (2).
```

- **Claimed** by `github:wstrinz/plane-jacobian-72-108`, Theorem A, and independently by
  Helali (Zenodo DOI `10.5281/zenodo.21479814`).
- **Unrefereed**. The repo's audit (`research/jacobian-wstrinz-audit.md`) found the argument
  substantive but the repository internally inconsistent.
- **If true**: the plane degree bound rises to 125, and the next region is `(75,125)` with
  445 unknowns on the full lattice.

**Why it is hard.** The proof is a case analysis over two Newton configurations, with a
valuation cascade and a finite polynomial identity. Every step must be checked; a single
gap invalidates the bound.

**Relation to our target.** This is the closest open sub-problem to the Jacobian terminal
object. A verified exclusion would not produce a counterexample, but it would move the
search to a region where a counterexample is still possible.

---

## The shape of the open set

The open analogues fall into three groups:

1. **A finite gap with a known method.** Kakeya in ℝ³ (now closed), Falconer (now closed),
   the 14-point minimum (closed). These were solved by a quantitative estimate plus a
   finite certificate. HN 6-vs-7 is not in this group: there is no candidate graph.
2. **A one-directional reduction.** JC(2) and the (72,108) exclusion. The finite object is
   known, but the reduction does not reconstruct the target. Progress is a degree bound,
   not a counterexample.
3. **A theorem that rules out the easy region.** BSD rank ≥ 2. GZK proves the rank ≤ 1
   case; the rank ≥ 2 case has no positive instance and no general method.

The corpus is organised so that a new claim can be placed in one of these groups before any
compute is spent. That placement is the main practical value of the collection.
