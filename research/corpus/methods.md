# Methods that transfer

Each entry is a technique, the corpus item that demonstrates it, and the failure mode that
makes it inapplicable outside its range. "Transfers" means the method has actually been
used on more than one of the three problems, or on a close analogue; it does not mean it is
a solution.

## 1. SAT + DRAT/LRAT

**What it does.** Encodes "this graph is not k-colourable" as a propositional formula,
runs a CDCL solver, and emits a clausal proof. `drat-trim` (or `cake_lpr`) checks the proof
independently. The graph is certified by a finite object; the colouring claim is certified
by a proof log.

**Worked pointers.**

- `arXiv:1805.12181` and `arXiv:1907.00929`: the 553 → 529 reduction is SAT plus clausal
  proof minimization.
- `arXiv:1610.06229`: the DRAT format and the checker.
- `github:HeliCorgi/fourteen-points-six-colors`: χ(G14) ≥ 6 is UNSAT under three
  independent solver families, with exact arithmetic on the edge lengths.

**Failure mode.** SAT proves a property of the *encoded* graph. If the coordinates are
floats, the encoding is not the geometric graph. Every SAT claim in this corpus is paired
with an exact-arithmetic check of the unit distances; the SAT proof alone is not enough.

## 2. Finite-graph compactness and Lean

**What it does.** De Bruijn–Erdős: a graph is k-colourable iff every finite subgraph is.
So an unrestricted lower bound can be proved by a compactness argument without exhibiting a
finite graph. Lean 4 checks the argument.

**Worked pointers.**

- `openai-math:euclidean-five`: the 2026 proof that χ(ℝ²) ≥ 6 is a covering-dimension
  argument, not a finite-graph argument.
- `openai-math:five-lean`: 131 lines, no `sorry`, theorem at line 124. The statement is
  about `ℂ → Fin 5`, not about a finite vertex set.
- `openai-math:unit-distances`: the same release formalizes a sub-4/3 exponent for the
  unit-distance counting function.

**Failure mode.** Compactness gives existence, not a witness. A Lean proof of
`¬ ∃ c : ℂ → Fin 5, ProperColoring 5 c` does not produce a `Fin N` graph, an edge list, or
coordinates. Do not schedule "extract the witness" after such a proof; the extraction is
not in the proof term.

## 3. Density and tiling bounds

**What it does.** Colours a large fraction of the plane with k colours, leaving an
uncoloured region of density δ. Any k-chromatic graph must then have at least `1/δ` edges;
the Erdős relation `e < v^{3/2}` converts that into a lower bound on the order. This turns
a geometric covering into a statement about how large a k-chromatic graph must be.

**Worked pointers.**

- `arXiv:2010.12668`: a 6-colouring of >99.985698% of the plane gives "order ≤ 6,992 ⇒
  6-colourable".
- `arXiv:2303.14714`: refined densities give `v₇ > 6,992` and `v₆ > 24`.
- `arXiv:2010.12665`, §6: the minimization pipeline that produced the 509-vertex graph is
  the complementary upper-bound technique.

**Failure mode.** The bound is on the *order*, not on the structure. It says a 7-chromatic
graph must be huge; it does not say one exists, and it gives no seed. The 6,993 floor is a
kill criterion for small compositions, not a construction.

## 4. Exact symbolic geometry

**What it does.** Represents coordinates in a number field (typically `ℚ(√3)` or a
multiquadratic extension), computes squared distances exactly, and verifies every edge is
exactly 1. No floating point enters the decision.

**Worked pointers.**

- `github:HeliCorgi/fourteen-points-six-colors`: all 91 vertex pairs of G14 are classified
  by exact arithmetic in `ℚ(√3)`; the SAT check is downstream of that.
- `research/hn-509-verify.md` in this repo: 509 vertices, 2,442 edges, worst
  `|d² − 1| = 1.5e-49` under exact/high-precision arithmetic.
- `arXiv:1805.00157`: the G40/G79 construction is given in `ℚ[√2, √3]`-style coordinates,
  which is why the 8/3 and √(11/3) distances are exact.

**Failure mode.** A near-miss is not a unit distance. The 509 graph's worst residual
`|d² − 1| = 1.5e-49` is evidence at the working precision, not a proof; only an exact or
interval argument closes the gap. A graph with `|d² − 1| = 1e-12` is not a unit-distance
graph, and no amount of SAT will fix it.

## 5. Newton polygon and the edge cascade

**What it does.** For a Keller pair, the Jacobian identity forces the top-degree block to
be determinantal. The Newton polygons of `P` and `Q` are similar; their edges encode the
valuations at infinity. The cascade solves the leading-form condition on the top edge, then
descends one edge at a time. At each step the unknowns are the coefficients on one edge, and
the equations are the coefficients of one degree of `det J`.

**Worked pointers.**

- `abhyankar:tifr57`, Theorem (19.4): JC is equivalent to `N(f)` being a triangle, and a
  Keller pair has at most two points at infinity per component.
- `doi:10.1515/crll.1983.340.140` (Moh): the degree-100 bound is a root-configuration
  analysis, the same cascade in a different language.
- `arXiv:2204.14178`, Proposition 4.3: the two surviving `(72,108)` configurations are the
  output of the cascade at the top edge.

**Failure mode.** The cascade produces a finite list of configurations and, if a
configuration survives, a solution in `L^(1) = K[x,x⁻¹,y]`. It does not by itself
reconstruct a polynomial map in `K[x,y]`. That is the reconstruction gap: the implication
`counterexample ⇒ normalised pair` is one-directional.

## 6. Gröbner bases and saturation

**What it does.** Writes the coefficient system as a polynomial ideal, computes a Gröbner
basis, and reads off whether the variety is empty (`[1]`) or positive-dimensional. Saturation
by a product of variables imposes "all four corners are present" without adding equations.

**Worked pointers.**

- `problems/jacobian/config2_system.py` (this repo): 72 unknowns / 92 equations on the full
  lattice; 8 / 7 on the vertices.
- The vertices-only system has basis `a1*b1 − 1, a2, a3, b2, b3`; saturating by
  `a2*a3*b2*b3` gives the unit ideal.
- `arXiv:1708.07936`: the corner enumerator that supplies the supports.

**Failure mode.** Gröbner elimination is doubly exponential in the worst case and the full
`(72,108)` lattice (72 unknowns) is out of reach. The vertices-only system is a truncation:
its `[1]` is consistent with a solution on the full lattice, not a proof that none exists.
A `[1]` from a truncated system is a negative result about the truncation, not about JC.

## 7. Degree reduction

**What it does.** Replaces "all polynomial maps" by a finite-dimensional normal form:
cubic-linear (Drużkowski), symmetric Jacobian (de Bondt–van den Essen), or a fixed degree
pair with a classified Newton polygon (Starkov, GGHV).

**Worked pointers.**

- `doi:10.1007/BF01459126`: JC(n) reduces to `F = x + (Ax)^{*3}`.
- `doi:10.1090/s0002-9939-05-07570-2`: JC reduces to symmetric Jacobian matrices.
- `doi:10.15393/j3.art.2016.3510` and `doi:10.15393/j3.art.2017.3870`: complete
  classifications for total degree 3 and 4.

**Failure mode.** A reduction changes the object, not the difficulty. The cubic-linear
reduction is equivalent to JC, so a search in that class is not a weakening. The symmetric
reduction is also equivalent. Do not treat a reduced class as "easier" unless the reduction
also lowers the degree or the number of unknowns.

## 8. Approximate functional equation and rigorous zero-counting

**What it does.** The AFE writes `L(E,s)` as a rapidly convergent sum of length about
`√N`. Evaluating it at several points gives a numerical order of vanishing. That is a
scout, not a certificate. To certify the order, count zeros in a rectangle with the
argument principle and interval arithmetic, and check that the count matches the expected
order.

**Worked pointers.**

- `arXiv:math/0412181` (Rubinstein): the practical reference for the AFE and zero-counting.
- `arXiv:math/0207280` (Dokchitser): computing special values and Taylor coefficients from
  a smoothed functional equation.
- `lmfdb:reliability`: the database's displayed analytic ranks are rigorously verified,
  including ranks > 1, by this method.
- `research/bsd-rank-smoke.md` (this repo): a PARI-free AFE estimator that returns 0, 1, 2,
  3 on 11a1, 37a1, 389a1, 5077a1, and matches `L(1)` to 29 digits. It is **numerical**, not
  certified.

**Failure mode.** A pointwise interval around a single value can prove it is non-zero; it
cannot prove an exact zero. `ellL1` on 389a1 returns `−5.4e−29` at default precision, which
is numerical noise, not a proof that `L(1) = 0`. The AFE ladder must be paired with a
zero-count, not a tolerance.

## 9. 2-descent

**What it does.** Computes the 2-Selmer group and a rank interval `r1 ≤ rank E(Q) ≤ r2`.
When `r1 = r2`, the rank is proved. The upper bound comes from the Selmer rank; the lower
bound comes from explicit points.

**Worked pointers.**

- `rubinstein:records`: the rank 4/5/6 calibration curves, all with proved ranks.
- `research/bsd-rank-tooling.md` (this repo): `ellrank` returns `[r1, r2, s, L]`; equality
  is exact; a gap is not a rank.
- `arXiv:1006.1002` (Bhargava–Shankar): the average size of the 2-Selmer group is 3, which
  is why the interval usually collapses for random curves and why the hard cases are
  special.

**Failure mode.** `r1 < r2` is an interval, not a rank. Sha can have 4-torsion, in which
case equality is impossible; the PARI documentation gives `E = [−113², 0]` with
`[0, 2, 0, []]` as an example. A pipeline that reports `r2` as the rank is reporting an
upper bound.

## 10. Heegner points and Euler systems

**What it does.** Constructs a rational point on `E` from a CM point on a modular curve,
computes its canonical height, and uses the Gross–Zagier formula to relate the height to
`L'(E,1)`. Kolyvagin's Euler system then bounds the Selmer group, giving
`r_an = r_MW` for `r_an ≤ 1`.

**Worked pointers.**

- `doi:10.1007/BF01388809`: the Gross–Zagier formula.
- `doi:10.1007/978-0-8176-4575-5_11`: Euler systems.
- `arXiv:0707.0032`: explicit Heegner points, Theorem 1.1, and the rank-2 Selmer corank.
- `arXiv:math/0208056`: the congruent-number family, where Heegner points are explicit.

**Failure mode.** The construction works for `r_an ≤ 1`. For `r_an ≥ 2` there is no
analogous point, and the Euler system bounds the Selmer group but does not compute the rank.
This is the exact boundary of the GZK gate.

## 11. Interval arithmetic

**What it does.** Replaces floating-point values by midpoint-radius intervals and propagates
rigorous error bounds through every operation. A result is certified if the final interval
has the required sign or contains the required count.

**Worked pointers.**

- `arXiv:1611.02831` (Arb): the reference implementation.
- `problems/bsd/certified_rank.py` (this repo): certifies `k = 0` (11a1) and `k = 1` (37a1),
  refuses `k = 3` (5077a1) because a pointwise interval cannot prove an exact zero.
- `github:HeliCorgi/fourteen-points-six-colors`: interval arithmetic certifies the
  non-embeddability of 85,431 candidate graphs.

**Failure mode.** Interval arithmetic is only as good as the formula it is applied to. A
pointwise interval around `L(1)` cannot certify a zero of order 2; the argument principle
is needed. An interval that contains zero is not a proof of vanishing, and an interval that
excludes zero is not a proof of the order.

## 12. Average-rank and Selmer distribution

**What it does.** Counts integral orbits of a representation to compute the average size of
the 2-Selmer group, then converts the average Selmer size into a bound on the average rank.

**Worked pointers.**

- `arXiv:1006.1002`: mean 2-Selmer size 3, so limsup average rank ≤ 1.5.
- `rubinstein:records`: the rank census, with 2.3 billion curves of rank ≥ 3 and 99.973%
  proved by explicit points plus a 2-descent upper bound.
- `lmfdb:completeness`: the complete region below 500,000.

**Failure mode.** An average bound says nothing about a particular curve. The existence of
a rank-4 curve is not excluded by an average-rank bound of 1.5; it is exactly the kind of
rare event that averages hide. Use this method to choose a search region, not to prove a
rank.

---

## Cross-problem transfers

| technique | HN | Jacobian | BSD |
|---|---|---|---|
| exact arithmetic / certificates | unit distances | polynomial identities | L-values |
| finite reduction | graph / gadget | Newton polygon | Selmer group |
| independent checker | DRAT | Gröbner / Lean | 2-descent / zero-count |
| a priori bound | 6,993 vertices | degree 104 / 125 | conductor 500,000 |
| missing piece | level-5 virtual edge | reconstruction map | rank ≥ 2 certificate |

The last row is the useful one. Each problem has a technique that works and a single missing
component; the corpus is arranged so that the missing component is visible before compute
is spent.


## 13. Second wave: the certificate and algebra toolbox

**LLL and PSLQ.** `doi:10.1007/bf01457454` (LLL) is polynomial-time lattice reduction;
`doi:10.1090/s0025-5718-99-00995-3` (PSLQ) finds integer relations. Together they turn a
high-precision numerical identity into an exact algebraic one. Worked pointer: recovering
the coefficients of a candidate L-value or a polynomial identity before asking a CAS to
prove it.

**Groebner complexity.** `doi:10.1016/0001-8708(82)90048-2` (Mayr–Meyer) is the
double-exponential lower bound; `doi:10.1007/978-3-319-16721-3` (Cox–Little–O'Shea) and
`doi:10.1090/cbms/097` (Sturmfels) are the algorithm references. Worked pointer: the
`(72,108)` full lattice has 72 unknowns, so the bound is a reason to use the edge cascade,
not a Groebner basis.

**Newton–Puiseux and real algebra.** `doi:10.1007/978-3-662-05355-3` (Basu–Pollack–Roy)
covers Newton–Puiseux expansions and quantifier elimination. Worked pointer: the
`(75,125)` edge cascade and the exact unit-distance checks both need this symbolic layer.

**SAT encodings and symmetry breaking.** `doi:10.3233/FAIA336` (Handbook of Satisfiability)
is the encoding reference; `doi:10.1109/dac.2003.1219135` (Shatter) constructs lex-leader
symmetry-breaking predicates. Worked pointer: a 6,993-vertex colouring instance is
symmetry-heavy, and the graph automorphism group is the first thing to break.

**Lean and mathlib.** `doi:10.1145/3372885.3373824` describes the library the OpenAI
colouring formalization builds on. Worked pointer: `openai-math:five-lean` imports
`OAI.Geometry.PlaneColoring.*`, which is mathlib-style finite combinatorics plus topology.

**Computer-assisted proof stack.** `doi:10.1090/s0002-9904-1976-14122-5` and
`doi:10.1007/978-3-540-87827-8_28` (four-colour), `doi:10.4007/annals.2005.162.1065` and
`doi:10.1017/fmp.2017.1` (Kepler/Flyspeck), `doi:10.1007/978-3-319-40970-2_15` (Pythagorean)
and `doi:10.1007/978-3-319-09284-3_17` (Erdős discrepancy). Worked pointer: the proof
obligation is a finite certificate plus an independent checker, not a hand argument.

**Negative results.** `doi:10.1038/nature16059` (spectral gap undecidable) and
`doi:10.1080/00029890.1973.11993265` (Hilbert's tenth unsolvable) bound what a general
algorithm can decide. Worked pointer: do not schedule a search for a decision procedure
where the problem is undecidable; schedule a certificate for the specific instance.
