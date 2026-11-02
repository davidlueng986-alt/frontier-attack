# Transfer: corpus → decisions

Only methods with an `index.json` item enter. "Wall" is the blocker in `STATUS.md`; EV is
the chance the method moves it, not solves the problem.

## 1. Method table

| method (item) | where it solved something | our wall? | cost | EV |
|---|---|---|---|---|
| SAT + DRAT/LRAT (`arXiv:1610.06229`, `arXiv:1805.12181`, `doi:10.1007/978-3-319-09284-3_31`) | 1,581 → 553 5-chromatic; 200 TB Pythagorean certificate | certifies; never generates | hours at 500 v; 200 TB precedent | terminal only |
| Clausal minimization (`arXiv:1907.00929`, `arXiv:2010.12665`) | 553 → 529 → 509; all χ=5 | needs a χ=6 seed; none | days per seed | 0 now |
| Symmetry breaking (`doi:10.1109/dac.2003.1219135`) | no colouring record in corpus | yes: Aut(G) in catalog and at 6,993 v | days | M |
| Catalog + exact arithmetic + SAT (`github:HeliCorgi/fourteen-points-six-colors`) | minimum χ=6 3-distance graph = 14 v | gadget catalog yes; full unit-distance catalog dead by `v6 > 24` | 2–3 weeks | H |
| Gadget composition (`arXiv:1804.02385`) | χ(ℝ²) ≥ 5 | only route to 6; level-5 primitive missing | weeks | H iff primitive |
| Density/tiling (`arXiv:2010.12668`, `arXiv:2303.14714`) | `v7 > 6,992`, `v6 > 24` | kill criterion only | hours | 0 |
| Edge cascade (`arXiv:2204.14178`, `abhyankar:tifr57`) | (72,108) → 2 configurations; degree 100/104 | 4 genuine edges survive; 408 dims left | weeks | M |
| Approximate roots (`doi:10.1515/crll.1983.340.140`) | degree ≤ 100, first planar bound | yes: couples the edges | 1–2 weeks | H |
| Gröbner + saturation (`doi:10.1007/978-3-319-16721-3`, `doi:10.1016/0001-8708(82)90048-2`) | 8-unknown vertices system | no: 445 unknowns | prohibitive | 0 |
| Frameworks (`doi:10.37236/9210`) | finite reformulation of 2D Keller maps | yes: the sufficiency side | 1 week to scope | H/M |
| Resultants (`doi:10.1090/cbms/097`, `tao:jacobian-2026`) | C³ counterexample is a resultant | reconstruction only | medium | M |
| Modular / char-p (`doi:10.1016/j.jpaa.2024.107797`, `arXiv:2608.02634`) | char-2 counterexample | filter only; no Q-transfer lemma in corpus | low | L |
| GZK / Euler systems (`arXiv:0707.0032`, `doi:10.1007/BF01388809`, `doi:10.1007/978-0-8176-4575-5_11`) | rank ≤ 1 | already closed | done | 0 |

## 2. HN: generate graphs, do not validate them

Our wall is the level-5 virtual edge: a unit-distance gadget forcing its terminals to differ
in every 5-colouring. With one, composing against `arXiv:2010.12656` (16 v, distances 1 and
φ) or `arXiv:1909.13177` (426 v, 1 and 2) gives a 6-chromatic unit-distance graph
(`plan/hadwiger-nelson.md` §5).

SAT+DRAT is not a generator. Every corpus SAT result starts from an existing graph:
1,581 → 553 (`arXiv:1805.12181`), 553 → 529 (`arXiv:1907.00929`), 529 → 509
(`arXiv:2010.12665`). A DRAT trace proves "this graph is not 5-colourable"; it cannot choose
coordinates. Local search has no corpus item and is the wrong objective: all
19,998 candidates tested are 4-colourable and admit a 5-colouring with `u = v` — the landscape
is flat. Symmetry breaking (`doi:10.1109/dac.2003.1219135`) is a throughput multiplier on
whatever enumeration we run, not a change of object.

The one corpus method that has produced a new minimum χ=6 3-distance graph is the three-layer
catalog of `github:HeliCorgi/fourteen-points-six-colors`: explicit object, structural lemma,
then an exhaustive catalog (1,497,264 edge-6-critical graphs on ≤13 vertices → 85,431
candidates → all non-embeddable by interval arithmetic). That is graph generation with a completeness
claim, and it matches the virtual edge: a finite graph with a ∀-over-colourings property.
**Move: generate, do not validate.**

Plan:

1. Enumerate abstract `G` on n = 8..16 with terminals u,v: `G` is 5-colourable, `G+uv` is
   not, edge-minimal. Canonical dedupe against `doi:10.1016/j.dam.2012.07.018`.
2. Geometric filter: exact embeddings with `|uv|²` equal to the replaced edge (3 or φ²); reject
   by interval arithmetic (`arXiv:1611.02831`). This is the bottleneck, not the colouring.
3. Certify survivors with two solvers plus DRAT (`arXiv:1610.06229`).
4. Compose with the 16-vertex φ graph; verify every unit edge exactly, then χ ≥ 6 by DRAT.
5. Gate: n ≤ 16 with no embeddable gadget is a publishable gadget-size lower bound; then
   n = 17..20.

Validate by reproducing G14 first. Cost 2–3 weeks. It yields the object or a size bound; no
other corpus method yields either.

## 3. Jacobian: approximate roots, not more cascade

Wall: `(75,125)` emptiness undecided. The full cascade finds 4 genuine edge conditions, not
16, all non-empty — dimensions 4 and 1, two vacuous — leaving 408 coefficients free
(`runs/jacobian_75_125_cascade_full.json`). A filter, not a decision. Reconstruction is
`PARTIAL`.

The item is Moh, `doi:10.1515/crll.1983.340.140`: the first planar bound, degree ≤ 100,
obtained by **approximate roots**. Our cascade solves each edge's leading form independently;
each non-vacuous edge already produces `P_e ∝ G³`, `Q_e ∝ G⁵`, the approximate-root shape,
but only as a leading form. Approximate roots impose the chain globally: `P` is a chain
indexed by the compact edges, and each edge's root is the truncation of the next. The cascade
becomes triangular, and edge i+1 either fixes or kills edge i's free parameters.

`(75,125)` is the right test because `gcd(75,125) = 25` (`arXiv:1401.1784` gives gcd ≥ 16 and
gcd ≠ 2p). So `P` has an approximate 25th root of degree 3 and `Q` one of degree 5: ten root
coefficients before the lower-edge corrections, not 408 free coefficients.

Plan:

1. Re-derive the `3Δ′` / `5Δ′` chain symbolically. Gate: reproduce our vertices-only support
   obstruction (min x-degree 8 → `[1]`) as the degenerate case.
2. Re-run the 4 genuine edge conditions in chain coordinates; the edge-1 family must survive
   edge 2 or die there.
3. Emptiness at any edge gives `max(deg P, deg Q) ≥ 125`, valid without reconstruction
   (`plan/jacobian.md` §3).
4. Only if the chain survives, attack sufficiency with frameworks (`doi:10.37236/9210`) or
   resultants (`tao:jacobian-2026`); not before.

Second choice: Newton–Puiseux (`doi:10.1007/978-3-662-05355-3`) expands each approximate root;
it is not a decider. Modular (`doi:10.1016/j.jpaa.2024.107797`) is a filter only:
`arXiv:2608.02634` shows the char-2 world is false, and the corpus has no Q-transfer lemma.
Degree bounds (`arXiv:1902.05923`) are the wrong tool: 100 → 104 in 36 years is the rate of
more cascade.

## 4. BSD: no corpus item suggests GZK is wrong

No. GZK is stated in `arXiv:0707.0032` Theorem 1.1, `doi:10.1007/BF01388809`, and
`doi:10.1007/978-0-8176-4575-5_11`, and every corpus computation agrees. Two caveats,
neither of which reopens the route:

- Miller (`doi:10.1112/s1461157011000180`) proves full BSD for 16,714 of 16,725 curves of
  conductor < 5,000 with k ≤ 1. Given GZK, the 11 missing curves cannot be rank mismatches;
  the gap is the full formula (Sha/regulator). "Closed" is true for our terminal object
  `(E, r, k)` with `r ≠ k`, not for every k ≤ 1 curve's full BSD.
- Our admissible-shape list uses parity `r ≡ k (mod 2)`. That is the parity conjecture, proved
  only under `Sha[2]` finiteness and local conditions at 2 and 3
  (`doi:10.4007/annals.2010.172.567`). GZK is unconditional; the shape list is not. Label it
  conditional or drop it: `(r,k) = (1,2)` is admissible and our list currently excludes it.

## 5. Three highest-value unused items

1. `doi:10.1515/crll.1983.340.140` — Moh. Action: build the `(75,125)` approximate-root
   chain, gate on the vertices-only obstruction, then descend. Owner: Jacobian. Largest
   reduction: 445 → ten root coefficients.
2. `doi:10.37236/9210` — Borisov frameworks. Action: one-week scoping pass — does a
   `(75,125)` normalized pair map to a framework? If yes, the reconstruction gap gets a
   sufficiency test instead of a rational-map inverse; if no, record it. Owner: Jacobian.
3. `doi:10.1109/dac.2003.1219135` — Shatter. Action: lex-leader predicates in the gadget
   catalog and in any 6,993-vertex UNSAT instance; measure candidates/min on the 509 seed
   before and after. Owner: HN.
