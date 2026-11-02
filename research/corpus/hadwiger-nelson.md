# Hadwiger–Nelson: annotated corpus

The question is whether the chromatic number of the plane is 6 or 7. The corpus splits
into five layers: the 2018 existence proof and its reductions, the 509-vertex record, the
two-distance graphs that are the only known 6-chromatic objects, the counting bounds that
make 7-chromatic graphs necessarily huge, and the 2026 proof that five colours are
impossible.

## The object table

| object | v | e | χ | distances | source |
|---|---:|---:|---:|---|---|
| Moser spindle | 7 | 11 | 4 | 1 | Polymath16 (`polymath16`) |
| de Grey G | 1,581 | 7,877 | 5 | 1 | `arXiv:1804.02385`, `polymath16` |
| de Grey G₀ | 1,585 | 7,909 | 5 | 1 | `polymath16` |
| Heule 2018 | 553 | 2,722 | 5 | 1 | `arXiv:1805.12181` |
| Heule 2019 | 529 | 2,670 | 5 | 1 | `arXiv:1907.00929` |
| Parts 2020 | **509** | **2,442** | 5 | 1 | `arXiv:2010.12665` |
| Parts G₁₆ | 16 | 56 | 6 | 1, φ | `arXiv:2010.12656` |
| Exoo–Ismailescu K | 426 | 2,901 | 6 | 1, 2 | `arXiv:1909.13177` |
| G₁₄ | 14 | 50 | 6 | 1, 1/√3, 2 | `github:HeliCorgi/fourteen-points-six-colors` |
| de Grey R³ | 59/60 | — | 6 | 1 | `mathworld:deGreyGraphs` |

The 509 record is still the smallest 5-chromatic unit-distance graph. Nothing in the
corpus exhibits a 6-chromatic **unit-distance** graph: the 6-chromatic examples use two or
three distances, and the 2026 theorem is non-constructive.

## 1. The 2018 construction

**`arXiv:1804.02385` (de Grey).** The proof is gadget-based, not a search over point sets.
The 7-vertex, 12-edge wheel `H` has four essential 4-colourings, two of which contain a
monochromatic √3-triangle. `L` contains 52 copies of `H` and forces at least one copy to
carry such a triple; `M` is a unit-distance graph containing `H` with no 4-colouring in
which that copy carries a triple. Assembling 52 copies of `M` against `L` gives a
non-4-colourable graph. The exhibited graph has 1,581 vertices; the edge count 7,877 is in
the Polymath16 table (`polymath16`), not in the paper's abstract.

**`polymath16`.** The wiki is the only consolidated source for the intermediate counts:
`G₀` 1,585/7,909, `L` 121/301, `M` 1,345/8,268, diamond 4/5, Moser spindle 7/11. Use it
when a paper omits a count, and cite it as a community ledger, not a referee.

**`arXiv:1805.12181` (Heule).** Replaces the hand-built colouring argument with a clausal
proof: the 1,581-vertex graph is encoded as SAT, the UNSAT proof is minimized, and the
surviving core is a 553-vertex graph. The point is not only the smaller graph but the
checkable certificate.

**`arXiv:1907.00929` (Heule).** Clausal proof optimization reduces 553/2,720 to 529/2,670.
This is the direct predecessor of the 509 record and the reference for how to shrink a
graph without losing the UNSAT core.

## 2. The 509 record

**`arXiv:2010.12665` (Parts).** A property-preserving minimization method, applied to
5-chromatic unit-distance graphs, produces 509 vertices and 2,442 edges. The graph
decomposes into a 374-vertex/1,860-edge part and a 136-vertex/564-edge part. The previous
510-vertex graph had 2,502 edges. This is the graph the repo verifies in
`tools/verify_509.py`; the corpus records the published counts, not the local
re-verification.

**`arXiv:2010.12668` (Parts).** The same author's tiling paper is the source of the
counting bound, not of a graph. It gives a 6-colouring of more than 99.985698% of the
plane, hence every unit-distance graph of order at most 6,992 is 6-colourable. A
7-chromatic graph therefore needs at least **6,993** vertices.

**`arXiv:2303.14714` (de Grey–Parts).** Refines the numerical bounds: `v₇ > 6,992` and
`v₆ > 24`. The phrasing matters: these are lower bounds on the *order* of a k-chromatic
unit-distance graph, obtained from tiling densities, not from a search.

## 3. Two-distance graphs: the 6-chromatic targets

A 6-chromatic unit-distance graph is not known. The known 6-chromatic planar graphs forbid
two distances, and that is exactly what makes them useful: replace every non-unit edge by
a unit-distance gadget that forces its endpoints to differ in every 5-colouring.

**`arXiv:2010.12656` (Parts).** A 16-vertex graph with distances 1 and d, 28 edges of each
length, symmetry group of order 48, and 5-cliques. It is the smallest known 6-chromatic
two-distance graph and the cheapest target for a virtual-edge composition.

**`arXiv:1909.13177` (Exoo–Ismailescu).** A 426-vertex graph with 2,009 unit edges and 892
edges of length 2; every 5-colouring forces three named vertices to split. This is the
large, explicit target.

**`arXiv:1805.00157` (Exoo–Ismailescu).** The human-readable 5-chromatic proof: `G₄₀` has
40 points, 82 unit edges and 59 pairs at distance 8/3, and forces an equal pair; `G₇₉` has
79 points and 165 unit edges. Useful as a sanity check on gadget vocabulary, not as a
record.

**`github:HeliCorgi/fourteen-points-six-colors`.** A 14-vertex graph for the three
distances `{1, 1/√3, 2}` with χ = 6, ω = 4, and a proof that this distance triple admits
no K₅ anywhere in the plane. It also certifies that 14 is minimum: no 13-point witness
exists. This is the smallest 6-chromatic planar distance graph in the corpus and a clean
calibration object for any pipeline that claims to find 6-chromatic graphs.

## 4. The size wall

Two independent counting arguments give the same order of magnitude:

- `arXiv:2010.12668`: order ≤ 6,992 ⇒ 6-colourable, so a 7-chromatic graph has ≥ 6,993
  vertices.
- `arXiv:2303.14714`: `v₇ > 6,992`, `v₆ > 24`.

Consequence for this repo: any construction that glues a gadget onto the 509-vertex graph
produces a graph of order ~516 and is provably 6-colourable. The 509-plus-gadget route to
χ = 7 is dead on size, not on compute. A 7-chromatic target must start from a seed with at
least 6,993 vertices, and no such seed is known.

## 5. The 2026 five-colouring theorem

**`openai-math:euclidean-five`.** Theorem 1.1: there is no proper 5-colouring of ℝ², even
with arbitrary (non-measurable) colour classes. Hence 6 ≤ χ(ℝ²) ≤ 7. The proof is in ZFC
and has two parts: a transfer from arbitrary colourings to a weak measurable condition
(Theorem 1.3), and a geometric obstruction to five labels under that condition (Theorem
1.4). It is a covering-dimension/measure argument, not a finite-graph argument.

**`openai-math:five-lean`.** The Lean 4 formalization is 131 lines with no `sorry`; the
final theorem `no_proper_five_coloring` is at line 124. It quantifies over `ℂ → Fin 5`,
not over a finite vertex set. The paper says the graph and the positive constant are
existential. **No finite witness can be extracted.**

This is the single most important negative fact in the HN corpus: the lower bound moved
from 5 to 6 without producing the object the repo's terminal check consumes. A 6-chromatic
unit-distance graph would be new, but it would only calibrate χ ≥ 6; it would not settle
6 versus 7.

**`openai-math:unit-distances`.** The same 2026 release formalizes a sub-4/3 exponent for
the planar unit-distance counting function: `∃ C β, 0 < C, 1 ≤ β < 4/3, ∀ n, u(n) ≤ C
n^β` (137 lines, no `sorry`). This is a power saving over the trivial `O(n^{3/2})`, not the
conjectured `n^{1+o(1)}`; it is the quantitative sibling of HN, not a solution of it.

## 6. Analogues

**`openai-math:falconer`.** The Falconer distance conjecture is resolved in every
dimension: for compact `E ⊂ ℝ^d` with `dim_H E > d/2`, the distance set has positive
Lebesgue measure. Same objects (distance sets in Euclidean space), opposite outcome: a
full resolution with a sharp threshold. The method is harmonic analysis and projection
theory, not graph colouring.

**`arXiv:2502.17655` (Wang–Zahl).** Every Kakeya set in ℝ³ has Minkowski and Hausdorff
dimension 3. This is the 2025 model of an incidence-geometry conjecture falling to tube
volume estimates. It does not transfer to colouring, but it is the current proof technique
for "the set is too small to avoid this configuration" arguments.

**`mathworld:deGreyGraphs`.** In ℝ³, 59- and 60-vertex unit-distance graphs have χ = 6, and
a 61-vertex triangle-free graph has χ = 5. Dimension 3 is strictly easier than dimension 2
for constructions and strictly less constrained for upper bounds; do not import R³ counts
into a plane argument.

**`doi:10.1007/978-0-387-74642-5_9` (Soifer).** The measurable chromatic number
`χ_m(ℝ²)` is the standard regularity variant. Falconer's 1981 bound gives `χ_m(ℝ²) ≥ 5`;
the exact value is open. The 2026 theorem deliberately drops regularity, which is why it
can prove 6 without producing a finite graph.

## 7. What the corpus implies for the search

1. **The missing primitive is a level-5 virtual edge**, not another 5-chromatic graph. A
   unit-distance gadget that forces its endpoints to differ in every 5-colouring, composed
   with `arXiv:2010.12656` or `arXiv:1909.13177`, would give an explicit 6-chromatic
   unit-distance graph. That is Target 0, not χ = 7.
2. **χ = 7 needs ≥ 6,993 vertices** (`arXiv:2010.12668`, `arXiv:2303.14714`). Any plan
   whose seed is the 509 graph plus a bounded gadget is provably 6-colourable.
3. **The 2026 proof is not a source of witnesses** (`openai-math:five-lean`). It is a
   proof of existence by compactness with no size bound. Extraction is dead; search is the
   only route.
4. **Calibrate on 14, not on 509.** `github:HeliCorgi/fourteen-points-six-colors` is the
   smallest explicit 6-chromatic planar distance graph and comes with exact arithmetic and
   SAT certificates. A pipeline that cannot reproduce it will not reproduce a 6-chromatic
   unit-distance graph.
