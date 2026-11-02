# `compact_open_cover_impossible`: finite-cover core, no finite graph

## Answer

**Finitary: yes only as a finite open cover; no as a finite point set, graph, or combinatorial certificate.** `OAI.PlanarCover.compact_open_cover_impossible` (`lean/OAI/Geometry/PlaneColoring/PlanarCovers.lean:374`) takes an arbitrary open cover `U : ι → Set Q` of a compact `Q ⊇ Metric.ball 0 5`, with point multiplicity ≤ 2 and each `U i` of diameter ≤ 2, and derives `False`. Its only finite-extraction step is `isCompact_univ.elim_finite_subcover U ho` at `PlanarCovers.lean:381`, returning `s : Finset ι`; line 384 applies the `[Fintype ι]` lemma `finite_open_cover_impossible` (`:348`) to the finite subtype `s`. So the obstruction reduces to a finite open cover — finitely many open subsets of `Q`, diameter ≤ 2, multiplicity ≤ 2 — not to a finite point set or unit-distance graph.

**Setup.** The compact space is `Q = Metric.closedBall (0:ℂ) 5` at the call site (`Five.lean:54`); the lemma is general over compact `Q ⊇ Metric.ball 0 5`. The cover is `W : I → Set Q`, where `I = {p : Fin 5 × Set Q // ∃ x, p.2 = connectedComponentIn (U p.1) x}` (`Five.lean:32–33`) and `U i = {z | 0 < (F z).val i}` for continuous `F : closedBall 0 5 → GraphSpace (⊤ : SimpleGraph (Fin 5))` (`Five.lean:23,29`). The given cover is arbitrary, possibly uncountable; it is not assumed finite or locally finite. Compactness extracts a finite subcover, which is then trivially locally finite.

**What is finite.** The only finite set introduced is `s` (plus `Finset` cardinal bookkeeping at 391–392). `p : ι → Q` at line 358 picks a point in each nonempty `U i`, but these points have no unit-distance edges and no colouring constraint. The closing lemma `finite_partition_impossible` (`:302`) is measure-theoretic: the barycentric map `F z = ∑ w_i(z) • p_i` has range containing a ball (`hit_center_subtype`, `:252`, via `hit_center`, `:180`), while `barycenter_two` (`:274`) and `affinePair_null` (`:242`) place that range in a finite union of lines of measure zero; 334–340 contradict positive measure. That is a covering-dimension (Lebesgue-type) argument, not a graph-colouring one. The lemma never mentions 5 colours; the colour count 5 enters only in `Five.lean`, through `GraphSpace (⊤ : SimpleGraph (Fin 5))` (`GraphRealization.lean:93–95`), the geometric realization of K5 in the probability simplex — colour weights, not a unit-distance graph.

**Graph verdict: no.** Nothing in the chain defines a vertex set, edge set, coordinates, or `SimpleGraph ℂ`; the only finite graph-like object is the colour-label `SimpleGraph (Fin 5)`. Within this development a finite graph would require a separate de Bruijn–Erdős compactness step on the infinite unit-distance graph of ℂ; that is not formalized, and the paper calls the resulting `H` existential (`build/sections/introduction.tex:310–323`). It would give only `Classical.choice`, with no size or coordinates.

### Chain, with file:line

1. `Five.lean:54` applies the lemma to `Q = closedBall 0 5` and `W : I → Set Q` (`Five.lean:29,32–33`).
2. `PlanarCovers.lean:380` installs `CompactSpace Q`; `:381` extracts `s : Finset ι`.
3. `PlanarCovers.lean:384` re-indexes by `s` and applies `finite_open_cover_impossible` (`:348`).
4. `PlanarCovers.lean:354` gets a subordinate partition of unity `w` via `PartitionOfUnity.exists_isSubordinate`; `:358` chooses the `p i`.
5. `PlanarCovers.lean:359` applies `finite_partition_impossible` (`:302`); `:334`, `:336`, `:339–340` close the contradiction.

Inputs not captured by the finite cover: the partition of unity (continuous functions), the choice of the `p i`, the winding-number argument, and Lebesgue measure. `finite_packed_cover` (`:9`) is a separate `Finset` extraction (totally bounded sets, used by `finite_disk_overlap_cover`, `:56`), not this chain.

**Size.** No upper bound and no graph vertex count; `s.card` is whatever the chosen subcover has, with no algorithm or estimate. For any hypothetical cover, `diam ≤ 2` gives `area(U_i ∩ ball 0 5) ≤ π`; since the sets cover the radius-5 ball of area `25π`, `∑ area ≥ 25π`, hence `s.card ≥ 25`. Multiplicity ≤ 2 only gives the consistent upper bound `∑ area ≤ 2·25π = 50π`. This is a hand area count (`research/scratch/check_planarcover_core.py`), not Lean-derived; it bounds cover elements, not vertices.

## Objects and numbers

| item | value | ran/cited |
|---|---|---|
| `compact_open_cover_impossible` | `PlanarCovers.lean:374` | cited |
| finite-subcover extraction | `elim_finite_subcover`, `PlanarCovers.lean:381`; `s : Finset ι` | cited |
| finite-index downstream lemma | `finite_open_cover_impossible`, `PlanarCovers.lean:348` | cited |
| measure closing lemma | `finite_partition_impossible`, `PlanarCovers.lean:302` | cited |
| topological input | `hit_center_subtype`, `PlanarCovers.lean:252`; `hit_center`, `PlanarCovers.lean:180` | cited |
| measure-zero input | `barycenter_two`, `PlanarCovers.lean:274`; `affinePair_null`, `PlanarCovers.lean:242` | cited |
| call site / cover | `Five.lean:54`; `Five.lean:29,32–33`; `GraphRealization.lean:93–95` | cited |
| finite vertex/edge data | none; only `SimpleGraph (Fin 5)` at `Five.lean:23` | cited |
| minimum cover cardinal (hypothetical) | `≥ 25` | ran |
| area radius-5 ball / max area diameter-2 set | `25π ≈ 78.5398` / `π ≈ 3.1416` | ran |

## Reproduced locally

- Clone `research/scratch/openai-math` (commit `ad22eb5`); local files are `PlaneColoring/*.lean`, upstream `lean/OAI/Geometry/PlaneColoring/*.lean`.
- Read `PlanarCovers.lean` (key ranges 1–120, 170–402) and `Five.lean:1–131`; only use of `compact_open_cover_impossible` is `Five.lean:54`.
- Ran `research/scratch/check_planarcover_core.py`: prints `Finset`/`Fintype`/`SimpleGraph` hits and computes `25π/π = 25` (multiplicity cap `2·25π = 50π`).
- No Lean build; `elan`/`lake`/`lean` absent.

## What I could not find

- No graph vertex set, edge list, coordinates, `SimpleGraph ℂ`, or `Fin N` vertex type in the chain. `finite_disk_overlap_cover` (`PlanarCovers.lean:56`) produces a `Finset ℂ` of disk centres, but is not a vertex set and not this chain.
- No upper bound on the finite subcover's cardinal and no algorithm producing it; `elim_finite_subcover` is classical existence.
- No formalized de Bruijn–Erdős graph-colouring compactness in the repository; `Compactness.lean` concerns weak-* subsequences/traces, not finite subgraphs. The paper's `H` is explicitly existential.
- No explicit finite 6-chromatic plane unit-distance graph in the sources checked (`research/hn-gadgets.md`); the 509-vertex Parts graph is 5-chromatic, hence not a witness.
