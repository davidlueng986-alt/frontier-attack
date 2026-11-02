# Public Lean source and finite-witness status for "The Euclidean plane is not five-colorable"

## Answer

**Public Lean source: yes.** OpenAI's `openai/math` repository is public and Lean-based (`"language":"Lean"`, 940,236 KB, created 2026-10-06). The proof is at `lean/OAI/Geometry/PlaneColoring/Five.lean` (6,231 bytes), theorem `OAI.EuclideanFiveColor.no_proper_five_coloring` at lines 124–127. It is registered in `lean/formalization.yaml` (lines 1227–1229) and imported by `lean/OAI.lean`. The file with the same name under `lean/ComparatorChallenges/` is only the challenge statement (`sorry`); its JSON names `OAI.Geometry.PlaneColoring.Five` as the solution module.

**Explicit finite unit-distance graph that is not 5-colourable: no.** `Five.lean` proves the infinitary statement `¬ ∃ c : ℂ → Fin 5, ProperColoring 5 c` by a measure-theoretic/topological contradiction. It never defines a finite vertex set, a finite edge set, coordinates, or a cardinality for a non-5-colourable graph. The only explicit finite unit-distance graph used as a colouring obstruction in the development is the 7-vertex Moser spindle in `Moser.lean` (11 unit edges, proved not 3-colourable); it is 4-chromatic and therefore 5-colourable, so it is not the requested witness. The paper's only finite non-5-colourable `H` appears in a corollary, obtained from de Bruijn–Erdős, and the paper says "The graph and the positive constant here are existential" (`introduction.tex:323`); `m = |E(H)|` is a variable, with no value or coordinates.

### What the Lean proof actually does

`proper_to_borel_weak` (`BorelTransfer.lean:436`) converts an arbitrary proper 5-colouring into a measurable weak 5-colouring. `no_measurable_weak_five_coloring` (`Five.lean:105`) derives a contradiction from `fullApprox_exists` and the planar-cover obstruction `PlanarCover.compact_open_cover_impossible` (`PlanarCovers.lean:374`). The `GraphSpace (⊤ : SimpleGraph (Fin 5))` objects in `Five.lean` are probability simplices on the five colour labels (`graphRealization`, `GraphRealization.lean:86–95`), not finite graphs embedded in the plane. The `Finset ℂ` objects elsewhere in `PlaneColoring/` are finite sets of disk centres used in the approximation (`FullApproximation.lean`, `HoleApproximation.lean`), not vertex sets of a unit-distance graph.

### What a witness extraction would take

Semantically, `no_proper_five_coloring` says the infinite unit-distance graph on `ℂ` is not 5-colourable, so de Bruijn–Erdős compactness implies some finite subgraph is not 5-colourable. That inference is not in the Lean development: there is no `SimpleGraph ℂ`, no finite-subgraph theorem, and no cardinality. Adding it would require (i) defining the unit-distance graph on `ℂ`, (ii) formalizing the de Bruijn–Erdős graph-colouring corollary (Mathlib has the compactness ingredient `Finset.rado_selection` in `Mathlib/Combinatorics/Compactness.lean`, but no packaged graph corollary was found), and (iii) applying it, which yields only a `Classical.choice` witness with no size or coordinates. Such a witness would be a 6-chromatic unit-distance graph; no explicit one is known (the 509-vertex record is 5-chromatic).

## Objects and numbers

| item | value | source |
|---|---|---|
| repo | `https://github.com/openai/math`; public; `"language":"Lean"`, `"size":940236` KB; created 2026-10-06, last push 2026-10-30 | `gh api repos/openai/math --jq '.private, .language, .size, .created_at, .pushed_at'` |
| repo root | `.github/`, `CONTENTS.md`, `LICENSE`, `README.md`, `history.md`, `lean/`, `overview.pdf`, `overview.tex`, `preprints/`, `reasoning_traces/` | `gh api repos/openai/math/contents --jq '.[].name'` |
| Lean entry point | `lean/OAI.lean` imports `OAI.Geometry.PlaneColoring.Five` and `OAI.Geometry.PlaneColoring.Main` | file |
| formalization registry | `lean/formalization.yaml` lines 1227–1229: comparator config `ComparatorChallenges/EuclideanFiveColor.json`; declaration `OAI.EuclideanFiveColor.no_proper_five_coloring`; file `OAI/Geometry/PlaneColoring/Five.lean` | file |
| solution file | `lean/OAI/Geometry/PlaneColoring/Five.lean`, 6,231 bytes, 131 lines; no `sorry`/`admit` | `gh api .../Five.lean --jq '.size'`; content inspection |
| theorem | `theorem no_proper_five_coloring : ¬ ∃ c : ℂ → Fin 5, ProperColoring 5 c` (lines 124–127) | `Five.lean` |
| challenge stub | `lean/ComparatorChallenges/EuclideanFiveColor.lean` (388 B) is `import Mathlib` plus `theorem ... := by sorry`; `EuclideanFiveColor.json` has `"solution_module": "OAI.Geometry.PlaneColoring.Five"` | files |
| proof route | `proper_to_borel_weak` (`BorelTransfer.lean:436`) → `no_measurable_weak_five_coloring` (`Five.lean:105`) → `PlanarCover.compact_open_cover_impossible` (`PlanarCovers.lean:374`) | files |
| finite witness explicit? | **No** — no `Finset ℂ` of plane vertices, no `SimpleGraph ℂ`, no `Fin N` vertex set in `Five.lean`; `GraphSpace` is the probability simplex on `Fin 5` labels (`GraphRealization.lean:86–95`) | `Five.lean`; `GraphRealization.lean` |
| only explicit finite unit-distance graph used as a colouring obstruction | Moser spindle: `Fin 7` vertices `{0, A, B, T, uA, uB, uT}`, `A=(√3+i)/2`, `B=(√3−i)/2`, `T=√3`, `u=(5+i√11)/6`; 11 unit edges; proved not 3-colourable | `Moser.lean:128–167` |
| cardinality found | 7 vertices, 11 edges for the Moser gadget; χ=4, hence 5-colourable, **not** a non-5-colourable witness | `Moser.lean`; standard fact |
| cardinality of a non-5-colourable witness | none given; paper writes `m = |E(H)|` for an existential `H` | `introduction.tex:310–323` |
| paper's finite `H` | "Theorem main and graph compactness ... give a finite unit-distance graph `H` that is not five-colorable. Set `m=|E(H)|` and `δ=1/m`." … "The graph and the positive constant here are existential." | `introduction.tex:310–323` |
| paper's stance on finite pictures | "an unrestricted lower bound can also be established without first displaying a finite graph witnessing it" | `introduction.tex:31–36` |
| known explicit lower bound | 509-vertex, 2,442-edge 5-chromatic unit-distance graph (Parts); no explicit 6-chromatic plane unit-distance graph known | `research/hn-gadgets.md:4`; arXiv:2010.12665 |
| build files / toolchain | `lean/lakefile.lean` (16,876 B), `lean/lake-manifest.json` (14,622 B), `lean/lean-toolchain` = `leanprover/lean4:v4.34.1`; `elan`, `lake`, `lean` all absent from PATH; build not attempted | `lean/` directory listing; `shutil.which`; `lean-toolchain` |
| tree size | 174,524 paths; 163,096 `.lean` files; 46 files in `lean/OAI/Geometry/PlaneColoring` | `git ls-tree -r --name-only HEAD`; `gh api .../PlaneColoring --jq '.[].name'` |

## Evidence for the two yes/no answers

### Yes: public Lean source

```text
$ gh api repos/openai/math --jq '.private, .language, .size, .created_at, .pushed_at'
false
Lean
940236
2026-10-06T21:47:02Z
2026-10-30T00:35:38Z

$ gh api repos/openai/math/contents --jq '.[].name'
.github
CONTENTS.md
LICENSE
README.md
history.md
lean
overview.pdf
overview.tex
preprints
reasoning_traces
```

The full solution is the 6,231-byte file listed by:

```text
$ gh api repos/openai/math/contents/lean/OAI/Geometry/PlaneColoring --jq '.[].name'
Accumulation.lean
...
Five.lean
...
WeakColoring.lean
```

`gh search code --owner openai "no_proper_five_coloring"` finds it in both the challenge stub and the solution:

```text
openai/math:lean/ComparatorChallenges/EuclideanFiveColor.lean: theorem no_proper_five_coloring : ¬ ∃ coloring : ℂ → Fin 5, ProperColoring 5 coloring := by
openai/math:lean/OAI/Geometry/PlaneColoring/Five.lean: theorem no_proper_five_coloring : ¬ ∃ c : ℂ → Fin 5, ProperColoring 5 c := by
openai/math:lean/ComparatorChallenges/EuclideanFiveColor.json: "OAI.EuclideanFiveColor.no_proper_five_coloring"
openai/math:lean/formalization.yaml: declaration: OAI.EuclideanFiveColor.no_proper_five_coloring
```

The paper and its Lean scope page are also public:

```text
$ gh search code --owner openai "five-colorable"
openai/math:preprints/The-Euclidean-plane-is-not-five-colorable-September-23-2026/README.md: # [The Euclidean plane is not five-colorable](paper.pdf)
openai/math:lean/docs/158.md: - [The Euclidean plane is not five-colorable](...)
openai/math:preprints/The-Euclidean-plane-is-not-five-colorable-September-23-2026/build/sections/introduction.tex: $H$ that is not five-colorable. Set $m=|E(H)|$ and $\delta=1/m$.
openai/math:lean/formalization.yaml: - title: "The Euclidean plane is not five-colorable"
```

### No: explicit finite witness

Content inspection of the complete 131-line `Five.lean` finds only `SimpleGraph (Fin 5)` (the five-colour target simplex) and `Finset (Fin 5)` (sums over colours), no `Fin 7`, no `SimpleGraph ℂ`, and no `Finset ℂ` vertex set:

```text
$ Select-String -Path lean\OAI\Geometry\PlaneColoring\Five.lean -Pattern 'Finset|SimpleGraph|Fin 7|Fin 6|card'
10: lemma graphRealization_positive_three {V : Type*} {G : SimpleGraph V} {p : V → ℝ}
23: (F : C(Metric.closedBall (0:ℂ) 5,GraphSpace (⊤ : SimpleGraph (Fin 5))))
59: have hsum := Finset.sum_nonpos (s := (Finset.univ : Finset (Fin 5))) (fun i _ => hh i)
95: (F : C(Metric.closedBall (0:ℂ) 5,GraphSpace (⊤ : SimpleGraph (Fin 5)))) :
```

The theorem's own statement is infinitary, and the proof ends at `no_measurable_weak_five_coloring`:

```lean
theorem no_proper_five_coloring : ¬ ∃ c : ℂ → Fin 5, ProperColoring 5 c := by
  rintro ⟨c,hc⟩
  obtain ⟨d,hd,hw⟩ := proper_to_borel_weak hc
  exact no_measurable_weak_five_coloring hd hw
```

The paper's only mention of a finite witness is explicitly existential:

```text
Theorem~\ref{thm:main} and graph compactness
\citep[Theorem~1]{deBruijnErdos1951} give a finite unit-distance graph
$H$ that is not five-colorable. Set $m=|E(H)|$ and $\delta=1/m$.
...
The graph and the positive constant here are existential.
```

— `preprints/.../build/sections/introduction.tex:310–323`.

The only fixed finite unit-distance graph with an explicit edge list in `PlaneColoring/` is the Moser spindle, and it is only 4-chromatic:

```lean
def moserPoint : Fin 7 → ℂ := ![0,moserA,moserB,moserT,moserU*moserA,moserU*moserB,moserU*moserT]
def moserEdges : Finset (Fin 7 × Fin 7) :=
  {(0,1),(0,2),(1,2),(3,1),(3,2),(0,4),(0,5),(4,5),(6,4),(6,5),(3,6)}
...
lemma moser_not_three_colorable (f : Fin 7 → Fin 3)
    (hf : ∀i j,(i,j)∈moserEdges → f i≠f j) : False := by
```

— `Moser.lean:132–167`. Being 4-chromatic, it admits a 5-colouring, so it cannot witness `¬ 5-colorable`.

## What I could not find

- No explicit finite unit-distance graph that is not 5-colourable: no vertex set, no edge list, no coordinates, no cardinality, and no `SimpleGraph ℂ` anywhere in `Five.lean`.
- No formalized de Bruijn–Erdős graph-colouring compactness step in the OpenAI repo. Mathlib has `Finset.rado_selection` (`Mathlib/Combinatorics/Compactness.lean`), but no packaged graph corollary was found; extracting a finite `H` would still be nonconstructive and give no size bound.
- No Lean build or `#print axioms` verification: `elan`, `lake`, and `lean` are not installed on this machine, and the task forbids a full build without them. The toolchain is pinned to `leanprover/lean4:v4.34.1`.
- GitHub code search returned HTTP 403 rate-limit errors for the later queries `Colorable`, `not_colorable`, and `five_coloring`; I used the complete 174,524-path tree listing and the 46-file `PlaneColoring` directory instead.
- No explicit 6-chromatic unit-distance graph in the plane anywhere in the literature searched for this note; the current explicit record is the 509-vertex 5-chromatic graph (Parts, arXiv:2010.12665), consistent with the finite witness being existentially guaranteed but not constructed.
