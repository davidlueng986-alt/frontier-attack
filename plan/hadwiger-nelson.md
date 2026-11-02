# Plan — Hadwiger–Nelson (6 vs 7)

Phase 0 output. Written after `research/hn-gadgets.md`.
Status: `ADMISSIBLE_ACTIONS_EXIST`.

## 1. Terminal object

**Target 1 (the question):** an explicit finite unit-distance graph `G ⊂ R²`
with `χ(G) = 7`, plus a machine-checkable proof of non-6-colourability
(DRAT/LRAT trace or equivalent).

**Target 0 (calibration, guaranteed to exist, no witness known):** an explicit
finite unit-distance graph with `χ(G) = 6`.

Target 0 does **not** advance the 6-vs-7 question — χ(ℝ²) ≥ 6 is already proved
(OpenAI, Sept 2026). It is the calibration object: the hardest construction that
is known to exist but has never been written down. A pipeline that cannot find
a 6 will not find a 7.

## 2. Terminal check

```
python problems/hadwiger_nelson/hn.py --verify <graph>.json     # exact unit distances
python problems/hadwiger_nelson/chi.py --lower-bound 7 <graph>  # DRAT / UNSAT
```

Two conditions, both machine-checkable:

1. every edge has squared length exactly 1 (symbolic, not float);
2. the graph is not 6-colourable (SAT returns UNSAT, trace checked by `drat-trim`).

Condition 2 is the expensive one. `hn.py`'s DSATUR is exact but only for small
graphs; at 500+ vertices the lower bound must come from a SAT proof, not from
branch and bound.

## 3. The wall

| fact | value | consequence |
|---|---|---|
| record 5-chromatic | Parts 509 v / 2,442 e | the seed graph; coordinates at `vasnesterov/HadwigerNelson`, `vtx/509_parts.vtx` |
| next best | Heule 510, 517, 529; de Grey 1,581 | all 5-chromatic; no 6 exists explicitly |
| 6-chromatic in the plane | **none known** | existence is non-constructive (de Bruijn–Erdős + the Sept 2026 not-5-colourable proof) |
| 7-chromatic | would settle χ(ℝ²) = 7 | no route known |
| de Grey's "H" | 7 v, 12 e wheel | forces a **monochromatic √3-triangle** in a 5-colouring — a 3-vertex pattern, *not* a two-terminal equality |
| diamond | 4 v, 5 e | forces its two tips **equal** in any 3-colouring |
| spindle operation | glue two copies at a vertex, rotate until a second pair is unit distance | if that pair is forced equal, the union is not k-colourable |

**This is why the previous attempt saturated at 4.** `glue2.py` performed rigid
gluing but never used the forced-equality step. Gluing without the equality
constraint is just union; the chromatic number cannot rise because the
constraint that would create the contradiction is never imposed.

## 4. Load-bearing decomposition

1. **Implement forced equality, not just gluing.** The diamond's tips are equal
   in every 3-colouring. A composition rule must carry that fact. *Delete it and
   the search is a union of 4-colourable graphs, which is 4-colourable.*
2. **Implement the spindle operation with the second-pair distance check.** Two
   copies, one shared vertex, rotate until a designated pair is exactly unit.
   *Delete it and no new constraint is created.*
3. **SAT lower bound with a DRAT trace.** At 500+ vertices, DSATUR is not a
   proof and is too slow. *Delete it and a claimed χ=7 has no certificate.*
4. **Exact coordinate verification.** *Delete it and near-misses pass as unit
   distances.*

Items 1–2 are the work. 3 is off-the-shelf (`cnp-5`, `CNP-SAT`). 4 already
exists in `hn.py`.

## 5. Admissible actions

**H3 result (measured, 18,889 pairs).** In the 509 graph at k = 5 there are
**zero forced-equal pairs**. The 721 diamonds in it force their apices equal at
k = 3 — as the gadget theory says — but those apices are **free at k = 4 and at
k = 5**.

H1 then confirmed the consequence directly: the diamond composition reproduces
the Moser spindle at χ = 4 and **saturates there**, across every composition
tried (`best_chi = 4`).

Consequence: equality gadgets are a **3-colouring tool**. Composing diamonds
therefore cannot produce anything above 4 colours, which is exactly what
`glue2.py` observed. Any plan built on `forced_relation == "equal"` is
calibrated to the wrong k.

### The levels are indexed, and we were measuring at the wrong one

| level | primitive | forces | status |
|---|---|---|---|
| k=3 | diamond | apices equal | known → Moser spindle, χ=4 |
| k=4 | H, and L = 52 copies of H | a monochromatic √3-triangle in some copy | known → χ=5 |
| k=5 | **virtual edge** | endpoints differ in **every** 5-colouring | **missing** |

H is a *4-colouring* gadget: it has four 4-colourings, two with a monochromatic
triple, and 52 copies (L) force the triple. At 5 colours H has a colouring with
no monochromatic √3-triangle. The earlier statement that a "7v/12e wheel forces
a monochromatic triangle at 5 colours" was a level error — three files in this
repo gave three different levels for the same object.

The 509 graph is 5-chromatic, i.e. a k=4 obstruction. Measuring its relations at
k=5 was therefore the wrong level, which is why the result was null.

**A size bound kills the 7-vertex approach outright.** A gadget that forces a
monochromatic triple in a 5-colouring needs at least **11 vertices**: every
7-vertex unit-distance graph is 4-colourable, and any 4-colouring of 7 vertices
refines to a 5-colouring with all colour classes of size ≤ 2 (5 classes of size
2 cover 10 vertices). A class of size ≤ 2 is not a monochromatic triple.

So the 7-vertex H cannot be a k=5 gadget at all, for a counting reason, and any
pattern-forcing search at k=5 must start at n ≥ 11. `six.py` found this and
proved it, rather than merely failing to find a gadget.

**The route to 6 is a level-5 virtual edge**: a unit-distance gadget forcing its
two endpoints to differ in every 5-colouring. Compose it with a 6-chromatic
*two-distance* graph — Parts' 16-vertex graph (distances 1 and φ), the 33-vertex
graph (distances 1 and √3), or Exoo–Ismailescu's 426-vertex graph — replacing the
non-unit edge by the gadget. That yields a 6-chromatic unit-distance graph.

**χ=7 is out of reach by composition.** Any 7-chromatic unit-distance graph has
at least **6,993 vertices** (Parts 2020, Table 2, cited). A 509-plus-gadget
composition is ~516 vertices and is provably 6-colourable. The old step 5 was
dead on size, not on compute.

**Extraction from the Lean proof is ruled out.** `research/hn-lean-witness.md`
checked the public source: `lean/OAI/Geometry/PlaneColoring/Five.lean` is 131
lines, no `sorry`, and proves the infinitary statement

```lean
theorem no_proper_five_coloring : ¬ ∃ c : ℂ → Fin 5, ProperColoring 5 c
```

It never defines a finite vertex set (`Fin N`), an edge list, or coordinates.
The paper says so itself: *"The graph and the positive constant here are
existential."* So the witness has to be **searched for**, not extracted.

The one remaining extraction lead is `PlanarCovers.compact_open_cover_impossible`
(`PlanarCovers.lean:374`). If that lemma is finitary — a finite open cover that
cannot exist — then the obstruction it encodes may be instantiable as an
explicit finite point set. That is a research question, not an assumption.

Ordered by cost, not by appeal:

1. ~~Read `PlanarCovers.lean` and decide whether the obstruction is finitary.~~
   **Resolved negative.** The obstruction is a finite *open cover*, not a finite
   graph: `elim_finite_subcover` (`PlanarCovers.lean:381`) produces a `Finset ι`
   of cover elements, and the closing lemma `finite_partition_impossible`
   (`:302`) is a Lebesgue covering-dimension argument that never mentions 5
   colours. There is no vertex set, no edge list, no coordinates, and no size
   bound on a graph. `s.card ≥ 25` is a bound on **cover elements**, not
   vertices — do not confuse the two. Extraction is dead.
2. **Load `509_parts.vtx` and verify it with `hn.py`.** Confirms our exact
   arithmetic on the current record and gives us a seed. Cheap, and it tests the
   gate on a real object rather than a 7-vertex toy.
3. **Re-run gluing with forced-equality bookkeeping.** The previous run had the
   right geometry and the wrong constraint. This is the smallest change that
   could move χ from 4 to 5.
4. **SAT-trim around the 509 seed** (Heule's pipeline) to look for a 6. This is
   calibration, not the target.
5. **Only if 1–4 produce a 6: push for a 7** by adding a pattern-forcing gadget
   for a monochromatic √3-triangle, which is the primitive de Grey actually used.

## 6. Kill criteria

- `509_parts.vtx` fails exact verification → the coordinate source is wrong;
  stop and re-source before any search.
- Forced-equality gluing still saturates at 4 after 50 trials → the composition
  rule is wrong, not the compute budget.
- SAT-trim at 509 produces no 6 within a fixed budget → this is expected; it
  does not justify scaling. Record it and move to a different seed family.

## What the research changed

Before: "find colour-forcing gadgets with two terminals (equal / different)."
After: the primitive is a **3-vertex pattern** (monochromatic √3-triangle), and
equality is a *property of the diamond*, not a general two-terminal gadget.

Before: "extract the finite witness from the Lean proof."
After: **impossible.** The 2026 theorem is proved by a covering-dimension
argument; it is not a finite-graph proof and never names 5 colours. The witness
exists only by de Bruijn–Erdős, non-constructively, with no size bound.

Consequence: an explicit 6-chromatic unit-distance graph would be genuinely new
— not a re-reading of an existing proof. That makes the search target worth the
compute, and it makes "we could not find one" an honest result rather than a
failure to transcribe.
