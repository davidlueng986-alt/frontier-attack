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

Ordered by cost, not by appeal:

1. **Load `509_parts.vtx` and verify it with `hn.py`.** Confirms our exact
   arithmetic on the current record and gives us a seed. Cheap, and it tests the
   gate on a real object rather than a 7-vertex toy.
2. **Re-run gluing with forced-equality bookkeeping.** The previous run had the
   right geometry and the wrong constraint. This is the smallest change that
   could move χ from 4 to 5.
3. **SAT-trim around the 509 seed** (Heule's pipeline) to look for a 6. This is
   calibration, not the target.
4. **Only if 1–3 produce a 6: push for a 7** by adding a pattern-forcing gadget
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

The earlier plan was looking for the wrong object. That is exactly the failure
the Think phase is supposed to prevent, and it was caught by Phase 1 for the
price of one research agent.
