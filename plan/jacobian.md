# Plan — planar Jacobian conjecture

Phase 0 output. Updated after `research/jacobian-enumeration.md` and
`research/jacobian-polygon-count.md`.
Status: `ADMISSIBLE_ACTIONS_EXIST` — and the admissible space is **tiny**.

## 1. Terminal object

A quadruple `(P, Q, p, q)` with

- `P, Q ∈ Q[x, y]`, `J(P, Q) ≡ 1`,
- `(deg P, deg Q) = (72, 108)` (or `(108, 72)`),
- `N(P)` and `N(Q)` one of the two configurations in GGHV 2022, Prop. 4.3,
- `p ≠ q ∈ C²`, `P(p) = P(q)`, `Q(p) = Q(q)`.

The Newton polygons are explicit:

```
config 1 (has a slope-1 edge, may be excluded by Cassou-Noguès):
  N(P) = {(0,0),(1,0),(8,14),(8,16),(0,8)}
  N(Q) = {(0,0),(2,1),(12,21),(12,24),(0,12)}

config 2 (no slope-1 edge):
  N(P) = {(0,0),(1,0),(8,14),(8,16)}
  N(Q) = {(0,0),(2,1),(12,21),(12,24)}
```

If config 1 is excluded by the no-bisector condition, **one configuration
remains**. That is the entire planar search space under the published
constraints.

**Normalisation, corrected by J1.** These are the *reduced* polygons:
`deg P = 24`, `deg Q = 36`. The `(72,108)` pair is these scaled by
`(m,n) = (3,2)`. The paper's normal form is `[P,Q] = x²`, so the system to
solve is

```
det J(P,Q) = x²
```

not `det J = 1`. On this support `det J` has **no constant term** (P has no
linear `y`-term and Q has no linear part), so `det J = 1` is infeasible by
inspection: the `(0,0)` equation reads `-1 = 0`.

Sizes, measured (`runs/jacobian_config2_size.json`):

| support | unknowns | equations | verdict |
|---|---|---|---|
| full lattice (25 + 47 points) | 72 | 92 | overdetermined by 22 |
| vertices only | 8 | 7 | small enough for a Groebner basis |

## 2. Terminal check

```
python problems/jacobian/jac.py --claim <file>.json
```

plus the degree filter `(72,108)` and the support check against config 2.

## 3. The wall

| constraint | value | consequence |
|---|---|---|
| total degree | ≤ 104 proved (Moh 1983; Nguyen 2025) | low-degree search is provably empty |
| total degree | < 125 forces `(72,108)` (GGHV 2022, **unrefereed**) | the band collapses to one degree pair |
| `(72,108)` configurations | 2 (Prop. 4.3); 1 after the literal no-bisector reading | the space is 1–2 objects, not `10^17` |
| gcd of degrees | `gcd(72,108) = 36 ≥ 16` | satisfied |
| points at infinity | exactly 2 per component | encoded by the unique top vertex |
| normal form | `m,n > 1`, `m < n`, no edge ∥ diagonal | encoded in the configurations |

The polygon enumeration also showed how badly the count depends on the normal
form: allowing equal degrees gives `~10^17` polygons, allowing `m = 1` gives
`~5×10^9`, and the published GGHV normal form gives **2**. The `10^9` gate is
therefore decided by a theorem, not by compute.

**Caveat that matters:** GGHV 2022 is an unrefereed preprint. If the `(72,108)`
reduction is wrong, the real space is larger and this plan is void. Treat the
reduction as a strong prior and re-check it before spending a week.

**Second caveat, added after the J2 audit.** `wstrinz/plane-jacobian-72-108`
claims a proof (Theorem A) that *both* configurations are empty, hence
`max(deg P, deg Q) ≥ 125`. The audit's verdict: unrefereed, internally
inconsistent (the repo's own `CURRENT_STATUS.md` still says `OPEN`, two cited
close-out documents are missing, the stated proof hash is absent from the public
history), but the argument is substantive and its checkers pass. If Theorem A
holds, **this exact target is dead**.

Consequence: `(72,108)` is no longer the only plan. The next region is the
24–27 degree pairs with `max_deg ≥ 125` still listed as open. Slot J4 is
transcribing and ranking them. J3 (Groebner on the 8-unknown system) doubles as
a test of Theorem A: a solution would refute it, `[1]` would agree with it.

**Pivot (J4 result).** The list is **24 pairs**, not 27 (the README heading is
stale; its own table and `corner_atlas.json` say 24). Only **9** have a Newton
configuration written down; the other 15 are unknown. Ranked by tractability:

| rank | `max deg` | `(deg P, deg Q)` | reduced degrees | config |
|---|---|---|---|---|
| 1 | 125 | (75,125) | (15,25) | known |
| 2 | 135 | (90,135) | (18,27) | known |
| 3 | 147 | (98,147) | (14,21) | known |

**New primary target: `(75,125)`.** It is the smallest `max_deg` on the list,
the only case that moves the 125 frontier, and it has infrastructure:

```
Delta' = {(0,0),(3,0),(4,1),(0,5)}
N(P) = 3*Delta' = {(0,0),(9,0),(12,3),(0,15)}
N(Q) = 5*Delta' = {(0,0),(15,0),(20,5),(0,25)}
[P,Q] = x^2 ; reduced degrees (15,25)
```

`(72,108)` stays as a cross-check: if J3 returns `[1]` on the vertices-only
system, that independently supports Theorem A and raises confidence that
`(75,125)` is the right place to look.

**Second result, same pattern.** The `(75,125)` corner (`N(P) = 3Δ'`,
`N(Q) = 5Δ'`) is also empty at the vertices-only level, and for a clean reason:

```
det J = 45*a1*b2*x^28*y^4 + 225*a1*b3*x^8*y^24 - 45*a2*b1*x^26*y^2
      + 300*a2*b3*x^11*y^27 - 225*a3*b1*x^14*y^14 - 300*a3*b2*x^19*y^19
min x-degree = 8
```

so the coefficient of `x^2` is identically zero and the equation reads
`-1 = 0`. Groebner returns `[1]` in 0.8 s, but no elimination is needed.

**The general lesson is the useful part: the interior lattice points are
load-bearing.** At both `(72,108)` and `(75,125)` the four-corner
specialisation fails, and it fails precisely because the low-degree monomials of
`det J` can only be produced by interior terms. A "minimal support" model is
therefore never a valid reduction — it deletes the terms that make the identity
possible.

Sizes for `(75,125)` (J5, measured):

| support | unknowns | note |
|---|---|---|
| vertices only | 8 | empty by support obstruction |
| full lattice | **445** | too large for a Groebner basis |

So the only viable method at this corner is the **edge cascade**: decompose the
Newton polygons into edges, solve the leading-form condition on the top edge,
then descend one edge at a time. That is the GGHV machinery, and it is the next
task.

The real target is the full-lattice system: 72 unknowns, 92 equations for
`(72,108)`. That is too large for a Groebner basis and needs the graded
structure instead: the top-degree block is determinantal, and the low-degree
block is where the obstruction lives.

**J3 result, for contrast.** On `(72,108)` the basis is

```
a1*b1 - 1,  a2,  a3,  b2,  b3
```

so every solution forces `a2 = a3 = b2 = b3 = 0`. The 3-parameter family is

```
P = a0 + a1*x,   Q = b0 + b1*x^2*y,   a1*b1 = 1,   det J = x^2
```

which does **not** have the required Newton polygon — the four genuine corners
are absent. It is also not a Keller map (`det J = x^2`, not constant), and its
only collision lies on the critical locus `x = 0`.

Adding the strict-support condition `a2*a3*b2*b3 ≠ 0` and saturating gives the
unit ideal:

> **no solution of `det J = x^2` on config 2 uses all four corners.**

That is a conditional negative for config 2 at the vertices-only level,
consistent with the wstrinz exclusion. It is not a proof: the full lattice (72
unknowns) is untouched.

### The reconstruction gap — the actual wall

Raised by the 2026-11-02 review, and it is decisive. The GGHV object lives in

```
L^(1) = K[x, x^-1, y]
```

and the reduction from a polynomial counterexample uses automorphisms of
`L^(1)` that are **not** automorphisms of `K[x,y]`. The implication is
one-directional:

```
counterexample  ⇒  normalized pair with [P,Q] = x²
```

It does not reverse. So a solution of the config-2 system is **not** a polynomial
map with `det J = 1`, and `jac.py` rejects it — correctly, because its
determinant is `x²`, not a constant.

Until an inverse normalisation is written and checked, **none** of the config-2
work can reach the terminal object. That makes reconstruction the top task,
ahead of any further elimination. If the inverse cannot be defined from the
published data, the whole Jacobian workstream is off-core and the plan must be
rewritten.

The review also corrected the "interior points are load-bearing" inference: the
vertices-only failure is a truncation artefact, not a theorem that minimal
supports are always invalid. The support of `det J` is the Minkowski difference
of the supports, and any truncation must be justified by the Newton-polygon
edge structure. That is what the cascade is for.

### Reconstruction verdict (measured): `PARTIAL`

`problems/jacobian/reconstruct.py` inverts the composite as a **rational** map
and restores `det J = 1` before clearing denominators. It is not a map to
`K[x,y]`: the inverse swap introduces `y^-1`.

For the J3 representative `P = x`, `Q = x^2 y`:

```
P0 = 1/y,  Q0 = x*y^2
cleared: P' = 1,  Q' = x*y^2,   det J(P',Q') = 0
```

No polynomial clearing helps. A polynomial representative is `(f*P0, g*Q0)` with
`y | f` and `y^2 | g`, and the determinant is then divisible by `y` unless it
vanishes.

The published data do not supply:

1. the exact shear constants `(c2,c3,c4)` or the branch of the proof;
2. the vanishing conditions for the negative powers of `y` (down to `y^-8` from
   the `(8,14),(8,16)` vertices);
3. a proof that a config-2 solution lies in the image `A(K[x,y])`.

**What this blocks, and what it does not.** The implication runs

```
counterexample  ⇒  normalized pair with [P,Q] = x²
```

so:

- **The success direction is blocked.** A normalized *solution* cannot be turned
  into a counterexample without an inverse. That is all `PARTIAL` means.
- **The emptiness direction is unaffected.** If the normalized system is empty,
  there is no counterexample of that degree. No reconstruction is needed,
  because the implication is used in the direction in which it actually holds.

Every computation run so far is an emptiness result, and emptiness is exactly
what a **degree-bound** proof consumes. The branch was not blocked; it was
mis-aimed.

The live target is therefore:

> **raise the planar degree bound from 104 to 125, and past it, by proving each
> normalized corner empty.**

That is a theorem we can finish with the tools already built, and it does not
require the inverse normalisation. A counterexample remains the long shot, gated
on reconstruction; a degree bound is the deliverable that is actually reachable.

The `(75,125)` cascade is now on-core again: an empty cascade at that corner is a
degree-bound step, not a dead end.

## 4. Load-bearing decomposition

1. **Write the coefficient system for config 2.** Unknown coefficients on the
   two supports; `det J = 1` is a polynomial identity, so it is a finite system
   of polynomial equations in the coefficients. *Delete it and there is no
   candidate at all.*
2. **Solve it over Q (or prove no rational solution).** *Delete it and the
   polygon is a shape, not a map.*
3. **Reject triangular solutions.** A composition of elementary maps is
   invertible by construction. *Delete it and a shear gets reported as a
   counterexample — exactly what happened in the degree-3 pilot.*
4. **Collision search + exact verification.** *Delete it and there is no
   terminal object.*

Item 1 is the whole job. Items 2–4 are standard.

## 5. Admissible actions

1. ~~Build the config-2 coefficient system.~~ **Done.** 72 unknowns / 92
   equations on the full lattice; 8 / 7 on the vertices. Constructed by
   `problems/jacobian/config2_system.py`, reconstructed and checked exactly.
2. **Solve the vertices-only system.** 8 quadratic equations, 8 unknowns. A
   Groebner basis either returns `[1]` (no solution — decisive) or a variety.
   This is the only part small enough to settle by elimination.
3. **Decide the full-lattice system without eliminating it.** The top-degree
   block is determinantal: the leading forms of `P` and `Q` must share a factor,
   or the top equations are inconsistent. That check is cheap and settles
   whether the generic support can work at all.
4. **Check `wstrinz/plane-jacobian-72-108`** (slot J2). If both configurations
   are excluded and the argument holds, this whole branch closes.
5. Only after 2–4: collision search and exact verification.

There is no broad sweep. There is no degree-105–110 enumeration.

## 6. Kill criteria

- The GGHV `(72,108)` reduction cannot be re-derived or is refuted → stop; the
  plan is void.
- `wstrinz/plane-jacobian-72-108` survives scrutiny → stop; both configurations
  are excluded.
- The config-2 system has only triangular solutions → stop. That is consistent
  with the conjecture and is a real (if conditional) result.

## What the research changed

Before: "enumerate degree-3 maps and look for collisions." That box is provably
empty.

Then: "enumerate admissible Newton polygons at degree ≥105." That produced
`10^17` polygons under weak normal forms, and **2** under the published one.

Now: "solve one coefficient system attached to one explicit pentagon/quadrilateral
pair." The space went from `10^17` to `2` in two research passes, without
writing a line of search code.
