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
