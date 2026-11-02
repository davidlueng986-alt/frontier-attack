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

1. **Build the config-2 coefficient system.** Five vertices in `N(P)`, four in
   `N(Q)`; the supports are fixed, so this is mechanical. Output: a polynomial
   system in `< 30` unknowns.
2. **Check the 2026 repo `wstrinz/plane-jacobian-72-108`**, which claims a
   conditional structural exclusion of both configurations. If that claim is
   correct and the argument survives review, this branch is closed and the
   correct output is `NO_ADMISSIBLE_ACTION` for the whole planar route.
3. **Only after 1–2:** solve the system, reject triangular solutions, test
   collisions.

There is no broad sweep. There is no degree-105–110 enumeration. That work was
done and it returned "two objects".

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
