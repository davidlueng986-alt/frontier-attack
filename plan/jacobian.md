# Plan — planar Jacobian conjecture

Phase 0 output. Written after `research/jacobian-enumeration.md`.
Status: `ADMISSIBLE_ACTIONS_EXIST` (conditional on the polygon count, see §5).

## 1. Terminal object

A quadruple `(P, Q, p, q)` with

- `P, Q ∈ Q[x, y]`, `J(P, Q) ≡ 1` (constant non-zero determinant),
- `deg_total(P, Q) ≥ 105`,
- `gcd(deg P, deg Q) ≥ 16`,
- `p ≠ q ∈ C²` and `P(p) = P(q)`, `Q(p) = Q(q)`.

That is the whole object. Not a proof, not a reduction, not a family: five
explicit things.

## 2. Terminal check

```
python problems/jacobian/jac.py --claim <file>.json
```

which runs, in order:

1. `verify_map([P,Q],[x,y],[p,q]) == COUNTEREXAMPLE` (exact symbolic, det must
   be a non-zero constant and the two images must coincide);
2. `deg_total ≥ 105`;
3. `gcd(deg P, deg Q) ≥ 16`.

All four must hold. The checker already certifies a real object of this shape in
dimension 3 (see `alpoge_counterexample()` in `jac.py`), so the acceptance path
is exercised, not hypothetical.

## 3. The wall

| constraint | value | consequence |
|---|---|---|
| total degree | ≤ 100 proved (Moh 1983); ≤ 104 (Nguyen 2025) | any counterexample has degree ≥ 105 |
| total degree | < 125 except (72, 108) (GGHV 2022, unrefereed) | narrows the band further |
| gcd of component degrees | ≥ 16 | both degrees are large; rules out most small supports |
| points at infinity | each component ≤ 2; one ⇒ invertible | forces exactly 2 per component |
| Newton polygons | similar trapezoids, top vertex (m, n), n > m > 1, no edge ∥ diagonal | the space is a thin arithmetic stratum, not "all polynomials" |

The important consequence: **our earlier degree-3 enumeration was not merely
degenerate, it was provably empty.** Everything below degree 105 is settled by
theorem. Any future search that samples low-degree coefficients is wasting
compute by construction.

## 4. Load-bearing decomposition

1. **Enumerate admissible Newton polygons.** Trapezoid, gcd ≥ 16, total degree
   ≥ 105, no edge parallel to the diagonal. *Delete it and the space has no
   index — coefficient sampling becomes meaningless.*
2. **Solve `J(P,Q) ≡ 1` on a fixed support.** A polynomial system in the
   coefficients, linear in the top-degree unknowns. *Delete it and there is no
   candidate pair.*
3. **Collision search** on each surviving pair. *Delete it and there is no
   terminal object, only a Jacobian pair.*
4. **Exact verification** of det and of both images. *Delete it and there is no
   proof.*

Items 1 and 2 are the work. Items 3 and 4 are cheap.

## 5. Admissible actions

- For each polygon, build the coefficient system and solve it.
- **Reject triangular solutions before the collision search.** A composition of
  elementary maps is invertible by construction and cannot be a counterexample.
  This is the filter whose absence produced the 3/3 shear result in
  `findings/2026-11-02-jacobian-cubic.md`.
- Only pairs that survive 1–3 are admissible.

**Gate:** if the number of admissible polygons at total degree 105–110 with
`gcd ≥ 16` exceeds `10^9`, write `NO_ADMISSIBLE_ACTION` and stop. The next move
is then a stronger invariant, not more compute.

## 6. Kill criteria

- Polygon count at degree 105–110 exceeds `10^9` → stop.
- 20 consecutive polygons yield only triangular solutions → stop, re-derive the
  support assumptions.
- 100 admissible pairs yield 0 collisions → stop. That is consistent with the
  conjecture; it is not evidence for it, and it is not progress.

## What the research changed

Before: "enumerate degree-3 maps and look for collisions."
After: "degree ≤ 104 is a theorem; the only admissible region is degree ≥ 105
with gcd ≥ 16 and a constrained Newton polygon."

One research pass deleted ~100% of the previous search space and replaced it
with a much smaller, structured one. That is Phase 1 doing its job.
