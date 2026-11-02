# 2026-11-02 — Jacobian, cubic sparse family: everything collapses to a shear

Agent: `gen-jac-cubic-01`. Checker: `problems/jacobian/jac.py`.
Run: `python problems/jacobian/cubic_search.py --samples 2000`.
Artifact: `runs/jacobian_cubic.json`.

## Result

| quantity | value |
|---|---|
| sampled maps | 2,000 |
| constant Jacobian determinant | 3 (0.15%) |
| counterexamples | **0** |
| all survivors triangular | **3 / 3** |

The three survivors:

```
F = (x, y + 3x^3 - x^2)        det J = 1
F = (x, y - x^3 - 2x^2)        det J = 1
F = (x, y + 2x^2)              det J = 1
```

All have the form `(x, y + h(x))`, i.e. a shear. Each is invertible in one
line: `(u, v - h(u))`.

## Why this is the interesting part

The family was specified as `f = x + (cubic terms)`, `g = y + (cubic terms)`.
But the only way to keep det J constant is for all of `f`'s non-linear terms to
vanish — otherwise `f_x` picks up an `x`-dependence that cannot be cancelled by
`g_y = 1 + ...` without forcing `g` to depend on `x`, which then re-introduces
a determinant term.

So the search space is **degenerate by construction**: 100% of its solutions
lie in the triangular subvariety, and every triangular map is invertible. The
agent cannot find a counterexample here for the same reason it cannot find one
by enumerating the identity map.

## What to change in the next iteration

1. Do not fix the linear part to `(x, y)`. Normalise by the *leading
   homogeneous part* instead, allowing a non-identity linear part.
2. Require both components to carry non-linear terms: reject any candidate with
   `f - x` or `g - y` identically zero before the det test.
3. Parameterise by Newton polygon and solve `det J = const` as a polynomial
   system in the coefficients, rather than sampling coefficient boxes.
4. The Cartographer must report the dimension of the triangular subvariety and
   subtract it from the space size. A space that is 99.85% empty and 0.15%
   trivial has an effective size of zero.

## Status

Negative result. No claim about the planar Jacobian conjecture is made or
implied. This only says: this particular family is the wrong place to look, and
now we know why.

## Addendum (after Phase 1 research)

`research/jacobian-enumeration.md` establishes that the planar case is **proved
for total degree ≤ 104** (Moh 1983 gives ≤ 100; Nguyen 2025 gives ≤ 104; the
unrefereed GGHV preprint gives < 125 except (72, 108)).

So this box was not merely degenerate — it was **provably empty by theorem**.
The 0 counterexamples were guaranteed before the first sample was drawn. The
correct search region is total degree ≥ 105 with `gcd(deg P, deg Q) ≥ 16` and a
constrained trapezoidal Newton polygon; see `plan/jacobian.md`.

This is a better outcome than the original framing: the negative result now has
a theorem behind it instead of an empirical observation.
