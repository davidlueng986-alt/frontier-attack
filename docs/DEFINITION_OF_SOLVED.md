# What counts as solved

This file is the contract. If a claim does not satisfy it, it is not a solution,
regardless of how interesting it is.

## The rule

A problem is solved iff **both** of these pass:

1. `python tools/accept.py --problem <p> --claim <file>` exits 0.
2. A second verifier, implemented by a **different method**, exits 0.

Verifier A may be the repo's checker. Verifier B must not import A's code.
Two runs of the same program is not a cross-check.

A solution also ships an **evidence bundle**:

```
claim.json          the object
check.log           exact command + stdout + exit code of A
check_b.log         exact command + stdout + exit code of B
env.txt             python version, package versions, commit hash
```

No bundle, no claim. This is not bureaucracy: the last three "results" in this
repo were all correct and all irrelevant, and a bundle is what made that visible.

## Hadwiger–Nelson: `chi(R^2) = 7`

**Claim object.** A finite set of points with algebraic coordinates and a list of
unit-distance edges.

**Accept:**

| # | condition | check |
|---|---|---|
| A1 | every declared edge has `|d² − 1| = 0` exactly, at ≥50 digits | `tools/verify_509.py` geometry path |
| A2 | no undeclared pair at distance 1 | same |
| A3 | the graph is **not 6-colourable** | SAT `k=6` returns UNSAT |
| A4 | the UNSAT result has a checkable DRAT/LRAT trace | `drat-trim` exits 0 |
| A5 | the graph is finite and the coordinates are explicit | schema check |

**Reject:**

- A 6-chromatic graph. It proves `chi >= 6`, which was already known in
  September 2026. It is a calibration object, not a solution.
- A 7-chromatic graph with no UNSAT certificate. A colouring heuristic saying
  "looks like 7" is not a proof.
- Floating-point coordinates with tolerance `1e-9`. Rounding is not geometry.
- A Lean proof of `chi >= 6` repackaged as `chi >= 7`.

**Current target.** No explicit 6-chromatic graph is known. That is Target 0
(calibration). Target 1 is 7.

## Planar Jacobian: a counterexample to JC

**Claim object.** `(P, Q, p, q)` with `P, Q ∈ Q[x,y]`, `p ≠ q ∈ C²`.

**Accept:**

| # | condition | check |
|---|---|---|
| J1 | `det J(P,Q)` is a non-zero constant, computed symbolically | `jac.py` |
| J2 | `P(p) = P(q)` and `Q(p) = Q(q)` exactly | `jac.py` |
| J3 | `max(deg P, deg Q) ≥ 125` | degree filter |
| J4 | the pair is **not triangular** | `tools/accept.py` rejects `(x + h(y), y)` and `(x, y + h(x))` forms |
| J5 | the Newton polygon is one of the open configurations | support filter |

**Reject:**

- A constant-Jacobian map with no collision. That is a candidate, not a
  counterexample.
- A collision found numerically. Exact substitution only.
- Degree ≤ 104. Proved invertible; the box is empty by theorem.
- A shear `(x, y + h(x))`. Invertible by construction. This is exactly what the
  degree-3 pilot produced, and it must never be reported again.

### Jacobian, secondary target: the degree bound

A different and finishable claim:

> every planar Keller map of total degree ≤ D is invertible, for a new D.

**Accept:**

| # | condition |
|---|---|
| D1 | every counterexample with `max(deg) < D` must normalise into one of finitely many corners (GGHV) |
| D2 | each of those corner systems is proved empty |
| D3 | the emptiness proof is machine-checked or symbolically certified |

**No reconstruction is required.** The implication is used only in the direction
`counterexample ⇒ normalized pair`, so emptiness of the normalized system
suffices.

This does not solve JC. It raises the bound. It is a theorem, and it is the
target the current tooling can actually reach. It must be labelled as a degree
bound, never as a solution to the conjecture.

## BSD: a counterexample

**Claim object.** `(E, r, k)` with `r` the algebraic rank and `k` the analytic
rank, both **proved**.

**Accept:**

| # | condition | check |
|---|---|---|
| B1 | discriminant ≠ 0 | `bsd.py` |
| B2 | algebraic rank proved: 2-descent interval `r1 = r2 = r` | `ellrank` output attached |
| B3 | analytic rank proved: certified order of vanishing, `k ≤ 3` | interval certificate |
| B4 | `r ≠ k` | arithmetic |
| B5 | root-number parity consistent | `bsd.py` |

**Reject:**

- A numerical near-zero of `L`. `ellL1` on 389a1 returns `−5.4e−29` at default
  precision and its analytic rank is 2.
- A 2-descent interval with `r1 < r2`. A gap is not a rank.
- `k ≥ 4`. The analytic order is not proved for any curve, per the Sage
  documentation cited in `research/bsd-rank-tooling.md`.
- A database field taken on faith. LMFDB is curated; that is the point of it.

## What "solved" is not

- It is **not** peer review, and it is **not** Clay Institute acceptance. The
  correct phrasing is "mathematically settled, pending community verification".
- It is **not** a Lean build of a related statement.
- It is **not** "no counterexample found". That is a negative result, and it goes
  in `findings/`, not in a solution claim.
- It is **not** a 6-chromatic graph, a triangular map, or a numerical rank.

## Why this file exists

Every failure in this repo has been a correct artifact that did not meet this
contract:

| artifact | correct | meets contract | why not |
|---|---|---|---|
| 3/2000 constant-Jacobian maps | yes | no | all triangular |
| 20/20 four-chromatic gluing | yes | no | target was 7 |
| 50-curve LMFDB scan | yes | no | curation forbids a mismatch |
| AFE rank estimator | yes | no | not certified |

The checker answers "is it true". This file answers "is it the thing".
