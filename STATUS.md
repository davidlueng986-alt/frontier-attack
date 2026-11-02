# STATUS

One page. What we are doing, where we are, and what is running. Updated every
time the pipeline moves.

## Mission

Solve three open problems. Not "make progress on" — produce the object that
settles each one.

| # | Problem | Terminal object | Terminal check |
|---|---|---|---|
| 1 | Hadwiger–Nelson | explicit finite unit-distance graph with χ = 7 | exact unit distances + DRAT UNSAT for 6-colouring |
| 2 | Planar Jacobian | `(P, Q, p, q)`, `J(P,Q) ≡ 1`, degrees `(72,108)`, `p ≠ q`, same image | `jac.py --claim` → `COUNTEREXAMPLE` |
| 3 | BSD | `(E, r, k)`, `r` proved, `k` proved, `r ≠ k` | `ellrank` interval collapsed + certified `k` |

## Honest odds, revised after Phase 1

| Problem | Odds 24 mo (now) | Was | Why it moved |
|---|---|---|---|
| HN χ=7 | **0.05** | 0.40 | no explicit χ=6 exists; the record is 509 for χ=5; the 2026 proof is a covering-dimension argument and yields no finite witness |
| Jacobian `(72,108)` | **0.10** | 0.20 | the space collapsed to 2 configurations, but they are exactly the hard ones; an unrefereed 2026 repo claims both are excluded |
| BSD counterexample | **0.05** | 0.35 | analytic rank ≥ 4 is uncertifiable; LMFDB is complete below 500,000; no positive instance has ever been observed |

The earlier numbers were guesses made before research. These are guesses made
after it. They are all much lower, and that is the correct direction.

## Where each problem stands

### 1. Hadwiger–Nelson

- ✅ Record graph verified independently: 509 vertices, 2,442 edges, worst
  `|d²−1| = 1.5e-49`, 4-colouring UNSAT in 134 s. `tools/verify_509.py`.
- ✅ The 2026 "not 5-colourable" Lean proof is public but **infinitary**; it
  constructs no finite witness. The finite graph is existential.
- ✅ The covering-dimension core is not a graph argument and cannot be
  instantiated into one.
- ❌ No explicit 6-chromatic graph. No explicit 7-chromatic graph.
- **Next:** the diamond's forced-equality relation is the missing primitive.
  `glue2.py` glued geometry without it and saturated at χ=4.

### 2. Planar Jacobian

- ✅ Degree ≤ 104 is a theorem, so the degree-3 pilot was provably empty.
- ✅ Below 125 only the degree pair `(72,108)` survives (GGHV 2022,
  **unrefereed**).
- ✅ That leaves **2 Newton configurations**, 1 under the literal no-bisector
  reading.
- ❌ The coefficient system has not been written.
- ❌ An unrefereed 2026 repo (`wstrinz/plane-jacobian-72-108`) claims both are
  excluded — must be checked before any solve.
- **Next:** write the config-2 system and size it.

### 3. BSD

- ✅ A working PARI-free analytic-rank estimator: ranks 0/1/2/3 correct on four
  curves; `L(1)` matches LMFDB to 29 digits.
- ✅ The obvious 20k truncation **fails** on rank 1 and rank 3 — recorded so
  nobody repeats it.
- ❌ The estimator is numerical, not certified. It cannot supply the `k` in the
  terminal object.
- ❌ No algebraic-rank backend without PARI/GP.
- **Next:** certify the AFE ladder with interval arithmetic; then a coverage
  sweep above `N = 500,000`.

## Pipeline

Max 10 concurrent agents; target 5–7. A slot is refilled as soon as any agent
returns — the pipeline does not wait for the slowest member.

| slot | problem | task | output | gate |
|---|---|---|---|---|
| J1 ✅ | Jacobian | config-2 system sized | `problems/jacobian/config2_system.py` | 72 unk / 92 eq (full), 8/7 (vertices) |
| J2 | Jacobian | audit `wstrinz/plane-jacobian-72-108` | `research/jacobian-wstrinz-audit.md` | states exclude yes/no with file:line |
| J3 | Jacobian | Groebner on the 8-unknown vertices system | `problems/jacobian/solve_vertices.py` | `[1]` or a variety, ≤700 s |
| H1 | HN | equality-gadget class + composition + search | `problems/hadwiger_nelson/gadgets.py` | reproduces Moser spindle from two diamonds |
| H2 | HN | SAT relation prober (equal / different / free in k-colourings) | `problems/hadwiger_nelson/relation.py` | agrees with brute force on ≤10-vertex graphs |
| B1 | BSD | certify the AFE ladder for k ≤ 1 | `problems/bsd/certified_rank.py` | interval bound contains the known value |

Refill rule: first return → immediately launch the next task for that problem,
not a copy of the finished one.

## Definition of done

A problem is done when the terminal checker accepts the object **and** the
cross-verifier, implemented by a different method, agrees. Everything else is
`NO_ADMISSIBLE_ACTION` or a negative result, both of which get written down.

## What we will not do

- Report a lemma as progress.
- Call a numerical rank a rank.
- Call a 6-chromatic graph a solution to 6-vs-7.
- Scale agents past `1/p` before `p > 0` has been measured.
