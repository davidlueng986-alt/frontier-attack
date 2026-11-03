# STATUS

## Current root status and closed-attempt plan (2026-11-03)

All three root problems remain OPEN in this project. The active root definitions, direct attempts,
and evidence requirements are in [docs/CLOSED_ATTEMPTS.md](docs/CLOSED_ATTEMPTS.md).
The prior snapshot below is historical: its pipeline slots are not currently running, its probabilities
are not validated estimates, and its checker shortcuts are not sufficient for a solved claim.
In particular, `tools/accept.py` does not yet actually replay HN proof files or BSD rank certificates;
the legacy degree125 filter is not a mathematical definition of a planar counterexample.
Old statements such as "CLOSED as a counterexample hunt" close only a method, not BSD itself.

Calibration executed: all affine planes through the three listed source collision pairs,
with all linear output projections, cannot yield a planar Keller map. Exact certificates and an
independent chain-rule verifier are in `runs/jc_affine_projection_certificate.json` and
`attempts/verify_jc_affine_projection.py`. The degree<=7 family was already excluded by known
planar low-degree results and should have been screened out first. It is retained as a verifier
calibration, not new main-goal progress or planar JC CLOSED. Quadratic embeddings are likewise
not an admissible new lead (output degree<=14); see the degree-screened next-task contract.

## Historical snapshot

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
- ✅ The config-2 and `(75,125)` vertices-only systems are both empty; for
  `(75,125)` by a support obstruction (min x-degree 8).
- ❌ **Reconstruction is `PARTIAL`.** The GGHV inverse exists only as a rational
  map; clearing denominators kills `det J = 1`. So a normalized solution is
  necessary, not sufficient, and cannot become a counterexample.
- **Status: BLOCKED at the published interface.** Missing: shear constants,
  vanishing conditions for negative powers of `y`, and a proof that a solution
  lies in `A(K[x,y])`.

### 3. BSD

- ✅ A working PARI-free analytic-rank estimator: ranks 0/1/2/3 correct on four
  curves; `L(1)` matches LMFDB to 29 digits.
- ✅ The obvious 20k truncation **fails** on rank 1 and rank 3 — recorded so
  nobody repeats it.
- ✅ Coverage check: 25 curves with conductor in `[500k, 520k]`, analytic ranks
  0–3, all reproduced exactly, 0 mismatches, 0 parity inconsistencies.
- ✅ **GZK closes the counterexample route.** Gross–Zagier + Kolyvagin: analytic
  rank ≤ 1 forces `r = k`. So any counterexample needs `k ≥ 2`.
- ❌ No algebraic-rank backend on this machine (no PARI, no C compiler).
- ❌ No rigorous `k ≥ 2` certificate (needs zero-counting, not an AFE interval).
- **Status: CLOSED as a counterexample hunt.** The remaining work is
  infrastructure, not a solution. Keep the tools; stop spending slots here.

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

The contract is `docs/DEFINITION_OF_SOLVED.md`, enforced by
`python tools/accept.py --problem <p> --claim <file>`.

A claim is accepted only if `accept.py` exits 0 **and** a second verifier,
written by a different method, exits 0. `accept.py` is deliberately stricter
than the individual checkers: the real 2026 C³ counterexample passes
`jac.verify_map` and is **rejected** here, because its degree is 7 and the planar
bound is 125. That is the gate working — it separates *true* from *the thing*.

Everything else is `NO_ADMISSIBLE_ACTION` or a negative finding, both of which
get written down.

## Workstreams, not single tasks

Each slot owns a **queue of 3 tasks**, executed one at a time. The agent reports
after each and continues without waiting for a new assignment. A slot is only
idle when its queue is empty; that is when it is refilled.

## What we will not do

- Report a lemma as progress.
- Call a numerical rank a rank.
- Call a 6-chromatic graph a solution to 6-vs-7.
- Scale agents past `1/p` before `p > 0` has been measured.
