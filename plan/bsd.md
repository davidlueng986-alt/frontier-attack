# Plan — BSD

Phase 0 output. Written after `research/bsd-rank-tooling.md`.
Status: **`NO_ADMISSIBLE_ACTION` for "find a counterexample".**
Status: `ADMISSIBLE` for "extend certified coverage", which is a different and
much smaller claim.

## 1. Terminal object

`(E, r, k)` where

- `E/Q` is an elliptic curve, non-singular;
- `r` is the algebraic rank, **proved**: `ellrank` returns `r1 = r2 = r`;
- `k` is the analytic rank, **proved**: `ord_{s=1} L(E,s) = k` with an explicit
  numerical error bound, and `k ≤ 3`;
- `r ≠ k`.

All four. A numerical near-zero of `L` is not `k`. A 2-descent interval with
`r1 < r2` is not `r`.

## 2. Terminal check

```
gp -q < claim.gp     # ellrank -> r1,r2 ; ellanalyticrank -> k, L^(k)(1)
python problems/bsd/bsd.py --claim claim.json
```

Accept only if `r1 == r2`, the analytic order is certified at the stated
precision, and `bsd.check_claim` returns `COUNTEREXAMPLE_CANDIDATE` (not
`INCOMPLETE`).

## 3. The wall

| fact | value | consequence |
|---|---|---|
| analytic rank ≥ 4 | not proved for any particular curve, **per the Sage documentation** cited in the research | any counterexample at rank ≥ 4 is currently uncertifiable — treat as a strong prior, not a theorem we verified |
| `ellrank` output | `[r1, r2, s, L]`, exact only when `r1 == r2` | a gap is not a rank |
| numerical noise | `ellL1` on 389a1 returns −5.4e−29 at default precision | tolerance is a trap |
| LMFDB coverage | complete for `N < 500,000` | below that, BSD is already checked |
| rank conjecture | verified for `N < 130,000` (Miller 2011) | the interesting band starts at `N > 500,000` |
| laptop 2-descent | routine `N ≲ 10^6`; `10^8–10^9` is hours–days | the sweep is real work, not a script |

Combining these: the only certifiable search region is `N > 500,000` with
analytic rank `k ≤ 3`. Inside that region BSD is expected to hold, and no
counterexample has ever been observed anywhere. **There is no known positive
instance to calibrate against, so `p` is unknown and probably 0.**

That is why "find a counterexample" is `NO_ADMISSIBLE_ACTION`: not because the
checker is missing, but because the object has never existed and the only
region where it could be certified is the region where the conjecture has
already survived every test.

## 4. Load-bearing decomposition

1. **A working rank backend on Windows** (PARI 2.19, native installer;
   `ellrank`, `ellanalyticrank`). *Delete it and there is no `r` at all.*
2. **Analytic certification with an error bound**, not a tolerance. *Delete it
   and `k` is a numerical guess.*
3. **A conductor sweep above 500,000**, ordered by feasibility. *Delete it and
   there is no new curve.*
4. **The claim gate** (`bsd.py`). *Delete it and interval ranks get reported as
   ranks.*

## 5. Admissible actions

1. Install PARI/GP 2.19 (native Windows, `Pari64-2-19-0.exe`) and verify
   `ellrank` on `234446.a1` reproduces rank 4 with `r1 = r2`. **Still needed:
   there is no algebraic-rank backend without it.**
2. ~~Build an analytic-rank estimator.~~ **Done as a prototype.**
   `research/bsd-rank-smoke.md` has a PARI-free AFE estimator that returns
   0, 1, 2, 3 on `11a1`, `37a1`, `389a1`, `5077a1`, and matches LMFDB's
   `L(1)` to 29 digits. It is **numerical, not certified**: the returned order
   has no interval bound, so it cannot yet serve as the `k` in the terminal
   object.
3. Sweep `500,000 < N < 10,000,000` for curves with `k ≤ 3`, recording
   `(r1, r2, k)` and flagging `r1 == r2 != k`.
4. Report the count and the first mismatch, or the absence of one.

The gap between 2 and the terminal object is now specific: turn the AFE ladder
into a **certified** order of vanishing (interval arithmetic on the ladder, or
modular-symbol methods), or accept that the search can only produce
`COUNTEREXAMPLE_CANDIDATE`, never `COUNTEREXAMPLE`.

Step 4 is the only one that could produce the terminal object. Steps 1–2 are
prerequisites, not progress.

## 6. Kill criteria

- `r1 < r2` on more than 20% of the sweep → the 2-descent budget is too low;
  widen before drawing any conclusion.
- 10^5 curves with `r1 == r2` and no mismatch → stop. That is a coverage
  statement, not a counterexample.
- Any claimed mismatch where `k ≥ 4` → reject immediately; `k` is unprovable.

## What the research changed

Before: "scan LMFDB for a rank mismatch."
After: LMFDB is complete to 500,000 and curated to be consistent; a mismatch
there is impossible by construction. The only certifiable region is
`N > 500,000, k ≤ 3`, and in that region the expected number of counterexamples
is unknown and probably zero.

The honest output for this branch is **coverage, not a solution**. That is a
much smaller claim than the one this project started with, and it is the one the
evidence supports.
