# BSD: annotated corpus

The terminal object is `(E, r, k)` with `r` proved by descent, `k` proved by rigorous
zero-counting, and `r ≠ k`. The corpus is organised around one gate (Gross–Zagier–
Kolyvagin), one coverage boundary (conductor 500,000), and one calibration ladder (the rank
records).

## 1. The conjecture

**`doi:10.1515/crll.1965.218.79` (Birch–Swinnerton-Dyer 1965).** The conjecture in the form
the terminal object must refute:

```
ord_{s=1} L(E,s) = rank E(Q)
```

with the leading coefficient predicted as

```
L^{(r)}(E,1)/r! = (Ω · Reg · |Sha| · ∏ c_p) / |E(Q)_tors|²
```

The rank statement and the leading-term statement are separate: a rank mismatch refutes
the first, a Tamagawa/Sha mismatch refutes the second without changing `r` or `k`. The
repo's terminal object only tests the first.

## 2. The GZK gate

**`doi:10.1007/BF01388809` (Gross–Zagier 1986).** For analytic rank 1, the first derivative
is an explicit non-zero multiple of the canonical height of a Heegner point:

```
L'(E,1) = c · ĥ(P)
```

The rank-0 case is also treated. This is the formula that turns a vanishing order into a
rational point.

**`doi:10.1007/978-0-8176-4575-5_11` (Kolyvagin 1990).** Euler systems. The Heegner-point
Euler system bounds the Selmer group; together with Gross–Zagier it gives:

```
r_an(E/Q) ≤ 1  ⇒  r_an(E/Q) = r_MW(E/Q),  and Sha(E/Q) is finite.
```

**`arXiv:0707.0032` (Jetchev–Lauter–Stein 2007), Theorem 1.1 (p. 3)** is the precise modern
citation:

> Let `E/Q` be an elliptic curve with `r_an(E/Q) ≤ 1`. Then
> `r_an(E/Q) = r_MW(E/Q)`.

The paper also proves the contrapositive direction used for rank ≥ 2: Kolyvagin's
conjecture implies that `r_an ≥ 2` forces `Z_p-corank Sel_p ≥ 2` (abstract). That is the
statement that makes a rank-2 curve the first place a mismatch can hide.

**Consequence.** A counterexample must have `k ≥ 2`. Parity (`r ≡ k mod 2`, from the sign
of the functional equation) then restricts the shapes. Under parity alone the smallest are:

```
(r,k) = (0,2), (4,2), (6,2), ...     k = 2, r even, r ≠ 2
(r,k) = (1,3), (5,3), (7,3), ...     k = 3, r odd,  r ≠ 3
```

The 2026-11-02 review quotes `(4,2)` as the smallest shape. That is the smallest shape with
`r > k`; it is not the smallest under `r ≠ k`, because `(0,2)` and `(1,3)` are allowed by
parity. If the gate accepts any `r ≠ k`, `(0,2)` is the cheapest object to certify: `r = 0`
by an upper bound of 0 from 2-descent, `k = 2` by zero-counting. The corpus records both
readings rather than silently adopting one.

**What is ruled out.** `(2,0)`, `(3,1)`, `(4,0)`, `(5,1)` are not "hard"; they are
impossible by GZK. Any plan whose certifiable region is `k ≤ 1` is aimed at the one region
where the theorem already proves `r = k`.

## 3. Verification landscape

**`doi:10.1112/s1461157011000180` (Miller 2011).** Proves the BSD formula for **16,714 of
the 16,725** elliptic curves of conductor **< 5,000** with analytic rank 0 or 1 (abstract).
This is the largest published pointwise verification in the GZK range. Note the correction:
the repo's `plan/bsd.md` says "verified for `N < 130,000`"; the paper's abstract says
`N < 5,000` and 16,714/16,725. The 130,000 figure is not supported by this source.

**`lmfdb:reliability`.** The LMFDB states that for classical modular form L-functions the
displayed analytic rank has been rigorously verified in all cases, **including rank > 1**
(page text). For rank ≤ 1 the displayed value "must be correct". This is the source that
invalidates the claim "k ≥ 2 is not certifiable". The method is rigorous zero-counting /
interval arithmetic, not a pointwise numerical tolerance.

**`lmfdb:completeness`.** All curves of conductor **< 500,000** are in the database;
3,824,372 curves in 2,917,287 isogeny classes, maximum conductor 299,996,953. It also
contains all 7-smooth conductors and all prime conductors ≤ 300,000,000. Below 500,000 a
rank mismatch is impossible by construction (the data is curated to be consistent), so the
only interesting region is above it.

**`rubinstein:records`.** Calibration curves:

| rank | conductor | curve |
|---:|---:|---|
| 4 | 234,446 | `[1,-1,0,-79,289]` |
| 5 | 19,047,851 | `[0,0,1,-79,342]` |
| 6 | 5,187,563,742 | `[1,1,0,-2582,48720]` |

**`arXiv:math/0208056` (Elkies 2002).** For `E_D : D y² = x³ − x`, proves positive rank for
every `|D| < 10^6` with `|D| ≡ 5, 6, 7 mod 8` (abstract). The proof is computational:
Heegner points from CM points on modular curves, then a height check. This is the largest
explicit odd-analytic-rank family in the corpus.

**`arXiv:1006.1002` (Bhargava–Shankar 2015).** Ordered by height, the mean size of the
2-Selmer group is 3, so the limsup of the average rank is at most 1.5 (abstract). A random
curve is not where a rank-4 object will be found; the search must be targeted, not uniform.

**What is actually certifiable.** `r` is exact when a 2-descent returns `r1 = r2`;
`k` is exact when zero-counting (argument principle + interval arithmetic) certifies the
order. The AFE pointwise interval method can prove `k = 0` and `k = 1` (non-vanishing) but
cannot prove an exact zero of higher order; that is a limitation of the method, not of
mathematics. LMFDB's rigorous ranks are the existence proof that `k ≥ 2` is certifiable.

## 4. Analogues

**`goldfeld:numdam-1976` (Goldfeld 1976).** Reduces the effective class-number problem for
imaginary quadratic fields to the existence of an elliptic curve with a zero of order ≥ 3
at `s = 1` (Annali SNS Pisa 3 (1976), 623–663). Gross–Zagier supplied the rank-3 curve, and
the class-number problem followed. This is the cleanest example of a BSD-type analytic
input settling a classical Diophantine problem.

**`arXiv:1210.8231` (Tian 2014).** Congruent numbers via Heegner points: constructs points
on `m y² = x³ − x` for integers `m` with many prime divisors, extending Heegner's
prime-conductor construction to composite `m`. The congruent-number problem is the standard
conditional corollary of BSD: Tunnell's criterion decides it **if** BSD holds.

**`doi:10.2977/prims/31` (Barnet-Lamb–Geraghty–Harris–Taylor 2011).** Potential automorphy
for a family of Calabi–Yau varieties, yielding Sato–Tate for elliptic curves over Q. A
local-statistics conjecture about every elliptic curve, proved by automorphy; the model for
"a big conjecture with a complete method".

**`doi:10.2307/2118559` (Wiles 1995).** Every semistable elliptic curve over Q is modular,
hence Fermat's Last Theorem for all exponents ≥ 3. The canonical case study of a theorem
that is bigger than its original target and whose method (modularity lifting) became a
tool.

## 5. What the corpus implies

1. **`k ≤ 1` is a dead region.** GZK proves `r = k` there. The certifiable region and the
   interesting region do not overlap.
2. **The smallest admissible mismatch is `(0,2)` under `r ≠ k`, or `(4,2)` under `r > k`.**
   Either way it needs a certified analytic rank of 2 or more and a certified algebraic
   rank that differs.
3. **LMFDB is complete to 500,000.** A sweep below that cannot find anything new. Above it,
   the cost is 2-descent, not L-function evaluation: `r` is the expensive certificate.
4. **The analytic side is solved in principle.** Zero-counting plus interval arithmetic can
   certify `k`; the AFE ladder is a scout, not a certificate. The open engineering problem
   is a certified `k ≥ 2` implementation, and the open mathematical problem is an algebraic
   rank that disagrees.
5. **No positive instance has ever been observed.** The honest output of this branch is
   coverage, not a counterexample. `plan/bsd.md` already says this; the corpus supplies the
   numbers behind it.
