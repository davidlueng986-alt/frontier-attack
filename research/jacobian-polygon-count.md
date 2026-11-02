# Admissible Newton polygons at total degree 105–110

## Answer

Under the full published constraints the band is essentially empty: GGHV 2022
(arXiv:2204.14178, Theorem 2.1) proves that a planar counterexample with
`max(deg P, deg Q) < 125` must have degree pair `(72,108)` or `(108,72)`, and
Proposition 4.3 reduces the `(72,108)` case to two normalized Newton-polygon
configurations (one if the displayed configuration 1 is rejected by the
Cassou-Noguès no-bisector condition). The resulting degree → count table for
105…110 is `0, 0, 0, 2, 0, 0` — the largest value is 2, so nothing exceeds
`10^9`. The five bullets alone are **not** a decision procedure: an exact
enumeration of a conservative subfamily gives 3,691,505 polygons at 105,
5,610,784 at 108 and 41,559 at 110, while dropping the GGHV normal-form
condition `m,n > 1` gives `3.36e9` at 106, `4.81e9` at 108 and `6.84e9` at 110,
and allowing equal component degrees gives `~10^17`. The `10^9` gate therefore
turns on the published normal-form/complete-chain theorems, not on the five
bullets being numerically small. Gate verdict: `ADMISSIBLE_ACTIONS_EXIST`, but
the only action is the `(72,108)` case, not a broad degree-105–110 sweep.

## Objects and numbers

| item | value | source | checked |
|---|---|---|---|
| degrees below 125 | only `(72,108)` / `(108,72)` | GGHV, arXiv:2204.14178, Thm. 2.1 | cited (PDF text read) |
| `(72,108)` Newton configurations | 2 in Prop. 4.3: config 1 has slope-1 edge, config 2 does not | GGHV, arXiv:2204.14178, Prop. 4.3 | cited; slope check ran |
| full published count, 105…110 | `0, 0, 0, 2, 0, 0` | Thm. 2.1 + Prop. 4.3 | cited |
| count with no-bisector imposed on the displayed configs | `0, 0, 0, 1, 0, 0` | this note, `check_config` | ran |
| strict normal-form arc, `D=105` | 3,691,505 polygons in the S3 subfamily | this note | ran |
| strict normal-form arc, `D=108` | 5,610,784 polygons in the S3 subfamily | this note | ran |
| strict normal-form arc, `D=110` | 41,559 polygons in the S3 subfamily | this note | ran |
| S3 at the GGHV corner `(a,b)=(8,28)`, `e=36` | 994,614 polygons | this note | ran |
| weak arc (`m=1` allowed), `D=106/108/110` | 3,356,627,123 / 4,807,171,343 / 6,839,140,446 | this note | ran |
| equal-degree arc, `D=105…110` | 5.57e16 … 2.36e17 (log10 16.7458 … 17.3731) | this note | ran |
| growth, inclusive arc | `1.3371, 1.3360, 1.3349, 1.3338, 1.3327` per degree | this note | ran |
| small-degree brute force, `d=3` | 0 polygons (all degree-3 Keller pairs triangular) | Starkov 2016, doi:10.15393/j3.art.2016.3510 | cited; brute force ran |
| small-degree brute force, `d=4` | 1 shape `conv{(0,0),(1,0),(1,3),(0,1)}`, but 0 strict minimal-counterexample data (`m=1`, gcd 4 < 16) | Starkov 2017, doi:10.15393/j3.art.2017.3870 | cited; brute force ran |
| GGHV 2017 complete-chain cases, `max degree <= 150` | 34 | arXiv:1708.07936, §6 | cited |
| GGHV chain enumerator source | C++ implementation described, no code URL | arXiv:1708.07936, §4 | cited |

### What was enumerated

Published inputs: (C1) Abhyankar's two-points-at-infinity restriction (TIFR
lectures, `mathweb.tifr.res.in/Documents/Publications/Lectures/tifr57.pdf`);
(C2) similarity of `N0(f)` and `N0(g)` with the origin as centre (Lang 1991,
doi:10.1016/0022-4049(91)90128-O; restated in Lee–Li 2024, Prop. 2.1);
(C3) the trapezoid `T(a,b)` and the Cassou-Noguès no-bisector condition
(Makar-Limanov, arXiv:2106.06869); (C4) `gcd(deg f,deg g) >= 16` (Heitmann
1990, doi:10.1016/0022-4049(90)90005-3); (C5) the degree wall `<= 104`
(Nguyen 2025, arXiv:1902.05923; Moh 1983, doi:10.1515/crll.1983.340.140).

The script `research/scratch/polygons.py` implements two layers.

**Layer 1 — degree arithmetic (exact).** For a total degree `D`, write
`e = a+b`, `D = n e`, `deg f = m e`, `deg g = n e`, with `1 <= a < b`,
`e >= 16`, `m,n` coprime and `2 <= m < n`. This is the GGHV minimal-
counterexample normal form (`m,n > 1` is part of their standard `(m,n)`-pair
definition). It gives the degree-data counts

| D | admissible `(e,a,b,m,n)` tuples |
|---|---|
| 105 | 47 |
| 106 | 0 |
| 107 | 0 |
| 108 | 38 |
| 109 | 0 |
| 110 | 30 |

The zeros are arithmetic, not a failure of the generator: 106 = 2·53 has only
`e=53, n=2`, so no `m >= 2`; 107 and 109 are prime.

**Layer 2 — primitive polygon S3 (exact for the subfamily, a lower bound for
all admissible polygons).** A primitive lattice polygon `P0` is built inside

```
T(a,b) = conv{(0,0), (a,0), (a,b), (0,b-a)}
```

with:

- `(0,0)` and `(a,b)` as vertices; `1 <= a < b`;
- bottom edge `(0,0)–(p,0)` with `1 <= p <= a`;
- a convex lower-right chain `(p,0) -> (a,k)`, `0 <= k < c`, all edge slopes
  positive, strictly increasing and `!= 1`;
- vertical right edge `(a,k)–(a,b)`;
- a convex upper-left chain `(a,b) -> (0,c)`, `1 <= c <= b-a-1`, all edge
  slopes strictly decreasing and `> 1`;
- left edge `(0,c)–(0,0)`.

The two chains are separated vertically (`k < c`), so the polygon is convex and
inside `T(a,b)`. It contains `(1,0)` and `(0,1)`, has the unique top vertex
`(a,b)` (hence two points at infinity), and has no edge of slope 1. The number
of such polygons is counted exactly by a 2-D unbounded knapsack over primitive
edge directions, with multiplicities representing one edge vector `t·(p,q)`.
`S1` is the `k=0` subfamily (maximal vertical right edge); `S2` adds exactly
one lower-right edge; `S3` is the main model with arbitrary lower-right chains.
`S3(e)` denotes the sum over all `a < b`, `a+b=e`.

**What is exact.**

- Lattice polygons, all arithmetic in Python integers (no floating point).
- Trapezoid containment, origin/top vertices, `(1,0)` and `(0,1)` membership.
- Unique top vertex and therefore the two-points-at-infinity condition for
  every polygon in the family.
- Similarity with the origin as centre: `N0(f)=m P0`, `N0(g)=n P0`.
- `gcd(deg f, deg g) = e >= 16`.
- No edge parallel to the bisector `y=x` (slope 1).
- Total degree `D = n e`.

**What is approximated or omitted.**

- `S3` is a **proper subfamily**: polygons whose lower chain rises above
  `c` but still stays below the upper chain are not counted. Thus all S3
  numbers are lower bounds on the five-bullet count.
- The GGHV complete-chain machinery (Algorithms 1–8 of arXiv:1708.07936:
  valid edges, children, `PLLC`, divisibility conditions, last lower corners)
  is **not** implemented. Those conditions are strong: at the GGHV corner
  `(a,b)=(8,28)` the S3 subfamily alone has 994,614 polygons, while the full
  chain conditions leave 2 (Prop. 4.3).
- The Makar-Limanov side conditions `deg_y(g) > deg_y(f)` and
  `deg_y(f) ∤ deg_y(g)` are not implemented.
- Heitmann's `gcd != 2p` condition is not implemented; only `gcd >= 16`.
- The GGHV 2022 `< 125` theorem is used only in the full-literature table, not
  in the five-bullet enumeration.
- "Exactly two points at infinity" is enforced by the unique top vertex; that
  is exact inside `T(a,b)`, but it is a sufficient, not a necessary, lattice
  condition in general.

### Table 1 — full published constraints (GGHV 2022 is unrefereed)

Applying GGHV 2022 Theorem 2.1 and Proposition 4.3:

| total degree D | admissible Newton-polygon configurations |
|---|---|
| 105 | 0 |
| 106 | 0 |
| 107 | 0 |
| 108 | 2 |
| 109 | 0 |
| 110 | 0 |

This is the table used by the JSON summary: `count_d105 = 0`,
`count_d110 = 0`, `exceeds_1e9 = false`. For context, the earlier GGHV
2017 complete-chain enumeration lists 34 possible cases with
`max(deg P,deg Q) <= 150`; the only one in the 105–110 band is the 108 case.
The two configurations (GGHV 2022, Prop. 4.3) are:

```
config 1:
  N(P) = {(0,0),(1,0),(8,14),(8,16),(0,8)}
  N(Q) = {(0,0),(2,1),(12,21),(12,24),(0,12)}
config 2:
  N(P) = {(0,0),(1,0),(8,14),(8,16)}
  N(Q) = {(0,0),(2,1),(12,21),(12,24)}
```

Config 1 contains the edge `(8,16)–(0,8)` of slope 1 (and the corresponding
`(12,24)–(0,12)` in `N(Q)`), so applying the Cassou-Noguès no-bisector
condition literally to the displayed polygons leaves only config 2, i.e. count
1 at degree 108. The paper's polygons live in a localized ring `L(1)` after
non-polynomial transformations, so the no-bisector theorem may not transfer
directly; the table reports both readings.

### Table 2 — five bullets only, three normal-form choices

`strict` is the GGHV normal-form arc (`m,n > 1`, `m < n`). `weak` allows
`m = 1` but keeps `m < n`. `inclusive` also allows equal degrees (`m=n=1`). All
entries are exact counts of the S3 subfamily, hence lower bounds on the full
five-bullet count. Applying the further Heitmann/GGHV condition `gcd != 2p`
removes the `D=110` entry, because `e=22=2·11`; the 105 and 108 entries are
unchanged (21=3·7, 35=5·7, 18=2·9, 27=3³, 36=2²·3²). Table 2 is a sensitivity
table, not a competing count: its 105 and 110 entries are killed by the
GGHV 2022 theorem in Table 1.

| D | strict (`m,n>1`) | weak (`m>=1`) | inclusive (equal degrees allowed) | log10(inclusive) |
|---|---|---|---|---|
| 105 | 3,691,505 | 7,365,510 | 55,692,298,242,837,849 | 16.7458 |
| 106 | 0 | 3,356,627,123 | 74,468,950,988,958,239 | 16.8720 |
| 107 | 0 | 0 | 99,492,100,865,201,547 | 16.9978 |
| 108 | 5,610,784 | 4,807,171,343 | 132,812,706,388,624,149 | 17.1232 |
| 109 | 0 | 0 | 177,146,457,564,382,321 | 17.2483 |
| 110 | 41,559 | 6,839,140,446 | 236,086,600,450,098,024 | 17.3731 |

The strict column is the main five-bullet result. It counts distinct
`N0(f)=m P0` supports; `n=D/e` is fixed by `D` and `e`, so this is also the
number of distinct `(N0(f),N0(g))` support pairs for that `D`. In the weak and
inclusive columns the `a=1` primitive shapes are removed when `m=1`, because
they have top `x`-coordinate 1 and violate `n>m>1`; the correction is `O(e)` and
does not change any exponent. The weak and inclusive columns show that the
`10^9` answer is model-dependent: dropping `m,n > 1`
puts 106, 108 and 110 above `10^9`; allowing equal degrees puts every degree
at `~10^16.7–10^17.4`. The published standard `(m,n)`-pair normal form is what
excludes both.

### Table 3 — strict breakdown by primitive degree

| D | e | S3(e) | allowed `m` | contribution |
|---|---|---|---|---|
| 105 | 21 | 8,725 | 2,3,4 | 26,175 |
| 105 | 35 | 3,665,330 | 2 | 3,665,330 |
| 108 | 18 | 2,100 | 5 | 2,100 |
| 108 | 27 | 129,087 | 3 | 129,087 |
| 108 | 36 | 5,479,597 | 2 | 5,479,597 |
| 110 | 22 | 13,853 | 2,3,4 | 41,559 (0 after `gcd != 2p`) |

### Growth ratios

For the strict column the consecutive ratios are `0, undefined, undefined, 0,
undefined`; the zeros are the divisor arithmetic above, not a decrease in the
underlying polygon count. The inclusive column, dominated by `S3(D)`, gives the
meaningful growth rate:

| transition | inclusive ratio |
|---|---|
| 105 → 106 | 1.3371 |
| 106 → 107 | 1.3360 |
| 107 → 108 | 1.3349 |
| 108 → 109 | 1.3338 |
| 109 → 110 | 1.3327 |

So the raw polygon count grows by about `1.33×` per degree (≈ `4.2×` over the
five-degree band).

### Does anything exceed `10^9`?

- **Strict five-bullet S3 subfamily:** no; the maximum is 5,610,784 at 108.
  This is a lower bound on the full five-bullet count, so it does not by
  itself prove that the full count is below `10^9`.
- **Weak model (`m=1` allowed):** yes. 106: `3.36e9` (exponent 9.53);
  108: `4.81e9` (9.68); 110: `6.84e9` (9.84).
- **Inclusive model (equal degrees):** yes, by a wide margin:
  `5.57e16` at 105 (exponent 16.75) to `2.36e17` at 110 (exponent 17.37).
- **Full published constraints:** no; the count is 2 at 108 and 0 elsewhere.

### Small-degree sanity check

The script includes an independent exhaustive check for small `(a,b)`: it
enumerates all subsets of lattice points of `T(a,b)`, takes convex hulls,
dedupes, and filters the same shape conditions. It does not call the S3
generator.

| e = a+b | exhaustive brute force | S1 | S2 | S3 |
|---|---|---|---|---|
| 3 | 0 | 0 | 0 | 0 |
| 4 | 1 | 1 | 1 | 1 |
| 5 | 2 | 2 | 3 | 2 |
| 6 | 7 | 4 | 7 | 4 |
| 7 | 14 | 6 | 13 | 6 |

- **Degree 3:** 0 shapes. This matches Starkov's complete classification
  (all degree-3 Keller pairs are triangular). The only pair with
  `n>m>0` and `m+n=3` is `(1,2)`, whose `m=1` fails `n>m>1`; the trapezoid
  `T(1,2)` also has a slope-1 top edge, so even the relaxed shape is empty.
- **Degree 4:** the exhaustive shape count is 1, the quadrilateral
  `conv{(0,0),(1,0),(1,3),(0,1)}`. It has top vertex `(1,3)`, so `m=1`, and
  the strict minimal-counterexample model gives 0 at degree 4; the known
  quartic non-triangular Keller families `Φ_{a,b,d}` and `F_{c,s}` (Starkov
  2017) have `deg f = 2`, `deg g = 4`, `gcd = 2 < 16`; `F_{c,s}` has
  `N(f)` with top vertex `(0,2)` on the `y`-axis, and `Φ_{a,b,d}` has the
  non-monomial top form `(ax+by)^2`, so neither is the minimal-counterexample
  trapezoid shape. The generator therefore reproduces the empty
  **counterexample** set at degree 4, but not the full list of non-triangular
  tame Keller maps, which it was never meant to enumerate.
- **Degrees 5–7:** exhaustive totals 2, 7, 14 versus S3 2, 4, 6. The S3
  generator is a proper subfamily as designed; it never overcounts.

## Reproduced locally

Exact command (from the repository root):

```text
python research/scratch/polygons.py --json research/scratch/polygons.json
```

Output (key lines):

```text
== small-degree brute force vs S3 ==
e  brute_total  S1  S2  S3
3 0 0 0 0
4 1 1 1 1
5 2 2 3 2
6 7 4 7 4
7 14 6 13 6

== strict GGHV model: total degree -> count ==
D=105  data=    47  strict=        3691505  weak=        7365510  inclusive=   55692298242837849
D=106  data=     0  strict=              0  weak=     3356627123  inclusive=   74468950988958239
D=107  data=     0  strict=              0  weak=              0  inclusive=   99492100865201547
D=108  data=    38  strict=        5610784  weak=     4807171343  inclusive=  132812706388624149
D=109  data=     0  strict=              0  weak=              0  inclusive=  177146457564382321
D=110  data=    30  strict=          41559  weak=     6839140446  inclusive=  236086600450098024

== growth ==
strict ratios [0.0, None, None, 0.0, None]
weak ratios [455.7222952653652, 0.0, None, 0.0, None]
inclusive ratios [1.3371498993316389, 1.336021248371735, 1.3349070452192737, 1.3338065489459487, 1.3327198505467965]

== GGHV 2022 (72,108) configurations ==
config_1 {'N(P)': {'top': (8, 16), 'n_gt_m_gt_1': True, 'has_slope_1_edge': True, ...}, ...}
config_2 {'N(P)': {'top': (8, 16), 'n_gt_m_gt_1': True, 'has_slope_1_edge': False, ...}, ...}
```

The JSON payload with the full S1/S2/S3 tables for every degree 2…110 is
`research/scratch/polygons.json`.

Two independent recursive chain enumerators (not the knapsack) were also run:
`verify_chains.py` agrees with `U[a][b-c]` for `(a,b) <= (5,8)`, and
`verify_lower.py` agrees with `Lcum[a][k]` for `(a,k) <= (6,4)`. Both print
`ok` on every tested case.

## What I could not find

- No released source for the GGHV C++/PostgreSQL enumerator (arXiv:1708.07936
  says only "website under development"), so their 34-case table for
  `max(deg P,deg Q) <= 150` cannot be regenerated; the two `(72,108)` entries
  are taken from arXiv:2204.14178.
- No machine-readable definition of "admissible Newton polygon" beyond the
  GGHV chain algorithms. The five bullets are necessary conditions, not a
  decision predicate; this is why the report separates the five-bullet count
  from the full-literature count.
- No exact count of **all** lattice polygons in `T(a,b)` for `e=36` (or any
  `e` in the band). The S3 numbers are lower bounds, and the brute-force
  comparison is only feasible up to `e=7`.
- No independent confirmation of the GGHV 2022 `< 125` theorem; it is an
  unrefereed preprint (as noted in `research/jacobian-enumeration.md`).
- No resolution of the config-1 slope-1 discrepancy: GGHV's configuration 1
  has a slope-1 edge in the displayed coordinates, but the polygons are
  normalized in a localized ring, so it is not clear whether the
  Cassou-Noguès condition applies to them directly. The table gives both the
  published count (2) and the literal no-bisector count (1).
- A 2026 repository, `wstrinz/plane-jacobian-72-108`, claims a conditional
  structural exclusion of both `(72,108)` configurations. It is not
  peer-reviewed and is not used in the count; if accepted, the full count in
  105–110 is 0 at every degree. The present note counts what the published
  constraints admit, not what has been subsequently excluded.
- No proof that the full five-bullet count (as opposed to the S3 subfamily)
  stays below `10^9`; the weak model already exceeds it at 106, 108 and 110.
  The only unconditional `10^9` statement is the full-literature one: 2 at
  108.
