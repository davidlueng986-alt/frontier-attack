# Degree pairs at `max(deg P, deg Q) >= 125`: transcription and triage

## Answer

**24 pairs, not 27.** `README.md:183` heads the section "27 cases still open", but its own
table (`README.md:190-191`) gives 10 settled + **24 open**; `d2/corner_atlas.json` has
`"open_frontier": 24` and exactly 24 rows at `max_deg >= 125`. The heading is stale.
**The list is split across two artefacts of different shape:** `d2/NEXT_CASES.md:149-169`
is the only prose table but has 19 lines compressed by multiplicity (`x2 chains`,
`x4 chains`); `corner_atlas.json` expands to one row per case. They agree exactly — 24.

**Only 9 of the 24 have a Newton-polygon configuration written down anywhere** in the repo
(the GGV5/GGHV22 `.tex` sources are not redistributed, `d2/paper_src/README.md`, so
"published" means the repo's line-cited transcription). Seven full, two partial; the
other 15 are **unknown**. I have not filled a gap by scaling a neighbouring polygon.

## Objects and numbers

`L` = chain length; `cfg` = Newton configuration; `src` = `d2/NEXT_CASES.md` line. Every
row is also an id in `d2/corner_atlas.json`.

```text
 max  (degP,degQ)   A0              (m,n)  L  cfg      src
 125  (75,125)      (5,20)          (3,5)  1  KNOWN    151
 126  (84,126)      (7,35)          (2,3)  1  unknown  152
 126  (84,126)      (12,30)         (3,2)  2  unknown  153
 128  (96,128)      (8,24)          (3,4)  2  unknown  154
 132  (88,132)      (11,33)         (2,3)  1  unknown  155
 135  (90,135)      (9,36)          (2,3)  1  KNOWN    156
 135  (90,135)      (9,36)          (3,2)  1  KNOWN    157
 135  (90,135)      (12,33)         (2,3)  1  partial  158
 135  (90,135)      (9,36)->(9,24)  (2,3)  2  unknown  159
 140  (56,140)      (7,21)          (2,5)  1  KNOWN    160
 140  (84,140)      (7,21)          (3,5)  1  KNOWN    161
 144  (96,144)      (12,36)         (2,3)  2  unknown  162 x4
 144  (96,144)      (8,40)          (3,2)  2  unknown  163
 144  (96,144)      (12,36)         (3,2)  3  unknown  164
 144  (108,144)     (8,28)          (3,4)  1  partial  165
 147  (42,147)      (6,15)          (2,7)  1  unknown  166
 147  (63,147)      (6,15)          (3,7)  1  unknown  167
 147  (98,147)      (7,42)          (2,3)  1  KNOWN    168
 147  (98,147)      (7,42)          (3,2)  1  KNOWN    168
 150  (100,150)     (10,40)         (3,2)  2  vtx only 169 x2
```

Rows 12-15 are four distinct cases (intermediates `(12,33)`, `(9,24)`, `(21/4,9)`,
`(21/4,9)`); rows 23-24 likewise. Rows 23-24 have one vertex (`j*=10`,
`SECOND_CORNER.md:342`) but no polygon, so they are **not** counted as config-known.

**Why the nine are known.** `(5,20)`: `PASSPORT_75_125_REPAIR.md:27-33`. `(7,21)`:
published in GGHV22, `SECOND_CORNER.md:19,45-53`. `(7,42)`, `(9,36)`: explicit `N(P)`,
`N(Q)` at `(2,3)`, "(3,2): the same two polygons, swapped" (`SECOND_CORNER.md:181-191`).
All follow the proved closed form `SECOND_CORNER.md:132-141`:
`Delta'={(0,0),(l-1,0),(l,1),(0,a0)}`, `kappa=l-2`, `N(P)=m*Delta'`, `N(Q)=n*Delta'`.
`(12,33)`: `Delta'={(0,0),(2,0),(3,11),(3,12)}`, `kappa=1`, but `N(C)` undefined
(`SECOND_CORNER.md:354-360`) — partial. `(8,28)/(3,4)`: reduced base edge `x^4 C4(y)`,
`C4=y^3(y+1)h4(y)`, `deg h4=4` (`CORNER_144_COMPARISON.md:53-63`), but no vertex set.

**Why the 15 are unknown** (repo says so, not for lack of search): `(6,15)` reduction
"performed in no paper", `Delta'` still a work item (`CORNER_ATLAS.md:196`, `PHI_F7.md:98`);
`(8,24)`, `(12,30)`, `(10,40)` have `A_0'=(2,0)` with one vertex undetermined
(`SECOND_CORNER.md:346-352,386-388`); `(12,36)`, `(8,40)` and the second `(9,36)` row are
length-2 chains with integer `A_1`, "out of scope by §5.1" (`SECOND_CORNER.md:348-352`);
`(7,35)` and `(11,33)` appear nowhere.

## Ranking

Cost proxy `(deg_red/25)^2 x (max_deg/200)`; config completeness is a gate. The proxy
ignores cascade depth `(t,kappa)` and tooling reuse, which decide ranks 2-3, so #3's lower
number is not a lower rank.

| rank | case | max | deg_red | cost |
|---|---|---|---|---|
| 1 | (75,125), F_2(3,5) | 125 | (15,25) | 1.00 |
| 2 | (90,135), (9,36)/(2,3) | 135 | (18,27) | 1.31 |
| 3 | (98,147), (7,42) | 147 | (14,21) | 0.86 |

**1. `(75,125)`, `F_2(3,5)`, `A_0=(5,20)`.** Smallest `max_deg` in the list; reduced degree
25, the second smallest among config-known rows; the only row with infrastructure: `C_SERIES_75_125.md`,
`G_SYSTEM_75_125.md`, `PHI_75_125.md`, `weight_lemma_75_125.py`, plus a proved obstruction
(`PROOF_72_108.md`, Thm 11.4). Its `j=0` sibling `(50,75)` is killed by GGV3 §5 at the same
corner, so a kill template exists one rung down. Unique row at the 125 bound: it alone
moves the frontier.

**2. `(90,135)`, `(9,36)/(2,3)`.** The only other corner with `N(P)`, `N(Q)` printed in
both orientations (`SECOND_CORNER.md:181-185`), and `(t,kappa)=(4,2)` matches `(72,108)`,
so cascade/`Phi` tooling transfers without re-deriving the chart. Four rows share the pair,
amortising setup; the fourth (row 9) has no polygon and would be dropped.

**3. `(98,147)`, `(7,42)`.** Reduced degrees `(14,21)`, the smallest coefficient system
here, both orientations printed. Penalty: `max=147` and `(t,kappa)=(6,4)` give the deepest
window cascade of the three.

Runners-up: `(56,140)` has the list's smallest total (196) and a *published* `(7,21)`
polygon, but `deg_red=35`; `(84,126)` has `max=210` and `deg_red=21`, but its `(7,35)`
polygon is class-inferred, not written.

### Top-1 corner data, verbatim

`d2/PASSPORT_75_125_REPAIR.md:29-33`:

```text
  Delta' = {(0,0),(3,0),(4,1),(0,5)}
  N(P) = 3*Delta' = {(0,0),(9,0),(12,3),(0,15)}      N(Q) = 5*Delta' = {(0,0),(15,0),(20,5),(0,25)}
  [P,Q] = x^2      reduced degrees (15,25)
```

`d2/case_dossier_F2_j1_75_125.json` (`case`): `A0=[5,20]`, `A0'=[1,0]`, `A1="(7\5,2)"`,
`a=3`, `b=5`, `m=3`, `n=5`, `k=1`, `v11=25`, diophantine
`"(m+n)*q*k - n*(q*l - p) = 1 = k  OK"`; `corner_signature`: `t=4`, `kappa=2`, `q=1`,
`dg=0`, `e=3`. The closed `(72,108)` pair for contrast: `N(P)=conv{(0,0),(1,0),(8,14),(8,16)}`,
`N(Q)=conv{(0,0),(2,1),(12,21),(12,24)}` (`README.md:79`).

## What I could not find

- Any Newton polygon at `(6,15)`, `(7,35)`, `(11,33)`, `(8,24)`, `(12,30)`, `(12,36)`,
  `(8,40)`: no vertex set, no `Delta'`, no published reduction.
- The GGV5/GGHV22 `.tex` files. `paper_src/` holds `next_cases.py`,
  `upstream_facts.json`, `upstream_quotes.json` only; every "published" claim above is
  second-hand from the repo's line-cited transcriptions.
- Any referee report for GGV5's 34-case enumeration, or for the (72,108) exclusion that
  makes 125 the bound.
- A justification of the "27" in `README.md:183`; nothing in the tree supports it.
