# Jacobian reconstruction literature: is there an inverse of the GGHV/GGV normal form?

Scope: the planar Jacobian normal form of Guccione–Guccione–Horruitiner–Valqui
(GGHV, arXiv:2204.14178) and the earlier Guccione–Guccione–Valqui (GGV) reductions
(arXiv:1401.1784, arXiv:1406.0886, arXiv:1605.09430, arXiv:1708.07936), plus the
standard monographs/surveys and the 2026 `(72,108)` artifacts.  All URLs checked
2026-11-02.

## Answer

No general inverse and no general image characterization is published.  The
GGHV/GGV reduction is one-directional: a polynomial counterexample is mapped into
a normal form in `L^(1) = K[x,x^-1,y]`, and every degree-bound paper uses only
that direction.  Partial reconstruction exists: GGHV22 §§5–6 reconstruct
polynomial pairs in the `(9,24)/(9,27)` and `(7,21)` cases, and GGV16 §8/Thm 8.12
returns to `L` in the `B=16` case, but neither covers `(8,28)`; GGV's separate
`C`-form (arXiv:1406.0886, Thm 1.9/Cor 1.17) is a genuine two-way
polynomial-representative theorem, but it is not the GGHV `L^(1)` normal form.
The standard practice is emptiness only: prove the normalized system empty
(wstrinz, Helali), never invert it.

## Objects and numbers

| item | value | source | checked |
|---|---|---|---|
| GGHV22 reduction | counterexample ⇒ `P,Q ∈ L^(1)`, `[P,Q]=x^2`, config (1) or (2); not an iff | arXiv:2204.14178, Prop. 4.3, p.10 | cited |
| GGHV22 last map | `φ(x)=x^-1`, `φ(y)=x^4 y` is an automorphism of `L^(1)` but not of `K[x,y]` | ibid. p.7; p.12 repeats "not an automorphism" | cited |
| GGHV22 §5 reconstruction | `D = Σ x^k D_k`, `P̃=D^2`, `Q̃=D^3+… ∈ K[x,y]`; `ψ(x)=x−D^2`; then `ψ(P̃),ψ(Q̃) ∈ K[x,y]` | ibid. §5, pp.17–18 | cited |
| GGHV22 §6 reconstruction | same `D`-trick; `P̃,Q̃ ∈ K[x,y]`, then `ψ` | ibid. §6, p.23 | cited |
| GGHV22 `(8,28)` | left open: "we couldn't solve the corresponding system of polynomial equations, thus it is left open" | ibid. p.1 | cited |
| GGHV22 polygons | every listed vertex lies in `N0^2`; with `N(P)` the convex hull of `Supp(P)` this forces the normalized `P,Q` into `K[x,y]` | ibid. p.10; GGV16 Def. of `H(P)`, p.4 | ran+cited |
| GGV16 method | "we remain all the time in `L` (or `L(l)`)" | arXiv:1401.1784, p.1 | cited |
| GGV16 Prop. 5.20 | standardizes an `(m,n)`-pair in `L(1)` by an automorphism of `L(1)`; "if `P,Q ∈ L`, then we can take `φ ∈ Aut(L)`" | ibid. p.29 | cited |
| GGV16 Thm 8.12 | after the `B=16` degree reduction, "we obtain `P6,Q6 ∈ L`" | ibid. §8, p.64 | cited |
| GGV14 Thm 1.9 | JC false iff `∃ P,Q∈K[x,y], C,F∈K[y]((x^-1))`, `C^n=P`, `Q=Σ λ_i C^{m−i}+F` | arXiv:1406.0886v3, Thm 1.9 | cited |
| GGV14 Cor 1.17 | counterexample with degrees `(m,n)` iff the finite system `S_t(n,m,(λ_i),y)` has a solution in `K[y]^{m+n−2}` | ibid. Cor 1.17 | cited |
| GGV14 polynomiality | `P=C^n∈K[x]` iff `(C^n)_{-k}=0` for all `k>0`; `Q∈K[x]` iff all its negative coefficients vanish | ibid. Rem. 1.13; arXiv:2506.05697 §1 | cited |
| van den Essen 2000 | Ch. 3 "Invertibility criteria and inversion formulae", pp.61–75: formula/algorithm for the inverse of a given Keller map | doi:10.1007/978-3-0348-8440-2 | cited |
| van den Essen 2021 | Ch. 4 "The Jacobian Conjecture: New Equivalences", pp.91–111; no normal-form reconstruction | doi:10.1007/978-3-030-60535-3 | cited |
| Moh 1983 | JC for `max(deg P,deg Q) ≤ 100` via approximate roots; four exceptional cases, discarded by hand-solving coefficient systems (detailed proof only for the smallest) | J. reine angew. Math. 340 (1983) 140–212; described in GGHV22 p.1 and GGV14 intro | cited |
| Abhyankar 1977 | Thm (18.13), p.131, and Cor. (18.14): each component of a Keller pair has at most two points at infinity; Thm (19.4), p.136: four equivalent formulations | TIFR Lectures 57, Part II, Ch. 6 | cited |
| Nowicki 1994 | derivation formulation (`k[X]^d`, rings of constants); no `L^(1)` normal form | doi:10.13140/2.1.3906.3844 | cited |
| Derksen 1994 | JC(n) equivalent to bounded inverse degree `BI(n)`; inverse degree, not inverse map | Comm. Algebra 22(12) (1994) 4793–4794 | cited |
| van den Essen Mathieu–Zhao | JC ⇒ Mathieu conjecture; Mathieu–Zhao spaces; no normal-form reconstruction | arXiv:1907.06107; doi:10.1007/978-3-030-42136-6_14 | cited |
| wstrinz Thm A | no pair `(P,Q) ∈ K[x,x^-1,y]^2` with `[P,Q]=x^2` and config (1) or (2) | github README:77–84; `d2/PROOF_72_108.md:38–42` | cited |
| wstrinz proof logic | `counterexample in config (1)/(2) ⇒ admissible germ ⇒ contradiction`; never the converse | `d2/PROOF_72_108.md:380–390` | cited |
| Helali 2026 | excludes both transcribed Laurent/Newton systems; conditional on GGHV22 Prop. 4.3 and transcription | doi:10.5281/zenodo.21479814 | cited |
| Ramírez–Valqui 2025 | Gröbner basis and complete solution set of the GGV system; restates the coefficient-wise polynomiality conditions | arXiv:2506.05697 | cited |
| Nguyen 2019 | complex JC through degree 104 by Newton-polygon technique; no inverse | arXiv:1902.05923 | cited |

## Source-by-source classification

Meaning of the columns: **(a)** constructs an inverse of the normalisation, or a
polynomial representative of a normalized pair; **(b)** characterizes the image
`A(K[x,y])` inside `L^(1)`, or states the polynomiality conditions in the normal
form; **(c)** neither, but uses the normal form one-directionally and proves
emptiness of the normalized set.

| source | (a) inverse? | (b) image? | (c) emptiness only? | verdict |
|---|---|---|---|---|
| GGHV22 (arXiv:2204.14178) | **partial / case-specific**: §§5–6, not `(8,28)` | no | yes for all cases except `(8,28)`, which is left open | necessary-only normal form; case-specific reconstruction devices |
| GGV16 (arXiv:1401.1784; J. Algebra 471 (2017) 13–74) | **case-specific**: §8/Thm 8.12 returns to `L` | no | yes | shape/degree constraints, one-directional |
| GGV14 (arXiv:1406.0886v3) | **yes, in the `C`-coordinate** (Thm 1.9/Cor 1.17) | **yes, in coefficient form** | no: it is a two-way parametrization of counterexamples | closest published "polynomial representative" theorem, but not the GGHV `L^(1)` normal form |
| GGV17 (arXiv:1708.07936) | no | no | yes: possible-corner enumeration | emptiness only |
| GGV16b (arXiv:1605.09430) | no | no | yes: Newton-polygon restrictions | emptiness only |
| GGV26 (Quaestiones Math. 2026, doi:10.2989/16073606.2026.2701437) | no | no | yes: lower-side restrictions, discards corners/families | emptiness only |
| van den Essen 2000, Ch. 3 | **yes for a given Keller map** (formal inverse formula/algorithm), not for the normalisation | no | neither for GGHV | different inverse problem |
| van den Essen 2021, Ch. 4 | no | no | neither | new equivalences, no reconstruction |
| Moh 1983 | no | no | yes: degree ≤100, coefficient systems | emptiness only |
| Abhyankar 1977 | no | no | yes: structural constraints (two points at infinity) | no Laurent normal form |
| Nowicki 1994 | no | no | neither: derivation formulation | different formulation |
| Derksen 1994 | no (inverse **degree** only) | no | neither | different inverse problem |
| van den Essen Mathieu–Zhao | no | no | neither | different reduction |
| wstrinz 2026 | no | no | **yes**: Theorem A over `K[x,x^-1,y]` | emptiness only, explicitly |
| Helali 2026 | no | no | **yes**: exact exclusion of the two Laurent systems | emptiness only |
| Ramírez–Valqui 2025 | **yes, in the `C`-coordinate** | **yes, coefficient form** | no | Gröbner analysis of the GGV system |
| Nguyen 2019 | no | no | yes: degree ≤104 | emptiness only |

### Notes on the decisive sources

**GGHV22.** Proposition 4.3 is explicitly one-directional: "*If there is a
counterexample … then there exist `P,Q ∈ L^(1)` with `[P,Q]=x^2` and one of the
following cases*" (p.10).  The last map of the chain is `φ(x)=x^-1`, `φ(y)=x^4y`;
the paper says at p.7 that this is "*an automorphism of `L^(1)=K[x,x^-1,y]` but
not of `K[x,y]`*", and at p.12 that it is "*not an automorphism*" (sc. of
`K[x,y]`).  No inverse is written down for the `(8,28)` case; the paper says it
"*couldn't solve the corresponding system of polynomial equations, thus it is
left open*" (p.1).  In §§5–6, for `(9,24)/(9,27)` and `(7,21)`, they do
reconstruct a polynomial pair in a special form (not necessarily the original
pair): "*However `C_k` in general is not a polynomial, so we will transform
the system into a system for certain `D_k`'s which are polynomials*" (p.17);
then `P̃:=D^2 ∈ K[x,y]`, `Q̃:=D^3+D^-1 λC_3^20+x^-4F_-4C_3^23+F̃ ∈ K[x,y]`
(p.17), and the automorphism `ψ(x)=x−D^2`, `ψ(y)=y` gives
`ψ(P̃),ψ(Q̃) ∈ K[x,y]` (p.18; §6 p.23 is the same pattern).  This is a
case-specific polynomial representative, not a general inverse of the
normalisation, and it does not cover `(8,28)`.

A useful correction to the repo's framing: every vertex in the two GGHV22
`(8,28)` configurations (and in the `(9,24)/(9,27)` configurations) lies in
`N0^2`.  Since `N(P)` is the convex hull of `Supp(P)` (GGV16 p.4), the normalized
`P,Q` are automatically in `K[x,y]`.  The obstruction is not polynomiality of the
normalized pair; it is membership in the proper subring `A(K[x,y]) ⊂ L^(1)`.
For the simplified T1–T5 composite recorded in `reconstruct.py`, `A` is an
automorphism of `L^(1)` (each stage is invertible on `L^(1)`), with
`A(x)=x^4y+s(x^-1)`, `A(y)=x^-1`; hence
`A(K[x,y]) = K[x^-1, x^4y+s(x^-1)]`.  What fails is that `A^-1` does not preserve
`K[x,y]`.  The full published chain has more stages, but the same phenomenon
occurs: its inverse does not preserve `K[x,y]`.  I found no paper that states a
membership criterion for `A(K[x,y]) ∩ K[x,y]`.

**GGV16.** The program deliberately stays in the Laurent ring: "*Opposed to most
of the authors working in this area … we remain all the time in `L` (or
`L(l)`)*" (p.1).  Proposition 5.20 standardizes an `(m,n)`-pair in `L(1)` using
an automorphism of `L(1)`, with the important caveat "*if `P,Q ∈ L`, then we can
take `φ ∈ Aut(L)`*" (p.29).  Section 8/Theorem 8.12 is the one place where the
degree-reduction chain is arranged to land back in `L`: "*performing on `(P,Q)` a
series of operations as described at the beginning of this section we obtain
`P6,Q6 ∈ L`*" (p.64).  That is a genuine published "return to `K[x,y]`", but it
is a case-specific `B=16` proof device, not an inverse of the GGHV normal form.

**GGV14 (the `C`-form).** Theorem 1.9 is a two-way criterion: JC is false iff
there exist `P,Q∈K[x,y]`, `C,F∈K[y]((x^-1))`, `n,m` with `n∤m`, `m∤n`, such that
`C=x+C_-1x^-1+…`, `C^n=P`, and `Q=Σ_{i=0}^{m+n-2} λ_i C^{m-i}+F`; "*furthermore,
under these conditions `(P,Q)` is a counterexample*".  Corollary 1.17 removes
`P,Q`: a counterexample with degrees `(m,n)` exists iff the finite system
`S_t(n,m,(λ_i),y)` has a solution in `K[y]^{m+n-2}`.  The associated pair is
reconstructed as `P=C^n`, `Q=…`; the image condition is exactly the vanishing of
the negative coefficients, `(C^n)_{-k}=0` and `(Σ λ_i C^{m-i}+F)_{-k}=0`
(Remark 1.13; restated in arXiv:2506.05697 §1).  This is a real published
"polynomial representative" theorem, but in the GGV `C`-coordinate, not the
GGHV `L^(1)` normal form.  It does not invert `A`.

**wstrinz 2026.** Theorem A is stated over the Laurent ring: "*There is no pair
`(P,Q)` in `K[x,x^-1,y]^2` with `[P,Q]=x^2` and*" either configuration
(`README.md:77–84`; `d2/PROOF_72_108.md:38–42`).  The proof is the two-step
`counterexample ⇒ admissible germ ⇒ contradiction` (`PROOF_72_108.md:380–390`);
it never constructs an inverse or a polynomial pair from a normal-form solution.
`README.md:233–239` explicitly says "Nothing here proves or refutes JC(2)" and
Corollary B only raises the degree bound.  This is the cleanest example of the
standard one-directional use.

**Helali 2026.** The Zenodo abstract says the computation "*excludes both
transcribed systems*" and that "*the conclusion … remains conditional on the
correctness and exhaustiveness of the published reduction and on its faithful
transcription*".  No reconstruction; emptiness of the two Laurent/Newton systems
is the whole claim.

**The enumeration/restriction papers (GGV17, GGV16b, GGV26, Moh, Nguyen).**
They compute possible corners, lower-side restrictions, or degree ranges for a
hypothetical counterexample.  None constructs an inverse or characterizes
`A(K[x,y])`; they all use the normal form (or a Newton-polygon restriction) as a
necessary condition and then exclude cases.  This is the standard degree-bound
method.

**Monographs and surveys.** Van den Essen (2000) Ch. 3 gives the formal inverse
formula and an algorithm that computes the inverse of a *given* polynomial map
(pp.61–75); that is a different problem from inverting the normalisation
`K[x,y] → L^(1)`.  The 2021 sequel Ch. 4 gives new equivalences of JC (pp.91–111)
and the Mathieu–Zhao chapter (pp.113–177) gives the Mathieu–Zhao reduction; neither
contains a normal-form reconstruction theorem.  Abhyankar (1977) works in
`A=k[x1,x2]`, not a Laurent ring; Theorems (18.13) and (19.4) are structural
constraints.  Moh (1983) proves the degree ≤100 bound by approximate roots and
coefficient systems.  Nowicki (1994), Derksen (1994), and van den Essen's
Mathieu–Zhao survey address different formulations (derivations/rings of constants;
inverse degree bounds; Mathieu–Zhao spaces) and do not address the GGHV normal
form.

### Standard practice

Across the degree-bound literature the normal form is used as a **necessary
condition**.  Moh, GGV17, GGHV22, Nguyen, wstrinz and Helali all assume a
polynomial counterexample, normalize it, and prove the normalized set empty or
exclude its Newton data.  None uses a normal-form solution to construct a
counterexample.  The only reconstructions I found are GGHV22 §§5–6 and GGV16
§8/Thm 8.12, and both are internal contradiction devices in specific cases, not a
general inverse and not the `(8,28)` case.  Therefore the degree-bound program is
**emptiness-only**; reconstruction is needed only for the opposite direction
(normal-form solution ⇒ counterexample), which no published degree-bound proof
attempts.

## Reproduced locally

`python problems/jacobian/reconstruct.py` (from the repo root):

```text
self-test: forward det = x^2; inverse(forward(x,y)) = (x,y);
           positive control (x+y^2, y) reconstructs with det = 1  [ok]
...
step 1: undo phi: (x,y) -> (x^-1, x^4 y)
  P = a0 + a1/x
  Q = b0 + b1*x**2*y
...
step 3: undo swap: (x,y) -> (y,x)
  P = a0 + a1/y
  Q = b0 - b1*c2 - b1*c3/y - b1*c4/y**2 + b1*x*y**2

exact inverse, before clearing:
  P0 = a0 + a1/y
  Q0 = b0 + b1*y**2*(-c2/y**2 - c3/y**3 - c4/y**4 + x)
  det J(P0,Q0) = a1*b1   (= 1 under a1*b1 = 1)

cleared pair (multiply P by y^1, Q by y^2):
  P' = a0*y + a1
  Q' = b0*y**2 - b1*c2*y**2 - b1*c3*y - b1*c4 + b1*x*y**4
  det J(P',Q') = -a0*b1*y**4
  constant? False
...
VERDICT: PARTIAL
```

This reproduces the `y^-1` obstruction exactly: `det J = 1` before clearing,
`det J = -a0 b1 y^4` after, so no polynomial clearing gives a nonzero constant.

`python research/scratch/jacobian-lit/checkA.py`:

```text
A(x)= x**2*(c2 + c3*x + c4*x**2 + x**2*y)
A(y)= 1/x
inv(A(x))= 0
inv(A(y))= 0
A(inv(x))= 0
A(inv(y))= 0
image of y under inv: -c2*y**2 - c3*y - c4 + x*y**4
det A: x**2
```

So `A` is an automorphism of `L^(1)` with `A^-1(u,v)=(v^-1, v^4(u-s(v)))`, and
`y ∉ A(K[x,y])`; the obstruction is the `v^-1` in the first component.

`python research/scratch/jacobian-lit/check_polygons.py`:

```text
GGHV22 (9,27) N(P): min x=0, min y=0; all in N0^2 = True
GGHV22 (9,27) N(Q): min x=0, min y=0; all in N0^2 = True
GGHV22 (9,24) N(P): min x=0, min y=0; all in N0^2 = True
GGHV22 (9,24) N(Q): min x=0, min y=0; all in N0^2 = True
GGHV22 (8,28) N(P): min x=0, min y=0; all in N0^2 = True
GGHV22 (8,28) N(Q): min x=0, min y=0; all in N0^2 = True
GGHV22 (8,28) N(P) config1: min x=0, min y=0; all in N0^2 = True
GGHV22 (8,28) N(Q) config1: min x=0, min y=0; all in N0^2 = True
```

This is the check behind the statement that the normalized pairs are already
polynomial; the reconstruction gap is membership in `A(K[x,y])`, not the ambient
Laurent ring.

## What I could not find

- No published **general** inverse of the GGHV/GGV `L^(1)` normalisation, and no
  published characterization of `A(K[x,y]) ⊂ L^(1)` or of `A(K[x,y]) ∩ K[x,y]`
  for the composite used in GGHV22 Prop. 4.3.
- No published "polynomial representative" or "reconstruction" theorem for the
  `(8,28)` normal form.  GGHV22 §§5–6 cover `(9,24)/(9,27)` and `(7,21)` only;
  GGV16 §8/Thm 8.12 covers `B=16` only.
- No paper that uses the GGHV normal form as a **sufficiency** or as a
  parametrization of counterexamples.  GGHV22 Prop. 4.3 is one-directional and
  the `(8,28)` case is explicitly left open (p.1).
- No peer-reviewed version of GGHV22, wstrinz, or Helali as of 2026-11-02.
  GGHV22 is arXiv:2204.14178v1 with no journal reference; wstrinz and Helali are
  public repository/Zenodo artifacts.
- No full text of Moh 1983 (De Gruyter paywall); I relied on the descriptions in
  GGHV22 p.1 and GGV14 intro.  No full text of Nowicki 1994 (DOI metadata and
  abstract only).
- No source explicitly says "the standard practice is emptiness only"; that
  statement is my inference from the one-directional logic in GGHV22, GGV16,
  wstrinz, and Helali, and from the absence of any sufficiency use.
