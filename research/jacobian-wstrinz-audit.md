# Jacobian (72,108): audit of `wstrinz/plane-jacobian-72-108`

## Answer

**Claim, not a conjecture and not a pure computation.** `d2/PROOF_72_108.md:38-42` (mirrored at `README.md:77-84`) states Theorem A: for every characteristic-0 field `K` there is no pair `(P,Q) in K[x,x^-1,y]^2` with `[P,Q]=x^2` and either GGHV22 Prop. 4.3 (8,28) Newton configuration — config (2) `N(P)=conv{(0,0),(1,0),(8,14),(8,16)}`, `N(Q)=conv{(0,0),(2,1),(12,21),(12,24)}`, or config (1) with the extra corners `(0,8)`, `(0,12)` — and Corollary B: conditional on Prop. 4.3 being exhaustive and on the alpha-strip WLOG `[QQ1]` (Prop. 2.1), every plane JC counterexample has `max(deg P,deg Q) >= 125`. It is a prose proof with exact sympy sub-lemmas; it uses the GGHV22 §5 normal form (`ell_{1,0}(P)=R^2`, `R=x^4 y^7(y+1)`, `P=C^2`; `PROOF_72_108.md:121-141`, `d2/paper_src/upstream_facts.json:33-35`) but **not** Cassou-Noguès/no-bisector (0 occurrences in the repo), so config (1) is treated directly. The Lean artifacts check an `f31` ideal identity and a corner-144 bridge (`lean/Cert/Certificate.lean`, `lean/Cert/BridgeOrder.lean`), not Theorem A; `lean/Challenge.lean:185,194` still contains the two challenge `sorry`s. **Verdict: No.** Even granting Theorem A, it excludes one degree pair and raises the plane lower bound to 125; it does not kill the planar Jacobian search — the repo says so itself (`README.md:233-236`), and 24 degree pairs with `max_deg >= 125` remain open in GGV5's `<=150` list (`README.md:185-191`). It does kill this project's `(72,108)` terminal object, if accepted.

## Objects and numbers

| item | value | source | checked |
|---|---|---|---|
| exact claim | no `(P,Q)` in `K[x,x^-1,y]^2`, `[P,Q]=x^2`, either (8,28) configuration | `d2/PROOF_72_108.md:38-40`; `README.md:77-82` | cited |
| conclusion | `max(deg P,deg Q) >= 125`, conditional | `d2/PROOF_72_108.md:42` | cited |
| assumption 1 | GGHV22 Prop. 4.3 exhaustive; audit 20/20; `L^(1)=K[x,x^-1,y]` | `PROOF_72_108.md:74-79`; `prop43_audit.py`; `HELALI_ADJUDICATION.md:72` | ran/cited |
| assumption 2 | `[QQ1]`: `Q=C^3+lambda C^-1+F`, `v(F)=-5`; uses GGV1 Props 1.13/2.1 | `PROOF_72_108.md:121-141`; `t6_premises_verify.py` | ran |
| normal form | `ell_{1,0}(P)=R^2`, `ell_{1,0}(Q)=R^3`, `R=x^4 y^7(y+1)` | `upstream_facts.json:33-35`; `PROOF_72_108.md:121` | cited |
| Cassou-Noguès | not used; zero hits for `Cassou` or `bisector` | repo-wide grep | ran |
| proof status | 23k-word prose proof, revision 2; no formal proof of Theorem A | `PROOF_72_108.md:1-5`; `lean/LEAN_CERTIFICATE.md:19-24` | cited |
| registry C0 | `closed: true`, level `claimed`, `subcases_closed: 5` | `d2/proof_dag.json:45124-45128` | ran |
| repo's stale status | `C0 (OPEN) ... Not proven` | `d2/CURRENT_STATUS.md:45-47` | cited |
| zero-margin inputs | 5, incl. `v_t(S)>=11`, `v_t(h6)>=11`, `deg R<=12` | `PROOF_72_108.md:90-103` | cited |
| upper bound | `a_t<=9` single-legged; second route "not established" | `PROOF_72_108.md:80-89`; `at_le9_audit.py` 76/0 | ran |
| Lean | no Theorem A: `f31` identity + corner-144 bridge; two `sorry`s in challenge statements | `lean/Cert/Certificate.lean`; `lean/Challenge.lean:185,194` | cited |
| missing docs | `C0_CLOSEOUT.md`, `AT_LE9_AUDIT.md` cited but absent | `PROOF_72_108.md:104-114,1027`; `Test-Path` | ran |
| proof hash | `3b87d82` absent from 44-commit public history | `git cat-file -t 3b87d82` | ran |
| repo checkers | 55/55, 48/48, 14/14, 20/20, 70/70, 81/81, 76/0, f37 closure, d3 all pass | `d2/*.py`, `d3/verify.py` | ran |
| independent math | ODE residual 0; K-syzygy 0; `deg d2=16`; `disc(q)=12837954459480883200` | `research/scratch/audit_*.py` | ran |
| open frontier | 24 pairs at `max_deg >= 125` (heading says 27) | `README.md:183-191` | cited |
| priority | Helali 2026-07-21, DOI 10.5281/zenodo.21479814; 89 MB certificate; conditional on Prop. 4.3 | Zenodo API record 21479814 | cited |
| archived wstrinz | v1.0.0 metadata: "target claim C0 remains OPEN" | Zenodo API record 21633408 | cited |

## Reproduced locally

```text
python research/scratch/audit_lemma22.py
  ODE residual: 0
  q irreducible; disc(q)=12837954459480883200; q(-1)=3315
  unique solution of degree <=14; deg f1=14

python research/scratch/audit_ksyz.py
  K-syzygy residual: 0
  G1,G2,G3 after the Prop. 4.4 chart: 0
  determinant of the (d0,d1,T) system: 27*e^3*(-R^2+S*e)/4

python research/scratch/audit_deg16.py
  first nonzero remainder coefficient at u^18 => deg d2 = 16
  both sign branches; coefficient 3976272238669575811/3746994889972252672

cd d2
python support_certificates.py   -> 55/55 checks pass
python verify_derivation.py     -> ALL 48 DERIVATION CHECKS PASSED
python t6_premises_verify.py    -> ALL 14 PREMISE CHECKS PASSED
python prop43_audit.py          -> 20/20 checks passed
PYTHONUTF8=1 python caps_audit.py -> 70/70 checks pass
python spine9_audit.py          -> 81/81 checks pass
python at_le9_audit.py          -> 76 checks, 0 failures
python f37_sat_verify.py        -> CONCLUSION: f31 in <G1,G2,G3,G5body+Phi>
cd ../d3; python verify.py      -> all 9 checks pass; JC(3) refuted
```

## What I could not find

- No peer-reviewed version and no external referee report. GGHV22 itself is arXiv:2204.14178v1, 25 pp, no journal reference on arXiv; the whole (72,108) exclusion is unrefereed.
- No independent verification of Theorem A. The in-repo "independent" auditors are same-project; the only external check on the same exclusion is Helali's certificate, not downloaded or replayed (89 MB, conditional on Prop. 4.3).
- No line-by-line verification of §§6-8. I ran the repo's checkers and independent spot checks (ODE, syzygy, caps, Cor. 4.5, degree 16), not a full re-derivation of the cascade.
- No authoritative resolution of the status split: `README.md` claims exclusion, while `d2/CURRENT_STATUS.md:45-47` and the archived Zenodo v1.0.0 record say `C0` is OPEN; the cited `C0_CLOSEOUT.md` and `AT_LE9_AUDIT.md` are absent from the public tree.
- No trace of the claimed revision `3b87d82` (proof header) in the 44-commit public history; `git cat-file -t 3b87d82` fails.
- No in-repo proof that the displayed Theorem A follows from its stated hypotheses alone: the transfer proof (`PROOF_72_108.md:378-391`) assumes the §5 normal form, which the theorem statement does not list. Corollary B is unaffected if every counterexample satisfies §5.
- No use of Cassou-Noguès/no-bisector (a negative finding, not a gap), and no machine-readable certificate for Theorem A as a whole.
