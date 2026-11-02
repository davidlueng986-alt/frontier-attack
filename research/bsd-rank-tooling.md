# Windows tooling for algebraic and analytic rank of elliptic curves over Q

## Answer

PARI/GP 2.19.0 is the native-Windows option: its installer provides `ellrank` (unconditional r1 ≤ rank ≤ r2, exact when equal), `ell2cover` (locally soluble 2-covers) and `ellanalyticrank` (tolerance-based numerical vanishing order); SageMath, mwrank and lcalc have no native Windows build. SageMath runs under WSL, where `E.rank()` uses mwrank/eclib by default and errors unless the rank is proved, while `E.analytic_rank()` uses PARI and is only “probably” the analytic rank. mwrank exists only as eclib source, and lcalc computes L-values/zeros plus analytic rank with `--with-pari`; neither ships Windows binaries. LMFDB is complete for N < 500,000 (plus prime N ≤ 3×10^8 and 7-smooth N), laptop 2-descent is routine to ~10^6 and feasible to 10^8–10^9 (rank 6 proved at N = 5,187,563,742), and numerical analytic rank is tested to N ≤ 10^9. No peer-reviewed algebraic-versus-analytic rank disagreement was found; for rank ≥ 4 the analytic rank is numerical only.

## Objects and numbers

| item | value | source |
|---|---|---|
| PARI/GP and Windows build | 2.19.0 (2026-09-28); `Pari64-2-19-0.exe`, 109.8 MB | [PARI download](https://pari.math.u-bordeaux.fr/download.html) |
| `ellrank` | `ellrank(E,{effort=0},{points})` → `[r1,r2,s,L]`; r1 ≤ rank E(Q) ≤ r2; r2 = C − T − s, computed unconditionally | [PARI docs](https://pari.math.u-bordeaux.fr/dochtml/ref-stable/Elliptic_curves.html) |
| `ellrank` and Sha | r1 = r2 proves the rank; equality is impossible when Sha has 4-torsion; E = [−113²,0] gives [0,2,0,[]], analytic rank 0 and Sha = (Z/4)² (assuming BSD) | [PARI docs](https://pari.math.u-bordeaux.fr/dochtml/ref-stable/Elliptic_curves.html) |
| `ell2cover` | basis of everywhere locally soluble 2-covers; each `[R,P]` is a quartic y² = R(x) plus a point on E(k) | [PARI docs](https://pari.math.u-bordeaux.fr/dochtml/ref-stable/Elliptic_curves.html) |
| `ellanalyticrank` | order of vanishing at s=1 plus first nonzero derivative; values < eps treated as zero; default eps = 2^(−b/2) | [PARI docs](https://pari.math.u-bordeaux.fr/dochtml/ref-stable/Elliptic_curves.html) |
| numerical noise | `ellL1` on 389a1 at default precision gives −5.384e−29, although its analytic rank is 2 | [PARI docs](https://pari.math.u-bordeaux.fr/dochtml/ref-stable/Elliptic_curves.html) |
| Sage `E.rank()` | default `algorithm='mwrank_lib'` (eclib C++); options `mwrank_shell`, `pari`; “assuming no conjectures”; RuntimeError if it cannot prove the rank | [Sage docs](https://doc.sagemath.org/html/en/reference/arithmetic_curves/sage/schemes/elliptic_curves/ell_rational_field.html) |
| Sage `E.analytic_rank()` | default `algorithm='pari'`; “probably” the analytic rank; also `sympow`, `rubinstein` (lcalc), `magma`, `zero_sum`, `all` | [Sage docs](https://doc.sagemath.org/html/en/reference/arithmetic_curves/sage/schemes/elliptic_curves/ell_rational_field.html) |
| analytic rank ≥ 4 | open problem to prove it for any particular elliptic curve | [Sage docs](https://doc.sagemath.org/html/en/reference/arithmetic_curves/sage/schemes/elliptic_curves/ell_rational_field.html) |
| Sage on Windows | WSL/Ubuntu (or Docker/VM); no native Windows build | [Sage installation](https://doc.sagemath.org/html/en/installation/index.html) |
| mwrank on Windows | only inside eclib source (in Sage); no Windows binary; conda-forge eclib has no win-64; source build needs NTL + PARI | [mwrank page](https://johncremona.github.io/mwrank/index.html), [conda-forge eclib](https://api.anaconda.org/package/conda-forge/eclib) |
| lcalc | zeros and values of L-functions; analytic rank if built `--with-pari`; latest 2.2.2; conda-forge has no win-64 | [README](https://gitlab.com/sagemath/lcalc/-/raw/master/README.md), [tags](https://gitlab.com/sagemath/lcalc/-/tags), [conda-forge lcalc](https://api.anaconda.org/package/conda-forge/lcalc) |
| LMFDB coverage | all N < 500,000; prime N ≤ 3×10^8; 7-smooth N; max N 299,996,953; 2,086 rank 4, 19 rank 5 | [completeness](https://www.lmfdb.org/EllipticCurve/Q/Completeness), [stats](https://www.lmfdb.org/EllipticCurve/Q/stats) |
| rank 4 calibration | N = 234,446; [1,−1,0,−79,289] | [rank 4 records](https://www.math.uwaterloo.ca/~mrubinst/elliptic_curves/records/rank_4.html) |
| rank 5 calibration | N = 19,047,851; [0,0,1,−79,342] | [rank 5 records](https://www.math.uwaterloo.ca/~mrubinst/elliptic_curves/records/rank_5.html) |
| rank 6 calibration | N = 5,187,563,742; [1,1,0,−2582,48720] | [rank 6 records](https://www.math.uwaterloo.ca/~mrubinst/elliptic_curves/records/rank_6.html) |
| rank census | 2,302,875,461 curves of rank ≥ 3 with \|Δ\| ≤ 10^13; 99.973% proved by explicit points plus a 2-descent upper bound | [rank census](https://www.math.uwaterloo.ca/~mrubinst/elliptic_curves/index.html) |
| Sha-conditional status | “S”: rank_lo = rank_hi − 1; equals rank_hi if Sha is finite (2-primary part enough), otherwise rank_lo or rank_hi | [rank census](https://www.math.uwaterloo.ca/~mrubinst/elliptic_curves/index.html) |
| 1997 mwrank timings | log\|Δ\| 8.52–8.83: mean 50 s, worst 7,206 s; log\|Δ\| 12.90–15.34: mean 2,053 s, worst 141,682 s | [Djabri–Smart HPL-97-137](https://shiftleft.com/mirrors/www.hpl.hp.com/techreports/97/HPL-97-137.pdf) |
| laptop 2-descent (estimate) | routine N ≲ 10^6; feasible 10^8–10^9 at hours–days; no sharp cutoff (Selmer rank/Sha dominate) | extrapolated from [Djabri–Smart](https://shiftleft.com/mirrors/www.hpl.hp.com/techreports/97/HPL-97-137.pdf) and the rank-6 record |
| analytic-rank cost | N^(1/2+ε) terms to estimate L′(E,1) | [Elkies math/0208056](https://ar5iv.labs.arxiv.org/html/math/0208056) |
| analytic-rank tested range | zero-sum bound correct for 99.75% of curves with N ≤ 10^9 | [Sage docs](https://doc.sagemath.org/html/en/reference/arithmetic_curves/sage/schemes/elliptic_curves/ell_rational_field.html) |
| rank conjecture verified | all E/Q of conductor N < 130,000 (analytic rank ≤ 3) | [Miller, LMS JCM 2011](https://doi.org/10.1112/s1461157011000180) |
| algebraic vs analytic disagreement | none found in peer-reviewed sources; LMFDB, the rank census and Miller’s checks agree | [rank census](https://www.math.uwaterloo.ca/~mrubinst/elliptic_curves/index.html), [Miller](https://doi.org/10.1112/s1461157011000180) |

## What I could not find

- No published laptop benchmark gives a sharp 2-descent cutoff; the bands above extrapolate from Djabri–Smart and record computations, not measurements here.
- No official Windows binaries for eclib/mwrank, lcalc, or native Sage; only WSL/MSYS2/Cygwin source builds or conda-forge packages without win-64.
- No peer-reviewed proved-rank disagreement; apparent conflicts are 2-descent intervals (r1 < r2, e.g. Sha = (Z/4)²) or numerical-tolerance artifacts (`ellL1` on 389a1).
- Rank-7 records differ: OEIS gives 382,623,908,456, the census 536,670,340,706; ranks 4–6 agree.
