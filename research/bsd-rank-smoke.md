# Analytic rank of E/Q without PARI or Sage: a 20,000-term smoke test

## Answer

With only brute-force point counts mod $p$ ($p \le 20000$) and mpmath, the requested sharp 20,000-term Dirichlet truncation fails as a rank estimator on these three curves: its ladder returns $r=0$ for all three (right for 11a1 only by accident), because at $s=1$ the truncated sum converges to a cutoff-dependent offset ($L_{\rm naive}(1)=0.2095367$ vs $0.2538419$ for 11a1, $0.0657963$ vs $0$ for 37a1, $-0.0465231$ vs $0$ for 5077a1), not to $L(1)$. Replacing the sharp cutoff by the standard smoothed approximate functional equation (same $a_n$, same brute-force point counts, no PARI/Sage) gives tail slopes $0.000010$, $1.000005$ and $2.999985$, i.e. ranks $0,1,3$, matching LMFDB's analytic ranks $0,1,3$ with gap $0$. The 11a1 value matches LMFDB's $L(1)=0.25384186085591068433775892335$ to all 29 shown digits; 37a1's $L'(1)$ and 5077a1's $L'''(1)/3!$ match to $1.1\times10^{-11}$ and $3.6\times10^{-9}$ relative. Runtime is 4.03 s, 6.18 s and 34.44 s including the derivative cross-check (AFE-only estimator: 0.82 s, 1.55 s, 13.16 s); extrapolating from measured per-term costs, the AFE should stay practical to $N\approx10^{10}$ and break around $N\approx10^{12}$, while the sharp truncation is already broken at $N=37$ for rank 1.

## Objects and numbers

### Method (one line each)

- Point count: for each odd $p$, complete the square to $(2y+a_1x+a_3)^2=D(x)$ and sum the Legendre symbol of $D(x)$ over $x\in\mathbb F_p$ with a quadratic-residue table; $p=2$ by brute force. Then $a_p=p+1-\#E(\mathbb F_p)$.
- $a_n$: multiplicative from $a_p$, with $a_{p^k}=a_pa_{p^{k-1}}-p\,a_{p^{k-2}}$ at good $p$ and $a_{p^k}=a_p^k$ at multiplicative bad primes (all three test conductors are prime); additive bad primes would give $a_{p^k}=0$.
- Sharp truncation (A): $L_N(s)=\sum_{n\le20000}a_nn^{-s}$ at $s=1+\epsilon$, exactly as requested.
- Smoothed AFE (B): $\Lambda(s)=\sum_n a_n[c_n^{-s}\Gamma(s,c_n)+w\,c_n^{-(2-s)}\Gamma(2-s,c_n)]$ with $c_n=2\pi n/\sqrt N$, $\Lambda(s)=N^{s/2}(2\pi)^{-s}\Gamma(s)L(s)$, and $w$ determined by matching at $s=2$; this is exact and exponentially convergent.
- Rank: slope of $\log|L(1+\epsilon)|$ against $\log\epsilon$ over the smallest four ladder points, with central differences of $L$ at $s=1$ as an independent check.

### Test curves and LMFDB targets

| item | value | source | checked |
|---|---|---|---|
| 11a1 model | `[0,-1,1,-10,-20]`, conductor 11 | task; LMFDB `11.a1` (Cremona 11a2) is isogenous | ran: point counts; $L(1)$ matches |
| 37a1 model | `[0,0,1,-1,0]`, conductor 37 | LMFDB `37.a1` (Cremona 37a1) | ran: `ainvs` match API |
| 5077a1 model | `[0,0,1,-7,6]`, conductor 5077 | LMFDB `5077.a1` (Cremona 5077a1) | ran: `ainvs` match API |
| 5077a1 note | smallest conductor of an elliptic curve of rank 3 over Q | <https://www.lmfdb.org/EllipticCurve/Q/5077/a/1> | cited |
| LMFDB analytic rank | 11a1: 0; 37a1: 1; 5077a1: 3 | <https://www.lmfdb.org/EllipticCurve/Q/11/a/1>, <https://www.lmfdb.org/EllipticCurve/Q/37/a/1>, <https://www.lmfdb.org/EllipticCurve/Q/5077/a/1> | cited |
| LMFDB special value, 11a1 | $L(E,1)=0.25384186085591068433775892335$ | <https://www.lmfdb.org/EllipticCurve/Q/11/a/1> | cited |
| LMFDB special value, 37a1 | $L'(E,1)=0.30599977383405230182048368332$ | <https://www.lmfdb.org/EllipticCurve/Q/37/a/1> | cited |
| LMFDB special value, 5077a1 | $L^{(3)}(E,1)/3!=1.7318499001193006897919750851$ | <https://www.lmfdb.org/EllipticCurve/Q/5077/a/1> | cited |
| 389a1 extra calibration | conductor 389, analytic rank 2, $L''(1)/2!=0.75931650028842677023019260790$ | <https://www.lmfdb.org/EllipticCurve/Q/389/a/1> | cited |

### Rank estimates

The gap is defined as $r_{\rm est}-r_{\rm LMFDB}$.

| curve | $N$ | LMFDB $r_{\rm an}$ | naive $r_{\rm est}$ | naive gap | AFE $r_{\rm est}$ | AFE gap |
|---|---|---|---|---|---|---|
| 11a1 | 11 | 0 | 0 | 0 | 0 | 0 |
| 37a1 | 37 | 1 | 0 | -1 | 1 | 0 |
| 5077a1 | 5077 | 3 | 0 | -3 | 3 | 0 |

The naive $r_{\rm est}$ column is the literal tail-fit output; as the next section shows, for 37a1 and 5077a1 it is a truncation artifact, not evidence of rank 0.

| item | value | source | checked |
|---|---|---|---|
| 11a1, naive 20k truncation | tail slope $0.000029$, $r_{\rm est}=0$, gap $0$ | `analytic_rank.py` | ran |
| 11a1, smoothed AFE | tail slope $0.000010$, $r_{\rm est}=0$, gap $0$ | `analytic_rank.py` | ran |
| 37a1, naive 20k truncation | tail slope $-0.000043$, $r_{\rm est}=0$, gap $-1$ | `analytic_rank.py` | ran |
| 37a1, smoothed AFE | tail slope $1.000005$, $r_{\rm est}=1$, gap $0$ | `analytic_rank.py` | ran |
| 5077a1, naive 20k truncation | tail slope $-0.000080$, $r_{\rm est}=0$, gap $-3$ | `analytic_rank.py` | ran |
| 5077a1, smoothed AFE | tail slope $2.999985$, $r_{\rm est}=3$, gap $0$ | `analytic_rank.py` | ran |
| 389a1, smoothed AFE (extra) | tail slope $1.999995$, $r_{\rm est}=2$, gap $0$ | `test_389.py` | ran |

### High-accuracy special values from the AFE

| item | value | source | checked |
|---|---|---|---|
| 11a1 $L(1)$ | `0.25384186085591068433775892335090946104` | `special_values.py`, mp.dps=45 | ran |
| 37a1 $L'(1)$ | `0.3059997738306325` (central difference, $h=10^{-5}$) | `special_values.py`, mp.dps=45 | ran |
| 5077a1 $L'''(1)/6$ | `1.7318499064195` (central difference, $h=3\times10^{-5}$) | `special_values.py`, mp.dps=45 | ran |
| 389a1 $L''(1)/2!$ | $7.5931607\times10^{-13}/(10^{-6})^2=0.75931607$ (ladder at $\epsilon=10^{-6}$) | `test_389.py`, mp.dps=30 | ran |

The two positive-rank values also satisfy the BSD quotient $\Omega\cdot\mathrm{Reg}$ from the LMFDB pages: $5.986917292463919\times0.051111408239969=0.305999773834052$ for 37a1 and $4.151687983086933\times0.417143558758384=1.731849900119301$ for 5077a1. The AFE values above match both products, so the agreement is not just a numerical coincidence with one quoted constant.

### Where the sharp truncation fails

| item | value | source | checked |
|---|---|---|---|
| $L_{\rm naive}(1)$, $N_{\rm cut}=20000$ | 11a1 `0.2095367245066492`; 37a1 `0.0657963393295491`; 5077a1 `-0.0465231381793374` | `naive_at_one.py` | ran |
| true $L(1)$ | 11a1 `0.2538418608559107`; 37a1 `0`; 5077a1 `0` | LMFDB / AFE | cited + ran |
| 11a1 naive error at $s=1$ | $-0.0443051363492615$ ($-17.45\%$) | `naive_at_one.py` | ran |
| 11a1 partial sums at $s=1$ | 10k `0.275394997691`; 20k `0.209536724507`; 50k `0.225233665909`; 100k `0.267508289991` | `naive_convergence.py` | ran |
| 37a1 partial sums at $s=1$ | 10k `-0.0507864537302`; 20k `0.0657963393295`; 50k `-0.0511055245358`; 100k `-0.0241301834944` | `naive_convergence.py` | ran |
| 37a1 naive local slopes | `0.5890, -0.0210, -0.0669, -0.0276, -0.0093, ...` | `analytic_rank.py` | ran |
| 5077a1 naive local slopes | `0.3721, -0.6350, -0.1781, -0.0566, -0.0177, ...` | `analytic_rank.py` | ran |
| 5077a1 naive ladder at $\epsilon=10^{-6}$ | `-0.046522685879905` (not $0$) | `analytic_rank.py` | ran |

### Runtime and scaling

| curve | point count to 20k | $a_n$ build + $w$ | naive ladder | AFE ladder | derivative check | total |
|---|---|---|---|---|---|---|
| 11a1 | 0.287 s | 0.085 s | 1.648 s | 0.742 s | 1.266 s | 4.03 s |
| 37a1 | 0.384 s | 0.128 s | 2.427 s | 1.420 s | 1.821 s | 6.18 s |
| 5077a1 | 0.364 s | 0.376 s | 2.342 s | 12.783 s | 18.572 s | 34.44 s |

| item | value | source | checked |
|---|---|---|---|
| point counts $p\le20000$ | 2262 primes, $\sum p=21{,}171{,}191$, $0.287$--$0.384$ s | `analytic_rank.py` | ran |
| point counts $p\le50000$ | 5133 primes, $\sum p=121{,}013{,}308$, $2.400$ s | `bench_scaling.py` | ran |
| AFE evaluation, $N=5077$, $n_{\max}=955$ | $1.07$ s per evaluation; $12.783$ s for 12 ladder points | `analytic_rank.py` | ran |
| AFE evaluation scaling, actual $a_n$ | $n_{\max}=1000$: 1.011 s; 5000: 1.531 s; 10000: 2.247 s | `bench_scaling.py` | ran |
| AFE evaluation scaling, dummy $a_n$ | $n_{\max}=100000$: 17.633 s; 200000: 37.212 s | `bench_afe2.py` | ran |
| AFE-only estimator runtime | 11a1 0.82 s; 37a1 1.55 s; 5077a1 13.16 s | `analytic_rank.py` timings + `bench_scaling.py` | ran |
| full smoke-test runtime | 11a1 4.03 s; 37a1 6.18 s; 5077a1 34.44 s | `analytic_rank.py` summary | ran |
| break, sharp truncation | already at $N=37$ for rank 1; at $N=11$ (rank 0) $L(1)$ off by $17.45\%$ | `analytic_rank.py` | ran |
| break, smoothed AFE | $\sim N=10^{10}$ (about 1 h) to $\sim N=10^{12}$ (about 1--1.5 days), extrapolated | measured per-term costs | extrap. |
| method assumptions | conductor $N$ supplied; bad primes multiplicative ($a_p=\pm1$) or additive with $a_p=0$; no Tate algorithm | `analytic_rank.py` docstring | -- |
| software | Python 3.14.2, numpy 2.4.4, mpmath 1.3.0 | `analytic_rank.py` header | ran |

## Reproduced locally

### Main estimator: `research/scratch/analytic_rank.py`

```powershell
PS C:\Users\daive\Downloads\frontier-attack> python research/scratch/analytic_rank.py
```

```text
analytic_rank.py -- minimal analytic-rank estimator (no PARI/Sage)
python 3.14.2  numpy 2.4.4  mpmath 1.3.0  mp.dps=30

Method A: sharp Dirichlet truncation at 20000 terms (as requested).
Method B: smoothed approximate functional equation, same a_n.

==============================================================================
11a1   ainvs=[0, -1, 1, -10, -20]   N=11   LMFDB 11.a1 (analytic_rank=0)
note: task model; isogenous to LMFDB 11.a1 (Cremona 11a2), same L-function
==============================================================================
self-test vs O(p^2) brute force, p<=50: OK
primes counted up to 20000 (2262 primes); point-count time 0.287s; a_n build 0.006s
root number: w=+1   |AFE_w(2)-Dirichlet(2)|: w=+1 2.220e-06, w=-1 3.336e-01
AFE truncation: nmax=200, c_max=2*pi*nmax/sqrt(N)=378.9

  eps        L_naive(1+eps)        L_AFE(1+eps)
  --------------------------------------------------------------
  3.0e-01         0.34428356986336         0.34655580006734
  1.0e-01         0.26832927891317         0.28479018025744
  3.0e-02         0.23019277366173         0.26311232822937
  1.0e-02         0.21680143981517         0.25693004238272
  3.0e-03         0.21175961027791         0.25476808741992
  1.0e-03         0.21028195052012         0.25415058068123
  3.0e-04         0.20976074404021         0.25393447443463
  1.0e-04         0.20961144083427         0.25387273182255
  3.0e-05          0.2095591439397         0.25385112212212
  1.0e-05         0.20954419808301         0.25384494794238
  3.0e-06         0.20953896662492         0.25384278698161
  1.0e-06         0.20953747188373         0.25384216956446

  naive local slopes: [0.2269, 0.1273, 0.0546, 0.0195, 0.0064, 0.0021, 0.0006, 0.0002, 0.0001, 0.0000, 0.0000]
  AFE   local slopes: [0.1787, 0.0658, 0.0216, 0.0070, 0.0022, 0.0007, 0.0002, 0.0001, 0.0000, 0.0000, 0.0000]
  naive tail-fit slope (last 4 eps): 0.000029 -> r_est=0
  AFE   tail-fit slope (last 4 eps): 0.000010 -> r_est=0

  central differences at s=1 (D_k ~ L^(k)(1) for k <= r):
    h=1e-2: D_0=0.2538418609  D_1=0.3087076163  D_2=0.02265787965  D_3=-0.2202401134
    h=1e-3: D_0=0.2538418609  D_1=0.3087085248  D_2=0.02265612615  D_3=-0.2202412542

  timings: naive ladder 1.648s, AFE ladder 0.742s, derivatives 1.266s, w-detect 0.079s

==============================================================================
37a1   ainvs=[0, 0, 1, -1, 0]   N=37   LMFDB 37.a1 (analytic_rank=1)
==============================================================================
self-test vs O(p^2) brute force, p<=50: OK
primes counted up to 20000 (2262 primes); point-count time 0.384s; a_n build 0.006s
root number: w=-1   |AFE_w(2)-Dirichlet(2)|: w=+1 1.922e-01, w=-1 3.300e-06
AFE truncation: nmax=200, c_max=2*pi*nmax/sqrt(N)=206.6

  eps        L_naive(1+eps)        L_AFE(1+eps)
  --------------------------------------------------------------
  3.0e-01         0.10843866543461          0.1050635994384
  1.0e-01        0.056777617800334        0.032330449021519
  3.0e-02         0.05823305197351        0.009344206352772
  1.0e-02         0.06267277368247       0.0030785158895185
  3.0e-03        0.064790501167731      0.00091967455961722
  1.0e-03        0.065454323047468      0.00030618618485598
  3.0e-04        0.065693021023257       9.1816717758731e-5
  1.0e-04        0.065761831755843       3.0601842724588e-5
  3.0e-05        0.065785979896813       9.1801611043458e-6
  1.0e-05        0.065792885503146       3.0600163929835e-6
  3.0e-06        0.065795303109995       9.1800100042864e-7
  1.0e-06        0.065795993916209       3.0599996038171e-7

  naive local slopes: [0.5890, -0.0210, -0.0669, -0.0276, -0.0093, -0.0030, -0.0010, -0.0003, -0.0001, -0.0000, -0.0000]
  AFE   local slopes: [1.0728, 1.0310, 1.0106, 1.0035, 1.0011, 1.0004, 1.0001, 1.0000, 1.0000, 1.0000, 1.0000]
  naive tail-fit slope (last 4 eps): -0.000043 -> r_est=0
  AFE   tail-fit slope (last 4 eps): 1.000005 -> r_est=1

  central differences at s=1 (D_k ~ L^(k)(1) for k <= r):
    h=1e-2: D_0=0.0  D_1=0.3059963541  D_2=0.3730988156  D_3=-0.8207208842
    h=1e-3: D_0=0.0  D_1=0.3059997396  D_2=0.3730956267  D_3=-0.8207484997

  timings: naive ladder 2.427s, AFE ladder 1.420s, derivatives 1.821s, w-detect 0.122s

==============================================================================
5077a1   ainvs=[0, 0, 1, -7, 6]   N=5077   LMFDB 5077.a1 (analytic_rank=3)
==============================================================================
self-test vs O(p^2) brute force, p<=50: OK
primes counted up to 20000 (2262 primes); point-count time 0.364s; a_n build 0.016s
root number: w=-1   |AFE_w(2)-Dirichlet(2)|: w=+1 2.129e-02, w=-1 2.784e-06
AFE truncation: nmax=955, c_max=2*pi*nmax/sqrt(N)=84.2

  eps        L_naive(1+eps)        L_AFE(1+eps)
  --------------------------------------------------------------
  3.0e-01        0.024318162791271        0.026835784337042
  1.0e-01        -0.01615897445193       0.0014383290758705
  3.0e-02       -0.034709534695897       4.4230576963765e-5
  1.0e-02       -0.042211334719261       1.7000700137143e-6
  3.0e-03       -0.045185783963849       4.6500948767032e-8
  1.0e-03        -0.04607302639713       1.7286467937945e-9
  3.0e-04        -0.04638764533267      4.6733986271576e-11
  1.0e-04       -0.046477929994197      1.7315293375606e-12
  3.0e-05       -0.046509571108972      4.6757350587736e-14
  1.0e-05       -0.046518615382882      1.7318178413434e-15
  3.0e-06       -0.046521781294231      4.6759687625549e-17
  1.0e-06       -0.046522685879905      1.7318466942164e-18

  naive local slopes: [0.3721, -0.6350, -0.1781, -0.0566, -0.0177, -0.0057, -0.0018, -0.0006, -0.0002, -0.0001, -0.0000]
  AFE   local slopes: [2.6636, 2.8919, 2.9662, 2.9892, 2.9966, 2.9989, 2.9997, 2.9999, 3.0000, 3.0000, 3.0000]
  naive tail-fit slope (last 4 eps): -0.000080 -> r_est=0
  AFE   tail-fit slope (last 4 eps): 2.999985 -> r_est=3

  central differences at s=1 (D_k ~ L^(k)(1) for k <= r):
    h=1e-2: D_0=0.0  D_1=4.329799754e-5  D_2=-0.0006411976634  D_3=10.39529911
    h=1e-3: D_0=0.0  D_1=4.3296265e-7  D_2=-6.411812832e-6  D_3=10.3911414

  timings: naive ladder 2.342s, AFE ladder 12.783s, derivatives 18.572s, w-detect 0.360s

==============================================================================
SUMMARY
==============================================================================
curve          N  r_LMFDB  r_naive   r_AFE  gap_naive  gap_AFE  t_total
-----------------------------------------------------------------------
11a1          11        0        0       0         +0       +0    4.03s
37a1          37        1        0       1         -1       +0    6.18s
5077a1      5077        3        0       3         -3       +0   34.44s

results written to C:\Users\daive\Downloads\frontier-attack\research\scratch\analytic_rank_results.json
```

The central-difference check is consistent with the ladder: for 37a1 the first stable nonzero derivative is $D_1$; for 5077a1, $D_1$ and $D_2$ shrink by a factor 100 when $h$ is reduced from $10^{-2}$ to $10^{-3}$ (pure truncation error), while $D_3$ stays at $10.39$, so the first stable nonzero derivative is $D_3$ (and $D_3/6$ is compared with LMFDB's $L'''(1)/3!$ below).

### High-accuracy special values (mp.dps = 45)

```powershell
PS C:\Users\daive\Downloads\frontier-attack> python research/scratch/special_values.py
```

```text
11a1 N 11 w 1 nmax 300
  L(1) = 0.25384186085591068433775892335090946104   LMFDB 0.25384186085591068433775892335
37a1 N 37 w -1 nmax 300
  h 1e-3 L'(1) = 0.305999739636187689743464338578   LMFDB 0.30599977383405230182048368332
  h 1e-4 L'(1) = 0.30599977349207364419373650224   LMFDB 0.30599977383405230182048368332
  h 1e-5 L'(1) = 0.305999773830632515243065613863   LMFDB 0.30599977383405230182048368332
5077a1 N 5077 w -1 nmax 1525
  h 1e-3 L'''(1)/6 = 1.73185690034313754755044272026   LMFDB 1.7318499001193006897919750851
  h 1e-4 L'''(1)/6 = 1.73184997012160927231633871248   LMFDB 1.7318499001193006897919750851
  h 3e-5 L'''(1)/6 = 1.73184990641950852030501502655   LMFDB 1.7318499001193006897919750851
```

### Naive partial sums at $s=1$ do not converge at 20k (or 100k)

```powershell
PS C:\Users\daive\Downloads\frontier-attack> python research/scratch/naive_convergence.py
```

```text
11a1 ap time 13.3
11a1 naive L(1) by cutoff:
    10000  0.275394997691  (0.02s)
    20000  0.209536724507  (0.04s)
    50000  0.225233665909  (0.09s)
   100000  0.267508289991  (0.18s)
37a1 ap time 11.75
37a1 naive L(1) by cutoff:
    10000  -0.0507864537302  (0.02s)
    20000  0.0657963393295  (0.03s)
    50000  -0.0511055245358  (0.10s)
   100000  -0.0241301834944  (0.16s)
```

The partial sums oscillate with the cutoff; 20,000 is not a special point, and increasing the cutoff does not monotonically approach $L(1)$.

### Scaling measurements

```powershell
PS C:\Users\daive\Downloads\frontier-attack> python research/scratch/bench_scaling.py
```

```text
point counting scaling (5077a1 coefficients):
  P=   500  nprimes=   95  sum_p=      21536  time=0.002 s
  P=  1000  nprimes=  168  sum_p=      76127  time=0.003 s
  P=  2000  nprimes=  303  sum_p=     277050  time=0.010 s
  P=  5000  nprimes=  669  sum_p=    1548136  time=0.033 s
  P= 10000  nprimes= 1229  sum_p=    5736396  time=0.094 s
  P= 20000  nprimes= 2262  sum_p=   21171191  time=0.341 s
  P= 50000  nprimes= 5133  sum_p=  121013308  time=2.400 s
AFE eval scaling (s=1.0001, N=5077, w=-1):
  nmax=   200  c_max=  17.64  time=0.132 s  L=1.5974063e-12
  nmax=   500  c_max=  44.09  time=0.351 s  L=1.7315293e-12
  nmax=  1000  c_max=  88.18  time=1.011 s  L=1.7315293e-12
  nmax=  2000  c_max= 176.36  time=1.138 s  L=1.7315293e-12
  nmax=  5000  c_max= 440.91  time=1.531 s  L=1.7315293e-12
  nmax= 10000  c_max= 881.81  time=2.247 s  L=1.7315293e-12
```

Dummy-$a_n$ run (isolates the incomplete-gamma cost at large $n_{\max}$):

```powershell
PS C:\Users\daive\Downloads\frontier-attack> python research/scratch/bench_afe2.py
```

```text
nmax=  10000 c_max=   881.8 time=2.666 s
nmax=  20000 c_max=  1763.6 time=3.712 s
nmax=  50000 c_max=  4409.1 time=7.721 s
nmax= 100000 c_max=  8818.1 time=17.633 s
nmax= 200000 c_max= 17636.3 time=37.212 s
```

Extrapolation used for the break estimate: at $N=5077$, $n_{\max}=955$ and 12 AFE evaluations take 12.783 s. The large-$n_{\max}$ dummy run is linear at about $1.9$--$2.3\times10^{-4}$ s per term. For $N\approx10^{10}$, $n_{\max}\approx12\sqrt N=1.2\times10^6$, so 12 evaluations are about 45--55 min plus about 17 min of point counting; for $N\approx10^{12}$, $n_{\max}\approx1.2\times10^7$, so 12 evaluations are about 7--9 h plus about 24 h of point counting. This is an extrapolation from measured per-term costs, not a run at those conductors.

### LMFDB records used

```powershell
PS C:\Users\daive\Downloads\frontier-attack> python research/scratch/lmfdb_api3.py
```

```text
11.a1 {'conductor': 11, 'ainvs': [0, -1, 1, -7820, -263580], 'analytic_rank': 0, 'rank': 0, 'regulator': 1, 'sha': 1, 'torsion': 1, 'lmfdb_label': '11.a1', 'Clabel': '11a2'}
37.a1 {'conductor': 37, 'ainvs': [0, 0, 1, -1, 0], 'analytic_rank': 1, 'rank': 1, 'regulator': {'__RealLiteral__': 0, 'data': '0.051111408239968840235886099757', 'prec': 97}, 'sha': 1, 'torsion': 1, 'lmfdb_label': '37.a1', 'Clabel': '37a1'}
5077.a1 {'conductor': 5077, 'ainvs': [0, 0, 1, -7, 6], 'analytic_rank': 3, 'rank': 3, 'regulator': {'__RealLiteral__': 0, 'data': '0.41714355875838396981711954462', 'prec': 97}, 'sha': 1, 'torsion': 1, 'lmfdb_label': '5077.a1', 'Clabel': '5077a1'}
```

The LMFDB `11.a1` record is the minimal model `[0,-1,1,-7820,-263580]` (Cremona 11a2), isogenous to the task's non-minimal `[0,-1,1,-10,-20]` (Cremona 11a1); isogenous curves share the same $L$-function, and the AFE's $L(1)$ agrees with LMFDB's $11.a1$ special value to 29 digits.

### Extra rank-2 calibration: 389a1

```powershell
PS C:\Users\daive\Downloads\frontier-attack> python research/scratch/test_389.py
```

```text
selftest []
w 1 nmax 336 cands {1: '0.3600928636', -1: '0.4305028659'} target 0.3600976311
  3.0e-01 0.0560379734
  1.0e-01 0.007147761132
  3.0e-02 0.000671620858
  1.0e-02 7.549945823e-5
  3.0e-03 6.822214777e-6
  1.0e-03 7.588860049e-7
  3.0e-04 6.83268653e-8
  1.0e-04 7.592734681e-9
  3.0e-05 6.833732319e-10
  1.0e-05 7.593121972e-11
  3.0e-06 6.833836884e-12
  1.0e-06 7.5931607e-13
slopes [1.8744, 1.9642, 1.9894, 1.9967, 1.999, 1.9997, 1.9999, 2.0, 2.0, 2.0, 2.0]
tail 1.9999953391648864
```

## What I could not find

- No reliable way was found to make the sharp 20,000-term Dirichlet truncation estimate rank $\ge1$ on these curves: at $s=1$ its partial sums oscillate with the cutoff (11a1: 0.2754, 0.2095, 0.2252, 0.2675 at 10k/20k/50k/100k; 37a1: -0.0508, 0.0658, -0.0511, -0.0241), so 20k is neither sufficient nor special. This is a property of conditional convergence at $s=1$, not an implementation bug.
- No implementation of Tate's algorithm: the script takes the conductor $N$ as an input from LMFDB. Computing $N$ from $(a_1,a_2,a_3,a_4,a_6)$, and detecting all local reduction types, is out of scope. The bad-prime treatment here is exact for the three test curves (prime conductor, multiplicative) and for additive primes with $a_p=0$, but not a general local-factor implementation.
- No direct run at $N>5077$: the break-conductor numbers ($N\approx10^{10}$--$10^{12}$) are extrapolated from measured point-count and incomplete-gamma costs, not measured at those conductors. The extra 389a1 check ($N=389$, rank 2, gap 0) is a calibration, not a large-conductor test.
- No independent proof of the AFE formula was attempted; it is derived from $\Lambda(s)=w\Lambda(2-s)$ in the script docstring and validated numerically against LMFDB's special values. The root number $w$ is determined numerically by matching the AFE to the Dirichlet sum at $s=2$, not from local root-number theory.
- LMFDB's `11.a1` is not the same Weierstrass model as the task's `11a1`; the comparison relies on the isogeny class (same $L$-function). The numerical $L(1)$ match to 29 digits is the check that this is the right $L$-function, but the exact `ainvs` differ (LMFDB: `[0,-1,1,-7820,-263580]`; task: `[0,-1,1,-10,-20]`).
