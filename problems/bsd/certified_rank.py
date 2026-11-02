#!/usr/bin/env python3
"""
certified_rank.py -- a numerical certificate for the analytic order of
vanishing of L(E, s) at s = 1, built on the smoothed approximate functional
equation (AFE).

The point of this script
------------------------
The AFE estimator in research/scratch/analytic_rank.py returns a numerical
order of vanishing (the slope of log|L(1+eps)| against log(eps)) with no error
bar.  A slope is not a proof.  Here the same AFE is equipped with

  (1) an explicit truncation-tail bound, derived from the incomplete gamma
      function and checked numerically by doubling the number of terms;
  (2) a rounding bound obtained by evaluating at two precisions;
  (3) for positive order, a finite-difference/Richardson estimate of the
      derivative and of its truncation error.

That is enough to certify k = 0 for 11a1 and k = 1 for 37a1, and to explain
exactly why k = 3 for 5077a1 is *not* certifiable by this method.

The AFE, and why it is the right object
---------------------------------------
For a weight-2 newform f = sum a_n q^n of level N and root number w,

    Lambda(s) = N^{s/2} (2 pi)^{-s} Gamma(s) L(s),
    Lambda(s) = w Lambda(2-s).

With c_n = 2 pi n / sqrt(N) and g(t) = f(i t / sqrt(N)), modularity gives
g(1/t) = w t^2 g(t), hence the exact representation

    Lambda(s) = sum_{n>=1} a_n [ c_n^{-s} Gamma(s, c_n)
                                 + w c_n^{-(2-s)} Gamma(2-s, c_n) ],   (AFE)

    L(s)      = (2 pi / sqrt(N))^s / Gamma(s) * Lambda(s),

where Gamma(s, c) is the upper incomplete gamma function.  This is the
standard smoothed approximate functional equation.  Every term decays like
exp(-c_n), so the sum converges exponentially and a finite truncation at
n <= M is a genuinely good approximation.

Truncation-tail bound
---------------------
Write R_M(s) for the tail n > M of (AFE).  For 0 <= sigma = Re(s) <= 2 and
c >= 1, integrating by parts gives

    Gamma(s, c) <= Gamma(sigma, c) <= 2 c^{sigma-1} e^{-c}.

(The constant 2 comes from Gamma(sigma,c) = c^{sigma-1} e^{-c} *
int_0^inf (1+u)^{sigma-1} e^{-c u} du and (1+u)^{sigma-1} <= 1+u for
sigma in [0,2].)  Both gamma factors in (AFE) therefore contribute at most
2 c_n^{-1} e^{-c_n}, and Deligne's bound |a_n| <= d(n) sqrt(n), together with
d(n) <= 2 sqrt(n), gives |a_n| <= 2 n.  Hence, with a = 2 pi / sqrt(N),

    |R_M(s)| <= (4 sqrt(N) / pi) * sum_{n>M} e^{-a n}
             = (4 sqrt(N) / pi) * e^{-a(M+1)} / (1 - e^{-a}).      (TAIL)

The bound is uniform in sigma in [0,2].  For L rather than Lambda it is
multiplied by max |(2 pi/sqrt(N))^s / Gamma(s)| over the s-range used.

Rounding and derivative bounds
------------------------------
mpmath does not provide interval arithmetic.  The rounding error is estimated
by evaluating the same quantity at two precisions (50 and 70 decimal digits)
and by an a priori condition-number estimate (the sum of absolute values of
the terms, scaled by 10^-dps).  This is a *heuristic* bound, not a proof.

For k >= 1 the derivatives are computed with central differences at three
step sizes h, h/2, h/4 and two Richardson extrapolations.  The difference
between the last two extrapolants estimates the remaining truncation error.
This is also heuristic (it is a standard numerical error estimate, not a
Cauchy-integral bound).

Honesty statement
-----------------
The result is therefore labelled "numerical with heuristic bound".  The tail
bound (TAIL) is rigorous given the AFE identity, modularity, and Deligne's
bound; the rounding and derivative bounds are not.  The root number w is
determined numerically by matching (AFE) at s = 3 to the absolutely
convergent Dirichlet series; it is not obtained from local root-number
theory.  These caveats are printed and written to the JSON output.

What is certified
-----------------
* 11a1, k = 0: the interval for L(1) excludes 0.  Since L(1) != 0, the order
  of vanishing is 0.  No lower-order information is needed.
* 37a1, k = 1: w = -1, so Lambda(1) = -Lambda(1) = 0 and L(1) = 0 exactly
  (given w = -1; w is determined numerically, see the honesty statement).
  The interval for L'(1) excludes 0, so the order is exactly 1.
* 5077a1, k = 3: NOT certified.  L'''(1) != 0 gives order <= 3 and w = -1
  gives odd order, so the order is in {1, 3}; but L'(1) = L''(1) = 0 cannot
  be established numerically.  Any correct error interval for a quantity that
  is actually zero must contain zero, so no amount of extra precision or
  extra terms can prove the vanishing.  A structural argument is required.

Run:
    python problems/bsd/certified_rank.py

Writes:
    runs/bsd_certified.json
"""

from __future__ import annotations

import json
import math
import platform
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import mpmath as mp
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "runs"

DPS_LO = 50            # low precision, for the two-precision rounding check
DPS_HI = 70            # production precision
DPS_LADDER = 45        # rank-ladder precision
DPS_TAIL_CHECK = 60    # tail-doubling verification precision
DERIV_H = "0.02"       # base step for the derivative stencils
DIR_NMAX = 2000        # terms in the Dirichlet cross-check at s = 3


# --------------------------------------------------------------------------
# curves
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Curve:
    label: str
    ainvs: tuple[int, int, int, int, int]
    conductor: int
    lmfdb_label: str
    lmfdb_rank: int
    special_name: str
    special_value: str
    target_derivative: int      # 0 -> L(1), 1 -> L'(1), 3 -> L'''(1)/3!
    expected_root_number: int
    note: str = ""


CURVES = (
    Curve(
        "11a1", (0, -1, 1, -10, -20), 11, "11.a1", 0,
        "L(1)", "0.25384186085591068433775892335", 0, +1,
        "task model; isogenous to LMFDB 11.a1 (Cremona 11a2), same L-function",
    ),
    Curve(
        "37a1", (0, 0, 1, -1, 0), 37, "37.a1", 1,
        "L'(1)", "0.30599977383405230182048368332", 1, -1,
    ),
    Curve(
        "5077a1", (0, 0, 1, -7, 6), 5077, "5077.a1", 3,
        "L'''(1)/3!", "1.7318499001193006897919750851", 3, -1,
        "smallest conductor of an elliptic curve of rank 3 over Q",
    ),
)


# --------------------------------------------------------------------------
# elementary arithmetic: point counts and a_n
# --------------------------------------------------------------------------

def primes_upto(n: int) -> list[int]:
    if n < 2:
        return []
    sieve = np.ones(n + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i:: i] = False
    return np.nonzero(sieve)[0].tolist()


def brute_count(ainvs, p: int) -> int:
    """O(p^2) affine count plus infinity; used only for the p <= 50 self-test."""
    a1, a2, a3, a4, a6 = ainvs
    cnt = 1
    for x in range(p):
        rhs = (x * x * x + a2 * x * x + a4 * x + a6) % p
        for y in range(p):
            if (y * y + a1 * x * y + a3 * y - rhs) % p == 0:
                cnt += 1
    return cnt


def count_points(ainvs, p: int) -> int:
    """#E(F_p) by a quadratic-residue table, O(p) per prime."""
    a1, a2, a3, a4, a6 = ainvs

    if p == 2:
        cnt = 1
        for x in range(2):
            for y in range(2):
                if (y * y + a1 * x * y + a3 * y
                        - (x * x * x + a2 * x * x + a4 * x + a6)) % 2 == 0:
                    cnt += 1
        return cnt

    B2 = (a1 * a1 + 4 * a2) % p
    B4 = (2 * a4 + a1 * a3) % p
    B6 = (a3 * a3 + 4 * a6) % p

    x = np.arange(p, dtype=np.int64)
    D = (4 * x * x * x + B2 * x * x + 2 * B4 * x + B6) % p

    leg = np.zeros(p, dtype=np.int8)
    y = np.arange((p + 1) // 2, dtype=np.int64)
    leg[(y * y) % p] = 1
    leg[leg == 0] = -1
    leg[0] = 0

    return p + 1 + int(leg[D].sum())


def selftest(ainvs, pmax: int = 50) -> list[tuple[int, int, int]]:
    bad = []
    for p in primes_upto(pmax):
        b = brute_count(ainvs, p)
        f = count_points(ainvs, p)
        if b != f:
            bad.append((p, b, f))
    return bad


def build_an(N: int, nmax: int, aps: dict[int, int]) -> list[int]:
    """Multiplicative a_n for n <= nmax from the a_p.

    Good p: a_{p^k} = a_p a_{p^{k-1}} - p a_{p^{k-2}}.
    Multiplicative bad p: a_{p^k} = a_p^k.  Additive bad p: a_p = 0 -> all 0.
    This is exact for the three test curves (prime conductor, multiplicative
    reduction); computing N itself (Tate's algorithm) is out of scope.
    """
    spf = list(range(nmax + 1))
    for i in range(2, int(nmax ** 0.5) + 1):
        if spf[i] == i:
            for j in range(i * i, nmax + 1, i):
                if spf[j] == j:
                    spf[j] = i

    a = [0] * (nmax + 1)
    a[1] = 1
    for n in range(2, nmax + 1):
        p = spf[n]
        m, k = n, 0
        while m % p == 0:
            m //= p
            k += 1
        if N % p == 0:
            val = 0 if aps[p] == 0 else aps[p] ** k
        else:
            ap = aps[p]
            v0, v1 = 1, ap
            for _ in range(2, k + 1):
                v0, v1 = v1, ap * v1 - p * v0
            val = v1
        a[n] = a[m] * val
    return a


def discriminant(ainvs) -> int:
    a1, a2, a3, a4, a6 = ainvs
    b2 = a1 * a1 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3 * a3 + 4 * a6
    b8 = a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4
    return -b2 * b2 * b8 - 8 * b4 ** 3 - 27 * b6 * b6 + 9 * b2 * b4 * b6


def prime_factors(n: int) -> dict[int, int]:
    n = abs(n)
    out: dict[int, int] = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


# --------------------------------------------------------------------------
# the smoothed AFE
# --------------------------------------------------------------------------

def afe_parts(s, a: list[int], N: int, M: int):
    """Return (S1, S2, Q1, Q2) with

        Lambda(s) = S1 + w S2,

    and Q1, Q2 the sums of absolute values of the two families of terms.
    Q1+Q2 (times the prefactor) is the condition number used for rounding.

    The caller is responsible for the mpmath precision (mp.workdps).
    """
    c0 = 2 * mp.pi / mp.sqrt(N)
    s1 = mp.mpf(0)
    s2 = mp.mpf(0)
    q1 = mp.mpf(0)
    q2 = mp.mpf(0)
    for n in range(1, M + 1):
        an = a[n]
        if an == 0:
            continue
        c = c0 * n
        t1 = c ** (-s) * mp.gammainc(s, c)
        t2 = c ** (-(2 - s)) * mp.gammainc(2 - s, c)
        s1 += an * t1
        s2 += an * t2
        q1 += abs(an) * abs(t1)
        q2 += abs(an) * abs(t2)
    pre = c0 ** s / mp.gamma(s)
    return pre * s1, pre * s2, abs(pre) * q1, abs(pre) * q2


def L_afe(s, a: list[int], N: int, w: int, M: int):
    """L(s) from the AFE truncated at n = M."""
    s1, s2, _, _ = afe_parts(s, a, N, M)
    return s1 + w * s2


def prefactor(s, N: int):
    """|(2 pi/sqrt(N))^s / Gamma(s)|, the Lambda -> L conversion factor."""
    c0 = 2 * mp.pi / mp.sqrt(N)
    return abs(c0 ** s / mp.gamma(s))


def tail_bound_lambda(N: int, M: int):
    """(TAIL): bound on |Lambda(s) - Lambda_M(s)| for 0 <= Re(s) <= 2."""
    a = 2 * mp.pi / mp.sqrt(N)
    return (4 * mp.sqrt(N) / mp.pi) * mp.e ** (-a * (M + 1)) / (1 - mp.e ** (-a))


def detect_root_number(a: list[int], N: int, M_afe: int, M_dir: int = DIR_NMAX,
                       dps: int = DPS_LADDER) -> dict:
    """Determine w in {+1,-1} by matching the AFE to the Dirichlet sum at s=3.

    At s = 3 the Dirichlet series converges absolutely, so the comparison is
    meaningful.  The two candidates are S1 + S2 and S1 - S2; the gap between
    them is reported as the discrimination margin.
    """
    with mp.workdps(dps):
        s = mp.mpf(3)
        s1, s2, _, _ = afe_parts(s, a, N, M_afe)
        target = mp.fsum(a[n] * mp.power(n, -3)
                         for n in range(1, M_dir + 1) if a[n])
        cand = {+1: s1 + s2, -1: s1 - s2}
        err = {w: abs(cand[w] - target) for w in cand}
        w = min(err, key=lambda key: err[key])
        gap = abs(cand[+1] - cand[-1])
        margin = gap / err[w] if err[w] != 0 else mp.inf
        return {
            "w": w,
            "dirichlet_target": mp.nstr(target, 25),
            "afe_plus": mp.nstr(cand[+1], 25),
            "afe_minus": mp.nstr(cand[-1], 25),
            "err_selected": mp.nstr(err[w], 8),
            "err_other": mp.nstr(err[-w], 8),
            "gap": mp.nstr(gap, 8),
            "margin": mp.nstr(margin, 8),
            "M_afe": M_afe,
            "M_dir": M_dir,
        }


# --------------------------------------------------------------------------
# derivatives by central differences + Richardson
# --------------------------------------------------------------------------

def derivative_table(a: list[int], N: int, w: int, M: int, h_str: str,
                     dps: int) -> dict:
    """L^(k)(1) for k = 0,1,2,3, with heuristic error components.

    The stencils are evaluated on a shared grid at steps h, h/2, h/4:

        k=1:  (f(h)-f(-h))/(2h)
        k=2:  (f(h)-2f(0)+f(-h))/h^2
        k=3:  (f(3h/2)-3f(h/2)+3f(-h/2)-f(-3h/2))/h^3

    then two Richardson extrapolations (h^2, then h^4 elimination).  The
    reported truncation error is 4*|R2 - R1(h/2)|, a conservative version of
    the usual Richardson difference estimate.

    All returned quantities are mpf values computed at the requested dps.
    """
    with mp.workdps(dps):
        h = mp.mpf(h_str)
        cache: dict = {}

        def F(off):
            if off not in cache:
                s = mp.mpf(1) + off
                s1, s2, q1, q2 = afe_parts(s, a, N, M)
                cache[off] = (s1 + w * s2, q1 + q2)
            return cache[off]

        def D(k: int, level: int):
            hl = h / 2 ** level
            if k == 1:
                return (F(hl)[0] - F(-hl)[0]) / (2 * hl)
            if k == 2:
                return (F(hl)[0] - 2 * F(mp.mpf(0))[0] + F(-hl)[0]) / hl ** 2
            if k == 3:
                p1 = hl / 2
                p3 = 3 * hl / 2
                return (F(p3)[0] - 3 * F(p1)[0]
                        + 3 * F(-p1)[0] - F(-p3)[0]) / hl ** 3
            raise ValueError(k)

        tail_lam = tail_bound_lambda(N, M)
        coefsum = {1: mp.mpf(1), 2: mp.mpf(4), 3: mp.mpf(8)}
        h_finest = h / 4

        raw = {}
        for k in (1, 2, 3):
            dvals = [D(k, level) for level in range(3)]
            r1 = [(4 * dvals[l + 1] - dvals[l]) / 3 for l in range(2)]
            r2 = (16 * r1[1] - r1[0]) / 15
            trunc = 4 * abs(r2 - r1[1])
            raw[k] = (r2, trunc)

        pref_max = max(prefactor(mp.mpf(1) + off, N) for off in cache)
        scale_max = max(v[1] for v in cache.values())

        out = {}
        for k in (1, 2, 3):
            r2, trunc = raw[k]
            # The R2 combination has coefficients (64, -20, 1)/45 on the
            # three levels; 4 is a safe overall factor for the propagation.
            round_ap = (mp.mpf(2000) * mp.mpf(10) ** (-dps) * scale_max
                        * coefsum[k] / h_finest ** k)
            tail_der = (pref_max * tail_lam * coefsum[k]
                        / h_finest ** k * 4)
            out[k] = {
                "value": r2,
                "trunc": trunc,
                "round_apriori": round_ap,
                "tail": tail_der,
            }

        v0, q0 = F(mp.mpf(0))
        out[0] = {
            "value": v0,
            "trunc": mp.mpf(0),
            "round_apriori": mp.mpf(200) * mp.mpf(10) ** (-dps) * q0,
            "tail": pref_max * tail_lam,
        }
        return out


# --------------------------------------------------------------------------
# rank ladder (the original estimator, for k_est)
# --------------------------------------------------------------------------

def ladder_slope(a: list[int], N: int, w: int, M: int,
                 dps: int = DPS_LADDER) -> dict:
    """Slope of log|L(1+eps)| against log(eps) over the last four ladder points."""
    with mp.workdps(dps):
        eps = [mp.mpf(s) for s in
               ("1e-1", "3e-2", "1e-2", "3e-3", "1e-3", "3e-4")]
        vals = [L_afe(1 + e, a, N, w, M) for e in eps]
        xs = [mp.log(e) for e in eps[-4:]]
        ys = [mp.log(abs(v)) for v in vals[-4:]]
        n = len(xs)
        mx = sum(xs) / n
        my = sum(ys) / n
        num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
        den = sum((x - mx) ** 2 for x in xs)
        slope = num / den
        return {
            "slope": float(slope),
            "eps": [float(e) for e in eps],
            "values": [mp.nstr(v, 18) for v in vals],
        }


# --------------------------------------------------------------------------
# tail-doubling verification
# --------------------------------------------------------------------------

def verify_tail(a: list[int], N: int, w: int, M_list: list[int],
                dps: int = DPS_TAIL_CHECK) -> list[dict]:
    """Check |L_M(s) - L_2M(s)| <= bound(M) at s = 1.05 (sigma in [0,2])."""
    with mp.workdps(dps):
        s = mp.mpf("1.05")
        floor = mp.mpf(10) ** (-dps + 2)
        rows = []
        for M in M_list:
            v1 = L_afe(s, a, N, w, M)
            v2 = L_afe(s, a, N, w, 2 * M)
            diff = abs(v1 - v2)
            bound = prefactor(s, N) * tail_bound_lambda(N, M)
            meaningful = bound > floor
            if meaningful and diff > bound:
                raise AssertionError(
                    f"tail bound violated at M={M}: |diff|={mp.nstr(diff, 8)} "
                    f"> bound={mp.nstr(bound, 8)}")
            rows.append({
                "M": M,
                "L_M": mp.nstr(v1, 20),
                "abs_diff": mp.nstr(diff, 8),
                "bound": mp.nstr(bound, 8),
                "bound_over_diff": (mp.nstr(bound / diff, 6)
                                    if diff != 0 else "inf"),
                "meaningful": bool(meaningful),
            })
        return rows


def derivative_label(k: int) -> str:
    return {0: "L(1)", 1: "L'(1)", 2: "L''(1)", 3: "L'''(1)"}[k]


# --------------------------------------------------------------------------
# one curve
# --------------------------------------------------------------------------

def run_curve(curve: Curve) -> dict:
    t_start = time.perf_counter()
    mp.mp.dps = 30

    M_afe = max(200, int(12 * math.sqrt(curve.conductor)) + 100)
    M_build = max(2 * M_afe, DIR_NMAX, curve.conductor)

    # --- point counts and a_n ---------------------------------------------
    t0 = time.perf_counter()
    bad = selftest(curve.ainvs, 50)
    if bad:
        raise AssertionError(f"{curve.label}: point-count self-test failed: {bad}")
    aps = {p: p + 1 - count_points(curve.ainvs, p)
           for p in primes_upto(M_build)}
    t_point = time.perf_counter() - t0

    t0 = time.perf_counter()
    a = build_an(curve.conductor, M_build, aps)
    t_build = time.perf_counter() - t0

    # --- sanity checks on the supplied conductor --------------------------
    disc = discriminant(curve.ainvs)
    fac = prime_factors(disc)
    if set(fac) != {curve.conductor}:
        raise AssertionError(
            f"{curve.label}: discriminant {disc} has bad primes {sorted(fac)}, "
            f"but the conductor was supplied as {curve.conductor}")
    if abs(a[curve.conductor]) != 1:
        raise AssertionError(
            f"{curve.label}: expected multiplicative reduction at "
            f"{curve.conductor}, got a_N = {a[curve.conductor]}")

    # --- root number ------------------------------------------------------
    t0 = time.perf_counter()
    wdiag = detect_root_number(a, curve.conductor, M_afe)
    t_w = time.perf_counter() - t0
    w = wdiag["w"]
    if w != curve.expected_root_number:
        raise AssertionError(
            f"{curve.label}: numerical root number {w:+d} disagrees with "
            f"the expected {curve.expected_root_number:+d}")

    # --- rank ladder ------------------------------------------------------
    t0 = time.perf_counter()
    ladder = ladder_slope(a, curve.conductor, w, M_afe)
    t_ladder = time.perf_counter() - t0
    k_est = int(round(ladder["slope"]))
    if not 0 <= k_est <= 4:
        k_est = None

    # --- derivative table at two precisions -------------------------------
    t0 = time.perf_counter()
    tab_lo = derivative_table(a, curve.conductor, w, M_afe, DERIV_H, DPS_LO)
    t_deriv_lo = time.perf_counter() - t0

    t0 = time.perf_counter()
    tab_hi = derivative_table(a, curve.conductor, w, M_afe, DERIV_H, DPS_HI)
    t_deriv_hi = time.perf_counter() - t0

    # --- combine -> derivative values and total bounds --------------------
    with mp.workdps(DPS_HI):
        deriv = {}
        for k in (0, 1, 2, 3):
            d_hi = tab_hi[k]
            d_lo = tab_lo[k]
            delta = abs(d_hi["value"] - d_lo["value"])
            round_bound = max(d_hi["round_apriori"], 100 * delta)
            total = d_hi["tail"] + round_bound + d_hi["trunc"]
            val = d_hi["value"]
            deriv[k] = {
                "value": val,
                "bound": total,
                "lo": val - total,
                "hi": val + total,
                "contains_zero": bool((val - total) <= 0 <= (val + total)),
                "exact_zero": bool(k == 0 and w == -1),
                "components": {
                    "tail": d_hi["tail"],
                    "rounding": round_bound,
                    "rounding_empirical": 100 * delta,
                    "rounding_apriori": d_hi["round_apriori"],
                    "derivative_truncation": d_hi["trunc"],
                },
            }

        k_first_nonzero = None
        for k in (0, 1, 2, 3):
            if not deriv[k]["contains_zero"]:
                k_first_nonzero = k
                break

        # --- tail-doubling verification -----------------------------------
        if curve.conductor <= 11:
            check_Ms = [10, 20, 40]
        elif curve.conductor <= 37:
            check_Ms = [25, 50, 100]
        else:
            check_Ms = [119, 238, 477, 955]
        tail_rows = verify_tail(a, curve.conductor, w, check_Ms)

        # --- certification logic ------------------------------------------
        target_k = curve.target_derivative
        exact_zeros = {0} if w == -1 else set()
        lower_ok = all(k in exact_zeros for k in range(target_k))
        target_excludes = not deriv[target_k]["contains_zero"]
        certified = bool(lower_ok and target_excludes
                         and k_est == target_k)

        if certified:
            if target_k == 0:
                basis = "L(1) interval excludes 0, hence ord = 0."
            else:
                basis = (
                    f"L(1) = 0 by the functional equation with w = {w:+d} "
                    f"(w determined numerically at s=3, margin "
                    f"{wdiag['margin']}); the interval for "
                    f"{derivative_label(target_k)} excludes 0, hence "
                    f"ord = {target_k}.")
        else:
            if target_k == 3:
                basis = (
                    "NOT certified: ord >= 3 requires L(1) = L'(1) = "
                    "L''(1) = 0.  Only L(1) = 0 follows from the functional "
                    "equation (w = -1, determined numerically with margin "
                    f"{wdiag['margin']}).  The "
                    "AFE intervals for L'(1) and L''(1) contain 0, which is "
                    "consistent with vanishing but is not a proof of it.  "
                    "A numerical error bound cannot prove an exact zero: any "
                    "correct interval around a true zero must contain zero.  "
                    "Hence only ord <= 3 (from L'''(1) != 0) and odd order "
                    "(from w = -1) are established, i.e. ord is in {1, 3}.")
            else:
                basis = (
                    f"NOT certified: lower-order vanishing missing "
                    f"(need {sorted(set(range(target_k)) - exact_zeros)} "
                    f"to vanish exactly) or the target interval contains 0.")

        # --- special-value interval ---------------------------------------
        special_scale = 6 if target_k == 3 else 1
        if target_k == 3:
            special_value = deriv[3]["value"] / special_scale
            special_bound = deriv[3]["bound"] / special_scale
            special_lo = special_value - special_bound
            special_hi = special_value + special_bound
        else:
            special_value = deriv[target_k]["value"]
            special_bound = deriv[target_k]["bound"]
            special_lo = deriv[target_k]["lo"]
            special_hi = deriv[target_k]["hi"]

        # Bound components for the *special value* (divide by 6 for
        # L'''(1)/3!), for reporting only.
        comp = deriv[target_k]["components"]
        special_components = {
            key: val / special_scale for key, val in comp.items()
        }

        # Does the quoted LMFDB value agree, allowing for the fact that the
        # quote is a truncation at its last printed digit?
        lmfdb_val = mp.mpf(curve.special_value)
        lmfdb_decimals = len(curve.special_value.split(".")[1])
        lmfdb_ulp = mp.mpf(10) ** (-lmfdb_decimals)
        lmfdb_inside = bool(special_lo - lmfdb_ulp <= lmfdb_val
                            <= special_hi + lmfdb_ulp)

        # order bracket, when it adds information
        if certified:
            order_candidates = [target_k]
        elif w == -1 and k_first_nonzero == 3:
            order_candidates = [1, 3]
        else:
            order_candidates = None

        der_json = {}
        for k in (0, 1, 2, 3):
            der_json[str(k)] = {
                "value": mp.nstr(deriv[k]["value"], 40),
                "bound": mp.nstr(deriv[k]["bound"], 8),
                "interval": [mp.nstr(deriv[k]["lo"], 40),
                             mp.nstr(deriv[k]["hi"], 40)],
                "contains_zero": deriv[k]["contains_zero"],
                "exact_zero": deriv[k]["exact_zero"],
                "bound_components": {
                    key: mp.nstr(val, 8)
                    for key, val in deriv[k]["components"].items()
                },
            }

        result = {
            "curve": curve.label,
            "ainvs": list(curve.ainvs),
            "conductor": curve.conductor,
            "discriminant": str(disc),
            "discriminant_factorisation": {str(p): e for p, e in fac.items()},
            "lmfdb_label": curve.lmfdb_label,
            "lmfdb_analytic_rank": curve.lmfdb_rank,
            "lmfdb_special_name": curve.special_name,
            "lmfdb_special_value": curve.special_value,
            "root_number": w,
            "root_number_diagnostics": wdiag,
            "nmax_afe": M_afe,
            "nmax_built": M_build,
            "k_est": k_est,
            "k_est_ladder_slope": ladder["slope"],
            "k_est_ladder_values": ladder["values"],
            "k_first_nonzero_derivative": k_first_nonzero,
            "special_value_name": curve.special_name,
            "special_value": mp.nstr(special_value, 50),
            "special_value_scale": special_scale,
            "interval": {
                "lo": mp.nstr(special_lo, 50),
                "hi": mp.nstr(special_hi, 50),
            },
            "bound": mp.nstr(special_bound, 8),
            "bound_components": {
                key: mp.nstr(val, 8)
                for key, val in special_components.items()
            },
            "lmfdb_value_in_interval": lmfdb_inside,
            "lmfdb_quote_ulp": mp.nstr(lmfdb_ulp, 3),
            "certified": certified,
            "certification_basis": basis,
            "order_candidates": order_candidates,
            "derivatives": der_json,
            "tail_checks": tail_rows,
            "timings_s": {
                "point_count": t_point,
                "build_an": t_build,
                "root_number": t_w,
                "ladder": t_ladder,
                "derivatives_dps_lo": t_deriv_lo,
                "derivatives_dps_hi": t_deriv_hi,
            },
        }

    # values needed for printing, outside the high-precision block
    result["_print"] = {
        "special_name": curve.special_name,
        "special_value": mp.nstr(special_value, 50),
        "lo": mp.nstr(special_lo, 50),
        "hi": mp.nstr(special_hi, 50),
        "bound": mp.nstr(special_bound, 8),
        "components": {
            key: mp.nstr(val, 8)
            for key, val in special_components.items()
        },
        "lmfdb_value_in_interval": lmfdb_inside,
        "lmfdb_quote_ulp": mp.nstr(lmfdb_ulp, 3),
        "deriv": {
            k: {
                "value": mp.nstr(deriv[k]["value"], 30),
                "bound": mp.nstr(deriv[k]["bound"], 8),
                "contains_zero": deriv[k]["contains_zero"],
                "exact_zero": deriv[k]["exact_zero"],
            } for k in (0, 1, 2, 3)
        },
        "total_time": time.perf_counter() - t_start,
    }
    return result


# --------------------------------------------------------------------------
# reporting
# --------------------------------------------------------------------------

def print_curve(res: dict) -> None:
    p = res["_print"]
    wd = res["root_number_diagnostics"]
    print()
    print("=" * 82)
    print(f"{res['curve']}   N={res['conductor']}   "
          f"ainvs={res['ainvs']}   "
          f"LMFDB {res['lmfdb_label']} (analytic rank "
          f"{res['lmfdb_analytic_rank']})")
    print("=" * 82)
    print(f"  discriminant {res['discriminant']} = "
          f"{res['discriminant_factorisation']}; only bad prime "
          f"{res['conductor']}")
    print(f"  root number w={res['root_number']:+d}  "
          f"(s=3 match: err={wd['err_selected']}, gap={wd['gap']}, "
          f"margin={wd['margin']})")
    print(f"  k_est = {res['k_est']}  "
          f"(ladder slope {res['k_est_ladder_slope']:.6f}; "
          f"first non-zero derivative k={res['k_first_nonzero_derivative']})")

    print()
    print("  derivative / order table  (h=0.02, dps 50 -> 70)")
    print("  k   L^(k)(1)                         bound        contains 0   exact 0")
    print("  " + "-" * 76)
    for k in ("0", "1", "2", "3"):
        d = p["deriv"][int(k)]
        print(f"  {k}   {d['value']:<30}  {d['bound']:<11}  "
              f"{'yes' if d['contains_zero'] else 'no':<11}  "
              f"{'yes' if d['exact_zero'] else 'no'}")

    print()
    print(f"  interval for {p['special_name']}:")
    print(f"    [{p['lo']}, {p['hi']}]")
    c = p["components"]
    print(f"  bound = {p['bound']}  =  tail {c['tail']} + "
          f"rounding {c['rounding']} + deriv.trunc "
          f"{c['derivative_truncation']}")
    if res["special_value_scale"] != 1:
        print(f"    (all three components are scaled by 1/"
              f"{res['special_value_scale']} for {p['special_name']})")
    print(f"  LMFDB value {res['lmfdb_special_value']} consistent with the "
          f"interval (quote ulp {p['lmfdb_quote_ulp']}): "
          f"{'yes' if p['lmfdb_value_in_interval'] else 'no'}")
    print(f"  certified = {str(res['certified']).lower()}")
    print(f"  basis: {res['certification_basis']}")
    if res["order_candidates"] is not None and not res["certified"]:
        print(f"  order candidates from the certified facts: "
              f"{res['order_candidates']}")

    print()
    print(f"  tail-doubling check at s=1.05 (dps=60):")
    print("    M      |L_M - L_2M|    bound        bound/diff   meaningful")
    print("    " + "-" * 64)
    for row in res["tail_checks"]:
        print(f"    {row['M']:<6} {row['abs_diff']:<14} {row['bound']:<12} "
              f"{row['bound_over_diff']:<12} {'yes' if row['meaningful'] else 'no'}")
    print(f"  time: {p['total_time']:.2f}s")


def main() -> int:
    t_all = time.perf_counter()
    print("certified_rank.py -- numerical certificate for the analytic order")
    print("of vanishing from the smoothed approximate functional equation")
    print(f"python {platform.python_version()}  numpy {np.__version__}  "
          f"mpmath {mp.__version__}")
    print()
    print("AFE: Lambda(s) = sum_n a_n [c_n^-s Gamma(s,c_n) + w c_n^-(2-s) "
          "Gamma(2-s,c_n)]")
    print("     c_n = 2 pi n / sqrt(N),  "
          "L(s) = (2 pi/sqrt(N))^s / Gamma(s) * Lambda(s)")
    print("Tail: |R_M(s)| <= (4 sqrt(N)/pi) exp(-a(M+1))/(1-exp(-a)),  "
          "a = 2 pi/sqrt(N)")
    print("Rigour: numerical with heuristic bound.")
    print("  * tail bound: explicit and rigorous given the AFE identity,")
    print("    modularity, and |a_n| <= d(n) sqrt(n) (Deligne);")
    print("  * rounding: heuristic (two-precision agreement + a priori scale,")
    print("    no interval arithmetic);")
    print("  * derivative truncation: heuristic Richardson difference estimate;")
    print("  * root number: numerical match at s=3, not local root-number theory.")

    results = [run_curve(c) for c in CURVES]

    for res in results:
        print_curve(res)

    print()
    print("=" * 82)
    print("SUMMARY")
    print("=" * 82)
    print(f"{'curve':<9} {'k_est':>5}  {'special value':<14} "
          f"{'interval (center +/- bound)':<47} {'bound':<10} certified")
    print("-" * 112)
    for res in results:
        p = res["_print"]
        center = mp.nstr(mp.mpf(p["special_value"]), 22)
        iv = f"{center} +/- {p['bound']}"
        print(f"{res['curve']:<9} {str(res['k_est']):>5}  "
              f"{p['special_name']:<14} {iv:<47} "
              f"{p['bound']:<10} {str(res['certified']).lower()}")

    print()
    print("Verdict:")
    for res in results:
        if res["certified"]:
            print(f"  {res['curve']}: certified k = {res['k_est']} "
                  f"(numerical with heuristic bound)")
        else:
            print(f"  {res['curve']}: k_est = {res['k_est']} but NOT "
                  f"certified; see basis above")

    # strip the printing-only helper before serialising
    out = {
        "script": "problems/bsd/certified_rank.py",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "method": ("smoothed approximate functional equation with an explicit "
                   "incomplete-gamma truncation bound, two-precision rounding "
                   "estimate, and Richardson finite-difference derivatives"),
        "rigour": "numerical with heuristic bound",
        "rigour_notes": {
            "tail": ("explicit; rigorous given the AFE identity, modularity, "
                     "and |a_n| <= d(n) sqrt(n).  "
                     "|R_M(s)| <= (4 sqrt(N)/pi) exp(-a(M+1))/(1-exp(-a)), "
                     "a = 2 pi/sqrt(N), for 0 <= Re(s) <= 2."),
            "rounding": ("heuristic: two-precision agreement (dps 50 vs 70) "
                         "plus an a priori condition-number estimate; mpmath "
                         "does not provide interval arithmetic"),
            "derivative": ("heuristic: central differences at h, h/2, h/4 with "
                           "two Richardson extrapolations; the reported "
                           "truncation bound is 4*|R2 - R1(h/2)|"),
            "root_number": ("numerical: matched to the absolutely convergent "
                            "Dirichlet series at s=3; not obtained from local "
                            "root-number theory"),
        },
        "assumptions": [
            "the conductor N is supplied, not computed from the model",
            "bad primes are multiplicative with |a_p| = 1 (true for these "
            "three prime-conductor curves; checked via the discriminant)",
            "the AFE identity is used as a standard theorem (modularity + "
            "functional equation); it is not re-proved here",
        ],
        "precision": {"dps_low": DPS_LO, "dps_high": DPS_HI,
                      "dps_ladder": DPS_LADDER, "dps_tail_check": DPS_TAIL_CHECK},
        "curves": [],
    }
    for res in results:
        clean = {k: v for k, v in res.items() if not k.startswith("_")}
        out["curves"].append(clean)

    RUNS.mkdir(parents=True, exist_ok=True)
    path = RUNS / "bsd_certified.json"
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
        fh.write("\n")

    print()
    print(f"results written to {path}")
    print(f"total time {time.perf_counter() - t_all:.2f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
