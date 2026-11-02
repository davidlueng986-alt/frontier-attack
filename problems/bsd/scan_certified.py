#!/usr/bin/env python3
"""
scan_certified.py -- certifying analytic-rank scanner for prime-conductor
elliptic curves above LMFDB's complete range.

What this does
--------------
For each curve in the sample, decide whether the analytic order of vanishing
k = ord_{s=1} L(E,s) is 0 or 1, using the smoothed approximate functional
equation (AFE) with an explicit truncation-tail bound.  Then compare the
decision with the LMFDB analytic_rank (a DECLARED input, not a proof) and with
the parity of the LMFDB algebraic rank.

Why only k = 0 or k = 1
-----------------------
A numerical interval can prove that a value is NON-ZERO; it can never prove
that a value is EXACTLY zero, because every correct interval around a true
zero contains zero.  Therefore:

  * k = 0 is certifiable: exhibit an interval for L(1) that excludes 0.
  * k = 1 is certifiable: for prime conductor N, the root number is exactly
    w = a_N (standard local root-number theory: w_inf = -1 and
    w_N = -a_N for multiplicative reduction).  If a_N = -1 then
    L(1) = 0 exactly by the functional equation, and an interval for L'(1)
    that excludes 0 proves the order is exactly 1.
  * k >= 2 is not certifiable by this method: it requires proving exact
    vanishing of L(1) and L'(1), which no error interval can do.

This is the structural fact recorded in plan/bsd.md (B1).  The scanner is
therefore a coverage tool, not a counterexample search: an interval mismatch
is a LEAD, never a result.

Generalisation of the certifier
-------------------------------
`certify_k(ainvs, conductor, kmax=1)` is the requested general entry point.
It accepts any prime-conductor curve (minimal or not), not just the three
hardcoded curves of `certified_rank.py`.  The existing file is left untouched
(task constraint: do not modify existing files); this module imports its
AFE primitives (point counting, a_n, root-number cross-check, tail bound,
prefactor) and supplies the generalised function and the scanner.

Rigour
------
The truncation tail bound is explicit and rigorous given the AFE identity,
modularity, and Deligne's |a_n| <= d(n) sqrt(n).  The rounding bound is
a priori (mpmath has no interval arithmetic) and is cross-checked by a second
run at lower precision.  The derivative truncation error is a Richardson
difference estimate, which is heuristic.  The overall label is therefore
"numerical with heuristic bound".

Run:
    python problems/bsd/scan_certified.py

Reads:
    research/scratch/bsd_stress.json    (15 LMFDB curves, cached sample)

Writes:
    runs/bsd_stress.json                (per-curve results + verdict)
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))

import certified_rank as base  # noqa: E402  (path inserted above)

DPS_HI = 40          # production precision: the interval is quoted at dps 40
DPS_LO = 30          # cross-check precision for the rounding estimate
DERIV_H = "0.02"     # base step for the central-difference stencils
SAMPLE_DEFAULT = ROOT / "research" / "scratch" / "bsd_stress.json"
OUT_DEFAULT = ROOT / "runs" / "bsd_stress.json"
SAMPLE_LIMIT = 15    # defensive cap on the number of curves
CONDUCTOR_LO = 500_000
CONDUCTOR_HI = 520_000

# Calibration curves.  The first two are the original certifier's positive
# controls; the third is a rank-2 prime-conductor curve that must be REFUSED
# at kmax=1 (L(1) is numerically zero but w=+1, so no k=0 certificate is
# possible, and k=2 needs an exact-zero argument).
#
# 234446.a1 (the rank-4 curve from plan/bsd.md) is NOT usable here: its
# conductor is 234446 = 2 * 117223, composite, and this certifier is
# restricted to prime conductor.  That is recorded in the tooling block.
CALIBRATION = (
    ("11a1", (0, -1, 1, -10, -20), 11, 0),
    ("37a1", (0, 0, 1, -1, 0), 37, 1),
    ("389a1", (0, 1, 1, -2, 0), 389, None),
)


# --------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------

def is_prime(n: int) -> bool:
    """Deterministic trial division; n is at most ~5.2e5 in this scan."""
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def dec(x, digits: int = 30) -> str:
    return mp.nstr(x, digits)


def excludes_zero(lo, hi) -> bool:
    return bool(lo > 0 or hi < 0)


def contains_zero(lo, hi) -> bool:
    return bool(lo <= 0 <= hi)


def nmax_for(conductor: int) -> int:
    """AFE truncation point.

    The derivative tail is roughly
        4 * prefactor * tail_lambda / h_finest
    with tail_lambda ~ (4 sqrt(N)/pi) exp(-a M), a = 2 pi/sqrt(N).
    Choosing a M >= 25 makes the derivative tail ~1e-5 at N ~ 5e5, far
    below the observed Richardson truncation (~1e-6) and far below any
    expected L'(1).  The floor keeps the small control curves accurate.
    """
    return max(400, int(math.ceil(4.0 * math.sqrt(conductor))))


# --------------------------------------------------------------------------
# AFE derivative pair: L(1) and L'(1) with bound components
# --------------------------------------------------------------------------

def _deriv_pair(a, N: int, w: int, M: int, h_str: str, dps: int) -> dict:
    """L(1) and L'(1) at precision dps, with heuristic error components.

    Central differences at h, h/2, h/4 and two Richardson extrapolations.
    The reported truncation error is 4*|R2 - R1(h/2)|, the same conservative
    estimate used by certified_rank.derivative_table.
    """
    with mp.workdps(dps):
        h = mp.mpf(h_str)
        hf = h / 4
        cache: dict = {}

        def F(off):
            if off not in cache:
                s1, s2, q1, q2 = base.afe_parts(mp.mpf(1) + off, a, N, M)
                cache[off] = (s1 + w * s2, q1 + q2)
            return cache[off]

        for off in (mp.mpf(0), h, -h, h / 2, -h / 2, h / 4, -h / 4):
            F(off)

        dvals = []
        for level in range(3):
            hl = h / 2 ** level
            dvals.append((F(hl)[0] - F(-hl)[0]) / (2 * hl))
        r1 = [(4 * dvals[1] - dvals[0]) / 3,
              (4 * dvals[2] - dvals[1]) / 3]
        r2 = (16 * r1[1] - r1[0]) / 15
        trunc1 = 4 * abs(r2 - r1[1])

        tail_lam = base.tail_bound_lambda(N, M)
        pref_max = max(base.prefactor(mp.mpf(1) + off, N) for off in cache)
        scale_max = max(v[1] for v in cache.values())

        v0, q0 = F(mp.mpf(0))
        # Same a-priori rounding constants as certified_rank.py.
        round0_ap = 200 * mp.mpf(10) ** (-dps) * q0
        round1_ap = 2000 * mp.mpf(10) ** (-dps) * scale_max / hf
        tail0 = pref_max * tail_lam
        tail1 = pref_max * tail_lam / hf * 4

        return {
            "v0": v0,
            "v1": r2,
            "trunc1": trunc1,
            "round0_ap": round0_ap,
            "round1_ap": round1_ap,
            "tail0": tail0,
            "tail1": tail1,
            "tail_lam": tail_lam,
            "pref_max": pref_max,
            "scale_max": scale_max,
            "h_finest": hf,
        }


# --------------------------------------------------------------------------
# the generalised certifier
# --------------------------------------------------------------------------

def _certify_k_full(ainvs, conductor: int, kmax: int = 1,
                    dps_hi: int = DPS_HI, dps_lo: int = DPS_LO,
                    nmax: int | None = None) -> dict:
    """Full diagnostic version of certify_k.  Returns a dict."""
    if kmax not in (0, 1):
        raise ValueError(
            "kmax must be 0 or 1: k >= 2 is not certifiable by an interval "
            "method (an interval can prove non-zero, never exact zero); "
            "see plan/bsd.md (B1)")
    if len(ainvs) != 5:
        raise ValueError("ainvs must have five coefficients")
    ainvs = tuple(int(x) for x in ainvs)
    N = int(conductor)
    if not is_prime(N):
        raise ValueError(f"conductor {N} is not prime; this certifier is "
                         f"restricted to prime-conductor curves")
    disc = base.discriminant(ainvs)
    if disc == 0:
        raise ValueError("singular Weierstrass model (discriminant 0)")
    fac = base.prime_factors(disc)
    if set(fac) != {N}:
        raise ValueError(
            f"discriminant {disc} has bad primes {sorted(fac)}, not exactly "
            f"the prime conductor {N}")

    t0 = time.perf_counter()
    M = int(nmax) if nmax is not None else nmax_for(N)

    # a_n for n <= M.  Only primes p <= M are needed: the AFE is truncated
    # at n = M, and for prime conductor N > M the bad prime N never enters
    # the sum.
    aps = {p: p + 1 - base.count_points(ainvs, p)
           for p in base.primes_upto(M)}
    a = base.build_an(N, M, aps)
    t_an = time.perf_counter() - t0

    # Exact local root number for prime conductor: w = w_inf * w_N
    # = (-1) * (-a_N) = a_N.  a_N is an integer point count, not a
    # numerical match.
    a_N = N + 1 - base.count_points(ainvs, N)
    if abs(a_N) != 1:
        raise ValueError(
            f"expected multiplicative reduction at {N} (|a_N| = 1), "
            f"got a_N = {a_N}")
    w = int(a_N)

    # Independent numerical root number from the s=3 AFE/Dirichlet match.
    t0 = time.perf_counter()
    wdiag = base.detect_root_number(
        a, N, M, M_dir=min(base.DIR_NMAX, M), dps=min(DPS_LO, dps_hi))
    t_w = time.perf_counter() - t0

    # Two precisions: the high one supplies the value, the low one the
    # empirical part of the rounding bound.
    t0 = time.perf_counter()
    hi = _deriv_pair(a, N, w, M, DERIV_H, dps_hi)
    t_hi = time.perf_counter() - t0
    t0 = time.perf_counter()
    lo = _deriv_pair(a, N, w, M, DERIV_H, dps_lo)
    t_lo = time.perf_counter() - t0

    with mp.workdps(dps_hi):
        delta0 = abs(hi["v0"] - lo["v0"])
        delta1 = abs(hi["v1"] - lo["v1"])
        round0 = max(hi["round0_ap"], 100 * delta0)
        round1 = max(hi["round1_ap"], 100 * delta1)
        total0 = hi["tail0"] + round0
        total1 = hi["tail1"] + round1 + hi["trunc1"]
        lo0, hi0 = hi["v0"] - total0, hi["v0"] + total0
        lo1, hi1 = hi["v1"] - total1, hi["v1"] + total1

        # Decision.  For w = -1, L(1) = 0 EXACTLY by the functional
        # equation, so the L(1) interval must contain 0; if it does not,
        # the interval is inconsistent with the exact identity and we
        # refuse to certify rather than report a spurious zero.
        inconsistency = False
        if w == -1:
            if not contains_zero(lo0, hi0):
                inconsistency = True
            decisive = 1
            k = 1 if (kmax >= 1 and excludes_zero(lo1, hi1)) else None
            value, lo_v, hi_v = hi["v1"], lo1, hi1
            components = {"tail": hi["tail1"], "rounding": round1,
                          "truncation": hi["trunc1"], "total": total1}
        else:
            decisive = 0
            k = 0 if excludes_zero(lo0, hi0) else None
            value, lo_v, hi_v = hi["v0"], lo0, hi0
            components = {"tail": hi["tail0"], "rounding": round0,
                          "truncation": mp.mpf(0), "total": total0}

        certified = k is not None
        return {
            "k": k,
            "certified": certified,
            "interval": (lo_v, hi_v),
            "value": value,
            "bound_components": components,
            "bound_components_all": {
                "L1": {"tail": hi["tail0"], "rounding": round0,
                       "truncation": mp.mpf(0), "total": total0,
                       "interval": (lo0, hi0),
                       "contains_zero": contains_zero(lo0, hi0)},
                "L1prime": {"tail": hi["tail1"], "rounding": round1,
                            "truncation": hi["trunc1"], "total": total1,
                            "interval": (lo1, hi1),
                            "contains_zero": contains_zero(lo1, hi1)},
            },
            "decisive_derivative": decisive,
            "w": w,
            "w_numerical": int(wdiag["w"]),
            "w_diagnostics": wdiag,
            "a_N": int(a_N),
            "nmax": M,
            "inconsistency": inconsistency,
            "timings_s": {
                "a_n": t_an,
                "root_number": t_w,
                "derivatives_hi": t_hi,
                "derivatives_lo": t_lo,
                "total": time.perf_counter() - t0,
            },
        }


def certify_k(ainvs, conductor: int, kmax: int = 1):
    """Certify the analytic order of vanishing k for a prime-conductor curve.

    Parameters
    ----------
    ainvs : sequence of 5 ints
        Weierstrass coefficients [a1, a2, a3, a4, a6].
    conductor : int
        Prime conductor N.
    kmax : int
        0 or 1.  Only k = 0 and k = 1 are structurally certifiable.

    Returns
    -------
    (k, interval, certified, bound_components)
        k : int or None
            The certified order (0 or 1), or None if it cannot be certified.
        interval : (mpf, mpf)
            The interval for the decisive derivative: L(1) for k = 0,
            L'(1) for k = 1 (and for w = -1 refusals).
        certified : bool
        bound_components : dict
            tail / rounding / truncation / total for the decisive derivative.
    """
    r = _certify_k_full(ainvs, conductor, kmax=kmax)
    return r["k"], r["interval"], r["certified"], r["bound_components"]


# --------------------------------------------------------------------------
# sample loading and validation
# --------------------------------------------------------------------------

def load_sample(path: Path, limit: int = SAMPLE_LIMIT) -> list[dict]:
    with open(path, "r", encoding="utf-8") as fh:
        doc = json.load(fh)
    curves = doc.get("curves", [])
    if not isinstance(curves, list) or not curves:
        raise ValueError(f"{path}: no curves list")
    if len(curves) > limit:
        curves = curves[:limit]
    out = []
    for c in curves:
        N = int(c["conductor"])
        if not (CONDUCTOR_LO <= N <= CONDUCTOR_HI):
            raise ValueError(
                f"{c.get('label')}: conductor {N} outside "
                f"[{CONDUCTOR_LO}, {CONDUCTOR_HI}]")
        if not is_prime(N):
            raise ValueError(f"{c.get('label')}: conductor {N} is not prime")
        ka = int(c["analytic_rank"])
        if ka not in (0, 1):
            raise ValueError(
                f"{c.get('label')}: analytic_rank {ka} not in {{0,1}}")
        if len(c["ainvs"]) != 5:
            raise ValueError(f"{c.get('label')}: bad ainvs")
        out.append({
            "label": str(c["label"]),
            "conductor": N,
            "ainvs": tuple(int(x) for x in c["ainvs"]),
            "analytic_rank": ka,
            "rank": int(c["rank"]),
        })
    return out


# --------------------------------------------------------------------------
# calibration
# --------------------------------------------------------------------------

def run_calibration() -> list[dict]:
    rows = []
    for label, ainvs, N, expected in CALIBRATION:
        t = time.perf_counter()
        try:
            r = _certify_k_full(ainvs, N, kmax=1, nmax=400 if N < 1000 else None)
            got = r["k"]
            err = None
        except Exception as exc:  # noqa: BLE001 - report, do not abort
            r = None
            got = None
            err = f"{type(exc).__name__}: {exc}"
        rows.append({
            "label": label,
            "conductor": N,
            "expected_k": expected,
            "k": got,
            "ok": (got == expected and err is None),
            "error": err,
            "root_number": None if r is None else r["w"],
            "nmax": None if r is None else r["nmax"],
            "interval": None if r is None else {
                "lo": dec(r["interval"][0]), "hi": dec(r["interval"][1]),
                "contains_zero": contains_zero(*r["interval"]),
            },
            "seconds": round(time.perf_counter() - t, 2),
            "note": (
                "expected refusal: rank 2, but kmax=1 can only certify "
                "k=0 or k=1; L(1) is numerically zero and w=+1, so no "
                "k=0 certificate is possible"
                if expected is None else
                "positive control"),
        })
    return rows


# --------------------------------------------------------------------------
# main scan
# --------------------------------------------------------------------------

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", default=str(SAMPLE_DEFAULT))
    ap.add_argument("--out", default=str(OUT_DEFAULT))
    ap.add_argument("--limit", type=int, default=SAMPLE_LIMIT)
    ap.add_argument("--skip-calibration", action="store_true")
    args = ap.parse_args(argv)

    t_all = time.perf_counter()
    sample_path = Path(args.sample)
    out_path = Path(args.out)

    print("scan_certified.py -- certifying analytic-rank scanner")
    print(f"sample: {sample_path}")
    print(f"dps: {DPS_HI} (cross-check {DPS_LO}); h = {DERIV_H}")
    print("rigour: numerical with heuristic bound; LMFDB ranks are declared "
          "inputs, not proofs")
    print()

    try:
        curves = load_sample(sample_path, limit=args.limit)
    except Exception as exc:  # noqa: BLE001
        print(f"FATAL: could not load sample: {type(exc).__name__}: {exc}")
        return 2

    calibration = [] if args.skip_calibration else run_calibration()
    if calibration:
        print("calibration:")
        for row in calibration:
            print(f"  {row['label']:<12} expected k={row['expected_k']!s:<5} "
                  f"got k={row['k']!s:<5} ok={row['ok']} "
                  f"({row['seconds']}s)"
                  + (f"  error={row['error']}" if row['error'] else ""))
        print()

    results = []
    anomalies = []
    uncertified = []

    for idx, c in enumerate(curves, 1):
        t0 = time.perf_counter()
        try:
            r = _certify_k_full(c["ainvs"], c["conductor"], kmax=1)
            err = None
        except Exception as exc:  # noqa: BLE001 - a failure is data
            r = None
            err = f"{type(exc).__name__}: {exc}"

        elapsed = time.perf_counter() - t0
        if r is None:
            row = {
                "label": c["label"],
                "conductor": c["conductor"],
                "ainvs": list(c["ainvs"]),
                "lmfdb_analytic_rank": c["analytic_rank"],
                "lmfdb_rank": c["rank"],
                "k": None,
                "certified": False,
                "error": err,
                "flag": "UNCERTIFIED_ERROR",
                "seconds": round(elapsed, 2),
            }
            results.append(row)
            uncertified.append({"label": c["label"], "reason": err})
            print(f"[{idx:>2}/{len(curves)}] {c['label']:<12} "
                  f"ERROR {err}  ({elapsed:.1f}s)")
            continue

        k = r["k"]
        ka = c["analytic_rank"]
        kr = c["rank"]
        flags = []
        if k is None:
            flags.append("UNCERTIFIED")
        else:
            if k != ka:
                flags.append("LEAD_NUMERICAL_MISMATCH")
            if (kr % 2) != (k % 2):
                flags.append("LEAD_PARITY_MISMATCH")
        if r["w_numerical"] != r["w"]:
            flags.append("LEAD_ROOT_NUMBER_DIAGNOSTIC")
        if r["inconsistency"]:
            flags.append("LEAD_INTERVAL_CONTRADICTS_EXACT_ZERO")

        flag = "+".join(flags) if flags else "OK"
        row = {
            "label": c["label"],
            "conductor": c["conductor"],
            "ainvs": list(c["ainvs"]),
            "lmfdb_analytic_rank": ka,
            "lmfdb_rank": kr,
            "a_N": r["a_N"],
            "root_number_exact": r["w"],
            "root_number_numerical": r["w_numerical"],
            "root_number_margin": r["w_diagnostics"].get("margin"),
            "nmax": r["nmax"],
            "k": k,
            "certified": r["certified"],
            "decisive_derivative": r["decisive_derivative"],
            "value": dec(r["value"], 30),
            "interval": {
                "lo": dec(r["interval"][0], 30),
                "hi": dec(r["interval"][1], 30),
                "contains_zero": contains_zero(*r["interval"]),
            },
            "bound_components": {
                key: dec(val, 12)
                for key, val in r["bound_components"].items()
            },
            "bound_components_all": {
                name: {
                    "tail": dec(d["tail"], 12),
                    "rounding": dec(d["rounding"], 12),
                    "truncation": dec(d["truncation"], 12),
                    "total": dec(d["total"], 12),
                    "interval": [dec(d["interval"][0], 30),
                                 dec(d["interval"][1], 30)],
                    "contains_zero": d["contains_zero"],
                }
                for name, d in r["bound_components_all"].items()
            },
            "k_match_lmfdb_analytic": bool(k is not None and k == ka),
            "parity_match_lmfdb_rank": bool(
                k is not None and (kr % 2) == (k % 2)),
            "flag": flag,
            "note": (
                "certified k=0: L(1) interval excludes 0"
                if k == 0 else
                "certified k=1: a_N=-1 gives L(1)=0 exactly; "
                "L'(1) interval excludes 0"
                if k == 1 else
                "not certifiable with kmax=1 (requires an exact-zero "
                "argument for L(1), or k>=2)"),
            "timings_s": {key: round(val, 3)
                          for key, val in r["timings_s"].items()},
            "seconds": round(elapsed, 2),
        }
        results.append(row)
        if k is None:
            uncertified.append({
                "label": c["label"],
                "reason": "k not certifiable with kmax=1 "
                          f"(w={r['w']:+d}, L(1) interval "
                          f"[{dec(r['bound_components_all']['L1']['interval'][0], 8)}, "
                          f"{dec(r['bound_components_all']['L1']['interval'][1], 8)}])",
            })
        if flags:
            anomalies.append({"label": c["label"], "flags": flags})

        print(f"[{idx:>2}/{len(curves)}] {c['label']:<12} "
              f"N={c['conductor']} a_N={r['a_N']:+d} w={r['w']:+d} "
              f"k={k!s:<4} LMFDB(k,r)=({ka},{kr}) "
              f"flag={flag}  ({elapsed:.1f}s)")
        print(f"      interval for {'L(1)' if r['decisive_derivative']==0 else chr(76)+chr(39)+'(1)'}: "
              f"[{dec(r['interval'][0], 20)}, {dec(r['interval'][1], 20)}]")

    certified_count = sum(1 for row in results if row.get("certified"))
    anomaly_count = len(anomalies)
    uncertified_count = len(uncertified)

    if anomalies:
        first = anomalies[0]
        verdict = (f"ANOMALY: {first['label']} "
                   f"(LEAD: {'+'.join(first['flags'])}; numerical, "
                   f"not a result)")
    elif uncertified_count:
        verdict = (f"ANOMALY: {uncertified[0]['label']} "
                   f"(UNCERTIFIED: coverage gap, not a mathematical "
                   f"counterexample)")
    else:
        verdict = (f"COVERAGE EXTENDED: {certified_count} curves certified, "
                   f"0 anomalies")

    out = {
        "script": "problems/bsd/scan_certified.py",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "method": (
            "smoothed approximate functional equation with an explicit "
            "incomplete-gamma truncation bound, a priori rounding bound "
            "cross-checked at two precisions, and Richardson "
            "finite-difference derivatives"),
        "rigour": "numerical with heuristic bound",
        "structural_note": (
            "A numerical interval can prove a value is non-zero, never that "
            "it is exactly zero.  Only k=0 (L(1) interval excludes 0) and "
            "k=1 (a_N=-1 forces L(1)=0 exactly; L'(1) interval excludes 0) "
            "are certifiable.  k>=2 is not certifiable by this method."),
        "tooling": {
            "gp": "absent (the PowerShell 'gp' command is an alias for "
                  "Get-ItemProperty; no PARI/GP binary found)",
            "sage": "absent",
            "cypari2": "absent (ModuleNotFoundError)",
            "consequence": (
                "no algebraic-rank backend; LMFDB 'rank' and 'analytic_rank' "
                "are DECLARED inputs, not proofs"),
        },
        "inputs": {
            "sample": str(sample_path),
            "sample_source": (
                "LMFDB ec_curvedata / EllipticCurve search, conductor "
                "[500000, 520000], prime conductor, analytic_rank in {0,1}"),
            "lmfdb_complete_range": "N <= 500000 (Cremona); above that the "
                                    "rows are Stein-Watkins / BGR "
                                    "prime-conductor data",
            "declared_not_proved": True,
        },
        "precision": {
            "dps": DPS_HI,
            "dps_cross_check": DPS_LO,
            "derivative_step_h": DERIV_H,
            "nmax_rule": "max(400, ceil(4*sqrt(N)))",
        },
        "calibration": calibration,
        "calibration_ok": all(row["ok"] for row in calibration),
        "curves": results,
        "certified_count": certified_count,
        "uncertified_count": uncertified_count,
        "anomaly_count": anomaly_count,
        "anomalies": anomalies,
        "uncertified": uncertified,
        "verdict": verdict,
        "total_seconds": round(time.perf_counter() - t_all, 2),
    }

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
        fh.write("\n")

    print()
    print("=" * 78)
    print(f"certified_count = {certified_count}")
    print(f"anomaly_count = {anomaly_count}")
    print(f"uncertified_count = {uncertified_count}")
    print(f"verdict: {verdict}")
    print(f"results written to {out_path}")
    print(f"total time {time.perf_counter() - t_all:.1f}s")
    print("=" * 78)

    if calibration and not out["calibration_ok"]:
        print("CALIBRATION FAILED")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
