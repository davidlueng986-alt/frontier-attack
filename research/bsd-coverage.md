# BSD coverage: no rank backend, and a 25-curve analytic-rank cross-check at N ~ 5x10^5

## Answer

There is no PARI/GP, Sage or mwrank on this machine: `gp` is only PowerShell's
`Get-ItemProperty` alias, and `pip install --user cypari2` dies at metadata
generation because no C compiler is installed, so the requested
`ellrank`/`ellanalyticrank` test for 37a1 could not be run at all. LMFDB's
documented range-filter API returns 404 on this deployment, so the sample was
taken from the server-rendered search page (410 curves with conductor in
[500000, 520000]) and each record was fetched by exact label; 25 were kept,
covering analytic ranks 0, 1, 2 and 3. Re-running the existing AFE estimator
from `research/scratch/analytic_rank.py` on all 25 reproduced LMFDB's
`analytic_rank` exactly: 25 evaluated, 0 mismatches, 0 root-parity
inconsistencies. Total wall time was 1190.6 s because the first pass was cut by
the 900 s ceiling after 17 curves and a second pass finished the remaining 8;
this is estimator coverage, not a certified analytic rank and not a BSD result.

## Objects and numbers

### Task 1 - tooling audit (10 min cap)

| item | value | source | checked |
|---|---|---|---|
| `gp` on PATH | **no**; `Get-Command gp` returns `Alias -> Get-ItemProperty` | local shell | ran |
| `gp.exe` on disk | none under `C:\Program Files`, `C:\Program Files (x86)`, `cygwin64`, `msys64`, `scoop`, `chocolatey`, `AppData\Local\Programs`; no `*pari*` directory under `C:\` or `C:\Users\daive` (depth 2) | local search | ran |
| `sage` | **not found** | `Get-Command sage` | ran |
| `mwrank` | **not found** | `Get-Command mwrank` | ran |
| `import cypari2` | **fails**: `ModuleNotFoundError: No module named 'cypari2'` | `python -c` | ran |
| `pip install --user cypari2` | **fails after 24.2 s** at metadata generation | pip log | ran |
| exact blocker | meson: `Unknown compiler(s): [['icl'], ['cl'], ['cc'], ['gcc'], ['clang'], ['clang-cl'], ['pgcc']]`; `vswhere.exe` absent, so no MSVC either. cypari2 2.2.4 additionally needs PARI headers/library. | pip log | ran |
| native Windows installer | <https://pari.math.u-bordeaux.fr/pub/pari/windows/Pari64-2-19-0.exe> | PARI download page | cited |
| installer version / size | PARI/GP **2.19.0**, 2026-09-28, 109.8 MB | PARI download page | cited |
| installer sha256 | `70ff123c371fa071e77632832a914292fb1c8cf98b47e25a57223640bf31f7f7` | PARI download page | cited |
| system-wide install | **not attempted** (per the task's instruction) | -- | -- |
| `ellrank`/`ellanalyticrank` on 37a1 | **not run**: `gp` does not exist; no output to report | -- | -- |

Consequence: the algebraic side of the terminal object is unavailable here.
Any rank statement below is LMFDB's `analytic_rank` field, not a proof, and no
`r1 = r2` interval was produced.

### Task 2 - conductor sample (15 min cap)

| item | value | source | checked |
|---|---|---|---|
| requested API URL | `.../api/ec_curvedata/?_format=json&_fields=lmfdb_label,ainvs,conductor,rank,analytic_rank&conductor=500000-520000&analytic_rank=0-3` | task | -- |
| its HTTP status | **404** (`LMFDB Page Not Found`) | local request + DevTools Chrome | ran |
| API behaviour | every non-underscore argument is an exact-equality column; range/list values (`500000-520000`, `500000..520000`, `500000,520000`, JSON `$gte/$lte`) all 404 | probes | ran |
| API pagination | `_limit` >= 500 works, but `_offset` fails at 20000 (404); deep paging to N~5x10^5 is impossible | probes | ran |
| working route | server-rendered search table `/EllipticCurve/Q/?conductor=500000-520000` + one API call per exact `lmfdb_label` | DevTools Chrome | ran |
| matches in conductor range, all ranks | **410** (`Results (1-50 of 410 matches)`) | search page | ran |
| labels scanned | 50 (page 1 of 9) | `bsd_sample.json` | ran |
| records fetched / in range | 50 / 50 | `bsd_sample.json` | ran |
| **records kept** | **25** (hard cap 25) | `bsd_sample.json` | ran |
| dropped by `analytic_rank <= 3` | 1 curve: `501029.a1`, `analytic_rank = 4` | API record | ran |
| selection rule | keep the two rank-3 curves, then earliest conductors, up to 25 | `finalize_sample.py` | ran |
| kept rank histogram | 0: 6, 1: 12, 2: 5, 3: 2 | `runs/bsd_coverage.json` | ran |
| kept conductor range | 500393 .. 505691 | `bsd_sample.json` | ran |
| output | `research/scratch/bsd_sample.json` | -- | -- |

The search-page URL's `analytic_rank=0-3` argument was silently ignored by the
page (the rank-4 curve appeared in its own results), so the filter was applied
locally after fetching each record.

### Task 3 - recompute and compare (900 s cap)

Estimator: `L_afe` + `tail_slope` imported from
`research/scratch/analytic_rank.py`, unchanged. Ladder points are the four
smallest of the smoke test - `1e-4, 3e-5, 1e-5, 1e-6` - and the fit is the same
`tail_slope(..., tail=4)`; the omitted larger-eps points cannot enter that fit.
`mp.dps = 30`. Conductor is supplied from LMFDB; the root number is detected
numerically by matching the AFE to the Dirichlet sum at `s = 2`, not taken from
the rank.

| item | value | source | checked |
|---|---|---|---|
| curves sampled | 25 | `bsd_sample.json` | ran |
| **curves evaluated** | **25** | `bsd_recompute.json` | ran |
| failures | 0 | `bsd_recompute.json` | ran |
| **mismatches (non-zero absolute difference)** | **0** | `bsd_recompute.json` | ran |
| first mismatch | **none** | -- | -- |
| root-parity inconsistencies (`w != (-1)^k`) | 0 | `bsd_recompute.json` | ran |
| point-count self-test mismatches (`p <= 50`) | 0 on all 25 | `bsd_recompute.json` | ran |
| total runtime | **1190.6 s** | `runs/bsd_coverage.json` | ran |
| run 1 | 17 curves, 863.0 s, stopped by the 900 s ceiling | stdout | ran |
| run 2 | remaining 8 curves, 327.6 s | stdout | ran |
| per-curve wall time | 20.7 s .. 112.6 s (load-dependent) | `runs/bsd_coverage.json` | ran |
| observed slopes, k = 0/1/2/3 | -0.000080 / 0.999920 / 1.999920 / 2.999920 | `bsd_recompute.json` | ran |

| label | N | LMFDB k | k_est | tail slope | w | wall |
|---|---:|---:|---:|---:|---:|---:|
| 505061.a1 | 505061 | 3 | 3 | 2.999920 | -1 | 20.7 s |
| 505691.a1 | 505691 | 3 | 3 | 2.999920 | -1 | 39.5 s |
| 500393.a1 | 500393 | 2 | 2 | 1.999920 | +1 | 33.4 s |
| 500719.a1 | 500719 | 0 | 0 | -0.000080 | +1 | 38.7 s |
| 500921.a1 | 500921 | 1 | 1 | 0.999920 | -1 | 36.5 s |
| 500921.b1 | 500921 | 1 | 1 | 0.999920 | -1 | 36.2 s |
| 501121.a1 | 501121 | 2 | 2 | 1.999920 | +1 | 47.0 s |
| 501139.a1 | 501139 | 1 | 1 | 0.999920 | -1 | 45.2 s |
| 501139.b1 | 501139 | 0 | 0 | -0.000080 | +1 | 44.7 s |
| 501317.a1 | 501317 | 0 | 0 | -0.000080 | +1 | 23.0 s |
| 501341.a1 | 501341 | 1 | 1 | 0.999920 | -1 | 47.9 s |
| 501367.a1 | 501367 | 1 | 1 | 0.999920 | -1 | 61.8 s |
| 501451.a1 | 501451 | 0 | 0 | -0.000080 | +1 | 61.4 s |
| 501563.a1 | 501563 | 2 | 2 | 1.999920 | +1 | 73.8 s |
| 502237.a1 | 502237 | 2 | 2 | 1.999920 | +1 | 82.5 s |
| 502237.b1 | 502237 | 1 | 1 | 0.999920 | -1 | 58.1 s |
| 502237.c1 | 502237 | 0 | 0 | -0.000080 | +1 | 112.6 s |
| 502507.a1 | 502507 | 1 | 1 | 0.999920 | -1 | 57.4 s |
| 502543.a1 | 502543 | 1 | 1 | 0.999920 | -1 | 38.4 s |
| 502543.b1 | 502543 | 1 | 1 | 0.999920 | -1 | 41.9 s |
| 502819.a1 | 502819 | 1 | 1 | 0.999920 | -1 | 36.4 s |
| 502829.a1 | 502829 | 1 | 1 | 0.999920 | -1 | 22.4 s |
| 503621.a1 | 503621 | 1 | 1 | 0.999920 | -1 | 36.1 s |
| 503803.a1 | 503803 | 2 | 2 | 1.999920 | +1 | 53.2 s |
| 503857.a1 | 503857 | 0 | 0 | -0.000080 | +1 | 41.6 s |

All 25 slopes sit at `k - 8.0e-5`, the expected `O(eps)` drift of a log-log
fit whose smallest point is `1e-6`; no slope is within rounding distance of a
half-integer, so the `round()` is not doing hidden work. The two rank-3 curves
are the load-bearing cases: at conductor ~5x10^5 the fit still lands on
2.999920, the same 8e-5 offset as every other rank, rather than drifting toward
2 or 4.

## Reproduced locally

```powershell
PS C:\Users\daive\Downloads\frontier-attack> Get-Command gp,sage,mwrank | Format-List Name,CommandType,Definition
Name        : gp
CommandType : Alias
Definition  : Get-ItemProperty
# sage, mwrank: not found

PS> python -m pip install --user cypari2
...
  meson.build:1:0: ERROR: Unknown compiler(s): [['icl'], ['cl'], ['cc'], ['gcc'], ['clang'], ['clang-cl'], ['pgcc']]
  WARNING: Failed to activate VS environment: Could not find ...vswhere.exe
error: metadata-generation-failed
# EXIT=1, elapsed 24.2 s
```

LMFDB's documented range URL, and the fallback that worked:

```text
GET .../api/ec_curvedata/?_format=json&_fields=...&conductor=500000-520000&analytic_rank=0-3
-> HTTP 404  LMFDB Page Not Found

GET .../EllipticCurve/Q/?conductor=500000-520000
-> "Results (1-50 of 410 matches)"
```

Sample construction (records already collected into `lmfdb_raw_result.json` by
an in-page script; this step is the filter and cap):

```powershell
PS> python research/scratch/finalize_sample.py
wrote ...\research\scratch\bsd_sample.json records 25
rank histogram: {0: 6, 1: 12, 2: 5, 3: 2}
conductors: 505061 .. 503857
# the second line is first/last record after the rank-3-first sort,
# not min/max; the kept range is 500393 .. 505691
```

Task 3, pass 1 (stopped by the 900 s tool timeout at 17/25):

```powershell
PS> python research/scratch/bsd_recompute.py
[ 1/25] 505061.a1  N=505061 k=3 k_est=3 slope=2.999920 w=-1 20.7s
[ 2/25] 505691.a1  N=505691 k=3 k_est=3 slope=2.999920 w=-1 39.5s
...
[17/25] 502237.c1  N=502237 k=0 k_est=0 slope=-0.000080 w=+1 112.6s
# killed at 900 s; no summary file written
```

Task 3, pass 2 (resume; the 17 rows above were transcribed from stdout, the
other 8 recomputed):

```powershell
PS> python research/scratch/bsd_finish.py
seeded 17 from run 1; evaluating 8 more
[ 1/8] 502507.a1  k_est=1 slope=0.999920 57.4s
...
[ 8/8] 503857.a1  k_est=0 slope=-0.000080 41.6s

evaluated=25/25 mismatches=0 parity_bad=0 run2=327.6s total=1190.6s
written ...\research\scratch\bsd_recompute.json
```

Machine-readable output: `runs/bsd_coverage.json`
(`{"task1": ..., "task2": ..., "task3": ...}`).

## What I could not find

- **No `gp` output at all.** The requested `E=ellinit([0,0,1,-1,0]); ellrank(E);
  ellanalyticrank(E)` for 37a1 was not executed; there is no PARI binary and no
  admissible way to install one without a system-wide installer run.
- **No cypari2.** The blocker is a missing C toolchain (and PARI headers), not a
  missing Python package; `pip` cannot build it in this environment.
- **No working LMFDB range API.** The exact URL in the task returns 404; the
  deployed API accepts only exact-equality filters, and `_offset` stops working
  before 20000, so the table cannot be paged to conductor 500000.
- **No random sample over the whole interval.** Only page 1 (50 of 410 labels)
  was scanned; kept conductors are 500393..505691, not spread to 520000.
- **No certified `k`.** `k_est` has no interval bound; the estimator's own
  docstring calls it a smoke test, and `mp.dps = 30` gives no error certificate.
  A mismatch would have been a lead; zero mismatches is a coverage statement.
- **No algebraic rank.** No `ellrank`, no `r1 = r2`, no `mwrank`; Task 3 says
  nothing about BSD beyond "the numerical analytic order agrees with LMFDB on
  these 25 curves".
- **No independent check of LMFDB's `analytic_rank` values.** They are used as
  the comparison target, i.e. taken on faith as a curated database field.
