# HN-509: independent verification of the 509-vertex 5-chromatic unit-distance graph

## Answer

The record checks out. `vtx/509_parts.vtx` from `vasnesterov/HadwigerNelson` (master, commit
`ad7843546da5f4db4fcf9827e0951f34295ad7cc`) is a **509-vertex unit-distance graph with 2442 unit-distance
edges**, connected, minimum degree 4. Its **4-colouring instance is UNSAT** and it **is 5-colourable**, so this
graph has chromatic number exactly 5. Since every edge is exactly unit length, any proper 4-colouring of the
plane would induce a proper 4-colouring of this graph; therefore **chi(R^2) >= 5**.

Three independent solvers return UNSAT (CaDiCaL 1.9.5 133.2 s, CaDiCaL 1.5.3 122.2 s, Glucose 4.2.1 103.2 s,
all inside a 600 s cap), and a structurally different binary encoding gives the same verdict. The file
contains **coordinates only** —there is no `p edge` header and no edge list; the 2442 edges are recovered from
exact unit-distance relations, and every one of them was proved symbolically to have squared length exactly 1.
The derived count matches the published record (Jaan Parts, arXiv:2010.12665: 509 vertices, 2442 edges).

## Objects and numbers

| object | value | how obtained | checked |
|---|---|---|---|
| vertices `N` | **509** | 509 non-empty lines, each a `{x, y}` pair; no header, no comments, no DIMACS records | yes |
| edges `M` | **2442** | all pairs with `abs(d^2 - 1) < 1e-50` at 100-digit precision; nearest non-edge has `abs(d^2-1) = 1.602e-4`, so the cut is unambiguous | yes |
| published value of `M` | **2442** | Jaan Parts, *Graph minimization...*, arXiv:2010.12665, abstract and Table 1 ("509 vertices and 2442 edges"); also MathWorld "Parts Graphs" table (509 / 2442) | matches derived `M` |
| coordinates present? | **yes, exact algebraic** | Mathematica syntax: rationals and `Sqrt[...]`, incl. nested radicals `Sqrt[(5*(7 ± Sqrt[33]))/2]` at 4 vertices | yes |
| explicit edge list present? | **no** | only `{x, y}` lines; graph is implied by geometry (repo's `ParseVtx.lean`/`BuildEdge.lean` do the same) | yes |
| max edge deviation `abs(d^2 - 1)` | **1.457e-99** | 100-digit `mpmath` arithmetic over all 2442 edges (numerically zero) | yes |
| exact edge proof | **2442 / 2442 edges have `d^2 - 1` simplifying to exactly 0** | `sympy.simplify` on exact radical expressions, 11.1 s, 0 nonzero | yes |
| min `abs(d^2 - 1)` over non-edges | **1.602e-4** | 100-digit arithmetic | yes |
| connected components | **1** (all 509 vertices) | BFS on derived edge list | yes |
| degree min / max / mean | **4 / 36 / 9.597** | derived edge list (deg-4 vertices: 6) | yes |
| duplicate / coincident vertices | **0**; closest distinct pair 0.021130 | 100-digit arithmetic | yes |
| 5-colouring | **SAT**, 5 colours used, 0 monochromatic edges | CaDiCaL 1.9.5, 0.011 s | yes |
| 4-colouring, one-hot CNF | **UNSAT**, 2036 vars, 13331 clauses | CaDiCaL 1.9.5, **133.208 s** | yes |
| 4-colouring, one-hot CNF (2nd solver) | **UNSAT** | CaDiCaL 1.5.3, **122.245 s** | yes |
| 4-colouring, one-hot CNF (3rd solver) | **UNSAT** | Glucose 4.2.1 (`glucose4`), **103.180 s** | yes |
| 4-colouring, binary CNF | **UNSAT**, 1018 vars, 9768 clauses | CaDiCaL 1.9.5, **132.400 s** | yes |
| encoding self-test | K4: UNSAT at k=1,2,3 and SAT at k=4; C5: UNSAT at k=1,2, SAT at k=3; K3 likewise | same `build()` code used for the 509 instance | yes |
| negative control, exact check | perturbing one coordinate by 1/1000 makes `simplify(d^2 - 1) != 0` (3/3) | guards against a vacuous exact check | yes |
| provenance | blob `854a7df2d17d6f3e0aeaab53db0dd6eefa7935c5`, 25,515 bytes, sha256 `770a585a6c1e1222355322707479cb826e9ada560279da904ef89c15c99ff0b5` | `git hash-object` of downloaded file equals GitHub blob SHA (byte-identical) | yes |

Interpretation of the SAT run: variables `x[v][c]` ("vertex v has colour c", one-hot, 4x509 = 2036 vars);
clauses are one at-least-one per vertex (509 x size 4), six at-most-one pairs per vertex (3054 binary), and
`not(x[u][c] and x[v][c])` for every edge and colour (2442 x 4 = 9768). Total 13,331 clauses. UNSAT means no
proper 4-colouring exists. The binary encoding uses 2 bits per vertex (1018 vars) and no one-hot clauses at
all —only the 4 clauses per edge forbidding equal 2-bit patterns —so it is an independent CNF for the same
combinatorial claim.

Environment: Windows, Python 3.14.2, `python-sat 1.9.dev17` (bundled CaDiCaL 1.9.5 = `cadical195` and CaDiCaL
1.5.3 = `cadical153`), `sympy 1.14.0`, `mpmath 1.3.0`.

## Reproduced locally

Working directory `C:\Users\daive\Downloads\frontier-attack`; all scratch files in `research/scratch/`.

### 1. Fetch (the `main` URL is 404; `master` is the default branch)

```
PS> Invoke-WebRequest https://raw.githubusercontent.com/vasnesterov/HadwigerNelson/main/vtx/509_parts.vtx
    -> 404
PS> Invoke-WebRequest https://raw.githubusercontent.com/vasnesterov/HadwigerNelson/master/vtx/509_parts.vtx -OutFile research\scratch\509_parts.vtx
    -> OK
PS> git hash-object research\scratch\509_parts.vtx
854a7df2d17d6f3e0aeaab53db0dd6eefa7935c5        # == GitHub blob SHA for vtx/509_parts.vtx
PS> (Get-FileHash research\scratch\509_parts.vtx -Algorithm SHA256).Hash
770A585A6C1E1222355322707479CB826E9ADA560279DA904EF89C15C99FF0B5
```

`gh api repos/vasnesterov/HadwigerNelson/contents/vtx` lists only `509_parts.vtx`, `510_heule.vtx`,
`moser.vtx`, `triangle.vtx`. The repo README states that formalizing the 509-vertex Parts graph is still TODO
(the `build_edge` tactic cannot handle `Sqrt[(5*(7 + Sqrt[33]))/2]`), so no edge list or machine-checked proof
for this graph ships with the repo.

### 2. Parse and build the graph (`research/scratch/parse_vtx.py`)

```
PS> python research\scratch\parse_vtx.py research\scratch\509_parts.vtx
{
  "n": 509,
  "m": 2442,
  "max_dev_edge": 1.4573089238849006e-99,
  "min_nonedge_gap": 0.00016021784676687766,
  "deg_min": 4,
  "deg_max": 36,
  "deg_hist": {"4": 6, "5": 42, "6": 56, "7": 62, "8": 45, "9": 40, "10": 69, "11": 50,
               "12": 61, "13": 28, "14": 15, "15": 14, "16": 6, "17": 4, "18": 4, "22": 6, "36": 1},
  "components": [509],
  "isolated": []
}
```

Parsing is exact: each coordinate is `sympy.sympify`ed (`Sqrt[x]` -> `sqrt(x)`), then evaluated at 100 decimal
digits. (A first naive parser used Python `1/6` float division and produced only 1654 edges with 94 isolated
vertices; that artefact is why the final script uses sympy rationals. The clean result above has 2442 edges and
no isolated vertices.)

### 3. Exact edge-length verification (`research/scratch/exact_check.py --all`)

```
PS> python research\scratch\exact_check.py --all
  checked 2420/2442 in 11.1s
  checked 2440/2442 in 11.1s
checked 2442 edges exactly in 11.1s; nonzero: 0
```

Negative control: for the first three edges, adding `1/1000` to one x-coordinate and re-running
`sympy.simplify(d^2 - 1)` returned nonzero each time, so the symbolic test is not vacuous
(`research/scratch/summary_stats.py`):

```
"negative_control_perturbed_edges_nonzero": [true, true, true]
```

### 4. CNF export (`research/scratch/write_cnf.py`)

```
research/scratch/509_4coloring.cnf 2036 13331     (187,153 bytes)
research/scratch/509_5coloring.cnf 2545 17809     (255,067 bytes)
```

### 5. SAT verdicts

```
PS> python research\scratch\solve_color.py --colors 4 --solver cadical195 --timeout 600
n=509 m=2442 k=4 vars=2036 clauses=13331 solver=cadical195
{"verdict": "UNSAT", "solver": "cadical195", "seconds": 133.208, "timeout": 600.0,
 "vars": 2036, "clauses": 13331}

PS> python research\scratch\solve_color.py --colors 4 --solver cadical153 --timeout 600
n=509 m=2442 k=4 vars=2036 clauses=13331 solver=cadical153
{"verdict": "UNSAT", "solver": "cadical153", "seconds": 122.245, "timeout": 600.0,
 "vars": 2036, "clauses": 13331}

PS> python research\scratch\solve_color.py --colors 4 --solver glucose4 --timeout 600
n=509 m=2442 k=4 vars=2036 clauses=13331 solver=glucose4
{"verdict": "UNSAT", "solver": "glucose4", "seconds": 103.18, "timeout": 600.0,
 "vars": 2036, "clauses": 13331}

PS> python research\scratch\binary_solve.py
binary clause semantics self-test: OK
n=509 m=2442 binary k=4 vars=1018 clauses=9768
{"verdict": "UNSAT", "solver": "cadical195", "seconds": 132.4,
 "vars": 1018, "clauses": 9768}

PS> python research\scratch\solve_color.py --colors 5 --solver cadical195 --timeout 600
n=509 m=2442 k=5 vars=2545 clauses=17809 solver=cadical195
{"verdict": "SAT", "solver": "cadical195", "seconds": 0.011, ...}
colours used: 5 bad edges: 0
```

The 600 s cap was enforced by a timer calling `solver.interrupt()`; all 4-colour runs finished well inside it,
so the verdicts are not timeout artefacts. Encoding self-test (`research/scratch/encoding_selftest.py`, same
`build()` function):

```
K4: k=1 UNSAT OK   k=2 UNSAT OK   k=3 UNSAT OK   k=4 SAT OK
C5: k=1 UNSAT OK   k=2 UNSAT OK   k=3 SAT OK
K3: k=1 UNSAT OK   k=2 UNSAT OK   k=3 SAT OK
```

## What I could not find

- **No explicit edge list or DIMACS graph header in the file.** `509_parts.vtx` is a bare list of 509
  Mathematica-syntax coordinate pairs; the 2442 edges are *derived* (all unit-distance pairs), not stated.
  Parts' paper (arXiv:2010.12665) also defines *strict* unit-distance graphs —every pair at distance exactly 1
  is an edge —so the derivation matches the published convention, and the derived `M = 2442` matches Table 1.
- **No machine-checked proof of non-4-colourability for the 509-vertex graph in the repository.** The README
  says the smallest known Parts graph is "currently beyond the capabilities of `build_edge`" and lists its
  formalization as TODO; only the 510-vertex Heule graph is formalized there.
- **No UNSAT certificate (DRAT/LRAT) was produced or checked.** The verdict rests on two CaDiCaL versions
  (1.9.5, 1.5.3) plus a structurally different binary encoding, all agreeing, and on encoding self-tests.
  A DRAT proof would make the result independently machine-checkable; `python-sat`'s CaDiCaL wrappers do not
  expose proof logging here.
- **One solver wrapper could not be used as an extra cross-check.** The `kissat404` wrapper in this
  `python-sat` build crashes the interpreter (no output, exit 1) on even a 2-clause CNF; it was dropped.
  (Glucose 4.2.1 did finish —UNSAT in 103.18 s —and is included above.)
- **No claim about minimality.** I verified that *this* graph is 5-chromatic; I did not verify that 509 is the
  smallest possible, nor re-derive the graph. Parts (2020) likewise does not claim vertex-minimality —the
  paper's method returns a minimum among the subgraphs it searches —and the file is not checked against the
  paper's own vertex list beyond the matching vertex/edge counts.
