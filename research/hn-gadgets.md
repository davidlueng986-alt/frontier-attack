# Hadwiger–Nelson gadgets and small 5-chromatic unit-distance graphs

## Answer
The smallest published 5-chromatic unit-distance graph in the Euclidean plane is Parts' 509-vertex, 2,442-edge graph (arXiv:2010.12665), still the record in August 2026 (https://mathworld.wolfram.com/Hadwiger-NelsonProblem.html). De Grey's original construction has 20,425 vertices/151,311 edges, reduced to 1,581/7,877, with coordinates in §5.1 of arXiv:1804.02385. Heule's SAT trimming reached 553/2,722 (arXiv:1805.12181) and 529/2,670 (arXiv:1907.00929), all `.vtx` public (https://github.com/marijnheule/CNP-SAT). De Grey's "H" is a 7-vertex, 12-edge wheel used in L (121 vertices, 52 copies of H) to force a monochromatic √3-triangle in at least one copy—a 3-vertex pattern, not a two-terminal equality; equality is forced by the diamond (4/5), and spindling two diamonds gives the Moser spindle (7/11) (https://web.archive.org/web/2018id_/http://michaelnielsen.org/polymath1/index.php?title=Hadwiger-Nelson_problem). No 6-chromatic unit-distance graph in the plane is known (MathWorld, above); 6-chromatic examples are two-distance (16 vertices, arXiv:2010.12656) or in R³ (47–60 vertices, https://mathworld.wolfram.com/deGreyGraphs.html).

## Objects and numbers
| item | value | source |
|---|---|---|
| de Grey G (corrected) | 1,581 v; 7,877 e; χ=5; coordinates in §5.1 | arXiv:1804.02385; Polymath16 wiki |
| de Grey G₀ | 1,585 v; 7,909 e; DIMACS + Sage vertices | Polymath16 wiki |
| Heule intermediates | 874/4,461; 826/4,273; 803/4,144; 633/3,166; 610/3,000 (v/e) | Polymath16 wiki; CNP-SAT |
| Heule 553 (2018 paper) | 553 v; 2,722 e (Heule 2019 abstract: 2,720) | arXiv:1805.12181; arXiv:1907.00929 |
| Heule 529 (2019) | 529 v; 2,670 e | arXiv:1907.00929 |
| Heule 517 | 517 v; 2,579 e | https://mathworld.wolfram.com/HeuleGraphs.html |
| 510-vertex pair | Heule 510/2,504; Parts 510/2,502 | MathWorld HeuleGraphs; arXiv:2010.12665 |
| Parts 509 (current record) | 509 v; 2,442 e; coordinates NOT in arXiv—Polymath16 comment 24942/Dropbox `JP/Large`; mirror `vtx/509_parts.vtx` | arXiv:2010.12665; https://github.com/vasnesterov/HadwigerNelson |
| Moser-spindle-free | 64,513 v (Voronov 2022); 1,441 v (Heule 2021); 2,131 v/12,530 e (Haugland 2026) | arXiv:2106.11824; arXiv:2608.04542; https://mathworld.wolfram.com/HauglandGraphs.html |
| Exoo–Ismailescu G₄₀ | 40 v; 82 e; forces an equal pair at distance 8/3 if no mono √(11/3)-edge | arXiv:1805.00157 |
| Exoo–Ismailescu G₇₉ | 79 v; 165 e; forces a mono √(11/3)-edge | arXiv:1805.00157 |
| Exoo–Ismailescu G₄₉ | 49 v; 180 e; mono √(11/3)-edge ⇒ mono 1/√3-triangle | arXiv:1805.00157 |
| Exoo–Ismailescu G₆₂₇ | 627 v; 2,982 e; kills that mono 1/√3-triangle | arXiv:1805.00157 |
| H (de Grey's wheel; "H-graph") | 7 v; 12 e; four essential 4-colourings, two with a mono √3-triangle; every 5-colouring has a mono √3- or 2-edge | Polymath16 wiki; arXiv:1804.02385 |
| L | 121 v; 301 e; 52 copies of H; every 4-colouring has an H with mono √3-triangle | Polymath16 wiki |
| M | 1,345 v; 8,268 e; no 4-colouring with mono triangle in central H | Polymath16 wiki |
| Moser spindle | 7 v; 11 e; χ=4; diamond tips forced equal in any 3-colouring | Polymath16 wiki; https://en.wikipedia.org/wiki/Moser_spindle |
| Diamond | 4 v; 5 e; tips at distance √3 forced equal in any 3-colouring | https://mathworld.wolfram.com/DiamondGraph.html; Wikipedia Moser spindle |
| Spindle operation | glue two copies at one vertex, rotate so a second corresponding pair is unit distance; if that pair is forced equal, result is not k-colourable | Polymath16 wiki |
| 6-chromatic two-distance (Parts) | 16 v | arXiv:2010.12656 |
| 6-chromatic two-distance (Exoo–Ismailescu) | 426 v; 2,009 unit + 892 length-2 edges | arXiv:1909.13177 |
| 6-chromatic unit-distance in R³ | 47/48 v (de Grey–Haugstrup); 59/60 v (de Grey) | https://mathworld.wolfram.com/deGrey-HaugstrupGraphs.html; MathWorld deGreyGraphs |
| 5-chromatic triangle-free in R³ | 61 v (de Grey 2026) | MathWorld deGreyGraphs |

### Code / data
| repo | what it does | source |
|---|---|---|
| CNP-SAT | Heule SAT trimming; `.vtx` for 510–874 and T721 | https://github.com/marijnheule/CNP-SAT |
| cnp-5 | Python/PySAT library; generates de Grey's N; UNSAT-core shrinking | https://github.com/oisincar/cnp-5 |
| HadwigerNelson | Lean 4 formalisation of the 510-vertex Heule graph; `vtx/509_parts.vtx` | https://github.com/vasnesterov/HadwigerNelson |
| wustep/maths 509 | independent checker: 509 v / 2,442 e + DRAT | https://github.com/wustep/maths/blob/main/problems/unit-distance-509/compute/README.md |
| haugland-2131-certificates | machine-checked DRAT/cake_lpr certificates; 1,501-v reduction | https://doi.org/10.5281/zenodo.22435777 |

Reference URLs: Polymath16 wiki = https://web.archive.org/web/2018id_/http://michaelnielsen.org/polymath1/index.php?title=Hadwiger-Nelson_problem · HeuleGraphs = https://mathworld.wolfram.com/HeuleGraphs.html · deGreyGraphs = https://mathworld.wolfram.com/deGreyGraphs.html · Moser spindle = https://en.wikipedia.org/wiki/Moser_spindle

## What I could not find
- No published source calls de Grey's H the "H-graph"; the term collides with the established H-graph class (arXiv:1706.00575). No two-terminal H gadget: H forces a 3-vertex monochromatic pattern.
- No single composed vertex/edge count in arXiv:1805.00157; only the G₄₀/G₇₉/G₄₉/G₆₂₇ sizes.
- MathWorld (Aug 2026) lists no 2024–2026 plane reduction below 509: 2026 has 2,131 (spindle-free) and 61 in R³ (triangle-free); the 1,501-v certificate is above Heule's 1,441.
- 509 coordinates are absent from arXiv:2010.12665 (only Polymath16 comment 24942/Dropbox and mirrors); Mixon's 1,577 edge count not found.
