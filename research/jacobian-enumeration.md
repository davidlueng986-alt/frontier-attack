# Jacobian pairs of degree 3–5: what is actually enumerated

## Answer

Nothing in the literature or in public code exhibits a non-triangular degree-3–5 **planar** Jacobian counterexample, and none can exist: the plane case is proved for total degree ≤104, and below 125 apart from (72,108) per the unrefereed GGHV preprint (Moh 1983; Nguyen 2025; GGHV 2022), and every Keller pair of total degree ≤4 is explicitly classified (Starkov 2016, 2017). The July 2026 Alpöge/Fable object is not a planar pair at all — it is a map ℂ³→ℂ³ of total degree 7 (component degrees 7, 6, 4; geometric degree 3) with det JF ≡ −2 and a three-point collision. The only genuinely new non-triangular small-degree families in print are Starkov's quartic families Φ_{a,b,d} and F_{c,s}; total degree 5 has no complete classification. Normalisation theory says a planar counterexample needs exactly two points at infinity per component, similar trapezoidal Newton polygons with top vertex (m,n), n>m>1, no edge parallel to the diagonal, and gcd(deg f, deg g) ≥ 16. Software-wise the only published Newton-polygon/degree enumerator is the Guccione–Guccione–Horruitiner–Valqui algorithm (arXiv:1708.07936), whose C++ implementation is described but not released, and no Sage/Magma/Singular/Maple package is dedicated to Jacobian pairs.

## Objects and numbers

| item | value | source |
|---|---|---|
| Alpöge–Fable map (dimension 3) | F = ((1+xy)³z + y²(1+xy)(4+3xy), y + 3x(1+xy)²z + 3xy²(4+3xy), 2x − 3x²y − x³z) | https://x.com/__alpoge__/status/2079028340955197566 |
| degree data | component degrees (7,6,4); total degree 7; generic fibre 3 points; det JF ≡ −2 | https://terrytao.wordpress.com/2026/07/21/a-digestion-of-the-jacobian-conjecture-counterexample/ |
| collision 1 | (0,0,−1/4), (1,−3/2,13/2), (−1,3/2,13/2) ↦ (−1/4,0,0) | https://sbseminar.wordpress.com/2026/07/20/the-new-counterexample-to-the-jacobian-conjecture/ |
| collision 2 | (1,−2,9), (−1/3,4,27), (−2/3,−1/2,−9/8) ↦ (−1,1,−1); hand-checked here | https://openproblem.space/wiki/Statement:alpoge-counterexample (unrefereed Fong–Fable ledger) |
| geometry | geometric degree 3; monodromy S₃; Im F = ℂ³∖M, M = {12a=b², 3bc=4} ≅ ℂ* | https://zenodo.org/records/21461572 |
| announcement | X post 19 Jul 2026 (02:19 UTC 20 Jul); no peer-reviewed article found, only tweet + expository notes | https://x.com/__alpoge__/status/2079028340955197566 |
| Wang 1980 | Keller maps of degree ≤2 are invertible, every dimension | https://doi.org/10.1016/0021-8693(80)90233-1 |
| total degree ≤3 | complete classification up to linear conjugacy: (x+α₂s²+α₃s³, y−α₂s²−α₃s³), s = x+y, det 1 | https://doi.org/10.15393/j3.art.2016.3510 |
| total degree ≤4 | complete classification; new quartic families Φ_{a,b,d} and F_{c,s} ≈ Φ_{1,1,1} | https://doi.org/10.15393/j3.art.2017.3870 |
| F_{c,s} (non-triangular, tame) | (x + cy², y − sx² − 2csxy² − c²sy⁴), c,s ≠ 0, det 1; inverse y = v+su², x = u−cy² | https://doi.org/10.15393/j3.art.2017.3870 |
| geometric-degree walls | no 3-sheeted Keller map (Orevkov 1987); no 4-sheeted with one dicritical component (Domrina–Orevkov 1998); topological degree ≤5 excluded (Żołądek 2008) | https://doi.org/10.1070/IM1987v029n03ABEH000984; https://doi.org/10.1007/BF02313031; https://doi.org/10.1016/j.top.2008.04.001 |
| total-degree walls | ≤100 (Moh 1983); ≤104 (Nguyen 2025); <125 except (72,108) (GGHV 2022, unrefereed preprint) | https://doi.org/10.1515/crll.1983.340.140; https://arxiv.org/abs/1902.05923; https://arxiv.org/abs/2204.14178 |
| gcd constraint | gcd(deg P, deg Q) ≥ 16; ≠ 2p for prime p | https://doi.org/10.1016/0022-4049(90)90005-3; https://arxiv.org/abs/1401.1784 |
| Newton polygons | N(f) and N(g) of a Jacobian pair are similar | https://doi.org/10.1016/0022-4049(91)90128-O |
| Abhyankar constraint | each component has ≤2 points at infinity; one point ⇒ invertible; JC ⇔ every component's Newton polygon is a triangle | https://mathweb.tifr.res.in/Documents/Publications/Lectures/tifr57.pdf; https://en.wikipedia.org/wiki/Jacobian_conjecture |
| minimal-counterexample shape | N(f) in a trapezoid with top vertex (m,n), n>m>1, edges ∥ y-axis and diagonal; no edge ∥ bisector (Cassou-Noguès) | https://arxiv.org/abs/2106.06869 |
| reduced Jacobian pair | P = xⁿ + P_{n−1}(y)x^{n−1}+…, Q = x^m + Q_{m−2}(y)x^{m−2}+…, J(P,Q)=1 | https://doi.org/10.1007/BF02808206 |
| Drużkowski | cubic-linear reduction proved 1983; A²=0 refinement 2001; symmetric cubic-linear case tame 2005; no explicit non-invertible Drużkowski map found | https://doi.org/10.1007/BF01459126; https://ruj.uj.edu.pl/xmlui/handle/item/80085; https://doi.org/10.4064/ap87-0-7 |
| GGHV enumerator | computes all admissible Newton-polygon corners, gcd ≤35 and max degree ≤150; C++/PostgreSQL/D3 pipeline described, no code URL given | https://arxiv.org/abs/1708.07936 |
| Windows support | Sage: no Windows binaries (WSL); Singular: Cygwin or WSL; Magma: Windows unsupported per FAQ; Maple: native Windows | https://www.sagemath.org/download-windows.html; https://singular.mathematik.uni-kl.de/index.php/singular-download/install-windows.html; https://magma-maths.org/faq/; https://maplesoft.com/products/system_requirements.aspx |
| public code (verification, not enumeration) | SymPy/Lean audit repos; Wolfram Language checker; Windows-native Python parts | https://github.com/wstrinz/plane-jacobian-72-108; https://github.com/Kakarottoooo/jacobian-2d-research; https://github.com/alok/jacobian-two; https://github.com/jzkay12/jacobian_conjecture |
| this repo's degree-3 search | 2,000 sampled maps, 3 constant-det, 3/3 triangular, 0 counterexamples | https://github.com/davidlueng986-alt/frontier-attack (findings/2026-11-02-jacobian-cubic.md) |

## What I could not find

- No explicit degree-3, -4 or -5 **planar counterexample** — provably none in that range (Moh/Nguyen; Starkov).
- No complete classification or coefficient-level table of total-degree-5 planar Keller pairs; Starkov stops at 4.
- No released source code for the GGHV corner/Newton-polygon enumerator (paper says "website under development").
- No explicit non-invertible Drużkowski / cubic-linear Keller map in print.
- No peer-reviewed article for the Alpöge map found; only the X post, blog expositions and unrefereed verification notes.
- No dedicated Jacobian-pair enumeration package in Sage, Magma, Singular or Maple.
