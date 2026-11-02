# Gaps: what was searched for and not verified

This file records the negative search results. It is deliberately specific: the query, the
tool, the outcome, and what was added instead (if anything). A failed search is a fact about
the literature, not a defect of the corpus, and it should stop the next agent from repeating
the query.

## How the searches were run

| tool | use |
|---|---|
| `https://api.openalex.org/works?search=...` | title/author discovery, then exact DOI lookup |
| `https://api.openalex.org/works/doi:<DOI>` | bibliographic verification (title, year, authors, abstract) |
| `https://arxiv.org/abs/<id>` | title, authors, abstract, submission date |
| direct `urllib` fetch | repository pages, PDFs, MathWorld, LMFDB, blogs |
| Exa web search | repository/thesis URLs that OpenAlex does not index |

`verified: true` in `index.json` means one of the first four returned a matching title. It
does **not** mean the paper's theorem was reproved here.

## Searches that failed

| # | query | hoped for | outcome | used instead |
|---:|---|---|---|---|
| 1 | `"The Hadwiger-Nelson problem" Coulson` (OpenAlex) | Coulson's survey of the problem | returned unrelated planar-graph colouring papers; no survey record | `doi:10.1016/s0012-365x(01)00183-2` (Coulson 2002, the 15-colouring of R³) |
| 2 | `Exoo unit distance graphs table chromatic number` (Exa) | Exoo's own table page | returned the 2004 epsilon-unit-distance paper, not a table | `doi:10.1007/s00454-004-1092-8` |
| 3 | `Polymath16 first thread` (Exa) | a single index of threads 1–20 | found thread 1 and the wiki; no consolidated thread index | `polymath16:thread1` plus the existing `polymath16` wiki item |
| 4 | `Keller's conjecture in positive characteristic` (OpenAlex) | a primary counterexample paper | returned Lang 2024 and unrelated papers; no primary counterexample verified | `doi:10.1016/j.jpaa.2024.107797` (Lang survey) |
| 5 | `Selmer group ranks survey` (OpenAlex) | a dedicated Selmer-rank survey | no clean survey record | `doi:10.4007/annals.2010.172.567` (parity) + existing Bhargava–Shankar |
| 6 | LMFDB reliability page, searched for a cited zero-counting paper | the exact reference LMFDB uses | the page describes the method but cites no single paper | `doi:10.1080/10586458.2006.10128976` (Booker) and `platt:thesis`; **not** claimed to be LMFDB's exact implementation |
| 7 | DOI `10.7551/mitpress/4316.001.0001` for Matiyasevich's book | the book record | resolved to an unrelated "Logic Programming" record | `doi:10.1080/00029890.1973.11993265` (Davis 1973) |
| 8 | arXiv `1704.05530` for cake_lpr | the cake_lpr paper | resolved to an unrelated martingales paper | `doi:10.1007/978-3-319-63046-5_14` |
| 9 | direct `repository.ubn.ru.nl/bitstream/.../75148.pdf` | de Bondt's thesis PDF | HTTP 403 from the bitstream; the handle page is 200 and gives the title | `https://repository.ubn.ru.nl/handle/2066/75148`; PDF fetched with a Referer header |
| 10 | `ams.org/notices/200811/tx081101382p.pdf` and the 1976 Bull AMS PDF | local copies of Gonthier and Appel–Haken | HTTP 403 from ams.org in this environment | DOI metadata stands; both are **not** in `manifest.json` |
| 11 | `"Matiyasevich Hilbert's tenth problem"` (OpenAlex) | the MIT Press book | only reviews and the Davis 1973 article were returned | Davis 1973 |
| 12 | `"Symmetry breaking for SAT"` (OpenAlex) | a dedicated SAT symmetry paper | first query rate-limited; the exact Shatter DOI was found on retry | `doi:10.1109/dac.2003.1219135` |
| 13 | `"Newton Puiseux algorithm"` (OpenAlex) | a dedicated algorithm paper | rate-limited; no clean paper record | `doi:10.1007/978-3-662-05355-3` (Basu–Pollack–Roy textbook) |
| 14 | `"SAT encoding graph coloring"` (OpenAlex) | a dedicated encoding paper | rate-limited; no clean paper record | `doi:10.3233/FAIA336` (Handbook of Satisfiability) |

## Items added with weaker-than-ideal evidence

These are in `index.json` with `verified: true` (the item was fetched), but the *fact* is
from a secondary source or a standard attribution, not from the fetched text. Treat the
fact as a pointer, not a quotation.

| id | fact | why it is weaker |
|---|---|---|
| `doi:10.1090/s0002-9904-1976-14122-5` | "1,936 reducible configurations" | the 1976 abstract states only the four-colour theorem; the count is from the book/secondary literature |
| `doi:10.1006/jctb.1997.1750` | "633 configurations" | the DOI page has no abstract; the count is the standard attribution, not read from the paper here |
| `doi:10.1007/bfb0062705` | "average rank over quadratic fields" | only the title was resolved; no abstract on the DOI record |
| `doi:10.1007/bf02571929` | "det J > 0 everywhere, not injective" | the theorem is standard, but the paper's abstract was not available; no degree was recorded |
| `doi:10.1016/j.dam.2012.07.018` | "canonical augmentation" | no abstract on the DOI record; the phrase is from the paper's known content |
| `doi:10.1016/s0012-365x(00)00406-4` | "chi(R^3) >= 6" | no abstract; the construction size (60 vertices) was **not** verified and is not recorded |
| `doi:10.1007/s00454-004-1092-8` | "chi(G_epsilon) = 7 for certain epsilon" | the abstract says "certain values" but does not list them; no numeric epsilon is recorded |
| `doi:10.1007/978-3-662-05355-3` | "Newton–Puiseux and quantifier elimination" | book contents, not a theorem read from the text |
| `doi:10.1007/978-3-319-16721-3` | "Buchberger algorithm" | book contents, not a theorem read from the text |
| `stein:bsd-book` | year 2005 | the PDF preface is dated 2005; the copyright page/imprint was not checked |

## Things deliberately not added

- **A "Coulson survey"**: no such record was found. The 2002 paper is added instead; a
  survey may exist in a book chapter without a DOI.
- **Exoo's table page**: no stable URL was verified. The 2004 paper is the citable object.
- **Polymath16 threads 2–20**: they exist as blog posts, but no single index was fetched.
  Thread 1 is in the corpus; the wiki is the consolidated record.
- **A primary positive-characteristic counterexample**: Lang's 2024 survey is in; the
  primary construction was not identified.
- **A primary "JC fails in char p" theorem**: the Mondello 2026 arXiv paper is already in
  the corpus; the older char-p literature is represented only by the survey.
- **Matiyasevich's book**: the DOI could not be verified. Davis 1973 carries the theorem.
- **A local copy of Gonthier/Appel–Haken**: ams.org returns 403 here. The bibliographic
  records are verified; the PDFs are not in the manifest.

## Next searches (queued, not done)

1. `site:ams.org "Formal Proof - The Four-Color Theorem"` for a mirror of the Notices PDF.
2. `Coulson "Hadwiger-Nelson" survey` in zbMATH/MathSciNet, not OpenAlex.
3. `Keller map positive characteristic counterexample` in zbMATH, filtering pre-2000.
4. `Selmer group rank distribution survey` in MathSciNet.
5. `LMFDB zero verification Booker Platt` in the LMFDB source repository, not the
   reliability page.
6. `Polymath16 thread 20` for the final thread and any consolidated index.
7. `Appel Haken 1977 Illinois J. Math Part I Part II` for the full papers (the 1976 Bull
   AMS announcement is not the full proof).
8. `Tate algorithm Kodaira types table` for a page-level reference rather than the title.
