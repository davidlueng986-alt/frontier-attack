# Research corpus

Permanent, machine-readable source base for the three problems in this repo:

| tag | problem | terminal object |
|---|---|---|
| `hn` | Hadwiger–Nelson | explicit finite unit-distance graph with χ = 7 + DRAT refutation of 6-colouring |
| `jac` | planar Jacobian | `(P,Q)`, `J ≡ 1`, degrees `(72,108)` (or a later pair), `p ≠ q` with `P(p)=P(q)`, `Q(p)=Q(q)` |
| `bsd` | Birch–Swinnerton-Dyer | `(E,r,k)` with both ranks proved and `r ≠ k` |
| `methods` | cross-cutting | techniques that transfer between the three |

**59 items** as of 2026-11-02. Every entry was resolved against a live URL; `verified` means the URL
(or the DOI resolver record for that exact DOI) was fetched and the title matched. No entry is
marked verified on the strength of a citation alone.

## Files

| file | what it is |
|---|---|
| `index.json` | the canonical list: one object per item, 12 fields, deduplicated by arXiv id / DOI / repo path |
| `hadwiger-nelson.md` | annotated bibliography, grouped by theme, with the object table |
| `jacobian.md` | annotated bibliography: degree ladder, Newton polygons, the 2026 C³ map, the (72,108) claim |
| `bsd.md` | annotated bibliography: GZK gate, verification landscape, the only admissible mismatch shapes |
| `solved-analogues.md` | problems of the same shape that are already settled, and the method that settled them |
| `open-analogues.md` | sibling problems that remain open, with the current bounds |
| `methods.md` | transferable techniques, each with a worked pointer into the corpus |
| `manifest.json` | the 20 downloaded PDFs: url, sha256, bytes, licence, local path |

## What is committed and what is not

The repo is public. **Do not commit the PDFs.** Most are not redistributable.

- **Committed** (this directory): `README.md`, `index.json`, the six `.md` notes,
  `manifest.json`, and links/identifiers. These are metadata and original prose.
- **Not committed**: `research/scratch/corpus/pdfs/*.pdf` and
  `research/scratch/corpus/txt/*.txt`. `research/scratch/` is already in `.gitignore`
  (line 6). The PDFs are 20 files, 18,155,732 bytes. Re-download with
  `research/scratch/corpus/download.py`; hashes in `manifest.json` let you check you got
  the same bytes.
- The extracted `.txt` files are working copies for `grep`; they are not part of the corpus.

Licences in `manifest.json` are conservative:

- `arXiv` = the arXiv non-exclusive distribution licence (no CC badge on the abs page).
- `CC-BY` = a Creative Commons Attribution 4.0 link was present on the abs page
  (`2303.14714`, `2204.14178`, `2502.17655`).
- `unknown` = no licence statement was found at the URL (the two OpenAI preprints, the
  TIFR lecture notes, and all publisher PDFs that were not downloaded).

`license` in `index.json` follows the same convention. Where the value is `unknown`,
link to it; do not redistribute it.

## index.json schema

```json
{"id": "arXiv:1804.02385",
 "title": "The chromatic number of the plane is at least 5",
 "authors": ["de Grey, Aubrey D. N. J."],
 "year": 2018,
 "url": "https://arxiv.org/abs/1804.02385",
 "kind": "paper",                       // paper | repo | dataset | note
 "problem": ["hn"],                     // hn | jac | bsd | methods
 "role": "core",                        // core | method | analogue-solved | analogue-open
 "why": "one sentence, specific",
 "facts": ["1581 vertices ... (abstract; Fig. 9)"],
 "license": "arXiv",                    // arXiv | CC-BY | MIT | Apache-2.0 | unknown
 "verified": true}
```

Rules used when building it:

- `facts` are numbers or theorem statements with a location in the source (`abstract`,
  `Theorem 1.1`, `Table 1`, `p. 3`), never adjectives.
- Every `core` item has at least one fact.
- One item per arXiv id / DOI / repository. The same work is never entered twice under
  a preprint and a journal DOI.
- `verified: true` requires an `https://` URL. There are no `http://` URLs in the file.
- `year` is the publication year of the version at `url` (arXiv submission year for
  preprints whose journal version is later).

## Validation

```bash
python tools/corpus_check.py
```

Exits 0 only if there are at least 40 items, every entry has all 12 fields, ids are
unique, every `core` item has a fact, and every verified URL is `https://`. It prints
counts per problem and per role.

Current output: 59 items — `hn` 18, `jac` 23, `bsd` 14, `methods` 4;
`core` 37, `method` 7, `analogue-solved` 10, `analogue-open` 5; 0 unverified.

## How to use it

1. Start with the problem note (`hadwiger-nelson.md`, `jacobian.md`, `bsd.md`). Each is
   organised by *what the item does for the terminal object*, not by author or date.
2. Pull the exact numbers from `index.json` rather than retyping them. `facts` is designed
   to be quotable in a plan or a claim file.
3. Before adding an item, check the id against `index.json`; the check rejects duplicates.
4. Before citing a PDF, confirm the hash in `manifest.json`. arXiv replaces PDFs; a
   changed hash means the local copy no longer matches what was read.

## What this corpus does not do

It does not contain the terminal objects. It contains the evidence needed to avoid
re-deriving dead ends: the 6993-vertex floor for 7-chromatic plane graphs, the GZK
exclusion of `k ≤ 1` for BSD, and the one-directional GGHV normalisation for the
Jacobian problem. Those three facts are the reason the three current plans look the way
they do.
