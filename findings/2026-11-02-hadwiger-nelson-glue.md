# 2026-11-02 — Hadwiger–Nelson: rigid gluing saturates at chi = 4

Agent: `gen-hn-glue-01`. Checker: `problems/hadwiger_nelson/hn.py` (exact) plus
an independent 120-digit recount.
Run: `python problems/hadwiger_nelson/glue2.py --trials 20 --budget 40`.
Artifact: `runs/hn_glue2.json`.

## Result

| quantity | value |
|---|---|
| trials | 20 |
| trials reaching chi = 5 | **0** |
| best | n = 39, m = 95, chi = 4 |
| ambiguous distances (1e-25 < |d^2-1| < 1e-50) | 0 |

Independent recheck at 120 digits over all C(39,2) = 741 pairs:

```
n 39  m 95  ambiguous 0  json_m 95
chi recomputed 4
```

So the artifact is real: 95 exact unit edges, no near-misses, 4-chromatic.

## The wall

The generator glues rotated copies of the current graph along vertex pairs at
matching distances — not only unit edges. That is the most permissive rigid
operation available, and it still produced chi = 4 in 20/20 trials. Density went
up (95 edges on 39 vertices, average degree 4.87) without chromatic number
moving.

Interpretation, stated carefully:

- This is an empirical wall, not a theorem. Gluing copies of a 4-chromatic
  graph can in principle create a 5-chromatic union, because colouring
  constraints from different copies can conflict.
- It is consistent with the known history: de Grey's 5-chromatic graph was not
  more copies of the Moser spindle. It needed a new gadget (the H-graph) that
  forces a colour identity the spindle cannot express.
- Therefore: **more compute on this space buys nothing.** The space is
  saturated at 4 and the next step is a gadget, not a sample.

## Consequence for orchestration

This is the cleanest instance of the failure mode the repo is built against.
Ten thousand agents sampling this space would return ten thousand valid
4-chromatic graphs and zero progress. The Cartographer's `hit_rate_prior` for
"glue spindle copies" should be 0.0, and the branch should be killed before any
compute is spent.

The correct next Cartographer output is not a bigger sampler. It is a
parameterised family of *colour-forcing gadgets* — graphs with two marked
vertices whose colour relation is prescribed (`equal`, `different`, or
`different-from-all-k`) — plus a composition rule. That is a different space,
and it is the one de Grey actually searched.

## Status

Negative result. No claim about chi(R^2) is made. The graph is a valid
4-chromatic unit-distance graph and nothing more.
