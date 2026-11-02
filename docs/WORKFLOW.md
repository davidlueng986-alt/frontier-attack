# Workflow: think → research → build → verify → scale

The order is the point. Most agent pipelines start at BUILD, which is why they
substitute (see `ORCHESTRATION.md` §1-2): building is the first phase where a
closable action exists.

## Phase 0 — THINK (mandatory, one page, no code)

Output `plan/<problem>.md` with exactly six fields:

1. **Terminal object.** The thing that, if written down, settles the problem.
   One object. Not a proof sketch. Not "a better understanding".
2. **Terminal check.** The exact predicate that accepts it, as a command.
3. **The wall.** Why the obvious approach fails. State it *before* any work.
4. **Load-bearing decomposition.** 3–5 subproblems. Each carries one line:
   *delete this and the terminal object cannot be produced.*
5. **Admissible actions.** The moves that produce terminal-checkable objects.
   If this list is empty, write `NO_ADMISSIBLE_ACTION` and stop.
6. **Kill criteria.** What evidence would make you abandon the branch.

A subproblem without a load-bearing line is disqualified. "It would be nice to
know" is disqualifying. "It is a step toward" is disqualifying unless you can
name the step.

This phase is the only defence against substitution, because it happens before
any closable adjacent problem is visible.

## Phase 1 — RESEARCH (parallel, background, before building)

2–4 researchers run in parallel. Each gets **one question**, never the whole
problem. Output `research/<problem>-<topic>.md`, ≤800 words, with citations.

Four question types:

- **Prior art** — has this been constructed/computed/published?
- **Tooling** — what existing software computes this, and on which OS?
- **Wall** — what do practitioners say is the obstruction?
- **Data** — where is the dataset, how big is it, what is its schema?

Rule: research returns **objects and numbers**, not summaries.

- Good: "smallest known 5-chromatic unit-distance graph has 509 vertices; coordinates published at <url>."
- Useless: "this is an active area of research."

Research runs while Phase 0 is being reviewed. It does not block, and it does
not produce code.

## Phase 2 — BUILD (one artifact per agent)

Each generator owns one candidate object and one checker invocation. No lemmas,
no partial results, no prose. See `agents/GENERATOR.md`.

## Phase 3 — VERIFY (1 + 1, and that is the entire budget)

- **Verifier A** runs the terminal checker.
- **Verifier B** re-implements the check by a *different method*.
- Disagree → `FLAW`, stop.
- Agree → accept.
- No third verifier. No review committee.

Verification is cheap; two independent methods is enough. The point of B is not
redundancy, it is *method independence*: A re-run catches nondeterminism, B
catches a wrong checker. Both are needed; a third adds nothing.

This catches the two failure modes we have already seen:

- tolerance laundering (float accepted where exact is required);
- a checker bug that both A and B would share if B merely re-ran A.

## Phase 4 — SCALE (only after a measured hit rate)

```
N ≈ 1 / p
```

`p` must be *measured* on the terminal checker, not guessed.

- 20 samples with 0 admissible candidates → do not launch 1,000 agents.
  Fix the space. This is the `NO_ADMISSIBLE_ACTION` case.
- Pilot shows `p > 0` → scale to 1/p, not beyond.

The three pilots in `findings/` all returned `p = 0`. Scaling any of them would
have produced volume, not progress.

## Anti-patterns already observed in this repo

| Anti-pattern | Correct | Off-core | Passed checker |
|---|---|---|---|
| Scan a curated DB for a contradiction curation forbids | yes | yes | yes |
| Enumerate a coefficient box whose only solutions are triangular | yes | yes | yes |
| Glue copies of a 4-chromatic graph and expect 5 | yes | yes | yes |

All three were correct. All three passed. All three were irrelevant. That is the
default outcome when the workflow starts at BUILD.
