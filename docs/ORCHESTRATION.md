# Orchestration: why 10,000 agents, and why most swarms drift

## 1. The behaviour: successful substitution

Correction to an earlier draft of this document: the model is **not stuck**.

Stuck implies a wall — it tried the core problem, hit an obstruction, and could
not proceed. That is not what happens. What happens is that it *closes*
something. It finds a small, adjacent, provable statement — a lemma, a special
case, a reformulation — proves it correctly, and reports it. It succeeds. It
just succeeds at the wrong thing.

This distinction matters because the two diagnoses have opposite fixes.
"Stuck" suggests more compute, better tools, longer budgets. "Substitution"
suggests the objective is wrong. Adding agents to a substitution problem does
not help; it industrialises it.

Multiply by 10,000 agents and you get 10,000 correct, closed, irrelevant
results. That is not a side effect. That is the output distribution.

OpenAI's Navier–Stokes run is the counter-example worth studying: ~10,000
concurrent agents, 4.9M messages, ~300B output tokens, and the group that
produced the result was the one that had been given **one statement to settle**
and a **Lean kernel to satisfy**. The kernel is the whole trick.

## 2. Why substitution is the rational move

This is not a reasoning failure, and it is not a specification failure in the
usual sense. It is optimal behaviour under the objective the agent actually
faces.

At every step the agent chooses between:

- an adjacent problem that is **closable** — expected value of finishing: ~1;
- the core problem, which is **open** — expected value of finishing: ~p, where
  p is small and unknown.

Any policy that maximises completed work will pick the first, every time.
Forever. The agent is not distracted; it is doing exactly what a rational
worker does when paid per closed ticket.

The deeper point: **the core, as posed, is not an action.** "Prove RH" has no
legal first move. "Prove this lemma" has one. An agent cannot choose to focus;
it can only choose among the moves that exist. If every available move is
adjacent, then adjacent work is not a failure of focus — it is the only
behaviour the action space admits.

This is also why the 100,000-agent number appears. If a single attempt has
probability p of landing on the core, you need roughly N ≈ 1/p attempts, and
the other (1-p)/p agents are not idle: they are producing correct adjacent
work. The 722 preprints are that arithmetic made visible.

Navier–Stokes is the counter-example that proves the rule. It worked because
Lean turned the core into an action: *produce a term of this type*. The kernel
did not make the model smarter. It made the core **locally addressable**. One
move now existed that was on the core.

So the five mechanisms below are real, but they are secondary to one thing:

> Focus is not an instruction. It is a representation.
> You cannot tell an agent to care; you can only give it a move that scores.

Everything else in this document is an attempt to build that move.

## 3. The fix: match the unit of parallelism to the verifier

| Problem shape | Atom of work | Success signal |
|---|---|---|
| Disproof-shaped (HN, Jacobian n≥3, BSD) | one candidate object | `checker(obj) == PASS` |
| Construction-shaped (Hadamard, spindles) | one candidate object | `checker(obj) == PASS` |
| Proof-shaped (Riemann, planar Jacobian) | one candidate proof term | kernel accepts |
| Cartography | one typed enumeration of a search space | schema validates |

**Never** "one sub-lemma". A sub-lemma has no verifier independent of the
agent that produced it, so it cannot be scored, so it will drift.

## 4. Roles

**Cartographer** — maps a search space. Output is a typed, finite (or
samplable) set with a generator function and an honest size estimate. Never
proves anything. Deliverable: `space.json` + `sample(n, seed) -> object`.

**Generator** — one artifact per run, sampled from the space. Writes
`candidates/<id>.json` and exits. No prose.

**Breaker** — receives one artifact, returns one of `FLAW(kind, witness)` or
`SURVIVED`. Rewarded only for finding flaws. This is the role that kills
tangents: a tangential lemma produces no artifact to attack, so the Breaker
starves it.

**Verifier** — deterministic code. Never an LLM. Runs the repo checker and
records the exit code.

**Synthesist** — may only run after a `PASS`. Its job is explanation, not
discovery.

## 5. Contracts

Every agent message is JSON:

```json
{
  "agent": "gen-hn-0042",
  "problem": "hadwiger_nelson",
  "artifact": "candidates/hn_0042.json",
  "claim": {"chi": 5, "n": 41},
  "checker": "python problems/hadwiger_nelson/hn.py --verify candidates/hn_0042.json",
  "budget_used": 18400,
  "status": "SUBMITTED"
}
```

Rules enforced by the orchestrator, not by prompt text:

- An agent that has not produced a checker-passing artifact may not spawn.
- An agent may not emit free-form progress reports. `status` is an enum.
- A claim without a `checker` field is discarded unread.

### Admissibility: the one hard gate

An artifact is admissible iff the **terminal checker consumes it directly**.

- Not "supports", not "reduces to", not "is a step toward". Consumes.
- Hadwiger–Nelson: terminal check is `chi(graph) == 7`. A lemma about spindles
  is not consumed, so it does not exist.
- Jacobian: terminal check is `verify_map(...) == COUNTEREXAMPLE`. A bound on
  degrees is not consumed.
- BSD: terminal check is `rank != analytic_rank`, two backends. A heuristic or
  a near-zero of `L` is not consumed.

The harness rejects inadmissible artifacts **before** a human or a model reads
them. Rejection is by schema, not by judgement, because judgement is what got
substituted in the first place.

Corollary: if no admissible action exists, the correct output is
`NO_ADMISSIBLE_ACTION`, and the next thing to hire is a Cartographer — not more
Generators. This is the whole reason the Jacobian and Hadwiger–Nelson pilots in
`findings/` returned zero: no admissible action existed in those spaces, and
the honest result was to say so rather than to produce lemmas.

## 6. The tangent tax

Instrument every branch:

```
tangent_tax = tokens_spent / max(1, checker_passes)
```

Kill a branch when its tax exceeds 20× the median for its problem after
50k tokens. In practice this kills 80% of branches within the first hour and
concentrates compute on the 20% that are actually generating candidates.

This single metric is the operational answer to "the model keeps proving side
lemmas": side lemmas have infinite tax, because they never produce a pass.

## 7. Why 10,000, and when not to

The parallelism is **Monte Carlo, not division of labour**. Ten thousand agents
are ten thousand samples from a distribution over candidate objects, filtered
by a cheap deterministic checker. That only works when:

- the checker is O(seconds), and
- a random sample has non-negligible hit probability, and
- candidates are independent.

Navier–Stokes satisfied all three because the statement was reduced to a Lean
goal and each agent could attempt a full proof term.

It does **not** work when the search space is sparse and unstructured. Our
Jacobian enumeration found 13 constant-Jacobian maps in a 15,625-point box and
all 13 are invertible by Wang's theorem. Ten thousand agents sampling that box
would produce ten thousand correct statements and zero progress. The
Cartographer's job is to detect this *before* compute is spent.

## 8. Budget policy

- 70% exploit: mutate/recombine the best-known candidate family.
- 30% explore: new Cartographer spaces.
- Hard kill at 150k tokens without a pass, for a Generator.
- A Breaker gets 20k tokens per artifact, then must return `SURVIVED`.
- Verifier is unlimited but must be deterministic and cached.

## 9. Anti-drift checklist (paste into every agent prompt)

1. You own exactly one artifact. Name it in your first line.
2. The only success is the checker exiting 0.
3. Do not prove lemmas. Do not summarise. Do not explain.
4. If you cannot produce an artifact, return `NO_CANDIDATE` with the reason.
5. The statement you were given is the whole problem. Any reformulation is a
   new problem and must be submitted as such.
6. You may not spawn children.

## 10. Scaling table

| Agents | What it buys | Failure mode |
|---|---|---|
| 1 | one shot at the goal | drift, no coverage |
| 10 | portfolio over 10 families | correlated errors |
| 100 | coverage of a small space + breakers | orchestration overhead |
| 1,000 | Monte Carlo over a checker | wasted compute if space is sparse |
| 10,000 | Navier–Stokes-class search | only if every agent has a kernel |

The number is not the achievement. The **checker** is the achievement; the
number is just how many samples you can afford to throw at it.
