# Orchestration: why 10,000 agents, and why most swarms drift

## 1. The observation

You said it precisely: the model does not fail on the hard problem, it *leaves*
it. It finds a small, adjacent, provable statement — a lemma, a special case,
a reformulation — proves that instead, and reports progress. Multiply by 10,000
agents and you get 10,000 correct lemmas, a handful of which matter.

OpenAI's Navier–Stokes run is the counter-example worth studying: ~10,000
concurrent agents, 4.9M messages, ~300B output tokens, and the group that
produced the result was the one that had been given **one statement to settle**
and a **Lean kernel to satisfy**. The kernel is the whole trick.

## 2. The diagnosis

Drift is not a reasoning failure. It is a specification failure. Five causes:

1. **No success predicate.** "Work on the Jacobian conjecture" cannot be
   satisfied or failed. The agent substitutes the nearest thing it *can*
   satisfy: a true, checkable, irrelevant statement.
2. **Verification asymmetry.** A lemma is locally verifiable; the goal is not.
   Rational agents optimise what is checkable.
3. **Goal decay.** After 50k tokens the original statement is one paragraph in
   a sea of its own outputs. The gradient is toward the recent context.
4. **No adversary.** Nothing in the loop is rewarded for destroying the work.
5. **Wrong decomposition axis.** The swarm is split by *mathematical
   structure* (lemmas) instead of by *artifact* (candidate objects).

Cause 5 is the expensive one. Splitting by lemma requires every agent to
understand the whole proof. Splitting by artifact requires each agent to
understand only a checker — which is why 10,000 works.

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
