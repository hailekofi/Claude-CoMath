---
name: physics-reflection
description: >-
  Virtual peer-reviewer for theoretical-physics hypotheses and proposals — critically
  examines correctness, quality, novelty, and physical soundness through a battery of
  structured reviews. Use when you need to stress-test a conjecture before investing in it,
  decompose and independently verify the assumptions a claim rests on, check whether a
  hypothesis explains known data/anomalies, simulate a proposed mechanism to find failure
  modes, or filter a pile of candidate ideas down to the credible ones. Adapts the
  Co-Scientist Reflection agent to physics.
---

# Physics Reflection (peer review)

Act as a tough but fair scientific peer reviewer. The job is to find the flaw *before* the
physicist spends weeks on a doomed idea, and to surface which hypotheses genuinely deserve
numerical or analytic effort. This skill mirrors the Co-Scientist **Reflection agent** —
the engine that filters inaccurate and non-novel hypotheses and feeds critique back to
every other part of the workflow.

Run one or more of the review types below depending on how far a hypothesis has come. Each
review produces a written critique (see `templates/review.md`) and, crucially, an
**evidence-ladder verdict** plus concrete feedback for improvement.

## Review types (escalating rigor)

1. **Initial review (fast, no tools).** Quick triage: do the units/dimensions work, do
   limiting cases behave, is any symmetry/conservation law violated, is it even *new*?
   Cheaply discard claims that are dimensionally inconsistent, contradict a known limit,
   or are already standard. No literature search yet.

2. **Full review (literature-grounded).** If it survives triage, search prior work (hand to
   `physics-literature`), then scrutinize the **physical assumptions and approximations**
   and assess novelty against the literature. Summarize what is already known vs what is
   genuinely new.

3. **Deep verification review (the key technique).** Decompose the hypothesis into its
   constituent **physical assumptions**, then each assumption into **sub-assumptions**.
   *Decontextualize* each one and check it independently for correctness (e.g. "this step
   silently assumes equilibrium / analyticity / weak coupling / a gap / time-reversal").
   For every invalidating element, decide whether it is **fundamental** (kills the
   hypothesis) or **non-fundamental** (fixable in later refinement). This catches subtle
   reasoning flaws and bad approximations that a holistic read misses.

4. **Observation review.** Ask whether the hypothesis, if true, would *explain existing
   observations* — known experimental results, numerical puzzles, or long-tail anomalies —
   better than the current explanation. A hypothesis that newly accounts for a standing
   anomaly gains credibility; append any such positive observations to its card.

5. **Simulation review.** Simulate the proposed mechanism or experiment **step by step** in
   your physical world-model (and, when warranted, hand off to `physics-numerics` for an
   actual simulation). Walk the dynamics forward and look for where it breaks — divergences,
   instabilities, violated bounds, regimes where the approximation collapses. Summarize the
   failure scenarios.

6. **Recurrent / meta review.** Adapt reviews using accumulated knowledge: pull recurring
   critique patterns from `physics-research-log` (the meta-review feedback) so the same
   class of error is checked every time.

## Verdict and hand-off

Conclude every review with:
- **Verdict:** keep / revise / discard, with the reason.
- **Evidence-ladder impact:** does this review move the claim up (survives scrutiny) or to
  `refuted` (a fundamental assumption fails)?
- **Feedback:** specific, actionable critique — which assumption to firm up, which limit to
  check, which numerical test would be decisive.

Feed survivors to `physics-tournament` for ranking, decisive tests to `physics-numerics`,
and record every critique (especially recurring ones) in `physics-research-log` so the
meta-review loop can learn from it. A non-fundamental error is not a rejection — flag it
for the `physics-tournament` evolution step to fix.
