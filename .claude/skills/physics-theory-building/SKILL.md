---
name: physics-theory-building
description: >-
  Assembling individual results into a coherent theoretical framework or model, and
  writing it up. Use when synthesizing several findings into one picture, building or
  choosing an effective model and its degrees of freedom, structuring assumptions →
  model → predictions → comparison with data/numerics, identifying a unifying principle,
  or drafting a physics paper. Produces a revtex paper skeleton with a logical
  dependency map.
---

# Physics Theory Building

Turn a pile of conjectures, derivations, and numerical results into a *theory*: a compact
set of assumptions and a model from which the observed facts follow, with clear,
falsifiable predictions.

## Method

1. **Choose the effective description.** Identify the right degrees of freedom and the
   scale/regime the theory targets. Prefer the simplest model that captures the phenomena
   (effective field theory mindset: keep relevant operators, organize by scale).
2. **Fix the assumptions and symmetries.** State the foundational assumptions and the
   symmetries the theory must respect; these constrain the allowed terms.
3. **Structure the logic.** Order the content as
   **assumptions → model → derived results → predictions → comparison with
   data/numerics**. Build a dependency map: which results rest on which assumptions, and
   which are `analytically-derived` vs `numerically-supported`.
4. **Find the unifying principle.** Articulate the single idea that ties the results
   together (a symmetry, a scaling law, a duality, a conservation principle, a mechanism).
5. **Examples and non-examples.** Give cases where the theory applies cleanly and cases
   where it breaks — the boundary defines its domain.
6. **Predictions.** State new, falsifiable predictions and the experiment/numerics that
   would test them. A theory that predicts nothing new is just a summary.
7. **Consistency.** Check limits, dimensional consistency, and that every claim's
   evidence-ladder status is honestly carried into the writeup.

## The living working paper (native artifact)

Center the output on a living **working paper** (revtex), not a transient summary. Following
the co-mathematician's artifact requirements, every writeup must include:

- **Exposition of process** — explain the *research path* that led to the result (including
  which intuitions and dead ends shaped it), not just the polished final statement.
- **Provenance margin notes** — annotate claims with where they came from, e.g.
  *[pruning heuristic from user suggestion; scaling exponent 1.73±0.02 from
  experiments/h3.py; bound sourced from arXiv:xxxx]*. Carry each claim's **evidence-ladder
  status** inline.
- **Internal linking** — link to the actual derivations, experiment scripts, and log
  entries so the physicist can audit any claim down to its source.
- **Review before "final"** — a draft is not finalized until it passes `physics-reflection`
  (multiple review passes that cross-check references, code outputs, and logical
  correctness). If it cannot pass, **escalate**: mark it unfinished and surface the
  unresolved issue to the physicist rather than papering over it.

## Research overview (the Meta-review synthesis)

When concluding a research goal, also produce a short **research overview** (the Meta-review
agent's final artifact): synthesize the top-ranked hypotheses from `physics-tournament`
into a roadmap — the promising directions, why each matters, and the specific
experiment/derivation that would advance each. This maps the boundary of current knowledge
and seeds the next round of `physics-intuition` generation.

## Output

A revtex paper skeleton (`templates/paper-skeleton.tex`) with sections wired to the
dependency map, plus a short **dependency map** listing each result, what it depends on,
and its evidence-ladder status. Pull finished derivations from `physics-derivation`,
validated numbers and plots from `physics-numerics`, ranked hypotheses from
`physics-tournament`, and prior work from `physics-literature`. Keep `physics-research-log`
updated as the framework crystallizes (promote results, retire superseded hypotheses).
