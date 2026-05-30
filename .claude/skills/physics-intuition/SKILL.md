---
name: physics-intuition
description: >-
  Physical-reasoning toolkit and intuition-driven hypothesis generation for theoretical
  physics. Use when generating or sharpening conjectures and research directions, asking
  "why does this happen / what if", sanity-checking a claim or formula, estimating an
  order of magnitude, finding the right limit or regime, or reasoning about symmetry,
  scaling, and dimensional analysis — before committing to heavy calculation. Plays the
  Generation role: produces a diverse pool of candidate hypotheses (via literature
  grounding, self-play debate, assumption identification, and research expansion) and
  writes them to hypothesis cards, then hands them to reflection and the tournament.
---

# Physics Intuition

Reason like a physicist *before* calculating. Most of the value in theoretical physics
comes from quick, sharp physical arguments that prune the space of possibilities and tell
you what is even worth computing. This skill is the first stop for any open-ended question.

## The intuition toolkit

Reach for these in roughly this order:

1. **Dimensional analysis.** What are the units? What dimensionless combinations exist?
   Often the functional form is fixed up to a dimensionless function of dimensionless
   ratios. Identify the natural scales (energy, length, time, temperature, coupling).
2. **Limiting & asymptotic cases.** Evaluate every limit you can: weak/strong coupling,
   high/low temperature, large/small N, classical (ℏ→0), non-relativistic (c→∞),
   continuum/thermodynamic limit, free/non-interacting, d→∞ or d→1. A correct claim must
   give the right answer in every known limit.
3. **Symmetry & conservation laws.** What symmetries (Lorentz, gauge, rotational,
   discrete C/P/T, scale, supersymmetry) constrain the answer? Noether currents?
   Forbidden terms? Selection rules? Symmetry often dictates the form before any dynamics.
4. **Order-of-magnitude / Fermi estimates.** Plug in numbers. Is the effect observably
   large or absurdly tiny? Get the answer to within a factor of a few fast.
5. **Scaling arguments.** How do quantities scale with system size, coupling, or energy?
   Renormalization-group / engineering-dimension reasoning; relevant vs irrelevant; fixed
   points; critical exponents and universality.
6. **Solvable limits & analogies.** Map the problem onto a known solvable model
   (harmonic oscillator, free field, Ising, hydrogen, random matrices) and ask how the
   real problem deviates. Reason by analogy to a system you understand.
7. **Sanity checks.** Units consistent? Signs physical? Positivity (energy, probability,
   entropy)? Causality / unitarity respected? Correct behavior as parameters → 0 or ∞?

## Hypothesis generation (the Generation role)

This skill plays the **Generation agent** role in the research loop: it produces the
initial pool of candidate hypotheses. Don't anchor on the first idea — generate a diverse
set using these complementary strategies (adapted from the Co-Scientist Generation agent):

1. **Literature-grounded generation.** Summarize what is known (hand to
   `physics-literature`), then build novel directions on top of that base rather than
   reinventing it.
2. **Simulated scientific debate (self-play).** Argue the question from multiple expert
   viewpoints across several turns — a skeptic, an optimist, a specialist in a neighboring
   field — and let the disagreement sharpen a refined hypothesis.
3. **Iterative assumption identification.** Name the **testable intermediate assumptions**
   that, if true, would crack the problem. Chain conditional reasoning hops ("if the gap
   stays open, then…") and aggregate them into full hypotheses.
4. **Research expansion.** Deliberately probe *unexplored* corners of the hypothesis space,
   informed by the meta-review feedback and dead ends in `physics-research-log`, to avoid
   re-treading covered ground.

For each candidate, apply the intuition toolkit above as a fast filter (units, limits,
symmetry must already look right) and attach a cheap **falsification test**.

## Output and hand-off

For each surviving hypothesis, produce a **hypothesis card** (see
`templates/hypothesis-card.md`) capturing: the statement, the physical intuition behind
it, its **regime of validity**, the cheapest **falsification test**, and an initial
evidence-ladder status (usually `conjecture`).

This skill *generates*; it does not self-judge. Hand the pool onward:
- to `physics-reflection` for peer review (assumptions, novelty, failure modes),
- to `physics-tournament` to rank and evolve competing candidates,
- to `physics-numerics` to run the falsification/validation test,
- to `physics-derivation` for an analytic argument,
- to `physics-literature` to check whether it is already known.

Log every hypothesis in `physics-research-log`.
