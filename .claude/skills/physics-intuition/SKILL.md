---
name: physics-intuition
description: >-
  Physical-reasoning toolkit and intuition-driven hypothesis generation for theoretical
  physics. Use when generating or sharpening conjectures and research directions, asking
  "why does this happen / what if", sanity-checking a claim or formula, estimating an
  order of magnitude, finding the right limit or regime, or reasoning about symmetry,
  scaling, and dimensional analysis — before committing to heavy calculation. Runs a
  lightweight generate–debate–evolve pass and writes hypotheses to a hypothesis card.
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

## Generate–debate–evolve (hypothesis loop)

For an open question, don't anchor on the first idea:

1. **Generate** ≥3 candidate hypotheses or explanations using different tools above.
2. **Debate** each: argue for and against; check correctness, novelty vs known results,
   and — crucially — **testability** (is there a cheap numerical or analytic test?).
3. **Rank** by *plausibility × testability*. Prefer hypotheses that are easy to kill.
4. **Evolve** the survivors: combine, simplify, sharpen, or take a cleaner limit.
5. **Meta-review**: state what this round established, what to test next, and which ideas
   to discard (and record the discards via `physics-research-log`).

## Output

For each surviving hypothesis, produce a **hypothesis card** (see
`templates/hypothesis-card.md`) capturing: the statement, the physical intuition behind
it, its **regime of validity**, the cheapest **falsification test**, and an initial
evidence-ladder status (usually `conjecture`).

Then hand off:
- to `physics-numerics` to run the falsification/validation test,
- to `physics-derivation` for an analytic argument,
- to `physics-literature` to check whether it is already known.

Log every hypothesis and verdict in `physics-research-log`.
