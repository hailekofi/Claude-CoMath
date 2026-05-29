---
name: physics-numerics
description: >-
  Designing well-constructed numerical experiments to validate or reject physics
  conjectures and promote them toward near-proof. Use when testing a conjecture
  numerically, running a simulation, computing a quantity, searching for a counterexample,
  fitting a scaling exponent, or checking a formula against data. Enforces convergence
  studies, error bars, independent cross-checks, reproduction of known limits, and
  artifact elimination. Backends: Python (NumPy/SciPy/SymPy/mpmath), SageMath, and
  Wolfram/Mathematica. This is how numerics earn the status of evidence.
---

# Physics Numerics

In theoretical physics, a careful numerical experiment is often the decisive arbiter. But
**numerics only count as evidence when they are constructed to be falsifiable and
artifact-free.** This skill turns "I ran a quick simulation" into a result that can sit at
`numerically-supported` (≈ near-proof) on the evidence ladder.

## Choosing a backend

- **Python — NumPy / SciPy / SymPy / mpmath**: default for general numerics, linear
  algebra, ODE/PDE, optimization, Monte Carlo, and **arbitrary-precision** checks
  (`mpmath`) when round-off is suspect.
- **SageMath**: algebra, number theory, combinatorics, exact/symbolic computation,
  algebraic geometry, special structures.
- **Wolfram / Mathematica** (`wolframscript`): symbolic integration, special functions,
  closed-form simplification, high-precision arithmetic, and quick exact cross-checks.

Pick the tool that makes the *cross-check* easy — using a second tool for an independent
method is itself a validation step.

## The validation checklist (mandatory before claiming evidence)

A conjecture is only `numerically-supported` once it survives:

1. **Convergence study.** Vary every discretization/truncation knob (grid spacing, time
   step, basis size, cutoff, Monte Carlo samples, working precision) and show the result
   converges. Extrapolate to the continuum/infinite limit where relevant.
2. **Error quantification.** Report uncertainties — statistical (error bars, bootstrap)
   and systematic (discretization, finite-size). A number without an error bar is not a
   result.
3. **Independent methods.** Reproduce the key number by ≥2 genuinely different methods or
   tools (e.g. direct sum vs. integral transform; SciPy vs. Mathematica; two algorithms).
4. **Known-limit reproduction.** Recover analytically known special cases, exactly
   solvable limits, or published values. If it gets the known cases wrong, stop.
5. **Parameter sweep.** Test across the regime of interest, not one lucky point. Look for
   where the claim breaks — that defines its regime of validity.
6. **Artifact elimination.** Rule out round-off (raise precision), finite-size effects
   (scale system size), discretization bias, RNG seeding issues, and silent overflow/NaN.

See `templates/validation-checklist.md` for the per-experiment form to fill in.

## Reproducibility (non-negotiable)

- Save the script as a durable artifact; record tool **versions** and the **random seed**.
- Use the `templates/experiment-template.py` scaffold (header states the hypothesis,
  method, parameters; results + plots saved to files).
- Outputs (numbers, tables, plots) are saved alongside the script, not just printed.

## Output and hand-off

Report: the quantity, its value **with uncertainty**, the checklist items passed/failed,
and the resulting **evidence-ladder level**. Be explicit and honest:
- Passed the full checklist → `numerically-supported` (near-proof); state the regime.
- Mixed/partial → `numerically-suggestive`; say what is still needed.
- Contradicts the conjecture → `refuted`; record *why* in `physics-research-log`.

Hand strong numerical evidence to `physics-derivation` (to seek an analytic explanation)
and log everything in `physics-research-log`.
