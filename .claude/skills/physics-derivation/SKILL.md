---
name: physics-derivation
description: >-
  Constructing controlled analytic derivations and arguments for theoretical physics —
  rigorous but physical, not formal proofs. Use when deriving a result by hand, doing
  perturbation theory or asymptotic analysis, applying variational or self-consistency
  methods, exploiting symmetry, or turning strong numerical evidence into an analytic
  explanation. Emphasizes distinguishing controlled vs uncontrolled approximations,
  tracking the regime of validity, and cross-checking against numerics. Outputs LaTeX.
---

# Physics Derivation

Produce derivations to the standard of a careful theoretical physicist: every step
physically motivated and checkable, with approximations made *explicit and controlled* —
but without demanding the full machinery of formal proof. The goal is an argument a
referee would accept, paired with numerical corroboration.

## Method

1. **Set up precisely.** State assumptions, the model/Hamiltonian/action, boundary and
   initial conditions, the **regime**, and the **small (or large) parameter** that
   organizes the calculation.
2. **Choose the strategy.** Perturbation theory, asymptotic/WKB analysis, scaling/RG,
   variational or saddle-point, self-consistency (mean field), symmetry/Ward identities,
   Green's functions, effective action. Name it and say *why* it suits the regime.
3. **Derive step by step.** Keep each step physically transparent. Carry units. Track the
   order of every approximation.
4. **Classify approximations.** For each, state whether it is **controlled** (there is a
   small parameter and an estimate of the neglected terms) or **uncontrolled** (an
   uncontrolled truncation/ansatz). Uncontrolled steps are allowed but must be flagged and
   later tested numerically.
5. **Self-audit.** Check every limit (does it recover known/free/classical results?),
   signs, units, symmetry constraints, and positivity/causality. Ask "where would a
   referee object?" and address it.
6. **Cross-check with `physics-numerics`.** A derivation and an independent numerical
   experiment that agree promote the result toward `established`; if they disagree, the
   derivation has a bug or a hidden assumption — find it.

## On rigor (the pivot)

This is **not** a formal-proof skill. Standard physics moves — interchanging limits,
assuming smoothness/analyticity, truncating series, dropping subleading terms — are
permitted *when physically justified and flagged*. State the **regime of validity** for
the final result rather than claiming universal truth. When a fully rigorous proof is
actually needed, say so and scope it separately.

## Output

A clean LaTeX/revtex derivation using `templates/derivation.tex`: assumptions → strategy →
steps → result → regime of validity → list of approximations (controlled/uncontrolled) →
consistency checks performed. Set the result's evidence-ladder status to
`analytically-derived`, or `established` if numerics independently confirm it. Log the
outcome (and any discovered hidden assumptions) in `physics-research-log`, and pass the
finished result to `physics-theory-building`.
