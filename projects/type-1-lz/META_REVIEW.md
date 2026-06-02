# Meta-review — recurring failure patterns (Type-1 N=3 LZ)

The Meta-review agent's standing critique list. **Propagate forward:** `physics-reflection` must
check each of these on *every* future hypothesis; `physics-intuition` (Generation) must avoid
re-introducing them. Updated 2026-06-01.

## M1 — Confident STRUCTURAL claims that fail an independent check
The single most recurring error this session. Instances (all caught only by independent
coordinator computation):
- WS-D: "rigidity / null-space dimension 0" obstruction — arithmetic wrong (null-space is 3); the
  real obstruction is the slope-free width lemma.
- WS-F: "the exact-crossing node is complex / real gaps ≥0.14" — wrong; the node is **real and
  universal**, and WS-F's *own* near_node sample had a real crossing at u≈−3 its bounded search missed.
- Coordinator (me): "the home run needs a *non-commuting* `Ê`" — false premise; Malikis–Cheianov's
  `Ê` is Abelian.
- WS-PA1: "Painlevé V via middle convolution" — refuted by the capstone (rank-3).
- H-R1: "triangular Stokes shear / `e^{2πic_i}` as a matrix factor" — wrong (cyclic; category error).
- WS-PV initially under-counted vs WS-CAP; WS-CAP's own band (T5) accepted on the agent's
  construction, not coordinator-reconstructed.
**RULE:** *Never promote a load-bearing structural claim on a sub-agent's (or one's own) assertion.
Independently recompute the decisive leg first.* This rule caught a real error every single time it
was applied. Corollary: a *correct verdict can ride on a wrong argument* (WS-D) — check the argument,
not just the conclusion.

## M2 — Naming a special-function class by analogy, not by computing OUR invariants
We labeled the answer "Kampé de Fériet," then "Painlevé V," then "confluent-Heun" — each imported
from an adjacent solved model (BBGY hyperbolic; the 2×2 PV) **without first computing the genus /
rank / monodromy of the actual Type-1 connection.** Each was the wrong rank (genus-1 vs genus-0;
rank-2 vs rank-3). Only direct computation (genus of `Σ`; the cyclic-Wronskian; the wild-character-
variety dimension 6 vs 2) settled it.
**RULE:** *Before naming the special-function class, compute the governing curve's genus, the
connection's rank/singularity structure, and the wild-character-variety dimension for the actual
model. Do not import the class from a neighbor.*

## M3 — Notation hygiene must precede claim-scaling
Three collisions corrupted communication (not the numbers): `Γ` (BE exponent vs Cauchy form factor),
`S₁₂` (conflated `P_{m→m}` diagonal probability, the amplitude `C_mm`, and the Stokes multiplier `σ`),
and channel ordering (ε-index vs slope vs energy). Each forced a cleanup pass.
**RULE:** *Pin conventions in `NOMENCLATURE.md` before scaling claims/agents. State, per quantity:
basis (ε-index), the physical role (slope), and the gauge status (probability vs amplitude).*

## M4 — Classificatory drift vs constructive goal
The deepest strategic miss: the goal was a *concise computable expression* (constructive), but the
workflow drifted into *classifying* the object (name-the-special-function), which by construction
yields a name + a no-go, not a formula. The classification reached "unpublished rank-3" and we
nearly mislabeled that as "no concise form" — when the *constructive* exact-WKB product (two explicit
Weber shears × a junction) was never attempted.
**RULE:** *Keep the deliverable form in view. If the goal is a computable expression, prefer
constructive workstreams (build from explicit local pieces) over classificatory ones (name the
transcendent). "Not a single published special function" ≠ "no concise constructive formula."*

## M5 — Random parameter search is the wrong tool for matched-constraint corroboration
Two attempts to find matched-window/different-χ pairs (the rank-3 signature) by random sampling found
nothing usable — such pairs are measure-rare.
**RULE:** *To test "does X depend on Y at fixed Z," use a constrained deformation (Newton-project to
hold Z fixed, sweep Y), not random search.*

## M6 — Distinguish evidence tiers honestly; do not let "verified on N samples" become "proven"
Several pillars are symbolic on 2–4 samples + a structural argument (node=accessory; the `{0,1,3}`
scheme; accessory algebraicity), and one decisive result is gold-gated-numerical (T5) with analytic
support not machine-proved (T1 dim count). These are strong but **not theorems.**
**RULE:** *Label sample-verified + structural-argument as `numerically/analytically-supported`, not
`established`; list the general proof as an explicit owed item.*
