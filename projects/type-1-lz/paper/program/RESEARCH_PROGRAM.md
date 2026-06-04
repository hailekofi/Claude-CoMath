# Type-1 N=3 Landau–Zener: an autonomous research program

**Status:** draft for review & discussion (2026-06-01). **Owner:** co-physicist working session.
**Branch:** `claude/comathematician-agent-skills-9HPSJ`. **Shared state:** `RESEARCH_LOG.md`.
**Companion artifacts:** `paper/drafts/type1_lz_working_paper.tex` (living paper), `paper/notes/gate_test_genus.md`,
`paper/notes/REVIEW_kz_line.md`, `experiments/{anchor_experiment,ring_structure}.{py,md}`.

This plan is written so that the workstreams below can be dispatched **autonomously** to the
co-physicist sub-skills (intuition / literature / numerics / derivation / reflection / tournament /
research-log), with explicit inputs, outputs, pass/fail gates, and hand-off rules. Each workstream
is a self-contained agent task; the research log is the single source of truth that lets parallel
tracks coordinate without colliding.

---

## 0. Goal and practicality bar (unchanged)

Closed-form, computationally practical `P` for the `3×3` Type-1 MLZ matrix in `{γ_i,ε_i,a_i}` — no
ODE matrizant, ≤ a handful of geometric integrals. The two extreme survivals are exact (BE); the
open content is the **middle survival `P_{m→m}`** and **one independent off-diagonal**.

**Honest target (revised by this session's findings):** a fully-elementary generic closed form is
unlikely; the achievable, publishable deliverable is a **structured exact result**: elementary
exponents (BE + residue data) × a **genus-0 confluent-Heun / Kampé de Fériet** connection-coefficient
prefactor, plus an **elementary exact locus** (wherever the connection factorizes or the node
reduces it), all benchmarked to `10⁻⁹`.

---

## 1. The theory in one paragraph (the unifying principle)

> **`P` is the squared modulus of the Stokes/connection data of a rank-2 irregular meromorphic
> connection living on the Type-1 *spectral curve*, which is genus-0 rational, marked by three
> regular poles `{ε_i}`, one irregular point at `∞`, and a structural node — the exact crossing.
> The commuting-ring (Abelian) structure fixes the elementary part: the adiabatic spectrum, the WKB
> phases, and the BE survivals (residue-collapsing periods). The remaining, *non-Abelian* content —
> the middle survival and off-diagonals — is the connection coefficient threading the turning points
> and the node, a genus-0 confluent-Heun / Kampé de Fériet object. All slope-dependence enters
> linearly, through the dynamical phases, on a fixed geometric connection.**

Everything in the program either (i) makes this statement precise, (ii) computes the prefactor, or
(iii) maps the locus where it degenerates to elementary.

---

## 2. Current state — evidence ladder (carried into every workstream)

**ESTABLISHED (proof / machine-precision numerics):**
- E1 BE extreme survivals `= ∏e^{−2πΓ_ij}`, realized as residue-collapsing window periods. [assay]
- E2 Commuting family is a **ring**; the time-quadratic partner reduces to the linear family with
  linear-in-`u` coefficients; the family **shares one eigenbasis**. ⇒ **Abelian ceiling**: the
  commuting partner gives only spectrum + `τ`-invariance, never the prefactor.
  [`experiments/ring_structure.py`]
- E3 **Spectral curve `Σ` is genus-0 rational**, via Gaudin parametrization `(E,u)=(m/p,n/p)`; a
  **structural node = the exact crossing** at real `u*` (two eigenvalues exactly degenerate) + 4
  simple branch points (avoided crossings). [`paper/notes/gate_test_genus.md`]
- E4 `P^(a) = |holonomy of a fixed algebraic connection `W_ij=⟨φ_i|φ_j'⟩` (a-independent) twisted by
  a-linear phases|²`. [`experiments/ring_structure.py`]
- E5 Middle-survival prefactor **dominates**: `P_mid/P_mid^{incoherent} = 2.5–104×`, growing with
  crossing overlap. [`experiments/anchor_experiment.py`]

**ANALYTICALLY-DERIVED / FRAMEWORK:**
- F1 `P` = Stokes/connection data of a rank-2 irregular ODE at `u=∞` (N=2 = Weber exactly).
- F2 C–S `τ`-deformations are **isomonodromic** (preserve `P`).
- F3 Function class is **genus-0 confluent-Heun / Kampé de Fériet**, NOT elliptic/Painlevé.

**CONJECTURE / OPEN (the program's targets):**
- O1 Closed form of the prefactor (the confluent-Heun/Kampé de Fériet connection coefficient).
- O2 Whether the **exact-crossing node** reduces accessory parameters / pins `P_{m→m}`.
- O3 Whether the 3rd-order `λ`-connection **factorizes** on Type-1 sub-loci (⇒ elementary there).
- O4 Whether a **non-commuting zero-curvature `Ê`** exists for Type-1 (⇒ exact product `S=∏S_ij`).

**FALSIFIED (do not revisit):** incoherent product; two-Dykhne; coherent three-rotation; algebraic
frame maps; point-local Airy; BE-extended diagonals.

---

## 3. Dependency map

```
                         [WS-F oracle: 1e-9 harness + strata]  (foundational, continuous)
                                          |
                    +---------------------+----------------------+
                    |                                            |
        [WS-A Riemann scheme of the connection]        [WS-D non-Abelian Ê search]  (high-risk)
                    |  (gate 1: 3-pt vs 4-pt / irregular type)            |  (gate 4)
        +-----------+-----------+                          exact product OR
        |                       |                          obstruction proof
 [WS-B node reduction]   [WS-C factorization test]
   (gate 3: P_{m→m}?)        (gate 2: elementary locus?)
        |                       |
        +-----------+-----------+
                    |
        [WS-E closed-form prefactor assembly]  --->  [working paper + tournament re-rank]
```

Rule: a workstream may start when its inputs reach the stated evidence level; results promote only
after a `physics-reflection` pass; the `physics-tournament` re-ranks O1–O4 after every gate.

---

## 4. Workstreams (autonomous task specs)

Each spec is dispatchable as one agent task. Format: **objective · method · skills · inputs ·
outputs · PASS/FAIL gate · parallel? · escalation**.

### WS-F — Benchmark oracle (foundational, continuous)
- **Objective:** a trusted `10⁻⁹` `P(γ,ε,a)` oracle and a stratified parameter suite (well-separated,
  overlapping, near-node, near-degenerate-slope, strong/weak coupling).
- **Method:** port `experiments/anchor_experiment.py` onto the project interaction-picture harness
  (`assay/ip.py`, `benchmark.py`); add the `2.5–104×` middle-survival cases and the exact-crossing
  `u*` locus as named strata; convergence + error bars per `physics-numerics`.
- **Skills:** physics-numerics. **Inputs:** assay package, E5 cases. **Outputs:** `experiments/
  oracle.py`, a strata table in the log.
- **PASS:** reproduces BE to `≤10⁻⁸` on all strata; middle survival stable to `≤10⁻⁶`. **FAIL:**
  cannot converge near the node ⇒ flag a genuine singular-stratum caveat.
- **Parallel:** yes (independent). **Escalation:** if harness disagrees with standalone by `>10⁻³`,
  stop and reconcile before any other workstream trusts numbers.

### WS-A — Riemann scheme of the connection (the class-fixing gate)
- **Objective:** the exact local data of the connection on `P¹_λ`: regular-singular exponents at each
  `ε_i`, the irregular/Stokes structure at `∞`, and the local type at the node `u*`. This decides
  hypergeometric (3-pt) vs confluent-Heun (4-pt) vs reducible.
- **Method:** write the scalar 3rd-order ODE (or the `3×3` first-order system) for the amplitudes;
  read off the indicial exponents at `ε_i` (expect `∝ Γ_ij`), the Poincaré rank and formal/Stokes
  data at `∞`, and the local monodromy/exponents at the node. Cross-check exponents `∝ Γ_ij`
  against the residue structure of E1.
- **Skills:** physics-derivation (lead) → physics-numerics (verify exponents numerically) →
  physics-reflection. **Inputs:** E3, E4, geometry.py. **Outputs:** a Riemann-scheme table; the
  named ODE class.
- **PASS (gate 1):** a definite Riemann scheme with counted accessory parameters. **FAIL:** if the
  scalar reduction introduces uncontrolled apparent singularities ⇒ stay with the `3×3` system and
  use the spectral-network description instead.
- **Parallel:** start after WS-F online; blocks WS-B/C/E. **Escalation:** if the exponents are *not*
  `∝ Γ_ij`, that contradicts E1 — halt and re-examine the dictionary.

### WS-B — Exact-crossing (node) reduction  [rigorous form of old H-B]
- **Objective:** decide O2 — does the node carry computable local connection data that pins `P_{m→m}`
  (or reduces the accessory-parameter count)?
- **Method:** construct the `2×2` local model at `u*` (the crossing pair whose coupling vanishes);
  compute its exact local connection/rotation; compose with the BE data + unitarity + double
  stochasticity; predict `P_{m→m}`. Test against WS-F including the strong-overlap (`104×`) stratum.
- **Skills:** physics-derivation → physics-numerics → physics-reflection. **Inputs:** WS-A scheme,
  WS-F oracle, E3 node location. **Outputs:** candidate `P_{m→m}(γ,ε,a)` + benchmark deltas.
- **PASS (gate 3):** prediction matches oracle to `≤10⁻⁶` on ≥1 stratum ⇒ promote to
  numerically-supported; pursue an analytic proof. **FAIL:** `O(1)` miss on the overlap stratum ⇒
  the node alone does not pin it; log and feed the residual to WS-E.
- **Parallel:** with WS-C after gate 1. **Escalation:** none; both outcomes are results.

### WS-C — Factorization / reducibility of the connection (gate-test #2)
- **Objective:** decide O3 — does the `3×3` `λ`-connection factor (a sub-rep splits off) on any
  Type-1 sub-locus, dropping confluent-Heun → hypergeometric → elementary there?
- **Method:** symbolic factorization of the connection matrix / the scalar ODE over `Q(γ,ε,a)` and
  on candidate sub-loci (e.g. the `ε_CS=0` slice, symmetric-slope `b₁=b₂`, weak-coupling). Map the
  elementary locus explicitly.
- **Skills:** physics-numerics (symbolic, sympy) → physics-derivation → physics-reflection.
  **Inputs:** WS-A. **Outputs:** the factorization locus (possibly empty) + the elementary `P` there.
- **PASS (gate 2):** a non-trivial elementary locus found and benchmarked. **FAIL:** irreducible
  everywhere ⇒ confluent-Heun is genuine; bank that as a structural result.
- **Parallel:** with WS-B. **Escalation:** none.

### WS-D — Non-commuting zero-curvature `Ê` (the home-run / its obstruction)
- **Objective:** decide O4 — does Type-1 admit a **non-commuting** `Ê` (`[Ê,H]≠0`) solving
  `∂_εH−∂_uÊ+i[Ê,H]=0` that separates crossings into an exact product `S=∏S_ij`? (The commuting
  partner cannot, by E2.)
- **Method:** literature pass on Lax/zero-curvature for Gaudin/Cauchy systems (is there a known
  non-Abelian `Ê` for rational Gaudin?); then attempt construction by ansatz on the rational spectral
  curve; if construction fails, attempt an **obstruction proof** (e.g. the absence of a spin-rep /
  Lie-algebra structure that the bow-tie uses).
- **Skills:** physics-literature → physics-derivation → physics-numerics (test path-deformation
  invariance) → physics-reflection. **Inputs:** E2, E3, Malikis–Cheianov. **Outputs:** either an
  exact factorized `S` (benchmark!) or a clean obstruction theorem.
- **PASS (gate 4):** exact `S=∏S_ij` matching oracle ⇒ supersedes WS-E (home run). **FAIL/obstruct:**
  a no-go that confirms the transcendental prefactor is irreducible ⇒ strengthens WS-E's framing.
- **Parallel:** independent high-risk track from the start. **Escalation:** if construction looks
  live but messy, ask the physicist before a multi-day push.

### WS-E — Closed-form prefactor assembly (synthesis)
- **Objective:** O1 — assemble the connection coefficient into a closed form (confluent-Heun /
  Kampé de Fériet of explicit arguments built from the window actions + node data), `P=|S|²`.
- **Method:** use WS-A's Riemann scheme to write the connection-problem solution in the named
  special functions; fix arguments from the geometry (window periods `I_X`, node data); validate
  against WS-F across all strata; recover BE and the WS-C elementary locus as degenerations.
- **Skills:** physics-derivation → physics-numerics → physics-literature (connection formulae) →
  physics-reflection. **Inputs:** WS-A,B,C (and WS-D outcome). **Outputs:** the closed-form `P`
  (in named functions) + the working-paper Results section.
- **PASS:** matches oracle to `≤10⁻⁶` on all strata. **FAIL:** if confluent-Heun connection
  coefficients are not available in closed form, deliver the **integral representation** + the
  controlled asymptotics as the honest result, and flag the genuinely-hard remainder.
- **Parallel:** after gates 1–3. **Escalation:** mark unfinished and surface if it cannot pass
  reflection.

---

## 5. Coordination protocol (how this runs autonomously)

1. **Shared state = `RESEARCH_LOG.md`.** Every workstream reads the log on start and appends its
   result (with evidence-ladder tag + artifact link) on finish. No result is "real" until logged.
2. **Reflection gate.** No claim promotes past `conjecture` → `numerically-supported` →
   `analytically-derived` → `established` without a `physics-reflection` pass recorded in the log.
   Recurring-error checklist applies (esp. *"compute the curve's genus/scheme for OUR model, never
   import from a neighbor"*).
3. **Tournament cadence.** After each gate (1–4), `physics-tournament` re-ranks O1–O4 and the
   workstreams; deprioritized tracks pause, not delete. Leaders get the next compute budget.
4. **Dispatch rule.** WS-F + WS-D + WS-A can launch immediately (WS-A waits only for WS-F online).
   WS-B/C launch on gate 1. WS-E on gates 1–3. A gate-4 PASS (home run) short-circuits WS-E.
5. **Stop conditions.** (i) Any workstream PASS that reproduces the oracle to `10⁻⁶` on all strata
   ⇒ promote to working paper + seek proof. (ii) WS-D obstruction + WS-C irreducible + WS-E integral
   form ⇒ the "structured-but-transcendental" theory is the deliverable; finalize. (iii) Numerics
   disagreement `>10⁻³` anywhere ⇒ global halt, reconcile.
6. **Escalation to physicist.** Any multi-day analytic push (WS-D construction, WS-E special-function
   identification) or any contradiction with E1–E5 pauses for review rather than proceeding.

---

## 6. Phased timeline (critical path)

- **Phase 0 (now):** WS-F online; launch WS-A and WS-D in parallel. *Deliverable: oracle + Riemann
  scheme draft + literature verdict on non-Abelian `Ê`.*
- **Phase 1 (gate 1 passed):** launch WS-B and WS-C. *Deliverable: `P_{m→m}` node-prediction vs oracle;
  the elementary locus map.*
- **Phase 2 (gates 2–3):** WS-E assembles the prefactor; tournament re-rank. *Deliverable: candidate
  closed-form `P` (named functions) + benchmark.*
- **Phase 3:** reflection hardening; working paper to "review-ready"; research overview + next-round
  intuition seeds. *Deliverable: the living paper passes `physics-reflection`.*

## 7. Risk register

| Risk | Likelihood | Mitigation |
|---|---|---|
| Confluent-Heun connection coeffs have no closed form | high | WS-E falls back to integral rep + asymptotics; still a result |
| Scalar reduction spawns apparent singularities (WS-A) | med | stay with `3×3` system / spectral-network description |
| Node does not reduce `P_{m→m}` (WS-B FAIL) | med | feed residual to WS-E; node still cuts accessory params |
| No non-Abelian `Ê` (WS-D obstruct) | med-high | obstruction is itself a publishable structural theorem |
| Harness/standalone mismatch near node | med | WS-F escalation; treat node stratum as singular |

## 8. Definition of done / deliverables

**Decisions (user, 2026-06-01):** (1) WS-F + WS-A + WS-D launched now as parallel background agents.
(2) **WS-D is the emphasized track** (run literature-first as a probe), weighted above the WS-A→E
spine. (3) Target paper = **honest but maximally elegant/parsimonious**; "elementary everywhere" is
too ambitious as a requirement, but hope for the best and present the cleanest true result. (4) WS-F
**ports onto the project interaction-picture harness** (`assay/ip.py`) and makes all improvements
needed so the gold-standard oracle is highly accurate and faithful to the Type-1 construction.


1. **Living working paper** (`paper/drafts/type1_lz_working_paper.tex`) — assumptions → model → results →
   predictions → comparison, with provenance margin notes and links to every derivation/script,
   passing `physics-reflection`.
2. **Closed-form `P`** in named special functions (or the integral rep + asymptotics + elementary
   locus), benchmarked to `10⁻⁹`.
3. **Two structural theorems regardless of O1:** (a) the Abelian ceiling (E2); (b) the genus-0
   rational spectral curve with exact-crossing node (E3) — both new, both publishable.
4. **Research overview** mapping the boundary of knowledge + next-round seeds.

## 9. Predictions already implied (falsifiable, for early testing)

- P-i: the middle-survival prefactor diverges (or peaks sharply) as parameters approach the
  exact-crossing node `u*` — test the `near-node` stratum (the `104×` case is the leading edge).
- P-ii: on any locus where the connection factorizes (WS-C), `P` collapses to a product of
  elementary LZ factors to `10⁻⁹` — a sharp, bimodal prediction (elementary vs transcendental).
- P-iii: the off-diagonal carries a genuine interference phase (relative period between the two
  windows) that is *not* a closed window period — measurable as a deviation from any product form.
