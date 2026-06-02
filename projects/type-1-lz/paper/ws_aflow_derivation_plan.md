# WS-AFLOW — derivation roadmap: the deformation flow of the commuting family

**Type:** `physics-derivation` scoping (a roadmap, **not** a finished proof). **Date:** 2026-06-02.
**Status of the premise it builds on:** *established, machine precision* (deformation_family_probe.py,
ring_structure.py). **Honest headline:** the make-or-break is Milestone 1 (does the flow *close*?); the
realistic best case is a *conditioning/transport* lever, not a free closed form.

---

## 0. The central principle this workstream operationalizes

> **All the transcendental content of the entire Type-1 commuting family is the holonomy of a single,
> `a`-independent, geometrically-computable connection `W`. Each member's open data `{P_mm, b}` is the
> holonomy of that fixed `W` against the member's algebraic adiabatic-energy profile `E^{(a)}(u)`.**

Established inputs: (i) the family `{H^{(a)}}` (fixed `γ,ε`, varying slope `a`) commutes and **shares one
`a`-independent eigenbasis** `φ_i(u)` (`||[H^a,H^{a'}]||~1e-16`, shared-eigvec defect `~1e-15`); (ii)
hence the derivative coupling `W_{ij}(u)=⟨φ_i|φ_j'⟩` is **`a`-independent** and algebraic in `(γ,ε,u)`;
(iii) `{P_mm(a), b(a)}` vary smoothly along the family (the deformation data). The goal: turn "the family
is `W`'s holonomy orbit" into a usable **flow** in `a`.

## 1. Precise setup

Adiabatic frame `ψ(u)=Φ(u)χ(u)`, `Φ=[φ_1\,φ_2\,φ_3]` (`a`-independent). The Schrödinger problem becomes
```
  i χ'(u) = [ D^{(a)}(u) − i W(u) ] χ(u),
  D^{(a)}(u) = diag(E_i^{(a)}(u))  (adiabatic energies; algebraic in u; a-DEPENDENT),
  W(u) = Φ^†(u) Φ'(u)             (anti-Hermitian, off-diagonal; a-INDEPENDENT; geometric).
```
Scattering `S^{(a)} = lim_{u→+∞} U_χ(u,−u)` (phase-regularized at the rank-1 irregular point `u=∞` whose
leading term is `D^{(a)}~u·diag(a)`); `P^{(a)}_{jk}=|S^{(a)}_{jk}|^2`; the open data are `{P_mm, b}`.
**Only `D^{(a)}` carries the `a`-dependence.** Deformation parameter: the slope vector `a` modulo the
2-dim affine gauge `a→αa+β` (so `N−1=2` essential directions; pick a 1-parameter path for the first pass).

**Regime / small parameter.** The BE window actions `δ_{ij}=s_{ij}^2|a_i−a_j|` organize the family:
`δ→0` (well separated, `W` weak vs the diagonal phases) is the natural *initial* end; `δ=O(1)`–large
(deep overlap, the rank-3 core) is the *target* end.

## 2. Strategy

Two ingredients, used together:
- **(a) Variation of holonomy (Duhamel).** Since `W` is `a`-fixed,
  ```
    ∂_a S^{(a)} = −i ∫_{−∞}^{∞} U(+∞,u)\, [∂_a D^{(a)}(u)]\, U(u,−∞)\, du,
    ∂_a D^{(a)}(u) = diag(∂_a E_i^{(a)}(u)),   ∂_a E_i = ⟨φ_i|(∂_a H_0 + u\,∂_a A)|φ_i⟩  (algebraic).
  ```
  The deformation *source* `∂_a D` is purely geometric. The question is whether the right side reduces to
  a function of `S` (and geometry) alone.
- **(b) The integrable zero-curvature (Lax) structure.** Yuzbashyan integrability = *a flat non-abelian
  gauge field in the space of system parameters*: there exists `M_a(u)` with the compatibility
  `∂_a U_u − ∂_u M_a + [U_u, M_a]=0`, `U_u = −i(D^{(a)}−iW)`. The commuting partner supplies `M_a`
  explicitly. This is the standard mechanism by which a parameter-deformation becomes **isomonodromic**
  (Schlesinger/Garnier) and the flow closes. NB WS-D found the Malikis–Cheianov `Ê` for Type-1 is
  **abelian** — so the relevant flat connection is abelian, which (see Milestone 1) is exactly what makes
  the *abelian* data flow closed while leaving the non-abelian seed as the conserved constant.

## 3. Milestones (Milestone 1 is the make-or-break — do it first)

### Milestone 1 — DOES THE `a`-FLOW CLOSE? (and is the seed `σ` conserved?)  [decisive]
Construct `M_a(u)` from the commuting partner; verify the zero-curvature compatibility in `(u,a)`
symbolically/numerically. Then decide among three outcomes:
- **(i) Closes, `σ` conserved (isomonodromic).** The monodromy/Stokes seed `σ` (= `W`'s Stokes data, the
  rank-3 connection constant WS-RH isolated) is **`a`-invariant**; only the abelian connection-normalization
  (`G±`, the phases) flows. Then `{P_mm(a),b(a)}` are `σ` dressed by a *closed, geometric* abelian flow.
  ⇒ the deformation/transport lever is viable (Milestones 2–4).
- **(ii) Closes into a finite-dim ODE for `{P_mm,b}` but `σ` not literally conserved.** Still a usable
  flow (a 2-component ODE in `a` with geometric coefficients). ⇒ Milestones 2–4 with that ODE.
- **(iii) Does NOT close** (the Duhamel integral genuinely needs the full propagator history, no `M_a`
  closing it). ⇒ the principle is *organizing but not computational*; record as a first-class negative and
  stop. **This is the most likely failure mode** and must be tested honestly first.

*Decisive test for Milestone 1:* compute `∂_a{P_mm,b}` two ways — (A) finite-difference along the family
(deformation_family_probe gives the data), (B) the proposed closed-form flow RHS — and check they agree
across the family. Agreement on a path ⇒ the flow closes (numerically-supported); disagreement ⇒ outcome (iii).

### Milestone 2 — derive the explicit flow ODE  [only if M1 = (i)/(ii)]
Write `d/da {P_mm, b} = F(P_mm, b; geometric data(a))` (or the `σ`-conserved + abelian-dressing form).
Identify the special-function class of the flow (expected: a confluent `GL_3`/Garnier / Schlesinger system
— the deformation-side avatar of the WS-RH `σ` / the W₃ frontier). Classify each step controlled vs
uncontrolled; flag the Abelian-ceiling boundary (the abelian part is controlled; `σ` is the conserved seed).

### Milestone 3 — the separated-regime initial condition
At `δ→0` the crossings are independent: `S` factorizes into commuting 2-level (Weber) events, so
`{P_mm,b}` are leading-order **elementary** (independent-crossing / WS-O3 separated limit, RMS ~1% there)
and the exact `σ` is well-conditioned (this is the regime WS-RH/the v-plane solve is best behaved).
Set the IC there. *Caveat (honest):* the exact `σ` at the IC is still the rank-3 transcendental — the
separated regime makes it *well-conditioned and near-elementary*, not elementary. The payoff is
**conditioning/transport**, not closed-form.

### Milestone 4 — integrate to deep overlap and gold-gate
Integrate the flow from the separated IC along a path to canonical and **sampleB (deep overlap)**; compare
`{P_mm,b}` to oracle.py (target ≤1e-6, incl. the deep-overlap regime where WS-O3 and double-precision
v-transport fail). This is the validation *and* the payoff test: does transporting along the closed flow
reach the hard regime that direct methods cannot?

## 4. What success / partial / failure each mean (evidence-ladder honest)
- **Full success (M1=(i), M4 passes):** a closed geometric `a`-flow that transports a well-conditioned IC
  into deep overlap = an independent route to the hard regime + the τ/Schlesinger structure for the named
  closed form. (Would promote the central principle to `analytically-derived` + `numerically-supported`.)
- **Partial (M1=(ii) or M4 only near the IC):** a flow ODE valid in a neighborhood; useful, bounded.
- **Failure (M1=(iii)):** the principle is structural/organizing only — the family cannot compute the
  transcendentals (consistent with the Abelian ceiling, R1). First-class negative; stop and record.

## 5. Risks / honest assessment
1. **Closure (M1) is the gamble.** Generic holonomy variation does NOT close; integrability *might* close
   it via `M_a`, but WS-D's abelian `Ê` warns that the closing structure may only govern the *abelian*
   data, leaving `σ` as a conserved constant that still must be supplied (no free lunch — the Abelian
   ceiling, restated as a flow statement).
2. **Even with closure, the IC is not elementary** — only well-conditioned. The realistic win is a
   conditioning/transport lever, not a closed form. Do not oversell.
3. **The flow's special-function class** is the same `GL_3`/W₃ confluent-Garnier object WS-RH isolated and
   the literature leaves open (Gavrylenko); deriving it in closed form is frontier (rung 3, out of scope —
   the goal here is the *flow*, and numerical transport, not the closed `σ`).
4. **Gauge/parametrization.** Track the `a→αa+β` gauge carefully; the essential deformation is `N−1=2`
   dimensional; do the 1-parameter path first.

## 6. Reusable assets
deformation_family_probe.py (the data + the foundation), ring_structure.py (the commuting Lax/zero-curvature
seed), ws_d_nonabelian_E.md (the abelian `Ê`), gl3_rh_problem.md / gl3_rh_solver.py (the `σ` = `W`-Stokes
isolation, the v-plane engine), oracle.py / num_S12.py (gold gate), skeleton_two_transcendentals.py (the
`{P_mm,b}` target). The whole chain — R1 → skeleton → WS-RH → this — is one statement: *characterize the
holonomy of the single geometric seed `W`.*
