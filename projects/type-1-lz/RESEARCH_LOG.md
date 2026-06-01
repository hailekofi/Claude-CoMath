# Research Log: Type 1 Landau–Zener

_Last updated: 2026-06-01_

## Current state (read me first)
Source materials ingested and oriented. This is a mature exact-WKB / spectral-network
program on the **Type-1, N=3 multistate Landau–Zener** model. The 3×3 doubly-stochastic
transition matrix P has **4 free real parameters**; **2 are solved exactly** (the
Brundobler–Elser extreme-level survivals, realized as the imaginary window periods I_X of
the phase one-form on the elliptic curve E: µ²=Q₄). The **open problem** is the remaining
**2 parameters** (middle-level survival + one independent off-diagonal), which are provably
*not* periods of the phase curve Σ_Y. Next: generate candidate attack lines (`physics-
intuition`), reflect, and tournament them. Awaiting user steer on which sub-line to open.

## Definitions (from uploads/)
- **Model**: H(u)=H₀+uA, A=diag(aᵢ), (H₀)ᵢⱼ=γᵢγⱼ(aᵢ−aⱼ)/(εᵢ−εⱼ) — the Cauchy off-diagonal
  is the **Type-1** property; Type-1 matrices form a **commuting family** (integrability).
  [LZ_summary.tex gives the augmented-ODE solver in the gauged adiabatic basis; 2N−1 ODEs.]
- **Observable**: P_{x→j}=|f_j(+∞)|², the N×N transition matrix (rows sum to 1).
- **Γ_ij = γᵢ²γⱼ²|aᵢ−aⱼ|/(εᵢ−εⱼ)²** (elementary, no elliptic integral).
- **Sextic gap polynomial** (Type1 Note): Δ⁶−2A(t)Δ⁴+A(t)²Δ²−D(t)=0, A=½(3 tr H²−(tr H)²),
  D=Disc χ_H — universal for any 3×3 pencil; explains the 16/28 monomial sparsity. Pairwise
  WKB lives on this explicit rational cover; **virtual turning points / global selector**
  remain the hard part.

## Key established results (per internal notes — provenance: uploads/, not yet re-derived here)
- BE extreme-level survivals = exp(−|Im I_X|), validated to ~1e-9 (Richardson-extrapolated).
- I_X is purely imaginary (Re I_X=0): residue collapse on E + conjugate-pair symmetry.
- Topology: g(E)=1 ⇒ b₁=2 ⇒ exactly two window periods ↔ two extreme levels. Middle level
  has no cycle on the genus-1 phase curve ⇒ "the open part is the part Σ_Y cannot see."
- Falsified ansätze (do NOT retry as stated): naive 2⊕1 embedded product; incoherent
  stochastic grid (only valid ~15% of param space, well-separated crossings); coherent
  three-rotation amplitude product; algebraic frame maps (any working frame map must be
  transcendental).

## HARD REQUIREMENT (user, 2026-06-01) — computational practicality
The closed form must be **reducible to the core parameters {γ, ε, a}** and must **avoid
numerical integration entirely**, or use **at most a handful of Abelian / geometric integrals
that derive directly from the established geometry** (window/transport periods on the curves
E:µ²=Q₄, ν²=W₄, the mixed cover). Exact-WKB / Stokes data has NOT yielded a computationally
practical implementation so far — so a strategy's rank is gated by how elementary/finite its
output is in {γ,ε,a}, not just whether it is "in principle" determined. Prefer: residues,
algebraic functions, finite period sums. Avoid: order-by-order Borel sums requiring numerical
resummation, ODE matrizants, dense Stokes-graph integration.

## Goal (refined 2026-06-01)
**Closed-form solution for the Type-1 N=3 transition amplitudes/probabilities P.** Four-part
program requested by the user:
1. **Leverage** the commuting-partner structure + the *guaranteed exact crossing*.
2. **Exhaustive assessment** of what's already achieved (proved vs validated vs target).
3. **Literature review** of ideas/innovations in multistate LZ theory.
4. **Rank-ordered recommendation** of strategies + near-term targets toward closed form.

Workstreams: (2) deep code+notes audit — RUNNING (Explore agent); (3) literature — RUNNING
(general-purpose agent); (1) commuting-triple/crossing leverage — driven directly; (4)
synthesis after 1–3 land.

## Updated open-problem structure (per PRL draft, dated 2026-04-18 — newer than the BE brief)
The local selector reduction has advanced. **Theorem-level (proved in the canonical exact
two-sheet reduction at a simple transport ramification point):**
- T1: exact local gauge is **torus-valued to all orders** (Ĝ_v=exp(Θ_v σ₃), scalar exponent).
- T2: odd phase coeffs D_{2m+1}(H) obey a **closed Lagrange–Bürmann recurrence** (Supp eq 12–18).
- T3: full member-dependent phase series is **affine-linear in L_H(v), L′_H(v)** across the
  commuting family: D_{2m+1}=A_m(v)L_H(v)+B_m(v)L′_H(v) (geometry-only A_m,B_m).
- T4: same-radius local selector multiplier is **scalar**: s_BS = −i·exp(2Θ_v(t;H)) (eq 29).

**Still open (sharply concentrated):**
- O3 **Transport recurrence**: all-orders recurrence for the *transport* exponent Θ^tr_v(t)
  — the analog of T2/T3 on the transport side (NAMED highest-leverage symbolic step).
- O4 **Nongeneric collisions**: degenerate transport ramification; common zeros of Q₄ and W₄.

**UPDATE (workstream-2 audit, 2026-06-01): O1 and O2 are now PROVED (generic simple regime).**
The PRL §8 phrasing lagged the notes. Per `selector_sufficiency_outer_descent_memorial.tex`
(Thm 1, 1302/1302 @60-digit, wall residual ≤1.3e-8) and `exact_global_insertion_rule_note.tex`
(Thm 2.1 + overlap-independence Lemma 2.2): selector sufficiency and the exact global insertion
rule are proved for the generic simple case. **H1/H2 are superseded** by a *bridge lemma*
(Thm 1.1, memorial): in matched half-form normalization the adiabatic and resolved-sheet frames
differ by a constant geometry-only G_v=I — no hypothesis needed. So H1/H2 are no longer the
frontier.

**⇒ The remaining symbolic frontier is essentially ONE object: O3, the transport exponent
recurrence Θ^tr_v(t).** Phase side is closed to all orders (D_{2m+1} affine-linear in L_H,L′_H);
transport side has only K₁=−A₂/(2A₁) and K₃ explicit (in `transport_phase_jump_note_augmented.tex`,
geometry-only, from Σ^(2) derivatives along inverse branches). No code module computes Θ^tr yet.
Once Θ^tr_v is closed, the full torus exponent Θ_v=Θ^tr_v+Θ^ph_v and the selector multiplier
s_BS=−i·exp(2Θ_v) become algorithmic closed-forms to all orders, locally.

## Hypotheses
| id | statement | status | evidence / links |
|----|-----------|--------|------------------|
| H1 | After spectator decoupling + curvature subtraction, no nondecaying torus-breaking off-diagonal survives in outer matcher | numerically-supported | 1200/1200 resolved-sheet probe; 1001/1001 torus-survival (Supp Tab.1) |
| H2 | Residual constant non-torus part is sector-independent + geometry-only (absorbable) | numerically-supported | 999/1001 base, 2 borderline pass 100-digit retest (Supp §7) |
| O3 | Θ^tr_v admits a closed all-orders recurrence analogous to D_{2m+1} | conjecture | not yet attempted in notes (per PRL §8) — TARGET |
| L1 | **The entire Type-1 commuting family shares ONE coupling matrix S_ij(u); the slope vector a enters the gauged ODE only through the diagonal phase diag(E_i−E_x), which is LINEAR in a.** So P(a) is generated by a fixed geometric transport S + an a-linear phase. | analytically-derived (from LZ_summary.tex) | λ_i, Σ²_i, Γ_i, S_ij depend only on (γ,ε) [eqs (1)-(5)]; E_i=Σ a_j γ_j²/(λ_i−ε_j) linear in a [eq (E)]. Coupling = "transport" (η=S du), phase = E. NOT exploited by integrability.py. |
| L2 | Because dD_a/da is explicit and S is a-independent, P(a) obeys first-order matrix ODEs in slope-space a, anchored by the exactly-known BE extreme entries — candidate route to pin the middle-level/off-diagonal. | conjecture | follows from L1 via variation-of-parameters on the matrizant M(+∞;a); UNTESTED — workstream (1) target |

**Reflection on L1 (cross-check):** L1 is *corroborated* by the proven T3 — the phase's
a-linearity is exactly why D_{2m+1}=A_m L_H+B_m L′_H is affine-linear (L_H is the
commuting-family member's linear factor, ∝ slope data). So my independent leverage reasoning
reproduces an established theorem — good. **Caution (from integrability.py audit):** the *naive*
commuting-partner exploitations (two-Dykhne strong product; BE-extended min/max/sum) are already
FALSIFIED (residual ~1e-1; inversion gives non-symmetric implied Γ′). So L2 must NOT collapse to a
product Ansatz; its only viable form is a genuine differential/transport statement in slope-space,
subordinate to the assembly picture below.

## Path to closed-form P (assembly view, post-audit)
P is assembled as: **(known BE extreme diagonals) × (selector/insertion product)**. The
insertion product is now PROVED closed (O1,O2) *except* it contains the local scalar block
s_BS=−i·exp(2Θ_v), and Θ_v=Θ^ph_v+Θ^tr_v. Θ^ph is closed; **Θ^tr (O3) is the one missing
piece.** ⇒ Closing O3 ⇒ closed-form local block ⇒ closed-form P (generic simple regime),
with nongeneric collisions (O4) as a separate finish. This makes **O3 the single highest-
leverage target.**

## Evidence ledger
- (empty)

## Dead ends (do not retry as stated)
- (none yet)

## Dead ends (do not retry as stated)
- 2⊕1 embedded product ansatz — forces structural zeros not present in generic P.
- Incoherent stochastic LZ grid — exact only in well-separated limit (~15% of param space).
- Coherent three-rotation amplitude product — fails for overlapping crossings.
- Algebraic frame maps — ruled out; any working frame map must be transcendental.
- (point-local) Airy collapse of projected Q₄ windows — falsified per graphify graph.

## Decisions & next steps
- 2026-06-01 — Created `projects/type-1-lz/`, ingested uploads, oriented on the program.
  Next: read prl_journal_main_text + supplement and the assay/ package to see exactly what
  is coded/validated, then run `physics-intuition` to generate attack lines for Q1/Q2 and
  tournament them. Awaiting user steer on sub-line + whether to prioritize Q1 (full P) or
  Q2 (prove the BE residue-bundling rule).
