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

## Open questions
- Q1 (the open problem): Determine the remaining 2 parameters of P (middle-level survival
  P_{1→1} and one independent off-diagonal) for the generic overlapping-crossing family.
- Q2: Pin the window↔level dictionary / prove the residue-bundling rule symbolically.

## Hypotheses
| id | statement | status | evidence / links |
|----|-----------|--------|------------------|
| — | (none yet) | | |

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
