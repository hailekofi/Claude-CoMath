# Research overview — Type-1 N=3 Landau–Zener (meta-review synthesis)

**As of 2026-06-01.** Companion to `type1_lz_working_paper.tex` (full results), `OPEN_PROBLEM.md`
(the open problem + plan), `NOMENCLATURE.md` (conventions), `RESEARCH_LOG.md` (chronology).
This is the roadmap: what is achieved, the boundary of knowledge, and the one capstone step.

## The problem
Closed-form Type-1 N=3 MLZ transition matrix `P` (`H(u)=H₀+uA`, Cauchy coupling
`(H₀)_ij=γ_iγ_j(a_i−a_j)/(ε_i−ε_j)`). Two extreme-slope survivals are Brundobler–Elser (BE), exact.
**Open:** the middle-slope survival `P₂→₂=|𝒮_{mid,mid}|²` (`mid=argsort(a)[1]`) and one off-diagonal.

## Headline
We went from "the middle survival is open and unnamed" to **"the middle survival
`P₂→₂=|𝒮_{mm}|²` (`m=slope-middle`) is `|C_{mm}|²`, built from a specific, named Stokes multiplier
`σ` of the Painlevé-V / confluent-Heun / `c=1` family, with monodromy data fixed *algebraically* by
the exact crossing, computable to `10⁻⁹`."** [Notation: the symbol `S₁₂` used in earlier notes is
DEPRECATED — it conflated `P₂→₂` (diagonal probability), the amplitude `𝒮_{mm}/C_{mm}`, and the
Stokes multiplier `σ`; see NOMENCLATURE.md.] A fully-elementary closed form does
not exist (ruled out, multiple ways); the deliverable is a *named* constant + a *computable* model.

## The theory (why) — established structural results [EST]
| # | result | evidence |
|---|--------|----------|
| R1 | Integrable structure is **Abelian** (commuting ring; shared eigenbasis); the quadratic partner reduces. ⇒ integrability fixes spectrum/BE, **never** the prefactor (the "Abelian ceiling"). | ring_structure.py |
| R2 | Spectral curve `Σ` is **genus-0 rational** (Gaudin parametrization `(E,u)=(m/p,n/p)`), with a structural **node** = a real exact level crossing. | gate_test_genus.md |
| R3 | The exact crossing is **universal** (every Type-1; rational `u_*,E_*`). | structural_crossing.py [NS] |
| R4 | **Slope-free width lemma** `w_ij=|2s_ij|`, `s_ij=γ_iγ_j/(ε_i−ε_j)` ⇒ crossings **permanently marginally overlap** ⇒ no exact factorization; Type-1 is **never tridiagonal**. | coscaling_derivation.md |
| R5 | **Classifier:** elementary ⟺ a level decouples ⟺ joint-free spectral network ⟺ rank-2 boundary. | ws_c_factorization_locus.md [NS] |

## The object (what) — the named connection constant
- `P=|𝒮|²` = Stokes data of a **single rank-2 irregular point** (3-level Weber). [F1, AD]
- Scalar reduction: one **apparent `{0,1,3}` singularity** (= the node, `v_*=E_*`) + the irregular
  point — confluent-Heun class. **Derived 3 independent ways** (WS-E, WS-CH, WS-PA1). [EST]
- **Structure:** `C_{mm} = ` (formal monodromy `c_i=Σ_j s_ij²(a_i−a_j)`, *known* = signed BE) ⋉ two
  Stokes shears in the `{mid,lo}`,`{mid,hi}` carrier spaces (the `12×13` joint). The middle is hard
  because `c_mid` is a **cancelling** sum, not a definite survival. [EST: WS-CH, WS-R]
- **Arguments:** `{two window actions (= BE exponents, scale), scale-invariant cross-ratio χ (shape)}`.
  On fixed shape, `P₂→₂` is a clean 1-D function. [WS-NUM, EST]
- **Named:** Painlevé-V / Lisovyy–Naidiuk confluent-Heun connection coefficient = quasiclassical
  `c=1` Virasoro conformal block = PV τ-ratio (WS-PV and WS-CH converge on this family). [literature]

## The pivot (resolved) — PA-2 = ALGEBRAIC [EST, coordinator-verified]
The accessory parameter (the free knob of a generic Heun connection problem) **is the node** —
the apparent singularity `v_*=E_*`, a *rational* double root of the spectral discriminant. The ODE
is built by rational operations from `H₀,1/a_j`, so the accessory parameter is algebraic in
`{γ,ε,a}`, **not a free transcendental modulus.** ⇒ `C_{mm}` is a *specific* named constant.
**The exact crossing IS the accessory parameter** — the program's arc closes here.

## Computable model [gold]
`P₂→₂(γ,ε,a)` validated vs the `10⁻⁹` oracle across sep/width 0.1→4 (canonical `0.2147`,
sampleB `0.0210`); BE survivals analytic. (num_S12.py, ws_ch.)

## Boundary of knowledge — what remains
1. **Capstone [in progress]:** evaluate the *published* PV/Lisovyy–Naidiuk connection formula at the
   algebraic monodromy data and benchmark vs the oracle. Match ⇒ PV (WS-PA1), `C_{mm}` written
   *explicitly*; no match ⇒ rank-3 (WS-PV), the higher (unpublished) constant. Settles the one open
   structural dispute (PV-via-middle-convolution vs rank-3).
2. **Three proofs to promote `[NS]→established`** (physics-derivation, no new agents): node
   *universality*; the `≈1.6` overlap bound; full accessory-parameter algebraicity (have: 2 samples +
   structural argument).
3. **Simplification hope:** the two shears are rank-2 (Heun) pieces; only their composition is
   rank-3. `C_{mm}` may reduce to (Γ-function dressing of the BE data) × (one joint constant). The
   algebraic accessory parameter is the lever (special monodromy can collapse the conformal block).
   Floor: cannot become elementary (residual curvature + finite coherent floor, WS-NUM).

## The contribution, stated plainly
**Type-1 N=3 is the minimal *full* (non-tridiagonal), confluent MLZ system; its middle survival is
an irreducible confluent-Heun / `c=1` connection constant whose accessory parameter is fixed
algebraically by a universal exact level crossing — and `elementary ⟺ a level decouples` organizes
the entire solvable-MLZ "constellation."** Even if the capstone returns rank-3, this stands.
