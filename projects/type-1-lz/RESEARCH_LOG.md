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

## Literature findings (workstream-3, 2026-06-01) — confidence flags in memo
Decisive reframing. Key references (metadata cross-checked; full texts NOT opened — fetcher 403):
- **Barik, Bakker, Gritsev, Yuzbashyan, "KZ equations and integrable (hyperbolic) LZ,"
  arXiv:2409.17053, SciPost Phys. 18, 212 (2025)** — solves integrable LZ for **N=2,3,4** via
  the **Gaudin/Knizhnik–Zamolodchikov** structure. MOST RELEVANT: same Type-1/Gaudin machinery,
  does N=3 directly. The transition matrix = **KZ monodromy/connection matrix** of a rational
  flat connection with Cauchy/Gaudin residues at {ε_i} — reducible to {γ,ε,a} + monodromy data.
- **Chernyak & Sinitsyn, time-quadratic commuting operators, arXiv:2006.15144 (JPA 54,115204)**
  — quadratic commuting partner gives **constraints on the N=3 S-matrix** + **asymptotically
  exact (adiabatic-limit) probabilities**. Closest generic-N=3 partial result; baseline to beat.
- **Malikis & Cheianov, arXiv:2505.06048 (PRA 2025)** — **Lax / non-Abelian zero-curvature**
  → algebraic S-matrices; importable recipe for the commuting family.
- **Lin & Sinitsyn, three-state LZ–Coulomb, arXiv:1310.7245 (JPA 47,015301)** — only published
  **exact full 3×3** LZ matrix (hypergeometric), but for a non-generic 1/t Hamiltonian; tells us
  the expected function class + gives a hard analytic check in a limit.
- **Patra & Yuzbashyan 1412.4926; Yuzbashyan 1802.01571** (Type-1/Type-2 classification);
  **Wang & Sun 2209.02888** (no-go rules). **No published uniformly-exact closed form for the
  generic Type-1 N=3 matrix** (~80% confidence) ⇒ the problem is genuinely open.

## Answer to "is Θ^tr high-leverage?" (recorded 2026-06-01)
**No, under the practicality requirement.** Closing Θ^tr (O3) yields one factor (the local
scalar block s_BS) at one ramification point, as a power SERIES in the local Langer coordinate;
extracting an amplitude still needs Borel resummation (the practicality wall), plus
insertion-rule assembly across all active selector points, plus sector-propagator products.
It does not deliver an elementary {γ,ε,a} form. ⇒ DEMOTE O3 from "single highest-leverage."
The literature points to a different top strategy: the **KZ/Gaudin monodromy route**, which
targets a closed amplitude in {γ,ε,a} + a handful of monodromy/period integrals directly. See
strategy memo `paper/type1_lz_strategy_memo.tex`.

## Theory-building: autonomous research program drafted — 2026-06-01
`paper/RESEARCH_PROGRAM.md` (the plan) + `paper/type1_lz_working_paper.tex` (living paper, DRAFT,
not finalized until physics-reflection passes). Unifying principle: P = |Stokes/connection data of a
rank-2 irregular connection on the genus-0 rational spectral curve marked by {eps_i}, irregular inf,
and the exact-crossing node|^2; Abelian ring fixes elementary exponents (BE), non-Abelian connection
coefficient = the confluent-Heun/Kampe de Feriet prefactor. Six workstreams WS-A..F with pass/fail
gates, skill routing, dependency map, coordination protocol (log=shared state; reflection gates;
tournament re-rank per gate). Targets O1 (prefactor closed form), O2 (node reduction), O3
(factorization locus), O4 (non-Abelian E). Awaiting user review/discussion before dispatch.

## WS-F completed + coordinator integration (one finding REJECTED) — 2026-06-01
Deliverables `experiments/oracle.py` (gold oracle wrapping assay ip.propagate_ad_ip + 7-stratum
suite) and `experiments/oracle_report.md`. Verification `experiments/structural_crossing.py`.
- **ACCEPTED:** oracle is faithful (assay builds exactly Type-1 Cauchy H; IP vs lab agree 1e-12).
  BE survivals ≤1e-8 @T=80, ~1e-9 @T>=120; double-stochasticity 1e-11-1e-13. **sampleB reconciled:
  converged P_mid=0.02102, ratio 128.7x** (my standalone 0.017/104x was UNDER-converged; oracle
  trusted). canonical P_mid=0.21472 (2.55x). FIX 2 (Richardson): tails are MIXED-order per entry
  (most 1/T^4, some off-diag ~1/T^3); 8:1 sub-optimal, oracle uses 16:1; persisted
  validation_dataset.pkl (T=120,8:1) carries ~1e-7 not 1e-9 — corrected honestly.
- **ACCEPTED PENDING VERIFY (labeling only):** FIX 1 endpoint-permutation. WS-F claims the assay's
  pi_in=argsort(a),pi_out=argsort(-a) is wrong; correct sheet->diabatic map is fixed PI_IN=(0,1,2),
  PI_OUT=(2,0,1). Affects WHICH transition each P-entry is (not the physics values). Not re-derived
  from ip.py internals here; WS-F's BE-landing test is reasonable evidence. Flag for formula-
  validation downstream.
- **REJECTED:** WS-F report Sec.7 "the node is COMPLEX / real gaps stay >=0.14." REFUTED: the real
  exact crossing is STRUCTURAL and universal — every Type-1 sample has real min gap ~1e-8 (canonical
  -0.249, sampleB -0.689, 8/8 randoms), INCLUDING WS-F's OWN near_node params (crossing at u~-3.02,
  gap 1.9e-8; its bounded search missed it). WS-F's λ-interlacing argument is a non-sequitur (E=m/p:
  two distinct interlacing λ can give the same E). γ-signs are gauge (eigenvalues invariant).
  Confirms gate_test_genus + WS-D + the user's "exact crossing" hint, now shown UNIVERSAL.
- **META (3rd time):** independently check a sub-agent's surprising structural claim before adopting
  — WS-F's complex-node finding would have mis-centered the near-node stratum and the local models.

## Co-scaling derivation (obstruction upgraded) — 2026-06-01
`paper/coscaling_derivation.md` + `experiments/coscaling.py`. Done by coordinator (not delegated;
the obstruction statement needed judgment).
- **WIDTH LEMMA [established, exact, symbolic]:** w_ij = 2|γ_iγ_j|/|ε_i−ε_j| — avoided-crossing
  width is SLOPE-INDEPENDENT (Δa cancels). Primary variable = w_ij; Γ_ij keeps its established
  meaning (BE/LZ exponent γ²γ²|Δa|/Δε²) — the slope-free width vs slope-carrying adiabaticity split.
- **Scale-invariance [analytic]:** both separation and width ~ γ²/Δε (slope scale cancels) ⇒
  sep/width is a dimensionless shape function; no overall scaling grows it.
- **THEOREM (corrected, sharp) [analytic backbone + num bound]:** on the genuine-LZ locus
  (all Γ_ij∈[0.2,5]) the closest crossing pair is within ~1.6 max-widths (35k samples; median 0.20;
  TIGHTER than the unconstrained 4.2). At least TWO crossings are ALWAYS marginally-overlapping ⇒
  MC all-isolated config does not exist in Type-1 ⇒ S=∏S_ij unreachable.
- **Loophole closed:** degenerate-slope limit sends one crossing to ∞ but Γ→0 there (trivial pass);
  the other two stay overlapping. So the earlier worry (ratio unbounded in degenerate limit) does
  NOT provide an escape with genuine crossings.
- Upgrades O4-FAIL obstruction from numerically-supported to analytically-derived backbone. Quotable:
  the Cauchy structure welds widths to level spacing (slope-free), so Type-1 crossings are
  permanently marginally-overlapping — the structural reason it's outside the factorizable class.
- Open: fully-analytic sup of the O(1) constant over shape space.

## WS-A completed (gate-1 PASS) + coordinator integration — 2026-06-01
Deliverable `paper/ws_a_riemann_scheme.md`. Claims structurally forced (H(u) polynomial), low risk.
- **Riemann scheme (established):** honest object = 3×3 diabatic system ψ'=−i(H0+uA)ψ; entire
  coefficients ⇒ EXACTLY ONE singular point: u=∞, rank-2 irregular (the 3-level Weber/parabolic-
  cylinder point). No finite singularities. Accessory parameters (Fuchsian) = 0.
- **Class verdict:** NOT hypergeometric (one irregular point, not 3 regular). The transcendental
  unknowns are PURELY the off-diagonal STOKES connection coefficients of this rank-2 irregular point
  — the 3-level generalization of the Weber/Zener connection problem. (Refines "confluent-Heun":
  Heun needs regular singular points; the linear model is the pure single-irregular confluent case.
  BBGY's Kampé de Fériet was the hyperbolic regular+irregular model.)
- **Γ-dictionary confirmed (corrected):** diagonal exponent ρ_i=−i·b1_i (signed Γ row sum);
  off-diagonal Weber/Stokes exponent = Γ_ij = BE rate. These are residues of the ACTION/Stokes
  one-form, NOT the geometric connection (adiabatic transport ~1/u² ⇒ no simple-pole residue).
  Strengthens E1, no contradiction.
- **CORRECTS my framing (candor):** the exact-crossing NODE is an ORDINARY point of the honest 3×3
  system — it carries NO local connection data and NO accessory parameter. My repeated "exact
  crossing = local reduction handle / rigorous H-B" was WRONG: it conflated the adiabatic-frame
  gauge singularity (eigenvectors rotate fast at the degeneracy) with a local ODE feature. In the
  diabatic frame u* is regular. ⇒ **WS-B as specced (model the node as a local 2×2 connection
  coefficient pinning P₂→₂) is CANCELLED.** The node's 2.5–104× enhancement is a GLOBAL Stokes
  near-degeneracy, not local data.
- **Net for program:** the target is now sharp and singular: compute the off-diagonal STOKES
  CONNECTION COEFFICIENT of the rank-2 irregular point of the 3×3 system. Elementary ceiling = Weber
  data (BE + quadratic/linear phases). This is exactly the project's ORIGINAL exact-WKB / global-
  selector machinery (Voros symbols / virtual turning points). WS-E ← this object. WS-C (factorization
  on sub-loci, e.g. ε₂=½(ε₁+ε₃)) is the remaining gate-2 question. WS-F (oracle) still running.

## WS-D completed + coordinator reflection — 2026-06-01
Deliverable `paper/ws_d_nonabelian_E.md`; verification `experiments/ws_d_verification.py`.
- **KEY CORRECTION (established):** Malikis–Cheianov's `Ê` is **Abelian** ([Ê,H]=0, Ê∈span{I,H,H²},
  machine precision) — it is a Chernyak–Sinitsyn time-quadratic commuting partner re-read as
  ε-translations. There is NO non-commuting generator for either model. ⇒ my earlier framing
  ("the home-run needs a non-commuting Ê outside the ring") was a FALSE PREMISE. Corrected in
  ring_structure.md / RESEARCH_PROGRAM.md narrative (to do at next consolidated integration).
- **WS-D verdict (Type-1 cannot MC-factorize) UPHELD, but its obstruction legs were FLAWED:**
  Coordinator reflection (physics-reflection role) caught it. c1 "rigidity/null-space 0" is WRONG
  (fixing couplings+slopes leaves a 3-dim family; crossings DO move, rank 6→7). c2 "node" is not a
  real blocker (a node = diabatic zero-coupling crossing = trivial 2-level event).
- **REAL obstruction (numerically-supported, 14k-sample scan):** Type-1 crossing separation and
  avoided-crossing width BOTH scale ~γ²/Δε ⇒ ratio bounded ~O(1) (median 0.5, 99%<2.3, max 4.2,
  none>5). Crossings are PERMANENTLY marginally-overlapping ⇒ MC isolated-2-level limit UNREACHABLE
  ⇒ exact S=∏S_ij obstructed (right reason). Also explains the empirical 15% (high-ratio tail) vs
  85% (overlapping). Open: ANALYTIC proof of the sep~width~γ²/Δε co-scaling.
- **Net for program:** O4 = FAIL (home run closed), now for a clean, correct reason. Strengthens
  WS-E: the prefactor is a genuine genus-0 confluent-Heun/Kampé de Fériet connection coefficient
  (BBGY: 3×3 = Kampé de Fériet F^{0:1;1}_{1:0;0}, reduces to ₁F₂ on ε₂=½(ε₁+ε₃) — a concrete
  WS-C elementary-locus candidate). The "marginally-overlapping crossings" fact is itself a clean
  structural theorem of Type-1 worth the paper.
- **META (recurring):** accept sub-agent VERDICTS only after independently checking the load-bearing
  leg; a correct conclusion can ride on a wrong argument.

## Gate-test #1 (genus of spectral curve) — function class CORRECTED — 2026-06-01
Full result: `paper/gate_test_genus.md`. Ran the review's decisive test; it OVERTURNS the review's
own A4 elliptic-refutation (in the constructive direction).
- **Spectral curve Σ: χ_H(E,u)=0 is GENUS-0 RATIONAL, universally for Type-1.** Proof: global Gaudin
  parametrization (E,u)=(m(λ)/p(λ), n(λ)/p(λ)) verified |·|≤9e-16. Riemann-Hurwitz agrees once the
  node is found: D(u)=Disc_E sextic has mults [1,1,1,1,2]; the double root has TRIVIAL monodromy
  [0,1,2] (a NODE, not a branch point), the 4 simple roots are transpositions → smooth 3:1 cover
  with 4 branch pts → g=0.
- **The node = the EXACT CROSSING** (user's hint): real u*≈−0.249 where two eigenvalues are EXACTLY
  degenerate (gap=0). Structural (all samples). Codim-2 degeneracy forced by Type-1 ⇒ crossing-pair
  coupling vanishes at u*. The 4 complex branch pts = avoided crossings at Re(Q4 pairs).
- **Reconciliation:** genus-1 μ²=Q4 is the WKB PHASE/action double-cover (residue-collapsing →
  elementary exponents), NOT the spectral curve. Earlier refutation conflated the two — that step
  was wrong.
- **Corrected class:** NOT elliptic/Painlevé. Genus-0 special-function family. Class set by marked
  points: 3 poles {ε_i} + irregular ∞ + node ⇒ confluent-Heun / Kampé de Fériet (genus-0), not ₂F₁.
- **Honest residual:** genus-0 ≠ elementary (confluent-Heun connection coeffs generically hard);
  but the exact-crossing node has ELEMENTARY (log) local data ⇒ a real reduction handle (rigorous
  form of old H-B). Exponents elementary (BE) confirmed consistent (genus-0 periods = residues).
- **Verdict:** theory-building NOW on correct footing — target genus-0 connection problem on P¹_λ
  with {ε_i}, irregular ∞, node u*; exploit node to cut accessory parameters; expect confluent-Heun/
  Kampé de Fériet prefactor. Gate-test #2 still open (does 3rd-order λ-connection factorize on a
  sub-locus → drop to hypergeometric/elementary there).
- **META:** gate-test caught that BOTH my optimistic (Kampé) and pessimistic (elliptic) priors were
  imprecise — compute the curve's genus+singularities for OUR model, never import from a neighbor.

## Deep-verification review of the KZ line — 2026-06-01
Full review: `paper/REVIEW_kz_line.md`. Peer-review BEFORE theory-building (user-requested).
**SURVIVES (keep, established/analytically-derived):** P = connection/Stokes coefficient of the
rank-2 irregular ODE at u=∞; C-S τ-deformations are isomonodromic (preserve P); commuting-partner
Abelian ceiling; "closed form in special functions, not elementary" at the existence level.
**REVISE/REFUTE (→ conjecture / one sub-claim refuted):**
- A2 "Fuchsian KZ reduction": naive Laplace of i ψ'=(H₀+uA)ψ → Aφ'=(H₀−iλ)φ, coeff LINEAR in λ ⇒
  still IRREGULAR at λ=∞, NOT Fuchsian. KZ identification borrowed from BBGY hyperbolic; doesn't
  transfer. Use "isomonodromy" (true) not "KZ" (unproven here).
- A3 "explicit product-form Euler integral": N=3 Laplace transform is 3rd-order; solution generically
  NOT ∏(λ−ε_i)^{α_i}. Product form only on a degenerate 2nd-order (hypergeometric) locus = the
  Lin-Sinitsyn ε₀=0,g₁₃=0 point already shown OUTSIDE Type-1. Writing the integrand is circular.
- A4 "Kampé de Fériet": **REFUTED as stated.** Phase curve μ²=Q₄ is GENUS-1 ELLIPTIC (4 distinct
  roots verified, both samples). Kampé/₂F₁ live on genus-0. Live class = Heun/confluent-Heun /
  elliptic-theta, or (worst case) a Painlevé-type isomonodromic τ-function. BBGY's Kampé came from
  the rational (genus-0) hyperbolic model; linear/Cauchy Type-1 is genuinely elliptic.
**VERDICT: do NOT build theory on "Kampé de Fériet."** Framework ready; function-class/representation
not. Gate theory-building on: (1) genus+Picard-Fuchs order of the OFF-diagonal governing curve
(sextic gap discriminant), (2) does the N=3 λ-ODE reduce to 2nd order on benchmark strata?,
(3) saddle↔Q₄-window check for any candidate integrand.
**META (recurring error):** committing to a special-function class by analogy to an adjacent solved
model without computing the monodromy/genus of OUR curve (same class as the falsified product
ansatze). Standing checklist item added.

## Ring structure → commuting-partner viability verdict — 2026-06-01
Code: `experiments/ring_structure.{md,py}`. User's lever verified (machine precision):
- Type-1 family {H^(a)(u)} is a COMMUTING RING (max||[H^a,H^b]||=9e-16); spans full 3-dim commutant.
- Time-quadratic partner REDUCES: H(u)²=Σ c_k(u)H^(a_k)(u), c_k LINEAR in u (residual 8e-14).
- Commuting ⇒ ENTIRE family shares one eigenbasis φ_i(u) (overlap defect 3e-16); members differ
  only in eigenvalues E_i^(a)(u), linear in a.
**Viability verdict (negative for naive S2):** a commuting partner carries ONLY Abelian data
(spectrum→WKB/Dykhne→BE survivals; τ-invariance). The quadratic partner adds NOTHING (reducible).
⇒ commuting-partner route has an ABELIAN CEILING and CANNOT reach the middle-survival prefactor —
which the anchor run shows DOMINATES (2.5–104×). This is a proof, not a difficulty. Explains C-S
asymptotic-only, "no Dykhne analog," our prefactor size.
**Constructive flip side:** ring ⇒ P^(a) = |holonomy of FIXED algebraic connection
W_ij=⟨φ_i|φ_j'⟩ (a-independent) twisted by a-LINEAR phases|² = isomonodromy/KZ problem → Kampé de
Fériet, NOT elementary. The non-Abelian holonomy of W IS the missing prefactor.
**Relocates the home run:** S1's zero-curvature Ê must be NON-COMMUTING ([Ê,H]≠0), outside the ring,
NOT subject to this reduction. Bow-tie gets Ê from su(2)/su(3) spin-rep structure; Type-1 is
Cauchy/Gaudin with no manifest spin rep ⇒ existence of a useful non-Abelian Ê is THE decisive open
question. Memo updated to v2.1.

## Day 1–2 anchor experiment (Type-1 vs C-S mismatch resolved) — 2026-06-01
Full write-up + reproducible code: `experiments/anchor_experiment.{md,py}`.
- **Model-mismatch resolved.** Type-1 → C-S Eq.(11) dictionary derived & computed (time-shift t₀
  + diagonal gauge flatten the middle level; ε_CS, b1,b2, couplings read off). KEY: the C-S
  exactly-solvable bow-tie point needs ε_CS=0 AND g_OUTER=0; Type-1's outer-outer coupling
  g_OUTER=γ_loγ_hi(a_lo−a_hi)/(ε_lo−ε_hi) is NEVER zero ⇒ **C-S Eq.49 anchor is OUTSIDE Type-1.**
  (Earlier proposal to reproduce Eq.49 was wrong for Type-1; corrected.) ε_CS=0 IS reachable
  (one condition) but g_OUTER stays nonzero there.
- **Corrected intrinsic experiment** (numpy+scipy, self-calibrated on EXACT BE extreme survivals):
  benchmark P by direct integration vs incoherent-product middle survival.
- **Headline result (numerically-supported, ~1% finite-T):** middle survival is 2.5×–104× LARGER
  than the incoherent product (canonical 0.212 vs 0.084 = 2.5×; sampleB 0.0169 vs 0.00016 = 104×),
  growing with crossing overlap. The interference/prefactor is the DOMINANT part of P_mid, not a
  small correction. Convergence = oscillatory ~1/T Stückelberg tail.
- **Strategic updates:** (i) S2 leading-Dykhne CANNOT reach P_mid (wrong order of magnitude when
  crossings overlap) — confine S2 to BE entries + deep-adiabatic limit. (ii) Re-confirms no
  product-of-real-tunneling-factors closed form (falsified class). (iii) sampleB (104×) is now the
  sharpest discriminator: any candidate must hit P_mid≈0.017. Next: S1 zero-curvature — does the
  path-deformed two-level product reproduce P_mid via Stückelberg PHASES?

## Primary-source read + strategy v2 re-rank — 2026-06-01
Read all four cited papers in full (were search-metadata only before). Memo updated to v2
(`paper/type1_lz_strategy_memo.tex`). Key findings:
- **Chernyak–Sinitsyn 2021 (2006.15144)** = OUR exact model class (linear 3-state, time-quadratic
  commuting partner, Eq.20). Verdict: integrability ⇒ τ-invariant P but "generally NOT expressible
  in known special functions"; for ε₀≠0 "likely no analytical solution." Their tool = Dykhne
  complex-turning-point formula = OUR Q4-window machinery. They state there is **"no general analog
  of the Dykhne formula" for N=3** (only "limited progress" on the prefactor η + subdominant
  exponents) — that prefactor IS our missing P_{2→2}/off-diagonal. Exact handles: ε₀=0 slice
  (confluent hypergeometric, P_{2→2}(0)=2e^{−πg²/b}/(1+e^{−πg²/b}), their Eq.49) + adiabatic-limit
  asymptotics via time-scale separation.
- **Malikis–Cheianov 2025 (2505.06048)** = the elementary route. Zero-curvature operator Ê ⇒
  path-deform in (t,ε) ⇒ exact S = product of 2-level LZ S-matrices (bow-tie: S=S13·S23, Eq.29).
  This is the rigorous "when is the incoherent product EXACT" statement → explains our 15%.
- **BBGY 2024 (2409.17053)** = KZ but HYPERBOLIC (A+B/t); N=3 amplitude = two-variable Kampé de
  Fériet. Tells transcendentality class, not elementary form. **Demoted from v1 #1.**
- **Lin–Sinitsyn 2013 (1310.7245)** = template: LZ-Coulomb N=3 has elementary entries out of the
  special level, hypergeometric between ordinary levels; "simplicity does NOT follow from
  integrability."

**Re-rank (v2):** S1 zero-curvature/Lax path-deformation (Malikis–Cheianov) → TOP (only elementary
output; uses commuting partner directly). S2 adiabatic Dykhne + ε₀=0 exact slice → foundation +
exact anchors. S3 contour/KZ off-diagonal (BBGY+Lin-Sinitsyn) → exact-but-transcendental. S4
KZ/Gaudin direct → demoted (hyperbolic, Kampé de Fériet). S5 transport recurrence → unchanged low.
**Honest reset:** fully-elementary generic (ε₀≠0) closed form likely does NOT exist (~70–80%);
realistic deliverables = BE+elementary entries, exact ε₀=0 slice, adiabatic asymptotics, and
(home run) exact factorisation IFF Ê exists & separates crossings. The pivot question: does the
constant middle level ε₀ obstruct path-separation? (Malikis-Cheianov yes for bow-tie;
Chernyak-Sinitsyn warn ε₀≠0 spoils clean reduction.)

This reconciles the earlier ideation cards: H-B/H-C (exact-crossing/Q4-interference) ARE the
multistate-Dykhne prefactor problem; H-C's "is φ geometric or transcendental?" test now has a
literature-predicted answer (transcendental). H-A slope-flow ⊂ S1 (zero-curvature is its closed form).

## Ideation pass (commuting partners / exact crossing / Q4 pairs) — 2026-06-01
Attack lines generated + reflected + tournamented in `paper/attack_lines.md`. Grounding facts
VERIFIED numerically (standalone, not their harness):
- Commuting family shares p,n,W₄,Q₄,κ_S **exactly** (max|Δ|=0); a enters only via L_H; L_H linear
  in a (residual 3e-17). ⇒ transport curve ν²=W₄ is family-invariant; all slope-dependence is the
  linear L_H in the phase form. [evidence ladder: numerically-supported, machine precision]
- The guaranteed Type-1 crossings are **avoided / phase-active** (min gap nonzero): canonical pair
  01 Δ≈2.3e-4 @ u≈−0.25 (strong avoided crossing), pair 12 Δ≈1.67 @ u≈+0.13; overlap sample
  similar. ⇒ crossings CARRY phase ⇒ H-B/H-C viability test PASSES.
- Q₄ has two conjugate pairs at distinct real parts (two window centers) — consistent with the
  two-window interference picture of H-C.

**Tournament leaders (see attack_lines.md):**
1. **H-B** exact (avoided) crossing → middle survival P_{2→2} via 2×2 local rotation + unitarity.
2. **H-C** avoided crossings at Re(Q₄ pairs) → missing off-diagonal as finite 2-window period sum;
   decisive cheap test = is the interference phase φ a GEOMETRIC period or a transcendental Stokes
   constant. H-B+H-C are complementary (target both missing params, both meet practicality bar).
3. H-A slope-flow/isomonodromy in a (principled backbone; must close on amplitudes not probs).
4–5. H-D KZ/Gaudin, H-E constraint closure (literature-import hedges).

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
- 2026-06-01 — Completed the 4-part program. (1) commuting-family leverage = slope-linearity
  of phase + guaranteed crossing as a regular singular point (but NOT a product ansatz —
  falsified). (2) audit: O1/O2 proved; lone local gap O3 is low-leverage for closed form.
  (3) literature: KZ/Gaudin (BBGY 2409.17053, does N=3) is the live lead. (4) Strategy memo
  written: `paper/type1_lz_strategy_memo.tex`. **Rank: S1 KZ/Gaudin monodromy > S2 constraint
  closure > S3 Lax > S4 contour integral > S5 transport recurrence (finishing only).**
- **Highest-conviction next move:** determine whether the BBGY N=3 KZ solution specializes to
  Type-1 H₀+uA. If yes, the open 2 params are likely KZ monodromy/period data in {γ,ε,a}.
  Verify directly in 2409.17053 + 2006.15144 (full texts not yet opened — fetcher 403).
- Two-week critical path in the memo (S1+S2 parallel days 1–3; candidate days 4–8; falsify on
  8 benchmark strata to 1e-9 days 9–14).
