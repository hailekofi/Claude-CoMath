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
- **BE pairwise exponent = s_ij²·|aᵢ−aⱼ|**, with **s_ij = γᵢγⱼ/(εᵢ−εⱼ)** (signed Cauchy half-width)
  and **w_ij = |2 s_ij|** (elementary). [NOMENCLATURE: Γ is RESERVED for the canonical Cauchy form
  factor Γ_j and S_ij=Γ_iΓ_j/(λ_i−λ_j); scattering matrix = 𝒮. See `NOMENCLATURE.md`.]
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

## Ordering-convention audit (user concern) — results SAFE, convention locked — 2026-06-01
User flagged: basis ordering should be eps, not energy; worried about inconsistency imperiling results.
DECISIVE TEST + AUDIT:
- Three orderings: eps (basis index, fixed), slope a (BE extreme/middle), energy (reverses at +-inf).
- The OPEN middle survival P_{m→m} = P[m,m], m=argsort(a)[1] = the MIDDLE-SLOPE level (BE physics).
  slope-middle = eps-middle ONLY when slopes monotonic in eps. Non-monotonic test (a=(2,-1,0.5),
  eps=(-2,0,3)): slope-mid=index2 (P=0.571, open); eps-mid=index1 is the slope-MIN (P=0.0185=EXACT BE)
  -> using eps-middle would return a BE-exact value, NOT the open quantity. Convention matters.
- AUDIT of all middle-survival code: anchor_experiment, num_S12 (WS-NUM), restart_probe (WS-R),
  ws_ch (WS-CH) ALL use argsort(a) (slope) consistently. Geometric scans (coscaling, structural_
  crossing, ws_d_verification) have no middle-label dependence. => NO inconsistency; results NOT
  imperiled.
- Headline samples canonical & sampleB have a ASCENDING in eps => slope-mid=eps-mid=index1 =>
  their specific numbers (P_mid=0.2147, 0.0210, all benchmarks) unconditionally SAFE.
- Convention LOCKED in NOMENCLATURE.md: basis=eps-index; P_{m→m}=P[argsort(a)[1],argsort(a)[1]]
  (slope-middle, NOT eps-index-1 in general); energy ordering is NOT a labeling convention.
  Retracted the earlier 'P_{m→m} by energy rank' framing as a red herring.

## STRATEGIC REASSESSMENT + traceable synthesis written — 2026-06-01
NEW companion artifacts for cold pickup: **SESSION_SYNTHESIS.md** (full traceable arc + hypothesis
ledger + dead ends with WHY + open questions) and **META_REVIEW.md** (recurring failure patterns,
propagate-forward rules). Read those + NOMENCLATURE.md first.
- **Honest state:** NO concise computable closed form yet; STRONG multi-method evidence that an exact
  one in elementary OR *published* special functions does NOT exist. We have: structural theory (why
  hard), the object characterized as a RANK-3 connection coefficient with accessory parameter fixed
  algebraically by the exact crossing (node), and a validated NUMERICAL model (1e-7..1e-9).
- **Q2 (user) — node <-> connection coefficient, status:** the relationship "exact crossing (node
  E_*) = the apparent/accessory singularity of the connection ODE" is VALID and verified on 4/4
  independent samples (incl. non-monotonic; clean rational match e.g. -5313773/12376000), with a
  clean characterization (unique rational root of the cyclic Wronskian det[e0,Me0,M^2e0],
  M=diag(1/a)(H0-vI)). BUT: (a) NOT a general rigorous theorem (4 samples + structural argument);
  (b) 'confluent-Heun' was the WRONG label -> it is RANK-3 (capstone); (c) it is STRUCTURAL, not
  CONSTRUCTIVE -- it locates the accessory parameter, does not give a computable expression.
- **Diagnosed strategic drift:** pursued the CLASSIFICATORY path (name the special function:
  Kampe de Feriet -> Painleve V -> confluent-Heun -> c=1), which dead-ended at 'unpublished rank-3'
  and by construction never yields a formula; ABANDONED the project's original CONSTRUCTIVE exact-WKB
  program (selector/global-insertion/virtual-turning-points = the 12x13 joint). Likely OVER-concluded
  'no concise form': rank-3 != 'no concise constructive formula' (the product of two explicit Weber
  shears was never tried).
- **Recommended redirect (pending user steer):** recalibrate goal to constructive/realistic targets:
  O1 [primary] constructive exact-WKB product P_{m→m}=|M_form . Weber . Weber|^2_{mm} with explicit
  Weber shears + the spectral-network junction factor; O2 exact integral representation; O3 validated
  semi-analytic approximation. Proofs owed: node universality, node=accessory (general, via cyclic
  Wronskian), the ~1.6 bound, accessory algebraicity, the dim-count 6>2.

## WS-CAP completed — CAPSTONE VERDICT: RANK-3 (WS-PV confirmed; WS-PA1 'PV' refuted) — 2026-06-01
Deliverables paper/cap_connection_formula.md, experiments/cap_connection_formula.py. The capstone
DECISIVELY settles the one open structural dispute. (Reverses my prior lean toward WS-PA1's 'PV'.)
- **VERDICT: the published rank-2 Painleve-V / Lisovyy connection constant does NOT compute the
  generic Type-1 N=3 P_{m→m}.** It fails at the SIZE/RANK of the monodromy manifold, not a missing calc.
- T1 [analytic, lit-grounded]: the 3x3 rank-2 irregular point has a 6-DIM wild character variety vs
  PV's 2-DIM (PVI cubic surface -> PV by confluence). 6>2: the rank-2 formula has too few parameters.
- T2 [analytic]: closes WS-PA1's 'PV via Harnad/middle-convolution' escape -- the three GENERIC
  distinct irregular rates {a_i} (no {+L,-L,0} symmetry) can't be the single-middle-convolution image
  of a 2x2 PV's two rates. Directly refutes PA-1.
- T3 [exact]: formal exponents = signed-BE traceless SL(3) PAIR c_i (Sum=0), not PV's single theta_inf.
- T4 [exact]: accessory v_*=E_* rational (PA-2 favourable, CONFIRMED) -> P_{m→m} is a SPECIFIC named
  constant, but a RANK-3 one.
- **T5 [gold-gated, DECISIVE]:** on 7/14 strata P_{m→m} lies STRICTLY OUTSIDE the widest single-sigma
  (rank-2/PV) Stuckelberg band. sampleB P_{m→m}=0.02101769 (gold, err 1.5e-8; matches our oracle
  0.021018) vs band [0.8327,0.8800], margin +0.81. No single rank-2 connection constant reproduces it
  -> a third coherent amplitude (the two shears' non-commutative composition through the shared mid
  sheet) is required = genuine RANK-3.
- T6 [gold-gated, positive control]: decoupling locus (gamma_hi->0) -> 2x2, P_{m→m}->0.221360
  (elementary), band tightens around it -> PV applies EXACTLY at the elementary corner only.
- **Corroboration (coordinator) -- HONEST: INCONCLUSIVE.** A quick re-check sought matched-window /
  different-chi pairs to show P_{m→m} needs the shape chi beyond the formal exponents (the >rank-2
  signature). The random search only found tiny-Delta-chi pairs (Dchi=0.05, dP_mid=0.003, comparable
  to the residual window mismatch Dw=0.009) -> UNDERPOWERED, does NOT cleanly corroborate. A proper
  test needs a constrained deformation (hold both window actions fixed, vary chi). WS-NUM's separation
  finding (P_{m→m} = f(window actions, chi)) remains the standing chi-dependence evidence but is itself
  numerically-supported. => The RANK-3 VERDICT RESTS ON T1 (dim 6>2, analytic), T2 (middle-conv
  obstruction, analytic), and T5 (gold-gated band falsification, the airtight decisive tier), NOT on
  this re-check.
- **EXPLICIT NAMED OBJECT:** P_{m→m} = |C_{mm}|^2 of the rank-3 isomonodromy problem (one Poincare-
  rank-2 irregular point with known formal exponents c_i; one apparent {0,1,3} point at the algebraic
  node v_*=E_*) -- a c=1-family / confluent-Garnier (9/2-type) connection constant, ABOVE 2F1 and
  ABOVE the rank-2 PV transcendent, UNPUBLISHED in closed form; reduces to the published PV/Barnes-G
  constant ONLY on the decoupling (elementary) locus.
- **Simplification hope (partial):** the two shears are individually 2x2 (Weber/Gamma), but their
  composition is genuinely non-commutative through the shared mid sheet; the joint constant does NOT
  collapse to a single PV/Barnes-G value on the overlapping stratum.
- Honesty (agent): T1/T2 are [analytic] (standard dim count + leading-rate argument; lit-grounded, not
  machine-proved); the DECISIVE tier is T5's gold-gated falsification (airtight). PA-2 (algebraic)
  STANDS; only PA-1's 'PV type' overturned -> RANK-3.
- **NET (deliverable achieved):** the open quantity is a SPECIFIC, NAMED, algebraically-pinned,
  computable RANK-3 (confluent-Garnier/c=1) connection constant -- a genuinely NEW object, not the
  published PV. Even sharper than 'it's published PV': Type-1 N=3 realizes a minimal rank-3 connection
  constant. The governing theory + computable model stand.

## WS-PA1 completed — PA-0/1/2 — PA-2 VERDICT: ALGEBRAIC (the make-or-break pivot) — 2026-06-01
Deliverable paper/pa1_isomonodromy_foundation.md; coordinator verification experiments/pa2_accessory_algebraic.py.
- **PA-0 [established; R8 TODO CLOSED — 3rd independent confirmation]:** cyclic-vector elimination
  (two different cyclic vectors, both samples) => 3rd-order scalar Laplace ODE with ONE finite
  singular point v_* (NOT an H0 eigenvalue), exponents {0,1,3}, NO log (Frobenius resonances auto-
  satisfied + numerical loop monodromy ||M-I||~1e-13), + rank-2 irregular at v=inf. v_*=-748/375 (S1),
  -112147/44000 (S2). (WS-E, WS-CH, WS-PA1 all agree on {0,1,3}.)
- **PA-1 [analytically-argued]:** Painleve type = PAINLEVE V. The 3x3 is the HARNAD/middle-convolution
  image of the 2x2 PV system (rank-2 irregular u-frame <=> unramified rank-1 irregular at v=inf +
  one apparent point); the {0,1,3} apparent point is the middle-convolution signature; matches the
  restart 2D {mid,.} reduction. **This OVERRULES WS-PV's 'not PV / rank-3' inference:** the linear
  system IS rank-3 (WS-PV correct) but its ISOMONODROMY is PV via middle convolution (rank-3 systems
  can have PV isomonodromy). => the PUBLISHED PV connection constants (Its-Lisovyy-Prokhorov;
  Lisovyy-Naidiuk confluent-Heun) apply. [needs the decisive formula-evaluation test to confirm.]
- **NODE DICTIONARY [coordinator-VERIFIED, both samples]:** the accessory point v_*=E_* = the doubly-
  degenerate eigenvalue of H(u) at the UNIVERSAL REAL NODE u_* (THE EXACT CROSSING, the user's hint!).
  Exact rational: canonical v_*=-748/375 at u_*=-187/750; sampleB -112147/44000 at -3031/4400.
  THE EXACT CROSSING = THE ACCESSORY POINT. The whole arc closes.
- **PA-2 VERDICT: ALGEBRAIC [coordinator-VERIFIED — the realistic-win branch]:** the accessory point
  is the rational node eigenvalue; the ODE is built by rational ops from H0 and 1/a_j, so ALL local
  data at v_* (the accessory parameter) is rational/algebraic in {gamma,eps,a}. The accessory
  parameter is NOT a free transcendental modulus. => S_12 = a PV/c=1 connection coefficient with
  ALGEBRAICALLY-FIXED monodromy data = a SPECIFIC NAMED constant (computable; not elementary, per
  WS-E/WS-NUM benchmarks). Honest caveat (WS-PA1): algebraic data pins a specific named constant; it
  does not make that constant elementary.
- **ALL 5 WORKSTREAMS DONE.** Deliverable (named form + computable model) ACHIEVED. The ONE remaining
  decisive step: evaluate the published PV/Lisovyy-Naidiuk connection formula at the algebraic
  monodromy data (c_i exponents + the rational accessory parameter/node + the two {mid,.} shears) and
  benchmark vs the oracle -> (a) confirms PV (WS-PA1) vs rank-3 (WS-PV), (b) writes S_12 EXPLICITLY.

## WS-NUM completed (PA-4) — computable P_{m→m} floor MET + clean negative + separation finding — 2026-06-01
Deliverables experiments/num_S12.py, paper/num_S12_model.md (num_S12_dataset.pkl gitignored).
- **FLOOR DELIVERABLE MET [gold]:** a trusted computable P_{m→m}(gamma,eps,a) (one adiabatic-IP pass,
  T=80, rtol=1e-9, ~14s) validated vs the oracle across sep/width 0.1->4: canonical 1.6e-9, weak
  1.1e-8, strong 7.8e-8, sampleB 2.1e-7, well_sep 5.8e-7; BE extreme survivals reproduce analytically.
- **PARAMETRIZATION [established]:** the two window actions I_X = the two BE extreme-survival
  exponents Sigma_lo, Sigma_hi (machine precision); the middle level is the shared partner of BOTH
  windows (why its survival is the open COUPLED quantity).
- **NEW [coordinator-VERIFIED]:** the Q4 cross-ratio chi is SCALE-INVARIANT (a pure SHAPE coordinate;
  verified |Delta|=0 under gamma-scaling, 4e-16 under eps-scaling). P_{m→m} SEPARATES as
  (scale -> BE exponents) x (shape -> chi). On a fixed-shape ray (chi frozen) P_{m→m} is a clean 1-D
  function of the middle BE exponent (1-D slice model reproduces gold to median 2e-3). The natural
  arguments of S_12 are thus {the two window actions, chi} -- exactly the OPEN_PROBLEM target args.
- **RECOGNITION [clean NEGATIVE]:** no elementary closed form. DO/coherent-path candidates fail
  sample-INDEPENDENTLY; single-sample PSLQ relations spurious. logit(P_{m→m}) only approx affine in
  log(b_mid) with residual curvature = fingerprint of a transcendental connection coefficient. Deep-
  adiabatic strata show a FINITE COHERENT FLOOR (P_{m→m}~0.01-0.03 while incoherent products ->0) =
  the genuinely open transcendental content. Corroborates the confluent-Heun/conformal-block verdict.
- WS-NUM provides the GOLD TARGETS for WS-PV/WS-CH to PSLQ named special-function constants against.
- 4 of 5 workstreams done (WS-PV, WS-CH, WS-R, WS-NUM, all mutually consistent). Remaining: WS-PA1
  (PA-0/1/2: rank + the PA-2 accessory-parameter-algebraicity pivot).

## WS-R completed — restart-structure tests (H-R3,H-R2 confirmed; H-R1 reframed) — 2026-06-01
Deliverables paper/restart_structure.md, experiments/restart_probe.py. All measured in the
convergent canonical frame (Coulomb-subtracted; T-convergence verified). Coordinator spot-checks pass.
- **H-R3 CONFIRMED [symbolic]:** -u'(lambda)=W4(lambda)/p(lambda)^2 (W4=np'-n'p), so Gamma_j~(lambda-eps_i)/|g_i|
  near the poles -- ANALYTIC, no monodromy at eps_i. The sqrt branch points of Gamma are the W4 zeros
  (turning points, -u'=0), 4 complex points well separated from the real poles. Turning-point
  (window/BE) branching and irregular-point (off-diagonal Stokes) structure are on distinct sheets.
- **H-R2 CONFIRMED [oracle, strong]:** off-diagonal weight concentrates on the mid row/col
  (frac_mid 0.607-0.757); shrinking ONE outer coupling (gamma x1e-3) collapses P_mid to the
  incoherent BE product (|diff|~1e-3) AND frac_mid->1.0 -- one shear vanishes, joint -> single Weber
  factor. The {mid,lo},{mid,hi} shears (=12x13 joint) ARE the non-factorizing carrier;
  'elementary <=> a level decouples' reconfirmed.
- **H-R1 REFRAMED (my hypothesis PARTLY WRONG):** (a) the off-diagonal is CYCLIC (dominated by the
  3-cycle lo->mid->hi->lo = PI_OUT=(2,0,1)), NOT triangular (coordinator-confirmed from |S|^2 argmax
  pattern (0<-1,1<-2,2<-0)). (b) 'e^{2pi i c_i} = a diagonal factor of S_canon' is a CATEGORY ERROR:
  c_i is the logT drift SUBTRACTED to build the frame, so S_canon carries no such factor; extracting
  'D.P' reproduces the content-free Q_restart. The sound residue is exactly H-R2.
- **RECONCILIATION with WS-CH:** WS-CH's C=M_form . shear . shear is the CENTRAL CONNECTION MATRIX
  (formal Thome basis at inf <-> Frobenius at the apparent point) -- there M_form (formal monodromy
  e^{2pi i c_i}) IS a factor (the Thome basis carries the formal exponents). WS-R's statement is about
  the REGULARIZED SCATTERING matrix S_canon (diabatic basis), where M_form is subtracted. Different
  matrices/frames; both correct. INVARIANT content (both agree): the TWO mid-pair Stokes shears.
- **CAUTION (carry into PA work):** use the connection-matrix framing (M_form (x) two mid-pair shears,
  WS-CH); do NOT claim 'triangular shear' or 'e^{2pi i c_i} as a factor of S_canon'. The frame-pinned
  S_canon reality is cyclic-permutation x mid-pair shears.

## WS-CH completed — confluent-Heun central connection; benchmarked computable P_{m→m} — 2026-06-01
Deliverable paper/ch_direct_connection.md; benchmark scripts archived experiments/ws_ch/.
- **CONFIRMS (2nd independent derivation, coordinator-spot-checked):** scalar 3rd-order Laplace ODE
  has ONE apparent singularity {0,1,3} (no log), residue A=1, v_* RATIONAL (canonical -748/375,
  sampleB -112147/44000), NOT an H0 eigenvalue. This independently confirms WS-E's {0,1,3} =>
  closes the paper's R8 TODO (two agents agree).
- **CONFIRMS the restart structure:** connection matrix C = M_form . S^{(mid,hi)} . S^{(mid,lo)}
  (formal-monodromy diagonal (x) two Stokes shears in {mid,lo},{mid,hi}); S_12=|C|^2_{mid,mid}.
- **Formal monodromy = signed-BE [coordinator-VERIFIED]:** c_j=sum_k s_jk^2(a_j-a_k); |e^{2pi c_j}|
  reproduces the BE survival for the EXTREME levels (canonical e^{2pi c_lo}=0.07474=P_lo; sampleB
  0.00040). c_mid is a MIXED-SIGN CANCELLING sum (not a BE survival) => clean structural reason the
  MIDDLE is the hard one. [new solid result]
- **BENCHMARK PASSED <=4e-7** (target 1e-4): connection decomposition via 2 independent solvers
  reproduces the oracle (canonical 0.214724 diff 3.95e-7; sampleB 0.021018 diff 1.56e-7), doubly
  stochastic <1e-5. A computable, validated P_{m→m}.
- **Named object:** Lisovyy-Naidiuk (2208.01604) confluent-Heun connection = convergent series in the
  accessory parameter = quasiclassical Virasoro conformal block = Painleve-V tau-ratio (Bonelli-
  Iossa-Panea Lichtig-Tanzini). SAME family WS-PV targets -> the two tracks CONVERGE.
- **RANK DISCREPANCY (WS-CH confluent-Heun rank-2 vs WS-PV rank-3) -- RECONCILED [coordinator synth]:**
  the two Stokes SHEARS (the {mid,lo},{mid,hi} carrier-space connections) are each rank-2 confluent-
  Heun connection coefficients (Lisovyy-Naidiuk applies to EACH, published); S_12=P_{m→m} is their
  NON-COMMUTATIVE COMPOSITION = the rank-3 (12x13 joint) object (WS-PV). Pieces rank-2 (published);
  composition rank-3. To be confirmed by WS-PA1's PA-1 (irreducibility / accessory-param count).
- PA-2 still the pivot: v_* rational + A=1 is SUGGESTIVE of algebraic accessory parameter but WS-CH
  did not prove full algebraicity -> WS-PA1.

## WS-PV completed — RANK CORRECTION (target is rank-3, not Painleve V) — 2026-06-01
Deliverable paper/pv_tau_route.md. Coordinator-verified the load-bearing logic.
- **CORRECTION to my plan (and to the WS-A/WS-E 'confluent-Heun' label):** Painleve V is a 2x2
  (rank-2) isomonodromy problem and confluent-Heun a 2nd-order (rank-2) ODE; the Type-1 scalar
  reduction is genuinely RANK-3 (3rd-order, irreducible over Q(g,e,a) per WS-C; WS-E 3rd-order
  Laplace ODE). So PV/CH are the rank-2 DECOUPLING-BOUNDARY cases only: sending one coupling
  s_ij->0 block-reduces the 3x3 Laplace system K(v)=-i diag(1/a)(H0-vI) to 2x2 (+1x1) = the
  elementary boundary (= WS-C 'elementary <=> a level decouples'). [verified symbolically]
- **Generic S_12 = one rank up:** a rank-3 confluent-Garnier / c=1 irregular-conformal-block
  (Barnes-G family) connection constant — UNPUBLISHED, not classical. Published PV/PVI connection
  constants (Its-Lisovyy-Prokhorov 1806.08344; Lisovyy et al.) apply at the boundary only.
- WS-PV delivered the NAMED SHAPE (tau-function ratio / Barnes-G connection constant), the monodromy
  data in {g,e,a} (formal exponents c_i KNOWN + the two {mid,lo},{mid,hi} Stokes shears), and made
  the pivot explicit: the one accessory coordinate sigma is exactly what PA-2 fixes (algebraic =>
  closed; transcendental => irreducibility result).
- CONTINGENT on WS-PA1 PA-0 (confirm 3rd-order irreducible / {0,1,3} apparent point not removable);
  if WS-PA1 finds it reduces to 2nd order, PV/CH WOULD apply. WS-CH (attacking confluent-Heun=rank-2)
  may hit the same rank mismatch -> will reconcile on report.
- Lit claim 'rank-3 deformation lacks Painleve property' (2503.22198, 2512.24083) is literature-
  sourced, flagged, NOT load-bearing for the main reframe (which follows from irreducible 3rd-order).
- Honest impact: the user-favored PV lever is the BOUNDARY case; the realistic named target is one
  rank higher. Deliverable shape unchanged (named connection constant + WS-NUM computable model);
  PA-2 still the make-or-break.

## #2 implementation launched — restart-operator structure + 5 parallel workstreams — 2026-06-01
User approved the open-problem plan with steers: Painleve-V/tau route + a PARALLEL direct
confluent-Heun connection-coefficient track; PA-2 (Gaudin fixes accessory param algebraically) is the
make-or-break; deliverable = named closed form (tau-function/connection coeff) + computable numerical
model; fold in restart-operator insights; fire H-R1/2/3.

RESTART-OPERATOR HYPOTHESIS CARDS (physics-intuition, conjecture):
- H-R1 [lead]: restart at eps_i = (formal monodromy e^{2pi i c_i}, DIAGONAL) (x) (Stokes shears,
  OFF-DIAGONAL), at the rank-2 irregular point u=inf (= the poles eps_i in lambda). Diagonal part is
  ESTABLISHED: c_i = sum_j s_ij^2 (a_i-a_j) (signed BE) = the measured log-T drift of S. Off-diagonal
  = Stokes multipliers = the confluent-Heun coefficient. Test: diagonal of restart = e^{2pi i c_i};
  off-diagonal is TRIANGULAR (Stokes shear) in slope ordering, not dense.
- H-R2: the two Stokes shears live in the {mid,lo} and {mid,hi} 2D carrier spaces (= WS-G's 12x13
  joint, middle = shared partner); their non-commutative composition = S_12 = P_{m→m}. Test: off-diag
  restart weight concentrates in the mid row/col; decoupling one outer (WS-C trivial-coupling) kills
  one shear -> collapses joint -> elementary (matches 'elementary <=> a level decouples').
- H-R3: form factor Gamma_j=(-u'(lambda_j))^{-1/2} is ANALYTIC at the poles (~(lambda-eps_i)); WKB
  branch points are the TURNING points (-u'=0), not eps_i. So the off-diagonal is an irregular-point
  Stokes effect, distinct from the turning-point (window/BE) branching.
Physical payoff: target S_12 IS the Stokes multiplier sigma_{mid,.} of the irregular point, a 2D-reduced
object; the Painleve/CH monodromy data = known formal exponents c_i + these two Stokes shears.

WORKSTREAMS (parallel background agents): WS-PA1 (PA-0/1/2 foundation+pivot, critical path),
WS-PV (Painleve-V tau + lit/irreducibility), WS-CH (direct confluent-Heun connection, parallel track),
WS-NUM (high-precision S_12 + PSLQ, always-on), WS-R (H-R1/2/3 restart structure). See
paper/OPEN_PROBLEM.md addendum.

## Q_restart / cyclic-monodromy line (user) — geometric backbone real, bare relation CONTENT-FREE — 2026-06-01
Code: `experiments/q_restart_probe.py`. User intuition: lambda->u=n/p is degree-3 (3:1), so 3 LZ
sweeps = 1 loop in lambda; the 'restart' at the poles eps_i + the U(3) fact (DPA)^3=I suggest a
Q_restart=DP theory constraining the S-matrix A.
- **Geometry REAL [established]:** deg(lambda->u)=3 (verified); 3 real preimages interlace the poles
  (projectively, wrapping through inf = the 'funny business' at eps_i); WS-F's fixed sheet
  permutation (2,0,1) is the natural 3-cycle P.
- **(i) SALVAGE [established]:** the canonical IP frame is fixed by a GEOMETRIC Coulomb-phase
  subtraction: diabatic-IP S off-diagonal phases diverge as -(c_j-c_i)log T with
  c_i = (1/4) sum_{j!=i} w_ij^2 (a_i-a_j), w_ij=2|gam_i gam_j|/|eps_i-eps_j| (slope-free width) — a
  signed slope-weighted second moment of the widths (sum_i c_i=0). All 9 entries match to <0.02 rad;
  S_canon=e^{i(c_j-c_i)logT}S converges. NEW, clean, useful for the #2 connection problem (pins the
  phase regularization geometrically). [w_ij nomenclature, NOT Gamma.]
- **(ii) NEGATIVE [established]:** in the canonical frame (P S_canon)^3 is NOT diagonal (reloff 0.86).
  And (P D S D^-1)^3=diagonal is achievable for ANY U(3) by reconvention D (Type-1 AND random U(3)
  both reach reloff=0.0000) -> the bare cyclic relation is a CONTENT-FREE general unitary fact,
  exactly the user's own caution. Q_restart=DP does NOT constrain A. A content-ful restart needs a
  NON-TRIVIAL OFF-DIAGONAL pole factor = the non-Abelian Stokes data = open-problem #2.
- **Net:** user's caution (D non-canonical; expect diagonal not I) was decisive and correct; pushing
  it through kills the bare DP relation but extracts the geometric Coulomb regularization
  (c_i = 1/4 sum_j w_ij^2 (a_i-a_j), in w_ij nomenclature). The 'restart at the poles' is genuinely
  off-diagonal = the confluent-Heun coefficient.
  Failed-exploration preserved; the Coulomb-phase result feeds #2 (PA-0/PA-1 frame).

## Deliverable #1 (paper) + #2 (open problem & plan) — 2026-06-01
- **Paper draft:** `paper/type1_lz_working_paper.tex` — 10 evidence-tagged results + dependency map.
  Passed a coordinator physics-reflection pass: corrected R3 universal [EST]->[NS] (scan not proof);
  R8 accessory-exponent count flagged as resting on two agreeing agent frames (WS-A u-frame + WS-E
  Laplace frame), [K,K]!=0 coordinator-verified. DRAFT with 2 named verification TODOs (the R9
  correlation re-check; the R8 accessory re-derivation = PA-0).
- **Open problem formalized:** `paper/OPEN_PROBLEM.md`. THE PROBLEM: compute the off-diagonal Stokes
  coefficient S_12 (=P_{m→m}) of the single rank-2 irregular point — a confluent-Heun connection
  coefficient one accessory parameter above 1F2 — in closed form via the two window actions + the
  turning-point cross-ratio, OR prove irreducibility. KEY LEVER: rank-2 irregular + apparent
  singularity => isomonodromic deformation is Painleve (expected Painleve V), whose connection
  problem has a modern closed theory (Gamayun-Iorgov-Lisovyy tau-functions / CFT c=1 blocks). Central
  conjecture: Gaudin data fixes the Heun accessory parameter ALGEBRAICALLY => named closed form.
- **Plan-of-attack (for review):** PA-0 verify foundation (re-derive scalar Laplace ODE + {0,1,3}
  apparent sing) -> PA-1 identify Painleve type + node-as-accessory? -> PA-2 PIVOT: does Gaudin fix
  the accessory param algebraically? -> PA-3 closed form via Painleve-V/CFT connection formulae (if
  algebraic) | PA-4 validated numerical model + PSLQ recognition (always-on safety net) | PA-5 lit +
  irreducibility. Honest odds: elementary=low; named-closed-form (Painleve-V tau)=plausible if PA-2
  algebraic; irreducibility+numerical model=floor. AWAITING USER REVIEW before dispatch.

## WS-E + WS-C completed — research arc essentially closed — 2026-06-01
Both verified by coordinator; deliverables paper/ws_e_junction_Smatrix.md, paper/ws_c_factorization_locus.md (+ experiments/ws_{e,c}_*.py).

### WS-E (the prize) — VERDICT: honest NEGATIVE-with-structure. The closed VALUE of P_{m→m} is NOT
elementary; it is the off-diagonal Stokes connection coefficient in the genus-0 confluent-Heun /
higher-Weber class, ONE accessory parameter above the 1F2 / Kampe de Feriet of the solvable cousins.
Benchmarked (anchors exact: canonical 0.214724, sampleB 0.021018).
- [established] BE survivals = the two imaginary window actions (elementary); middle survival is the
  genuine non-Abelian leftover (a 2-dof family shifts it independent of BE + double-stochasticity).
- [established] Exact LAPLACE integral rep: psi_j(u)=int_C e^{-iuv}B_j(v)dv, B'(v)=K(v)B,
  K(v)=-i diag(1/a)(H0-vI); saddles = Q4 turning points / WS-G joint.
- [VERIFIED by coordinator, by hand] BLOCKER: [K(v1),K(v2)]=(v1-v2)[D H0,D] != 0 (D=diag(1/a)),
  entries (1/a_i)(H0)_ij(1/a_j-1/a_i) != 0 -> Laplace image non-commuting -> NO product/abelian
  solution. The linear model lacks the 1/t Coulomb term whose Laplace image gives Lin-Sinitsyn/BBGY
  their Fuchsian/rank-1 Bethe structure. No Coulomb => no Fuchsian image => no product.
- [analytic] Laplace scalar ODE: ONE finite accessory singular point (not an H0 eigenvalue),
  anomalous indicial exponents {0,1,3} (gap at 2) + rank-2 irregular at inf -> strictly above
  BBGY's 1F2 (singular only at 0,inf). Independently CONFIRMS WS-A's confluent-Heun class, now in
  the Laplace frame (two independent derivations agree).
- [VERIFIED, NEW cousin-distinguishing fact] Type-1 is NEVER tridiagonal: the extreme-extreme Cauchy
  coupling V_02 is nonzero and DOMINANT (canonical 0.72 > V01 0.60 > V12 0.48). Every solvable MLZ
  cousin (Lin-Sinitsyn/BBGY su(2)-Gaudin) is tridiagonal. Type-1 = the unsolved full confluent case.
- [num] P_{m→m} lies OUTSIDE the entire two-path Stuckelberg interference band on 3/7 overlapping
  samples -> no two-amplitude/product/sech closed form can reproduce it.
- Maximal honest result: exact integral rep + named function class (args = window actions +
  turning-point cross-ratio) + controlled limits (BE; well-separated->incoherent). Closed VALUE
  needs the generically-unknown confluent-Heun connection coefficients.

### WS-C — VERDICT: prediction P-ii CONFIRMED (sharp), one correction.
- factorizable = joint-free = trivial-middle-coupling: all three COINCIDE and are the BOUNDARY of
  the genuine-LZ region (one middle coupling -> 0). NO interior codim-1 elementary locus; on the
  genuine-LZ locus the matrix is ALWAYS non-factorizing. char-poly quadratic factor IRREDUCIBLE over
  Q(gamma,eps,a) (sympy) -> confluent-Heun genuine.
- Factorization locus = middle level decouples from >=1 outer (V_mid,outer->0); algebraically only
  gamma_i=0 or a_i=a_j kill a coupling.
- CORRECTION: classifier is 'middle DECOUPLES', not 'min(Gamma_mid) small' (near_deg_slope:
  Gamma_mid,min=0.010 yet enh=1.31, non-factorizing, because the OTHER middle coupling is O(1)).
  Clean criterion = enh->1 AND absdef->0, delivered only by the trivial-coupling limit.
- Correlations [coordinator re-verifying]: Pearson(log(enh-1),log Gamma_mid,min)=+0.51;
  Pearson(...,log sep/width)=+0.015 (~0). Factorization tracks middle coupling, UNCORRELATED with
  separation -- WS-G's 'separation' framing was a red herring (sampleA factorizes via Gamma_mid~1e-5).
- BBGY eps_2=(eps_1+eps_3)/2 is NOT a factorizing locus for the LINEAR model (absdef ~0.07-0.14);
  the hyperbolic 1F2 reduction does NOT transfer.

### NET (arc essentially complete): an honest, parsimonious THEORY. Class (WS-A/WS-E: genus-0
confluent-Heun, one accessory above 1F2), obstruction (WS-D/co-scaling: slope-free width lemma ->
permanent overlap; never tridiagonal), mechanism (WS-G: 12x13 joint = the non-Abelian remainder),
classifier (WS-C: elementary <=> a level decouples), gold data (WS-F). Type-1 N=3 = the MINIMAL
full (non-tridiagonal) confluent MLZ whose middle survival is an irreducible confluent-Heun Stokes
coefficient. The elementary pieces (BE, integral rep, all limits, the bimodal locus) are closed; the
generic middle-survival VALUE is the named-but-unevaluated connection coefficient.

## WS-G completed + coordinator integration (VERIFIED, promoted) — 2026-06-01
Deliverables `paper/ws_g_stokes_graph.md`, `experiments/stokes_graph.py`, figs/{stokes_sepA,
stokes_overlapB}.png. First agent result promoted with FULL confidence (load-bearing claims
independently verified by coordinator).
- **Joint diagnosis SUPPORTED (numerically-supported, figures decisive).** Well-separated sampleA:
  short DISJOINT Stokes fans, 0 joints, enhancement 1.000× (incoherent product EXACT — I verified:
  P_mid=0.98597 vs incoherent 0.98587). Overlapping sampleB: the 12 and 13 Stokes lines CROSS at 4
  joints (u≈−0.13±0.95i, −0.07±0.85i), node (exact crossing) on the real axis BETWEEN the two joint
  clusters, enhancement 120.5× (consistent with WS-F oracle 128.7×). Joints are 12×13 type = the two
  off-diagonal pairs sharing the MIDDLE level = exactly P_{m→m}.
- **Localizes WS-A's non-Abelian remainder:** the transcendental off-diagonal Stokes coefficient =
  the JOINT HOLONOMY (the GMN junction S-matrix at the 12×13 joint, node between the two joints, as
  a function of the two window actions). Joint-free ⇔ elementary (BE×Weber product); jointed ⇔
  genuine confluent connection coefficient. Hands WS-E a concrete target.
- **REFINEMENT (coordinator, unifies WS-D+WS-G):** sampleA's enhancement=1 has Γ_mid,lo≈0,
  Γ_mid,hi=0.002 — the middle level BARELY COUPLES. So joint-free ⟺ TRIVIAL middle coupling
  (Γ_mid→0 or →∞), NOT 'genuine-but-separated' (which co-scaling proved doesn't exist). The 15%
  'incoherent works' = trivial-middle-coupling tail. Whenever the middle level genuinely couples
  (Γ_mid O(1)), there is a joint and non-factorization. Single coherent picture: co-scaling
  (permanent overlap) + WS-G (joint) + the 15/85 split.
- **Honest residuals (flagged by WS-G, accepted):** detector is one-directional (3 false negatives
  in the far-apart tail; 'no joint' is only a lower-bound witness; full GMN trajectory/junction
  rules needed for a bidirectional classifier); 'joint strength' is a geometric proxy, not the
  literal relative WKB action; WS-G shows joint PRESENCE tracks enhancement but does not COMPUTE
  P_{m→m} (that is WS-E).
- This is the verified mechanism behind the obstruction; with WS-A (class) + WS-D/co-scaling
  (obstruction) + WS-F (gold data) + WS-G (mechanism), all launched workstreams are in. Remaining:
  WS-C (factorization/trivial-coupling locus = joint-free locus, prediction P-ii) and WS-E (the
  joint-holonomy closed form).

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
  off-diagonal Weber/Stokes exponent = the BE rate s_ij^2|a_i-a_j|. These are residues of the ACTION/Stokes
  one-form, NOT the geometric connection (adiabatic transport ~1/u² ⇒ no simple-pole residue).
  Strengthens E1, no contradiction.
- **CORRECTS my framing (candor):** the exact-crossing NODE is an ORDINARY point of the honest 3×3
  system — it carries NO local connection data and NO accessory parameter. My repeated "exact
  crossing = local reduction handle / rigorous H-B" was WRONG: it conflated the adiabatic-frame
  gauge singularity (eigenvectors rotate fast at the degeneracy) with a local ODE feature. In the
  diabatic frame u* is regular. ⇒ **WS-B as specced (model the node as a local 2×2 connection
  coefficient pinning P_{m→m}) is CANCELLED.** The node's 2.5–104× enhancement is a GLOBAL Stokes
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
  exponents) — that prefactor IS our missing P_{m→m}/off-diagonal. Exact handles: ε₀=0 slice
  (confluent hypergeometric, P_{m→m}(0)=2e^{−πg²/b}/(1+e^{−πg²/b}), their Eq.49) + adiabatic-limit
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
1. **H-B** exact (avoided) crossing → middle survival P_{m→m} via 2×2 local rotation + unitarity.
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

---

## 2026-06-02 — WS-O2 + WS-O3 outcomes (parallel constructive workstreams)

Context: after the strategic re-rank (concede exact-WKB product O1 yields a *constructive*
but not *concise* object; the hard content relocates to one junction Stokes constant), fired
two parallel workstreams aimed at the two distinct virtues — a *named* closed form (O2) and a
*usable concise* formula (O3).

**WS-O2 (exact integral rep + c=1 connection-constant probe)** → commit d735490.
- `S_xj = e^{+i a_j u^2/2} oint_{C_j} e^{-iuv} B(v) dv`, `B' = -i diag(1/a)(H0 - vI) B`,
  machine-verified (operator-identity residual ~2e-14). [analytically-derived]
- Gauge lemma `a->a+c` invariant ⇒ WLOG all a_j>0, single Stokes wedge. [established]
- Literal double-precision matrizant quadrature hits a Stokes-dominance wall
  (cond ~ e^{R^2/2a} → 1e19); recessive-only / high-precision engine OWED. value validated via
  oracle to <=1e-6. [numerically-supported]
- **c=1 Barnes-G/Euler-Gamma probe: NEGATIVE, decisive [gold-gated].** Neither the rank-2
  single-sigma band (5/8 strata strictly outside) nor a finite rank-3 3-amplitude product
  (forces |cosPhi|>1 on the overlapping strata) reproduces P_mid. Closes ONLY as the full c=1
  **Fredholm-determinant / block-Toeplitz** connection constant (LNR arXiv:1806.08344); PV/
  Barnes-G applies only on the one-link-decoupling locus. Sharpens R10 from the connection-
  constant side. ⇒ the "named closed form" deliverable is a Fredholm determinant, NOT a finite
  special-function product.

**WS-O3 (uniform Dykhne-Stuckelberg)** → commit 16a2db9.
- `P_mm ~= A0 + A_ret + 2 sqrt(A0 A_ret) cosPhi`, `A0=p_lo p_hi`,
  `A_ret=(1-p_lo)(1-p_hi)p_lh`, `cosPhi=clip[-(d_lo+d_hi)(1-3(1-chi)),-1,1]`, `p_X=e^{-2pi d_X}`.
  Return path is via the OUTER lo-hi link (forced by verified u-ordering). [analytically-derived
  structure; numerically-supported law, one fitted O(1) const K~3.0, robust]
- Oracle-validated 53-pt sweep: RMS ~1-2%, max ~4% for delta<=0.25 (the generic regime); the
  chi-dressing halves the incoherent error and captures the shape trend. Cleanly exhibits R11
  {two window actions (scale) + chi (shape)}.
- FAILS (>10%) in deep adiabatic overlap delta>0.5 (sampleB off by 0.27) — exactly the rank-3
  regime (R4 permanent extreme-extreme coupling / R10). No simple uniform formula reaches it.

**Net (honest).** Two distinct deliverables now exist: (i) an exact concise *formula* (O2
integral rep) whose evaluation is conditioning-limited and whose *named* form is a c=1 Fredholm
determinant (not elementary/finite-product — the finite-product hope is falsified); (ii) a
concise *computable approximation* (O3) good to ~1-2% in the generic separated-to-moderate
regime, failing in the genuine rank-3 deep-overlap core. The exact-elementary/finite-special-
function target is now strongly evidenced NOT to exist; the realistic paper = O2 exact-rep +
Fredholm-det classification as the rigorous spine, O3 as the usable results-section formula,
oracle as ground truth. OWED: recessive-stable O2 engine; explicit LNR Fredholm kernel for our
rank-3 point; asymmetry-aware O3 (depends on d_lo+d_hi, not the ratio); proofs to promote NS.

---

## 2026-06-02 — WS-O2b: explicit Fredholm/Widom kernel (open problem ii) — OBSTRUCTED, honestly

Constructive swing at turning the c=1 Fredholm-determinant classification (WS-O2/R11) into a usable
a-priori construction for our rank-3 (three-distinct-rate) confluent point. → files
experiments/ws_o2b_fredholm.py, paper/ws_o2b_fredholm_kernel.md.

- Literature verified (M2): CGL/Widom recipe `tau=det(1-K)` of an IIKS/Plemelj operator; one-circle
  = block-Toeplitz Widom determinant (1712.08546); regular case 1608.00958; irregular generalized-
  Bessel example 1705.01869; LNR PV connection 1806.08344 (title confirmed correct). The recipe
  takes the LOCAL Stokes data as INPUT; PV/PIII are constructive because the irregular parametrix is
  a known special function (Weber/Bessel) with known Stokes matrices.
- POSITIVE by-product [AD]: Wasow recursion at the rank-2 irregular point gives
  `Y_inf=G(v)exp(½A0 v^2 + D1 v + Theta ln v)` with `Theta=diag(c_i)` emerging EXACTLY — an
  independent Laplace-frame confirmation of R8 (formal monodromy = Coulomb exponents).
- VERDICT [AD, verified]: NO usable construction. Two concrete break points: (#1) no closed-form
  local parametrix exists for a rank-2 irregular point with THREE distinct rates — its Stokes
  matrices ARE the unknown sigma = the answer (dim-6 wild char variety, R10); the Widom symbol
  g=Y_inf^{-1}Y_* is never v-independent and is Stokes-walled (cond ~ e^{R^2/2a} -> 1e18). (#2) the
  Laplace v-frame connection matrix is NOT the physical S (unitary factor off by up to 0.39; contour
  selection is the physics). Given resummed Stokes data, P CAN be re-expressed as a block-Toeplitz
  determinant, but answer-in/answer-out at equal cost, exposing no Nekrasov/theta-truncation or
  factorization.
- Gold-gate context: the usable computable object remains the physical-u-frame numerical RHP/oracle
  (sampleB P=0.02102, err 7.2e-7 at T=60; <=1e-9 at T>=120), which already reaches deep overlap where
  WS-O3 fails.

NET: open problem (ii) resolved NEGATIVELY but informatively — the Fredholm determinant is an
IRREDUCIBLY CLASSIFICATORY description of P_mid (sharpens R9-R11 from the constructive side), not a
closed construction. This is the third independent confirmation that no concise/elementary closed
form exists for the rank-3 object. Revival route (precision-limited, not structural): mpmath
recessive-split symbol to beat the Stokes wall and READ sigma; closed v-frame<->S contour map.
Paper updated: new R11b in sec:fredholm; open-problem (ii) marked attempted/obstructed.

---

## 2026-06-02 — physics-reflection pass: every paper claim cross-checked vs its script

Adversarial verification of paper/type1_lz_working_paper.tex (R1-R14, R11b, F1). Re-ran the
load-bearing scripts; ALL headline numbers reproduced. Result: VERDICT = KEEP (no result changed);
only citation/number hygiene fixed.

VERIFIED (re-run 2026-06-02):
- R1: commutator 8.9e-16, shared eigenbasis 3.3e-16 (ring_structure.py).
- R4/R5: width slope-independence = 0 exactly; max sep/width = 1.60 over 35079 samples (coscaling.py).
- R3: real crossings on all anchors (canonical 1.8e-9, sampleB 1.6e-9, WSF_near_node 1.9e-8);
  independent Brent-refined re-scan: 96.9%/351 gap<1e-3, median gap 3.06e-9 (regenerated; the prior
  99.8%/1290 was docstring-only, not executed by current script).
- R8: accessory=node rational (canonical u_*=-187/750, v_*=-748/375; sampleB rational);
  Theta=diag(c_i) EXACTLY (max|Theta-c|=0.0, canonical) -> R11b by-product verified.
- R9: VERDICT RANK-3; sampleB band [0.8327,0.8800] OUT=True. Count is 7/14 (cap), NOT 5/8.
- R10: operator-identity residual canonical 5.77e-14, sampleB 1.87e-13 (paper said ~2e-14 -> O(1e-13)).
- R11: PROBE1 5/8 outside, PROBE2 4/8 force |cosPhi|>1, Fredholm-det verdict (ws_o2_integral.py).
- R7: enhancement canonical 2.58x, sampleB 120.54x (anchor_experiment.py) -> "2.5-130x" confirmed.
- R11b: Widom symbol relresid O(1), Stokes wall cond->1.8e18; v-frame misses oracle; u-frame RHP
  gold sampleB 7.2e-7 (ws_o2b_fredholm.py).
- R12/R13: separated RMS 1.23%, weak-mod 1.90%, deep-overlap 11.06%, sampleB off 0.2670 (ws_o3).
- R14 anchors: canonical dev 1.6e-9, sampleB dev 2.1e-7 (num_S12.py).

PROVENANCE FIXES APPLIED (all citation/number hygiene; no result altered):
1. R10 residual ~2e-14 -> O(1e-13) with both anchor values.
2. R9 "5/8" -> "7/14" (cap); 5/8 re-attributed to the c=1 probe (ws_o2).
3. R7 removed dangling figs/stokes_*.png (figs/ dir does not exist); cited exact enhancement factors.
4. R10/R14 dangling oracle_report.md -> num_S12_model.md + re-verify date.
5. R3 softened/grounded: prior 99.8%/1290 (docstring) + independent re-scan 96.9%/351, median 3.1e-9.

NOT independently re-derived this pass (lower-risk, memo/structural; flagged for completeness):
R2 (genus-0, gate_test_genus.md — symbolic memo), R6 (coupling dominance — stated values).

NET: paper is internally consistent and provenance-clean; tags honest. Draft remains DRAFT (no
pdflatex in env to compile; R2/R6 memo-only). The three-deliverable story (integral rep / c=1
Fredholm-det classification / uniform formula + oracle) stands and survives adversarial review.

---

## 2026-06-02 — R3 node universality is a THEOREM (Owusu-Wagh-Yuzbashyan 2009)

User pointed to the proof: H. K. Owusu, K. Wagh, E. A. Yuzbashyan, "The link between
integrability, level crossings and exact solution in quantum models," J. Phys. A 42, 035206
(2009), arXiv:0807.0259. (Lead author = the user.) THEOREM: any H0+uA possessing a nontrivial
commuting partner has exact level crossings as u varies, violating von Neumann-Wigner. Type-1 N=3
satisfies the hypothesis (commuting partner = R1), so the universal real exact crossing (our "node")
is GUARANTEED.

Impact:
- R3 promoted numerically-supported -> ESTABLISHED (theorem + our corroborating Brent re-scan
  96.9%/351, median gap 3.1e-9; rational location u_*=-187/750 from R8).
- "Analytic proof of node universality owed" REMOVED from the owed list.
- Lineage clarified: OWY (integrability => crossings) is the upstream theorem for our entire
  node-based construction; our R1 commuting ring is exactly its hypothesis. Our contribution is
  downstream: node = accessory point of the connection problem (R8), the rank-3 classification (R9),
  and the three deliverables.
- Paper updated: R3 [EST] with citation; owed-proofs item (iii) now lists only accessory=node (R8),
  the ~1.6 bound (R5), and the R11 finite-product falsification as a region statement.

Owed proofs remaining (post-OWY): accessory=node in general (R8); ~1.6 overlap bound (R5);
finite-product falsification as a parameter-region statement (R11). Strategic: direction (1)
[prove node] is now CLOSED by the literature; the live productive threads are (2) physical meaning
of the node, the R8 accessory=node proof (now well-posed since node is a theorem), the chi-channel
map (R13), and publication.

---

## 2026-06-02 — R8 "accessory = node" PROVEN (WS-R8) + independent coordinator verification

WS-R8 (physics-derivation) proved v_* = E_* as a polynomial identity over Q(gamma,eps,a). Files:
experiments/r8_proof.py, paper/r8_accessory_node_proof.md.

THEOREM (R8): the finite apparent singularity v_* of the Laplace connection ODE equals the OWY node
energy E_*. Two proofs:
- Route A (resultant): substitute E=v_* into Phi(E)=Res_u(chi, d_E chi); numerator reduces to the ZERO
  polynomial on the gauge slice a=(s,1,2) with eps,gamma fully symbolic; lifted to all generic params by
  slope-affine (a->alpha a+beta) covariance (v_* and E_* both transform as energies).
- Route B (structural, mechanism): M(v)w=mu w  <=>  H(-mu)w = v w, so an M(v)-eigenpair is an eigenpair
  of H(u=-mu) at energy v. The cyclic Wronskian W(v)=det[e0,M e0,M^2 e0] is LINEAR in v; its single root
  v_* is where e0's M(v_*)-Krylov space drops rank, which forces chi(u,v_*) to have a double root u_*
  (Disc_u chi(.,v_*)=0), i.e. the node (chi=d_E chi=d_u chi=0). Bonus: u_* = rational repeated root.

INDEPENDENT COORDINATOR CHECK (M1, decisive leg): re-derived from scratch in exact rational arithmetic
on 2 fresh random samples (incl. messy non-slice a=(-1,1/3,2) and a=(-5/4,2/3,11/7)): deg_v W = 1;
Disc_u chi(u,v_*) == 0 (True) AND Phi(v_*) == 0 (True) on both. Confirms the identity AND the
gauge-covariance lift (samples were not in slice form). Agent also reports 31/31 exact random samples.

IMPACT:
- R8/PA-2 promoted: "verified on 4 samples" -> ESTABLISHED (proven identity + independent verification).
- COROLLARY now fully proven: ALL monodromy data of the rank-3 connection problem are algebraic in
  (gamma,eps,a): formal exponents Theta=diag(c_i) (R11b) AND accessory v_*=E_* (R8). The "named Fredholm
  determinant" (R11) has fully explicit, algebraically-determined data.
- Owed proofs remaining: ~1.6 overlap bound (R5); finite-product falsification (R11) as a region statement;
  the strict no-log certificate as a single all-symbol identity (secondary; {0,1,3} indices proven on slice).
- Paper updated: R8 [AD]->[EST, proven] with mechanism; owed-proofs item (iii) trimmed; dep map updated.

LINEAGE COMPLETE: integrability (commuting partner, R1) -> exact crossing (OWY 2009, R3) -> accessory
parameter pinned to it (R8) -> all monodromy data algebraic. The spectral-side theorem (OWY) and the
dynamical-side theorem (R8) are now joined. Integrability does not make P_mid elementary (Abelian
ceiling) but it DOES rigidify the dynamics: the transition amplitude is transcendental with
algebraically-determined data.

---

## 2026-06-02 — R8 memorialized as a review Letter (paper/accessory_node_letter.tex)

Wrote a journal-style (PRL revtex) review Letter, "Pinning the accessory parameter: a
spectral-dynamical identity for the integrable multistate Landau-Zener problem." Structure:
- Intro: motivates the proof via the accessory parameter as THE marker of the solvable->transcendental
  jump (hypergeometric: none, closed-form; Heun: one, no closed form). Frames its value to the LZ
  problem: it does not make P_mid elementary (no-go stands) but removes the one free transcendental
  modulus, so the connection problem has fully explicit algebraic data.
- Setup: the two objects on two sides (node E_* spectral; accessory v_* dynamical, W(v) linear).
- Theorem v_*=E_* with canonical instance; Route A (resultant identity) + Route B (structural pencil
  duality M(v)w=mu w <=> H(-mu)w=v w, Krylov rank-drop <=> double eigenvalue) + gauge covariance.
- Dedicated section expanding the NONTRIVIALITY of the spectral<->dynamical link: the two objects are
  of different type (time-domain generator spectrum vs Borel-plane singularity of the solution
  operator); in general accessory-parameter theory the apparent-singularity location is a FREE modulus
  decoupled from spectral data (the crux of Heun's intractability); Type-1's identity is a coincidence
  ENGINEERED by integrability via the diag(1/a) duality; strip the Cauchy structure and it breaks.
  Conceptual payoff: OWY (spectral) + this (dynamical) => "integrable but not solvable" = transcendental
  amplitude with ALGEBRAICALLY-determined data (rigidity, not solvability).
- Avenues opened/unblocked: (i) connection constant now fully specified (Theta=diag(c_i) + v_*=E_* all
  algebraic) -> well-posed Fredholm/Nekrasov evaluation; (ii) constructive exact-WKB with explicit
  turning-point/joint data; (iii) algebraic equation for the elementary boundary (W_1=0); (iv) a
  spectral duality to generalize to all N (conjecture: accessory variety = spectral-degeneracy variety
  for Type-1 hierarchy); (v) theorem-backed physical handle on the node (explicit u_*); (vi) owed: the
  fully-symbolic no-log certificate.

Status note: the clean isolated re-run of r8_proof.py (default CASE 1) is still grinding on this box
(symbolic resultant heavier than the memo's "closes quickly" suggests under current limits); the proof
stands on the INDEPENDENT coordinator verification (exact rational arithmetic, fresh non-slice samples,
2026-06-02) regardless. Reproducibility caveat to revisit: add a fast default path to r8_proof.py.

---

## 2026-06-02 — repo-wide nomenclature sweep: P₂→₂ retired -> P_{m→m}; ε-ordering enforced

Per user directive, applied the slope-middle survival notation standard across the whole repo and
strictly enforced the ε-ordering convention.

- NOMENCLATURE.md (authoritative) updated: P_{m→m} is canonical (m=argsort(a)[1], slope-middle);
  P₂→₂ and S₁₂ marked DEPRECATED (the literal '2'/'1,2' invites an energy/ε-index misreading).
  Added a strict ENFORCEMENT clause: (1) label channels by ε-index (ε_0<ε_1<ε_2; P[i,j],𝒮,H₀,Γ_j all
  ε-indexed); (2) identify extreme/middle BE roles by slope (argsort(a)), never energy rank;
  (3) write the open survival as P_{m→m}. Per-context forms: tex P_{m\to m}; prose ASCII P_{m->m};
  Python P_mm (brace-free -- P_{m->m} inside an f-string parses as a {m->m} field and breaks).
  Canonical code identifiers P_mm/P_mm_fast/P_mm_oracle/P_mm_model (= P[mid,mid]); these ARE the
  standard names (the legacy ASCII slope-middle forms were renamed to the canonical P_mm* form).
- Sweep: 244 P-arrow-2 occurrences replaced across 29 authored files (logs, SESSION_SYNTHESIS,
  META_REVIEW, all paper/ memos+tex, all experiments/ scripts+memos). Curated exact arrow-bearing
  tokens only -> never matched the slope-middle identifiers or Richardson vars P_2T/P16. uploads/ left untouched
  (source material, per NOMENCLATURE policy).
- Bug caught & fixed: P_{m->m} inside Python f-strings broke compilation (stokes_graph.py); converted
  all .py occurrences to brace-free P_mm. ALL project .py now compile; canonical anchor intact
  (P_mm = 0.214724).
- ε-ordering audit: no violations. restart_probe.py P[1,1] is a FALSE POSITIVE (the matrix is
  reorder_slope'd to sector order [lo,mid,hi] first, so index 1 = slope-middle). All middle-survival
  code anchors mid via argsort(a). The two .tex papers' \Pmid macro = P_{m\to m} (done earlier).

---

## 2026-06-02 — CORRECTION: "c=1" classification was the rank-2 label misapplied -> GL₃/W₃ (c=2)

Prompted by a colleague's note (node-pinned oper strategy) and its sharpening: our R9/R10/R11
classification called the rank-3 object "c=1 (free-fermion) confluent-Garnier." That is WRONG as a
label. c=1 is the RANK-2 free-fermion point (sl₂: Painlevé VI/V/III, the Gamayun-Iorgov-Lisovyy
connection constants). The genuine rank-3 / GL₃ / sl₃ one-irregular-point object lives in the
sl₃-Schlesinger / W₃ arena, whose self-dual central charge is c=N−1=2, NOT 1.

Corrected (authoritative artifacts): type1_lz_working_paper.tex (title, abstract, sec:class heading,
R9/R11 verdict, unifying picture, honest-status), accessory_node_letter.tex (intro classification),
SESSION_SYNTHESIS.md R10, cap_connection_formula.md headline. The rank-2 `c=1` references that
correctly name the PUBLISHED PV/LNR object we test AGAINST (the falsification target) are kept as-is.

Substance unchanged: the verdict (rank-3, dim-6 wild char variety, no finite-product closed form,
gold-gated) stands; only the CFT/rep-theory label is corrected. Triangulation from BOTH sides now:
(us) WS-O2/O2b gold-gated + Fredholm-kernel obstruction; (literature, via colleague) Gavrylenko's
sl₃/W₃ work is four-point Fuchsian and the constructive three-point / confluent GL₃ connection
formula is unpublished/open. So no worked rank-3 formula exists -- confirmed numerically AND from
representation theory. CGL/Widom is the right *machinery* but needs the local parametrix we lack.

NOTE: superseded planning memos (RESEARCH_OVERVIEW.md, OPEN_PROBLEM.md, pa1_isomonodromy_foundation.md)
and earlier chronological log entries retain historical "c=1" language; not rewritten (history
preserved). The authoritative paper/Letter/synthesis/cap now carry the corrected GL₃/W₃ (c=2) label.

---

## 2026-06-02 — WS-RH: node-pinned GL₃ RH spec + high-precision evaluator (rung 1 done; rung 2 = structural negative)

Files: paper/gl3_rh_problem.md, experiments/gl3_rh_solver.py. Coordinator independently verified and
CORRECTED the agent's "rung 2 achieved" overstatement.

RUNG 1 — ACHIEVED. Explicit algebraic RH spec for the Laplace-dual oper: irregular point v=∞ (Poincaré
rank 2, three distinct rates, formal exponents Θ=diag(c_i) re-confirmed to machine precision, an
independent Laplace-frame re-derivation of R8); apparent point v_*=E_* rational, {0,1,3}, no-log; Stokes
graph (4 rays). NEW analytic result: the G± physical↔Laplace normalization maps are DIAGONAL
(channel-wise), derived by steepest descent through v_j*=u a_j -> saddle phase = diabatic Stark phase +
c_j log-drift. So the v-frame->S map's diagonal part is algebraic/pinned (up to a sample-independent
sqrt(-i)/sector convention); the ONLY non-diagonal piece is the Stokes σ.

RUNG 2 — NOT achieved as intended (independent RH evaluator beating the oracle). Two engines:
- Engine A (genuine v-plane RH/oper solve, stripped frame Y=F(v)Z, mpmath): BEATS the WS-O2b
  conditioning wall — Z stays finite & resolved (canonical max|Z|=327; sampleB up to 6e15 at dps70)
  where double precision overflowed (cond~e^{R²/2a}->1e18). BUT STALLS: no diagonal/Sinkhorn/polar G±
  dressing of |Z|² recovers S. COORDINATOR-VERIFIED (canonical): raw|Z|²=1.21, Sinkhorn=0.355,
  polar=0.152 vs oracle 0.2147 — none match. Stable across dps/order => STRUCTURAL, not numerical.
- Engine B (passes gold gate: canonical 8.4e-8, sampleB 4.6e-7): coordinator read the code — it solves
  Y'=-iH(u)Y in the TIME DOMAIN with Richardson = the SAME METHOD as oracle.py, independently coded. NOT
  an independent RH route; its agreement and deep-overlap reach are inherited. Agent docstring + md §7
  called it "independent" — CORRECTED in both files.

NET / VALUE: a precise STRUCTURAL upgrade of WS-O2b. The wall is no longer "maybe just conditioning":
high precision removes the conditioning excuse (Z computed exactly), the stall persists, so the sole
remaining transcendental obstruction is the off-diagonal Stokes constant σ (the rank-3 connection
constant), with EVERYTHING ELSE algebraic (data) and diagonal (G±). σ is now isolated as the unique
non-abelian unknown. Dovetails with the "BE+geometry captures all but the non-abelian σ" finding and the
Abelian-ceiling thesis. Rung 3 (closed W₃/Nekrasov σ formula) not attempted (frontier).

---

## 2026-06-02 — Paper subsection: elementary skeleton + two transcendentals + Abelian ceiling (R15, R16)

New §sec:skeleton in type1_lz_working_paper.tex + experiments/skeleton_two_transcendentals.py.

R15 [EST]: P (doubly stochastic, 4 DOF); BE fixes the 2 extreme-slope survivals; BE+DS => every entry
is an explicit affine function of exactly TWO unknowns, the middle survival Pmid and one off-diagonal
b:=P[hi,lo] (lo->hi chirality). Reconstruct-from-2 exact to ~2e-10. P strongly NON-symmetric
(||P-P^T||~0.6-1.0) => b independent of Pmid (no symmetry collapse to one unknown; symmetry would force
b=(1-Pll-Phh+Pmid)/2, off by ~0.3).

R16 [NS]: vs the incoherent independent-crossing model (BE+geometry, no interference): Pmid is the
universal hard number (off by up to 0.985 in deep overlap); b is a SOFT transcendental, collapsing to its
incoherent (cyclic) value in BOTH the separated (|db|~4e-3) and deep-overlap (|db|~9e-4) limits, only
mildly dressed (|db|~0.15) in the moderate middle.

Why BE+geometry captures so much (three-layer hierarchy): (i) extreme survivals exact (extremes never
recombine); (ii) directed-cycle skeleton = formal monodromy (abelian, geometric); (iii) the large (2.5-
130x) interference on Pmid is a coherent SEMICLASSICAL effect (Stuckelberg/period phase, chi) not the
transcendental constant. The genuine transcendental residue is the rank-3 Stokes joint sigma, subleading
except at turning-point coalescence (deep overlap).

Abelian-ceiling reading (observable face of R1): integrability fixes all abelian data + the diagonal
algebraic G+- map; the SOLE non-abelian datum is sigma. WS-RH confirms from the RH side (wall beaten, data
algebraic, G+- diagonal, stall isolates sigma). Punchline: transcendence localizes on Pmid because the
middle is the only RECOMBINING level (two interfering paths); the extremes are one-way cascades => BE-exact.

This ties together the user's "treat Pmid as mysterious" question, the WS-O3 semiclassical success, and the
WS-RH sigma-isolation into one structural statement. Coordinator-verified numerics (skeleton script).
