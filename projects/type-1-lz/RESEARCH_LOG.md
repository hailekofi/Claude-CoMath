# Research Log: Type 1 Landau–Zener

_Last updated: 2026-06-03 (T4 selector resolved: orientation = energy-position of the real node)_

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
strategy memo `paper/drafts/type1_lz_strategy_memo.tex`.

## Theory-building: autonomous research program drafted — 2026-06-01
`paper/program/RESEARCH_PROGRAM.md` (the plan) + `paper/drafts/type1_lz_working_paper.tex` (living paper, DRAFT,
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
Deliverables paper/notes/cap_connection_formula.md, experiments/cap_connection_formula.py. The capstone
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
Deliverable paper/notes/pa1_isomonodromy_foundation.md; coordinator verification experiments/pa2_accessory_algebraic.py.
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
Deliverables experiments/num_S12.py, paper/notes/num_S12_model.md (num_S12_dataset.pkl gitignored).
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
Deliverables paper/program/restart_structure.md, experiments/restart_probe.py. All measured in the
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
Deliverable paper/notes/ch_direct_connection.md; benchmark scripts archived experiments/ws_ch/.
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
Deliverable paper/notes/pv_tau_route.md. Coordinator-verified the load-bearing logic.
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
paper/program/OPEN_PROBLEM.md addendum.

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
- **Paper draft:** `paper/drafts/type1_lz_working_paper.tex` — 10 evidence-tagged results + dependency map.
  Passed a coordinator physics-reflection pass: corrected R3 universal [EST]->[NS] (scan not proof);
  R8 accessory-exponent count flagged as resting on two agreeing agent frames (WS-A u-frame + WS-E
  Laplace frame), [K,K]!=0 coordinator-verified. DRAFT with 2 named verification TODOs (the R9
  correlation re-check; the R8 accessory re-derivation = PA-0).
- **Open problem formalized:** `paper/program/OPEN_PROBLEM.md`. THE PROBLEM: compute the off-diagonal Stokes
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
Both verified by coordinator; deliverables paper/working_sessions/ws_e_junction_Smatrix.md, paper/working_sessions/ws_c_factorization_locus.md (+ experiments/ws_{e,c}_*.py).

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
Deliverables `paper/working_sessions/ws_g_stokes_graph.md`, `experiments/stokes_graph.py`, figs/{stokes_sepA,
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
`paper/notes/coscaling_derivation.md` + `experiments/coscaling.py`. Done by coordinator (not delegated;
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
Deliverable `paper/working_sessions/ws_a_riemann_scheme.md`. Claims structurally forced (H(u) polynomial), low risk.
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
Deliverable `paper/working_sessions/ws_d_nonabelian_E.md`; verification `experiments/ws_d_verification.py`.
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
Full result: `paper/notes/gate_test_genus.md`. Ran the review's decisive test; it OVERTURNS the review's
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
Full review: `paper/notes/REVIEW_kz_line.md`. Peer-review BEFORE theory-building (user-requested).
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
(`paper/drafts/type1_lz_strategy_memo.tex`). Key findings:
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
Attack lines generated + reflected + tournamented in `paper/program/attack_lines.md`. Grounding facts
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
  written: `paper/drafts/type1_lz_strategy_memo.tex`. **Rank: S1 KZ/Gaudin monodromy > S2 constraint
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
experiments/ws_o2b_fredholm.py, paper/working_sessions/ws_o2b_fredholm_kernel.md.

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

Adversarial verification of paper/drafts/type1_lz_working_paper.tex (R1-R14, R11b, F1). Re-ran the
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
experiments/r8_proof.py, paper/notes/r8_accessory_node_proof.md.

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

## 2026-06-02 — R8 memorialized as a review Letter (paper/drafts/accessory_node_letter.tex)

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

Files: paper/notes/gl3_rh_problem.md, experiments/gl3_rh_solver.py. Coordinator independently verified and
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

---

## 2026-06-02 — Commuting family as isomonodromic deformation (user idea: map family -> 2 transcendentals)

Probe: experiments/deformation_family_probe.py. Tested whether the 2-parameter commuting family maps
usefully onto the two transcendentals {P_mm, b}.

F1 [established, machine prec]: different slope vectors a (same gamma,eps) COMMUTE (||[H^a,H^a']||~1e-16)
and SHARE one a-INDEPENDENT eigenbasis (shared-eigvec defect ~1e-15) at all u tested. => the
derivative-coupling SEED W_ij=<phi_i|phi_j'> is a-independent: the ENTIRE transcendental content of the
WHOLE 2-param family is ONE geometric object W. (Strong confirmation/sharpening of R1.)

F2 [NS]: along a path a(t) in the family (separated->deep), {P_mm(a), b(a)} vary SMOOTHLY (P_mm 0.565->
0.193, b 0.293->0.867) along with the elementary BE actions. Each member's {P_mm,b} = holonomy of the
FIXED W against that member's a-linear phases.

VERDICT (honest): the map is the ISOMONODROMY structure -- W = monodromy/Stokes seed (a-invariant);
a = deformation times; {P_mm(a),b(a)} = tau-data flowing along the deformation. FRUITFUL as STRUCTURE
(one rigid shared seed for the whole family; characterizing W's Stokes sigma solves the entire family)
and as the correct SETUP for the Schlesinger/Garnier tau-route (= the named-closed-form deliverable;
matches WS-RH's sigma-isolation and the GL3/W3 frontier, approached from the deformation side with a clean
separated-regime initial condition). NOT a shortcut around the Abelian ceiling:
 - the a-INVARIANT is the SEED W, NOT the two transcendentals (they vary) -- the family preserves/transports
   the content, never produces it (R1);
 - no free deep-overlap shortcut: the holonomy of W with deep-overlap phases is still the hard connection
   problem; the genuine lever (derive the a-flow / Schlesinger ODE with the separated elementary IC) is
   itself the frontier tau-derivation.
 - dimension match 2=2 is N=3-special ((N-1)^2-2 transcendentals vs N-1 family params); the real link is
   the RH deformation<->monodromy correspondence, not the numerology.
Actionable hard-but-right next step (if pursued): derive/integrate the isomonodromic a-flow (Schlesinger
for our node-pinned GL3 point) from the separated-regime IC to deep overlap.

---

## 2026-06-02 — Core principle elevated + WS-AFLOW derivation roadmap scoped

User directive: make the single-seed-W holonomy insight the CORE of the project.
- Elevated in type1_lz_working_paper.tex: new "Organizing principle" clause in the abstract + a
  "The central object" paragraph in the intro threading R1/skeleton/WS-RH/deformation into one statement
  (characterize the Stokes data sigma of the single a-independent geometric seed W). SESSION_SYNTHESIS
  §1 now leads with the central principle.
- Scoped paper/working_sessions/ws_aflow_derivation_plan.md (physics-derivation roadmap, NOT a proof): derive the
  deformation flow of {P_mm(a),b(a)} along the commuting family using W a-independent + the integrable
  zero-curvature (Lax) M_a. MILESTONE 1 = make-or-break: does the a-flow CLOSE into a finite-dim ODE
  (and is sigma conserved)? Decisive test: finite-difference d_a{P_mm,b} (deformation_family_probe) vs
  the proposed closed RHS. Three outcomes: (i) closes + sigma conserved (isomonodromic) -> transport
  lever; (ii) closes into a 2-component ODE -> usable; (iii) doesn't close -> organizing-only (likely
  failure mode, consistent with Abelian ceiling). Then M2 derive flow, M3 separated IC (near-elementary,
  well-conditioned, NOT elementary), M4 integrate to deep-overlap + gold-gate vs oracle (incl sampleB).
  Honest: realistic best case is a conditioning/transport lever, not a free closed form; the flow's class
  is the same GL3/W3 confluent-Garnier sigma (frontier). WS-D's abelian Ê warns the closing structure may
  govern only the abelian data, leaving sigma as the conserved IC (no free lunch).

Status: NOT executed yet (scoped only). Next action when authorized: run Milestone 1 (construct M_a from
the commuting partner; test closure via the d_a finite-difference-vs-RHS decisive check).

---

## 2026-06-03 — WS-AFLOW Milestone 1: a-flow does NOT close (OUTCOME iii) — coordinator-verified

The make-or-break of the central programme. Files: experiments/aflow_closure_test.py,
paper/working_sessions/ws_aflow_milestone1.md. Relaunched after the first background agent died in an idle container
reclaim.

VERDICT [numerically-supported]: OUTCOME (iii) — the a-flow of {P_mm,b} does NOT close into a
finite-dimensional geometric ODE. Both candidate closed forms FAIL the decisive test:
- Mechanism 2 (algebraic M_a / Q_k via the homological eq d_u Q_k - i[H,Q_k] = d_{a_k}H): the polynomial
  residual PLATEAUS (k=0..2, degrees 1->6: ~0.76 -> 0.10 -> 0.06 -> 0.06 -> 0.026..0.035; never collapses
  to machine zero as a true algebraic Q_k would, cf. the commuting-partner reconstruction at 1e-14). The
  WS-D null signature. => no algebraic Q_k.
- Mechanism 3 (KZ/Gaudin linear R_k=(d_a S)S^{-1}): fits at a point (1e-16) but does NOT transport
  (transport_err 0.40-0.54 = O(1); |R_k|~1500 = divergent Stark content). => not a flat geometric
  connection; carries propagator history.
- Ground truth: converged d{P_mm}/da_k = (+0.0620, -0.0352, -0.0267) (sum 0 = a-shift gauge), reproduced
  ONLY by the full-history Duhamel, by NEITHER history-free mechanism.

COORDINATOR VERIFICATION (M1 discipline, this is a load-bearing negative):
- Independently reproduced the STARK SUBTLETY: naive d_a of a finite-T propagator oscillates and does NOT
  converge (my crude Richardson-at-T=80 gave dP_mm/da=(-1.82,+1.68,+0.27), a point on the agent's
  oscillation -4.4/+1.4/+0.78/+0.70) -- confirms d_a and lim_T do not commute (Stark phase a_k u^2/2).
- Re-ran the agent script: reproduced the residual plateau and the R_k transport failure.
- (Oracle-path converged-derivative cross-check via num_S12.P_mm_oracle launched as extra confirmation.)
The closure logic is robust: a closing flow needs an algebraic Q_k (ruled out) OR a geometric R_k (ruled
out); both fail, so the only "closure" is the full Duhamel history = not a finite-dim ODE. Not a false
negative.

MEANING: the central principle (one a-independent geometric seed W organizes the whole 2-param family)
SURVIVES as organizing structure but yields NO computational shortcut. The deformation/transport method
(compute sigma cheap, transport to deep overlap) is a CLOSED DEAD-END. This is the Abelian ceiling (R1)
restated as a flow statement, consistent with WS-D and deformation_family_probe. Milestones 2-4 of
ws_aflow_derivation_plan.md are NOT viable as a closed transport flow; sigma-conservation is moot (no
closed abelian dressing flow to carry it). The irreducible transcendental sigma (= W's Stokes data, the
rank-3 GL3/W3 connection constant, WS-RH) remains, with no closed flow to it.

---

## 2026-06-02 — physics-reflection pass: add R17 + Conclusion (the honest reframe); coherence fixes

Final review pass over type1_lz_working_paper.tex.
- Added R17 [NS] (sec:skeleton): the commuting-family deformation flow does NOT close (Abelian ceiling
  as a flow statement) -- plateau (no algebraic Q_k), R_k non-transport (O(1)), converged
  d_a P_mm=(+0.062,-0.035,-0.027) reproduced only by full-history Duhamel; Stark d_a/lim_T subtlety.
  Coordinator-verified this session.
- Added Conclusion section "a governed account, and why the answer is the ODE": the honest reframe so the
  reader finds insight in the dead end. Corollary: for the exact value, structured numerical integration of
  the governing ODE is the most efficient method (evidence-backed). NOT naive integration -- the theory
  delivers the right (adiabatic) frame (cond 1e18 -> 1e-9), reduces to 2 numbers (BE+DS), gives the uniform
  formula for the generic regime, and certifies the residue irreducible. "What the dead ends bought":
  the governing account for an ungoverned constellation.
- REVIEW FINDINGS FIXED (M6 honesty): conclusion overstated "theorem-backed"/"provably no algebraic
  object"/"we proved there is none" -> softened to "evidence-backed (gold-gated, short of a closure
  theorem)". Coherence: intro central-object pnote presented the deformation programme as an open hope ->
  reconciled with R17 (the seed W organizes but does not compute).
- Cross-check: R15/R16/R17 numbers match their scripts (verified this session). Structure: braces balanced,
  all envs OK, R1-R17 present, no stray c=1 mislabels, no P_{2->2} leftover. No pdflatex in env (no compile).

VERDICT: paper is internally coherent, provenance-clean, evidence-tags honest. The session's full arc
(structural ceiling -> classification -> 3 deliverables -> skeleton/2-transcendentals -> central seed W ->
deformation no-go -> conclusion) is one consistent governing account. Remaining for "final": a real
pdflatex compile (no toolchain here) and the owed proofs already listed (R5/R11 region statements).

---

## 2026-06-02 — New direction scoped: the GEOMETRIC selection question (WS-GEOM) + Magnus/Feynman primer

User pivot beyond closed form: how does the Type-1 Hamiltonian select its slice of U(3)? Two facets
(both sigma-FREE): (A) the image region of Phi:(g,e,a)->S in reduced-U(3) (directed-cycle bias, boundary
= classifier locus); (B) the holonomy/topological selection (the a-independent seed W's holonomy on the
spectral-curve double cover; spinor signs delta_j).

- paper/working_sessions/ws_geom_scope.md: scoping roadmap (physics-derivation + physics-numerics). Part I: the adiabatic-W
  Magnus expansion to 2nd order -- order 0 = directed-cycle permutation (topological skeleton), Ω_1/Ω_2 =
  first hop + interference; validate term-by-term vs oracle in the generic regime; Feynman-graph dictionary
  (sheets=lines, W=vertices, dynamical phase=propagators; W a-independent => geometric vertices, elementary
  phase edges). Part II: map the image region with cheap data + a CAGED symbolic-regression invariant hunt
  (structure-detector ONLY, not surrogate/formula-finder). Milestones M1-M4; honest limits (no resum in
  deep overlap; sigma still the irreducible all-orders datum; ML may find nothing = first-class negative).
- paper/drafts/primer_magnus_feynman_geometry.md: pedagogical primer -- (i) Magnus = exp of nested-commutator
  series, unitary at each order, converges in weak coupling; (ii) -> Feynman graphs on 3 adiabatic sheets;
  (iii) -> geometric theory: U(3) selection = (topological permutation/spinor sector from the eigenbundle
  monodromy) x (smooth W-holonomy dressing), boundary = decoupling locus, sigma = the quarantined
  irreducible resummation.

Status: SCOPED, not executed. Re-approaches the demoted WKB/selector "common cover & spinor lift"
geometry with the new seed-W understanding. Entry point: WS-GEOM M1 (write Ω_1,Ω_2; order-0=directed
cycle; term-by-term oracle validation in the generic regime).

---

## 2026-06-02 — WS-GEOM M1 DONE & verified: the adiabatic-W Magnus skeleton (order 0 = node-selected directed cycle)

WS-GEOM Milestone 1 complete. Files: experiments/ws_geom_magnus.py, paper/working_sessions/ws_geom_m1.md.

ORDER 0 = DIRECTED CYCLE (coordinator-VERIFIED independently): the overlap-continued adiabatic-following
permutation |S_0|^2 = the directed 3-cycle [1,2,0], = the oracle's dominant-entry pattern (MATCH). KEY
NEW GEOMETRIC INSIGHT (from coordinator verification): the directed cycle is SELECTED BY CONTINUING THE
EIGENFRAME THROUGH THE UNIVERSAL NODE (R3) -- energy-sorting (ignoring the node) gives the WRONG
transposition (0 2); continuing through the exact crossing turns it into the 3-cycle. So the node (R3) is
what makes the leading U(3) structure a directed cycle. This links R3 <-> R16 geometrically.

ORDER-BY-ORDER GATE PASSED (adiabatic regime, agent table, T-averaged to cancel a Stark-tail oscillation):
e0>e1>e2 at every sample; all decrease as adiabaticity grows (gamma up): e.g. gamma 1.0/1.3/1.6/2.0 ->
e0 .325/.129/.065/.027, e1 .073/.0146/.0111/.0067, e2 .059/.0104/.0105/.0066. Order 1 cuts error ~5-10x
(leading non-adiabatic correction); order 2 (interference) adds a smaller gain to a ~5e-3 floor.
Unitarity exact at every order (exp(Omega) unitary to 4e-17).

TWO HONEST CORRECTIONS to my scope (now fixed in ws_geom_scope.md + primer .tex):
1. The directed cycle is the ADIABATIC order-0, NOT a small-delta object. Small delta = DIABATIC corner
   (order 0 = IDENTITY). Consistent with R16 (small delta -> near-identity; deep overlap -> cycle).
2. The Magnus convergence parameter Lambda=∫||Wtil||≈3≈pi is MARGINAL EVERYWHERE (nearly slope-invariant)
   = the FLOW-SIDE FACE OF R4/R5 permanent marginal overlap. So the series is an order-improving SKELETON,
   never a fast delta-series; stops improving in deepest overlap (delta≳5) where Lambda>pi and sigma lives.

Stark-tail subtlety (agent, honest): Wtil has a Fresnel tail e^{icu^2}/u^2 -> Omega_1,2 pick up
O(sin(cT^2)/T^2) endpoint terms; single-T errors oscillate. Resolved by T-averaging (NOT faked); full fix
= canonical log-T Stark subtraction (beyond M1).

NET: M1 establishes the geometric SKELETON of S in U(3): order 0 = node-selected directed-cycle
permutation (topological), orders 1-2 = the W-holonomy dressing (analytic), with marginal (R4/R5)
convergence so the resummation = sigma stays irreducible. Next: M2 (image-region map) / M3 (invariant hunt).

---

## 2026-06-03 — WS-GEOM M2 DONE & coordinator-verified: the image region (facet A)

M2 complete (agent reclaimed mid-session -- transcript severed, no completion notice -- but the full
deliverables were written: experiments/ws_geom_m2_image.py, paper/working_sessions/ws_geom_m2.md, 289 lines, all gates).
Coordinator INDEPENDENTLY VERIFIED the decisive geometric claims (4 cases, fast solve):
 - weak coupling -> {P_mm,b}=(0.996,0.002) = IDENTITY vertex (1,0);
 - strong coupling -> (0.034,1.000) = directed CYCLE vertex (0,1);
 - mid decouple (g_mid->0) -> P_mm=0.9998 (P_mm=1 edge);
 - extreme decouple (g_lo->0) -> b=0.0002 (b=0 edge).

RESULTS:
 - M2a PASS: image closure touches EXACTLY two of six Birkhoff vertices -- identity (diabatic/small-delta)
   and the node-selected directed 3-cycle (adiabatic/large-delta). Reverse-cycle and transposition corners
   NOT reached (min-dist ~0.27): the cycle ORIENTATION is node-fixed.
 - M2b: monotone log-in-action bias toward the cycle (Spearman(dist,delta)~-0.77), chi-modulated -- the
   same {action x shape} split as WS-O3/R11.
 - M2c CONFIRMED (headline): image boundary = the decoupling locus, with a SHARP EDGE-SELECTION rule --
   extreme-level decoupling -> b=0; middle-level decoupling -> P_mm=1; global bounds b>=0, P_mm<=1 with no
   violations. (Honest: the agent flagged that convex-hull distance on a sparse random cloud is misleading;
   the load-bearing claim is edge-selection + bounds, which I reproduced.)
 - M2d: effective dimension FULL -- PCA spectrum [1,.70,.65,.49,.20,.07,0,0,0]; the 3 exact zeros are the
   KNOWN elementary relations among {delta_ij,c_i} (incl sum c_i=0). NO new analytic constraint.

NET (honest "partial" outcome, as anticipated): the section's ANALYTIC content is exhausted by BE+DS+chi
(no hidden invariant); the genuinely NEW geometric content is TOPOLOGICAL -- two node-selected vertices +
the edge-selection boundary. IMPLICATION FOR M3: M2d pre-empts the analytic-invariant hunt (full dimension
=> none to find); M3 reduces to FORMALIZING the topological label (the two-vertex / spinor sector, the
node-fixed cycle orientation), not searching for an analytic relation.

---

## 2026-06-03 — M3 reframed as FORMALIZATION (M2d killed the invariant hunt) + geometric-theory section (R18,R19)

- paper/working_sessions/ws_geom_m3_scope.md: M3 is NO LONGER a symbolic-regression invariant hunt -- M2d showed FULL
  effective dimension (no hidden analytic invariant; only the known elementary relations). M3 reframed as a
  physics-DERIVATION of the topological statements M1+M2 found numerically, all sigma-free:
  T1 [keystone] node->directed-cycle theorem (eigenbundle monodromy on Sigma: extreme-swap transposition
  o node-swap = 3-cycle); T2 two-vertex reachability (diabatic limit=identity, adiabatic limit=cycle);
  T3 edge-selection boundary = decoupling locus (extreme->b=0, middle->P_mm=1, from the decoupling block
  structure, R9); T4 [stretch] the delta_j spinor double-cover monodromy. Milestones M3a-d; honest ladder
  (T1 closed = proven topological core; T4 may stay partial). Does NOT touch sigma, claims no closed form.

- type1_lz_working_paper.tex: new section "Geometric structure of the U(3) section" with:
  R18 [AD/NS] the adiabatic-W Magnus skeleton -- order 0 = node-selected directed cycle (coordinator-
  verified), marginally-convergent dressing (Lambda~pi <-> R4/R5), Feynman-graph reading (W vertices
  a-independent/geometric, phase edges elementary).
  R19 [NS] the image = two-vertex reachable set (identity + node-selected cycle; other 4 forbidden),
  boundary = decoupling locus edge-selected (extreme->b=0, middle->P_mm=1), FULL effective dimension =>
  no hidden analytic invariant; analytic content = BE+DS+chi.
  Plus "the geometric account": U(3) section = topological permutation/spinor sector (node-selected) x
  W-holonomy dressing, bounded by the decoupling locus; sole irreducible datum = sigma. Dep map updated
  (R18,R19); fixed a stray \ii -> \mathrm{i}.

Status: M3 scoped (not executed). The geometric theory's RIGID core is now stated (R18-R19) and its
proof is scoped (M3 T1-T4). sigma remains the explicitly-quarantined irreducible remainder.

---

## 2026-06-03 — M3 T1 CLOSED: the node->directed-cycle theorem [analytically-derived]

Files: paper/working_sessions/ws_geom_m3.md, experiments/ws_geom_m3_t1.py (verified: T1 True on canonical + 4 random).

THEOREM T1: the order-0 adiabatic-following permutation, CONTINUED THROUGH the unique real node (R3),
is a DIRECTED 3-CYCLE (not the energy-sorted extreme-swap transposition). Proof (3 steps + verification):
 1. Bare reordering = extreme-SLOPE swap tau_ext=(lo hi): E_k~a_k u reverses ordering between -+inf;
    N=3 full reversal = swap extremes, fix middle. Universal (any distinct slopes). [verified (2,1,0)]
 2. Node = adjacent-rank transposition tau_node: R3's unique real crossing is between energy-adjacent
    levels; continuation through it swaps that adjacent pair. [verified: continued differs by 1 adj swap]
 3. Composition: pi_cont = tau_ext o tau_node; two transpositions sharing one index = a 3-cycle
    ((a b)(a c)=(a c b)); (lo hi) o (adjacent) = directed 3-cycle. [verified pi_cont=(1,2,0)=tau_ext o tau_node]
 4. Orientation set by which adjacent pair the node swaps (spectral data); uniform on the family (verified,
    not separately proved -- minor owed item, ties to T4).

Caveat (M1): order-0 = the full ORACLE permutation only in the adiabatic regime; for diabatic/small-action
samples the oracle is near identity. T1 is about the order-0 skeleton.

EPS-ORDERING (answers a user question): the extreme-swap is the diabatic image of the asymptotic
energy-order REVERSAL; it swaps the extreme-SLOPE pair argsort(a)[0]<->argsort(a)[2] expressed in the
EPS-indexed basis. In the samples slopes were eps-monotone so it read (2,1,0); generically it swaps
eps-indices argsort(a)[0],argsort(a)[2]. This is exactly why strict eps-ordering matters: basis=eps (fixed),
roles=slope (which states swap), energy=derived & end-reversing (hence never a channel label).

IMPACT: T1 [keystone] closed => the geometric theory (R18) has a DERIVED topological core (the leading
U(3) permutation is a node-selected directed cycle). Remaining M3: T2 (two-vertex), T3 (edge-selection
boundary), T4 (delta_j double-cover orientation/label).

---

## 2026-06-03 — M3 T2 + T3 CLOSED (analytically-derived + verified)

Files: paper/working_sessions/ws_geom_m3.md (T2,T3 sections added), experiments/ws_geom_m3_t2t3.py (verified).

T2 (two-vertex reachability) [AD+verified]: the {P_mm,b} image touches EXACTLY identity (1,0) and the
node-selected directed cycle (0,1); the other 4 permutations map to (0,0) or (1,1) and are unreachable.
Proof: diabatic limit g->0 => H->uA diagonal => S=I => (1,0); adiabatic limit g->inf => follow eigenstates
=> node-continued cycle (T1) => (0,1); no other vertex is a limit point (extreme-swap is node-converted to
the cycle; reverse cycle excluded by orientation; rest are limit points of neither). Verified: gscale 0.03
->(1.000,0.000), 2.0->(0.024,1.000); (0,0)&(1,1) never approached (>0.5 always).

T3 (edge-selection boundary = decoupling locus, R9) [AD+verified]: middle decouple (g_mid->0) => P_mm->1
edge (middle is a spectator, survives w.p. 1); extreme decouple (g_lo or g_hi ->0) => b->0 edge (the
directed cycle needs all 3 links; removing an extreme breaks it). Global bounds b>=0,P_mm<=1 => genuine
boundary edges; the decoupling locus traces them, edge-selected by WHICH level decouples. Geometric
realization of R9 (elementary <=> a level decouples). Verified: g_mid->0 P_mm=1.0000; g_lo->0,g_hi->0 b=0.0000.

STATUS: M3 T1+T2+T3 CLOSED. The rigid (sigma-free) half of the U(3)-selection theory is now DERIVED:
node-selected directed cycle (T1) + two-vertex reachable set (T2) + decoupling-locus boundary (T3). Only
T4 (the delta_j spinor double-cover ORIENTATION label) remains (a stretch, possibly partial). sigma (the
analytic dressing) is the sole irreducible remainder. Paper geometric "account" updated to cite the derived
T1-T3.

---

## 2026-06-03 — M3 T4 PARTIAL: cycle-orientation Z2 label ESTABLISHED, selector OPEN (with a retraction)

Files: paper/working_sessions/ws_geom_m3.md (T4 section + T1 Step-4 correction), experiments/ws_geom_m3_t4.py, paper update.

ESTABLISHED [numerically-supported, deterministic]: the directed-cycle ORIENTATION is a genuine Z2
topological label = the eigenframe spinor/double-cover sector (delta_j=+-1). The dynamics-free
overlap-continuation of the order-0 eigenframe gives BOTH orientations over broad samples: REV (2,0,1) x7,
FWD (1,2,0) x4 (clean). So the U(3) section splits into TWO MIRROR SECTORS. Group rule (T1 Step 4):
node swaps (mid,hi) => FWD, (lo,mid) => REV.

RETRACTION (M1 honesty): orientation = -sign(u_*) is REFUTED. It held ~92% on seed 7 (N=13) but 33%
(worse than chance) on seed 11 (N=9), with misses at LARGE |u_*| -> a seed-dependent fluke, not a rule.
I had stated the 92% in an interim; retracted. No validated selector/dividing-locus is known.

CORRECTION to T1: T1 Step-4's "uniform orientation (lo->hi->mid->lo on the family)" was a 5-sample
artifact (seeds 0,3 happened to be FWD). Broad sampling flips it -> Z2. T1's CORE (order-0 = a directed
3-CYCLE, not a transposition; Steps 1-3) is UNAFFECTED and stands. Corrected in ws_geom_m3.md + paper.

METHOD note: reading orientation from the adiabatic-limit P is UNRELIABLE -- the limit is the deep-overlap
regime where P does not cleanly reach a vertex (clean-defect ~0.4 at scale 8). The deterministic
overlap-continuation is the reliable measurement and is what establishes the Z2.

OWED (T4): the orientation SELECTOR -- the spectral-flow combinatorics (which of lo-mid / hi-mid is the
real node vs the avoided complex branch point; only one real crossing exists though the asymptotic order
reverses) + the delta_j monodromy on the genus-0 curve Sigma. T4 is the partial piece flagged at outset.

STATUS: M3 T1-T3 closed (rigid skeleton derived); T4 partial (Z2 label established, selector open). The
sigma-free topological theory is essentially complete: node-selected directed 3-cycle (T1) + two-vertex
reachable set (T2) + decoupling-locus boundary (T3) + two Z2 mirror sectors (T4-label). Two open items
remain, neither touching the derived skeleton: sigma (analytic dressing) and the Z2 orientation selector.

---

## 2026-06-03 — T4 selector: systematic predictor screen => NOT a simple parameter sign (genuine monodromy object)

Pushed the orientation selector empirically (against the reliable adaptive-resolution overlap-continuation
orientation; node properly resolved with a dense grid near u_*). Screened candidate predictors; ALL fail
(~chance on clean samples, N~9-11, FWD/REV balanced):
 - -sign(u_*): 92% one seed, 33% another -> fluke (already retracted).
 - node nearer hi-mid vs lo-mid diabatic crossing: 4/11 (0.36).
 - sign(s2_mid,hi - s2_mid,lo): ~0.56;  sign(a_mid-(a_lo+a_hi)/2): ~0.56.
 - sign(E_*-(H0)_mid,mid): ~0.67;  sign(ucr_lm+ucr_mh): ~0.78 (small-N, did not hold up);
   sign(ucr_lm*ucr_mh): ~0.56.
CONCLUSION: the Z2 orientation selector is NOT a simple parameter sign -- it is the genuine spectral-flow /
delta_j-monodromy object on the genus-0 curve Sigma (which complex branch point reconnects which sheets;
which mid-crossing the real node realizes). This SHARPENS "selector open": it's not "we didn't look" but
"the natural elementary predictors are refuted." Remaining path = the analytic spectral-flow/monodromy
derivation (hard; beyond predictor screening). STOPPED the screen (diminishing returns).

Measurement note: orientation must be read by ADAPTIVE-RESOLUTION continuation (dense grid near u_*);
coarse continuation gives many ambiguous (o==0) reads, and adiabatic-limit-P reads are unreliable
(deep-overlap, P not at a vertex). With proper resolution both orientations occur (FWD 5, REV 6).

STATUS unchanged at the theory level: M3 T1-T3 closed (rigid skeleton derived); T4 = Z2 label established,
selector OPEN and now confirmed to be the genuine monodromy object (not elementary). Two open items remain
(sigma; the Z2 selector), neither touching the derived sigma-free skeleton.

---

## 2026-06-03 — T4 selector RESOLVED: orientation = energy-position of the real node  [near-proof]

BIG SWING (user: "take a big swing... don't give up at the first dead end"). The predictor screen had
correctly REFUTED every *local* candidate; the lesson was that the selector is *global*. Acted on that:
asked which adjacent energy-rank pair the **real node** (R3) degenerates, computed from the explicit
spectral curve.

RESULT (the selector):
    orientation = sign(E_3 - E_*) = sign(tr H(u_*) - 3 E_*),
  with (u_*, E_*) the unique real node = real double root of the E-discriminant Disc_E(u) of the
  spectral curve, E_* the doubly-degenerate eigenvalue (= R8 accessory v_*), E_3 = tr H(u_*) - 2E_* the
  spectator. LOWER pair degenerate (E_3 > E_*) <=> FWD (1,2,0); UPPER pair (E_3 < E_*) <=> REV (2,0,1).

DERIVATION (spectral flow): pi_cont = tau_ext o tau_k (T1). tau_ext = (lo hi) is fixed (BE / asymptotic
order-reversal), so the Z2 lives entirely in k = which adjacent energy-rank pair the node swaps. k is read
off the spectator's energy-position: lower pair degenerate <=> spectator above <=> E_3 > E_*. Both u_*, E_*
are explicit algebraic data (R3 real node; R8 v_*=E_*) => closed-form algebraic sign in (gamma,eps,a).

THE CRUX (why the screen missed it): the controlling point is the REAL node (an exact crossing = real root
of Disc_E), NOT the nearest COMPLEX branch point (the dominant avoided crossing). My first attempt used the
nearest complex branch point and failed exactly on the small-|u_*| REV cases (9/11) -- the avoided crossing
and the protected crossing can sit on different level-pairs. Switching to the real root fixed it.

EVIDENCE (near-proof): ws_geom_m3_t4_selector.py.
 - 154/154 = 100.00% vs deterministic overlap-continuation Z2 (clean random samples; FWD 85, REV 69; zero
   mismatch).
 - +33/33 on an independent seed block; canonical anchored (u_*=-0.2493, lower pair, FWD).
 - TWO INDEPENDENT ALGORITHMS agree on all 154: (a) Hermitian eigh of H(u_*); (b) companion-matrix roots of
   the characteristic polynomial p(E,u_*) -- no Hermitian solver at all.
 Ladder: analytically-derived (the spectral-flow group identity) + numerically-supported/near-proof (the
 identification of the real node as the controlling point, and uniqueness of the real node).

DEAD END THAT BECAME THE ANSWER: the refuted local-predictor screen (prior entry) was not wasted -- its
uniform failure DIAGNOSED the selector as a global object, which is exactly what pointed at the real-node
energy-position. Recorded as the methodological win.

STATUS: M3 T1-T2-T3-T4 ALL CLOSED. The sigma-free topological theory of the U(3) section is COMPLETE:
node-selected directed 3-cycle (T1) + two-vertex reachable set (T2) + decoupling-locus boundary (T3) +
Z2 mirror sectors with the real-node energy-position SELECTOR (T4). The SOLE remaining open item is sigma
(the irreducible analytic Fredholm/Widom dressing) -- which does not touch the topological skeleton.

Artifacts: experiments/ws_geom_m3_t4_selector.py (new); paper/working_sessions/ws_geom_m3.md (T4 RESOLVED + derivation);
paper/drafts/type1_lz_working_paper.tex (R18/R19 geometric account updated).

---

## 2026-06-03 — T4 selector SHARPENED to a parity law + M4 synthesis (Geometric Selection Theory)

While building M4 (the synthesis of M1-M3), tested the selector's predictions and FOUND a much cleaner
form than the energy-position sign.

PREDICTION TEST (physics-numerics discipline): conjectured the orientation flips through a TRIPLE
degeneracy (E_3=E_* at the node). REFUTED by a 1-param scan (a_mid through a_lo): |E_3-E_*| is LARGE on
both sides; the node instead runs u_*->+-inf at a SLOPE COLLISION. The flip is a parity-change at an
ordering wall, not a triple point. This refutation POINTED at the right answer:

PARITY LAW (the clean selector): orientation = sgn(sigma), the PARITY of the permutation sigma mapping the
eps-order to the slope-order (slope-rank at each eps-index). EVEN => FWD (1,2,0); ODD => REV (2,0,1). It is
the sign character S3->{+-1} (the eigenframe Z2/spinor sector). PURELY combinatorial:
  - INDEPENDENT of gamma and of all spacings (depends only on the relative ORDER of eps & slopes);
  - chamber-constant: no interior flip in eps_1 or gamma_mid scans;
  - all 6 slope-to-eps assignments give exactly the even/odd split;
  - gamma-INDEPENDENT over 80 random draws (odd stays REV, even stays FWD);
  - matches cont_orient on 119/119 clean samples, AND equals the energy-position selector on every sample.
DERIVATION: each adjacent slope-transposition swaps which adjacent pair the real node degenerates (node ->
inf at the slope collision), flipping the orientation; so orientation = sgn(sigma)*orientation(identity),
and the eps-monotone identity = FWD => orientation = sgn(sigma). [AD + near-proof.]

So the selector now has TWO equivalent closed forms: (1) combinatorial parity sgn(sigma) [clean, gamma-free];
(2) spectral sign(tr H(u_*)-3E_*) [the geometric realization on the real node]. Why the local screen failed
is now obvious: a permutation PARITY is a discrete global invariant, invisible to any continuous local probe.

M4 SYNTHESIS (ws_geom_m4.md): the Geometric Selection Principle. S = (rigid topological sector fixed by Sigma
+ its real node) x (W-holonomy dressing). Shape fixed by geometry: node-selected directed 3-cycle (T1),
two-vertex reachable set (T2), decoupling-locus boundary (T3), parity-selected orientation (T4). Position =
holonomy of the single a-independent seed W (M1), summing to ONE irreducible Stokes constant sigma. Geometry
selects the section; dynamics only dresses it. No further analytic invariant (M2d). Dependency map +
predictions (P1 refuted-kept, P2 parity supported, P3/P4 N>3 conjectures, P5 sigma-uniqueness consistency).

STATUS: M1-M4 COMPLETE. The sigma-free Geometric Selection Theory of the U(3) section is finished. Sole open
frontier = sigma (the Fredholm/Widom constant of the rank-3 W3 c=2 problem). Roadmap: (1) sigma as a Fredholm
determinant (WS-O2b); (2) N=4 test of P3/P4; (3) sigma-free corollaries already cashed (BE + uniform law +
two-component IP integration).

Artifacts: paper/working_sessions/ws_geom_m4.md (new synthesis); paper/working_sessions/ws_geom_m3.md T4 (parity law + equivalence + derivation);
experiments/ws_geom_m3_t4_selector.py (perm_parity_selector + two-form check); paper/drafts/type1_lz_working_paper.tex
(R20 Geometric Selection Principle + selector parity form in the geometric account).

---

## 2026-06-03 — M5 (Feynman-graph calculator + validity domain) and M6 tournament (improve Dykhne-Stuckelberg?)

M5 [physics-numerics + derivation]: the adiabatic-W Magnus/Feynman series is a practical CLOSED-FORM
calculator for S in the adiabatic window, with a MAPPED validity domain (ws_geom_m5_graphs.py,
ws_geom_m5.md). Findings:
 - VALIDITY (Task A): error is NOT monotone in Lambda (which barely moves, ~pi everywhere); there is a
   SWEET SPOT at gamma-scale~1.8 (e2~1.3e-3) and it degrades on BOTH sides (too diabatic: e0 large; too
   strong: order-2 a wash e2/e1~1, even worsening e0). The bare truncation CANNOT be driven below ~1e-3.
 - DDP/uniform-law BRIDGE (Task B): theta_ij'=E_i-E_j has NO real zero (adiabatic levels never cross) ->
   graph integrals dominated by COMPLEX turning points -> steepest descent = DDP. Order-1 -> single-crossing
   survivals p_x=exp(-2pi delta_x) (BE extremes); order-2 commutator -> Stuckelberg interference -> the
   uniform law's 2 sqrt(A0 Aret) cosPhi. The graphs ARE the microscopic derivation of the uniform law.
 - BOUNDARY (Task C): at Lambda~pi the graphs do NOT resum; the remainder is the irreducible sigma; the
   practical engine beyond the domain is the exact two-component adiabatic-IP integration.

M6 question (user): can the graph machinery be ADAPTED to IMPROVE Dykhne-Stuckelberg (not just truncated)?
TOURNAMENT (ws_geom_m6_tournament.md) over H1-H6 + evolved H7=H1(+)H4:
 - Leaderboard: H1 (first-principles Stuckelberg phase) wins; H7=H1+H4 (graph-derived uniform phase) is the
   target architecture; H4 (Weber uniformization), H5 (Borel/resurgence of sigma), H6 (graph residual), H2
   (Zhu-Nakamura prefactors), H3 (multi-path, opaque) follow.
 - DECISIVE FIRST-CUT for H1 (ws_geom_m6_h1_phase.py) [numerically-suggestive]: extract true cosPhi =
   (P_mm-A0-Aret)/2sqrt(A0 Aret) in the interference regime (N=17); at EQUAL parameter count (1 const each),
   the first-principles phase phi_ml+phi_stokes (mid-lo gap integral between the two diabatic crossings +
   the Stokes phase ALREADY in ws_o3_uniform but unused) BEATS the fitted K: mean|cos-cos_true| 0.161 vs
   0.189. DIRECTIONALLY H1 wins.
 - HONEST MAGNITUDE: cosPhi is near-quadrature/small in the interference regime, so net P_mm improvement
   ~ +0.006 (SUB-1%). VERDICT: H1 is a STRUCTURAL win (a parameter-free uniform phase of similar ~1%
   accuracy, removing the fitted K), NOT an accuracy breakthrough. The residual is the sigma-adjacent
   non-period content -- consistent with the hard ceiling (the open content is provably not a period).
 - ANSWER to "can it be done?": YES but bounded -- the machinery improves the *period part* of the phase
   (removes the fitted constant) and gives a few orders of systematic correction before the marginal-Lambda
   wall and sigma take over; it does NOT yield an exact closed form. Recommended M6 build: derive phi_dyn
   exactly + the offset (parameter-free law) and extend with H4 (Weber uniformization); ESCALATION FLAG: stop
   if phi_dyn turns into the opaque exact-WKB matrix product.

Artifacts: experiments/ws_geom_m5_graphs.py, experiments/ws_geom_m6_h1_phase.py; paper/working_sessions/ws_geom_m5.md,
paper/working_sessions/ws_geom_m6_tournament.md.

---

## 2026-06-03 — M6/H1 CORRECTION: the "structural win" was a sample-size artifact; H1 FAILS (robust negative)

Following the user's "build full M6 (H1 then H7)": pursued H1 to a parameter-free phase and ran a ROBUSTNESS
check. The earlier suggestive-positive is RETRACTED.

 - Earlier N=17 sample: parameter-free phi_dyn+3*stokes ~ 0.189 = fitted-K 0.189; best-offset 0.161 < fitted.
   Looked like a structural win (parameter-free at equal/better accuracy).
 - ROBUST N=42 sample (ws_geom_m6_h1_phase.py, larger scale x shift grid): fitted-K 0.090; parameter-free
   first-principles 0.281 (~3x WORSE); first-principles + best offset 0.179 (~2x worse, even WITH a fitted
   offset). The fitted K DECISIVELY BEATS the elementary first-principles phase.

VERDICT (corrected, numerically-supported NEGATIVE): H1 FAILS. The elementary turning-point phase (mid-lo
adiabatic gap action + the three crossings' Stokes phases) does NOT reproduce cosPhi; the fitted constant K
is absorbing genuinely NON-elementary, sigma-adjacent (non-period) content. Beating K requires the EXACT
turning-point connection = sigma (the Fredholm/Widom constant of the rank-3 W3 c=2 problem) = the opaque
exact-WKB matrix product the program rejected.

H7/H4 ESCALATION (do not build): the H4 Weber/parabolic-cylinder uniformization IS that exact connection; in
deep overlap it is MORE sigma-dominated. The escalation flag set at M6 scoping has fired -- there is no
non-opaque graph adaptation that beats the fitted K. Stopped rather than grind out a Weber computation of
sigma. (Honest dead-end, preserved.)

VALUE (sharpens sigma): sigma is NOT a small remainder on a good semiclassical estimate -- it is LOAD-BEARING
even in the interference PHASE, large enough that a single fitted constant outperforms first-principles
elementary semiclassics. This explains why no elementary closed form for the interference exists, and recasts
the fitted K in the uniform law as the cheapest honest stand-in for sigma-content (not a removable blemish).

ANSWER to the user's question "can the graph machinery improve Dykhne-Stueckelberg?": NO, not with the
elementary (non-opaque) machinery. Truncation is marginal (M5); the elementary first-principles phase is
worse than fitted-K (this M6); beating K requires sigma itself. The graph engine remains valuable as the
microscopic DERIVATION of the uniform law (M5 Task B) and a transparent adiabatic-window calculator, but it
does not sharpen the estimate. M6 closed as a first-class NEGATIVE.

Process note: the N=17 false positive was caught precisely by the convergence/robustness discipline
(physics-numerics) -- a good example of why single-sample "wins" must be re-tested before promotion.

---

## 2026-06-03 — Omega_3 directly computed: NO significant refinement (closes a real gap)

User flagged that the marginality verdict (M5) only used Omega_0,1,2 -- we never formed Omega_3. Closed the
gap: implemented Omega_3 = 1/6 int_{t1>t2>t3}([A1,[A2,A3]]+[A3,[A2,A1]]), A=-Wtil, via nested cumulative
integrals (ws_geom_m5_omega3.py).
 - CORRECTNESS GATE: Om1,Om2 reproduce magnus_terms() to 0.0; ||Om3+Om3^H||=0 (anti-Hermitian) => correct.
 - TERM NORMS shrink geometrically: canonical ||Om1,2,3||=1.03,0.51,0.12; sweet spot sc=1.8: 0.295,0.041,
   0.004 (ratio ~0.1). So Om3 is genuinely small, not a hidden large term.
 - REFINEMENT: e3 ~ e2 across the sweep, raw AND cleanly T-averaged. Sweet spot sc=1.8 (T-avg vs 1e-13
   oracle): e1,e2,e3 = 1.64e-3, 1.45e-3, 1.45e-3. ||Om3||=0.004 exceeds the 1.45e-3 residual yet does not
   reduce it.
 - VERDICT: Omega_3 gives NO significant refinement. The residual floor is OUTSIDE the perturbative Magnus
   tower (non-perturbative sigma-scale remainder + finite-T endpoint tail), not a missing finite order;
   Om4... are smaller still and cannot reach it. The M5 marginality verdict HOLDS and is now unambiguous.
 - Side product (pedagogy): documented WHY steepest descent through the COMPLEX turning points gives BE/DDP
   (imaginary action to u_c = 2*pi*delta_ij) and Stueckelberg (two-saddle interference = the uniform-law
   cosPhi); see where_we_are.tex sec 3.

---

## 2026-06-03 — H3 reopened (M7): multi-path / transfer-matrix sum -- NEGATIVE for Type-1, but DEMYSTIFIED

User (correctly) flagged that H3 (Sinitsyn-style multi-path interference) was deprioritized in the M6
tournament but never tested, and that the Omega_3 negative does NOT refute it (multi-path = a semiclassical
TRAJECTORY sum, not the perturbative Magnus tower -- a distinct object). Reopened it properly (test + explain
+ literature).

LITERATURE (WebSearch, grounded): Sinitsyn et al. obtain EXACT transition probabilities via path
interference for SPECIAL solvable MLZ models (4-state arXiv:..., 6-state arXiv:1501.06083), where
trajectories connecting the same endpoints interfere coherently -- explicitly NOT the incoherent LZ product.
The independent-crossing approximation factorizes only "if path interference and accidental crossings are
absent." Type-1 N=3 general is not such a solvable model (P_mm is sigma-transcendental).

TEST (ws_geom_m7_multipath.py): built the independent-crossing transfer-matrix (ICTM) cleanly -- each
crossing = exact isolated 2-level integration (DOP853), spectator gets its diabatic phase, product in
u-order (convention-safe; Stokes phases automatic; coherent multi-path interference included).
 - Type-1 P_mm: ICTM is MUCH WORSE than the calibrated 2-path uniform law. |ICTM-ref| vs |uni-ref|:
   sc=0.5: 0.101 vs 0.001; sc=0.8: 0.342 vs 0.018; sc=1.2: 0.049 vs 0.073; sc=1.8: 0.078 vs 0.039.
   Uniform wins 3/4; ICTM edges only at sc=1.2. Mean ICTM err ~0.14 >> uniform ~0.03.
 - Separable generic toy (b_i free): ICTM error DECREASES with crossing separation (0.16 -> 0.043 as
   +-5 -> +-40) but slowly; no clean exact-in-the-limit validation obtained -> ICTM numbers indicative,
   not definitive (finite-separation ICA residual +/- phase-bookkeeping).

STRUCTURAL RESULT (the real contribution): Type-1's Cauchy weld (R4/R5) PREVENTS crossing separation.
Since H0 ~ gamma^2, small coupling sends all (H0)_ii -> 0, COLLAPSING every crossing to u=0 (maximal
overlap); large coupling widens them. There is NO limit where Type-1 crossings are isolated, so the
independent-crossing/multi-path factorization can NEVER reach its exact regime; its error is the overlap
content = sigma. This DEMYSTIFIES Sinitsyn's effectiveness (his models permit separation or exactly-summed
interference) and bounds it: Type-1 is the structural worst case for multi-path.

VERDICT: H3 = NEGATIVE for Type-1 N=3 (multi-path is not a refinement; worse than the uniform law), but
with a clean explanation. Consistent with every prior result: the uniform law is already the best cheap
estimate, and the irreducible remainder under permanent overlap is sigma. M7 closed.

Process note: the M6 tournament's low ranking of H3 was a too-quick dismissal (conflated with the Magnus
tower); reopening on the user's prompt produced the structural R4 insight. Good example of scientist-in-the-
loop catching an under-examined branch.

---

## 2026-06-03 — PARALLEL ANALYSIS: perturb Type-1, remove the protected crossing — PERTURBATIVE

New goal (user): perturb the Type-1 matrix to remove the crossing; is the impact on the LZ amplitudes
perturbative, or is there a non-perturbative component?

FRAMING (engaging the user's claims): the protected real node (R3) is a codim-2 degeneracy stabilized by
integrability. For real-symmetric H(u)=H0+uA a degeneracy is codim 2, so a generic 1-parameter family has
NO real crossing -- a real crossing requires either fine-tuning (unprotected, perturbs away) or a structure
forcing it (commuting => Type-1, protected). So "remove the crossing" = break the commuting structure; the
crossing and integrability are the same feature. (User's "crossing => Type-1" holds for PROTECTED crossings.)

EXPERIMENT (ws_pert_crossing.py): H_eps(u) = H0 + eps*V + u*A, V fixed generic symmetric (||V||=1, seed 7),
A unchanged (so it stays an MLZ model). Sweep eps in [1e-3, 0.2]; measure min_gap(eps), dP_mm, db,
d(extreme survivals), scaling exponents. Reference: T-averaged DOP853.

RESULT (numerically-supported): PERTURBATIVE.
 - gap ∝ eps^0.98 (≈0.42 eps): the protected node opens LINEARLY.
 - |dP_mm| ∝ eps^0.96, |d_surv_lo| ∝ eps^1.06, |db| ~ eps^1 (noisier) -- ALL exponents ~1.
 - The slope dP/eps is CONSTANT in eps (the eps^1 signature) -- robust across two numerical configs.
 - CAVEAT: the precise dP_mm SLOPE COEFFICIENT is endpoint-tail-limited (~1e-4): +0.11 vs -0.022 across
   configs (sign not resolved). The EXPONENT ~1 is robust; the coefficient is not. Does not affect the
   verdict.
 - No fractional power, no eps^2 log eps at leading order, no exp(-c/eps). The crossing-specific new LZ
   channel (diabatic traversal of the opened gap) is exp(-c eps^2) -> O(eps^2), subdominant + analytic.

WHY (analytic): S(eps) is differentiable at eps=0 -- dS/deps = -i int U0^dag V U0 du is FINITE despite the
degeneracy, because the degeneracy is a single point (measure zero) in the u-integral. So first-order
perturbation theory for the S-matrix is NON-SINGULAR -> O(eps).

CONSISTENCY: since A is unchanged, H_eps is still an MLZ model, so the EXTREME survivals remain exactly
Brundobler-Elser (with eps-shifted couplings) -- which is why d_surv ~ eps is clean and analytic. Only
integrability + the protected crossing are removed, not the LZ structure.

INTERPRETATION (the answer): the crossing is TOPOLOGICAL for the eigenframe/integrability (it carries the
node-swap T1, the parity selector T4, and is where sigma is defined), but its removal is an ANALYTIC
perturbation of the OBSERVABLE. The non-perturbative object sigma lives AT the integrable point, not in the
response to leaving it -- sigma is a property of the limit, not of the perturbation.

Open follow-ups (not pursued): (i) isolate the crossing-specific O(eps^2) channel from the generic O(eps)
global response (localized perturbation in the node's 2D subspace); (ii) higher-accuracy propagator
(Richardson) to pin the dP_mm coefficient and probe for a subleading non-analytic (eps^2 log eps) term.

---

## 2026-06-03 — Perturbation analysis, broadened: multi-direction + node-isolated (ws_pert_crossing_isolate.py)

Followed up the perturbation goal with (1) multiple V directions and (2) a node-isolated comparison, per
user request.

PART 1 (direction-genericity): seeds 1-4, generic dense symmetric eps*V on H0. |dP_mm| ~ eps^(~1) each
(exponents scatter 0.65-1.76 in a lean single-T / narrow-eps config -- numerical noise, NOT physical; none
non-perturbative). With the clean seed-7 run (eps^1.0 to 1e-3), the PERTURBATIVE scaling is direction-generic.

PART 2 (isolate the crossing): built V_open (generic -> opens node gap O(eps)) vs V_keep (same M but with the
node's 2x2 block SCALARIZED in the degenerate eigenbasis at u_*, so the pair stays degenerate to first order
-> crossing preserved, gap opens only O(eps^2)), equal norm. Result:
 - gap_keep stays at the grid floor (~4e-4 ≈ crossing preserved) while gap_open grows -> construction works.
 - ROBUST: dPmm_open ~= dPmm_keep (differ ~30-40%; their DIFFERENCE is ~10x smaller than either). So WHETHER
   OR NOT the crossing is opened, dP_mm is nearly the same -- the observable responds to the GENERIC
   integrability-breaking; the crossing-removal per se is a SUBDOMINANT contribution.
 - The predicted O(eps^2) scaling of the isolated crossing channel is NOT cleanly resolved (diff exponent
   ~1.1-1.4, crossover+noise-limited). Honest caveat; would need T-averaging + finer gap + smaller eps.

NET (answers the broadening): the impact of removing the protected crossing is PERTURBATIVE, direction-generic,
AND the crossing-removal is not even the dominant part of a generic perturbation's effect on the observable --
the bulk is ordinary first-order integrability-breaking. Reinforces: sigma is a property of the integrable
LIMIT, not of the response to leaving it; the crossing's topological role is in the eigenframe/integrability,
not a special driver of the observable.

Compute note: propP (DOP853 over u in [-T,T] with oscillatory Stark phases) is the bottleneck; several runs
timed out at full accuracy. The clean exponents live in ws_pert_crossing.py (seed 7, T-averaged); this
broadened run is qualitative (direction-genericity + crossing-subdominance), honestly noise-limited on the
finer exponents.

---

## 2026-06-03 — CONCEPTUAL NOTE: why "non-perturbative in the eigenbasis" does NOT imply non-perturbative P

Resolving an intuition that the perturbation analysis usefully undermined. The intuition: because opening the
crossing is dramatic in the EIGENFRAME evolution, that drama should surface in the LZ probabilities as a
non-perturbative effect. It does not. Two clean distinctions explain why:

1. EIGENFRAME-SINGULAR vs DYNAMICALLY-SINGULAR (a gauge artifact).
   The intuition's correct half: opening a gap eps reorganizes the adiabatic frame violently -- the
   eigenvectors (ill-defined at the exact degeneracy) rotate over a u-window of width ~eps, so the derivative
   coupling W_ij = <phi_i|d_u phi_j> ~ 1/gap spikes to height ~1/eps; the eigenvector ordering and the
   node-swap (T1) reorganize. BUT this is the price of the ADIABATIC frame being a bad coordinate system at a
   degeneracy. The observable lives in the DIABATIC frame, where the generator H0+eps*V+uA is smooth and
   LINEAR in eps over a finite u-measure, so S(eps)=T exp(-i int H) is analytic (dS/deps finite). The 1/eps
   spike is a COORDINATE singularity, not a dynamical one; it is tall-but-narrow, so int W stays O(1) and the
   integrated effect is bounded.

2. GAP-eps (diabatic, perturbative) vs RATE-eps (adiabatic, non-perturbative).
   Non-perturbative LZ -- exp(-c/eps) -- is the regime where eps sets the SWEEP RATE (slow/adiabatic passage).
   Here eps sets the GAP. Opening a crossing is the SMALL-GAP / DIABATIC direction of LZ, where stay-diabatic
   = exp(-c eps^2) ~ 1 - O(eps^2) (analytic) and the new adiabatic channel is 1 - exp(-c eps^2) = O(eps^2).
   So LZ DOES capture the crossing-opening -- and tells you it is perturbative. The intuition attached
   "non-perturbative" to the wrong limit of LZ (slow passage), whereas opening a gap is the fast/diabatic limit.

THROUGH-LINE (recurring project lesson): the ADIABATIC frame is where the topological/structural drama lives
(node-swap + parity T1/T4, the 1/gap coupling, the Magnus cumulants), and the OBSERVABLE is largely blind to
it. The genuine non-perturbative content of THIS problem (sigma) sits at the integrable point and in the
permanent-overlap regime (Lambda~pi), NOT in the gauge-singular response to opening one crossing. A
crossing-removing perturbation stresses the FRAME, not the dynamics.

---

## 2026-06-03 — /goal: the fiber bundle over the ring -- how a WITHIN-RING perturbation perturbs the amplitude

New goal: Type-1 is a commuting ring; for each member there is a transition amplitude (mod phase). No
elementary closed form (sigma), but characterize the fiber bundle: how does a within-ring perturbation
(fixed gamma,eps; vary slopes a) perturb the amplitude? Non-transcendentally?

Builds on R17/WS-AFLOW (the autonomous a-flow does NOT close into a finite elementary ODE) and R15/R16
(P affine in two transcendentals {P_mm,b}), R11/R16 (base invariants delta_lo,delta_hi,chi).

STRUCTURE (the bundle): base = ring = a mod shift (2-dim with fixed gamma,eps); the map a -> (delta_lo,
delta_hi,chi) is ELEMENTARY. Fiber = {P_mm,b}. Connection = abelian-flat (formal monodromy, elementary) +
sigma-twist.

EXPERIMENT (ws_ring_bundle.py): finite within-ring steps (|da|=0.1), T-averaged DOP853.
 NUMERICAL LESSON (logged honestly): a first tiny-h (1e-3) central-difference gradient gave GARBAGE
 (~ -2.4) because the single-T propagator's endpoint Fresnel tail varies with a and is amplified by 1/h.
 The direct finite-step (h=0.05) gave the true ~ -0.08. RETRACTED the noisy gradient numbers
 (cos -0.994 / 26x). Fix: T-averaged P + finite steps.
 CORRECTED FINDINGS:
  * SHIFT-NULL exact: a->a+(1,1,1) gives dP_mm = -0.0000 (global-phase invariance).
  * RANK-2 TANGENT: any within-ring dP reproduced from {dP[lo,lo],dP[hi,hi],dP_mm,db} via
    double-stochasticity to 1e-6 -> the tangent map is ELEMENTARY except for the rank-2 transcendental
    {dP_mm,db}. (Tangent-level R15.)
  * ELEMENTARY CAPTURE: the elementary uniform law (R12) TRACKS the exact within-ring response of P_mm:
    cosine +0.92, magnitude ratio ~1.3 (per-direction 1.50,1.01,1.08 in substantial directions). It flips
    sign only in near-null directions (true response ~1e-3, swamped by the law's ~2% value-error).

ANSWER (non-transcendental characterization): YES, mostly. (i) The base deformation a->(delta,chi) is
elementary. (ii) The tangent map is elementary EXCEPT for the rank-2 pair {dP_mm,db}. (iii) Even that pair's
response is ~elementary: the uniform law captures the within-ring response direction (cos 0.92) and magnitude
(~30%); the transcendental residual is SUBDOMINANT, dominant only where the elementary response nearly
vanishes. So the perturbation is mostly characterizable non-transcendentally, with a subdominant sigma-residual.
No EXACT elementary flow (R17 stands), but a good APPROXIMATE one -- and the bundle's whole non-elementary
content is rank-2.

Open: (a) decompose the response into the parameter-free incoherent part (A0+Aret, elementary) vs the
interference (cosPhi, transcendental) to see which drives the ~0.92 capture; (b) characterize the near-null
locus where the elementary response vanishes (likely an extremum ridge of P_mm on the base); (c) holonomicity
of P_mm on the base (is it D-finite/period-like or genuinely wild/Painleve -- predicted wild from W_3 class).

---

## 2026-06-03 — THEORY-BUILDING: the commuting family as a section; the bundle over the transport cover

User broadened out (geometry/topology, NOT solving amplitudes): (1) theory of sections of the vector bundle
over the transport cover; (2) selecting a commuting family = a choice of section; (3) how this shapes the
space of amplitudes over the family. Built paper/working_sessions/ws_geom_bundle_theory.md.

GEOMETRIC HOME (lit-confirmed: Boalch/Biquard-Boalch/Sabbah wild nonabelian Hodge; Hitchin; BNR; JMU
isomonodromy): the LZ problem is a WILD HITCHIN SYSTEM on P^1 with one rank-2 irregular point at u=inf.
Dictionary: H(u)=H0+uA <-> irregular connection; E_i(u) <-> genus-0 spectral curve Sigma (R2); real crossing
<-> node of Sigma (R3); eigenvectors <-> spectral line bundle L (BNR) / eigenbundle E with a-independent
connection W (M4); commuting ring <-> Hitchin base/abelian algebra (R1); member a <-> Hitchin-section point
(phase differential on Sigma); amplitude S/P <-> Stokes data = point of the wild character variety (dim 6,
R9/R11); map a->S <-> wild Riemann-Hilbert restricted to the Hitchin section; sigma <-> RH transcendence
(non-algebraicity on the Hitchin base; R11b/R17). "Nonabelian Hodge encodes two nonlinear systems: integrable
(the family) + isomonodromy (the amplitude flow)" -- exactly our two sides.

THE THREE BULLETS:
 (1) Transport cover = Sigma + its Z2 spinor double cover (eigenframe single-valued on Sigma, phases on the
     double cover). Sections = global eigen-frames, constrained by: abelian formal monodromy (F1), node
     gluing (R3), spinor Z2 sign. Canonical section = adiabatic frame Phi; bundle marginally non-flat
     (Lambda~pi); holonomy = the one seed W = the whole transcendence (M4).
 (2) A commuting family is the SECTION (Hitchin section / polarization) that makes the multivalued eigenframe
     coherent with an a-INDEPENDENT W; Type-1 is one such section; a-independence of W = the section
     trivializes the eigenbundle connection over the base. Abelian ceiling (R1) = the section is abelian =>
     it fixes the BASE (spectral data), never the FIBER (Stokes data).
 (3) Space of amplitudes = wild-RH image of the Hitchin section = the bundle over the ring: base (delta,chi)
     elementary; fiber {P_mm,b} rank-2 transcendental; reachable = two-vertex region + decoupling boundary
     (M2/T2/T3); orientation Z2 = sgn(eps<->slope) (T4, the only section-dependent topology); connection =
     abelian-flat + rank-2 sigma-twist (ws_ring_bundle); transcendence = sigma = RH non-algebraicity.

UNIFYING PRINCIPLE (bundle form of Geometric Selection): one intrinsic datum (E,W on the transport cover);
a family = the section trivializing W; the amplitude = the section's wild-RH image; the whole non-elementary
content = one isomonodromy constant sigma. Geometry selects which slice of the wild cv the family sweeps;
it never computes the slice.

INTRINSIC vs CHOSEN vs IRREDUCIBLE ledger: intrinsic (Sigma, node, irregular type, formal monodromy, cycle
structure, W, wild-cv ambient); member-dependent (orientation Z2, position {P_mm,b}, periods delta,chi);
irreducible (sigma = a property of the RH map).

STATUS: mostly a REFRAMING (all quantitative claims are established results re-expressed); the new content
is the organizing identification (family<->Hitchin section, amplitude<->wild RH, sigma<->RH transcendence),
tagged 'framing' throughout -- the natural, consistent math home, NOT theorems proved here. Roadmap (S9):
(1) pin the wild-RH dictionary against Boalch/Sabbah; (2) identify the BNR line bundle (why is W a-indep?);
(3) two-vertex set as a wild-cv stratum; (4) N=4 test of the structural predictions. No amplitude solved --
the value is the map of where the geometry lives and what is intrinsic vs chosen vs irreducible.

Refs: Hitchin 1987; Beauville-Narasimhan-Ramanan 1989; Jimbo-Miwa-Ueno 1981; Sabbah; Biquard-Boalch; Boalch
(wild character varieties / meromorphic Hitchin systems, e.g. arXiv:1703.10376, 1512.08091, 1203.6607).

---

## 2026-06-03 — Tame neighbor built: the tanh (Demkov-Kunike) analog vs the linear LZ sweep

Demonstration (ws_tanh_demkov_kunike.py) that LZ's transcendental sigma is the price of the UNBOUNDED linear
ramp (the rank-2 irregular point at u=inf), not of the level structure. Compared, same canonical (gamma,eps,a):
LINEAR (LZ, wild) H=H0+u*A  vs  TANH (DK, tame) H=H0+tanh(u)*A.
 (A) tanh eigenvalues SATURATE to eig(H0+-A) by u~5 (regular singular); linear diverge ~a_i*u (irregular).
     Dynamical phase scaling: theta_tanh ~ T^1.05 (e^{i*lambda*u}, regular) vs theta_linear ~ T^1.99
     (e^{i*lambda*u^2}, Fresnel/irregular). The exact regular-vs-irregular signature, measured.
 (B) Channel amplitude (read in the eigenbasis of H(+-inf), NOT diabatic): P_mm^tanh converges
     EXPONENTIALLY to ~1e-12 (increments -2.7e-5,-1.1e-7,+2.2e-8,+8e-12 at T=6,8,12,16); P_mm^linear only
     oscillates with a slow ~1e-4 Fresnel tail. Tame=well-conditioned, wild=Stokes-tailed.
     [Subtlety fixed: the tanh diabatic-basis P oscillates indefinitely because H(+-inf)=H0+-A is
      constant-but-not-diagonal -- itself the regular-singular e^{i*lambda*u} hallmark; the physical
      amplitude is the H(+-inf)-eigenbasis projection.]
 (C) Closed-form anchor: sech/Rosen-Zener P = sin^2(Omega*pi/2) reproduced to 1e-12 -> tame => exact
     elementary amplitude.
 (D) P_mm^tanh = 0.20383 (tame, trivially converged) vs P_mm^linear = 0.21472 (= sigma, wild).

CONCLUSION (near-proof of the structural claim): the wildness/sigma is the rank-2 IRREGULAR point from the
unbounded ramp. Tame neighbors (bounded sweep -> Fuchsian/hypergeometric; or finite-pole Gaudin/KZ) have
closed-form/holonomic amplitudes. Second axis: N=2 LZ is elementary even when wild; N>=3 (W3) is what makes
the connection constant transcendental. Folded into ws_geom_bundle_theory.md S8b (tame neighbors).

---

## 2026-06-03 — CORRECTION: Demkov-Kunike does NOT restore closed form for N=3 (rigidity)

User asked: does DK place us unambiguously back in closed-form territory? Answer: NO -- and this corrects an
overgenerous phrasing ("hypergeometric-class") in the prior tanh writeup.

TAME != CLOSED-FORM. Tame (regular-singular, no Stokes) = HOLONOMIC (D-finite, a linear ODE). Holonomic
splits by RIGIDITY (Katz). In z=(1+tanh u)/2 the tanh model is a rank-N Fuchsian system with 3 regular
singular points (z=0,1,inf; residues (H0-+A)/2i and -iA, regular-semisimple). Rigidity index
rig=(2-3)N^2+3N=N(3-N):
 - N=2: rig=2 => RIGID => Gauss 2F1 => Gamma-ratios, ELEMENTARY closed form (the (C) sech anchor, 1e-12).
 - N=3 (OURS): rig=0 => NON-rigid, 1 accessory parameter => HEUN class: holonomic but NOT elementary.

So for N=3, Demkov-Kunike moves us WILD (sigma; Stokes; nonlinear isomonodromy/Fredholm; beyond holonomic)
-> HOLONOMIC-but-non-rigid (Heun; a linear ODE + 1 accessory parameter). A real, large gain (linear ODE
replaces wild isomonodromy) but NOT closed form. Closed form needs BOTH a bounded ramp AND rigidity (N=2,
or a reducible/degenerate config).

TAXONOMY: elementary (rigid tame: 2F1, N=2 DK) C holonomic (tame: Heun, N=3 DK) -- both below wild (sigma,
N=3 LZ). REFRAMING: the N=3 LZ is the irregular/confluent limit of the N=3 tanh Heun system (z=0,1 confluence
to the rank-2 Weber point); under confluence the Heun ACCESSORY PARAMETER becomes the Stokes constant sigma.
=> sigma is the isomonodromy avatar of a Heun accessory parameter (consistent with rank-3/W3/confluent-Garnier).

Artifacts corrected: ws_geom_bundle_theory.md S8b (rigidity caveat + reframing); ws_tanh_demkov_kunike.py
(D + docstring: "Heun-class, holonomic, non-rigid", not "hypergeometric-class").
