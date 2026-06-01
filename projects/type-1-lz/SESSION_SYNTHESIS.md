# Session synthesis — Type-1 N=3 Landau–Zener (traceable research record)

**Purpose.** A self-contained, honest trace of the whole investigation (2026-06-01 session) for a
researcher picking this up cold. Chronological arc + hypothesis ledger (evidence-ladder status) +
**dead ends with WHY** + decisions/next steps. Companion files: `RESEARCH_LOG.md` (raw chronology),
`NOMENCLATURE.md` (symbol conventions — READ FIRST), `META_REVIEW.md` (recurring failure patterns),
`paper/` (deliverables), `experiments/` (reproducible code).

---

## 0. The problem and regime
Type-1 N=3 MLZ: `i dψ/du = (H₀+uA)ψ`, `A=diag(a)`, Cauchy coupling
`(H₀)_ij = γ_iγ_j(a_i−a_j)/(ε_i−ε_j)`, `(H₀)_ii=−Σ_{k≠i}γ_k²(a_i−a_k)/(ε_i−ε_k)`; `ε_0<ε_1<ε_2`.
Observable: doubly-stochastic `P_{x→j}=|𝒮_{xj}|²`, 4 free real params. Two **extreme-slope**
survivals are exact (Brundobler–Elser). **GOAL:** a *concise, computable* expression for the open
content — the **middle-slope** survival `P₂→₂=|𝒮_{mm}|²` (`m=argsort(a)[1]`) + one off-diagonal.
Nomenclature (locked): `s_ij=γ_iγ_j/(ε_i−ε_j)`, `w_ij=|2s_ij|`, BE exponent `=s_ij²|a_i−a_j|`;
`𝒮`=scattering matrix; `Γ_j`=Cauchy form factor (NOT the BE exponent); `S₁₂` DEPRECATED.

## 1. Honest bottom line (current state)
We have **not** produced a concise computable closed form, and we have **strong, multi-method
evidence that an exact one in elementary or *published* special functions does not exist.** What we
DO have: (i) a complete *structural theory* of why Type-1 N=3 is hard; (ii) the open quantity
*identified and characterized* as a **rank-3** connection coefficient with monodromy data fixed
*algebraically* by the exact crossing; (iii) a validated *numerical* model to 1e-7..1e-9. The
realistic remaining targets for "concise computable" are the **constructive exact-WKB product**
(two Weber shears × an explicit spectral-network junction), the **exact integral representation**,
or a **validated semi-analytic approximation** — see §6.

## 2. Established / supported results (evidence-ladder tagged)
| id | result | status | artifact |
|----|--------|--------|----------|
| R1 | Commuting Type-1 family is a **ring**; quadratic partner reduces to linear w/ linear-in-u coeffs; **shared eigenbasis** ⇒ integrability is **Abelian** and fixes only spectrum/BE, **never the prefactor** ("Abelian ceiling"). | **established** (machine-prec) | experiments/ring_structure.py |
| R2 | Spectral curve `Σ:det(E−H)=0` is **genus-0 rational** (Gaudin param `(E,u)=(m/p,n/p)`) with a structural **node** (real exact crossing). | **established** | paper/gate_test_genus.md |
| R3 | The real exact crossing (node) is **universal** (rational `u_*,E_*`). | **numerically-supported** (99.8% of 1290) | experiments/structural_crossing.py |
| R4 | **Width lemma** `w_ij=2|γ_iγ_j|/|ε_i−ε_j|` (slope-free, exact); crossings **permanently marginally overlap** ⇒ no exact `𝒮=∏𝒮_ij`; **never tridiagonal** (extreme–extreme coupling dominant). | **established** (symbolic) + NS (bound) | paper/coscaling_derivation.md, ws_d_verification.py |
| R5 | **Classifier:** elementary ⟺ a level decouples ⟺ joint-free spectral network ⟺ rank-2 boundary. | **numerically-supported** | paper/ws_c_factorization_locus.md |
| R6 | `P=|𝒮|²` = Stokes data of a **single rank-2 irregular point** at u=∞ (3-level Weber). | analytically-derived | paper/kz_isomonodromy_picture.md |
| R7 | Scalar Laplace ODE: one **apparent `{0,1,3}`** singularity + rank-2 irregular point. | **established** (3 independent derivations: WS-E, WS-CH, WS-PA1; symbolic, 2 samples each) | paper/{ws_e,ch_direct,pa1}*.md |
| R8 | **Formal monodromy = signed BE:** `c_i=Σ_j s_ij²(a_i−a_j)` (Σ=0); `|e^{2πc_i}|`=BE survival for extremes; `c_mid` is a **cancelling** sum ⇒ middle is the hard one. | **established** (verified) | ws_ch, ws_e |
| R9 | **Accessory point = node:** `v_*=E_*`=rational double root of discriminant; = unique rational root of cyclic Wronskian `det[e₀,Me₀,M²e₀]`. ⇒ **PA-2 = ALGEBRAIC** (accessory not a free modulus). | **numerically/symbolically-supported** (4 samples incl. non-monotonic) — NOT a general theorem | experiments/pa2_accessory_algebraic.py |
| R10 | **RANK-3 verdict:** the object is a rank-3 (confluent-Garnier / `c=1`) connection constant, **NOT** the published rank-2 Painlevé-V/confluent-Heun. dim(wild char variety) 6 vs PV's 2; gold-gated falsification `P₂→₂` outside the single-σ band (sampleB 0.021 vs [0.83,0.88]). PV/Heun apply **only at the decoupling (elementary) boundary**. | **numerically-supported (decisive tier T5 gold-gated)** + analytic (T1/T2, not machine-proved) | paper/cap_connection_formula.md |
| R11 | Computable `P₂→₂(γ,ε,a)` validated vs 1e-9 oracle across sep/width 0.1→4; args separate as `{BE window actions (scale), scale-invariant cross-ratio χ (shape)}`. | **gold / numerically-supported** | experiments/num_S12.py, oracle.py |

Anchor numbers: canonical `P₂→₂=0.214724`, sampleB `0.021018`; enhancement over incoherent 2.5–130×.

## 3. Hypothesis ledger (incl. demoted/refined)
- **Zero-curvature factorization `𝒮=∏𝒮_ij` (Malikis–Cheianov home run):** REFUTED for Type-1 — no MC `Ê` is non-Abelian (MC's is Abelian); the real obstruction is permanent overlap (R4). Status: **refuted**.
- **Painlevé-V (rank-2) via middle convolution (WS-PA1, PA-1):** REFUTED by capstone (R10): three generic distinct irregular rates can't be a 2×2 PV middle-convolution image; dim 6>2. Status: **refuted** (PV is the rank-2 boundary only).
- **Accessory parameter algebraic (PA-2):** **supported** (R9), the favorable branch.
- **`P₂→₂` = product of two Weber shears × explicit junction (constructive exact-WKB):** **UNEXPLORED** — the recommended redirect (§6).

## 4. DEAD ENDS (with WHY — the most valuable record)
1. **Commuting/quadratic partner → the prefactor.** WHY failed: the family is an Abelian ring (R1); a commuting partner carries only spectral data. *Could work only* for the BE/adiabatic part, never the non-Abelian prefactor. (Salvage: the framing led to R1, the Abelian ceiling.)
2. **Chernyak–Sinitsyn `ε₀=0` exact slice as an anchor.** WHY: it needs `g_OUTER=0`, but Type-1's outer–outer coupling is *never* zero (would need a zero coupling or equal slopes). The C–S solvable point is **outside** Type-1. (experiments/anchor_experiment.py)
3. **`Q_restart = D·P` cyclic relation (`(DPA)³=I` / `(P𝒮)³` diagonal).** WHY: `(DPDSD⁻¹)³`-type relations are achievable for *any* U(3) (verified: random U(3) hits it too) — a **content-free** general unitary fact; `D` is also non-canonical (cutoff/log-divergent). Salvage: extracted the geometric Coulomb phase `c_i` (= R8). (experiments/q_restart_probe.py)
4. **`(P𝒮)³` = formal-monodromy diagonal ⋉ *triangular* Stokes shear (H-R1 as stated).** WHY: the off-diagonal of `𝒮_canon` is **cyclic** (the 3-cycle sheet map), not triangular; and `e^{2πic_i}` is not a factor of the regularized `𝒮` (it's the subtracted log-drift — a category error). Refined to: cyclic perm × the two `{mid,·}` shears (= H-R2, which holds). (experiments/restart_probe.py)
5. **"Name it as a *published* special function" (Painlevé-V/Kampé de Fériet/confluent-Heun).** WHY: each was the wrong rank — Kampé de Fériet/Heun are genus-0/rank-2; the object is rank-3. The classification path proved *what it is not* but by construction never yields a concise formula. Strategic dead end (§5).
6. **Random-search corroboration of the χ-dependence (rank-3 signature).** WHY: matched-window/large-Δχ pairs are too rare to hit by chance (both 60- and 120-sample searches found none/only tiny-Δχ). Underpowered; the right tool is a **constrained deformation** (hold window actions fixed, sweep χ). Inconclusive, not corroborating.
7. **Coordinator/agent over-reaches caught by independent check (see META_REVIEW):** WS-D's "rigidity/null-space-0" obstruction (arithmetic wrong; real obstruction is R4); WS-F's "complex node" (node is real+universal, R3 — its own near_node sample had a real crossing at u≈−3 it missed); my "non-commuting Ê needed" framing (MC's Ê is Abelian).

## 5. Strategic arc & the diagnosed drift
Path taken: BE/window geometry → commuting-ring (Abelian ceiling) → genus-0 curve + universal node →
slope-free width lemma + classifier → `P` as rank-2-irregular Stokes data → scalar ODE `{0,1,3}` →
"name the special function" (KZ→Painlevé-V→confluent-Heun→`c=1`) → capstone **rank-3** verdict.
**Diagnosed error (2026-06-01 reassessment):** we pursued the *classificatory* (name-the-function)
path and reached "unpublished rank-3," **abandoning the project's original *constructive* exact-WKB
program** (the old "selector / global insertion / transport-recurrence / virtual-turning-points"
frontier). That original open piece — the "global selector" — is *exactly* the `12×13` joint we
characterized. We **characterized** it instead of **computing** it. We may also have **over-concluded
"no concise form"**: rank-3 means it's not a *single* rank-2 constant, but the *product* of two
explicit Weber shears is still a finite, concise, computable object — never tried.

## 6. Open questions & next steps
- **O1 [primary redirect]:** the **constructive exact-WKB product** —
  `P₂→₂=|M_form · 𝒮_Weber^{(mid,hi)} · 𝒮_Weber^{(mid,lo)}|²_{mm}` with each shear explicit
  (Γ-functions of `{s_ij,a_i}`) and the spectral-network **junction factor** computed. Candidate
  concise computable formula; the project's original aim. **Unexplored.**
- **O2:** the **exact integral representation** (WS-E skeleton) `P=|∮ e^{−iuv}B(v)dv|²` — a closed
  form (integrand solves the ODE).
- **O3:** a **validated semi-analytic approximation** (Dykhne leading + Stückelberg joint-phase).
- **Proofs owed (promote NS→established):** node universality (R3); node=accessory in general (R9 —
  cyclic-Wronskian path now in reach); the `≈1.6` overlap bound (R4); accessory algebraicity in
  general; T1's dimension count (R10).
- **Goal recalibration (pending user steer):** concede exact-elementary doesn't exist; aim at O1/O2/O3.

## 7. Reusable assets
Oracle `experiments/oracle.py` / `num_S12.py` (1e-9 `P`); the geometry/ODE builders in every
`experiments/*.py`; the four uploaded papers (BBGY 2409.17053, Chernyak–Sinitsyn 2006.15144,
Malikis–Cheianov 2505.06048, Lin–Sinitsyn 1310.7245) at `/root/.claude/uploads/...`; the assay
package `uploads/assay/`.
