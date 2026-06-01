# Open problem: the Type-1 N=3 middle-survival Stokes coefficient — statement & plan-of-attack

**Status:** for review. Everything below the "## Plan" line is a proposal to be critiqued before
dispatch. Companion: `paper/type1_lz_working_paper.tex` (the closed structural theory),
`paper/ws_e_junction_Smatrix.md` (the Laplace representation), `paper/ws_a_riemann_scheme.md`.

---

## 1. Formal statement

**Setup.** Type-1 N=3: `i ψ'(u) = (H₀+uA)ψ`, `A=diag(a₁,a₂,a₃)` (distinct), Cauchy couplings
`(H₀)_ij = γ_iγ_j(a_i−a_j)/(ε_i−ε_j)` (all nonzero), `(H₀)_ii = −Σ_{k≠i}γ_k²(a_i−a_k)/(ε_i−ε_k)`.
The transition matrix `P_{x→j}=|S_{xj}|²` is the modulus-squared of the connection (Stokes) matrix
`S` of this ODE across its single rank-2 irregular singular point at `u=∞`.

**Established reduction (this program).** `S` is determined by a scalar connection problem of
**confluent-Heun type, one accessory parameter above ₁F₂**. Concretely (WS-E, Laplace frame): with
`ψ_j(u)=∫_C e^{−iuv}B_j(v)dv`, `B'(v)=K(v)B`, `K(v)=−i·diag(1/a)·(H₀−vI)`, the eliminated scalar
3rd-order ODE `L[b]=0` has **one finite apparent (accessory) singular point** `v_*` with indicial
exponents `{0,1,3}` (gap at 2) plus the **rank-2 irregular point** at `v=∞`. Equivalently (WS-A,
`u`-frame) the same class appears directly. The two extreme survivals (BE) are the elementary
residue-collapsing periods; the **off-diagonal Stokes coefficient `S₁₂` (≡ middle survival `P₂→₂`)
is the open quantity**, carried by the single `12×13` spectral-network joint (WS-G).

**THE PROBLEM.** Compute `S₁₂` — equivalently `P₂→₂` and one independent off-diagonal — in closed
form as a function of `{γ,ε,a}`, with the natural arguments being the two window actions
`I_X = ∮_X √(Q₄) L_H W₄/p³ dλ` (the BE exponents) and the cross-ratio of the four complex turning
points; **or** prove that `S₁₂` is not reducible to classical special functions (a rigorous
irreducibility / transcendence statement), thereby pinning the named object exactly.

**Why it is hard (localized).** (i) No product/Bethe form: the Laplace-image matrices do not commute,
`[K(v₁),K(v₂)]≠0` (no `1/t` Coulomb term → no Fuchsian image). (ii) The accessory parameter of a
Heun-class equation is generically a transcendental function of the coefficients; a closed `S₁₂`
requires the Cauchy/Gaudin data to *fix the accessory parameter algebraically* (the make-or-break).
(iii) Confluent-Heun connection coefficients are not known in closed form in general.

**Key structural lever (the reason to be hopeful).** A rank-2 irregular point + regular/apparent
points on `P¹` is exactly the linear (isomonodromy) problem whose deformation is a **Painlevé
equation (expected: Painlevé V)**. The connection problem for the Painlevé V linear system *has* a
modern closed theory — isomonodromic τ-functions as Fourier series of `c=1` Virasoro conformal
blocks / Fredholm determinants (Gamayun–Iorgov–Lisovyy and successors). If the Type-1 `S₁₂` maps to
a Painlevé-V connection coefficient with the Gaudin data fixing the monodromy exponents and the
accessory parameter, then `S₁₂` has a *named* closed form (a τ-function ratio / a Nekrasov-type
sum), even though it is not elementary. **This is the central conjecture to test.**

---

## Plan (proposal — review before dispatch)

Ordered; each is a dispatchable workstream with a decisive gate. Numerics throughout are gated on
the `10⁻⁹` oracle (`experiments/oracle.py`).

**PA-0 — Verify the foundation (coordinator, fast).** Independently re-derive the scalar Laplace ODE
`L` and confirm the Riemann scheme: the apparent singularity `v_*`, its exponents `{0,1,3}` (no
log), and the rank-2 irregular point. *Gate:* the confluent-Heun-above-₁F₂ class is coordinator-
verified (currently rests on two agreeing agent frames). Blocks PA-1. *(This closes the paper's R8
TODO.)*

**PA-1 — Identify the Painlevé/isomonodromy type.** Put `L` (or the 3×3 system) in standard
isomonodromic form; determine which Painlevé equation governs its deformation (test Painlevé V vs
III/IV by the irregular Katz invariant / the singularity signature). Locate `v_*` relative to the
**universal real node** (is the accessory point the exact crossing?). *Gate:* a definite Painlevé
type + the monodromy/Stokes exponents expressed in `{γ,ε,a}` (expected `∝` the BE exponents `s_ij²|a_i−a_j|`).

**PA-2 — Does the Gaudin data fix the accessory parameter?** The make-or-break. Test whether the
Cauchy/Gaudin structure (rational spectral curve `(E,u)=(m/p,n/p)`; the Abelian ring; the BE +
double-stochasticity constraints) pins the Heun accessory parameter *algebraically* in `{γ,ε,a}`
(rather than transcendentally). *Method:* symbolic — compute the accessory parameter from `L` and
test for rational/algebraic dependence; check the node-as-apparent-singularity hypothesis. *Gate:*
algebraic ⇒ closed form is in reach (PA-3); transcendental ⇒ the answer is a genuine
Painlevé-V/confluent-Heun transcendent (PA-4 names it, PA-5 proves irreducibility).

**PA-3 — Closed-form connection coefficient (if PA-2 is algebraic).** Apply the Painlevé-V / CFT
τ-function connection formulae (or the confluent-Heun central-connection results) with the PA-1/PA-2
data to write `S₁₂` explicitly via the two window actions and the cross-ratio; recover BE and the
decoupling (incoherent) limit. *Gate:* matches the oracle to `≤10⁻⁶` across the ratio range
`0.1→4`, including the deep-overlap (sampleB, `~128×`) stratum.

**PA-4 — Validated numerical model + symbolic recognition (parallel safety net).** Build a
high-precision `S₁₂(I_{X1},I_{X2},\text{cross-ratio})` (mpmath/oracle); attempt inverse-symbolic /
PSLQ recognition against confluent-Heun and Painlevé-V connection constants. *Deliverable either
way:* a *computable* `P₂→₂` (a finite formula in named functions evaluated numerically), even if the
elementary-closed-form fails — this satisfies the practicality bar pragmatically.

**PA-5 — Literature + irreducibility.** Map the confluent-Heun / Painlevé-V connection-problem
literature (Lisovyy et al.; Jimbo's asymptotics; recent confluent-Heun connection results); determine
whether the relevant coefficient is published. If PA-2 is transcendental, assemble the rigorous
statement: `S₁₂` is the Painlevé-V/confluent-Heun connection coefficient with such-and-such
monodromy, provably not classical — the clean "named the object" result.

**Dispatch / coordination.** PA-0 then PA-1 first (foundation + type). PA-2 is the pivot. PA-3
(algebraic branch) and PA-4 (numerical/recognition, always-on) run after PA-1; PA-5 in parallel.
Same loop as before: research-log = shared state, `physics-reflection` gate before any promotion,
independent coordinator check of every load-bearing structural claim, tournament re-rank at PA-2.

**Honest odds (for review).** Elementary closed form: low (the structural results argue against it).
A *named* closed form (Painlevé-V τ / confluent-Heun connection coefficient with Gaudin-fixed data):
plausible, contingent on PA-2 being algebraic — this is the realistic "win." A clean irreducibility
theorem + a validated numerical model: the floor, and already a complete result.

---

## Addendum (user-approved, 2026-06-01): restart-operator structure + parallel tracks

**Decisions:** (1) Painlevé-V/τ route favored, **plus a parallel direct confluent-Heun
connection-coefficient track**. (2) PA-2 (Gaudin fixes the accessory parameter algebraically) is the
agreed make-or-break. (3) Acceptable deliverable = a **named closed form** (τ-function / connection
coefficient) **+ a computable numerical model**. (4) Fold in the restart-operator insights below;
fire H-R1/H-R2/H-R3 as parallel workstreams.

**Restart-operator structure (physical handle on `S₁₂`).** The single rank-2 irregular point at
`u=∞` is, in `λ`, split across the three poles `ε_i` (`λ→ε_i ⟺ u→∞`). Its local data factorizes
canonically as **(formal monodromy, diagonal) ⋉ (Stokes matrices, off-diagonal)**:
- the **diagonal/formal** part is already identified — `e^{2πi c_i}`, `c_i=Σ_{j≠i}s_{ij}²(a_i−a_j)`
  (the signed-BE Coulomb coefficient measured as the `log T` drift of `𝒮`); [established]
- the **off-diagonal** part is the **Stokes shears in the 2D carrier-space pairs `{mid,lo}` and
  `{mid,hi}`** — these *are* WS-G's `12×13` joint, and their non-commutative composition is `S₁₂≡P₂→₂`.
  The middle level is the shared subdominant partner; that is *why* the joint is `12×13` and why
  decoupling one outer link (WS-C trivial-coupling limit) collapses one shear → elementary. [conj]
- the form factor `Γ_j=(−u'(λ_j))^{-1/2}` is analytic at the poles (`∝(λ−ε_i)`); the WKB
  branch points are the turning points (`−u'=0`), not the `ε_i` — so the off-diagonal lives in the
  **Stokes sectors at `∞`**, not the `Γ_j` branch. [conj]

So the target `S₁₂` *is* the Stokes multiplier `σ_{mid,·}` of the irregular point — a concrete,
2D-reduced object. PA-1/PA-2/PA-3 inherit this: the Painlevé-V/confluent-Heun **monodromy data are
the formal exponents `c_i` (known) + these two Stokes shears**, and the connection coefficient is
their composition.

## Rank correction (WS-PV, 2026-06-01, coordinator-verified logic) — IMPORTANT
The "expected: Painlevé V" lever (and the "confluent-Heun" class label) are the **rank-2
decoupling-boundary** cases, **not** the generic answer. The Type-1 scalar reduction is genuinely
**rank-3** (3rd-order, irreducible over `Q(γ,ε,a)` — WS-C; WS-E's 3rd-order Laplace ODE). Painlevé V
is a 2×2 (rank-2) isomonodromy problem; confluent-Heun is a 2nd-order (rank-2) ODE — both are reached
**only on the decoupling locus** where the 3×3 Laplace system `K(v)=−i diag(1/a)(H₀−vI)` block-reduces
to `2×2 ⊕ 1×1` (verified: send one coupling `s_ij→0` ⇒ a level splits off). That locus is exactly the
*elementary* boundary (WS-C's "elementary ⟺ a level decouples"). **The generic `S₁₂` is therefore one
rank up: a rank-3 confluent-Garnier / `c=1` irregular-conformal-block connection constant** (Barnes-G
family), unpublished and not a classical special function. Published PV/PVI connection constants
(Lisovyy et al.) apply at the boundary only. [analytically-grounded; contingent on WS-PA1's PA-0
confirming the 3rd-order is irreducible / the `{0,1,3}` apparent point is not removable.]
**Unchanged:** PA-2 (does the Gaudin data fix the one accessory coordinate `σ` algebraically?) remains
the pivot — `σ` is exactly the Barnes-G argument; WS-NUM still guarantees a computable `P₂→₂`; the
deliverable is still a *named* (now rank-3) connection constant + the numerical model.

## Workstreams launched (parallel, 2026-06-01)
- **WS-PA1** (critical path): PA-0+PA-1+PA-2 — re-derive the scalar Laplace ODE, verify the `{0,1,3}`
  apparent singularity (closes R8 TODO), identify the Painlevé type, locate the accessory point vs
  the node, and **test whether Gaudin fixes the accessory parameter algebraically** (the pivot).
- **WS-PV** (PA-3+PA-5): Painlevé-V τ-function / `c=1` CFT connection route + literature/irreducibility.
- **WS-CH** (parallel direct track): the confluent-Heun **central connection problem** attacked
  directly (local solutions at the apparent singularity + irregular point; known CH connection results).
- **WS-NUM** (PA-4, always-on): high-precision numerical `S₁₂(`window actions, cross-ratio`)` +
  PSLQ/inverse-symbolic recognition → a computable `P₂→₂` regardless of closed form.
- **WS-R** (H-R1/H-R2/H-R3): extract the restart factor from `𝒮_canon`; test the
  (formal-monodromy diagonal `e^{2πi c_i}`) ⋉ (triangular Stokes shear in `{mid,·}`) structure, the
  mid-pair concentration + decoupling collapse, and the form-factor/branch-point loci.
