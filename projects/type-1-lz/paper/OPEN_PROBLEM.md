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
type + the monodromy/Stokes exponents expressed in `{γ,ε,a}` (expected `∝ Γ_ij`).

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
