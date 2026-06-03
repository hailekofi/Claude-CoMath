# WS-GEOM Milestone 4 — synthesis: the Geometric Selection Theory of the U(3) section

**Type:** `physics-theory-building` (synthesis of M1–M3). **Date:** 2026-06-03.
**Builds on:** M1 (`ws_geom_m1.md` — adiabatic-W Magnus skeleton), M2 (`ws_geom_m2.md` — image map /
no-hidden-invariant), M3 (`ws_geom_m3.md` — T1–T4, the rigid skeleton, now complete with the orientation
selector). **Spectral inputs:** R1 (Abelian ceiling / a-independent eigenbasis), R2 (genus-0 curve Σ),
R3 (universal real node, OWY 2009 arXiv:0807.0259), R8 (accessory = node, v_*=E_*), R15/R16 (the two
transcendentals {P_mm,b}). **Conventions:** NOMENCLATURE.md — ε-indexed states (`ε_0<ε_1<ε_2`); slope roles
`lo,mid,hi=argsort(a)`; `P[x,j]=prob(x→j)`; energy-rank is derived/end-reversing, never a channel label.

---

## 0. What M4 establishes (one paragraph)

The U(3) scattering section of a Type-1 N=3 MLZ model is a **dressed topological sector**. Its *shape* — the
sheet of U(3) it occupies — is a **rigid skeleton completely fixed by the geometry of the genus-0 spectral
curve Σ and its one real node** (R3): a node-selected directed 3-cycle (T1), reachable between exactly two
Birkhoff vertices (T2), bounded by the decoupling locus (T3), and carrying a Z₂ orientation that is the
**parity of the permutation relating the diabatic (ε) and adiabatic (slope) orderings** (T4). Its
*position on that sheet* — the two free numbers {P_mm,b} (R15) — is the **holonomy of a single
a-independent geometric seed W** against algebraic dynamical phases (M1), whose all-orders Magnus sum
contributes exactly **one irreducible Stokes constant σ** on the one recombining channel. **The geometry
selects the section; the dynamics only dresses it.** There is no further analytic invariant to find (M2),
and the only transcendental datum is σ. This is a *governing account* of the selection: every discrete and
abelian feature is derived; the entire residual transcendence is localized to σ.

---

## 1. The effective description (degrees of freedom and regime)

- **Object:** `S∈U(3)`, equivalently the doubly-stochastic `P_{xj}=|S_{xj}|²` (4 real DOF).
- **Frame:** the adiabatic frame `ψ=Φ(u)χ`, `iχ'=[D^{(a)}(u)−iW(u)]χ`, with `D=diag(E_i(u))` and the
  **geometric seed** `W_ij=⟨φ_i|∂_uφ_j⟩` (anti-Hermitian, off-diagonal). By R1, `Φ` and `W` are
  **a-independent** — the single connection whose holonomy is the whole commuting family.
- **Two layers (the central split):**
  1. a **topological / spectral skeleton** — discrete data read off Σ and its real node;
  2. an **analytic dressing** — the W-holonomy, i.e. the Magnus series, summing to σ.
- **Regime:** generic Type-1 N=3 (distinct `ε`, distinct slopes `a`, generic `γ`). The skeleton is
  exact everywhere; the *dressing's* leading skeleton equals the full S only in the adiabatic regime, and
  the dressing's all-orders content (σ) lives in the deep-overlap core (the M1/WS-O3 boundary).

## 2. Assumptions and symmetries (what the theory rests on)

| # | Assumption / symmetry | Used for |
|---|---|---|
| A1 | Type-1 Cauchy/Gaudin structure `(H_0)_{ij}=γ_iγ_j(a_i−a_j)/(ε_i−ε_j)` | the whole model; R4 width–spacing weld |
| A2 | Integrability / commuting ring (R1) ⟹ **a-independent eigenbasis** `Φ,W` | the single-seed principle; M1 frame |
| A3 | Genus-0 rational spectral curve Σ (R2) | one real node; explicit branch data |
| A4 | Universal real node (R3, OWY 2009): a unique real exact crossing `(u_*,E_*)` | T1–T4 skeleton; the node-swap |
| A5 | Accessory = node (R8): `v_*=E_*`, explicit algebraic | the spectral selector form |
| A6 | Brundobler–Elser law (extreme survivals) + double-stochasticity | R15 reduction to {P_mm,b}; τ_ext fixed |

Symmetry that organizes T4: the **two natural orderings** (ε-diabatic, slope-adiabatic) and the relabeling
group `S_3` acting between them; the orientation is the **sign character** of that relabeling.

## 3. The dependency map (each result, what it rests on, evidence-ladder)

```
A1 Cauchy ─┬─ R4/R5 width–spacing weld, permanent overlap ........ [AD]  (Λ≈π marginal everywhere)
           │
A2 ring ───┼─ R1 a-independent Φ,W  ── central principle (seed W) . [AD]
           │        │
A3 Σ g=0 ──┼────────┤
A4 node ───┼─ R3 unique real node ──┐
A5 R8 ─────┘                        │
                                    ▼
   M1: adiabatic-W Magnus skeleton; ORDER-0 = node-selected directed 3-cycle ......... [AD/NS]
        │  Ω_1→BE factors, Ω_2→interference; Λ≈π ⟹ order-improving, never resumming
        ▼
   M3-T1: node→directed-cycle   π_cont = τ_ext∘τ_node (3-cycle) ...................... [AD + verified]
   M3-T2: two-vertex reachable set {identity, node cycle}; others forbidden .......... [AD + NS]
   M3-T3: image boundary = decoupling locus (extreme→b=0, middle→P_mm=1) ............. [AD + NS]
   M3-T4: orientation = sgn(σ) [parity of ε→slope] = sign(tr H(u_*)−3E_*) ............ [AD + near-proof]
        │
   M2: image = {BE + double-stochasticity + χ}; NO hidden analytic invariant ......... [NS, first-class neg.]
   R15/R16: P affine in two transcendentals {P_mm,b} ................................. [AD + NS]
        ▼
   σ: the all-orders W-holonomy = ONE irreducible Fredholm/Widom Stokes constant ..... [OPEN — the sole datum]
```

Legend: `[AD]` analytically-derived, `[NS]` numerically-supported (near-proof), `[AD + near-proof]` both.

## 4. The unifying principle

> **Geometric Selection Principle.** For a Type-1 N=3 MLZ model the scattering operator factorizes as
> `S = (rigid topological sector) × (W-holonomy dressing)`. The rigid sector is **completely determined by
> the genus-0 spectral curve Σ and its unique real node**: it is the node-selected directed 3-cycle whose
> **orientation is the parity `sgn(σ)` of the diabatic↔adiabatic relabeling**, living between the two
> decoupling-bounded Birkhoff vertices. The dressing is the holonomy of the **single a-independent seed W**;
> its all-orders sum is **one irreducible Stokes constant σ** on the one recombining channel. *Geometry
> selects which section of U(3); dynamics only places the family within it.*

Three faces of the same statement: (i) *spectral* — the real node's sheet-pair fixes the discrete data;
(ii) *combinatorial* — `sgn(σ)` is the sign rep of `S_3`; (iii) *analytic* — everything not fixed by (i)/(ii)
is the marginally-convergent (Λ≈π) W-holonomy, irreducibly σ.

## 5. Examples and non-examples (the domain boundary)

- **Clean example (interior):** generic ε-monotone Type-1 (canonical `ε=[-2,0,3], γ=[1,.8,1.2], a=[-1,.5,2]`)
  — even `σ` ⟹ FWD cycle; node on the lower pair (`u_*=-0.249`); {P_mm,b}=(0.215,0.021) sits strictly
  inside the two-vertex region. The full theory applies.
- **Example (other sector):** any **odd** `σ` (e.g. `a` assigned `(0,2,1)` across ε) ⟹ REV cycle, node on the
  upper pair. Same theory, mirror sector.
- **Non-example / boundary (decoupling locus, T3):** drive a coupling to zero — the point leaves the
  interior and lands on a Birkhoff edge (`extreme→b=0`, `middle→P_mm=1`). This is the classifier boundary
  (R9): beyond it a level decouples and the model is elementary (rank-2). The theory *describes* this
  boundary; it is where the genuine-LZ region ends.
- **Non-example (orientation-flip walls):** the codimension-1 walls where `sgn(σ)` changes are exactly the
  **slope collisions** `a_i=a_j` (and ε collisions). There the real node runs to `u→±∞` and the slope
  labels relabel — the cycle is ill-defined on the wall. (Excluded by genericity.)
- **Outside the domain:** N>3, or Type-2/3 (non-Cauchy) networks — see §6 predictions.

## 6. Predictions (falsifiable, with status)

- **P1 [REFUTED — kept honestly].** *The orientation flips through a triple-degeneracy (`E_3=E_*` at the
  node).* **False.** A 1-parameter scan (`a_mid` through `a_lo`) shows the flip occurs with `|E_3−E_*|`
  *large* on both sides; the node instead runs `u_*→±∞` at a **slope collision**. The flip is a
  parity-change at an ordering wall, not a triple point. *(This refutation is what led to the correct
  combinatorial selector P2.)*
- **P2 [numerically-supported / near-proof].** *The orientation is `γ`-independent and equals `sgn(σ)`.*
  Confirmed: chamber-constant in `ε_1` and `γ_mid` scans; both sectors reproduced by the 6 slope-to-ε
  assignments; `γ`-independent over 80 draws; `sgn(σ)` matches on 119/119 + the spectral form on 154/154.
- **P3 [conjecture → testable].** *N>3 generalization.* The order-0 permutation is the **product of
  adjacent node-swaps over the real nodes of Σ**, dressed by W-holonomy; the orientation/parity data
  generalize to the sign of the diabatic↔adiabatic relabeling in `S_N`. Decisive test: N=4 Type-1, compare
  the continued order-0 permutation to the node-word and `sgn(σ)`.
- **P4 [conjecture → testable].** *Reachable-vertex count.* For general N the image touches exactly the
  **two** Birkhoff vertices `{identity, node-word permutation}` and no other; boundary = decoupling locus.
- **P5 [consistency, supported].** *σ is the unique transcendental.* Any independent representation (gauge,
  integral rep, Fredholm/Widom kernel, WS-O2b) must reproduce the **same** σ; no second analytic invariant
  appears (M2d). A mismatch would falsify the factorization.

## 7. Consistency checks

- **Limits.** Diabatic (`γ→0`): `S→1`, {P_mm,b}→(1,0) — identity vertex (T2). Adiabatic (`Λ` large):
  → node-selected cycle vertex (0,1) (T2, M1). Decoupling: edge (T3). All three consistent with the
  skeleton.
- **Unitarity / double-stochasticity.** Preserved at every Magnus order (M1), and BE+stochasticity is
  exactly the reduction to {P_mm,b} (R15). No over-counting: 4 DOF = 2 BE survivals + {P_mm,b}.
- **Dimensional/structural.** `Λ=∫‖W̃‖≈π` marginal *everywhere* (R4/R5) — the flow-side face of permanent
  overlap; guarantees the dressing never trivializes and never fast-resums, consistent with σ being
  irreducible.
- **Anchor.** Canonical numbers reproduced: `P_mm=0.214724`, sample-`b=0.021018`; node `u_*=-0.2493`,
  `E_*≈-1.9947`, lower pair, FWD.

## 8. Research overview / roadmap (meta-review)

The σ-free theory is **complete**. The single open frontier is **σ** itself — the irreducible
Fredholm/Widom connection constant of the rank-3 (W₃, c=N−1=2) confluent problem (WS-O2b). Ranked next steps:
1. **σ as a Fredholm determinant** (WS-O2b): nail the kernel; this is the only route to a closed form for
   {P_mm,b}, and the theory says it is the *whole* remaining content.
2. **N=4 test of P3/P4** (`physics-numerics`): cheapest decisive probe of the generalization; would promote
   the Geometric Selection Principle from N=3 to a structural law.
3. **σ-free corollaries** (already cashed): BE + the uniform semiclassical law (R12–R13) + two-component
   adiabatic-IP integration is the efficient recipe (Conclusion of the working paper).

**Bottom line for the physicist:** M4 closes the *geometry* question — *how these Hamiltonians select their
section of U(3)* — completely and σ-free. What remains is not geometric but analytic: the one Stokes
constant σ. The honest payoff of the geometric programme is a *map of why*, not a new closed form; the
closed form, if it exists, is σ's Fredholm determinant.
