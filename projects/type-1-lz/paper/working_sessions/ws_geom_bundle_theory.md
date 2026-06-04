# The geometry of the commuting family: the transport cover, sections, and the space of amplitudes

**Type:** `physics-theory-building` (geometric/topological synthesis — NOT solving amplitudes).
**Date:** 2026-06-03. **Builds on / re-expresses:** R1 (Abelian ring), R2 (genus-0 Σ), R3 (real node),
R8/F1 (rank-2 irregular point at u=∞), R9/R11 (dim-6 wild character variety, W₃), M4 (central principle: one
a-independent seed W), R15/R16 (rank-2 amplitude skeleton), M2/T1–T4 (reachable set, cycle, orientation),
R17 (a-flow does not close), ws_ring_bundle (within-ring response is rank-2, mostly elementary). **Math
home:** wild Hitchin systems / wild nonabelian Hodge (Hitchin 1987; Beauville–Narasimhan–Ramanan 1989;
Jimbo–Miwa–Ueno 1981; Sabbah; Biquard–Boalch; Boalch).

---

## 0. Scope

We are **not** solving for transition amplitudes here. We build the geometric/topological frame in which the
Type-1 N=3 LZ problem lives, answer three structural questions — (1) what the sections of the vector bundle
over the transport cover are; (2) how *selecting a commuting family* is a *choice of section*; (3) how that
*shapes the space of transition amplitudes over the family* — and place every established result of the
project into its geometric slot. The payoff is a single picture: **the LZ problem is a wild Hitchin system on
ℙ¹; a commuting family is a section (the Hitchin section); the amplitude is the wild Riemann–Hilbert image of
that section; and the lone transcendental datum σ is the non-algebraicity of that map.**

## 1. The LZ ↔ geometry dictionary

The Schrödinger system `i ψ'(u) = (H₀+uA)ψ` is a **meromorphic connection** `∇ = d/du − i(H₀+uA)` on the
trivial rank-3 bundle over `ℙ¹_u`, with a single **irregular singularity at u=∞ of Poincaré rank 2**
(R8/F1). Its classifying data are:

| LZ object | Geometric object | status |
|---|---|---|
| sweep line `u` | base curve `ℙ¹` | `EST` |
| `H(u)=H₀+uA` | Higgs/Lax field; `∇` an irregular connection, rank-2 wild type at `∞` | `EST` (R8/F1) |
| eigenvalues `E_i(u)` | **spectral curve** `Σ ⊂ T*ℙ¹`, a 3-sheeted cover of `ℙ¹`, **genus 0** | `EST` (R2) |
| real exact crossing | a **real node** of `Σ` (nodal degeneration) | `EST` (R3, OWY) |
| eigenvectors `φ_i(u)` | **spectral (eigen-)line bundle** `L → Σ` (BNR datum); the **eigenbundle** `E=π_*L` | `EST` (constructed) |
| `W_ij=⟨φ_i|∂_uφ_j⟩` | the (a-independent) connection on `E` / Gauss–Manin transport | `EST` (M4) |
| commuting ring (fixed `γ,ε`) | the **abelian algebra** of the integrable system = the **Hitchin base** direction | `EST`→`framing` (R1) |
| a member `a` (slopes) | a point of the **Hitchin section**: a Higgs field of fixed irregular type ↔ the phase one-form `(E_i−E_j)du` on `Σ` | `framing` |
| dynamical phase `∫(E_i−E_j)` | periods of the phase (Seiberg–Witten-type) differential on `Σ` | `EST` (framework) |
| BE exponents `δ_ij` / Stückelberg phases | **imaginary** / **real** periods of that differential | `EST` (R12, framework) |
| transition matrix `S`, `P=\|S\|²` | **Stokes data** / wild monodromy = a point of the **wild character variety** (dim 6) | `EST` (R9/R11) |
| map `a ↦ S` | the **wild Riemann–Hilbert / nonabelian-Hodge map**, restricted to the Hitchin section | `framing` |
| `σ` (irreducible Stokes constant) | the **transcendence of RH** — `S` is not algebraic along the Hitchin base | `EST`(σ) + `framing` |
| "phase ambiguity" of `S` | the `U(1)³` Cartan torus = framing/structure group of `L` | `EST` |
| transport cover | `Σ` (where the eigenframe is single-valued) **+ its spinor double cover** (`δ_j=±1`, T4) | `AD`/`NS` |

`EST` = established (a proven/near-proof project result, re-expressed); `framing` = organizing identification
(the natural math home; consistent with all results but not itself a theorem here); `AD`/`NS` = analytically
derived / numerically supported.

## 2. Sections of the eigenbundle over the transport cover  [bullet 1]

**The transport cover** is the locus on which the adiabatic transport is single-valued. On `ℙ¹_u` the
eigenframe is **multivalued** — the eigenvectors permute around the turning points (the branch points of
`Σ`). They become single-valued on `Σ` itself, and their *phase/sign* becomes single-valued only on the
**spinor double cover** `Σ̃ → Σ` (the `δ_j=±1` sector, T4). So:

> **transport cover** `= Σ̃` (the genus-0 spectral curve, with its `ℤ₂` spinor double cover);
> **vector bundle** `= ` the eigenbundle `E` (rank 1 per sheet) with the a-independent connection `W`.

A **section** is a global choice of eigen-frame over `Σ̃`. Sections are constrained by three topological
data, all **intrinsic** to the commuting family (section-of-the-member-independent):
- the **formal monodromy** at the irregular point — *abelian* (F1): the diagonal exponential/Stokes
  structure that the frame picks up around `u=∞`;
- the **node gluing** (R3): a section must respect the two-sheet reconnection at the real node — this is the
  topological seam that later produces the directed cycle (T1) and the orientation `ℤ₂` (T4);
- the **spinor sign sector** `δ_j∈{±1}` — sections live on the double cover; the global sign class is a
  `ℤ₂`.
The **canonical section** is the adiabatic frame `Φ(u)`. The space of sections is
`Φ × (U(1)³` gauge `× ℤ₂` spinor `×` flat deformations of `W)`. Because the dressed action `Λ=∫‖W̃‖≈π`
everywhere (R4/R5), the bundle is **marginally non-flat**: its holonomy is `O(1)`, and that holonomy is the
entire transcendental content (the single seed `W`, M4). *This is the bundle-theoretic restatement of the
central principle.*

## 3. Selecting a commuting family is a choice of section  [bullet 2]

A *generic* LZ connection `∇` does not single out a frame. An **integrable structure** — a commuting family
— is precisely the extra datum that makes the multivalued eigenframe into a **coherent single-valued section
over the transport cover**, with an `a`-INDEPENDENT connection `W`. In Hitchin language:

> **A commuting family is a section of the Higgs-bundle fibration over the Hitchin base** — the *Hitchin
> section* — and the Type-1 ring is one such section. Selecting it fixes the spectral line bundle `L` (the
> BNR datum) coherently over the whole base; its defining property, the **a-independence of `W`** (M4), is
> exactly the statement that this section **trivializes the eigenbundle's connection over the base** (the
> seed is shared by all members).

Two layers of "choice" must be kept distinct:
- **choosing the family** = choosing the *section/polarization* (which coherent diagonalization) — Type-1 is
  one; other integrable types (Type-2, …) are other sections. *[bullet 2]*
- **choosing a member `a`** within the family = sliding along the *Hitchin base* (choosing the phase
  differential on `Σ`). *[feeds bullet 3]*
The Abelian ceiling (R1) is the statement that the section is **abelian** — the spectral/Hitchin data is a
torus (Jacobian/Prym of `Σ`), so the family fixes *spectral* data only, never the *Stokes* data. That is why
selecting the section pins the geometry but leaves the amplitude (a Stokes datum) un-fixed: **the section
chooses the base, not the fiber.**

## 4. The space of transition amplitudes over the family  [bullet 3]

Compose the Hitchin section (member `a ↦ ` spectral data) with the wild **Riemann–Hilbert map** (spectral
data `↦` Stokes data). The image is the **space of amplitudes over the family** — the "bundle over the ring"
of `ws_ring_bundle`:

- **Base** = the Hitchin base = elementary invariants `(δ_lo, δ_hi, χ)` (with `γ,ε` fixed, the meaningful
  base is 2-dim, `a` mod shift); the map `a ↦ (δ,χ)` is **elementary**. `[EST]`
- **Fiber** = the amplitude mod its `U(1)³` framing → after BE reduction, the two numbers `{P_mm,b}`
  (R15). `[EST]`
- **Reachable set** = exactly **two** Birkhoff vertices (identity, node-selected directed cycle), bounded by
  the **decoupling locus** (M2/T2/T3); the section's image is this two-vertex region. `[AD+NS]`
- **Discrete topology** = the node-selected directed **3-cycle** (T1) and its **`ℤ₂` orientation** =
  `sgn(σ_{ε↔slope})` (T4) — the latter the only *section-dependent* topological datum (it reads the member's
  ε↔slope ordering). `[AD+NS]`
- **Connection over the base** = **abelian-flat (formal monodromy) + a rank-2 σ-twist**: the tangent map is
  elementary except for `∂{P_mm,b}` (rank-2, `ws_ring_bundle`, residual 1e-6), and even that response is
  *mostly* captured by the elementary BE/uniform structure (cosine 0.92), with a subdominant transcendental
  residual. `[NS]`
- **Transcendence** = `σ` = the RH map's **non-algebraicity along the Hitchin base** (isomonodromy is
  transcendental); the `a`-flow does not close into a finite elementary system (R17). `[EST]`

So the space of amplitudes is a **rank-2 transcendental fibration over an elementary base**, topologically the
two-vertex region with `ℤ₂` orientation sectors, dressed by one irreducible isomonodromy constant `σ`.

## 5. Unifying principle

> **Geometric Selection, bundle form.** The LZ problem furnishes a single intrinsic datum — the eigenbundle
> `E` with its `a`-independent connection `W` on the transport cover `Σ̃`. A *commuting family* is the
> *section* (the Hitchin/polarization choice) that makes the eigenframe coherent and trivializes `W` over the
> base; a *member* is a point of the Hitchin base; the *amplitude* is the wild Riemann–Hilbert image of the
> section. **The section fixes the base (elementary spectral data); the holonomy of the one shared connection
> `W` supplies the fiber; and the whole non-elementary content is the single isomonodromy constant `σ`.**
> Geometry selects *which slice* of the wild character variety the family sweeps; it never computes the slice.

## 6. Intrinsic vs section-dependent (the topological ledger)

- **Intrinsic** (fixed by the family, independent of the member): genus-0 `Σ`; the real node `R3`; the rank-2
  irregular type at `∞`; the abelian formal monodromy (F1); the directed-cycle *cycle structure* (T1); the
  seed `W` (M4); the wild-character-variety ambient (dim 6, R9).
- **Section-of-the-member-dependent**: the cycle **orientation** `ℤ₂` (T4, `=sgn` of the ε↔slope
  permutation); the position `{P_mm,b}` within the two-vertex region; the dynamical periods `(δ,χ)`.
- **Neither (the irreducible)**: `σ` — a property of the RH map itself, not of base or fiber separately.

## 7. Dependency map (geometric slot ← established result)

```
ℙ¹ + irregular connection (rank-2 at ∞) .......... R8/F1            [EST]
spectral curve Σ, genus 0 ........................ R2               [EST]
real node of Σ ................................... R3 (OWY)         [EST]
eigenbundle E + a-independent connection W ....... M4 / R1          [EST]
commuting family = Hitchin section ............... R1 (+ framing)   [EST→framing]
member a = Hitchin-base point .................... ring (ws_ring)   [EST]
amplitude = wild RH / Stokes data ................ R9/R11           [EST]
  target: wild character variety, dim 6 (W₃) ..... R9/R11           [EST]
RH image over the family:
  base = (δ_lo,δ_hi,χ), elementary ............... R11/R16          [EST]
  fiber = {P_mm,b}, rank-2 transcendental ........ R15/R16          [EST]
  reachable = two vertices, decoupling boundary .. M2/T2/T3         [AD+NS]
  orientation ℤ₂ = sgn(ε↔slope) .................. T4               [AD+NS]
  connection = abelian-flat + rank-2 σ-twist ..... ws_ring_bundle   [NS]
  transcendence σ = RH non-algebraicity .......... R11b/R17         [EST]
```

## 8. Predictions and domain

- **P1 (other types = other sections).** Type-2/3 integrable MLZ families are *other Hitchin sections* over
  the same wild character variety; their amplitude spaces are other slices. Testable: build a Type-2 N=3
  family and check it is the *same* wild-character-variety ambient with a *different* reachable slice.
  `[conjecture]`
- **P2 (intrinsic topology is member-independent).** The cycle structure, formal monodromy, and node
  reconnection do **not** change as `a` varies within the family (only the orientation `ℤ₂` can flip, and
  only at ordering walls). Already supported (T1/T4 chamber-constancy). `[NS→testable for N=4]`
- **P3 (dimension count).** The wild character variety has dim 6 (R9); after BE + double-stochasticity +
  phase framing, the *observable* fiber is 2-dim `{P_mm,b}`. The reduction `6 → 2` is the bundle's "visible"
  rank; predicts the analogous N-count (wild-cv dim vs observable fiber dim) for N>3. `[conjecture]`
- **Domain / non-examples.** The frame is *clean* on the genuine-LZ interior (permanent overlap, `Λ≈π`); it
  **degenerates** on the decoupling locus (a level decouples → `Σ` reducible → rank drops → elementary,
  R9 classifier) — that boundary is exactly where a section/polarization becomes singular. Breaking
  integrability (leaving the ring) destroys the section entirely (opens the node; the perturbation analysis,
  ws_pert_crossing) — that is *off* this geometry.

## 9. Research overview / roadmap (meta-review)

The geometry is now organized; the productive next moves, ranked by leverage:
1. **Pin the wild-RH dictionary** (make the `framing` rows `EST`): identify explicitly the irregular type,
   the Stokes structure, and confirm the dim-6 wild character variety is the Boalch wild-cv for this type
   (`physics-literature` + `physics-derivation`). Highest conceptual payoff; mostly bookkeeping against
   Boalch/Sabbah.
2. **Characterize the section as a point in `Jac(Σ)`/Prym** (the BNR line bundle): does the Type-1 family
   correspond to a distinguished (e.g., theta-characteristic / canonical) line bundle? Would explain *why*
   `W` is a-independent. `[conjecture]`
3. **The two-vertex reachable set as a wild-cv stratum**: is `{identity, cycle}` a natural cell of the wild
   character variety (Boalch's cluster/quiver structure)? Connect M2/T2 to the cluster stratification.
4. **N=4 test of P1–P3** (`physics-numerics`): the cheapest decisive probe that the picture is structural,
   not N=3-special.

## 10. Honest status

This note is **mostly a reframing**: every quantitative claim is an established project result
(R1–R20, M-series, ws_ring_bundle) re-expressed in the wild-Hitchin/RH language. The **new content is the
organizing identification** — commuting family ↔ Hitchin section, amplitude ↔ wild RH image, `σ` ↔ RH
transcendence — tagged `framing` throughout. These identifications are the natural and (to the evidence)
consistent mathematical home, but they are **not theorems proved here**; promoting them (§9.1) is the work.
No transition amplitude is solved; the value of the note is the map of *where the geometry lives and what is
intrinsic vs chosen vs irreducible*. The lone irreducible `σ` is, in this language, the statement that the
**wild Riemann–Hilbert map is transcendental on the Hitchin base** — which is exactly what isomonodromy
theory expects, and what R11b/R17 establish directly.
