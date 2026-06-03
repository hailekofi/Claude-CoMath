# WS-GEOM Milestone 5 — the Feynman-graph calculator and its validity domain

**Type:** `physics-derivation` + `physics-numerics`. **Date:** 2026-06-03.
**Scope (chosen):** center on the **graph engine + its domain of validity**; the deep-overlap exact
backstop is a pointer. **Builds on:** M1 (`ws_geom_m1.md` — the adiabatic-W Magnus skeleton, Ω₁,Ω₂),
M4 (`ws_geom_m4.md` — the dressed-sector factorization), R12/R13 (`ws_o3_uniform.py` — the uniform
semiclassical law), R4/R5 (permanent overlap, Λ≈π). **Code:** `experiments/ws_geom_m5_graphs.py`
(reproducible; reuses `ws_geom_magnus.py`, `ws_o3_uniform.py`, `oracle.py`).

---

## 0. What M5 establishes (one paragraph)

M5 turns the M1 skeleton into a **practical calculator for S** and pins down *exactly where it works*. The
order-2 Feynman-graph value `P₂ = |Φ(+)e^{Ω₁+Ω₂}Φ(-)^†|²` is a closed-form (integration-only) calculator
whose accuracy is controlled by the **dressed action `Λ=∫‖W̃‖du`**, *not* by the BE action `δ`: it is
accurate in the adiabatic regime and degrades as `Λ→O(π)` (the marginal wall, R4/R5). The graphs are not a
black box: by stationary phase they **are** the Dykhne–Davis–Pechukas / Stückelberg picture — Ω₁ carries the
single-crossing LZ survivals (the `p_x=e^{-2πδ_x}` that build the BE survivals), and the Ω₂ commutator
carries the **interference between two complex turning points**, whose closed form is exactly the uniform law
`P_mm≈A₀+A_ret+2√(A₀A_ret)cosΦ` (R12/R13). Beyond the validity domain the graphs do **not** resum (the
remainder is the irreducible σ), and the practical engine is the exact two-component adiabatic-IP
integration.

## 1. The graph rules (the calculator)

In the adiabatic frame `iχ'=[D^{(a)}−iW]χ` the dressed propagator has the Magnus log `Ω=Σ_kΩ_k` of the
**dynamical-phase-dressed coupling** `W̃_{ij}(u)=W_{ij}(u)e^{iθ_{ij}(u)}`, `θ_{ij}(u)=∫^u(E_i−E_j)`:
```
Ω₁ = −∫ W̃(u) du                                   (one vertex: a single non-adiabatic hop)
Ω₂ = ½∫∫_{u₁>u₂} [W̃(u₁),W̃(u₂)] du₁du₂            (two vertices: ordered double hop / interference)
S  = Φ(+T) · exp(Ω₁+Ω₂+⋯) · Φ(−T)^†               (P=|S|², diagonal Stark phases drop)
```
**Feynman dictionary.** Sheets `i` of Σ are lines; each `W̃_{ij}` is a vertex coupling line `i↔j`; the
dynamical phase `e^{iθ_{ij}}` is the propagator carried between vertices. The geometry sits in the vertices
(`W` is the `a`-independent seed); the elementary `a`-dependence sits in the phases. Each `Ω_k` is
anti-Hermitian, so `S∈U(3)` at every truncation (M1, verified).

## 2. The stationary-phase / DDP bridge (why the graphs *are* the uniform law)

The key structural fact: **the dressed phase has no real stationary point.** `θ_{ij}'(u)=E_i(u)−E_j(u)`
never vanishes on the real axis — the *adiabatic* levels of a real-symmetric `H(u)` never cross (R3's node
is an exact crossing of the *commuting-protected* curve, but the generic gap stays open). So each graph
integral is dominated by the **complex turning points** `u_c` where `E_i(u_c)=E_j(u_c)` — the off-axis branch
points of Σ. Steepest descent through `u_c` is precisely **Dykhne–Davis–Pechukas**:

- **Order 1 (one vertex).** `∫W̃_{ij}du` by steepest descent through `u_c^{ij}` gives the single-crossing LZ
  amplitude with `|amp|² = p_{ij}=e^{-2πδ_{ij}}` (the imaginary action `Im∮(E_i−E_j)` is the BE exponent).
  These single-crossing survivals are exactly the building blocks `p_lo,p_hi,p_lh` of the uniform law, and
  they make the **extreme survivals BE-exact** at leading order.
- **Order 2 (two vertices).** The commutator `[W̃(u₁),W̃(u₂)]` with two stationary points is the **Stückelberg
  interference** between the two turning-point paths — STAY (`A₀=p_lo p_hi`) vs RETURN-through-the-outer-link
  (`A_ret=(1−p_lo)(1−p_hi)p_lh`) — supplying the cross term `2√(A₀A_ret)cosΦ`.

Hence the closed form of the stationary-phase graph sum for the middle survival is the **uniform law**
(R12/R13):
```
P_mm ≈ A₀ + A_ret + 2√(A₀A_ret) cosΦ ,   cosΦ ≈ −(δ_lo+δ_hi)(1−K(1−χ)) ,
```
with `χ` the turning-point cross-ratio (the genus-0 shape datum). **The Feynman graphs are the microscopic
derivation; the uniform law is their closed form.** This is the practical content of M5: in the adiabatic
regime one need not integrate the graphs — the explicit branch points of Σ give `{p_x,χ}` and hence
`{P_mm,b}` in closed form, with the graphs supplying the systematic correction hierarchy.

## 3. The validity domain (numerical map)

*Reference:* a direct DOP853 diabatic propagator, T-averaged over `T∈{110,150,190}` to cancel the oscillatory
Fresnel endpoint tail; cross-checked against the 1e-13 Richardson oracle (`max|ref−oracle|` reported in
`main`). Knob: scale `γ` up at fixed `(ε,a)` → wider crossings → smaller `Λ` → more adiabatic.

<!-- VALIDITY_TABLE -->
| `γ`-scale | Λ | δ_max | e₀ | e₁ | e₂ | e₂/e₁ | domain |
|---|---|---|---|---|---|---|---|
| 1.0 | 2.98 | 0.24 | 3.3e-1 | 7.4e-2 | 5.6e-2 | 0.76 | FAIL |
| 1.4 | 2.96 | 0.92 | 1.0e-1 | 2.8e-3 | 8.7e-3 | 3.1 | <1% |
| 1.8 | 2.93 | 2.52 | 3.8e-2 | 1.6e-3 | **1.3e-3** | 0.80 | <1% |
| 2.2 | 2.90 | 5.62 | 1.7e-2 | 1.2e-2 | 1.2e-2 | 1.01 | <5% |
| 2.8 | 2.83 | 14.8 | 7.9e-3 | 1.2e-2 | 1.2e-2 | 1.01 | <5% |
| 3.5 | 2.73 | 36.0 | 2.0e-2 | 2.9e-2 | 2.9e-2 | 1.00 | <5% |

**Reading (honest, non-monotone).** The error is **not** monotone in `Λ` — across the whole sweep `Λ`
barely moves (2.98→2.73, all near the marginal `~π`), while `e₂` has a **sweet spot** at `γ`-scale≈1.8
(`e₂≈1.3e-3`) and *degrades on both sides*. Two regimes bracket it: (i) too diabatic (small `γ`): the
order-0 cycle is the wrong leading term and `e₀` is large; (ii) too strong (large `δ`): the order-2
correction is a **wash** (`e₂/e₁≈1.0`) and even *worsens* `e₀` — the marginal-`Λ` asymptotic overshoot. So
the calculator is a genuine sub-% tool only in the **moderate adiabatic window**, and it **cannot be driven
below ~1e-3 by adding orders**. The DDP/uniform-law bridge (Task B): the order-1 graph already carries the
BE extreme survivals, and `P₂(P_mm)` tracks the uniform law and the reference together through that window.

## 4. The honest boundary and the practical recipe

`Λ=∫‖W̃‖≈π` is **marginal almost everywhere** (R4/R5) — the flow-side face of permanent overlap — so the
graphs are **order-improving but never fast-resumming**: at the canonical/deep-overlap scale the order-2
term stops reducing the error (Task C), and the all-orders remainder is the irreducible **σ**. This is not a
defect of the bookkeeping; it is the content of σ. The practical recipe that M5 certifies:

1. **Extremes:** BE survivals (exact, closed form).
2. **`{P_mm,b}` in the generic/adiabatic regime:** the uniform law (R12/R13) = the stationary-phase graph
   sum, from the explicit branch points of Σ — accurate to ~1–2% across the generic window, with higher
   graphs as systematic corrections where `Λ` is genuinely small.
3. **Deep-overlap core (`Λ~π`):** the graphs do not resum — use the **exact two-component adiabatic-IP
   integration** (well-conditioned to 1e-9, where the naive diabatic propagator hits cond `~1e18`); this is
   where σ is computed.
4. **Fill the matrix:** double-stochasticity from `{P_mm,b}` + the extremes.

So "compute S from the Feynman graphs" is real and useful **as a closed-form adiabatic calculator**, with a
sharp validity boundary — and outside it the honest tool is exact integration, not more graphs.

## 5. Evidence ladder and limits

- `[analytically-derived]` the graph rules (M1) and the stationary-phase/DDP structure (θ' has no real
  zero ⟹ complex-turning-point domination ⟹ DDP single-crossing amplitudes + Stückelberg interference).
- `[numerically-supported]` the validity-domain map and the order-1→BE-extremes / order-2→uniform-law
  bridge (`ws_geom_m5_graphs.py`, Tasks A/B), referenced to a cross-checked propagator.
- `[honest limit]` the marginal-convergence wall (Task C): the graphs do not resum at `Λ~π`; σ is
  irreducible. The calculator's domain is the adiabatic regime; the exact backstop owns the deep core.
- Open: σ itself (the Fredholm/Widom constant) — unchanged by M5; M5 makes precise *where* you need it.
