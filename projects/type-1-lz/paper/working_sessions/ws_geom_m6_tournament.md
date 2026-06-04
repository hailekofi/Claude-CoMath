# WS-GEOM M6 — tournament: can the Feynman-graph machinery improve on Dykhne–Stückelberg?

**Type:** `physics-tournament` (rank & evolve) → feeding a `physics-derivation`/`physics-numerics` M6.
**Date:** 2026-06-03. **Question (user):** establish the Feynman-graph engine *not* for efficiency but to
see whether it can be **adapted to improve on the Dykhne–Stückelberg / uniform semiclassical law** (R12/R13).
**Baseline to beat:** `uniform_P_mm` (ws_o3_uniform.py): `P_mm ≈ A₀+A_ret+2√(A₀A_ret)cosΦ` with the
**fitted** proxy `cosΦ ≈ −(δ_lo+δ_hi)(1−K(1−χ))`, `K=3.0` calibrated on the gold set; ~1–2% accurate.

## Backdrop (measured, M5 Task A): truncation alone cannot do it
The bare Magnus/graph series is asymptotic (Λ≈π everywhere, R4/R5): order-2 is a wash vs order-1
(`e₂/e₁≈1.0` for sc≥2.2) and *worsens* at strong coupling. So "add graph orders" is refuted as an
improvement route. The question is whether *smarter use* of the same machinery (resummation /
uniformization / better saddle data) can improve the closed-form estimate.

## Candidates (Generation)
- **H1 — first-principles Stückelberg phase.** Replace the fitted `cosΦ` with the true phase
  `Φ = φ_dyn + φ_stokes`: `φ_dyn` = the action (period integral) between the two turning points on the
  explicit genus-0 curve Σ; `φ_stokes(δ_lo)+φ_stokes(δ_hi)` already in the code (computed, **unused**).
- **H2 — Zhu–Nakamura prefactors.** Subleading-saddle corrections to the single-crossing
  `p_x=e^{−2πδ_x}`. Sharpens the overlap regime; extremes are BE-exact already.
- **H3 — multi-path interference.** Beyond STAY/RETURN: higher graphs add multiple-return paths.
- **H4 — uniform Airy/Weber connection.** Uniformize across *coalescing* turning points (deep overlap,
  where DS fails). The genuine uniformization; closest to σ.
- **H5 — optimal-truncation / Borel–resurgence.** Treat the marginal Λ≈π series at optimal truncation;
  the non-perturbative ambiguity *is* σ (a resurgence reframing of σ).
- **H6 — graph residual.** Let the graphs compute `δP = P_mm − P_uniform` rather than `P_mm`.

## Proximity (clusters)
- *better phase:* **H1**  | *better amplitude:* H2 | *more paths/orders:* H3, H5 |
  *uniformization:* H4 (the coalescing-saddle generalization of H1) | *meta/residual:* H6.

## Ranking (Elo from pairwise debates on improvability × cost × correctness × non-opacity × σ-honesty)

| Rank | Hyp | Elo | One-line rationale |
|---|---|---|---|
| 1 | **H1** | **1322** | Highest improvability×(low cost)×transparency; targets the *one fitted constant*; the φ_stokes term is already in the code, unused; a single explicit period on a genus-0 curve (not opaque). |
| 2 | **H7=H1⊕H4** *(evolved)* | 1305 | Graph-derived **uniform** phase: isolated-saddle period (H1) continued by the Weber connection across coalescence (H4). Highest ceiling; **H1 is literally its first step**. |
| 3 | H4 | 1256 | High value exactly where DS degrades (deep overlap), but costlier and closest to the σ wall (full uniformization = σ, won't close). |
| 4 | H5 | 1232 | Deep structural insight (σ as a Borel/resurgence ambiguity) but needs many Magnus orders we don't have; not a cheap accuracy win. |
| 5 | H6 | 1212 | Cheap to test, synergistic with H1, but the marginal wall applies to δP too — uncertain payoff. |
| 6 | H2 | 1183 | Incremental; the exponents dominate and the extremes are already BE-exact. |
| 7 | H3 | 1141 | This is exactly the opaque matrix-product proliferation the program wanted to escape; marginal series. |

**Key debate outcomes.** H1 ≻ H3 (H1 fixes a *leading*-order defect — the phase — with one transparent
integral; H3 chases marginal higher-order paths and reintroduces opacity). H1 ≻ H2 (phase error, not
amplitude error, is what the fitted K is patching). H4 ≻ H3 (uniformization is the principled deep-overlap
fix; multi-path summation is its uncontrolled cousin). H1 vs H4: H1 wins on cost/decisiveness *now*; H4 is
the natural follow-on, hence the evolved **H7** that contains both.

## Evolution (new competitor)
**H7 (graph-derived uniform phase):** one formula `P_mm = A₀+A_ret+2√(A₀A_ret)cos Φ(Σ)`, with `Φ(Σ)` the
turning-point phase computed *uniformly* — the isolated-saddle period (H1) in the separated regime, smoothly
connected by the Weber/parabolic-cylinder form (H4) as the turning points coalesce — **no fitted constant**.
Its decisive first experiment **is H1**.

## The hard ceiling (carried into M6, non-negotiable honesty)
σ is irreducible (the open content is provably *not* a period). Best attainable: a **parameter-free** uniform
phase that **beats the fitted ~1–2%** in the moderate regime and degrades gracefully into the σ-core. H1
improves the *period part* of the phase; it cannot remove σ. If the construction turns opaque (the
exact-WKB trap the program already hit), that is the stop signal.

## Decision & decisive-test result
**Winner: H1** (with **H7 = H1⊕H4** as the target architecture). Ran the decisive first-cut
(`ws_geom_m6_h1_phase.py`): extract the *true* `cosΦ = (P_mm−A₀−A_ret)/2√(A₀A_ret)` from a reference solve
in the **interference regime** (moderate coupling, both A₀,A_ret appreciable; N=17), and compare, at *equal*
free-parameter count (one constant each), the fitted-K proxy vs a first-principles `cos(φ_dyn+φ_stokes)`.

**Result `[numerically-supported — a NEGATIVE]`.** *(Process note: an initial N=17 sample showed the
first-principles phase marginally **beating** the fitted K, 0.161 vs 0.189 — a tempting "structural win."
The robustness check overturned it.)* On a robust **N=42** interference sample:

| phase model | mean \|cosΦ − cosΦ_true\| |
|---|---|
| **fitted-K law** (1 fitted const) | **0.090** |
| parameter-free `φ_dyn + Σ₃ φ_stokes` (0 const) | 0.281 (~3× worse) |
| `φ_dyn + 2 φ_stokes`, best offset (1 const) | 0.179 (~2× worse) |

**The fitted K decisively beats the elementary first-principles phase** — even when the latter is given its
own fitted offset. The elementary turning-point phase (gap action + Stokes phases) does **not** reproduce
`cosΦ`; the fitted K is absorbing genuinely **non-elementary, σ-adjacent (non-period)** content. So **H1
fails**: the graph/turning-point machinery cannot improve Dykhne–Stückelberg with an *elementary* phase.

**Consequence for H7/H4 (escalation — do not build).** Beating K requires the **exact** turning-point
connection, which for this rank-3 (W₃, c=2) problem **is σ** — the Fredholm/Widom constant, i.e. the opaque
exact-WKB matrix product the programme already rejected. The H4 Weber/parabolic-cylinder uniformization is
exactly that exact connection; in deep overlap it is *more* σ-dominated, not less. **The escalation flag has
fired**: there is no non-opaque graph adaptation that beats the fitted constant. Per the hard constraint, we
stop rather than grind out a Weber computation of σ.

**The valuable takeaway (sharpens the σ story).** σ is not a small remainder bolted onto a good
semiclassical estimate: it is **load-bearing even in the interference *phase*** — large enough that a single
*fitted* constant outperforms a first-principles elementary phase. This is a sharp, first-class negative: it
tells you exactly why no elementary closed form for the interference exists, and that the fitted K in the
uniform law is the cheapest honest stand-in for σ-content, not a removable blemish.

## Answer to the original question
**Can the Feynman-graph machinery be adapted to improve on Dykhne–Stückelberg? — No, not with the
elementary (non-opaque) machinery.** Truncation is marginal (M5); the elementary first-principles phase is
*worse* than the fitted K (this M6); and the only route that could beat K is computing σ itself (the opaque
object). The graph engine remains valuable as the *microscopic derivation* of the uniform law (M5 Task B)
and as a transparent adiabatic-window calculator — but it does not sharpen the semiclassical estimate.
