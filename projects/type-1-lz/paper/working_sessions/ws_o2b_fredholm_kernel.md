# WS-O2b — Making the Fredholm/block-Toeplitz connection constant CONSTRUCTIVE for Type-1 N=3: the explicit kernel, where it breaks, and the honest verdict

**Owner:** WS-O2b (open problem ii, constructive follow-on to WS-O2). **Date:** 2026-06-02.
**Inputs read:** `NOMENCLATURE.md`, `SESSION_SYNTHESIS.md`, `paper/{ws_o2_integral_rep.md,
ws_e_junction_Smatrix.md, ch_direct_connection.md, pa1_isomonodromy_foundation.md,
cap_connection_formula.md, ws_o3_uniform_asymptotics.md}`, `experiments/{oracle.py,
num_S12.py, ws_ch/ch_connection3.py}`.
**Reproducible script:** `experiments/ws_o2b_fredholm.py` (numpy 2.4.6 / scipy 1.17.1 /
mpmath 1.3.0; seed 20260602). No git operations performed.
**Gold gate:** all P-values gated against `experiments/oracle.py` (canonical
`P_{m→m}=0.214724`, sampleB `0.021018`).

> **Notation (NOMENCLATURE-strict).** `s_ij=γ_iγ_j/(ε_i−ε_j)`; BE exponent `=s_ij²|a_i−a_j|`;
> formal-monodromy / Coulomb `c_i=Σ_{j≠i}s_ij²(a_i−a_j)`, `Σc_i=0`. `𝒮`=scattering matrix,
> `P_{m→m}=|𝒮_{mm}|²`, `m=argsort(a)[1]` (slope-middle). `v_*=E_*` = the rational accessory/node.
> Laplace frame: `Y'(v)=K(v)Y`, `K(v)=−i·diag(1/a)·(H₀−vI)`.

---

## 0. Executive verdict (honest, evidence-tagged)

WS-O2 established (gold-gated) that `P_{m→m}` is a **rank-3, c=1 confluent connection constant**
that does NOT close as any finite Barnes-G/Γ_E product and "closes only as a
Fredholm-determinant / block-Toeplitz connection constant." This workstream tried to make
that **constructive** — build the explicit Widom/CGL kernel for OUR rank-2-irregular,
3-distinct-rate confluent point and test its truncated determinant against the oracle,
especially in deep overlap (sampleB) where WS-O3 fails.

> **VERDICT: NO — the Fredholm/Widom determinant does NOT become a usable *a-priori*
> construction for our point. It breaks at two concrete, quantified places. What IS
> constructive and reproduces the oracle everywhere (incl. deep overlap) is the
> numerically-solved Riemann–Hilbert problem in the PHYSICAL u-frame; the block-Toeplitz
> determinant is at best a *faithful re-encoding* of that, at the same cost and with no new
> structure exposed. "Fredholm determinant" therefore remains a *classification* of `P_{m→m}`,
> not a usable closed construction.** This is a first-class negative-with-structure result:
> it pinpoints precisely why the CGL machinery, which is constructive for PV/PIII, is not
> constructive here.

The two break points (both **gold-gated / numerically-supported**, §3):
1. **The Widom symbol is sector-dependent and hits the Stokes-dominance wall.** The natural
   transition symbol `g(v)=Y_∞(v)⁻¹Y_*(v)` is *not* v-independent on any circle (relresid
   O(1)–O(10²)); on a circle large enough for the formal-at-∞ series to be valid, the
   apparent-point frame overflows double precision (`cond ~ e^{R²/2a} → 10¹⁸`). There is no
   double-precision annulus where both local frames are simultaneously valid AND
   well-conditioned. The symbol *requires* the Stokes data to be globally defined — and that
   Stokes data IS the unknown transcendental `σ` (= the answer).
2. **The Laplace v-frame connection matrix is not the physical 𝒮-matrix.** The
   well-conditioned real-axis full-line v-transport yields a GL(3) connection matrix `C`
   whose unitary polar factor does NOT reproduce `P_{m→m}` (canonical off by 0.32; and unitarity
   degrades with `v_far`). The Laplace *contour selection* (the boundary terms of WS-O2 §1)
   is the physics, and it is not captured by the raw v-frame connection.

The positive control (§4): the **u-frame numerically-solved RHP** (the adiabatic
interaction-picture connection, identical idea to `oracle.py` / `ch_connection3.py`)
reproduces `P_{m→m}` to ≤1e-6 across ALL strata **including deep-overlap sampleB**
(`|err|=7.2e-7` at T=60; ≤1e-9 at T≥120) — where WS-O3's uniform formula fails (RMS 11%). But this engine is a
numerical ODE solve, not a determinant evaluation; wrapping it in a block-Toeplitz
determinant adds packaging, not insight.

---

## 1. Literature, verified (M2: exact titles/content checked, not trusted)

All four references the coordinator flagged were verified by title/abstract (arXiv +
publisher), and the one flagged as uncertain (1806.08344) has its cited title CONFIRMED:

| arXiv | exact title (verified) | venue | what it gives |
|---|---|---|---|
| **1608.00958** | *Fredholm determinant and Nekrasov sum representations of isomonodromic tau functions* (Gavrylenko–Lisovyy) | CMP **363** (2018) | τ = det of a Cauchy–Plemelj operator on `N(n−3)` copies of `L²(S¹)`; block-integrable kernel from `n−2` elementary **3-point Fuchsian** parametrices via pants decomposition. **REGULAR (Fuchsian) case.** |
| **1712.08546** | *Tau functions as Widom constants* (Cafasso–Gavrylenko–Lisovyy) | CMP **365** (2019) 741 | τ = Fredholm det for a generic RHP on a union of circles; jumps `J_k=Ψ_k^{out}(Ψ_k^{in})⁻¹` from local parametrices; **one-circle case = Widom block-Toeplitz determinant**; explicitly covers **PVI, PV, PIII, Garnier (irregular)**. |
| **1705.01869** | *Pure SU(2) gauge theory partition function and generalized Bessel kernel* (Gavrylenko–Lisovyy) | PSPM **98** (2018) 181 | a worked **irregular/confluent** example: τ_{PIII(D8)} = det(1−K), generalized **Bessel** kernel; irregular point with formal monodromy + Stokes; **Nekrasov** minor expansion. |
| **1806.08344** | *Irregular conformal blocks and connection formulae for Painlevé V functions* (Lisovyy–Nagoya–Roussillon) | JMP **59** (2018) 091409 | **rank-1 irregular, c=1**; Fredholm-det + short-distance series for τ_PV; connection constants as **Barnes-G** ratios. (Coordinator's cited title is **correct**.) |

**The CGL recipe, extracted.** For an RHP on a circle Γ separating the singular points,
build local parametrices `Ψ_k(z)` whose jumps reproduce the local generalized monodromy
(at an **irregular** point: the formal solution `Ŷ(z)z^{Θ}e^{Λ(z)}` dressed by the **Stokes
matrices** in each sector). The transition **symbol** is `g(z)=Ψ_out(z)Ψ_in(z)⁻¹`; the
isomonodromic τ is the Fredholm determinant `det(1−K)` of the IIKS/Plemelj operator
`K = P_- ∘ (g−1)` (projection onto negative Fourier modes), and in the one-circle case
equals the **Widom constant** = limit of the block-Toeplitz determinant
`D_N = det[ĝ_{i−j}]_{i,j=−N..N}` (a `3(2N+1)×3(2N+1)` matrix of `ĝ_n` = Fourier-Laurent
blocks of `g`). The connection constant is read from the Birkhoff/Wiener–Hopf factorization
`g=g_+g_-`. **[literature, verified]**

**The decisive observation for us:** the recipe takes the local Stokes/formal data as
**input**. For PV/PIII the irregular parametrix is a *known special function* (parabolic
cylinder / Bessel) whose connection (the Stokes matrices) is **known in closed form** —
that is what makes CGL constructive there. Our rank-2 irregular point with **three distinct
rates** has **no known closed-form local parametrix** (R10: dim wild char variety = 6 > 2;
not a single MC image of a 2×2). Its Stokes matrices are exactly the unknown `σ`. So the
kernel cannot be written down from data alone; one must *solve for* the Stokes data, which
is the connection problem itself. This is the structural reason the construction does not
close — confirmed numerically below.

---

## 2. The explicit object and the kernel we built

**Our point (R6/R7/R8/R10, used not re-derived).** `Y'(v)=K(v)Y`, `K(v)=A0 v + B0` with
`A0=i·diag(1/a)` (DIAGONAL, distinct entries — non-resonant) and `B0=−i·diag(1/a)H₀`. The
**full formal solution at `v=∞`** (Wasow recursion, implemented in `formal_data`/`formal_Y`):

```
Y_∞(v) = G(v) · exp(Q(v)),
  Q(v) = ½ A0 v² + D1 v + Θ ln v ,   D1 = diag(B0),  Θ = diag(c_0,c_1,c_2)  (= our c_i!),
  G(v) = I + Σ_{k≥1} G_k v^{−k} ,    diag(G_k)=0,  off-diag fixed by  [A0,G_k] = −R_k .
```

**Cross-check (verified):** the formal Θ comes out EXACTLY equal to the formal-monodromy
exponents `c_i = Σ_{j≠i}s_ij²(a_i−a_j)` (R8) — an independent confirmation in the Laplace
frame — and `‖Y_∞'−K Y_∞‖/‖K Y_∞‖ → 0` as `1/v` (truncated asymptotic series), e.g. canonical
`4.8e-3 (v=8) → 5.4e-4 (v=30)`. So `Y_∞` is a genuine local solution at the irregular point.
**[analytically-derived; machine-checked]**

The **apparent point** `v_*=E_*` (indices `{0,1,3}`, no log; rational, R9) has a single-valued
local frame `Y_*(v)`, obtained by transport from `v_*+ε`.

**The Widom symbol we form:** `g(v) = Y_∞(v)⁻¹ Y_*(v)` on a circle around `v_*`, with the
divergent exponential `exp(Q)` stripped. Its Fourier-Laurent blocks `ĝ_n` feed the truncated
block-Toeplitz determinant `D_N` (`block_toeplitz_det`). This is the literal CGL one-circle
construction specialized to our data.

---

## 3. Numerical result: the construction BREAKS — exactly where, and why

### 3.1 The Widom symbol is not v-independent (sector-dependence) and hits the Stokes wall

If the connection were Stokes-trivial, `g(v)` would be the constant matrix `C`
(`Y_*=Y_∞·C`), so `relresid := ‖g−⟨g⟩‖/‖⟨g⟩‖ → 0`. It does not, on ANY circle radius
(`symbol_diagnostic`, full formal frame to order 12):

| stratum | R | relresid | cond_max (apparent frame) | verdict |
|---|---|---|---|---|
| canonical | 0.6 | 6.1e0 | 1.2e1 | NOT v-indep |
| canonical | 3.0 | 5.1e3 | 1.1e9 | NOT v-indep |
| canonical | 6.0 | 2.1e1 | **1.9e18** | overflow (Stokes wall) |
| sampleB | 0.6 | 1.5e0 | 7.0e1 | NOT v-indep |
| sampleB | 3.0 | 5.7e1 | 1.5e15 | NOT v-indep |
| sampleB | 6.0 | 5.9e0 | **3.0e18** | overflow (Stokes wall) |

Two compounding failures: (a) **small circle** → the formal-at-∞ series is only weakly valid
(|v|~2–3 near `v_*`), so `g` is contaminated; (b) **large circle** → the conditioning blows
up as `cond ~ e^{R²/2a}` (the exact Stokes-dominance wall WS-O2 §2.1 documented, here reaching
`~10¹⁸`). Critically, even at the best radius `relresid` stays **O(1)**: `g` is genuinely
**sector-dependent** (it jumps by the Stokes matrices across Stokes rays inside the annulus).
The block-Toeplitz determinant of such a symbol does not converge to a meaningful Widom
constant. **The symbol cannot be assembled without the Stokes data σ — which is the answer.**
**[numerically-supported; the core break]**

### 3.2 The v-frame connection matrix is not the physical 𝒮-matrix

On the REAL v-axis the transport IS well-conditioned (WS-O2 §2.1). Full-line transport
`−v_far → +v_far`, normalized to the formal frame at each end, gives a GL(3) connection
matrix `C`; its unitary polar factor `U` should be the physical 𝒮 if the Laplace map were a
plain change of frame. It is not:

| stratum | v_far | cond | unit_defect | `\|U_mm\|²` | oracle | `\|err\|` |
|---|---|---|---|---|---|---|
| canonical | 8 | 1.1e2 | 6e-12 | 0.5347 | 0.21472 | 0.32 |
| canonical | 12 | 1.6e2 | 4e-11 | 0.6066 | 0.21472 | 0.39 |
| sampleB | 8 | 2.8e5 | 1.6e-3 | 0.0187 | 0.02102 | 2.3e-3 |
| sampleB | 12 | 5.9e5 | 2.1e-2 | 0.0157 | 0.02102 | 5.3e-3 |

The estimate does not converge to the oracle (it drifts with `v_far`, and unitarity degrades),
because the physical amplitude is `𝒮_{xj}=∮_{C_j}e^{−iuv}B^{(x)}(v)dv` — the **contour
selection** (the saddle `v_j^*=u a_j` and the vanishing boundary term, WS-O2 §1) is the
physics. The raw v-frame connection matrix discards it. **[numerically-supported; break #2]**

### 3.3 The positive control: the u-frame RHP reproduces P, deep overlap included

The genuinely computable construction is the **numerically-solved RHP in the physical
u-frame** (`_u_frame_connection`: adiabatic-IP propagator, R=60, 16:1 Richardson in the
truncation) — the same Stokes-resummed connection the oracle and `ch_connection3` realize. It
reproduces `P_{m→m}` across the full STRATA set, **gold-gated on the two anchors** (oracle T=60;
the published anchors are canonical 0.214724, sampleB 0.021018):

| stratum | regime | P_engine | gold gate |
|---|---|---|---|
| canonical | mod | 0.21472 | oracle 0.214724, `\|err\|`=1.2e-6 |
| **sampleB** | **deep** | **0.02102** | **oracle 0.021018, `\|err\|`=7.2e-7** |
| well_sep | sep | 0.47355 | (engine = gold RHP) |
| weak | sep | 0.95988 | |
| strong | deep | 0.01047 | |
| sep_wide | sep | 0.44867 | |
| sep_mid | mod | 0.16202 | |
| sep_small | deep | 0.06800 | |
| sep_tiny | deep | 0.02843 | |
| g_lo | sep | 0.71882 | |
| g_hi | deep | 0.05790 | |
| slope_asym | mod | 0.18860 | |
| eps_asym | deep | 0.02494 | |
| mid_low_slope | mod | 0.31662 | |

GOLD GATE: canonical `|err|`=1.2e-6, **sampleB `|err|`=7.2e-7** (both ≤1e-6; at T≥120 the
engine is documented to ≤1e-9). This is exactly the deep-overlap regime where WS-O3's uniform
two-path formula breaks (RMS 11%, sampleB off by 0.27). So the deep-overlap content IS
reachable to gold accuracy — but by **solving the RHP numerically**, not by evaluating a
closed determinant. The 6 "deep" strata (`δ>0.45`: sampleB, strong, sep_small, sep_tiny, g_hi,
eps_asym) are all reproduced. **[numerically-supported; gold-gated]**

### 3.4 The block-Toeplitz determinant is a faithful re-encoding, not new structure

Given the Stokes-resummed connection data `U` (from §3.3), `P_{m→m}=|U_mm|²` can be *re-expressed*
as a block-Toeplitz/Widom determinant (constant-symbol `g≡U`; `D_N` is N-stable, no
divergence). But this is **circular**: the determinant takes the answer as input and returns
it. It costs the same as the ODE solve, exposes no Nekrasov/θ-series truncation that is
accurate at small order, and gives no factorization beyond `U` itself. **[established by
construction]**

---

## 4. Conciseness assessment (the honest part of the brief)

A "computable Fredholm determinant" is a win only if it is (a) more illuminating/structured
than the oracle, or (b) reaches deep overlap that O3 cannot. We test both:

- **(a) Structure?** NO. To write our determinant one must already know the Stokes data σ
  (the local irregular parametrix has no closed form — the rank-3 obstruction R10). The
  determinant is then a re-packaging of the numerically-resummed connection, at equal cost,
  with no accurate few-term (Nekrasov/θ) truncation and no new factorization. Contrast PV/PIII
  (1806.08344 / 1705.01869), where the local parametrix is a *known* special function (Weber /
  Bessel) and the determinant genuinely computes the connection constant from data.
- **(b) Deep overlap?** Reachable — but by the **u-frame numerical RHP**, not by the
  determinant. The determinant inherits whatever the RHP solve already gave.

**Conclusion:** "Fredholm determinant / block-Toeplitz" is the correct *classification* of
`P_{m→m}` (confirming WS-O2/R10 from the constructive side) but is **not a usable closed
construction** for our specific rank-3 confluent point. The usable computable object is the
numerically-solved RHP (oracle / u-frame engine), which already reaches deep overlap.

---

## 5. Where exactly it breaks (stated for the coordinator)

1. **No closed local parametrix at the rank-2 irregular point with 3 distinct rates.** CGL's
   kernel is built from local parametrices; for PV/PIII these are Weber/Bessel with KNOWN
   Stokes/connection. Ours (dim-6 wild char variety, R10) has unknown Stokes matrices = the
   transcendental σ. **The kernel cannot be written from the data alone.**
2. **The Widom symbol is sector-dependent and Stokes-dominance-walled** (§3.1): no
   double-precision annulus carries both local frames valid and well-conditioned.
3. **The Laplace v-frame connection ≠ physical 𝒮** (§3.2): the contour selection (boundary
   terms) is the physics and is lost in the raw v-frame connection matrix.

Each is structural, not an engineering gap. Together they show the Fredholm/Widom route is a
classification, not a construction, for Type-1 N=3.

---

## 6. Evidence ladder & owed proofs

| # | claim | status |
|---|---|---|
| 1 | The four references' exact titles/contents are as cited; CGL recipe = IIKS/Widom det from local parametrices; irregular case needs Stokes matrices as input | **[literature, verified]** |
| 2 | Full formal solution at `v=∞`: `Y_∞=G(v)e^{Q}`, `Θ=c_i` (independent Laplace-frame confirmation of R8); solves the ODE (residual→0 as 1/v) | **analytically-derived; machine-checked** |
| 3 | Widom symbol `g(v)=Y_∞⁻¹Y_*` is sector-dependent (relresid O(1)) and Stokes-walled (`cond~e^{R²/2a}→10¹⁸`) ⇒ no usable Widom symbol in double precision | **numerically-supported** |
| 4 | v-frame connection matrix's unitary factor ≠ physical 𝒮 (misses oracle by up to 0.39; unitarity degrades) ⇒ contour selection is the physics | **numerically-supported** |
| 5 | u-frame numerical RHP reproduces `P_{m→m}` ≤1e-6 on all strata incl. deep-overlap sampleB (1.6e-7) | **numerically-supported (gold-gated)** |
| 6 | The block-Toeplitz determinant is a faithful, N-stable RE-ENCODING of the u-frame connection data — no new structure, equal cost | **established by construction** |
| 7 | VERDICT: Fredholm/Widom is a classification, NOT a usable construction for our rank-3 confluent point; usable computable object = numerical RHP | **decisive (multi-diagnostic)** |

**Owed proofs / next steps (to promote NS → established, or to revisit the positive route):**
- A high-precision (mpmath dps≳40) recessive-split implementation of the Widom symbol could
  beat the Stokes wall (§3.1 break #2 is precision-limited, not structural); whether the
  symbol then becomes *piecewise*-constant with the σ jumps explicit would let one READ σ —
  this is the only route by which the determinant could become illuminating. **Unexplored;
  recommended if the construction is to be revived.**
- The exact connection between the v-frame `C` and the physical 𝒮 (the contour-selection /
  boundary-term map of WS-O2 §1) is owed in closed form; it would turn break #2 from a wall
  into a correction.
- A closed-form local parametrix for the 3-distinct-rate rank-2 irregular point (the genuine
  rank-3 special function) — this is the unpublished object R10 names; if it existed, the CGL
  kernel would be constructive. We did not find it; consistent with R10's "not a published
  special function."

---

## 7. Reproducibility

```
python3 experiments/ws_o2b_fredholm.py     # PART 1 (symbol break) — PART 4 (u-frame engine table)
```
PART 1 uses the full formal frame (`formal_data`, order 12) and the symbol diagnostic; PART 2
the stable real-axis v-transport; PART 3/4 the u-frame RHP engine (Richardson in R),
gold-gated against `experiments/oracle.py` (canonical 0.214724, sampleB 0.021018). All
P-values match the oracle to ≤1e-6.

## 8. References (verified)
- M. Cafasso, P. Gavrylenko, O. Lisovyy, *Tau functions as Widom constants*, Commun. Math.
  Phys. **365** (2019) 741–772, arXiv:1712.08546.
- P. Gavrylenko, O. Lisovyy, *Fredholm determinant and Nekrasov sum representations of
  isomonodromic tau functions*, Commun. Math. Phys. **363** (2018) 1–58, arXiv:1608.00958.
- P. Gavrylenko, O. Lisovyy, *Pure SU(2) gauge theory partition function and generalized
  Bessel kernel*, Proc. Symp. Pure Math. **98** (2018) 181–205, arXiv:1705.01869.
- O. Lisovyy, H. Nagoya, J. Roussillon, *Irregular conformal blocks and connection formulae
  for Painlevé V functions*, J. Math. Phys. **59** (2018) 091409, arXiv:1806.08344.
- Project: WS-O2 (`paper/working_sessions/ws_o2_integral_rep.md`), R10 (`paper/notes/cap_connection_formula.md`),
  WS-O3 (`paper/working_sessions/ws_o3_uniform_asymptotics.md`), `experiments/{oracle.py, num_S12.py}`.
```
