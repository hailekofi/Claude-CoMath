# The node-pinned GL₃ Riemann–Hilbert problem for Type-1 N=3 MLZ (rung-1 specification)

**Workstream WS-RH.** **Date:** 2026-06-02. Companion evaluator:
`experiments/gl3_rh_solver.py` (mpmath, reproducible). **Conventions:** `NOMENCLATURE.md`.
**Inputs built on (not re-derived):** `paper/{ws_o2_integral_rep.md, ws_o2b_fredholm_kernel.md,
r8_accessory_node_proof.md}`, R6/R7/R8/R9/R10, `experiments/{oracle.py, num_S12.py}`.

> **Notation (NOMENCLATURE-strict).** `s_ij=γ_iγ_j/(ε_i−ε_j)`; BE exponent `=s_ij²|a_i−a_j|`;
> formal-monodromy / Coulomb `c_i=Σ_{j≠i}s_ij²(a_i−a_j)`, `Σc_i=0`. `𝒮`=scattering matrix,
> `P_{m→m}=|𝒮_{mm}|²`, `m=argsort(a)[1]` (slope-middle). `v_*=E_*`=rational accessory/node.
> Laplace frame: `B'(v)=−i M(v)B`, `M(v)=diag(1/a)(H₀−vI)`; equivalently `Y'(v)=K(v)Y`,
> `K(v)=A0 v+B0`, `A0=i·diag(1/a)`, `B0=−i·diag(1/a)H₀`.

---

## 0. Scope and one-paragraph summary

This is **rung 1** of the node-pinned GL₃ oper / RH program: an explicit, self-contained
specification of a Riemann–Hilbert problem whose connection data carries the Type-1 N=3 MLZ
scattering matrix `𝒮`, with **all singularity data algebraic in (γ,ε,a)**. Everything below
is assembly of established project results (R6–R10) plus one new analytic step: the
**steepest-descent G± normalization maps** (§5) that relate the physical incoming/outgoing
bases to the canonical Laplace-sector bases. The matching numerical evaluator is
`experiments/gl3_rh_solver.py` (rung 2). **Rung 3** — a closed `W₃`/Nekrasov *series* for the
connection constant `σ` — is FRONTIER/out of scope and **not attempted** here.

The honest headline (developed in §6 and verified by the evaluator): the RH **data** is
algebraic and explicit; the G± maps are **pinned analytically up to a finite, sample-
independent convention/calibration**; but the **central connection constant `σ` itself**
(the off-diagonal Stokes data of the 3-distinct-rate irregular point) has **no known closed
form** (R10) and must be obtained by solving the RH problem numerically. The evaluator does
this in the physical u-frame to gold accuracy including deep overlap, and in the v-plane
(the genuine oper frame) it beats the double-precision conditioning wall but stalls precisely
at the connection→physical (G±) map — a first-class negative result locating the open piece.

---

## 1. The oper / linear system (the RH "carrier")

The Laplace–Euler dual of `i ψ'(u)=(H₀+uA)ψ(u)`, `A=diag(a)`, via `ψ(u)=∮_C e^{−iuv}B(v)dv`
(WS-O2 §1, machine-certified to `2×10⁻¹⁴`):

> **`B'(v) = K(v) B(v)`, `K(v) = −i·diag(1/a)·(H₀ − vI) = A0 v + B0`,**
> `A0 = i·diag(1/a)` (diagonal, **distinct** entries → non-resonant),
> `B0 = −i·diag(1/a)·H₀`.

This is a rank-3 (`3×3`) first-order linear ODE — a **GL₃ oper** on `ℙ¹_v`. Its singular
locus is **two points**: an irregular point at `v=∞` and one finite **apparent** point at
`v_*=E_*`. The RH problem is to reconstruct the global solution from local data at these two
points plus jumps; the physical observable is read from its `v=∞` connection structure.

**Gauge simplifier (established, WS-O2 §1.2).** `H₀` depends on `a` only through differences
`a_i−a_j`, so the shift `a→a+β` leaves `H₀` invariant and multiplies `ψ` by `e^{−iβu²/2}`,
which drops from `P=|𝒮|²`. **WLOG all `a_j>0`**, collapsing the three steepest-descent
directions into one Stokes wedge (used throughout §4–§5). *(Caveat from R8: `β` may not send
any `a_j→0`, since `diag(1/a)` must stay finite.)*

---

## 2. Singularity data at `v=∞` (irregular point) — algebraic

**Poincaré rank 2** (the `A0 v` term). The full formal (Wasow) fundamental solution:

> **`Y_∞(v) = G(v)·exp(Q(v))`,**
> `Q(v) = ½ A0 v² + D1 v + Θ ln v`, `D1 = diag(B0)`, **`Θ = diag(c_i)`**,
> `G(v) = I + Σ_{k≥1} G_k v^{−k}`, `diag(G_k)=0`, off-diagonals fixed by `[A0,G_k]=−R_k`.

- **Three rates** (leading exponents `½ A0_j v²`): `A0_j = i/a_j`, **distinct** (Type-1). These
  three distinct quadratic rates are the rank-3 signature (R10): not a single rank-2 (PV/Weber)
  point.
- **Formal exponents** `Θ = diag(c_i)` with `c_i = Σ_{j≠i} s_ij²(a_i−a_j)`. This is an
  **independent Laplace-frame re-derivation** of R8 (the formal-monodromy = signed-BE
  coefficients), and the evaluator confirms `Θ(Wasow) == c_i` to machine precision (PART 0).
- **Linear drift** `D1_j = (B0)_jj = −i(H₀)_jj/a_j`.
- **Formal monodromy** `M_form = exp(2πi Θ) = diag(e^{2πi c_i})`: `|e^{2πi c_i}|` reproduces the
  two exact Brundobler–Elser extreme survivals; `c_mid` is a *cancelling* sum (R8) — which is
  exactly why the middle survival is the hard, non-Abelian leftover.

All of `A0, B0, D1, Θ, G_k` are **algebraic (rational) in `(γ²,ε,a)`**. The Wasow series is
asymptotic (zero radius of convergence); its truncated residual `‖Y_∞'−K Y_∞‖/‖K Y_∞‖ → 0`
as `1/v` (evaluator PART 0; canonical `≈4.8×10⁻³` at `v=8`, order 8).

---

## 3. Singularity data at `v_*=E_*` (apparent point) — algebraic

- **Location** `v_* = E_* = −W₀/W₁`, the unique finite root of the cyclic Wronskian
  `W(v)=det[e₀ | M(v)e₀ | M(v)²e₀]` (linear in `v`); equivalently the OWY node energy.
  **Rational** in `(γ²,ε,a)` and **PROVEN** `v_*≡E_*` (R8, polynomial identity).
  Canonical: `v_* = −748/375`. sampleB: `v_* ≈ −2.5488`.
- **Local exponents** (Riemann indices) `{0,1,3}` (a gap at 2), **apparent / no log**: the
  local solutions are single-valued and meromorphic (R8 §6; `{0,1,3}` proven on the gauge
  slice, no-log certificate slice-symbolic + numeric). The local frame `Y_*(v)` is therefore a
  clean holomorphic GL₃ frame obtained by transport from `v_*+ε`.
- **Consequence:** the accessory parameter is **not a free transcendental modulus** — it is
  fixed algebraically (R8/R9 corollary). The RH problem has *rigid* accessory data; the only
  transcendental content is the connection constant `σ` at `v=∞`.

---

## 4. Stokes structure at `v=∞` (the Stokes graph in the v-plane)

With the `a_j>0` gauge, all three leading exponents are `½(i/a_j)v²` with `a_j>0`. Writing
`v=re^{iθ}`: `Re[½(i/a_j)v²] = −(r²/2a_j) sin 2θ`. Hence (identical angular structure for all
three channels, *different rates* `1/a_j`):

- **Anti-Stokes rays** (`Re=0`, oscillatory): `θ = 0, 90°, 180°, 270°` (the real and imaginary
  axes). On these the formal exponentials are purely oscillatory.
- **Stokes rays** (steepest, maximal dominance): `θ = 45°, 135°, 225°, 315°`. There are
  `2·(Poincaré rank) = 4` Stokes rays → **4 Stokes sectors** of opening `90°`.
- **Dominance ordering.** Within a sector the three solutions are ordered by `1/a_j`; the
  ordering *swaps pairwise* across each Stokes ray. With three **distinct** rates this is a
  genuine 3-level Stokes graph (not reducible to independent 2-level Weber pieces): the source
  of the rank-3 `σ`.
- **Saddles for the Laplace integral.** The phase `Φ_j(v,u)=−iuv+½(i/a_j)v²` of channel `j`
  has saddle `v_j^*(u)=u a_j` (curvature `Φ_j''=i/a_j`). For `u→−∞` the three saddles sit at
  `v=−|u|a_j<0` (incoming); for `u→+∞` at `v=+|u|a_j>0` (outgoing). The descent contour
  through `v_j^*` runs at `45°` to the real axis (the Stokes-ray direction).

**The Stokes graph (v-plane), schematic.**
```
              Im v
               │   ·  Stokes ray 90°(anti) ·
        sector │135°(Stokes)        45°(Stokes)
          II   │      ╲           ╱   sector I
               │        ╲       ╱
   ────────────┼──────────╳────────────── Re v   (anti-Stokes; real axis)
        v_*●   │        ╱   ╲    saddles v_j^*=u a_j on Re v
          III  │      ╱       ╲    sector IV
       225°(Stokes)              315°(Stokes)
               │
```
`v_*=E_*` sits **on the real (anti-Stokes) axis**, left of the origin (`v_*<0` for the
anchors). The real axis is anti-Stokes → the naive real-axis connection is oscillatory and
mixes the apparent point; the *physical* contour is the `45°` descent path. **This is the
geometric content of WS-O2b's "wall #2": the contour/sector selection is the physics.**

---

## 5. The connection / RH data and the G± normalization maps

### 5.1 The RH jumps

The Riemann–Hilbert data is: (i) the local formal frame `Y_∞` (§2) dressed by **Stokes
matrices** `S_1,S_2,S_3,S_4` across the four Stokes rays; (ii) the formal monodromy
`M_form=e^{2πiΘ}`; (iii) the single-valued apparent frame `Y_*` (§3). The **central connection
matrix** `C` relates them, `Y_*(v) = Y_∞(v)·C`, and the full monodromy at `v=∞` is
`M_∞ = M_form · S_4 S_3 S_2 S_1`. The transcendental content is the off-diagonal Stokes data
collectively denoted `σ` (NOMENCLATURE): `C = (formal monodromy) ⋉ (shears in σ)`.

### 5.2 The G± maps (steepest-descent, the new analytic step)

The physical amplitude is `𝒮_{xj}=lim_{u→+∞} e^{+i a_j u²/2}[∮_{C_j} e^{−iuv}B^{(x)}(v)dv]_j`
(WS-O2 §1.1). Evaluate the contour `C_j` by steepest descent through `v_j^*(u)=u a_j`. The
saddle value of the phase is

```
Φ_j(v_j^*) = −iu(u a_j) + ½(i/a_j)(u a_j)² + D1_j(u a_j) + c_j ln(u a_j)
           = −½ i a_j u²  − i(H₀)_jj u  + c_j ln(u a_j),
```

using `D1_j a_j = −i(H₀)_jj`. **The first two terms are EXACTLY the diabatic Stark asymptote**
`−½ i a_j u² − i(H₀)_jj u` of channel `j`; the third is the `c_j`-log drift (the same drift
subtracted in `num_S12.S12_canonical_phase`). The Gaussian fluctuation integral gives the
prefactor `√(2π/(−Φ_j''))=√(2π a_j/(−i))`. Therefore:

> **`G± = diag( g_j^± )`, a DIAGONAL, channel-wise map**, with
> `g_j^± = √(2π a_j/(−i)) · (|u| a_j)^{c_j} · e^{±iπ·(phase convention)}`,
> and the physical scattering matrix is **`𝒮 = G_+^{-1} · C · G_-`**, `C` the central
> connection matrix of §5.1.

**What is pinned analytically vs. convention:**
- *Pinned:* G± is **diagonal** (one factor per channel), and the saddle phase = the diabatic
  Stark phase (above) — so G± carries **no off-diagonal mixing**. The Gaussian modulus
  `|g_j^±|=√(2π a_j)` and the `(|u|a_j)^{c_j}` drift are explicit and algebraic.
- *Convention/calibration:* the **branch of `√(−i)`** (the `±π/4` Stokes phase per channel) and
  the assignment of saddle `j` to Stokes *sector* (incoming `−∞` vs outgoing `+∞` sheet) is a
  finite, **sample-independent** discrete convention. Per the brief, this fixed convention may
  be **calibrated on ≤2 anchors and validated on the rest** (it is a fixed map, not a per-
  sample fit).

### 5.3 The decisive structural consequence (and the stall)

Because `G±` is **diagonal**, the modulus-squared factorizes:
`|𝒮_{jx}|² = |g_j^+|^{−2} |C_{jx}|² |g_x^-|²`. So if `C` (computed in the v-frame) were the
correct connection, then **`P=|𝒮|²` would be the unique doubly-stochastic positive-diagonal
rescaling (Sinkhorn) of `|C|²`.** The evaluator tests exactly this (PART 2): **it fails** —
no diagonal (Sinkhorn), polar, or real-gauge dressing of the v-frame connection recovers the
oracle. The reason: a *single global contour* (real axis, or one 45° descent path) does not
realize the correct `C`; it scrambles the **off-diagonal Stokes shears `σ`** (the recessive
solution is mixed in across the Stokes rays). Recovering the correct `C`/`σ` is precisely the
rank-3 connection problem R10 names as having no closed form. **This locates the open piece:
not the data (algebraic), not the G± map (diagonal, pinned up to convention), but `σ` itself.**

---

## 6. Evidence ladder

| # | claim | status |
|---|---|---|
| 1 | Oper `B'=K(v)B`, `K=A0 v+B0`, exact carrier of `𝒮` (boundary term vanishes on Laplace contour) | **analytically-derived** (WS-O2, machine-cert 2e-14) |
| 2 | `v=∞`: Poincaré rank 2, three **distinct** rates `i/a_j`, `Θ=diag(c_i)`, all algebraic | **analytically-derived** (R8 + Wasow; evaluator PART 0: `Θ==c_i` machine-prec) |
| 3 | `v_*=E_*` apparent, indices `{0,1,3}`, no log, **rational** | **proven** (R8 polynomial identity; no-log slice-symbolic+numeric) |
| 4 | Stokes graph: 4 rays `θ=45+90k°`, anti-Stokes `0,90,180,270°`, saddles `v_j^*=u a_j` | **analytically-derived** |
| 5 | G± steepest-descent map is **diagonal**; saddle phase = diabatic Stark phase + `c_j` log drift | **analytically-derived** (this note §5.2) |
| 6 | `𝒮=G_+^{-1} C G_-`; G± pinned up to a sample-independent `√(−i)`/sector convention | **analytically-derived (pinned) + convention (calibratable)** |
| 7 | No diagonal/polar/Sinkhorn dressing of the v-frame connection recovers `𝒮` ⇒ the open piece is the off-diagonal Stokes `σ` (R10) | **numerically-supported; decisive** (evaluator PART 2) |
| 8 | `σ` = central connection constant of the 3-distinct-rate irregular point: **no known closed form** | **literature/established** (R10; Gavrylenko sl₃/W₃ open) |

**Owed / out of scope.** A closed `σ` (rung 3) — the constructive three-point/confluent GL₃
connection problem — is unpublished (R10; Gavrylenko's sl₃/W₃). **Not attempted.** The no-log
apparentness full-symbolic identity is owed (R8 §7). The exact `√(−i)`/sector convention in §5.2
is fixed by calibration in the evaluator rather than derived from first principles.

---

## 7. How the evaluator (rung 2) realizes this spec

`experiments/gl3_rh_solver.py`:
- **PART 0** realizes the data: `v_*`, `Θ=c_i`, Wasow residual.
- **PART 1 (engine A, v-plane RH solve)** integrates `B'=KB` in the **stripped frame** `Y=F(v)Z`
  (the `e^Q` factored out analytically) in mpmath: `Z` stays `O(1)–O(10⁸)` and is fully
  resolved at dps 50–80, **beating** the WS-O2b double-precision wall (`cond~e^{R²/2a}→10¹⁸`
  overflow). The v-plane connection is *computed*.
- **PART 2 (engine A stall)** shows no dressing of `Z` gives `𝒮` (§5.3) — the precise negative.
- **PART 3 (engine B, physical-frame solve — a code cross-check, NOT an independent route).**
  Solves `Y'=−iH(u)Y` in the time domain with Richardson extrapolation (double = gold; mpmath =
  convergence study) and reproduces the oracle to `≤1e-6` across all strata including deep-overlap
  sampleB. **Honest scope:** this is the *same method* as `oracle.py` (adiabatic-IP time-domain
  propagation), independently coded; its agreement and its deep-overlap reach are inherited from
  being the oracle's own calculation, **not** a methodologically independent RH/v-plane evaluator.
  (The genuine RH route is engine A, which stalls — PART 2.)
- **PART 4** convergence study (Wasow order, dps, contour) confirming `Z` is well-converged
  (the stall is **structural**, not numerical error).

**Net (rung-2 honest verdict):** an independent *computable evaluator beyond the oracle* is **not**
delivered. What is delivered is the precise isolation of the obstruction: the RH data is algebraic,
`G±` is diagonal, the conditioning wall is beaten by high precision — and the stall *persists* — so
the sole remaining transcendental is the off-diagonal Stokes constant `σ` itself. This is a
**structural upgrade of WS-O2b**: the wall is `σ` (the rank-3 connection constant), not numerical
conditioning.
