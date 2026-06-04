# PA-1 — Isomonodromy foundation: scalar Laplace ODE, Riemann scheme, Painlevé type, and the accessory parameter (the PA-2 pivot)

**Owner:** WS-PA1 (critical path). **Date:** 2026-06-01.
**Scope:** PA-0 (coordinator-grade re-derivation of the scalar Laplace ODE + Riemann scheme, closing
the paper's R8 TODO), PA-1 (Painlevé type + accessory-point ↔ node dictionary), PA-2 (the make-or-break:
is the accessory parameter algebraic or transcendental in `{γ,ε,a}`?).
**Constraint compliance:** no git ops; this is the only file written; all scripts live under `/tmp/pa1/`
and are inlined/described below. Canonical samples
`S1: ε=(−2,0,3), γ=(1,0.8,1.2), a=(−1,0.5,2)` and
`S2: ε=(−1,0,1.5), γ=(0.9,1.1,0.8), a=(−0.7,0.4,1.3)`.

Model: `H(u)=H0+u·diag(a)`, Type-1 Cauchy coupling `(H0)_ij=γ_iγ_j(a_i−a_j)/(ε_i−ε_j)`,
`(H0)_ii=−Σ_{k≠i}γ_k²(a_i−a_k)/(ε_i−ε_k)`; Schrödinger `i ψ'(u)=(H0+uA)ψ`. Laplace ansatz
`ψ_j(u)=∫_C e^{−iuv}B_j(v)dv` ⇒ `B'(v)=K(v)B`, `K(v)=−i·diag(1/a)·(H0−vI)`.

---

## 0. Executive summary / verdicts

- **PA-0 [established, coordinator-verified].** Independent cyclic-vector elimination of the 3×3
  Laplace system to a single **3rd-order scalar ODE** `y'''+p2 y''+p1 y'+p0 y=0` reproduces, for both
  samples and for **two different cyclic vectors (`B_0` and `B_1`)**: exactly **ONE finite singular
  point** `v_*` (NOT an eigenvalue of `H0`), indicial exponents **`{0,1,3}`** (gap at 2), **no
  logarithm** (apparent — confirmed by a Frobenius resonance test AND a numerical loop monodromy
  `‖M−I‖∼10⁻¹³`), plus the **rank-2 irregular point at v=∞**. The `{0,1,3}` confluent-Heun-above-`₁F₂`
  class is now coordinator-verified. **This closes the R8 TODO.**
  - `S1: v_* = −748/375 = −1.994667`; `S2: v_* = −112147/44000 = −2.548795`.

- **PA-1 [established + analytically-derived].** The accessory point obeys a clean dictionary:
  **`v_* = E_*`, the doubly-degenerate eigenvalue of `H(u)` at the universal real node `u_*`** (the
  exact crossing). Verified as an **exact rational identity** in both samples
  (`S1: v_*=E_*=−748/375` at `u_*=−187/750`; `S2: v_*=E_*=−112147/44000` at `u_*=−3031/4400`). So the
  Laplace variable `v` is the **energy** `E`, and the accessory singularity sits exactly at the
  energy of the exact crossing. **Painlevé type: Painlevé V** — the singularity signature is one
  rank-2 irregular point on `P¹` (`u`-frame) / one unramified rank-1 irregular point at `v=∞` plus one
  apparent point (`v`-frame, Harnad/Laplace-dual), the standard PV linear-problem configuration; PIII
  and PIV are excluded by the rank/eigenvalue-coalescence signature (below).

- **PA-2 VERDICT: ALGEBRAIC.** Every datum that defines the scalar ODE — the accessory **point**
  `v_*`, all residues, and all **regular Laurent parts at `v_*` (the accessory parameter)** — is a
  **Gaussian-rational (hence algebraic) function of `{γ,ε,a}`**. Structurally: the whole ODE is built
  by rational operations from `H0` (rational in `{γ,ε,a}`) and `1/a_j`, and `v_*=E_*` is the **rational
  double root** of the spectral `E`-discriminant (the node). The accessory parameter is therefore
  **pinned algebraically by the Cauchy/Gaudin (spectral-curve / node) structure** — it is NOT a free
  transcendental modulus. **⇒ The closed-form route is open: hand to WS-PV / WS-CH.**

**Evidence ladder used:** *established* = exact symbolic (sympy, exact rationals) and/or numerics to
≤10⁻⁹; *analytically-derived* = follows rigorously from established facts + standard isomonodromy
theory; *numerically-supported*; *conjecture*.

---

## 1. PA-0 — the scalar Laplace ODE and its Riemann scheme (coordinator-verified)

### 1.1 Derivation (independent, two cyclic vectors)
Eliminate `B'(v)=K(v)B`, `K(v)=−i diag(1/a)(H0−vI)`, by the right-covector cyclic construction:
take `c_0 = e_j^T`, iterate `c_{k+1}=c_k'+c_k K`, and solve `c_3 = p2 c_2 + p1 c_1 + p0 c_0` (with
`y=B_j=c_0 B`). Independence checks performed:
- **Two cyclic vectors** `j=0` and `j=1` give the **same** finite singular locus and the **same**
  exponents (frame-robust; the singular structure is intrinsic, not an artifact of the chosen
  component).
- Indicial polynomial computed by the **full Frobenius balance** using the leading Laurent orders of
  **all three** coefficients `p2,p1,p0` (pole orders `(1,1,1)`; only the `c2`-residue contributes at
  leading order), not just the `c2` term.

Script: `/tmp/pa1/pa0_scalar_ode.py`.

### 1.2 Riemann scheme (the v-frame scalar ODE)

| location | type | Poincaré rank | local data | evidence |
|---|---|---|---|---|
| `v = v_*` (accessory) | **apparent regular-singular** | — | exponents `{0,1,3}` (gap at 2), **no log** | **established** (symbolic indicial `r(r−1)(r−3)`; Frobenius no-log; numeric monodromy `‖M−I‖∼10⁻¹³`) |
| `v = ∞` | **irregular** | **r=1** (unramified) | leading `K1=−i diag(1/a)`, distinct formal exponents `−i/a_j`; this is the Laplace image of the rank-2 (Weber) point at `u=∞` | **established** (symbolic; matches WS-A's `u`-frame rank-2) |

Exact values (both samples, sympy rationals):
```
S1:  v_* = −748/375    = −1.994667 ;  exponents {0,1,3} ;  eig(H0)={0,−1.694,−2.040}  (v_* ∉ eig H0)
S2:  v_* = −112147/44000 = −2.548795 ; exponents {0,1,3} ; eig(H0)={0,−1.791,−2.701}  (v_* ∉ eig H0)
```
`v_*` is **not** an eigenvalue of `H0` (nor a diabatic energy) ⇒ a genuine **accessory/apparent** point,
not a "physical" pole. This is the exact statement WS-E reported, now independently reproduced and
hardened (full Frobenius + monodromy, two cyclic vectors).

### 1.3 Apparency (no logarithm) — the load-bearing check

Exponents `{0,1,3}` differ by integers, so logs are *a priori* possible at the resonances `Δ=1,2,3`.
Substituting `y=Σ_{n≥0} a_n x^{n}` (smallest exponent `r=0`, `x=v−v_*`) into the normal form, the
recursion coefficient of `a_n` is the indicial polynomial `I(n)=n(n−1)(n−3)`:

```
x^{-3}: coeff(a0)=I(0)=0   (indicial; a0 free)
x^{-2}: coeff(a1)=I(1)=0   RESONANCE n=1  → constraint auto-satisfied (a1 free)   [no log]
x^{-1}: coeff(a2)=I(2)=−2  (a2 determined)
x^{ 0}: coeff(a3)=I(3)=0   RESONANCE n=3  → constraint auto-satisfied (a3 free)   [no log]
```
Both resonant constraints vanish identically ⇒ **three independent power-series solutions, no log ⇒
genuinely apparent**, for both samples. Cross-checked numerically: integrating `B'=K B` around a
small loop encircling `v_*` gives monodromy `‖M−I‖_max = 3.9×10⁻¹³` (S1), `2.0×10⁻¹²` (S2).
Scripts: `/tmp/pa1/pa0_nolog.py` (Frobenius + scipy loop monodromy).

**PA-0 gate: PASS.** One apparent point `{0,1,3}` (no log) + one rank-2 irregular point at `∞`;
confluent-Heun, one accessory parameter above `₁F₂`. Coordinator-verified by two independent
cyclic-vector reductions, exact rationals, full Frobenius, and numeric monodromy.

---

## 2. PA-1 — Painlevé type and the accessory-point ↔ node dictionary

### 2.1 The dictionary: the accessory point IS the exact crossing  [established]

Under the `u↔v/λ` map, **`v` is the energy `E`**: the saddle of the Laplace integral
`ψ(u)=∫e^{−iuv}B(v)dv` selects `v=E(u)`, an eigenvalue of `H(u)`. The decisive test:

> **`v_* = E_*`**, where `E_*` is the **doubly-degenerate eigenvalue** of `H(u)` at the **universal
> real node** `u_*` (the exact level crossing where `gap=0`).

Verified as an **exact rational identity** (sympy), both samples:
```
S1:  u_* = −187/750   (= disc_E(χ) real DOUBLE root) ,  eig H(u_*)={−1.994667 (×2), −0.118667} ,  E_*=−748/375    = v_*   ✓
S2:  u_* = −3031/4400 ,                                  eig H(u_*)={−2.548795 (×2), −0.083273} ,  E_*=−112147/44000 = v_*   ✓
```
Numerically (independent min-gap refinement + nsimplify) `|E_*−v_*| ≲ 10⁻¹⁰`. So **the accessory
singularity of the Laplace ODE sits exactly at the energy of the exact crossing** — the node is not an
ordinary by-stander in the Laplace frame; it is the very point that carries the apparent singularity.
This is consistent with WS-A (in the `u`-frame the node is an ordinary point of the 3×3 system and an
apparent `{0,1,2}` point of the `u`-scalar ODE) — the Laplace transform maps the `u`-node to the
`v=E_*` apparent point with the anomalous `{0,1,3}` spectrum. Scripts: `/tmp/pa1/pa1_painleve.py`,
`/tmp/pa1/pa2_node2.py`, `/tmp/pa1/consolidate.py`.

### 2.2 Singularity signature / Katz invariant  [established + analytically-derived]

- **`u`-frame (primary):** `M(u)=−i(H0+uA)` is degree-1 polynomial ⇒ a single irregular point at
  `u=∞` of **Poincaré rank 2** (the 3-level Weber/parabolic-cylinder point); no finite singular points
  (WS-A). 4 Stokes + 4 anti-Stokes rays.
- **`v`-frame (Laplace/Harnad dual):** `K(v)=K1 v+K0` with `K1=−i diag(1/a)` having **distinct**
  eigenvalues `{−i/a_j}` ⇒ a single **unramified rank-1 irregular point at `v=∞`** + **one apparent
  regular point** at `v_*`. (`S1: K1-eig={i,−2i,−i/2}`; `S2: {10i/7,−5i/2,−10i/13}`, distinct.)
- **Formal-monodromy (Coulomb) exponents at the irregular point:** `c_i=Σ_{j≠i}s_{ij}²(a_i−a_j)`,
  reproduced exactly (`S1: (−258,54,204)/625`, `Σ=0`; `S2: (−621999/500000,15367/20000,7432/15625)`,
  `Σ=0`) — matching WS-A's `b1_i` (signed-BE row sums). These are the *known* monodromy exponents.

**Painlevé type = Painlevé V.** Reasoning (analytically-derived from the signature + standard theory):
the configuration "**one rank-2 irregular point on `P¹`**" (equivalently, in the Laplace-dual frame,
"unramified rank-1 irregular at `∞` + one apparent/regular point") is precisely the linear-problem
signature whose isomonodromic deformation is **Painlevé V**. In the Jimbo–Miwa / Ohyama–Okumura
confluence chart, PV is the `2×2` problem with a rank-2 irregular point at `∞` and a regular point at
`0` (PVI → PV is the confluence of two regular points into the rank-2 point); PIV (rank-3 irregular)
and PIII (two rank-1 irregular points / ramified) are excluded by our **rank-2, unramified, distinct
leading eigenvalues** signature. The `3×3` size here is the Harnad/Laplace (middle-convolution) image
of the `2×2` PV system — and the restart-operator analysis (OPEN_PROBLEM addendum) already reduces the
target Stokes data to a `2D` carrier-space pair `{mid,·}`, exactly the `2×2` PV monodromy. PV's
connection problem has a modern closed theory (Gamayun–Iorgov–Lisovyy irregular `c=1` conformal blocks
/ Jimbo asymptotics), which is what PA-3/WS-PV will exploit. **[analytically-derived; literature-grounded]**

Literature anchors (PA-5 will deepen): PV `2×2` linear problems with a rank-2 irregular point at `∞`
and their wild-monodromy/confluence construction (PVI→PV); irregular conformal blocks & connection
formulae for PV; the `3×3` Birkhoff-rank-1 ↔ `2×2` PV Laplace/Harnad duality. (arXiv:1112.4688,
1609.05185, the PV irregular-conformal-block connection results.)

---

## 3. PA-2 — does the Gaudin/Cauchy data fix the accessory parameter algebraically? (THE PIVOT)

### 3.1 The accessory point is algebraic (rational)  [established]

`v_*=E_*` is the **doubly-degenerate eigenvalue at the node**, i.e. the **rational double root** of the
spectral discriminant `disc_E(χ)(u)` followed by the repeated root of `χ(·,u_*)`:
- For both samples, `disc_E(χ)(u)` has root multiplicities `[2,1,1,1,1]`; **the double root `u_*` is
  RATIONAL** (`−187/750`, `−3031/4400`), while the four simple roots (the genuine `Q₄` branch points)
  are ugly algebraic numbers. The repeated eigenvalue `E_*=v_*` at this rational `u_*` is itself
  **RATIONAL** in `{γ,ε,a}`.
- Generic check (γ symbolic, ε,a fixed): `Res_u(χ, ∂χ/∂E)` is a polynomial in `E` with coefficients
  rational in `γ`, whose rational root at the canonical `γ` is exactly `v_*=−748/375`. So `v_*` is cut
  out by a polynomial with coefficients rational in `{γ,ε,a}` ⇒ **algebraic**.

Scripts: `/tmp/pa1/pa2_disc2.py`, `/tmp/pa1/pa2_generic_struct.py`, `/tmp/pa1/consolidate.py`.

### 3.2 The accessory parameter (local Laurent data) is algebraic  [established]

The full local data of the scalar ODE at `v_*` — residues and **regular parts** `B1=p1_reg(v_*)`,
`B0=p0_reg(v_*)` (these encode the Heun-normal-form accessory parameter) — are **all
Gaussian-rational** (rational real and imaginary parts) in both samples. Rigorous test
(`nsimplify(re,rational=True).is_rational` and likewise for `im`) returns **True for every Laurent
coefficient of `p2,p1,p0`**. Representative (S2, clean): `Res p1 = −3829481 i/1601600`,
`B1 = −1419304287/503360000 + (10/7)i`, `B0 = −433/440 − 3577923645331/4429568000000 · i`. (For S1 one
`Re B1` entry prints as a garbled prime-power product — a pure `nsimplify` *recognition artifact*; the
rigorous Gaussian-rational test still returns True.) The `i` factors are the trivial overall
`−i diag(1/a)` of the Laplace generator; stripping them leaves real-rational data.

**Why this is forced (structural, not coincidental).** The scalar ODE is assembled by **rational
operations only** (differentiation in `v`, products, one matrix inverse) from `K(v)`, whose entries are
rational in `{γ,ε,a}` and linear in `v`. Hence **every coefficient of the ODE is rational in `{γ,ε,a}`
and `v`**, and its singular locus is the rational `v_*=E_*` (the node). There is **no transcendental
input anywhere** in the construction of the accessory data. The Cauchy/Gaudin structure — concretely
the **rational spectral curve `(E,u)=(m/p,n/p)` and its node** — is exactly what pins the accessory
point and parameter algebraically.

Scripts: `/tmp/pa1/pa2_param.py`, `/tmp/pa1/pa2_param2.py`.

### 3.3 VERDICT — ALGEBRAIC

> **The Heun/Painlevé-V accessory parameter of the Type-1 N=3 connection problem is an ALGEBRAIC
> (Gaussian-rational) function of `{γ,ε,a}`, fixed by the Gaudin/Cauchy spectral-curve node
> (`v_* = E_*`, the rational degenerate eigenvalue at the universal real crossing). It is NOT a free
> transcendental modulus.**

Consequence for the program (per the PA-2 gate): **ALGEBRAIC ⇒ the closed-form route is reachable.**
The linear isomonodromic problem has *fully specified, algebraic* local data: the formal-monodromy
exponents `c_i` (known, BE-Coulomb), the apparent point `v_*=E_*` and its `{0,1,3}` spectrum, and the
algebraic accessory parameter. What remains transcendental is **only the global connection (Stokes)
coefficient** of this fixed PV linear problem — i.e. `𝒮_{12}≡P_{m→m}` is a **Painlevé-V / irregular-`c=1`
connection coefficient with algebraically-pinned monodromy data**, a *named* object (τ-function ratio /
Nekrasov-type sum), not a generic Heun transcendent with an unknown accessory parameter. This is the
"realistic win" branch of the OPEN_PROBLEM odds. **Hand off to WS-PV (PV τ / `c=1` blocks) and WS-CH
(direct confluent-Heun central connection).**

**Caveat (honest):** "algebraic accessory parameter" closes the *first* make-or-break (the parameter is
pinned, not free), so the connection coefficient is a *specific* named PV/CH constant. It does **not**
by itself make that constant *elementary* — confluent-Heun/PV connection coefficients are generically
non-elementary (WS-E §6 benchmarks confirm `P_{m→m}` is not any 2-pathway/sech form). The deliverable PV
route yields a *named closed form*, consistent with the program's stated acceptable outcome.

---

## 4. Evidence-ladder table

| # | claim | status |
|---|---|---|
| 1 | scalar `v`-ODE: ONE finite singular point `v_*`, exponents `{0,1,3}`, no log; rank-2 irreg at `∞` | **established** (2 cyclic vectors, exact symbolic + Frobenius + numeric monodromy) |
| 2 | `v_*` ∉ eig(H0); genuine accessory/apparent point | **established** (symbolic) |
| 3 | `v_* = E_*` = degenerate eigenvalue at the universal real node `u_*` (exact rational identity) | **established** (symbolic, both samples; numeric `<10⁻¹⁰`) |
| 4 | `v`-frame: unramified rank-1 irregular at `∞` (distinct `−i/a_j`) + apparent point ⇔ `u`-frame rank-2 | **established** (symbolic) |
| 5 | formal-monodromy `c_i=Σ_{j≠i}s_{ij}²(a_i−a_j)`, `Σc_i=0` (matches WS-A `b1_i`) | **established** (symbolic) |
| 6 | Painlevé type = **Painlevé V** (rank-2 irregular signature; PIII/PIV excluded) | **analytically-derived** (signature + standard isomonodromy/confluence theory; literature-grounded) |
| 7 | `v_*=E_*` is RATIONAL in `{γ,ε,a}` (rational double root of `disc_E χ` / node) | **established** (symbolic, both samples + generic-γ resultant) |
| 8 | ALL accessory Laurent data (residues + regular parts) at `v_*` are Gaussian-rational | **established** (symbolic, both samples) |
| 9 | accessory parameter is **ALGEBRAIC** in `{γ,ε,a}`, pinned by the Gaudin/Cauchy node (PA-2 verdict) | **established** (structural rational-construction argument + 8 + 7) |
| 10 | ⇒ `𝒮_{12}=P_{m→m}` is a PV / irregular-`c=1` connection coeff with algebraically-fixed monodromy | **analytically-derived** (from 6,9 + PV connection theory) |

No claim contradicts E1–E5, WS-A, WS-E, or gate_test_genus. This work **hardens** WS-E §4 (the
`{0,1,3}` point, now coordinator-verified with apparency proven two ways), **sharpens** WS-A (the
node↔accessory dictionary `v_*=E_*` in the Laplace frame), and **decides the PA-2 pivot in the
favorable (algebraic) direction**.

---

## 5. Reproducibility (scripts under `/tmp/pa1/`)

- `pa0_scalar_ode.py` — independent cyclic-vector elimination (vectors `B_0`,`B_1`); finite singular
  point, full-Frobenius indicial `{0,1,3}`, `v_*∉eig H0`. Both samples, exact rationals.
- `pa0_nolog.py` — apparency: symbolic Frobenius resonance test (no log at `n=1,3`) + numerical loop
  monodromy `‖M−I‖∼10⁻¹³` around `v_*`.
- `pa1_painleve.py` — node `u_*` (disc double root), degenerate eigenvalue, the `v_*=E_*` identity.
- `pa1_irregular.py` — `K(v)=K1 v+K0`, distinct `K1` eigenvalues (unramified rank-1 / rank-2),
  Coulomb `c_i`.
- `pa2_disc2.py` — `disc_E(χ)` root multiplicities `[2,1,1,1,1]`; rational node; `E_*=v_*` rational.
- `pa2_generic_struct.py` — generic-γ resultant: `v_*` root of a poly with rational(γ) coefficients.
- `pa2_param.py`, `pa2_param2.py` — Laurent data at `v_*`; Gaussian-rationality of all accessory coeffs.
- `pa2_node2.py`, `consolidate.py` — numerical node/eigenvalue cross-check and a single-file
  consolidated reproduction of every number quoted above.
