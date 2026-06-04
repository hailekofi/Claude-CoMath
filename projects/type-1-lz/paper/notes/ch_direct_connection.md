# WS-CH — The direct confluent-Heun central connection problem for `S₁₂`

> **Notation note (coordinator):** this deliverable's `S₁₂` is DEPRECATED — read it as the middle
> survival `P_{m→m}=|𝒮_{mm}|²` (`m=slope-middle`) `=|C_{mm}|²`, built from the off-diagonal Stokes
> multiplier `σ`. See `NOMENCLATURE.md`. (The benchmarked numbers are `|𝒮_{mm}|²`, correct.)


**Owner:** WS-CH (parallel direct confluent-Heun connection-coefficient track, independent of
the Painlevé-V/τ route WS-PV). **Date:** 2026-06-01.
**Inputs read:** `NOMENCLATURE.md`, `paper/program/OPEN_PROBLEM.md`, `paper/working_sessions/ws_e_junction_Smatrix.md`,
`paper/working_sessions/ws_a_riemann_scheme.md`, `experiments/oracle.py`, `experiments/ws_e_laplace_class.py`.
**Constraint compliance:** no git ops; the only file written under `paper/` is this one; all
scratch lives under `/tmp/wsch/`.

Target (NOMENCLATURE-strict): the off-diagonal Stokes coefficient `S₁₂ = P_{m→m}` (middle-slope
survival) — the off-diagonal Stokes multiplier `σ_{mid,·}` (or its composition) of the **single
rank-2 irregular point** of the Type-1 N=3 connection problem. `s_ij = γ_iγ_j/(ε_i−ε_j)`;
BE exponent `= s_ij²|a_i−a_j|`; signed-BE Coulomb coefficient `c_i = Σ_{j≠i} s_ij²(a_i−a_j)`.

---

## 0. Executive verdict (honest, evidence-tagged)

> **`S₁₂` is the central connection coefficient of a genus-0 confluent-Heun-class problem: one
> finite apparent singularity (Frobenius exponents `{0,1,3}`, gap at 2) plus one rank-2 irregular
> point. The local solution bases at both ends are written below in closed form (formal Thomé /
> WKB at ∞, Frobenius at the apparent point). The CONNECTION COEFFICIENT itself has NO classical
> closed form — this is a theorem-level fact for confluent-Heun connection problems (Bühring;
> Wolf; Lisovyy–Naidiuk) — but it DOES have a named convergent representation: the
> Lisovyy–Naidiuk (2022) perturbative connection formula, a convergent series in the accessory
> parameter, equivalent to quasiclassical Virasoro conformal blocks (the same object WS-PV reaches
> from the Painlevé-V side). The practical deliverable is the central-connection matrix computed
> directly (lifted to the 3×3 system in the flux-normalised / adiabatic frame), benchmarked to the
> gold oracle.**

This **agrees with and sharpens** WS-A and WS-E: the negative-with-structure verdict (no elementary
/ no ₁F₂ closed form) is now backed by the *connection-theoretic* reason (confluent-Heun connection
coefficients are provably non-elementary), and is paired with the *positive* outcome that the
relevant coefficient is a **named, convergent, computable** object (Lisovyy–Naidiuk series ↔ WS-PV's
τ-function), not merely "unknown."

Evidence ladder used: **[established]** = exact symbolic / reproduced numerically; **[derived]** =
analytic with controlled approximation; **[literature]** = cited published result; **[conjecture]**.

---

## 1. The two local solution bases (closed form)

### 1a. Frame and ODE

Two equivalent scalar reductions carry the same connection data:

- **u-frame** (WS-A §1b): scalar 3rd-order ODE from `i ψ'=(H0+uA)ψ`, one finite **apparent**
  singularity at the node `u* = −(node)`, exponents `{0,1,2}` (no log), + rank-2 irregular at `u=∞`.
- **v-frame / Laplace** (WS-E §4, the brief's frame): `ψ_j(u)=∫_C e^{−iuv}B_j(v)dv`, `B'(v)=K(v)B`,
  `K(v) = −i·diag(1/a)·(H0 − vI)`. Eliminating to a scalar 3rd-order ODE `L[b]=0` gives **one finite
  apparent singularity** `v_*` with Frobenius exponents **`{0,1,3}`** (gap at 2) + rank-2 irregular
  point at `v=∞`. **[established, reproduced both anchors]** (`/tmp/wsch`, `ws_e_laplace_class.py`):

  | anchor | `v_*` | indicial exponents | residue `A` |
  |---|---|---|---|
  | canonical | `−748/375 = −1.99467` | `{0,1,3}` | `1` |
  | sampleB | `−112147/44000 = −2.54880` | `{0,1,3}` | `1` |

  `v_*` is NOT an eigenvalue of `H0` (canonical `eig H0={−51/25,−847/500,0}`) ⇒ a genuine accessory
  point. The exponents `{0,1,3}` are the "one accessory above ₁F₂" signature: a 3rd-order ODE with
  exponents `{0,1,2}` and no accessory would be ₁F₂/Kampé-de-Fériet (BBGY); the **gap at 2** is the
  extra accessory parameter.

The brief specifies the v-frame `{0,1,3}`. (The u-frame's `{0,1,2}` apparent point is the same
accessory content in the dual variable; both are "one accessory above ₁F₂." For the connection
*computation* we lift to the 3×3 system, where the apparent point is an ordinary point and the only
singularity is the rank-2 irregular one — WS-A's recommended frame.)

### 1b. Frobenius basis at the apparent singularity `v_*` (exponents `{0,1,3}`)

Local solutions `b(v) = (v−v_*)^{ρ} Σ_{n≥0} a_n (v−v_*)^n`, `a_0=1`, with `ρ ∈ {0,1,3}`:

- `ρ=0`: `b^{(0)} = 1 + a_1(v−v_*) + …`  (holomorphic).
- `ρ=1`: `b^{(1)} = (v−v_*) + …`.
- `ρ=3`: `b^{(3)} = (v−v_*)³ + …`.

The gap at `ρ=2` means the second solution can in principle carry a logarithm; **WS-E/WS-A report no
log** (the singularity is apparent — trivial local monodromy, the recurrence closes with integer
exponents and no resonance obstruction). [established, both anchors — `ws_e_laplace_class.py` finds
exponents `{0,1,3}` and apparent (no-log) character.] The `a_n` are fixed by the recurrence from
`L`; they are **rational in `{γ,ε,a}`** and computable to any order — this is the data that feeds
the convergent connection series (§3).

### 1c. Thomé / WKB formal solutions at the rank-2 irregular point (closed form)

At `u=∞` (equivalently `v=∞`) the formal fundamental solutions are (WS-A §3, re-derived here):

```
ψ_j(u) ~ T_j(u) · exp(−i ∫^u E_j du')
       = [v_j(u)+O(1/u)] · exp(−i[ a_j u²/2 + (H0)_jj u + c_j log u ]),
```
per branch `j`, with **the algebraic-prefactor exponent = the signed-BE Coulomb coefficient:**
```
c_j = Σ_{k≠j} (H0)_jk²/(a_j−a_k) = Σ_{k≠j} s_jk²(a_j−a_k)   (signed BE row sum).
```
**Formal monodromy — precise statement (convention-fixed).** In the `exp(−i∫E)` convention the
algebraic prefactor is `u^{−i c_j}`, so the monodromy `u→e^{2πi}u` multiplies branch `j` by
`e^{2π c_j}` (real, since `c_j` is real). Because `c_j` is the *signed* BE row sum, `|e^{2π c_j}|`
**reproduces the BE survivals directly**: canonical `e^{2π c_lo}=e^{2π(−0.4128)}=0.07474 = P_lo`
(the exact extreme survival), etc. — the diagonal/formal part of the restart operator. The
OPEN_PROBLEM shorthand "`e^{2πi c_i}`" is this object with the `−i` from the `exp(−i∫E)` phase
absorbed (equivalently, the formal-monodromy *exponent* is `c_j`, and `Σ_j c_j = 0` is the trace /
double-stochasticity condition). **[established]** Verified independently in the v-frame by formal
diagonalisation of `B'=K(v)B` (`/tmp/wsch/formal_infty.py`): the 2nd-order perturbative diagonal
shift gives `c_j^{(v)} = −i·(signed BE)`, i.e. the v-frame prefactor `v^{−i·c_j}` — the SAME
monodromy. Numerically (exact rationals, `c_j` = signed-BE row sum = WS-A's `b1_j`):

| anchor | `c_j` = signed-BE row sum (`= b1_j` of WS-A) |
|---|---|
| canonical | `(−258/625, 54/625, 204/625) = (−0.4128, 0.0864, 0.3264)` |
| sampleB | `(−621999/500000, 15367/20000, 7432/15625) = (−1.2440, 0.7684, 0.4757)` |

`Σ_j c_j = 0` (the double-stochasticity / trace condition). [established, exact symbolic]

**Stokes geometry** (WS-A §3, rank-2 ⇒ Weber): 4 anti-Stokes rays at `arg u ∈ {0,π/2,π,3π/2}` and
4 Stokes rays at `arg u ∈ {π/4,3π/4,5π/4,7π/4}`; `2r=4` sectors. The **real axis is an anti-Stokes
line** (balanced/oscillatory) — the natural connection line `u:−∞→+∞`. The two carrier 2-spaces
`{mid,lo}` and `{mid,hi}` (WS-G's 12×13 joint) host the two Stokes shears whose non-commutative
composition is `S₁₂`. [established / conj per OPEN_PROBLEM]

---

## 2. The connection-matrix setup

The central connection problem is: express the incoming Thomé/WKB basis (`u→−∞`, sector containing
`arg u = π`) in terms of the outgoing basis (`u→+∞`, sector `arg u = 0`):
```
ψ^{in}_j(u)  =  Σ_k  C_{jk} ψ^{out}_k(u),       P_{j→k} = |C_{jk}|².
```
The connection matrix factors canonically (the restart structure of OPEN_PROBLEM §Addendum):
```
C  =  (formal monodromy, diagonal)  ⋉  (Stokes matrices, off-diagonal)
    =  M_form  ·  S^{(mid,hi)}  ·  S^{(mid,lo)}   (schematically),  M_form = diag(formal monodromy),
```
where each `S^{(mid,·)}` is a unipotent Stokes shear in the 2-D carrier space sharing the middle
sheet. **`S₁₂ = P_{m→m}` is the (mid,mid) entry of `|C|²`**. The formal-monodromy diagonal `M_form`
sets the BE *magnitudes* (its `|·|²` on the extreme rows = the exact BE survivals); the off-diagonal
content of `P_{m→m}` comes from the two Stokes shears `σ_{mid,·}` and, crucially, their non-commutative
*interference* through the shared middle sheet — this is what makes `P_{m→m}` a genuine three-crossing
coherence, not a product `p₁p₂` (WS-E §6).

**Flux normalisation (why `|C|²` is doubly-stochastic).** The raw *diabatic* Thomé basis is NOT
flux-normalised (its `T(u)=I+T₁/u+…` prefactor mixes channels), so a naive diabatic connection
matrix is not stochastic. The correct basis is the **adiabatic/WKB one** `v_j(u)exp(−i∫E_j)`, in
which the connection matrix is unitary and `|C|²` is automatically doubly-stochastic. We compute
`C` in that frame (§4); it is the central connection coefficient of the confluent-Heun-class
problem, lifted to the 3×3 system where the apparent point is ordinary.

---

## 3. Is there a closed / convergent form? — known confluent-Heun connection results

**Theorem-level fact (the honest blocker).** Confluent-Heun central connection coefficients are
**generically not elementary / not classical special functions.** They are represented as
[literature]:
- a **finite determinant** (size 2×2/3×3/4×4) whose entries are Taylor series with recursively
  available coefficients (Wolf 1998, *Math. Nachr.*; Bühring, double-confluent Heun characteristic
  exponent & connection formulae), or
- an **asymptotic / convergent series** in the recursively-known coefficients of the formal power
  series at the irregular point.

So a closed *value* of `S₁₂` in elementary functions of `{γ,ε,a}` does **not** exist (consistent
with WS-A/WS-E and the PA-2-transcendental branch of OPEN_PROBLEM).

**The named convergent representation (the positive result).** Lisovyy–Naidiuk, *Perturbative
connection formulas for Heun equations*, J. Phys. A (2022), arXiv:2208.01604 [literature]: for the
**confluent and reduced-confluent Heun equation**, the connection matrix is given as a **convergent
series in the accessory parameter**, computable to arbitrary order from the Frobenius/Thomé data of
§1 — and it **equals the quasiclassical (`c→∞`) Virasoro 4-point conformal block / the Painlevé-V
isomonodromic τ-function ratio** (the Bonelli–Iossa–Panea Lichtig–Tanzini conjecture, confirmed
there). This is **exactly the same named object WS-PV reaches from the Painlevé-V/τ side** — the two
tracks converge on it. The accessory parameter here is the residue/coefficient data of the `{0,1,3}`
apparent point (§1b), which the Type-1 Cauchy/Gaudin structure fixes *algebraically* in `{γ,ε,a}`
(this is the WS-PA1/PA-2 pivot; the apparent point `v_*` and its residue `A=1` are rational — a
strong indication the accessory data is algebraic, supporting the "named closed form is in reach"
branch).

**Concrete form of `S₁₂`.** Putting §1–§2 together, `S₁₂ = P_{m→m}` is the `(mid,mid)` modulus-squared
of the central connection matrix `C`, whose off-diagonal Stokes multipliers `σ_{mid,·}` are given by
the Lisovyy–Naidiuk convergent series in the accessory parameter with:
- monodromy exponents = the formal `c_j = Σ_{k≠j} s_jk²(a_j−a_k)` (§1c, KNOWN);
- accessory parameter = the `{0,1,3}` apparent-point data of `L` (rational `v_*`, residue `A=1`);
- argument data = the two window actions (BE exponents) + cross-ratio of the four turning points
  (WS-E §3,§5).

This is the **maximal closed-form statement**: a named (CFT/τ) convergent representation, not an
elementary value.

---

## 4. Benchmark — the direct central-connection computation vs the gold oracle

Computed `C` directly by transporting the flux-normalised (adiabatic/WKB) Thomé frame from `u=−R`
to `u=+R` along the real (anti-Stokes) axis and projecting onto the outgoing frame — the central
connection problem lifted to the 3×3 system — with Richardson extrapolation in `R`. Script:
`/tmp/wsch/ch_connection2.py`. Oracle: `experiments/oracle.py` (T=120 gold).

**Result (two independent runs agree, both reproduce the oracle to ≤4×10⁻⁷):**

| anchor | `P_{m→m}` (WS-CH connection) | `P_{m→m}` (oracle gold) | `|diff|` |
|---|---|---|---|
| canonical | `0.214724` | `0.214724` | `3.95×10⁻⁷` |
| sampleB | `0.021018` | `0.021018` | `1.56×10⁻⁷` |

Both runs give doubly-stochastic `|C|²` (row/col sums `=1` to `<10⁻⁵`) and the full diabatic matrix:

```
canonical  P =                      sampleB  P =
[[0.128628 0.056195 0.815177]       [[0.050359 0.002868 0.946773]
 [0.675196 0.214724 0.110080]        [0.926158 0.021018 0.052824]
 [0.196176 0.729081 0.074743]]       [0.023483 0.976114 0.000403]]
```
(rows = incoming diabatic slope channel, cols = outgoing; `[mid,mid]` = `P_{m→m}`.) **[established]**

Two independent solvers cross-check: `ch_connection3.py` (fast adiabatic interaction-picture,
Richardson `R∈{40,80}`, ~15 s) and `ch_connection2.py` (direct fundamental-matrix transport,
Richardson `R∈{120,240}`) agree to `<10⁻⁵` with each other and with the gold oracle. This confirms
the §1–§2 local-solution / connection-matrix construction *is* the object the oracle computes: the
central connection coefficient of the rank-2 irregular point, with `P_{m→m}` its `(mid,mid)` modulus-
squared. The deep-overlap sampleB stratum (`~10×` enhancement over the incoherent `p₁p₂≈0.000163`)
is reproduced — the three-crossing coherence is captured by the connection matrix, not by any
product/two-path form (WS-E §6). **[established]**

> Caveat: this benchmark validates the connection-matrix *framework and its numerical value*; it does
> NOT supply a closed-form symbolic `S₁₂`. The symbolic content is §3 (the named convergent
> Lisovyy–Naidiuk representation); the benchmark is the "computable `P_{m→m}`" deliverable (PA-4).

---

## 5. Status ladder

| # | claim | status |
|---|---|---|
| 1 | v-frame scalar ODE: 1 apparent sing. `v_*`, exponents `{0,1,3}` (gap at 2), no log, + rank-2 irregular at ∞ | **established** (symbolic, both anchors) |
| 2 | Frobenius basis at `v_*` (exponents `{0,1,3}`), coeffs rational in `{γ,ε,a}` | **established** |
| 3 | Thomé/WKB basis at ∞; formal-monodromy exponent `c_j=Σ s_jk²(a_j−a_k)` (signed BE); prefactor `u^{−ic_j}`, `\|e^{2πc_j}\|`=BE survivals; `Σc_j=0` | **established** (exact, both anchors; v-frame cross-check) |
| 4 | connection-matrix factorisation `C = diag(e^{2πi c_j})⋉ S^{(mid,hi)}S^{(mid,lo)}`; `S₁₂=|C|²_{mid,mid}` | **derived** (restart structure, OPEN_PROBLEM) |
| 5 | confluent-Heun connection coeff is NOT elementary/classical (determinant/series only) | **literature** (Wolf 1998; Bühring) |
| 6 | `S₁₂` HAS a named convergent representation = Lisovyy–Naidiuk series ↔ quasiclassical Virasoro block ↔ Painlevé-V τ | **literature** (arXiv:2208.01604) — *the WS-CH/WS-PV convergence point* |
| 7 | direct central-connection computation reproduces the oracle `P_{m→m}` | *(see §4 benchmark)* |

No claim contradicts WS-A/WS-E/E1–E5. This note **adds the connection-theoretic backing** for the
non-elementary verdict and **names the convergent form** of `S₁₂` (Lisovyy–Naidiuk / quasiclassical
conformal block), the same object WS-PV targets.

---

## 6. References

- O. Lisovyy, A. Naidiuk, *Perturbative connection formulas for Heun equations*, J. Phys. A **55**
  (2022) 434005, arXiv:2208.01604. — confluent-Heun connection matrix as a convergent series in the
  accessory parameter = quasiclassical Virasoro conformal block (≡ Painlevé-V τ).
- G. Wolf, *On the Central Connection Problem for the Double Confluent Heun Equation*, Math. Nachr.
  **195** (1998) 267. — connection coefficients as finite determinants of Taylor series.
- W. Bühring, *The double confluent Heun equation: characteristic exponent and connection formulae*.
- G. Bonelli, C. Iossa, D. Panea Lichtig, A. Tanzini — Heun connection matrix ↔ quasiclassical
  Virasoro conformal blocks (the conjecture confirmed by Lisovyy–Naidiuk).
- Litvinov et al. / *Accessory parameters in confluent Heun equations and classical irregular
  conformal blocks*, Lett. Math. Phys. (2021), arXiv (s11005-021-01400-6).

## 7. Reproducibility (scratch under `/tmp/wsch/`)

- `formal_infty.py` — formal diagonalisation at the rank-2 irregular point; `c_j = signed-BE`
  (v-frame), `Σc_j=0`, both anchors.
- `ws_e_laplace_class.py` (in `experiments/`) — v-frame scalar ODE, `v_*`, `{0,1,3}` exponents.
- `ch_connection2.py` — direct central-connection matrix (adiabatic/WKB frame), Richardson in `R`,
  benchmark vs `experiments/oracle.py`.
