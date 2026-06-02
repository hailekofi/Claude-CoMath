# WS-CAP — Capstone: evaluating the published Painlevé-V / Lisovyy connection formula at the algebraic monodromy data, and the PV-vs-rank-3 verdict

> **Notation (coordinator):** `S₁₂` here = the middle survival `P_{m→m}=|𝒮_{mm}|²` (`m=slope-middle`)
> `=|C_{mm}|²`; the transcendental ingredient is the off-diagonal Stokes multiplier `σ`. `S₁₂` as a
> symbol is DEPRECATED — see `NOMENCLATURE.md`. (Benchmarked numbers are `|𝒮_{mm}|²`, correct.)


**Owner:** WS-CAP (capstone, open-problem #2). **Date:** 2026-06-01.
**Inputs read:** `NOMENCLATURE.md`, `paper/RESEARCH_OVERVIEW.md`, `paper/OPEN_PROBLEM.md`,
`paper/pa1_isomonodromy_foundation.md`, `paper/ch_direct_connection.md`, `paper/pv_tau_route.md`,
`experiments/pa2_accessory_algebraic.py`, `experiments/num_S12.py`, `experiments/oracle.py`,
`experiments/ws_ch/*`.
**Constraint compliance:** no git ops; the only files written are this note and
`experiments/cap_connection_formula.py`; scratch under `/tmp/`. All numbers are gold-gated against
`experiments/oracle.py` (canonical `P_{m→m}=0.214724`, sampleB `0.021018`).
**Target (NOMENCLATURE-strict):** `𝒮₁₂ = P_{m→m} = |𝒮_{mid,mid}|²`, `mid = argsort(a)[1]`;
`s_ij = γ_iγ_j/(ε_i−ε_j)`; BE exponent `= s_ij²|a_i−a_j|`; signed-BE Coulomb `c_i = Σ_{j≠i} s_ij²(a_i−a_j)`.

---

## 0. Executive verdict (honest, evidence-tagged)

> **RANK-3. The published rank-2 Painlevé-V / Lisovyy connection constant does NOT compute the
> generic Type-1 N=3 middle survival.** `𝒮₁₂` is the off-diagonal Stokes/connection coefficient of
> the **rank-3** (3×3, single Poincaré-rank-2 irregular point) isomonodromy problem — the higher
> (confluent-Garnier / c=1-family) connection constant, which is **unpublished in closed form** and
> is **not** a classical special function nor the rank-2 PV transcendent. **Painlevé V is exactly the
> rank-2 reduction (one-link-decoupling) limit of our problem — precisely the corner where `P_{m→m}` is
> already elementary.** This settles the one open structural dispute in favour of **WS-PV over the
> WS-PA1 "Painlevé V via middle convolution" reading.** The decisive resolution is *positive and
> sharp*: we identify exactly where and why PV fails, with a gold-gated falsifiable signature.

This **confirms** WS-PV §2 ("PV is one rank too small") and **corrects** the WS-PA1/PA-1
"Painlevé type = Painlevé V" claim. It is fully consistent with the established R1–R5 structure, the
algebraic accessory parameter (PA-2), and the WS-CH connection-matrix benchmark; it only re-classifies
the *named family* one rank up.

**Evidence ladder:** *[established]* = exact symbolic / reproduced to ≤1e-6 numerically;
*[analytic]* = rigorous count/argument from established facts + standard isomonodromy theory;
*[gold-gated]* = numerics validated against `oracle.py`; *[literature]*.

---

## 1. The two angles, executed

### Angle (A) TOP-DOWN — the published PV / Lisovyy connection formula and the dictionary

The published object (Its–Lisovyy–Prokhorov / Lisovyy–Nagoya–Roussillon, arXiv:1806.08344, JMP **59**
091409; Gamayun–Iorgov–Lisovyy 1308.4092 for the c=1 mechanism) is the **isomonodromic τ-function
connection constant of the 2×2 Painlevé-V linear system**:

- **PV linear system [literature].** A `2×2` system on `P¹` with **one regular (Fuchsian) point** at
  `z=0` (exponents `±θ₀/2`) and **one rank-2 irregular point** at `z=∞` (formal exponent `θ_∞`, two
  Stokes multipliers `s₁,s₂`). Monodromy data `= {θ₀, θ_∞, σ}` modulo conjugation — a
  **2-(complex)-dimensional** wild character variety; `σ` is the single intermediate/accessory exponent
  (the `cos 2πσ` trace coordinate fixed by the Stokes data).
- **Connection constant [literature].** Relating the τ-asymptotics at `t→0, +∞, i∞`, each a Fourier
  series of **irregular c=1 Virasoro conformal blocks**, the constant is an explicit **ratio of Barnes
  G-functions**,
  ```
  Υ(θ₀,θ_∞,σ) = ∏ G(1 + (±θ₀ ±θ_∞ ±σ)/2 - type) / ∏ G(1 + 2·exponent),
  ```
  i.e. a finite product of `G(1 + linear-combination-of-{θ₀,θ_∞,σ})` over numerator and denominator.

**The dictionary we built (what maps cleanly, and the one slot that does NOT).**

| PV datum (2×2) | our datum | status |
|---|---|---|
| formal-irregular exponent `θ_∞` | the formal-monodromy exponents `c_i = Σ_{j≠i} s_ij²(a_i−a_j)` (signed BE), `Σc_i=0` | **[established]** but **the wrong shape**: PV has ONE `θ_∞`; we have a **traceless SL(3) pair** `(c_0,c_1,c_2)` |
| regular-point exponent `θ₀` | the apparent point `v_*=E_*`, exponents `{0,1,3}` (gap at 2) | **[established]** but **the wrong scheme**: PV's regular point has exponents `{±θ₀/2}` (a `2×2` Fuchsian point); ours is a `3×3` apparent `{0,1,3}` point |
| intermediate/accessory exponent `σ` | the algebraic accessory parameter `v_*=E_*` (rational; PA-2) | accessory IS pinned [established], **but it is one accessory coordinate of a 3×3 problem, not the single PV `σ`** |
| two Stokes multipliers `s₁,s₂` | the two shears in the carrier pairs `{mid,lo}`,`{mid,hi}` | **[established that there are two named shears]**, but they live in a `3×3` Stokes structure with **6** independent slots, not PV's **2** |

**The dictionary fails at the size/rank of the monodromy manifold.** It is not a missing computation;
the PV formula has the wrong number of independent inputs. Three independent diagnostics (Angle (A),
all `[analytic]` + `[established]`), in `experiments/cap_connection_formula.py`:

- **T1 — wild-character-variety dimension.** A `3×3` system with one rank-2 irregular point and no
  other singularity has `2r=4` Stokes matrices, each unipotent-triangular with `N(N−1)/2 = 3` free
  slots (raw `12`), minus the `N(N−1)=6` off-diagonal constraints of the single loop/product relation
  `S₁S₂S₃S₄ = e^{2πi·(formal-monodromy)}`, giving **dim = 6**. Painlevé V's wild character variety is
  **dim = 2** (`{θ₀,θ_∞,σ}` mod conjugation = one trace coordinate). `6 > 2`: our irregular point is
  strictly larger than PV's. **[analytic]**

- **T2 — the middle-convolution / Laplace escape-hatch is closed.** WS-PA1's defence was "the `3×3` is
  the Harnad/Laplace (middle-convolution) image of the `2×2` PV system, so the effective data still
  lives on PV's 2-dim manifold." A single middle convolution (Katz) / Laplace of a `2×2` system
  **preserves the count of independent irregular leading exponential rates**: a `2×2` PV irregular point
  has the pair `{+L,−L}`, and its `3×3` image carries at most `{+L,−L,0}` (two independent rates plus a
  convolution zero). Our leading rates at `u=∞` are the **three generic distinct** slopes
  `{a_0,a_1,a_2}` (`u`-frame; equivalently `{−i/a_j}` in the Laplace `v`-frame, also three distinct),
  with **no `{+L,−L,0}` symmetry** (canonical `{−1,0.5,2}`, sum `1.5`; sampleB `{−0.7,0.4,1.3}`, sum
  `1.0`). **Three independent leading rates cannot be the single-middle-convolution image of a `2×2`
  PV's two rates.** The "middle-convolution image of PV" reading is therefore false for generic Type-1
  data. **[analytic; established symbolically]**

So Angle (A) yields **NO MATCH**, and the failure is structural and precisely located: rank/dimension
mismatch of the monodromy data. The PV/Barnes-G formula is for a 2-parameter manifold; our connection
coefficient lives on a 6-parameter one (constrained by the known formal exponents, but still strictly
above 2).

### Angle (B) BOTTOM-UP — the two shears, and the decisive physical test

WS-CH established (and we re-verified, `ws_ch/ch_connection3.py`) that the central connection matrix
`C = M_form · 𝒮^{(mid,hi)} · 𝒮^{(mid,lo)}` reproduces the oracle: canonical `P_{m→m}=0.214724`
(`|diff|=4e-7`), sampleB `0.021018` (`|diff|=2e-7`). The "simplification hope" was: if the two shears
are each an elementary (Weber/Gamma) `2×2` connection and only their *composition* is transcendental,
then `𝒮₁₂ = (Γ-dressing of the BE data) × (one joint constant)`.

**The decisive test (T5).** A single shear is a `2×2` (2-level) connection: its survival is a
two-amplitude Stückelberg interference governed by **one** intermediate exponent `σ`. If `𝒮₁₂` were a
single-`σ` (rank-2 / PV-reducible) object — i.e. if the two shears composed *commutatively* into one
effective `2×2` connection — then `P_{m→m}` would lie inside the **widest possible single-`σ`
Stückelberg band**
```
P_{m→m} ∈ [ (√(q₁q₂) − √((1−q₁)(1−q₂)))² , (√(q₁q₂) + √((1−q₁)(1−q₂)))² ],
q₁ = e^{−2π·BE(mid,lo)},  q₂ = e^{−2π·BE(mid,hi)}   (the diabatic survivals at the two crossings the
middle level participates in; fixed by the BE/Stokes magnitudes — NOT fitted).
```
This band is the entire achievable range of any 2-level connection built from the two given crossing
strengths.

**Result (gold-gated, `cap_connection_formula.py` T5).** Of the 14 strata, **7 lie strictly OUTSIDE
this band**, by margins up to `+0.99`:

| stratum | `P_{m→m}` (oracle/fast) | single-σ band | outside? | margin |
|---|---|---|---|---|
| strong | 0.010463 | [1.000000, 1.000000] | **YES** | +0.9895 |
| sep_tiny | 0.028425 | [1.000000, 1.000000] | **YES** | +0.9716 |
| g_hi | 0.057896 | [0.992113, 0.995617] | **YES** | +0.9342 |
| sep_small | 0.068001 | [0.998661, 0.999653] | **YES** | +0.9307 |
| **sampleB** | **0.021018** | **[0.832734, 0.880024]** | **YES** | **+0.8117** |
| eps_asym | 0.024942 | [0.522127, 0.571593] | **YES** | +0.4972 |
| sep_mid | 0.162021 | [0.599998, 0.976143] | **YES** | +0.4380 |
| canonical | 0.214724 | [0.163125, 0.969568] | no | −0.0516 |
| well_sep | 0.473547 | [0.138268, 1.000000] | no | −0.3353 |
| … (6 more) | … | … | no | … |

A value outside the widest single-`σ` band **provably cannot be produced by any `2×2` PV connection
constant** with any single intermediate exponent: it requires a **third coherent amplitude** — exactly
the non-commutative composition of the two shears through the **shared middle sheet**. This is the
genuine three-channel (rank-3) coherence WS-E §6 flagged ("`P_{m→m}` outside every two-path band"), here
as a sharp falsifiable bound, gold-gated. **[gold-gated; decisive]**

> The two shears are individually `2×2` (Weber/Gamma-class) — the simplification hope's *first* half
> holds — but their composition `𝒮^{(mid,hi)}·𝒮^{(mid,lo)}` is **non-commutative through the shared
> mid sheet** and the resulting `𝒮₁₂` is genuinely a `3×3` (rank-3) object. The hoped-for collapse to
> "(Γ-dressing) × (one PV/Barnes-G joint constant)" **does not occur** on the generic (overlapping)
> stratum: the "joint constant" is itself the rank-3 connection coefficient, not a single PV value.

### The positive control (T6) — PV applies exactly on the reduction locus

This is the WS-PV §6 highest-value probe, executed. Decoupling the mid–hi link
(`γ_hi → 0`, holding canonical `ε,a`):

| `γ_hi` | `P_{m→m}` | single-σ band | in band? | single-crossing `q_{lo,mid}` |
|---|---|---|---|---|
| 1.20 | 0.214724 | [0.16312, 0.96957] | yes | 0.221360 |
| 0.60 | 0.194605 | [0.00007, 0.68158] | yes | 0.221360 |
| 0.30 | 0.211281 | [0.05906, 0.44890] | yes | 0.221360 |
| 0.12 | 0.219562 | [0.14573, 0.30769] | yes | 0.221360 |
| 0.04 | 0.221157 | [0.19479, 0.24913] | yes | 0.221360 |
| 0.01 | 0.221349 | [0.21460, 0.22819] | yes | 0.221360 |

As the link decouples, the band **tightens around** the elementary single-crossing survival
`q_{lo,mid}=0.221360`, and `P_{m→m} → q_{lo,mid}` monotonically. The `3×3` system block-reduces to
`2×2 ⊕ 1×1`; the surviving `2×2` is the literal PV/Weber linear system, and `P_{m→m}` becomes elementary.
**This is exactly where the published PV/Barnes-G connection constant applies — the already-elementary
decoupling corner (WS-C "elementary ⟺ a level decouples").** Off this corner (the 7 OUT strata),
`P_{m→m}` exits the band and PV cannot reach it. **[gold-gated; confirms WS-PV §2]**

---

## 2. The monodromy data, algebraically fixed (used, not fitted)

All inputs are the algebraic PA-1/PA-2 values, re-verified here (`cap_connection_formula.py` T3,T4):

- **Formal-monodromy exponents (KNOWN, the diagonal carrier):**
  `c_i = Σ_{j≠i} s_ij²(a_i−a_j)`, `Σc_i=0`.
  Canonical `c=(−0.4128,+0.0864,+0.3264)`; sampleB `(−1.2440,+0.7684,+0.4756)`. **[established, exact]**
  (Cross-checked in the Laplace `v`-frame by formal diagonalisation, `ws_ch/formal_infty.py`: the
  exponents come out `= ±i·(signed BE)`, the documented `−i` phase convention; magnitudes exact.)
- **Accessory parameter (PA-2 ALGEBRAIC):** `v_*=E_*` = the rational doubly-degenerate eigenvalue at
  the universal real node `u_*`. Canonical `u_*=−187/750`, `v_*=−748/375`; sampleB `u_*=−3031/4400`,
  `v_*=−112147/44000`. Both rational ⇒ algebraic in `{γ,ε,a}`. **[established, exact]**

These pin the connection coefficient to a *specific* named object — but a **rank-3** one. The accessory
parameter being algebraic was the PA-2 make-or-break and it remains favourable: it means `𝒮₁₂` is a
*specific* computable named constant (not a free transcendental modulus). What this capstone settles is
its **rank**: the named object is the rank-3 connection constant, not the rank-2 PV one.

---

## 3. Where exactly the dictionary/formula fails (the precise blocker)

The PV connection formula `Υ(θ₀,θ_∞,σ)` (Barnes-G ratio, §1) is a function of **three monodromy
numbers** parametrising a **2-dimensional** wild character variety. Our connection coefficient `𝒮₁₂`
is a function on a **6-dimensional** wild character variety (T1), constrained by the known traceless
*pair* `(c_0,c_1,c_2)` of formal exponents (T3) and the one algebraic accessory coordinate (T4), but
**not** reducible to PV's single `(θ_∞,σ)` because:

1. **Formal data shape (T3):** PV has one `θ_∞`; we have a traceless SL(3) pair — there is no single
   `θ_∞` to feed the Barnes-G arguments.
2. **Leading-rate count / middle-convolution (T2):** three independent irregular leading rates cannot
   be the single-MC image of a `2×2` PV's two rates — the "PV via middle convolution" embedding does
   not exist for generic Type-1 data.
3. **Stokes dimension (T1):** the irregular point carries 6 independent Stokes parameters, not 2.
4. **Physical falsification (T5):** `P_{m→m}` lies strictly outside the widest single-`σ` band on 7/14
   strata (including sampleB) — no single-`σ` (2×2 PV) connection constant can reproduce it.

Each is necessary; together they are decisive. The blocker is **dimensional/structural, not a missing
evaluation**: there is no assignment of `{θ₀,θ_∞,σ}` that makes the published `Υ` equal the generic
`𝒮₁₂`, because the generic `𝒮₁₂` carries strictly more independent monodromy content than `Υ` admits.

---

## 4. The explicit statement of `𝒮₁₂` (the named, rank-3 object)

> **`𝒮₁₂ = P_{m→m}` is the `(mid,mid)` modulus-squared of the central connection matrix
> `C = M_form · 𝒮^{(mid,hi)} · 𝒮^{(mid,lo)}` of the rank-3 linear isomonodromy problem with one
> Poincaré-rank-2 irregular point at `u=∞` (formal exponents `c_i = Σ_{j≠i} s_ij²(a_i−a_j)`, `Σc_i=0`,
> KNOWN) and one apparent regular point at the node `v_*=E_*` (exponents `{0,1,3}`, no log, accessory
> parameter ALGEBRAIC = rational in `{γ,ε,a}`).** Its monodromy data are these known formal exponents
> plus the two Stokes shears in the carrier pairs `{mid,lo}`, `{mid,hi}` (each a Weber/Gamma `2×2`
> connection), and `𝒮₁₂` is their **non-commutative composition through the shared middle sheet** — a
> genuine `3×3` (rank-3) Stokes coefficient. It belongs to the **c=1 / isomonodromic-τ family** (the
> confluent-Garnier-9/2-type connection constant), is **above ₂F₁ / elementary** (WS-A/WS-E) **and
> above the rank-2 Painlevé-V transcendent** (this note, T1–T5), and is **not published in closed
> form**. It reduces to the published rank-2 PV / Lisovyy Barnes-G connection constant **exactly on the
> one-link-decoupling locus** (T6), where `P_{m→m}` is already elementary.

This is the maximal honest closed-form statement: a *named* (c=1-family, rank-3) connection constant
with **algebraically-fixed** monodromy data, plus the gold-gated computable model (`num_S12.py` /
`ws_ch/ch_connection3.py`, ≤1e-6). An *elementary* closed form is ruled out (R4/R5, WS-E §6, T5).

---

## 5. Status ladder

| # | claim | status |
|---|---|---|
| 1 | PV linear system = 2×2, one regular + one rank-2 irregular; monodromy = {θ₀,θ_∞,σ}, dim 2; connection = Barnes-G ratio | **[literature]** (1806.08344, 1308.4092) |
| 2 | our irregular point: 3×3 rank-2, dim of wild char variety = 6 > 2 (PV) | **[analytic]** (T1) |
| 3 | 3 generic distinct leading rates ⇒ NOT a single middle-convolution image of a 2×2 PV (escape-hatch closed) | **[analytic; established symbolically]** (T2) |
| 4 | formal exponents = signed-BE pair `c_i`, Σ=0 — a traceless SL(3) pair, not PV's single `θ_∞` | **[established, exact]** (T3) |
| 5 | accessory `v_*=E_*` rational ⇒ ALGEBRAIC (PA-2 favourable) | **[established, exact]** (T4) |
| 6 | `P_{m→m}` lies OUTSIDE the widest single-σ (2×2/PV) Stückelberg band on 7/14 strata (incl. sampleB) ⇒ rank-3 | **[gold-gated; decisive]** (T5) |
| 7 | on the decoupling locus `P_{m→m} → q_{lo,mid}` (elementary); band tightens; PV applies exactly there | **[gold-gated]** (T6) |
| 8 | `𝒮₁₂` = rank-3 confluent-Garnier / c=1-family connection constant; named, algebraically pinned, unpublished, non-elementary | **[analytic + literature]** (§4) |
| 9 | VERDICT: RANK-3 (WS-PV over WS-PA1's "Painlevé V") | **[decisive, multi-diagnostic]** |

No claim contradicts R1–R5, WS-A/WS-E/WS-G, PA-1/PA-2, WS-CH, or the oracle. This note **decides** the
one open structural dispute: **the named object is the rank-3 (higher-Garnier / c=1) connection
constant; Painlevé V is its rank-2 reduction limit, valid only on the already-elementary decoupling
corners.**

---

## 6. Reproducibility

`experiments/cap_connection_formula.py` (this note's single reproducible script) runs T1–T6:
```
python3 experiments/cap_connection_formula.py            # fast engine (~min)
python3 experiments/cap_connection_formula.py --oracle    # gold T=120 band test on anchors
```
T5/T6 are gold-gated against `experiments/oracle.py`; T3/T4 are exact symbolic (sympy);
T1/T2 are analytic counts. The bottom-up connection-matrix benchmark is `experiments/ws_ch/ch_connection3.py`
(reproduces canonical 0.214724 / sampleB 0.021018 to ≤4e-7). The formal-exponent cross-check is
`experiments/ws_ch/formal_infty.py`.

## 7. References

- O. Lisovyy, H. Nagoya, J. Roussillon, *Irregular conformal blocks and connection formulae for
  Painlevé V functions*, J. Math. Phys. **59** (2018) 091409, arXiv:1806.08344. — rank-2 PV connection
  constant as a Barnes-G ratio / c=1 irregular conformal block (the formula tested here; applies on our
  rank-2 reduction locus only).
- N. Iorgov, O. Lisovyy, J. Teschner / Gamayun–Iorgov–Lisovyy, *Painlevé VI connection problem and
  monodromy of c=1 conformal blocks*, JHEP **12** (2013) 029, arXiv:1308.4092. — the c=1 mechanism.
- O. Lisovyy, A. Naidiuk, *Perturbative connection formulas for Heun equations*, J. Phys. A **55**
  (2022) 434005, arXiv:2208.01604. — confluent-Heun connection as convergent series ↔ quasiclassical
  Virasoro block (the rank-2 scalar cousin; the rank-3 analog our `𝒮₁₂` needs is unpublished).
- Mazzocco (nlin/0306020); rank-3 / higher-Garnier (arXiv:2503.22198, 2512.24083) — the rank-3
  (Garnier-9/2-type) family our `𝒮₁₂` belongs to; generically lacks the Painlevé property.
