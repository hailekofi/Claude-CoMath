# WS-E — Closed form of the 12×13 junction S-matrix / the middle survival `P₂→₂`

**Owner:** WS-E (prize workstream / closed-form prefactor assembly). **Date:** 2026-06-01.
**Inputs read (in order):** `paper/RESEARCH_PROGRAM.md`, `paper/ws_a_riemann_scheme.md`,
`paper/ws_g_stokes_graph.md`, `paper/gate_test_genus.md`, `paper/coscaling_derivation.md`,
`experiments/oracle_report.md`; plus `uploads/assay/{geometry,actions}.py`, `experiments/oracle.py`.
**Literature:** Lin–Sinitsyn (arXiv:1310.7245, the LZ-Coulomb 3-level exact matrix) and
Barik–Bakker–Gritsev–Yuzbashyan "BBGY" (arXiv:2409.17053, KZ ↔ hyperbolic-LZ, the Kampé de Fériet /
₁F₂ solvable 3×3). **Constraint compliance:** no git ops; the only files written are this note and
`experiments/ws_e_*.py`; throwaway probes lived under `/tmp/wse/`.

Model: `H(u)=H0+u·diag(a)`, Type-1 Cauchy coupling, `i dψ/du=H(u)ψ`, `u∈(−∞,∞)`. Observable: the
3×3 doubly-stochastic `P[x,j]`. Slope order `lo,mid,hi=argsort(a)`; the two extreme survivals are
exact BE, the open targets are the **middle survival `P₂→₂≡P[mid,mid]`** and one independent
off-diagonal. Anchors: canonical `eps=(−2,0,3),gam=(1,0.8,1.2),a=(−1,0.5,2) → P_mid=0.214724`;
sampleB `eps=(−1,0,1.5),gam=(0.9,1.1,0.8),a=(−0.7,0.4,1.3) → P_mid=0.021018`.

---

## 0. Executive verdict (honest)

**This is a controlled NEGATIVE-with-structure result, exactly the "clean negative" the brief allows.**

> **The 12×13 junction S-matrix / the middle survival `P₂→₂` of generic Type-1 N=3 MLZ is the
> connection (Stokes) coefficient of a rank-2 irregular point of a genuinely 3-level system whose
> Laplace-transformed scalar ODE carries ONE extra accessory singularity beyond the ₁F₂/Kampé de
> Fériet level. Generic Type-1 is therefore NOT in the solvable (tridiagonal / su(2)-Gaudin) class of
> Lin–Sinitsyn or BBGY — it is precisely the "full 3×3 LZ problem" those papers state is unsolved in
> terms of known special functions. There is no elementary or ₁F₂ closed form for generic `P₂→₂`; the
> honest deliverable is (i) the exact contour-integral representation for our model, (ii) the named
> function class (confluent-Heun / higher-Weber, one accessory parameter above Kampé de Fériet),
> (iii) the exact elementary limits, and (iv) a sharp structural obstruction that explains WHY the
> literature's product/₁F₂ closed forms cannot apply.**

What IS delivered in closed form:
- **The two BE survivals** = elementary, `P[lo,lo]=∏ₘe^{−2πΓ}`, realized as the two imaginary window
  actions `I_X` (verified, §3). [established]
- **The exact integral representation** for all amplitudes (a 3×3 first-order irregular Laplace
  system; the amplitudes are its Stokes data), with the saddle points = the Q₄ turning points / the
  WS-G joint. [analytically-derived, §2]
- **The function class, pinned to one accessory parameter** via the explicit scalar 3rd-order ODE in
  the Laplace variable: a single finite regular-singular point with exponents `{0,1,3}` plus the
  rank-2 irregular point at ∞ — strictly above the ₁F₂/Kampé de Fériet class of BBGY. [established, §4]
- **The precise blocker**: the product integrand does NOT close (the transformed coefficient matrix is
  non-abelian, `[K(v₁),K(v₂)]≠0`), because the linear model lacks the `1/t` Coulomb term whose Laplace
  image is the Fuchsian (rank-1 Bethe-ansatz) structure that makes Lin–Sinitsyn/BBGY solvable. [§2, §5]

What is NOT delivered: a generic elementary / named-special-function value of `P₂→₂`. The benchmark
(§6) shows `P₂→₂` is genuinely independent of the BE survivals and is NOT any of the natural
two-pathway / sech closed forms.

---

## 1. The two solvable templates and why ours is the confluent, non-tridiagonal cousin

**Lin–Sinitsyn (LZC, 1310.7245).** Their solvable 3-level model has one "central" level with energy
`k²/τ` coupled (`g₁,g₂`) to two linear levels `β₁τ,β₂τ` that do NOT couple to each other — a
**bow-tie / Demkov–Osherov** (rank-1-coupling) structure. The contour ansatz `b_j(t)=∫_A e^{−iut}B_j(u)du`
collapses the system to a **first-order** ODE for `B₀(u)`, with **product** solution
`B₀∝(−u)^α(−u+β₁)^{ξ₁}(−u+β₂)^{ξ₂}`, exponents `α,ξ₁,ξ₂` = BE/coupling data (their eq 7–8). The
contour integral of a product is a Beta × ₂F₁ — hence elementary `P` after taking |·|².

**BBGY (2409.17053).** Their solvable 3×3 HLZ (eq 52 / C.11) is
`H=diag(p/t,q/t,0)+offdiag(a₁,a₂)` — **tridiagonal**: the middle level couples to both extremes, the
two extremes do NOT couple. Its KZ/off-shell-Bethe solution is `₁F₂` (a special Kampé de Fériet
`F^{0:1;1}_{1:0;0}`, eq 66, A.27). For equal spacing the middle survival is the elementary
`P₂→₂ = 1 − 1/cosh(π/ν)` (eq 73). BBGY state plainly (p.14, and Conclusion): the **general 3×3 HLZ /
3×3 LZ is "as far as we know, not solvable in terms of known special functions"**, and "all HLZ models
in this work are represented by **tridiagonal matrices**."

**Our Type-1 model is neither.** The Cauchy form gives all three off-diagonals
`(H0)_ij=γᵢγⱼ(aᵢ−aⱼ)/(εᵢ−εⱼ)`; the extreme–extreme coupling `(H0)_{lo,hi}` vanishes only on the
degenerate-slope (trivial-crossing) locus `a_lo=a_hi`. So **Type-1 is never tridiagonal**, and in fact
the "missing in BBGY" coupling is the *dominant* one on the anchors:

| sample | `|H0_(lo,hi)| / max(adjacent)` |
|---|---|
| canonical | **1.20** (extreme–extreme is the LARGEST coupling) |
| sampleB | 0.53 (comparable to adjacent) |

Moreover Type-1 has **no `1/t` Coulomb term** — it is the *linear* model. In the
isomonodromy/confluence picture (WS-A, gate-test §1d) the regular point at `t=0` and the irregular
point at `t=∞` of the hyperbolic model **merge** into the single rank-2 irregular point at `u=∞`. So
Type-1 stands to BBGY as the **confluent, fully-coupled** limit: it inherits the *elementary* part
(the BE/Weber exponents = window actions) but loses the tridiagonal/Coulomb structure that makes the
middle survival elementary.

---

## 2. Exact integral representation for OUR model — and the precise non-closure (analytically-derived)

Apply the same Laplace/Euler ansatz Lin–Sinitsyn and BBGY use, to our linear model:
`ψ_j(u)=∫_C e^{−iuv}B_j(v)dv`. Using `u e^{−iuv}=i ∂_v e^{−iuv}` and integrating by parts (boundary
terms vanish on the contour), `i ψ'=(H0+uA)ψ` transforms to the **first-order 3×3 system**

```
        B'(v) = K(v) B,        K(v) = −i · diag(1/a) · ( H0 − v·I ).
```

(`a_j≠0` required; the degenerate-slope locus is the confluent boundary.) Explicitly
`K_jj = −i(v − H0_jj)/a_j`, `K_ij = −i H0_ij/a_i` (i≠j), i.e. **linear in `v`** with a diagonal
`−iv/a_j` term. This is the exact analog of Lin–Sinitsyn eq (7): the amplitudes are contour integrals
whose saddle points `v_*(u)` solve `det(H0 − v I + i a v'(...))`-type conditions and coincide with the
Q₄ turning points (the WS-G complex avoided crossings / the joint). **[analytically-derived]**

**Why it does NOT close to a product integrand (the precise blocker).** For Lin–Sinitsyn/BBGY the
`1/t` term makes the transformed *scalar* equation first-order, so `B₀` is a single power-product. For
our linear model:

1. The transformed system `B'=K(v)B` is **still a 3×3 first-order system, irregular of the same
   Poincaré rank at `v=∞`** — a Laplace transform of a linear (no-`1/t`) LZ system stays irregular.
2. `K(v)` is **non-abelian along the flow**: `[K(v₁),K(v₂)]≠0` (computed symbolically, nonzero;
   reproduced by `experiments/ws_e_laplace_class.py`).
   Hence there is no `v`-independent frame diagonalizing it, and the integrand is **not** a product
   `∏(v−b_k)^{ρ_k}`. The rank-1 (Bethe-ansatz `∏L̂⁺(λ_α)`) product structure that BBGY's su(2)-Gaudin
   model enjoys is absent.

This is exactly the caution flagged in the brief and the deep-verification review: *"a naive Laplace
transform of the linear model stays irregular; the N=3 transformed system is 3rd-order."* Confirmed,
and we now know the mechanism: **no Coulomb term ⇒ no Fuchsian image ⇒ no product integrand.**

---

## 3. The BE survivals ARE the two window actions (elementary part, established)

The two imaginary window actions `I_X` from `uploads/assay/actions.py` reproduce the two **extreme**
BE survival exponents exactly (each is a SUM over the two crossings that extreme level undergoes):

| sample | `|Im I_X|` (two windows) | matched to `2πΓ`-sums |
|---|---|---|
| canonical | `1.508`, `2.051` | `2πΓ_{01}=1.508`; `2πΓ_{02}+2πΓ_{12}=1.086+0.965=2.051` |
| sampleB | `2.989`, `7.816` | sums of the appropriate `2πΓ` pairs |
| well_sep | `0.566`, `0.566` | symmetric, `2πΓ`-sums |

So the windows deliver `P[lo,lo]`, `P[hi,hi]` (the elementary BE/Weber part), confirming WS-A §2 and
E1. **The middle survival is NOT a window action** — it is the genuine non-Abelian leftover. **[established]**

---

## 4. Function class, pinned by the scalar Laplace ODE (established)

Eliminating `B'=K(v)B` to a single 3rd-order scalar ODE for `B₀(v)` (cyclic-vector reduction, exact
sympy) gives `y''' = c₂(v)y'' + c₁(v)y' + c₀(v)y` with **a single finite singular point** at the root
of the common denominator (canonical: `v* = −748/375 = −1.9947`; sampleB: `v* = −112147/44000`), plus
the rank-2 irregular point at `v=∞`. The finite point is a genuine accessory point: it is NOT a
diabatic energy nor an eigenvalue of `H0` (canonical `eig H0 = {0,−2.04,−1.69}`), and its **indicial
exponents are `{0,1,3}`** — integer but with a GAP at 2.

Compare BBGY eq 64: their scalar ODE has singular points only at `0` and `∞` and reduces to `₁F₂`
(Kampé de Fériet `F^{0:1;1}_{1:0;0}`). **Our equation carries ONE extra finite regular-singular point**
with the anomalous `{0,1,3}` indicial set ⇒ one extra accessory parameter ⇒ it lies **strictly above
the ₁F₂ / Kampé de Fériet class**, in the **genus-0 confluent-Heun / higher-Weber** class — exactly
WS-A's verdict, now confirmed independently in the Laplace frame. **[established]**

This is the rigorous statement of the obstruction to a BBGY-style ₁F₂ closed form: the middle
survival is a connection coefficient of a confluent equation with an accessory parameter that ₁F₂ does
not have.

---

## 5. Why the node / cross-ratio does not rescue closure (numerically-supported)

WS-G identifies the transcendental content with the 12×13 joint holonomy, conjecturally a function of
the two window actions plus the node / the cross-ratio of the four complex turning points. We tested
this directly (§6, `experiments/ws_e_geometry.py`). The four complex turning points form two
conjugate pairs, so their **cross-ratio is real** — a single real geometric invariant per sample
(canonical `≈ 64.5`, sampleB `≈ 0.195`, well-sep `≈ 0.933`). The node sits on the real axis between
the joint clusters (WS-A/WS-G), and as an ORDINARY point of the `u`-system it contributes **no** local
connection data (WS-A §6) — confirmed: it cannot pin `P₂→₂`.

Crucially, even *with* the two window actions and the cross-ratio in hand, there is no elementary map
to `P₂→₂`: the value is the **connection coefficient** of the §4 confluent-Heun-class equation, whose
connection coefficients are generically not known in closed form (a standard fact for confluent Heun).
The geometric invariants fix the *arguments* of the transcendental object, not its value.

---

## 6. Benchmark: `P₂→₂` is genuinely transcendental (numerically-supported)

(Oracle = `experiments/oracle.py`, T-converged; T=40 and T=120 agree to ≤1e-5 on the anchors.)

### 6a. `P₂→₂` is independent of the BE survivals + double stochasticity
A 3×3 doubly-stochastic matrix has 4 dof; fixing the two BE survivals leaves a **2-dimensional**
family along which `P[mid,mid]` varies (shifts of −0.211 and +0.789 per unit; `experiments/ws_e_geometry.py`).
So `P₂→₂` is a genuine independent dynamical quantity — no algebraic shortcut. **[established]**

### 6b. Natural elementary / two-pathway closed forms FAIL
With `δ₁=2πΓ_{mid,lo}`, `δ₂=2πΓ_{mid,hi}`, `pᵢ=e^{−δᵢ}` (diabatic stay probs):

| candidate | canonical (true 0.21472) | sampleB (true 0.02102) |
|---|---|---|
| incoherent `p₁p₂` | 0.08432 (the WS-G baseline) | 0.000163 |
| `1−sech((δ₁+δ₂)/2)` (BBGY-shape) | 0.4644 ✗ | 0.9745 ✗ |
| `p₁p₂+(1−p₁)(1−p₂)` (max-coherence loop) | 0.5663 ✗ | 0.8564 ✗ |
| 2-path Stückelberg band `[lo,hi]` | `[0.163,0.970]` ∋ true (cosφ=−0.872) | `[0.833,0.880]` **∌** true ✗ |

The decisive line is the last: for **sampleB the true `P₂→₂` lies OUTSIDE** the entire two-pathway
interference band — no choice of a single Stückelberg phase between the two BE crossings can produce
it. The middle survival is a genuine **three-crossing (joint) coherence**, not a product or a
two-path form. This is the numerical signature of the §4 extra accessory parameter. **[established]**

### 6c. The independent off-diagonal carries the SAME transcendence
The brief also asks for one independent off-diagonal. By double stochasticity + the two exact BE
survivals, the middle row is `(P[mid,lo], P[mid,mid], P[mid,hi])` with `P[mid,lo]+P[mid,hi]=1−P[mid,mid]`;
the split between the two off-diagonals is the second independent dof (§6a). Both are entries of the
SAME rank-2 irregular connection matrix and inherit the §4 confluent-Heun-class transcendence (they are
the off-diagonal Stokes coefficients linking the middle sheet to lo/hi). They are NOT given by the
Demkov–Osherov / incoherent product either: e.g. for sampleB the converged off-diagonals
`P[mid,lo]=0.9761, P[mid,hi]=0.00287` are nowhere near the incoherent semiclassical split. So the
"one independent off-diagonal" is the partner connection coefficient of the same junction, equally
non-elementary. [established]

### 6d. Band-containment statistic (the sharp falsification)
7-sample scan (oracle T=40, anchors verified T=120). `δᵢ=2πΓ` for the two crossings the middle level
sees; `pᵢ=e^{−δᵢ}`; the two-path Stückelberg band is `[(√p₁p₂−√q₁q₂)², (√p₁p₂+√q₁q₂)²]`, `qᵢ=1−pᵢ` —
the FULL envelope of every two-amplitude (stay-stay vs flip-return) interference model. `P_mid` lying
outside it falsifies ALL such forms (no choice of a single Stückelberg phase fits).

| sample | ratio | `P_mid` (oracle) | `inc=p₁p₂` | `1−sech` | 2-path band | in band? |
|---|---|---|---|---|---|---|
| canonical | 0.28 | 0.21472 | 0.08433 | 0.4644 | [0.163, 0.970] | yes |
| sampleB | 0.14 | **0.02102** | 0.00016 | 0.9745 | [0.833, 0.880] | **NO** |
| well_sep | 0.00 | 0.47355 | 0.47049 | 0.0671 | [0.138, 1.000] | yes (≈incoherent) |
| weak | 0.20 | 0.95988 | 0.95977 | 0.0002 | [0.922, 0.999] | yes (≈incoherent) |
| near_deg | 0.35 | 0.28518 | 0.21827 | 0.2330 | [0.062, 0.470] | yes |
| ov1 | 0.15 | **0.11624** | 0.01022 | 0.7998 | [0.625, 0.986] | **NO** |
| ov2 | 0.04 | **0.09731** | 0.05072 | 0.5713 | [0.244, 0.891] | **NO** |

(The `1−sech` column uses the placeholder argument `(δ₁+δ₂)/2`; BBGY's exact `1−sech(π/ν)` is in a
different parametrization, but the point is only that no simple sech of the BE data fits — it is off
by O(1) everywhere.)

**`P_mid` is OUTSIDE the two-path band on 3 of 7 samples — all three strongly overlapping
(sampleB, ov1, ov2).** No two-amplitude / product / sech closed form can reproduce these. The
well-separated/weak samples sit at the incoherent edge (`P_mid≈p₁p₂`, the factorizing locus), exactly
as WS-G predicts. Conclusion: the overlapping-regime middle survival is a genuine three-crossing
(joint) coherence — the §4 confluent-Heun connection coefficient, not an elementary form. **[established]**

---

## 7. Status ladder

| # | claim | status |
|---|---|---|
| 1 | BE survivals = two imaginary window actions `I_X` (elementary part) | **established** |
| 2 | Exact Laplace integral rep `ψ_j=∫e^{−iuv}B_j dv`, `B'=K(v)B`; saddles = Q₄ turning pts/joint | **analytically-derived** |
| 3 | `[K(v₁),K(v₂)]≠0` ⇒ no product integrand ⇒ no Bethe-ansatz closure (no `1/t` ⇒ no Fuchsian image) | **established** (symbolic) |
| 4 | Scalar Laplace ODE: 1 finite singular pt, exponents `{0,1,3}`, + rank-2 irregular at ∞ | **established** (symbolic) |
| 5 | ⇒ class is confluent-Heun/higher-Weber, ONE accessory param ABOVE ₁F₂/Kampé de Fériet | **analytically-derived** |
| 6 | Type-1 is never tridiagonal (extreme–extreme coupling ≠0; dominant on canonical) | **established** |
| 7 | `P₂→₂` independent of BE + double-stochasticity (2-dof family) | **established** |
| 8 | `P₂→₂` ≠ any elementary 2-pathway / sech form (sampleB outside the interference band) | **established** (numeric) |
| 9 | generic `P₂→₂` = confluent-Heun connection coeff, not elementary / not ₁F₂ | **numerically-supported + analytically-framed** |

No claim contradicts E1–E5 or WS-A/D/G; this work **sharpens** WS-A (the extra accessory point, now
also in the Laplace frame) and **explains** WS-G/co-scaling (the three-crossing joint coherence is why
no product/two-path form works) and confirms WS-D's obstruction from the connection side.

---

## 8. The precise blocker, stated for the coordinator

The middle survival is the off-diagonal Stokes/connection coefficient of the rank-2 irregular point of
a 3×3 system whose Laplace image is a confluent-Heun-class scalar ODE with an accessory parameter
(`{0,1,3}` finite point). Closing it elementarily would require either (a) a `1/t` Coulomb term (to
make the Laplace image Fuchsian, giving a product/Bethe-ansatz integrand) — absent in the linear
model; or (b) a tridiagonal structure (su(2) Gaudin, BBGY) — absent because Type-1's extreme–extreme
Cauchy coupling never vanishes. Both are structural, not removable by a change of frame. Therefore the
clean result is: **`P₂→₂` is exactly a genus-0 confluent-Heun connection coefficient (one accessory
parameter above the ₁F₂/Kampé de Fériet of the solvable tridiagonal cousins), benchmarked, not
reducible to elementary functions or to ₁F₂ on the generic Type-1 locus.** A closed *value* would
require the (generically unknown) confluent-Heun connection coefficients; the exact integral
representation of §2 + the arguments fixed by the window actions/cross-ratio (§3,§5) is the honest
maximal result.

---

## 9. Reproducibility

- `experiments/ws_e_laplace_class.py` — symbolic: the transformed system `B'=K(v)B`, the commutator
  `[K(v₁),K(v₂)]≠0`, the scalar 3rd-order ODE in `v`, its single finite singular point and `{0,1,3}`
  indicial exponents (both anchors).
- `experiments/ws_e_geometry.py` — turning points (exact discriminant), cross-ratio, window actions =
  BE exponents, the double-stochastic dof count, and the candidate-form benchmark vs the oracle.
- Oracle: `experiments/oracle.py` (imported). Anchors reproduced: canonical `P_mid=0.214724`,
  sampleB `P_mid=0.021018` (matches `oracle_report.md`).
