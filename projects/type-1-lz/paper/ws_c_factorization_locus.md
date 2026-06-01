# WS-C — The factorization locus of the Type-1 N=3 MLZ transition matrix (gate-test #2 / P-ii)

**Owner:** WS-C (factorization / reducibility of the connection). **Date:** 2026-06-01.
**Inputs used:** WS-F gold oracle (`experiments/oracle.py`), WS-A Riemann scheme
(`paper/ws_a_riemann_scheme.md`), WS-G Stokes graph (`paper/ws_g_stokes_graph.md`),
co-scaling theorem (`paper/coscaling_derivation.md`), WS-D obstruction
(`experiments/ws_d_verification.py`). **Constraint compliance:** no git ops; the only files
written are this note and `experiments/ws_c_*.py`; throwaway probes lived under `/tmp/`.

Model: `H(u)=H0+u·diag(a)`, Type-1 Cauchy coupling
`V_ij=(H0)_ij=γ_iγ_j(a_i−a_j)/(ε_i−ε_j)`, `(H0)_ii=−Σ_{k≠i}γ_k²(a_i−a_k)/(ε_i−ε_k)`.
Observable: the `3×3` doubly-stochastic `P[x,j]=prob(x→j)`. The open quantity is the
**middle-slope survival** `P_mid≡P[mid,mid]`, `mid=argsort(a)[1]`. "Factorizes" ≡ the exact
`P_mid` equals the **incoherent product of the two two-level LZ factors**
`P_mid^{inc}=exp(−2π(Γ_{mid,lo}+Γ_{mid,hi}))`, `Γ_ij=γ_i²γ_j²|a_i−a_j|/(ε_i−ε_j)²`.

---

## 0. Executive summary / VERDICT on P-ii

**P-ii is CONFIRMED in its sharp form, with one precise correction.** The factorization locus
is exactly the **trivial-middle-coupling boundary**: `P_mid` reduces to the elementary
incoherent product `e^{−2π(Γ_{mid,lo}+Γ_{mid,hi})}` **iff the middle level effectively
decouples from (at least) one outer level over the passage**, i.e. one of its two LZ events
becomes trivial. There is **NO genuine codimension-1 elementary locus in the interior**
(all `γ≠0`, all `a` distinct, all `ε` distinct, all `Γ_ij` of order 1): on the genuine-LZ
locus the matrix is *always* non-factorizing.

Concretely, three things coincide (numerically-supported, and analytically explained):

- **factorizable** (`|P_mid − P_mid^{inc}| → 0`)
- **⟺ joint-free** (WS-G: no `12×13` Stokes joint)
- **⟺ trivial-middle-coupling** (the off-diagonal Cauchy coupling joining the middle level to
  one outer level → 0, equivalently `Γ_{mid,outer} → 0`).

The single sharpest number: across a random Type-1 scan, the factorization defect is
correlated with the middle coupling (`Pearson log₁₀ absdef vs log₁₀ Γ_mid,min = +0.51`) and
**uncorrelated with the separation ratio** (`vs log₁₀(sep/width) = +0.015`, §5). A controlled
weak-link sweep (§3) collapses the defect smoothly `0.13 → 3.7×10⁻⁴` as `γ_mid→0` with no
special surface crossed. Symbolically (§2) the only loci that kill an off-diagonal coupling are
`γ_i=0` (decoupling) or `a_i=a_j` (degenerate slope, where `Γ_ij→0` too); the connection is
irreducible everywhere else.

**The one correction to the naive P-ii.** Factorization is **NOT** controlled by the
`sep/width` separation ratio, and it is **NOT** controlled by `min(Γ_{mid,lo},Γ_{mid,hi})`
alone either. A single small middle coupling is *necessary but not sufficient*: if the
**other** middle coupling is `O(1)` **and** its crossing overlaps the third crossing, the
middle level still mixes and `P_mid` is enhanced (the `near_deg_slope` stratum, §1, is the
witness: `Γ_{mid,min}=0.010` yet enh `=1.31`). The clean classifier is "the middle level
decouples" — i.e. the *product* `Γ_{mid,lo}·Γ_{mid,hi}` together with the requirement that the
**non-trivial** middle crossing not overlap a third. Stated as the established structural
result: **elementary ⟺ a level decouples** (one off-diagonal Cauchy coupling is negligible
across the whole passage); the genuine-LZ interior is *entirely* transcendental.

This also reconciles WS-G's sampleA: it factorizes (`enh=1.00`) **not** because it is
"well-separated" (its `sep/width` ratio is 3.97) but because its middle coupling is trivial
(`Γ_{mid,lo}=1.2×10⁻⁵`, §1). The "separation" was a red herring; the joint-free-ness WS-G
saw there is a *consequence* of the vanishing middle coupling (no Stokes line of one pair has
any amplitude to joint the other), not of geometric separation.

**Evidence-ladder tags** are attached to every claim below.

---

## 1. Enhancement vs trivial-middle-coupling, on the named strata — [numerically-supported]

Gold oracle (`experiments/ws_c_factorization_scan.py`; `enh=P_mid/P_mid^{inc}`,
`absdef=|P_mid−P_mid^{inc}|`; `Gmid_lo/hi=Γ_{mid,lo/hi}`; `ratio=sep/width` from co-scaling).
The two BE extreme survivals are *always* exactly the incoherent product (E1) — only the
**middle** survival is at issue.

| stratum | enh (×) | absdef | Γ_mid,lo | Γ_mid,hi | Γ_mid,min | sep/width |
|---|---:|---:|---:|---:|---:|---:|
| canonical        |   2.546 | 1.3e−1 | 0.2400 | 0.1536 | 0.1536 | 0.282 |
| sampleB_overlap  | 128.747 | 2.1e−2 | 1.0781 | 0.3098 | 0.3098 | 0.142 |
| **wsg_sampleA**  |   1.000 | 6.9e−5 | **1.2e−5** | 0.0023 | **1.2e−5** | 3.970 |
| well_separated   |   1.006 | 3.1e−3 | 0.0600 | 0.0600 | 0.0600 | 0.000 |
| weak_coupling    |   1.000 | 1.0e−4 | 0.0041 | 0.0024 | 0.0024 | 0.196 |
| near_deg_slope   | **1.307** | 6.7e−2 | 0.2320 | **0.0102** | 0.0102 | 0.349 |

**Readings.**
1. **The two clean factorizing strata** (`wsg_sampleA`, `weak_coupling`) both have
   `Γ_mid,min ≲ 2×10⁻³` — a (near-)decoupled middle level. `well_separated` is at `0.06`
   (all `Γ` small, weak-coupling regime) and factorizes to `3×10⁻³`.
2. **`wsg_sampleA` factorizes by trivial coupling, not by separation.** Its `sep/width=3.97`
   (the "well-separated" tail), but the decisive quantity is `Γ_{mid,lo}=1.2×10⁻⁵`: the middle
   level is **decoupled** from the lo level. *Separation is neither necessary nor sufficient.*
3. **`well_separated` has `sep/width=0` yet factorizes** — two of its crossings *coincide*
   (ratio 0), the opposite of "separated", and it still factorizes because all `Γ` are small.
   Again: separation is the wrong variable; coupling strength is the right one.
4. **`near_deg_slope` is the sharp counterexample to "min middle coupling ⇒ factorizes".**
   `Γ_{mid,min}=0.010` is small, but `Γ_{mid,lo}=0.232` is `O(1)` **and** the `mid–hi`
   crossing (the weak one) sits *on top of* the `mid–lo` crossing (degenerate slopes
   `a_mid≈a_hi`). The middle level still mixes strongly through its `O(1)` `lo` coupling, so
   `P_mid` is enhanced 1.31× — **not factorizing**. ⇒ the correct condition is "the middle
   level **decouples**", i.e. its *non-trivial* coupling must not feed an overlapping crossing,
   not merely "some middle `Γ` is small."

---

## 2. Symbolic factorization test (gate-test #2 proper) — [analytically-derived (exact)]

`experiments/ws_c_symbolic_factor.py` (sympy, over `Q(γ,ε,a)`):

- **The off-diagonal Cauchy coupling vanishes on an algebraic locus iff a level decouples:**
  `V_ij = γ_iγ_j(a_i−a_j)/(ε_i−ε_j) = 0 ⟺ γ_i=0 ∨ γ_j=0 ∨ a_i=a_j`.
  There is **no interior locus** (all `γ≠0`, all `a` distinct, all `ε` distinct) on which a
  single coupling vanishes.
- **For the middle level to fully decouple** (BOTH `V_{mid,lo}=V_{mid,hi}=0`) with all `γ≠0`
  requires `a_lo=a_mid=a_hi` — all slopes equal — whereupon `H(u)=H0+u·a·I` is a **constant
  matrix up to a scalar phase**: there is no Landau–Zener problem at all (no crossings). So the
  *only* interior, all-slopes-distinct way to decouple the middle level is `γ_mid→0`
  (a trivial, degenerate limit), and the *only* single-coupling way is the
  degenerate-slope boundary `a_mid=a_outer` (where simultaneously `Γ_{mid,outer}→0`).
- **`H0` does not block-diagonalize, and its characteristic polynomial does not factor over
  `Q(γ,ε,a)`** (beyond the trivial `λ` factor coming from `tr`-bookkeeping); the genuine
  quadratic factor is irreducible. **The connection is irreducible on every interior locus.**

**Gate-2 verdict: FAIL in the constructive sense (no non-trivial interior elementary locus).**
The `3×3` `λ`-connection is irreducible everywhere except on the trivial-coupling boundary
(`γ_mid→0`) and the degenerate-slope boundary (`a_i=a_j`, where the runaway crossing is
diabatic and `Γ→0`). This is the **structural result**: confluent-Heun is *genuine* in the
Type-1 interior; the elementary tail is exactly the decoupling boundary. (It also *confirms*
WS-G's geometric criterion — joint-free ⟺ a coupling that could joint is absent ⟺ a level
decouples — from the algebra side.)

---

## 3. Controlled weak-link sweep (the decisive test) — [numerically-supported]

`experiments/ws_c_random_scan.py A`. Base point = canonical (`ε=(−2,0,3)`, `a=(−1,0.5,2)`,
`γ=(1,0.8,1.2)`; `mid=1`). Scale **only the middle coupling** `γ_mid→s·γ_mid`
(so `Γ_{mid,lo},Γ_{mid,hi}∝s²`), holding everything else fixed, and watch the factorization
defect collapse. The fast (`rtol=10⁻¹⁰,T=40`) and gold (`rtol=10⁻¹³,T=50`) settings agree to
all printed digits on the points computed at both.

| s | Γ_mid,min | Γ_mid,lo·Γ_mid,hi | absdef | enh (×) | sep/width |
|---:|---:|---:|---:|---:|---:|
| 1.00 | 1.536e−1 | 3.69e−2 | 1.304e−1 | 2.546 | 0.282 |
| 0.70 | 7.526e−2 | 8.85e−3 | 1.692e−1 | 1.568 | 0.645 |
| 0.50 | 3.840e−2 | 2.30e−3 | 1.386e−1 | 1.257 | 0.886 |
| 0.30 | 1.382e−2 | 2.99e−4 | 6.864e−2 | 1.086 | 0.975 |
| 0.20 | 6.144e−3 | 5.90e−5 | 3.372e−2 | 1.037 | 1.003 |
| 0.10 | 1.536e−3 | 3.69e−6 | 8.952e−3 | 1.009 | 1.019 |
| 0.05 | 3.840e−4 | 2.30e−7 | 2.272e−3 | 1.002 | 1.024 |
| 0.02 | 6.144e−5 | 5.90e−9 | **3.65e−4** | **1.000** | 1.025 |

**The defect `absdef→0` monotonically (past the first hump) as `Γ_mid→0`**, crossing the
`absdef<10⁻³` factorization threshold near `s≈0.04` (`Γ_mid,min≈3×10⁻⁴`). At `s=0.02`
(`Γ_mid,min=6×10⁻⁵`) the incoherent product is exact to `3.6×10⁻⁴` and `enh=1.000`. This is
the clean controlled realization of the trivial-coupling boundary: **the middle level
decoupling drives `P_mid` smoothly to the elementary incoherent product**, with no special
codimension-1 surface crossed along the way — the limit is reached only *at* `Γ_mid→0`.

(Note the non-monotone bump at `s=0.7`: weakening one coupling first slightly *raises* the
relative interference before the eventual collapse. The monotone variable is `Γ_mid,min`, and
the asymptotic `absdef ∝ Γ_mid` law is clear from the last four rows: `absdef/Γ_mid,min ≈
5.8, 5.5, 5.9, 5.9` — a clean linear vanishing.)

---

## 4. The BBGY symmetric slice `ε_2=(ε_1+ε_3)/2` is NOT a factorizing locus — [num. + symbolic]

BBGY-class hyperbolic 3-state models reduce their `3×3` Kampé de Fériet to a `₁F₂` on the
symmetric line (middle pole at the midpoint). We tested whether the **linear** Type-1 model
inherits an elementary middle survival there. `experiments/ws_c_bbgy_locus.py`: for two base
`(γ,a)` sets, put `ε=(−2, ε_2, 3)` with `ε_2=(−2+3)/2=0.5` (symmetric) vs off-symmetric
controls (`ε_2=−0.4, 1.4`), and measure `absdef`.

| case | ε_2 | P_mid | P_mid^{inc} | absdef | sep/width |
|---|---:|---:|---:|---:|---:|
| gA **SYMMETRIC** | 0.500 | 0.2300 | 0.0949 | **1.35e−1** | 0.302 |
| gA off | −0.400 | 0.1561 | 0.0447 | 1.11e−1 | 0.266 |
| gA off | 1.400 | 0.1043 | 0.0199 | 8.43e−2 | 0.343 |
| gB **SYMMETRIC** | 0.500 | 0.1970 | 0.1273 | **6.97e−2** | 0.141 |
| gB off | −0.400 | 0.1026 | 0.0418 | 6.08e−2 | 0.169 |
| gB off | 1.400 | 0.1132 | 0.0512 | 6.19e−2 | 0.154 |

**The symmetric slice is fully non-factorizing**, with `absdef` on the symmetric line *no
smaller* (indeed slightly larger) than at nearby off-symmetric points. **The BBGY `₁F₂`
reduction does NOT transfer to the linear Type-1 model.** This is consistent with the
symbolic finding (§2 method, run on the slice): on `ε_2=(ε_1+ε_3)/2` the couplings become
`V_{01}=2γ_0γ_1(a_0−a_1)/(ε_1−ε_3)`, `V_{12}=2γ_1γ_2(a_1−a_2)/(ε_1−ε_3)`,
`V_{02}=γ_0γ_2(a_0−a_2)/(ε_1−ε_3)` — **none vanishes**, and the characteristic polynomial's
quadratic factor remains **irreducible** over `Q`. The symmetric pole position is a
convenient parametrization, not an algebraic reduction, for the linear model. **BBGY's
symmetric elementary locus is a feature of the hyperbolic model's special-function structure,
not a universal MLZ factorization locus.**

---

## 5. Random scan: the defect tracks the middle coupling, NOT the separation — [numerically-supported]

`experiments/ws_c_random_scan.py B` (14 random Type-1 samples, gold oracle, `rtol=10⁻¹⁰`).
Recorded `absdef`, `enh`, `Γ_mid,min`, `Γ_mid,lo·Γ_mid,hi`, `sep/width` per sample.

**The headline correlations (Pearson vs `log₁₀ absdef`):**
```
  log10(absdef) vs log10(Γ_mid,min) = +0.514     (defect grows with the middle coupling)
  log10(absdef) vs log10(sep/width) = +0.015     (≈ ZERO: separation is irrelevant)
```
This is the quantitative core of P-ii: **the factorization defect is correlated with the
middle coupling and essentially uncorrelated with the separation ratio.** The samples with
`absdef<10⁻³` (factorizing) cluster at small `Γ_mid`; e.g. sample 1 (`Γ_mid,min=2.1×10⁻³`,
`absdef=5.6×10⁻⁴`, `sep/width=1.05`) and sample 8 (`Γ_mid,min=1.5×10⁻³`,
`Γ_mid,prod=8.9×10⁻⁶`, `absdef=8.0×10⁻⁴`). The strongly non-factorizing samples
(7: `absdef=0.215`; 14: `absdef=0.165`) have `Γ_mid` of order `0.04–0.09` and overlapping
crossings — and crucially, samples with **large** `sep/width` (2: ratio 1.10; 11: ratio 0.78;
5: ratio 1.29) are **non-factorizing** because their middle coupling is `O(10⁻²)`. Separation
does not buy factorization.

**Honest caveat (the small-`P_mid` subtlety).** `absdef=|P_mid−P_mid^{inc}|` can be small not
only when the middle level decouples but also when `P_mid` is *intrinsically tiny* (a strongly
adiabatic middle level): samples 13 (`Γ_mid,min=0.060`, `absdef=5.2×10⁻⁴`, **`enh=5×10⁶`**) and
3 (`absdef=5.9×10⁻⁴`, `enh=1.56`) have small *absolute* defect with `O(1)` middle coupling and
a large *ratio* `enh`, because `P_mid^{inc}` underflowed. These are **not** physical
factorizations — they are the deep-adiabatic corner where both `P_mid` and `P_mid^{inc}` are
`≈0`. The clean factorization criterion is therefore the **conjunction** `enh→1` **and**
`absdef→0` (incoherent product exact in *both* absolute and relative terms), which is exactly
what the trivial-coupling limit delivers (the §3 sweep: `enh→1.000`, `absdef→3.7×10⁻⁴`
*together*). With that conjunction the picture is clean; the `Γ_mid`-ranges overlap mildly
because of these deep-adiabatic and overlap edge cases, which is why the *operational/geometric*
("middle level decouples / joint-free", §6) statement — not a single scalar threshold — is the
precise classifier.

---

## 6. Reconciliation with co-scaling (genuine-LZ ⇒ always non-factorizing) — [analytic + numeric]

Co-scaling (`paper/coscaling_derivation.md`) proved that on the **genuine-LZ locus** — all
`Γ_ij∈[Γ_min,Γ_max]⊂(0,∞)` — the crossings are *permanently marginally-overlapping*
(`sep/width ≲ 1.6`). WS-C closes the loop:

- **Genuine-LZ ⇒ non-factorizing, always.** If all three `Γ_ij` are `O(1)` then in particular
  both middle couplings `Γ_{mid,lo},Γ_{mid,hi}` are `O(1)`, the middle level is *coupled* to
  both outer levels, and (by co-scaling) the relevant crossings overlap. Both the algebra (§2:
  the connection is irreducible) and the geometry (WS-G: a `12×13` joint must form) then force a
  non-trivial off-diagonal Stokes coefficient. **The genuine-LZ interior is entirely
  jointed/transcendental** — consistent with the established ~85% enhanced fraction.
- **The factorizing tail is exactly the trivial-coupling boundary.** The only way to leave the
  genuine-LZ band is to send some `Γ_ij→0`. If it is a *middle* coupling, the middle level
  (partially) decouples and `P_mid` heads to the incoherent product. If it is the *outer–outer*
  coupling, the middle level still couples to both and `P_mid` need **not** factorize unless one
  of *its* couplings is also small — which is why `min(Γ_mid)` small is necessary but the
  *decoupling* (the non-trivial middle crossing not overlapping a third) is what is sufficient
  (§1, `near_deg_slope`). Either way the factorizing set is a **boundary** (codimension ≥ 1 and
  attained only as `Γ_{mid,outer}→0`), **not** an interior codimension-1 elementary surface.

**Net:** `factorizable = joint-free = trivial-middle-coupling`, and all three are the
**boundary of the genuine-LZ region** reached as a middle coupling → 0 — *not* a generic
separation effect. The bimodal classifier "elementary ⟺ a level decouples" is sharp.

---

## 7. Honest residuals

1. **"Decouples" needs the overlap qualifier.** The cleanest one-line classifier is
   `Γ_{mid,lo}·Γ_{mid,hi} → 0`, but `near_deg_slope` shows a small *product* is not by itself
   sufficient when the surviving non-trivial middle crossing overlaps the third crossing
   (degenerate slopes). The fully precise statement is geometric (WS-G's joint-free criterion)
   or operational (`absdef→0`); the `Γ`-product is a good *scalar proxy* away from the
   degenerate-slope corner. Flagged.
2. **Boundary, not interior surface.** "Trivial-coupling boundary" means the factorizing set
   has empty interior in the genuine-LZ region; it is reached only in the limit
   `Γ_{mid,outer}→0`. There is no finite-codimension elementary *slice through the interior*.
   This is the precise (negative) answer to O3.
3. **Threshold dependence.** "Factorizes" is operationalized as `absdef<10⁻³` (the oracle is
   gold to ≲10⁻⁷, so this is safe). The strata and sweep show a clean monotone collapse, so the
   verdict is threshold-robust; the precise locus is the limit, not a level set.

---

## 8. Reproducibility

- `experiments/ws_c_symbolic_factor.py` — symbolic (sympy) irreducibility / coupling-vanishing
  test over `Q(γ,ε,a)`; the §2 results.
- `experiments/ws_c_factorization_scan.py` — named-strata enhancement-vs-trivial-coupling table
  (§1), using the WS-F gold oracle.
- `experiments/ws_c_random_scan.py` — the controlled weak-link sweep (§3) and the random scan
  (§5).
- `experiments/ws_c_bbgy_locus.py` — the BBGY symmetric-slice test (§4).
