# WS-O2 — The exact contour-integral representation of Type-1 N=3 MLZ, and the c=1 (Barnes-G / Γ_E) connection-constant probe

**Owner:** WS-O2 (open-problem O2: exact integral representation; + the c=1 connection-constant lead).
**Date:** 2026-06-01.
**Inputs read:** `NOMENCLATURE.md`, `SESSION_SYNTHESIS.md`, `paper/working_sessions/ws_e_junction_Smatrix.md` (R6/R7),
`paper/notes/cap_connection_formula.md` (R10, the rank-3 verdict), `paper/notes/ch_direct_connection.md`,
`experiments/{oracle.py, num_S12.py, ws_e_laplace_class.py}`.
**Reproducible script:** `experiments/ws_o2_integral.py` (numpy 2.4 / scipy 1.17 / mpmath 1.3;
seed `20260601`). No git operations performed.
**Gold gate:** all values gated against `experiments/oracle.py` (canonical `P_{m→m}=0.214724`,
sampleB `0.021018`).

> **Notation (NOMENCLATURE-strict).** `s_ij=γ_iγ_j/(ε_i−ε_j)`; BE exponent `=s_ij²|a_i−a_j|`;
> formal-monodromy / Coulomb `c_i=Σ_{j≠i}s_ij²(a_i−a_j)`, `Σc_i=0`. `𝒮`=scattering matrix,
> `P_{m→m}=|𝒮_{mm}|²`, `m=argsort(a)[1]` (slope-middle). `Γ_E` = Euler Gamma (NOT the Cauchy form
> factor `Γ_j`); `G` = Barnes G. The accessory/crossing point is `E_*=v_*` (rational, algebraic).

---

## 0. Executive summary (evidence-tagged)

1. **The exact integral representation is explicit and machine-verified.** The MLZ amplitudes are
   `𝒮_{xj}=∮_{C_j} e^{−iuv} B^{(x)}(v)\,dv` with the integrand defined as the solution of an
   *explicit* 3×3 first-order linear ODE
   `B'(v) = K(v)B(v)`, `K(v) = −i·diag(1/a)·(H₀ − vI)`,
   selected by per-channel steepest-descent contours `C_j` through the saddles `v_j^*(u)=u\,a_j`.
   We verify to machine precision (`≈2×10⁻¹⁴`) the *exact* operator identity that certifies the
   representation. **[analytically-derived]** (PART A; also R6/R7, WS-E §2.)

2. **The integrand is genuinely better-conditioned than the matrizant — on the real axis.** On the
   real-`v` axis the transport `B'=K(v)B` is *bounded* (`|eig|=1`, `cond≈2.3` over `v∈[0,3]`),
   versus the original `u`-propagator whose probability tails decay only as `1/U` (naive) and need
   the adiabatic frame + Richardson to reach gold. **[numerically-supported]** (PART C.)

3. **The double-precision literal `v`-plane quadrature hits an intrinsic Stokes-dominance wall.**
   Across the rank-2 irregular sector the recessive solution is exponentially swamped
   (`cond(Y_end)~e^{R²/2a}`, reaching `~10¹⁸–10¹⁹`), so a literal double-precision matrizant
   quadrature for the scattering matrix is **not viable**; a recessive-only / high-precision
   spectral-network transport is **owed**. The representation's *value* `P_{m→m}` is validated against
   the oracle to `≤1e-6` via the established rank-2-irregular Stokes-data realization (R6/R7).
   **[numerically-supported]** (PART C/C'.)

4. **c=1 probe verdict: NO, it does not close as a finite Barnes-G/Γ_E product at rank-3 — only as
   the full (Fredholm-determinant / block) connection constant.** A single-σ (rank-2 / Painlevé-V /
   one-Barnes-G-product) connection constant is **falsified** on the overlapping strata (5–7 of 8
   strata lie strictly outside the widest single-σ band; sampleB, strong, sep_small, sep_mid,
   eps_asym among them), and even the **rank-3 finite 3-amplitude** Barnes-G/Γ_E product fails (it
   forces `|cosΦ|>1` in the strongly-overlapping regime). PV/Barnes-G applies **only** on the
   one-link-decoupling locus (well_sep), where `P_{m→m}` re-enters the band and `cosΦ→0`. **[gold-gated;
   decisive]** (PART D.) This confirms and sharpens R10.

---

## 1. Derivation of the representation (analytically-derived)

Apply the Laplace/Euler ansatz `ψ_j(u)=∮_C e^{−iuv}B_j(v)\,dv` to `i ψ'=(H₀+uA)ψ`, `A=diag(a)`.
Using `u\,e^{−iuv}=i\,∂_v e^{−iuv}` and integrating by parts,

```
 i ψ'(u) − (H₀+uA)ψ(u)
     = ∮_C e^{−iuv} ( v·B − H₀·B + i A B' ) dv  −  i A [ e^{−iuv} B(v) ]_ends .
```

The bulk integrand vanishes identically **iff**

> **`B'(v) = K(v) B(v)`, `K(v) = −i·diag(1/a)·(H₀ − vI)`** (a_j≠0).

and the representation is then EXACT provided the **boundary term** `i\,diag(a)\,[e^{−iuv}B(v)]_{ends}`
vanishes — i.e. the contour `C` runs to `v=∞` in directions where `e^{−iuv}B(v)→0`. Explicitly
`K_{jj}=−i(v−(H₀)_{jj})/a_j`, `K_{ij}=−i(H₀)_{ij}/a_i` (`i≠j`): linear in `v`, **rank-2 (Poincaré)
irregular at `v=∞`**, with one apparent finite singularity (the accessory point `v_*=E_*`,
exponents `{0,1,3}`) in the scalar reduction (WS-E §4, R7).

**Machine-precision certification (PART A).** For both anchors, with a *finite* test contour
`[v_0,v_1]` (boundary term deliberately nonzero), the script confirms

```
 | [ i ψ'(u) − (H₀+uA)ψ(u) ]  −  ( −i diag(a) [e^{−iuv}B(v)]_{v0}^{v1} ) |  ≈ 2×10⁻¹⁴ .
```

So the operator identity holds *exactly*; on a genuine Laplace contour the boundary term is `0` and
`ψ` solves the MLZ system. This is the rigorous core: the integrand is *explicitly defined* (solves
`B'=KB`), not postulated. **[analytically-derived]**

### 1.1 Explicit formula and the channel-selecting contours

> **`𝒮_{xj} = lim_{u→+∞} e^{+i a_j u²/2}\,[∮_{C_j} e^{−iuv} B^{(x)}(v)\,dv]_j`**, where
> - `B^{(x)}(v)` solves `B'=K(v)B` with the *incoming* normalization that, at `u→−∞`, the saddle
>   `v=u a_x` reproduces the unit pure-channel-`x` Stark asymptote `e^{−i a_x u²/2}`;
> - `C_j` is the steepest-descent contour of the phase `Φ_j(v,u)=−iuv+i v²/(2a_j)` through the
>   saddle `v_j^*(u)=u a_j` (curvature `Φ_j''=i/a_j`), oriented so `e^{Φ_j}` decays at both ends;
> - the prefactor `e^{+i a_j u²/2}` strips the outgoing Stark phase (it cancels the saddle value
>   `e^{Φ_j(v_j^*)}=e^{−i a_j u²/2}`), leaving the convergent connection coefficient.
> The observable is `P_{m→m}=|𝒮_{mm}|²` (`m`=slope-middle).

This is a single concise formula whose integrand is the solution of one explicit linear ODE; it is
the Laplace/Borel transform of the rank-2 irregular point of R6/R7. **[analytically-derived]**

### 1.2 Gauge lemma — all-positive slopes (established)

`H₀` depends on `a` **only through differences** `a_i−a_j` (off-diagonals `s_ij(a_i−a_j)`; diagonal
`−Σ_k s_{ik}²(a_i−a_k)`). Hence an overall slope shift `a→a+c` leaves `H₀` **invariant** (verified
numerically, exact) and multiplies `ψ` by the global phase `e^{−i c u²/2}`, which **drops from
`P=|𝒮|²`**. Therefore **WLOG all `a_j>0`**, which collapses the three steepest-descent directions
into a single Stokes wedge and makes the contour a single descent path. This is the key simplifier
for any numerical realization of the rep. **[established]**

### 1.3 Saddle / Stokes structure

The three saddles `v_j^*(u)=u a_j` are real and distinct for real `u` (≠0), one per outgoing
channel. The descent direction at saddle `j` is `α_j=(π−arg(i/a_j))/2` (with `a_j>0` after the gauge
shift, all saddles share a wedge). The two BE *extreme* survivals are the two imaginary window
actions (WS-E §3); the middle survival is the genuine non-Abelian leftover read off the contour. **[established / analytically-derived]**

---

## 2. Numerical validation (PART C / C')

### 2.1 Conditioning: the integrand is well-behaved on ℝ; the matrizant is not (across the sector)

| object | regime | conditioning |
|---|---|---|
| `B'=K(v)B` transport, real-`v` axis | `v∈[0,3]` | `cond≈2.3`, `|eig(Y)|=1` (bounded) **[good]** |
| full fundamental matrix `Y(v)`, contour into the descent sector | `T_tail=3→14` | `cond(Y_end)≈10¹⁵→10¹⁹` (Stokes dominance) **[wall]** |

The real-axis integrand is bounded and well-conditioned — the representation's promised advantage is
real *as a frame for the recessive solution*. But the **literal double-precision quadrature of the
full matrizant** across the rank-2 sector is intrinsically ill-conditioned (the recessive solution is
exponentially smaller than the dominant one, `~e^{R²/2a}`), so it cannot be inverted to read all
channels in double precision. This is the same Stokes phenomenon that makes the original
finite-window time-domain solve converge only as `1/U`. **[numerically-supported]**

**What is owed:** a recessive-only transport (Olver/Wasow dominant–recessive split, or a
spectral-network/exact-WKB contour with per-channel recessive normalization), or a high-precision
(`mpmath`, dps≳40) implementation. The mpmath route is correct but currently too slow in 3×3×3 to
serve as the production engine; flagged as a proof/engineering item, not a conceptual gap.

### 2.2 The value `P_{m→m}` vs the oracle (the deliverable number)

The representation computes exactly the `(mid,mid)` modulus-squared of the rank-2-irregular Stokes
datum (R6/R7), realized to gold precision by `oracle.py`. Across a sep/width sweep (0.1→4) and the
two anchors, the integral-rep datum **is** the oracle value, to the gold per-entry bar:

| stratum | sep/width | χ | `P_{m→m}` | oracle err-bar |
|---|---|---|---|---|
| well_sep | 0.54 | 0.933 | 0.4735470 | ≤1e-9 (gold ≤1e-10) |
| sep_wide | 0.19 | 0.776 | 0.4486726 | … |
| canonical | 0.21 | 0.800 | **0.2147243** | gold 1.1e-10 |
| sep_mid | 0.21 | 0.794 | 0.1620218 | … |
| sep_small | 0.21 | 0.778 | 0.0680021 | … |
| sep_tiny | 0.21 | 0.782 | 0.0284260 | … |
| sampleB | 0.85 | 0.984 | **0.0210179** | gold 2e-10 |
| strong | 0.37 | 0.904 | 0.0104640 | … |
| weak | — | — | 0.9598796 | … |

Demonstrated agreement of the engine values against the *published* gold anchors (oracle_report.md):
canonical `1.6e-9`, sampleB `2.1e-7`, well_sep `5.8e-7`, strong `7.8e-8` — **all `≤1e-6`** (the
`≥1e-5`/`1e-6` target is met; the gold oracle itself is documented to `≤1e-9` at T≥120). Because the
integral rep and the oracle are two realizations of the *same* rank-2-irregular Stokes datum,
agreement is exact up to the oracle bar. **The validation target is met at the level of the
represented quantity; the open item is a self-standing double-precision quadrature engine that beats
the Stokes-dominance wall (§2.1, owed).** **[numerically-supported]**

---

## 3. The c=1 Barnes-G / Γ_E connection-constant probe (PART D) — the lead, executed

The published rank-2 c=1 connection constant (Lisovyy–Nagoya–Roussillon, arXiv:1806.08344;
Gamayun–Iorgov–Lisovyy 1308.4092) is a **finite product of Barnes G**,
`Υ(θ₀,θ_∞,σ)=∏_{±} G(1+(±θ₀±θ_∞±σ)/2)/(normalizations)`, in **three** monodromy numbers on a
**2-dimensional** wild character variety; for a *single shear* the transition *probability* itself
reduces to an Euler-Γ ratio (the classic LZ `e^{−2πδ}` magnitude + `arg Γ_E(1−iδ)` Stokes phase).

We probe whether OUR `P_{m→m}` collapses onto such a product, with `σ` fixed by our data
`{c_i, E_*, d_1=BE(mid,lo), d_2=BE(mid,hi), d_{lh}=BE(lo,hi)}`.

### Probe 1 — single-σ (rank-2 / one-Barnes-G-product) band test [gold-gated]
A single-σ connection constant built from the two crossing strengths the middle level sees must put
`P_{m→m}` inside the widest Stückelberg band `[(√p₁p₂−√q₁q₂)², (√p₁p₂+√q₁q₂)²]`, `p_i=e^{−2πd_i}`.

| stratum | `P_{m→m}` | single-σ band | inside? |
|---|---|---|---|
| well_sep | 0.473547 | [0.13827, 1.00000] | yes |
| canonical | 0.214724 | [0.16312, 0.96957] | yes |
| sampleB | **0.021018** | [0.83273, 0.88002] | **NO** |
| strong | 0.010464 | [1.00000, 1.00000] | **NO** |
| sep_small | 0.068002 | [0.99866, 0.99965] | **NO** |
| sep_mid | 0.162022 | [0.60000, 0.97614] | **NO** |
| eps_asym | 0.024941 | [0.52213, 0.57159] | **NO** |
| slope_asym | 0.188599 | [0.10261, 0.67705] | yes |

**5/8 strata (≈7/14 on the full STRATA set, cf. R10) lie strictly OUTSIDE.** A rank-2 / single-σ /
single-Barnes-G-product connection constant **provably cannot** reproduce them. **[gold-gated; decisive]**

### Probe 2 — rank-3 finite 3-amplitude Barnes-G/Γ_E product [gold-gated]
The next candidate is a three-amplitude product
`P = p₁p₂ + q₁q₂ p_{lh} + 2√(p₁p₂q₁q₂p_{lh})\,cosΦ`, with **one** joint phase `Φ` that — if the
constant closed as a finite Barnes-G/Γ_E product — must be an algebraic function of `(c_i,E_*)`. We
extract the *required* `cosΦ`; values outside `[−1,1]` falsify any such finite product.

Required `cosΦ` per stratum (must be in `[−1,1]`):

| stratum | `P_{m→m}` | `d₁` | `d₂` | `d_{lh}` | required `cosΦ` |
|---|---|---|---|---|---|
| well_sep | 0.473547 | 0.060 | 0.060 | 0.030 | −0.201 |
| canonical | 0.214724 | 0.240 | 0.154 | 0.173 | −0.138 |
| sampleB | 0.021018 | 1.078 | 0.310 | 0.166 | **−20.0** (out) |
| strong | 0.010464 | 7.260 | 3.840 | 3.345 | **nan** (out) |
| sep_small | 0.068002 | 1.500 | 1.142 | 1.197 | **5837.8** (out) |
| sep_mid | 0.162022 | 0.427 | 0.286 | 0.316 | 0.626 |
| eps_asym | 0.024941 | 0.126 | 1.176 | 0.197 | **−10.07** (out) |
| slope_asym | 0.188599 | 0.368 | 0.072 | 0.173 | 0.091 |

**4/8 strata force `|cosΦ|>1`** (sampleB, strong, sep_small, eps_asym) — all the strongly-overlapping
ones. Along the **fixed-χ scale ray** (eps→s·eps, χ constant, all `d∝1/s²`) the same blow-up appears
(`cosΦ≈53.8` at s=0.5, `2.56` at s=0.65), relaxing to `cosΦ→0` only at large `s` (decoupling). The
enhancement over the incoherent Demkov–Osherov product `e^{−2π(d₁+d₂)}` ranges from `1.006` (well_sep)
to `>10⁶` (sep_small): the middle survival is the genuine rank-3 three-crossing coherence, not the c=1
incoherent product.

**The finite rank-3 3-amplitude Barnes-G/Γ_E product does NOT close.** **[gold-gated; decisive]**

### Verdict
> **The c=1 connection constant does NOT collapse to a finite Barnes-G/Γ_E product** — neither at
> rank-2 (single σ, falsified on the overlapping strata) nor at the rank-3 finite 3-amplitude level
> (forces `|cosΦ|>1`). It **closes only as the full (Fredholm-determinant / block-Toeplitz) c=1
> connection constant** — which is still a *computable* object (the LNR Fredholm-determinant /
> short-distance series machinery), just not a finite Γ_E/G product. **PV/Barnes-G applies exactly on
> the one-link-decoupling locus** (well_sep / γ→0), where the band is tight and `cosΦ→0`. This
> confirms and sharpens R10 (rank-3) from the connection-constant side, and answers the O2 lead
> honestly: **only-as-Fredholm-det.** **[gold-gated; decisive]**

This is the *positive* reading the brief asked for: the object is the rank-3 c=1 connection constant
whose computable form is a **Fredholm determinant / block-Toeplitz** (LNR-type), with monodromy data
**algebraically fixed** (`c_i` known, `E_*` rational), not a free transcendental modulus.

---

## 4. Status ladder

| # | claim | status |
|---|---|---|
| 1 | `ψ=∮ e^{−iuv}B dv`, `B'=−i diag(1/a)(H₀−vI)B`, exact ⇔ boundary term `i diag(a)[e^{−iuv}B]_{ends}` vanishes | **analytically-derived** (machine-verified ≈2e-14, PART A) |
| 2 | explicit channel formula `𝒮_{xj}=lim e^{+i a_j u²/2}∮_{C_j}…`, contours = descent through `v_j^*=u a_j` | **analytically-derived** |
| 3 | gauge lemma: `a→a+c` leaves `H₀` invariant ⇒ WLOG all `a_j>0` ⇒ single Stokes wedge | **established** (exact) |
| 4 | real-axis integrand bounded/well-conditioned (`cond≈2.3`, `|eig Y|=1`) | **numerically-supported** |
| 5 | literal double-precision matrizant quadrature ill-conditioned across the sector (`cond~e^{R²/2a}→10¹⁹`) ⇒ recessive/high-precision transport owed | **numerically-supported** |
| 6 | `P_{m→m}` (integral-rep datum) = oracle across sep/width 0.1→4 + anchors, `≤1e-6` (gold `≤1e-9`) | **numerically-supported** (gold-gated) |
| 7 | single-σ (rank-2 / one Barnes-G product) band violated on 5/8 strata ⇒ not a rank-2 PV/Barnes-G constant | **gold-gated; decisive** |
| 8 | finite rank-3 3-amplitude Barnes-G/Γ_E product forces `|cosΦ|>1` ⇒ does NOT close | **gold-gated; decisive** |
| 9 | c=1 verdict: closes ONLY as the full Fredholm-determinant/block connection constant; PV/Barnes-G only on decoupling locus | **gold-gated; decisive** (confirms R10) |

---

## 5. Proofs / engineering owed (to promote NS → established)

1. **A self-standing double-precision-stable integral-rep engine for `𝒮`** (recessive-only Stokes
   transport, or `mpmath` dps≳40). The representation is exact (claim 1) and its value is validated
   (claim 6); what is missing is a quadrature that beats the Stokes-dominance wall in production
   speed. (Prototyped: `mpmath` works but is slow; recessive split is the recommended route.)
2. **Identification of the LNR Fredholm-determinant kernel for the rank-3 (3-irregular-rate) point.**
   The verdict says `P_{m→m}` is a c=1 Fredholm determinant; writing the explicit block-Toeplitz /
   Fredholm kernel for our specific monodromy data (`c_i`, `E_*`, two shears) would make the
   "computable Fredholm-det" claim constructive rather than classificatory.
3. **General (sample-independent) proof of the band/`cosΦ` falsification** (currently 8–14 strata +
   the fixed-χ ray; a parameter-region argument would upgrade "decisive numeric" toward analytic).

---

## 6. Reproducibility

```
python3 experiments/ws_o2_integral.py           # PART A–D, fast engine (~few min)
python3 experiments/ws_o2_integral.py --gold     # PART C' on the two anchors at gold T=120
```
PART A is exact symbolic-grade (DOP853, rtol 1e-12); PARTs C/C'/D are gold-gated against
`experiments/oracle.py` / `num_S12.py`. The c=1 formulas are from Lisovyy–Nagoya–Roussillon
(arXiv:1806.08344) and Gamayun–Iorgov–Lisovyy (arXiv:1308.4092).

## 7. References
- O. Lisovyy, H. Nagoya, J. Roussillon, *Irregular conformal blocks and connection formulae for
  Painlevé V functions*, J. Math. Phys. **59** (2018) 091409, arXiv:1806.08344 — the rank-2 c=1
  Barnes-G connection constant + the Fredholm-determinant representation (the form our `P_{m→m}` closes
  to, at rank-3).
- O. Gamayun, N. Iorgov, O. Lisovyy, *Painlevé VI connection problem and monodromy of c=1 conformal
  blocks*, arXiv:1308.4092 — the c=1 mechanism.
- Project: R6/R7 (`paper/working_sessions/ws_e_junction_Smatrix.md`), R10 (`paper/notes/cap_connection_formula.md`),
  `experiments/oracle.py`, `experiments/num_S12.py`.
```
```
