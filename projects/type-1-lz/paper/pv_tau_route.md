# WS-PV — The Painlevé-V / isomonodromic-τ route to `S₁₂` (≡ `P_{m→m}`) + literature & irreducibility

**Owner:** WS-PV (PA-3 + PA-5). **Date:** 2026-06-01.
**Inputs read (in order):** `NOMENCLATURE.md`, `paper/OPEN_PROBLEM.md`, `paper/ws_e_junction_Smatrix.md`,
`paper/ws_a_riemann_scheme.md`, `paper/kz_isomonodromy_picture.md`, `paper/ws_g_stokes_graph.md`,
`experiments/oracle.py`. **Literature:** Gamayun–Iorgov–Lisovyy (PVI, arXiv:1308.4092); Its–Lisovyy–
Prokhorov / "Irregular conformal blocks and connection formulae for Painlevé V" (arXiv:1806.08344);
Iorgov–Lisovyy–Teschner (sine-Gordon/PIII, arXiv:1403.1235); Mazzocco (rank-2 irregular Garnier,
arXiv:nlin/0306020); Bertola–Mazzocco and the rank-3 / higher-Garnier line (arXiv:2503.22198,
2512.24083). **Constraint compliance:** no git ops; the only file written is this note; oracle anchors
re-verified (canonical `P_mid=0.214724`, sampleB `P_mid=0.021018`, both at T=120 to ≤1e-5).

---

## 0. Executive verdict (honest — this is a *structured blocker*, with a named target and a precise gap)

> **The PV/τ route gives a fully-formed *named target* for `S₁₂`: it is a c=1 irregular-conformal-block /
> Barnes-G connection coefficient of an isomonodromic τ-function. But there is a hard, *load-bearing*
> dimensional obstruction between our problem and the published Painlevé-V theory: the published PV
> connection formula (arXiv:1806.08344) is for a `2×2` (rank-2) linear system, whereas the Type-1 N=3
> connection problem is the genuinely `3×3` (rank-3) system whose scalar reduction is *3rd*-order with
> an apparent `{0,1,3}` point (WS-A/WS-E, established). A rank-3 system with one rank-2 irregular point +
> one apparent regular point is NOT classical Painlevé V — its isomonodromic deformation is a higher
> (Garnier-type / rank-3) system that *generically lacks the Painlevé property* (arXiv:2503.22198). So
> the clean "`S₁₂` = the Lisovyy PV connection constant, benchmarked" is NOT available off the shelf. The
> realistic, defensible result is: `S₁₂` is the off-diagonal Stokes/connection coefficient of a rank-3,
> rank-2-irregular isomonodromy problem — the same *family* the c=1 CFT machinery names, but one rank
> above the published closed PV formula. The PV formula applies *exactly* only on the sub-loci where the
> rank-3 system reduces to rank-2 (a decoupling / the apparent point becoming non-apparent), which are
> the BBGY/Demkov–Osherov-type elementary corners — i.e. exactly where `P_{m→m}` is already elementary.**

What IS delivered:
1. The PV connection-formula machinery reconstructed and matched against our monodromy data: the formal
   exponents are the **known** `c_i = Σ_{j≠i} s_ij²(a_i−a_j)` (verified numerically, §3); the role of the
   two Stokes shears `{mid,lo}`, `{mid,hi}` is pinned to the off-diagonal Stokes multipliers of the
   irregular point. [analytically-framed]
2. The candidate closed form *shape*: `S₁₂` = a ratio of τ-functions = a c=1 irregular-conformal-block /
   Barnes-G connection constant of the monodromy data — with the precise statement of why the
   **published PV instance does not directly apply** (rank mismatch, §2). [analytically-framed]
3. The literature map (§4): the PV connection coefficient IS published (Lisovyy et al.) — *for rank 2*.
   The rank-3 confluent-Garnier connection coefficient our problem needs is **not** published in closed
   form. [established from literature]
4. The irreducibility framing (§5): `S₁₂` is the rank-3 confluent (Garnier-9/2-type) Stokes coefficient,
   provably above both ₂F₁/elementary (WS-A) and the rank-2 PV transcendent (this note). [analytically-framed]

What is NOT delivered: a benchmarked closed numerical value of `S₁₂` from a PV formula. The blocker is
the rank-3↔rank-2 gap, stated precisely in §2; it is structural, not a missing computation.

---

## 1. The Painlevé-V linear system and its monodromy data (the published object)

**The PV linear system (Jimbo–Miwa; Its–Lisovyy–Prokhorov).** Painlevé V is the isomonodromy
deformation of a **2×2** linear system on `P¹` with
- one **regular** (Fuchsian, simple-pole) singular point — conventionally at `z=0`, exponents `±θ₀/2`;
- one **rank-2 irregular** singular point at `z=∞`, with formal solution
  `Y ~ (I + O(1/z)) z^{Θ} e^{(t z/2 + …)σ₃}` carrying a formal-monodromy exponent `Θ` (the `z^{Θ}`
  power) and **two Stokes matrices** (rank 2 ⇒ `2r=4` Stokes sectors at ∞, hence 4 Stokes rays; the
  independent Stokes data reduces to **two** Stokes multipliers `s₁,s₂` by the cyclic relation).

So PV's monodromy data = `{θ₀ (regular exponent), θ_∞≡Θ (formal-irregular exponent), s₁, s₂ (two Stokes
multipliers)}`, modulo conjugation — a 2-(complex-)dimensional monodromy manifold (the PV wild character
variety). **[established, literature]**

**The connection formula (arXiv:1806.08344, Its–Lisovyy–Prokhorov / Lisovyy et al.).** The PV
τ-function has expansions at the three critical points `t→0`, `t→+∞`, `t→i∞`. Each expansion is a
**Fourier series of irregular c=1 Virasoro conformal blocks** (a Nekrasov-type sum) whose summand
parameters are the monodromy data; the **connection constant** relating two expansions is an explicit
**ratio of products of Barnes G-functions**,
```
   χ  =  ∏ G(1 + …) / ∏ G(1 + …) ,     arguments  =  (θ₀ ± θ_∞ ± σ)/2 - type combinations,
```
with `σ` the "intermediate" exponent fixed by the Stokes data (the `cos 2πσ` / `e^{2πiσ}` of the
trace-coordinate on the character variety). The Stokes multipliers `s₁,s₂` are *closed-form* functions
of `(θ₀,θ_∞,σ)` (the parametrization of the wild character variety), and `σ` is the one remaining
"accessory" coordinate — fixed for a *specific* solution by its asymptotic boundary condition.
**[established, literature]**

**Why this is the right *shape* for `S₁₂`.** WS-A established `𝒮` = Stokes/connection data of a single
rank-2 irregular point; WS-E/OPEN_PROBLEM established the scalar reduction has exactly one rank-2
irregular point + one apparent regular point. The restart-operator structure (OPEN_PROBLEM addendum)
says the monodromy data = (diagonal formal exponents `c_i`, KNOWN) ⋉ (two Stokes shears in `{mid,lo}`,
`{mid,hi}`), and `S₁₂` is their non-commutative composition. **This is structurally a PV-type connection
constant**: known formal exponents + two Stokes multipliers + one intermediate accessory coordinate,
composed into a Barnes-G ratio. The *form* of the answer is therefore named. **[analytically-framed]**

---

## 2. THE BLOCKER, stated precisely: rank-3 ≠ rank-2 (Painlevé V is one rank too small)

This is the make-or-break gap, and it is dimensional, not a missing calculation.

**Fact A (established, WS-A §1b / WS-E §4).** The Type-1 N=3 connection problem is a **3×3** first-order
system; its scalar reduction is a **3rd-order** ODE with one apparent regular point (exponents
`{0,1,3}`, gap at 2) plus the rank-2 irregular point at ∞.

**Fact B (established, literature §1).** The Painlevé-V linear system is **2×2**; its scalar reduction
is **2nd-order** (one regular + one rank-2 irregular point). Its published connection formula is for
*this* 2×2 system.

**Consequence.** The two are NOT the same isomonodromy problem. The deformation of a 3×3 system with
{one rank-2 irregular + one apparent} is a **higher / rank-3 (Garnier-type) isomonodromy system**, not
PV. The relevant references (Mazzocco arXiv:nlin/0306020; the Garnier-9/2 and 5/2+3/2 reductions,
arXiv:2503.22198; rank-3 Painlevé representations arXiv:2512.24083) make the rank-3 case explicit and —
critically — note that these higher systems **generically lack the Painlevé property** and their fourth-
order deformation equations are NOT PV. So:

> **OPEN_PROBLEM's central conjecture ("expected: Painlevé V") is, on the nose, FALSE for the generic
> rank-3 reduction — the deformation is a higher Garnier-type system, one rank above PV. The published
> Lisovyy PV connection formula therefore does NOT compute the generic Type-1 `S₁₂`.** [analytically-framed,
> literature-grounded]

**Where PV *does* apply (the reduction loci).** The 2×2 PV formula applies *exactly* when the rank-3
system reduces to rank-2. Two mechanisms:
- the apparent `{0,1,3}` point becomes a genuine ordinary point / the 3rd-order ODE factors a 1st-order
  piece off — i.e. one diabatic channel decouples (WS-C trivial-coupling limit; one Stokes shear
  collapses, OPEN_PROBLEM addendum). Then the surviving 2×2 connection problem is literally PV (or
  Weber, for full decoupling), and `S₁₂` is elementary/`₂F₁` — *exactly the corner WS-E found elementary*
  (well-separated / weak strata, `P_{m→m}≈p₁p₂`, the incoherent edge).
- the BBGY tridiagonal / Demkov–Osherov sub-locus (extreme–extreme coupling `s_lo,hi→0`): again a
  reduction to the solvable `₁F₂` cousin. Type-1 is generically NOT on this locus (WS-E §1: extreme–
  extreme coupling is the *largest* on canonical), so this corner is non-generic.

**Net:** PV names the answer on precisely the loci where `P_{m→m}` is *already elementary*; on the generic
overlapping locus (sampleB, ov1, ov2 — where `P_{m→m}` is genuinely enhanced and outside every two-path
band, WS-E §6d) the object is the **rank-3 confluent-Garnier Stokes coefficient**, strictly above PV.
This is consistent with — and sharpens — WS-E's "confluent-Heun, one accessory parameter above ₁F₂".

---

## 3. The monodromy data in model parameters `{γ,ε,a}` (what IS pinned)

The diagonal/formal half of the monodromy data is fully explicit and re-verified numerically here.

**Formal monodromy exponents (KNOWN, the diagonal carrier).** The `z^{ρ_i}` powers at the irregular
point are the signed-BE Coulomb coefficients (WS-A claim 4; NOMENCLATURE):
```
   c_i = Σ_{j≠i} s_ij² (a_i − a_j),     s_ij = γ_iγ_j/(ε_i−ε_j),     Σ_i c_i = 0.
```
Verified (this note, /tmp probe): canonical `c = (−0.4128, +0.0864, +0.3264)`, Σ=0 to 1e-16;
sampleB `c = (−1.2440, +0.7684, +0.4757)`, Σ=0. These are the rank-3 analog of PV's `θ_∞`: a
**traceless diagonal** of formal exponents (PV has one `θ_∞`; we have the SL(3)-type pair `c_i`).
**[established, numeric+symbolic]**

**Pairwise BE/Stokes exponents (KNOWN).** The pairwise weights `s_ij²|a_i−a_j|` set the magnitudes of
the Stokes shears (the BE survivals are their exponentials). Canonical
`(BE_01,BE_02,BE_12)=(0.240,0.173,0.154)`; sampleB `(1.078,0.166,0.310)`. These are the rank-3 analog
of PV's `θ₀` (the regular-point exponent governing the Stokes-multiplier magnitude). **[established]**

**The unknowns `S₁₂` is built from.** Exactly the **two Stokes multipliers / shears** in the carrier
pairs `{mid,lo}` and `{mid,hi}` (WS-G's `12×13` joint; the middle level is the shared subdominant
partner). In PV language these are `s₁,s₂`; here they are the two off-diagonal Stokes coefficients of a
rank-3 irregular point, and `S₁₂≡P_{m→m}` is their non-commutative composition. The **intermediate
exponent `σ`** (the one accessory coordinate that the PV Barnes-G formula needs) is the rank-3 analog of
WS-E's "one accessory parameter above ₁F₂" — and **this is exactly the quantity WS-PA1 is trying to fix
algebraically (PA-2).** If WS-PA1 finds it algebraic in `{γ,ε,a}`, the connection-coefficient program
closes (modulo writing the rank-3 Barnes-G ratio); if transcendental, `S₁₂` is a genuine rank-3
isomonodromy transcendent. **[analytically-framed; the PA-2 gate is the pivot]**

---

## 4. Literature map (PA-5): is the relevant coefficient published?

| object | published closed form? | reference | applies to our `S₁₂`? |
|---|---|---|---|
| PVI connection constant (c=1 conformal blocks, Barnes-G) | **yes** | Gamayun–Iorgov–Lisovyy 1308.4092 | no (3 reg. pts, 2×2; wrong scheme) |
| PIII / sine-Gordon connection constant | **yes** | Iorgov–Lisovyy–Teschner 1403.1235 | no (different confluent 2×2) |
| **PV connection constant** (irregular c=1 blocks; `t→0,+∞,i∞`; Barnes-G ratio) | **yes** | **Its–Lisovyy–Prokhorov / 1806.08344** | **only on the rank-2 reduction loci** (§2) |
| rank-2 irregular Garnier (1 simple pole + Poincaré-rank-2 ∞) deformation | structure known | Mazzocco nlin/0306020 | partial — still 2×2 / `N=2`-type |
| **rank-3 system, 1 rank-2 irregular + 1 apparent** (OUR case) connection coeff | **NO closed form** | Garnier-9/2 etc. 2503.22198, 2512.24083 (note: generically NO Painlevé property) | **this is the object; unpublished** |
| confluent-Heun central connection coefficient (general) | **NO closed form** | standard (Ronveaux; heun.xyz) | the 2nd-order scalar cousin; also open |

**Reading.** The c=1-CFT / τ-function connection program is real, modern, and exactly the right toolbox
— *for rank 2*. Lisovyy et al. solved PVI/PV/PIII/PII connection problems. **Our problem sits one rank
above the published PV result**: a rank-3 confluent (Garnier-9/2-type) connection coefficient, which is
NOT in the literature in closed form, and for which the underlying deformation generically lacks the
Painlevé property. So the honest literature verdict is: **the *named family* (isomonodromic τ /
c=1-type connection constant) is established; the *specific coefficient* our `S₁₂` equals is not
published and is not a classical special function.** [established from literature]

---

## 5. Irreducibility statement (PA-5, the clean "named the object" result)

Assembled from WS-A (rank-3, rank-2-irregular, apparent `{0,1,3}`), WS-E (one accessory parameter above
₁F₂; benchmarked-transcendental, outside every two-path band on the overlapping strata), the §2 rank
argument, and the §4 literature map:

> **`S₁₂` (≡ middle survival `P_{m→m}`) of generic Type-1 N=3 MLZ is the off-diagonal Stokes/connection
> coefficient of the rank-3 linear isomonodromy problem with one Poincaré-rank-2 irregular point at ∞
> (formal exponents `c_i = Σ_{j≠i}s_ij²(a_i−a_j)`, known) and one apparent regular point (exponents
> `{0,1,3}`). Its monodromy data are these known formal exponents plus the two Stokes shears in the
> carrier pairs `{mid,lo}`, `{mid,hi}`, of which `S₁₂` is the non-commutative composition. This object
> is:
> (i) **above ₂F₁ / elementary** — there is no Fuchsian image (no `1/t`, `[K(v₁),K(v₂)]≠0`, WS-E §2)
>     and no tridiagonal/Gaudin structure (extreme–extreme coupling ≠0, WS-E §1); benchmarked outside
>     every two-amplitude interference band on the overlapping strata (WS-E §6d);
> (ii) **above the published rank-2 Painlevé-V transcendent** — the linear system is rank 3, not 2; its
>      deformation is a higher Garnier-type (9/2-type) system that generically lacks the Painlevé
>      property (§2), so the closed Lisovyy PV/Barnes-G connection formula does not compute it (it
>      computes only the rank-2 *reduction* loci, which are the already-elementary decoupling/BBGY
>      corners).
> It therefore equals the rank-3 confluent-Garnier connection coefficient — a named isomonodromic object
> in the c=1-type τ-function family, but not reducible to classical special functions (₂F₁, confluent
> hypergeometric, or the rank-2 PV transcendent), on the generic (overlapping) Type-1 locus.**

Evidence ladder: (i) **established** (WS-A/WS-E symbolic + numeric); (ii) **analytically-framed,
literature-grounded** (the rank count is established symbolically; "deformation lacks Painlevé property
⇒ not PV" rests on the cited rank-3 Garnier literature, which should be coordinator-checked against the
exact `{0,1,3}` scheme by WS-PA1).

---

## 6. What would unblock a benchmarked closed form (handoffs)

1. **PA-2 / WS-PA1 (the pivot):** does the Gaudin/Cauchy data fix the rank-3 accessory coordinate `σ`
   algebraically in `{γ,ε,a}`? If yes, the rank-3 Barnes-G/τ ratio can in principle be assembled and
   benchmarked (target ≤1e-4 on canonical 0.214724, sampleB 0.021018). If no, §5(ii) is the final answer.
2. **Reduction test (decisive, cheap, hands to WS-CH/WS-PA1):** verify symbolically that on the
   decoupling locus (one `s_ij→0`) the 3rd-order ODE factors → 2nd-order with one regular + rank-2
   irregular point = the *literal* PV linear system; then the Lisovyy PV connection constant SHOULD
   reproduce the (now elementary) `P_{m→m}` there. This is the one place the published PV formula is
   directly benchmarkable against our oracle, and confirms the §2 reduction picture. *(I did not run this
   symbolic factorization here — flagged as the highest-value next probe; it is the concrete bridge
   between the published rank-2 theory and our rank-3 object.)*
3. **Rank-3 confluent-Garnier connection coefficient:** the genuinely new mathematics. Not in the
   literature; the honest "named object, unpublished, provably non-classical" result of §5 stands unless
   someone derives it.

---

## 7. Status ladder

| # | claim | status |
|---|---|---|
| 1 | PV linear system = 2×2, one regular + one rank-2 irregular; data = {θ₀,θ_∞,s₁,s₂} | **established** (literature) |
| 2 | PV connection constant = Barnes-G ratio of c=1 irregular-block data (1806.08344) | **established** (literature) |
| 3 | Type-1 N=3 reduction is rank-3 (3rd-order scalar, apparent `{0,1,3}`) — NOT 2×2 | **established** (WS-A/WS-E, symbolic) |
| 4 | ⇒ deformation is higher Garnier-type, generically NOT PV; published PV formula does not apply generically | **analytically-framed** (literature-grounded; coordinator-check vs `{0,1,3}` scheme advised) |
| 5 | formal exponents `c_i=Σs_ij²(a_i−a_j)` known, Σ=0; pairwise BE weights known | **established** (numeric+symbolic, this note) |
| 6 | unknowns = two Stokes shears `{mid,lo},{mid,hi}`; `S₁₂` = their composition; accessory `σ` = PA-2 pivot | **analytically-framed** |
| 7 | PV/Barnes-G formula applies exactly only on rank-2 reduction (decoupling/BBGY) loci = already-elementary corners | **analytically-framed** |
| 8 | irreducibility: `S₁₂` = rank-3 confluent-Garnier connection coeff, above ₂F₁ AND above rank-2 PV transcendent | **analytically-framed + literature** |
| 9 | rank-3 confluent-Garnier connection coefficient is unpublished in closed form | **established** (literature gap) |

No claim contradicts WS-A/WS-E/WS-G or the oracle. This note **corrects** OPEN_PROBLEM's "expected:
Painlevé V" to "expected: a rank-3 *higher*-Garnier system; PV is the rank-2 reduction" and **names**
the target object precisely, while honestly reporting that the off-the-shelf PV closed form does not
reach the generic `S₁₂`.
</content>
</invoke>
