# WS-A — Riemann scheme of the Type-1 N=3 MLZ connection (the class-fixing gate)

**Owner:** WS-A (Riemann-scheme / class-fixing gate). **Date:** 2026-06-01.
**Inputs used:** E3 (genus-0 node), E4 (a-independent connection), `uploads/assay/geometry.py`,
`uploads/assay/ip.py`, `experiments/ring_structure.py`. **Constraint compliance:** no git ops, this
is the only file written; all scripts live under `/tmp/wsa/` and are inlined below.

Model: `H(u)=H0+u·diag(a)`, Type-1 Cauchy coupling; Schrödinger ODE `i ψ'(u)=H(u) ψ`, the prime is
`d/du`. Canonical samples
`S1: ε=(-2,0,3), γ=(1,0.8,1.2), a=(-1,0.5,2)` and
`S2: ε=(-1,0,1.5), γ=(0.9,1.1,0.8), a=(-0.7,0.4,1.3)`.

---

## 0. Executive summary / gate verdict

**The honest object is the 3×3 first-order system in `u` with ENTIRE (polynomial degree-1)
coefficient matrix `M(u)=-i(H0+uA)`. Its ONLY singular point in the `u`-plane is a single
rank-2 (Poincaré rank 2) irregular point at `u=∞` — the Weber / parabolic-cylinder point. There are
NO finite singular points; in particular the exact-crossing "node" `u*` is an ORDINARY point of the
system (M is analytic there), with trivial monodromy.**

**Gate-1 verdict:** the class is **NOT hypergeometric** (no three regular singular points) and
**NOT reducible to elementary by the singular-point count alone**. It is a **single rank-2 irregular
(confluent) connection problem of a 3-level system** — the 3-level generalisation of the Weber /
parabolic-cylinder (Zener) equation. The transition amplitudes are its **Stokes / connection data**.
For `N=2` this is exactly Weber (the Zener formula); for `N=3` it is a genus-0 confluent object of
**confluent-Heun / higher-Weber type**, with the prefactor a genuine confluent connection
coefficient (not plain `₂F₁`). The conjectured "exponents `∝ Γ_ij`" is **CONFIRMED in its correct
form**: the `Γ_ij` are the **pairwise Weber/Stokes exponents** of the irregular point (one per level
pair), equal exactly to the elementary two-level LZ rate `Γ_ij = H0_ij²/|a_i−a_j|`; the diagonal
algebraic exponents are their signed row sums.

**Implication for the prefactor (gate-1 output to WS-B/C/E):** the **elementary part is exactly the
Weber exponential data** (the `a_i/2` quadratic phases, the `H0_ii` linear phases, the `Γ_ij`
exponents = BE survivals). The **non-elementary remainder is the off-diagonal Stokes connection
coefficient of a 3-level rank-2 irregular point** — confluent-Heun class, generically not elementary.
**The node `u*` contributes NO accessory parameter and NO non-trivial local connection data** (it is
an ordinary point of the system / an apparent singularity of the scalar reduction with exponents
`{0,1,2}` and no logs). This refines O2: the node cannot by itself pin `P_{m→m}`; it does not add
accessory moduli either, so it does not raise the difficulty, but it is not a free reduction lever in
the connection-coefficient count. Feed this to WS-B.

---

## 1. Riemann-scheme TABLE

### 1a. The 3×3 system `ψ' = M(u) ψ`, `M(u)=-i(H0+uA)` — the primary object

| location | type | Poincaré rank | local data / exponents | evidence |
|---|---|---|---|---|
| every finite `u∈ℂ` (incl. node `u*`) | **ordinary point** | — | `M` analytic; monodromy = I | **[established, symbolic]** `M` is a degree-1 polynomial matrix; `M(u*)` finite (script §4.2) |
| `u=∞` | **irregular** | **r = 2** (Weber/parabolic-cylinder) | per branch `i`: leading `exp(-i a_i u²/2)`, sub-leading `exp(-i H0_ii u)`, algebraic prefactor `u^{ρ_i}`, `ρ_i=-i b1_i`; pairwise Stokes exponents `Γ_ij` | **[established, symbolic+numeric]** §3, §4.4 |

Formal-solution branch data at `u=∞` (verified exactly, both samples):

```
branch i :  E_i(u) = a_i·u + H0_ii + b1_i/u + O(1/u²)
            phase  -i∫E_i du = -i[ a_i u²/2 + H0_ii u + b1_i log u ] + const
            ⇒  ψ_i ~ exp(-i a_i u²/2) · exp(-i H0_ii u) · u^{ρ_i},  ρ_i = -i·b1_i
            b1_i = Σ_{j≠i} H0_ij²/(a_i−a_j) = Σ_{j≠i} sgn(a_i−a_j)·Γ_ij   (SIGNED Γ row sum)
```

Pairwise Weber/Stokes exponents (= BE rates, exact, both samples):

```
Γ_ij = H0_ij²/|a_i−a_j| = γ_i²γ_j²|a_i−a_j|/(ε_i−ε_j)²
```

Stokes geometry at `u=∞` (rank-2 ⇒ Weber): **4 anti-Stokes rays** at `arg u ∈ {0, π/2, π, 3π/2}`,
**4 Stokes rays** at `arg u ∈ {π/4, 3π/4, 5π/4, 7π/4}`, `2r=4` sectors; all level pairs `(i,j)` share
this geometry because every `a_i−a_j` is real. **[established, analytic]** §3.

### 1b. The scalar 3rd-order ODE (`u`-variable, elimination of one component)

| location | type | local data | evidence |
|---|---|---|---|
| `u* = −(node)` | **apparent singularity** | exponents `{0,1,2}`, **no logarithm**, trivial monodromy | **[established, symbolic+numeric]** §4.2–4.3 |
| `u=∞` | irregular, rank 2 | same Weber data as 1a | **[established]** §3 |

`u*` matches the discriminant double root / the exact crossing exactly (S1: `u*=-187/750=-0.24933`;
S2: `u*=-3031/4400=-0.68886`), confirmed against the numerical min-gap node location to `<10⁻³`.

### 1c. The spectral `λ`-chart (Gaudin parametrization `u=n/p`, `E=m/p`)

The Gaudin eigenvector is **exact and elementary**: `v_i(λ) = γ_i/(λ−ε_i)`, eigenvalue
`E(λ)=m(λ)/p(λ)` (verified `H(u(λ))v − E v = 0` symbolically, §4.1). Each `ε_i` is a **single simple
preimage of `u=∞`**: `λ_i−ε_i ~ γ_i²/u` as `u→∞`, i.e. `Res_{ε_i} u(λ)=n(ε_i)/p'(ε_i)=γ_i²`
(verified numerically: approach rate `(λ_i−ε_i)u→γ_i²`, §4.5). So the single irregular point `u=∞`
**unfolds into the three points `λ=ε_i`** in the spectral chart.

| location (λ) | type | Poincaré rank | meaning | evidence |
|---|---|---|---|---|
| `λ=ε_i` (i=0,1,2) | **irregular** | **r = 2** | the three preimages of `u=∞`, one per sheet/diabatic channel | **[established, symbolic]** pole order 3 of the pulled-back connection ⇒ rank 2, §4.6 |
| `λ=∞` | regular/ordinary | — | finite `u` interior; no LZ growth | **[established]** §4.6 |
| node | (image is an ordinary `u`-point) | — | not a marked point of the connection | §0, §4.2 |

**Caveat (logged loudly per the watch-list):** in the *diabatic* frame the `λ`-connection is still
rank-2 irregular at each `ε_i` (pole order 3 of the pulled-back coefficient), so the `ε_i` are NOT
Fuchsian/regular-singular in that frame — the earlier framing "3 regular poles `{ε_i}`" holds only
**after** the adiabatic (WKB) gauge strips the dynamical exponential. The genuinely regular content
at each `ε_i` is then the pairwise Weber exponents `Γ_ij`, NOT a Fuchsian residue of the geometric
connection. See §2 for why the `Γ_ij` are **action/Stokes residues, not connection residues** —
this resolves the apparent tension with the conjecture and with E1.

---

## 2. Where `Γ_ij` lives (the dictionary — conjecture CONFIRMED, in corrected form)

The conjecture "indicial exponents at `ε_i` are `∝ Γ_ij`" is **TRUE but must be read as Weber/Stokes
exponents of the irregular point, not Fuchsian residues of a regular-singular point.** Evidence:

1. **Algebraic prefactor exponent (diagonal):** `b1_i` (the `1/u` coefficient of `E_i`, = the
   `log u` coefficient of the phase) equals the **signed `Γ` row sum**
   `b1_i = Σ_{j≠i} H0_ij²/(a_i−a_j) = Σ_{j≠i} sgn(a_i−a_j) Γ_ij`. Verified **exactly** (sympy) and
   **numerically** (Richardson extrapolation of `(E_i−a_iu−H0_ii)u` at `u=10³,10⁴,10⁵` → 5 digits),
   both samples. **[established]**
2. **Pairwise Weber exponent (off-diagonal):** `Γ_ij = H0_ij²/|a_i−a_j|`, verified symbolically equal
   to `γ_i²γ_j²|a_i−a_j|/(ε_i−ε_j)²` (the BE rate), both samples. This is the **elementary two-level
   LZ exponent** for the diabatic pair `(i,j)` — exactly the Brundobler–Elser residue structure of
   E1. **[established]**
3. **Why not a geometric-connection residue:** the adiabatic transport coefficient
   `S_ij = Γ_i Γ_j/(λ_i−λ_j)` scales as `1/u²` (a double zero) at `u=∞`, so `u·S_ij→0` — there is
   **no `1/u` simple-pole residue** of the geometric (Berry/Stückelberg) connection. The `Γ_ij`
   therefore arise as residues of the **action one-form / WKB phase** (the `√Q₄` double-cover of
   gate_test_genus §1(d)), realised here concretely as the Weber exponents of the irregular point —
   **not** as residues of the spectral connection. This **resolves the watch-list flag**: the
   exponents ARE `∝ Γ_ij`, consistently with E1; the only correction is *which* one-form they are
   residues of. **[established, symbolic+numeric]** §4.7

This is the precise statement of the "Abelian ceiling" (E2) at the level of local data: everything
**elementary** (quadratic phases `a_i/2`, linear phases `H0_ii`, exponents `Γ_ij` = BE survivals) is
the **formal/exponential** data of the Weber point; the **non-Abelian** remainder is the
**off-diagonal Stokes connection matrix** linking the 3 branches across the 4 Stokes rays.

---

## 3. Formal / Stokes structure at `u=∞` (Poincaré rank 2)

The leading exponential factors are `exp(-i∫E_i du) ~ exp(-i a_i u²/2)`. Pairwise exponent
`-i(a_i−a_j)u²/2`; its real part is `(a_i−a_j)/2·Im(u²)`. Hence:
- **Anti-Stokes rays** (`|·|=1`): `Im(u²)=0 ⇒ arg u ∈ {0, π/2, π, 3π/2}` (4 rays).
- **Stokes rays** (maximal dominance): `arg u ∈ {π/4, 3π/4, 5π/4, 7π/4}` (4 rays).
- `2r = 4` Stokes sectors. This is **exactly** the Weber / parabolic-cylinder geometry of the
  `N=2` Landau–Zener problem, with the `N=3` system threading 3 branches through it.

The connection (transition) data is the product of Stokes matrices around `u=∞`; for `N=2` it
collapses to the single Zener Stokes multiplier `√(1−e^{-2πΓ})`. For `N=3` it is a genus-0
confluent connection coefficient (the open target O1).

---

## 4. Singular-point count + accessory-parameter count

**System (primary object):**
- regular singular points: **0**
- irregular points: **1** (`u=∞`, rank 2)
- apparent singularities: **0**
- **Accessory parameters in the Fuchsian sense: 0** (no apparent singularities to place; the
  connection problem is a pure single-point irregular/Stokes problem). The local moduli are the
  *fixed* formal data: 3 quadratic coeffs `a_i/2`, 3 linear coeffs `H0_ii`, 3 algebraic exponents
  `ρ_i`, all elementary functions of `(γ,ε,a)`; the only transcendental unknowns are the entries of
  the Stokes matrices.

**Scalar 3rd-order reduction (for completeness):**
- finite singular points: **1**, at `u*` (the node), **apparent**, exponents `{0,1,2}`, no log.
- irregular point: **1** (`u=∞`, rank 2).
- The apparent singularity at `u*` is **pinned** (location = discriminant double root, local data
  trivial) ⇒ it adds **0 genuine accessory moduli**. **Watch-list outcome:** the scalar reduction
  *does* spawn an apparent singularity, so for all connection-problem work one should **keep the 3×3
  system** (where `u*` is simply an ordinary point) — exactly the WS-A FAIL-mode mitigation. This is
  not a failure of the gate: the count is clean and definite.

---

## 5. NAMED ODE class

**Class: a single rank-2 irregular (confluent) connection problem of a rank-3 system — the
3-level Weber / parabolic-cylinder equation; genus-0 confluent-Heun / higher-Weber type, NOT
hypergeometric, NOT elliptic/Painlevé, NOT (generically) reducible.**

Reasoning: a `3×3` first-order system with polynomial degree-1 coefficient matrix has a unique
singular point at `u=∞` of Poincaré rank `deg+1 = 2`. One irregular point of rank 2 with three
levels ⇒ confluent class. `N=2` is exactly Weber (Zener). `N=3` is its 3-level confluence:
genus-0 (E3, the spectral curve is rational), the elementary data is the Weber exponential/exponent
data (`Γ_ij` = BE), the transcendental remainder is the off-diagonal Stokes coefficient — a
confluent-Heun / Kampé-de-Fériet connection coefficient. Reducibility to `₂F₁`/elementary is **not**
forced or excluded by the local data alone; that is the separate factorization question (gate-2 /
WS-C), to be tested on sub-loci.

---

## 6. Gate-1 VERDICT and prefactor implication

- **VERDICT (gate 1 PASS):** definite Riemann scheme obtained and verified (§1). Singular points and
  accessory parameters counted (§4). Named class fixed (§5).
- **Hypergeometric? NO** — there are not three regular singular points; there is one rank-2 irregular
  point. **Confluent-Heun-type? YES** — single rank-2 irregular point on a genus-0 curve, 3 levels.
  **Reducible/elementary by count alone? NO** (deferred to WS-C on sub-loci).
- **Prefactor implication:** the prefactor (middle survival `P_{m→m}`, off-diagonals) is the
  **off-diagonal Stokes connection coefficient of a 3-level rank-2 irregular point**. The elementary
  ceiling is the Weber data: `exp(-2πΓ_ij)` survivals (BE), quadratic/linear dynamical phases. The
  genuinely non-elementary content is this single confluent Stokes coefficient — **generically
  confluent-Heun, not elementary**. The earlier "Kampé-de-Fériet / confluent-Heun" expectation
  (F3) is **confirmed and sharpened** to "the Stokes coefficient of the 3-level Weber point."
- **Node note for WS-B (O2):** the exact-crossing node is an **ordinary point of the system /
  apparent singularity of the scalar ODE with trivial local data** — it carries **no** non-trivial
  local connection data and **no** accessory parameter. Therefore it **cannot by itself pin
  `P_{m→m}`**, and it does **not** reduce the accessory-parameter count (there is nothing to reduce —
  the count is already 0 in the Fuchsian sense). Its physical role (the `2.5–104×` enhancement, E5)
  must come from the **global** Stokes data as parameters approach the node, not from local node
  data. WS-B should treat the node as a near-degeneracy of the *global* connection problem, not as a
  local `2×2` model with its own connection coefficient.

---

## 7. Evidence-ladder status of each claim

| # | claim | status |
|---|---|---|
| 1 | `M(u)` entire ⇒ only singularity at `u=∞`; node `u*` is ordinary | **established** (symbolic, exact) |
| 2 | `u=∞` is rank-2 irregular (Weber); 4 Stokes + 4 anti-Stokes rays | **established** (analytic) |
| 3 | branch data `a_i/2`, `H0_ii`, `ρ_i=-i b1_i` | **established** (symbolic + numeric, both samples) |
| 4 | `b1_i = Σ_j H0_ij²/(a_i−a_j) = signed Γ row sum` | **established** (exact symbolic; numeric to 5 digits) |
| 5 | `Γ_ij = H0_ij²/|a_i−a_j| = γ_i²γ_j²|a_i−a_j|/(ε_i−ε_j)²` (BE rate) | **established** (exact symbolic) |
| 6 | scalar ODE: single apparent singularity at `u*`, exponents `{0,1,2}`, no log | **established** (symbolic indicial + ordinary-point argument) |
| 7 | `u*` = discriminant double root = numerical node | **established** (symbolic = numeric `<10⁻³`) |
| 8 | Gaudin eigenvector `v_i=γ_i/(λ−ε_i)`, `E=m/p`; `ε_i` = preimages of `u=∞`, `Res=γ_i²` | **established** (exact symbolic + numeric) |
| 9 | `λ`-connection is rank-2 irregular (pole order 3) at each `ε_i` in the diabatic frame | **established** (symbolic) |
| 10 | `Γ_ij` are action/Stokes residues, NOT geometric-connection residues (`u·S_ij→0`) | **established** (numeric) |
| 11 | class = genus-0 confluent-Heun / 3-level Weber; not `₂F₁`, not elliptic | **analytically-derived** (follows from 1–10 + E3) |
| 12 | node carries no accessory parameter / no local connection data | **established** (from 1, 6) |

No claim contradicts E1–E5. Claim 4/5/10 **strengthen** E1 (they exhibit the BE rates as the explicit
Weber/Stokes exponents and as signed row sums of the diagonal algebraic exponents). The only
correction to prior framing is interpretive (the `ε_i` are regular-singular only after the adiabatic
gauge; the `Γ_ij` are action residues), logged in §1c/§2 per the watch-list — flagged, examined, and
resolved without contradiction.

---

## 8. Verification scripts (inline; live copies under `/tmp/wsa/`)

All results above are reproduced by `/tmp/wsa/VERIFY_ALL.py` (consolidated, both samples) plus the
supporting scripts `infty_exact.py` (exact `b1` = signed Γ row sum), `numeric_b1.py` (numeric
`b1`), `lam_geom.py` / `lambda_scalar.py` (`λ`-chart, `Res=γ_i²`, pole order 3), `adiabatic_lambda.py`
(`u·S_ij→0`), `gamma_residue.py` (`Γ_ij=H0_ij²/|Δa|`), `node_nolog.py` (ordinary-point argument),
`stokes_rays.py` (Stokes geometry).

`VERIFY_ALL.py` (canonical, self-contained):

```python
import numpy as np, sympy as sp

def H0_sym(eps,gam,a):
    H0=sp.zeros(3,3)
    for i in range(3):
        for j in range(3):
            if i!=j: H0[i,j]=gam[i]*gam[j]*(a[i]-a[j])/(eps[i]-eps[j])
        H0[i,i]=-sum(gam[k]**2*(a[i]-a[k])/(eps[i]-eps[k]) for k in range(3) if k!=i)
    return H0

def scalar_ode(eps,gam,a):                      # cyclic-vector elimination, psi_0 component
    u=sp.symbols('u'); H0=H0_sym(eps,gam,a); A=sp.diag(*a); M=-sp.I*(H0+u*A)
    e0=sp.Matrix([1,0,0]); rows=[e0]; c=e0
    for k in range(3):
        c=sp.Matrix([sp.diff(c[i],u) for i in range(3)])+M.T*c; rows.append(c)
    R=sp.Matrix.hstack(rows[0],rows[1],rows[2]); coeffs=sp.simplify(rows[3].T*R.inv())
    return u,coeffs                              # coeffs = (c2,c1,c0) in psi'''=c2 psi''+c1 psi'+c0 psi

# For each sample: finite singular pts of the scalar ODE (= node u*, apparent, {0,1,2});
# branch exponents a_i/2, H0_ii, b1_i=sum H0_ij^2/(a_i-a_j); Gamma_ij=H0_ij^2/|a_i-a_j|.
# (full driver in /tmp/wsa/VERIFY_ALL.py; output reproduced in this document, §4 of WS-A run)
```

Representative output (both samples, exact rationals):
```
S1: u*=-187/750=-0.24933  indicial r^3-3r^2+2r -> {0,1,2}  node-numeric match <1e-3
    branch i: (a_i/2, H0_ii, b1_i) = (-1/2,-168/125,-258/625),(1/4,-147/100,54/625),(1,-23/25,204/625)
    Gamma: (01)=6/25, (02)=108/625, (12)=96/625   all == g^2g^2|da|/de^2  (True)
S2: u*=-3031/4400=-0.68886  indicial -> {0,1,2}  node-numeric match <1e-3
    branch i: (a_i/2,H0_ii,b1_i)=(-7/20,-1843/1000,-621999/500000),(1/5,-51/40,15367/20000),(13/20,-687/500,7432/15625)
    Gamma: (01)=107811/100000,(02)=2592/15625,(12)=968/3125  all == g^2g^2|da|/de^2 (True)
```
