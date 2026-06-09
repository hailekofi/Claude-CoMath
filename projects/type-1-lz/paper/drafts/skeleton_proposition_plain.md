# Algebraic skeleton ⊕ two-dimensional Stokes core

*Plain-Markdown twin of the new Proposition + Remark folded into
`type1_n3_wildness_proof.tex` (§ between the Stokes corollary and the
Non-rigidity section). Faithful to the LaTeX; no compiler needed.*

---

## Setup (one-line recap)

The Type-1 N=3 MLZ model `i ψ'(u) = (H₀ + uA)ψ`, `A = diag(a₀,a₁,a₂)`, Cauchy/Gaudin
couplings `(H₀)_ij = γ_iγ_j(a_i−a_j)/(ε_i−ε_j)`. Observable: the doubly-stochastic
transition matrix `P_{x→j} = |S_xj|²`, where `S` is the scattering (connection) matrix
between the asymptotic frames at `u = ∓∞`. Channels are ordered by slope:
**lo < mid < hi** (by `a`).

---

## Proposition (Algebraic skeleton and the two-dimensional Stokes core)

The scattering matrix factorizes as

> **S = (formal monodromy) · (exponential torus) · (Stokes).**

The first two factors are fixed by the **local formal data** at `u=∞` (Theorem 1,
Wildness) and depend **rationally / Liouvillian-ly** on `(ε, γ, a)`:

- **(a) Exponential torus** `diag(e^{q_i})`, with `q_i = (i/2) a_i u²`.
  This is a transparent **grading of the three channels by slope** — *linear in the
  commuting-member datum `a_i`*, which **is** the irregular type itself.

- **(b) Formal monodromy.** Read on the doubly-stochastic observable `P`, it is a single
  **oriented 3-cycle** (a circulant permutation). Its orientation is the algebraic ℤ₂
  datum `sgn(π)` — the **parity of the slope-ordering permutation** `π : ε ↦ a`
  (the orientation selector, [OWY]).

Consequently `P` has the **exact normal form**

> **P = Φ_aff( P_lo, P_hi ; σ, b )**
> &nbsp;&nbsp;&nbsp;&nbsp;`└─ algebraic (BE) ─┘  └ Stokes core ┘`

where:

- the two **extreme-slope survivals** are pinned to the **Brundobler–Elser values**
  `P_lo, P_hi = exp(−2π Σ δ)` **exactly** (with `δ_ij = γ_i²γ_j²|a_i−a_j|/(ε_i−ε_j)²`);
- `Φ_aff` is the **affine map** that fixes the remaining seven entries by double
  stochasticity (row/column sums = 1);
- **only the middle (recombining) channel** carries the Stokes constants `{σ, b}`.

The skeleton is realized **concretely** as the **time-ordered product**

> **P_skeleton = M₃ M₂ M₁**

of three elementary **2-level Landau–Zener stochastic factors**, composed in the order
of the **rational diabatic crossing times**

> `u_ij = ( (H₀)_jj − (H₀)_ii ) / (a_i − a_j)`.

This product reproduces `P_lo`, `P_hi`, **and the oriented cycle exactly**.

### Proof (sketch, as in the LaTeX)

- The factorization and **(a)** are the **unramified formal decomposition** of Theorem 1
  ([Wasow], [JMU]): the formal solution at `∞` is `F̂(u) · u^L · e^{Q(u)}` with
  `Q = diag(q_i)` and `F̂` a formal gauge — giving the exponential torus `e^Q` and the
  formal monodromy `e^{2πi L}`.
- **(b)** is the value of that formal monodromy on `|S|²`; the orientation equals the
  **sign of the slope-ordering permutation**, an algebraic function of `(ε, a)`,
  established by the selector computation ([OWY]).
- The **exact pinning** of the two extreme survivals is the **Brundobler–Elser theorem**
  ([BE]). In Stokes-graph terms: the extreme rate `q_lo` (resp. `q_hi`) is **recessive**
  (resp. **dominant**) throughout, so its diagonal entry receives **no Stokes
  correction** — the extreme channel is **diagonally isolated**.
- **Double stochasticity** then determines the remaining seven entries affinely from
  `{P_lo, P_hi, σ, b}` — that is `Φ_aff`. The ordered-product realization is the standard
  independent-crossing (incoherent) composition of pairwise LZ factors in crossing-time
  order. ∎

---

## Remark (The Stokes core is two-dimensional)

The core `{σ, b}` has **transcendence degree two** over the elementary field
`Liouv(ε, γ, a)`: **`b` is not a function of `(P_lo, P_hi, σ)`.**

This is **numerically supported (not claimed as a theorem):**

- across the family, the Jacobian of `(ε,γ,a) ↦ (P_lo, P_hi, σ, b)` attains **full rank 4**;
- the best polynomial fit of `b` on `(P_lo, P_hi, σ)` leaves an **irreducible residual
  ≈ 0.26** (against `std(b) ≈ 0.35`), far above the `10⁻³` integration floor. Were `b`
  determined, the residual would collapse to that floor — it does not.

Hence **σ is the *universally hard* datum** — the only one not elementary in *any* limit —
but **not the sole transcendental**: the wild core of the Type-1 N=3 amplitude is genuinely
**two-dimensional, `{σ, b}`**, dressing the algebraic skeleton.

Proving `σ, b` genuinely **non-Liouvillian** (as opposed to "no closed form found") is a
question of **parameterized differential Galois theory**, beyond the present scope.

---

## How to read the slogan, in one line

> **"S = algebraic skeleton ⊕ 2-D Stokes core"** means:
> *2 continuous algebraic numbers* (the BE survivals `P_lo, P_hi`)
> **+** *2 discrete algebraic data* (the **grading** = which two are protected,
> and the **orientation** = the ℤ₂ cycle direction)
> **+** *the unitarity glue* `Φ_aff`
> **⊕** *2 transcendental Stokes numbers* `{σ, b}` dressing the one recombining channel.

The skeleton is not "two numbers + unitarity" — it carries a **combinatorial backbone**
(grading + oriented 3-cycle) and a **specific factorized matrix form** (the ordered LZ
product `M₃M₂M₁`), all explicit in `(ε,γ,a)`. It even explains *why* the transcendence is
exactly 2-dimensional and *where* it lives (the middle, recombining channel).

---

### References (keys as in the .tex)

- **[Wasow]** W. Wasow, *Asymptotic Expansions for Ordinary Differential Equations*, 1965.
- **[JMU]** Jimbo–Miwa–Ueno, *Monodromy preserving deformation…*, 1981.
- **[BE]** Brundobler–Elser, *S-matrix for generalized Landau–Zener problem*, J. Phys. A.
- **[OWY]** Owusu–Wagh–Yuzbashyan, *The link between integrability, level crossings…*.
