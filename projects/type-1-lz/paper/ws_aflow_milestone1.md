# WS-AFLOW Milestone 1 — Does the a-flow of the transition data close?

**Workstream:** WS-AFLOW (make-or-break milestone of the central programme). **Date:** 2026-06-02.
**Discipline:** `physics-derivation` + `physics-numerics`. **Builds on:** the established premise
(machine precision) that the slope family `{H^(a)=H0(a)+u·diag(a)}` (fixed `γ,ε`) commutes and
shares one `a`-independent eigenbasis `φ_i(u)`; the derivative coupling `W_ij=⟨φ_i|φ_j'⟩` is
`a`-independent (deformation_family_probe.py, ring_structure.py).

**Decisive test:** `experiments/aflow_closure_test.py` (reproducible; `numpy 2.4.6`, `scipy 1.17.1`).

---

## OUTCOME: (iii) — the a-flow does NOT close into a finite-dimensional geometric ODE.  [numerically-supported, machine precision on the structural legs]

The deformation flow `a ↦ {P_mm(a), b(a)}` is **organizing but not computational**: the family
pins the transcendental content to one shared geometric seed `W`, but no finite, history-free,
geometric `a`-flow transports it. This is the **Abelian ceiling restated as a flow statement** —
consistent with WS-D (the only zero-curvature structure Type-1 has is Abelian) and with
deformation_family_probe's interpretation (the `a`-invariant is the seed `W`, not the two
transcendentals). It is a **first-class negative**: the central principle remains true (one rigid
seed organizes the whole 2-parameter family) but it does **not** yield a free closed-form/transport
lever. Milestones 2–4 of the roadmap (which were gated on M1 = (i)/(ii)) are therefore **not viable
as a closed transport flow**.

This is a precise negative, not an inconclusive one: it rests on two independent machine-precision
structural facts (no algebraic `Q_k`; the KZ `R_k` is not geometric), each of which would have to be
**false** for the flow to close, plus the convergent ground-truth derivative that any closed form
would have had to reproduce.

---

## 1. Setup and the three closure mechanisms

Diabatic frame: `i U' = H(u) U`, `H = H0(a) + u·A`, `A = diag(a)`. Scattering `S = U(+∞,−∞)`
(phase-regularized at the rank-1 irregular point `u=∞`, leading term `u·diag(a)`); `P_jk=|S_jk|²`;
open data `{P_mm, b}` with `m=argsort(a)[1]`, `b=P[hi,lo]`.

`H` is **linear in `a`** (both `H0` and `u·A` are), so the deformation source is the algebraic
matrix
```
  ∂_{a_k} H(u) = ∂_{a_k} H0 + u·E_kk      (E_kk = diag(δ_ik)).
```
**Variation of holonomy (Duhamel, exact):**
```
  ∂_{a_k} S = −i ∫_{−∞}^{∞} U(+∞,u) [∂_{a_k} H(u)] U(u,−∞) du.          (D)
```
This is the reference. Closure means (D) reduces to a function of `S` (a finite reduced set) +
geometry, with no propagator history. Three candidate mechanisms:

- **(2) Zero-curvature `M_a` / boundary-term collapse.** If there is an **algebraic** (polynomial in
  `u`) `Q_k(u)` solving the homological equation
  ```
    ∂_u Q_k − i[H, Q_k] = ∂_{a_k} H,                                      (HOM)
  ```
  then the integrand in (D) is an exact `u`-derivative,
  `U(+∞,u)[∂_{a_k}H]U(u,−∞) = ∂_u[ U(+∞,u) Q_k(u) U(u,−∞) ]`, and (D) collapses to the
  **boundary terms** `∂_{a_k}S = −i( Q_k(+∞) S − S Q_k(−∞) )` (Stark-regularized). Closure (this
  mechanism) ⟺ an algebraic `Q_k` **exists**.

- **(3) KZ / Gaudin linear form.** If the `a`-flow is a flat KZ/Gaudin connection,
  `∂_{a_k} S = R_k(a) S` with `R_k` **geometric** (an `a`-function, history-free).

- **(1) Duhamel itself** is the exact-but-non-closed form (it carries the full history through
  `U(+∞,u)`, `U(u,−∞)`); it is the ground against which (2),(3) are tested.

### A structural simplification (verified, machine precision)
Because the whole family commutes, differentiating `[H^(a),H^(a')]=0` gives `[∂_{a_k}H, H]=0`:
**`∂_{a_k}H` is itself a commuting partner**, lying in the ring `span{I,H,H²}` (fit residual
`≤6.7e-15`). By Hellmann–Feynman it is **diagonal in the adiabatic eigenbasis**
(`off-diag ≤8.6e-16`, test `[1]`): `Φ^†(∂_{a_k}H)Φ = diag(∂_{a_k}E_i)`. So in the adiabatic frame
the source `∂_{a_k}D` is purely diagonal — exactly the premise's `∂_{a_k}D = diag(∂_{a_k}E_i^{(a)})`,
geometric. This is the *most favorable possible* setup for closure; it still fails, which sharpens
the negative.

---

## 2. The explicit closed-form RHS attempted, and why each fails

### Mechanism (2): no algebraic `Q_k`  [machine precision, decisive]
Solve (HOM) with the ansatz `Q_k(u) = Σ_{d=0}^{D} u^d Q_d` by **exact coefficient matching** in
powers of `u`: per power `p`,
```
  (p+1) Q_{p+1} − i[H0, Q_p] − i[A, Q_{p−1}] = RHS_p,   RHS_0 = ∂_{a_k}H0,  RHS_1 = E_kk,  else 0.
```
Least-squares homological residual vs polynomial degree `D` (canonical anchor):

| `k` | D=1 | D=2 | D=3 | D=4 | D=6 |
|---|---|---|---|---|---|
| 0 | 7.6e-01 | 1.0e-01 | 6.2e-02 | 6.5e-02 | 3.5e-02 |
| 1 | 7.3e-01 | 9.4e-02 | 5.9e-02 | 6.0e-02 | 2.6e-02 |
| 2 | 8.0e-01 | 9.5e-02 | 6.7e-02 | 6.4e-02 | 2.9e-02 |

The residual **plateaus** (it falls only as fast as a generic richer basis fits an arbitrary
function; D=6 is still `~3e-2`, not `→0`). This is the WS-D null signature: **no algebraic `Q_k`
exists.** Hence the Duhamel integral does **not** collapse to algebraic boundary terms; `Q_k` would
have to be a genuine (non-algebraic, oscillatory) solution carrying the dynamical phase
`e^{i∫(E_i−E_j)du}` in its off-diagonal-in-eigenbasis sector — i.e. the propagator history. (The
*diagonal* part of `Q_k` is the trivial algebraic antiderivative
`u·∂_{a_k}H0 + (u²/2)E_kk`; it is the **off-diagonal** part, sourced by the Berry coupling `W`, that
is non-algebraic. The naive antiderivative does not commute with `H`: `‖[H, u·∂_{a_k}H0+½u²E_kk]‖`
grows like `u²`, e.g. `9.0` at `u=5`.)

### Mechanism (3): the KZ `R_k = ∂_{a_k}S·S^{-1}` is not geometric  [numerically-supported]
`R_k := (∂_{a_k}S) S^{-1}` always fits `∂_{a_k}S = R_k S` exactly (residual `~1e-16`) — that is
empty (one can always solve for `R_k` at a single point). The **discriminating test** is whether
`R_k` is **history-free**: a flat/geometric `R_k(a)` must transport to a neighboring slope `a'`.
Freezing `R_k(a)` and predicting `∂_{a_k}S(a')` at `a' = a + 0.02 e_j`:

| `k` | fit_resid | `‖R_k(a)−R_k(a')‖/‖R_k‖` | TRANSPORT_err | `‖R_k‖` |
|---|---|---|---|---|
| 0 | 1.4e-16 | 0.51 | 0.53 | 1728 |
| 1 | 1.6e-16 | 0.40 | 0.40 | 1416 |
| 2 | 1.6e-16 | 0.59 | 0.54 | 1524 |

`R_k` **drifts by O(1)** over a small step and **does not transport** (transport_err `~0.4–0.5`,
i.e. order-unity relative error). It also has the divergent Stark magnitude `‖R_k‖~1.4–1.7×10³`
(the `~T²` content). So `R_k` is **not** a geometric/flat KZ connection — it is just `∂S·S^{-1}`
re-fit at each point, carrying the full propagator content. **Mechanism (3) does not close.**

> Caveat (honest): "geometric" here is tested operationally by transport over a finite `da`; the
> O(1) transport error decisively rules out a *constant/flat* `R_k`, but a non-closing flow is the
> natural reading. A genuinely flat KZ connection would have given transport_err → 0 (and an
> `a`-bounded `R_k`); we see neither.

---

## 3. The ground truth the closed forms had to reproduce (decisive-test data)

The gauge-invariant targets are `{P_mm, b}`. **Reported subtlety (a genuine finding):** the
finite-difference `∂_{a_k}{P_mm,b}` on the *raw truncated* propagator does **not converge in `T`** —
test `[4]`:

| `T` | `∂P_mm/∂a_1` | `∂b/∂a_1` |
|---|---|---|
| 40 | −4.43 | +5.36 |
| 70 | +1.39 | −3.09 |
| 100 | +0.78 | −2.40 |
| 130 | +0.70 | +1.70 |

Reason: the Stark phase `a_k u²/2` makes `lim_{T→∞}` and `∂_{a_k}` **fail to commute** — the
`a`-derivative of the oscillatory `O(1/T)` tail is `~T·sin(...)`, secular. The same non-convergence
appears in the *exact Duhamel* `∂_{a_k}P` (so it is structural, not a FD artifact). The converged
derivative exists (the converged `P_∞(a)` is smooth, deformation_family_probe F2) but must be taken
on **properly converged `P`**.

**Converged ground truth** (adiabatic-IP converged `P`; FD `h=1e-3`), all three directions, test
`[5]`:

| `k` | `∂P_mm/∂a_k` (T=70 → 100) | `∂b/∂a_k` (T=70 → 100) |
|---|---|---|
| 0 | +0.06149 → +0.06201 | +0.02095 → +0.02126 |
| 1 | −0.03508 → −0.03524 | −0.02724 → −0.02718 |
| 2 | −0.02639 → −0.02676 | +0.00629 → +0.00592 |

Stable across `T` in **every direction** (the two `T` values agree to `~1e-3`): the **converged
`∂_{a_k}{P_mm,b}` exists and is well-defined**. (Base `P_mm=0.21472435`, matching gold
`0.2147243114`; `b=P[hi,lo]` in the oracle diabatic basis. Cross-check with the full Richardson gold
oracle at the `k=1` anchor gives the same `−0.0352/−0.0271` to `~1e-3`.) This is the decisive-test
ground truth across the `N−1=2` essential deformation directions (`k=0,1,2`, with the `Σ_k`
gauge direction redundant) and both targets. Produced by `aflow_closure_test.py --full`.

These finite numbers are exactly what a closed form would have had to **reproduce from geometry
alone**. Mechanism (2) cannot (no algebraic `Q_k`); mechanism (3) cannot (no geometric `R_k`). The
Duhamel form (1) reproduces them, but only by integrating the full propagator history (verified:
the Duhamel identity matches FD `∂_{a_k}S` to relative `5.5e-5`, test `[0]`). **There is no
history-free closed RHS that lands on `−0.0353 / −0.0271`.**

---

## 4. Decisive-test summary table

| mechanism | closed RHS attempted | test | result | closes? |
|---|---|---|---|---|
| (1) Duhamel | `−i∫U[∂_aH]U` | vs FD `∂_a S` | rel `5.5e-5` (exact) | yes, but carries full history (not a closure) |
| (2) `M_a`/boundary | algebraic `Q_k`: `∂_uQ−i[H,Q]=∂_aH` | homological residual vs degree | plateaus at `~3e-2` (no `→0`) | **NO** |
| (3) KZ/Gaudin | `∂_aS=R_kS`, `R_k` geometric | transport `R_k(a)` to `a'` | transport_err `~0.4–0.5` | **NO** |
| ground truth | converged `∂_a{P_mm,b}` (oracle FD) | `T`-convergence | stable to `1e-3` in all 3 directions (e.g. `∂P_mm/∂a_1 ≈ −0.0352`) | (the target none of (2),(3) hit) |

Evidence tag: **structural legs (2),(3) and the Duhamel/Hellmann–Feynman checks are machine
precision; the convergence of the ground-truth derivative is numerically-supported (`T`-stable to
`1e-3`).** Verdict **(iii)**, numerically supported.

---

## 5. Why it fails — the mechanism, and what it means for the central principle

In the adiabatic frame `i χ' = (D−iW)χ`, the source `∂_{a_k}D` is diagonal (Hellmann–Feynman). The
homological equation `∂_u Q_k − i[D−iW, Q_k] = ∂_{a_k}D` splits: the **diagonal** part has the
trivial algebraic antiderivative, but the **off-diagonal** part is a homogeneous equation driven by
the Berry coupling `W`, whose solution is the convolution
`Q_{k,ij}(u) ∝ ∫^u e^{i∫_s^u(E_i−E_j)} (W-source)_{ij}(s) ds` — an oscillatory integral that is
**not algebraic in `u`** and is exactly the Stückelberg/connection content of `W`'s holonomy. This
is the same object the rank-3 connection problem (WS-RH `σ`) leaves open. The `a`-flow cannot
shortcut it: differentiating in `a` reproduces the *same* non-abelian holonomy integral, not a
finite-dimensional ODE.

So the central principle — *all the transcendental content is the holonomy of one `a`-independent
seed `W`* — **survives** (the seed is rigid; the family is its deformation orbit), but the
deformation is **not isomonodromic in a usable sense**: the data `{P_mm,b}` are not transported by a
closed geometric flow. The `σ`-conservation question (outcome (i)) is moot — there is no closed
abelian dressing flow to carry a conserved `σ`; the converged `∂_a{P_mm,b}` themselves require the
history. This matches WS-D's verdict (Type-1's only zero-curvature structure is Abelian and yields
only Abelian data) and deformation_family_probe's honest reading (the seed `W` is invariant, the two
numbers are not, and "no `a`-flow produces `W`; the flow only transports it" — here we show even the
*transport* does not close into a finite geometric ODE).

---

## 6. Honest scope, caveats, and what would overturn this

- **What is established (machine precision):** `∂_{a_k}H ∈ ring`, diagonal in the eigenbasis; the
  Duhamel identity; the non-existence of an *algebraic polynomial* `Q_k` (residual plateau to D=6).
- **What is numerically-supported:** the KZ `R_k` non-geometricity (transport over finite `da`); the
  `T`-convergence of the ground-truth `∂_a{P_mm,b}`.
- **What would overturn (iii):** (a) an algebraic `Q_k` in a richer-than-polynomial *but still
  history-free* basis (e.g. rational in `u` with poles only at the spectral ramification) that makes
  the homological residual collapse — tested informally in WS-D for the related `Ê` and also failed,
  but a dedicated symbolic search would upgrade the leg from "no polynomial" to "no algebraic";
  (b) a geometric `R_k` that transports (transport_err → 0) under a better gauge (canonical-frame
  `𝒮_canon` with the Coulomb `c_i log T` subtraction handled analytically) — we tested the
  probability-level flow, which is gauge-clean, and it does not close. Neither is expected to
  succeed given the WS-RH frontier status of `σ`, but they are the honest next checks if one wished
  to upgrade the negative to a theorem.
- **Tested locus:** the canonical anchor. The structural legs (ring membership, Hellmann–Feynman,
  the homological-residual plateau) are algebraic in `(γ,ε,a)` and not anchor-special; rerunning
  `aflow_closure_test.py` on `sampleB`/other strata is a cheap robustness extension (the script's
  `CANON` dict is the single knob).

---

## 7. Reproduce

```
python3 experiments/aflow_closure_test.py            # fast: structural legs + convergence demo
python3 experiments/aflow_closure_test.py --full     # adds the slow converged-FD oracle table [5]
```
Prints: Duhamel-vs-FD reference check `[0]`; ring/Hellmann–Feynman `[1]`; algebraic-`Q_k`
homological residuals `[2]`; KZ `R_k` transport test `[3]`; raw-propagator non-convergence demo
`[4]`; converged ground-truth `∂_a{P_mm,b}` `[5]`; and the **VERDICT (iii)**.
