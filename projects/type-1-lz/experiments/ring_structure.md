# The Type-1 ring structure and the viability of the commuting-partner approach

**User's lever:** "Type-1 matrices are a ring — all time-quadratic commuting partners are just
linear-`u` commuting partners with a linear-in-`u` coefficient. This bears on viability."

This note verifies the claim and works out exactly what it does (and does not) buy. Reproduce with
`python3 ring_structure.py` (numpy only). All checks pass to machine precision.

---

## 1. The claim is correct (verified, ≤10⁻¹³)

The family `{H^(a)(u) = H₀(a) + u·diag(a)}` (fixed `γ,ε`; varying slope `a`) satisfies:
1. **Commuting ring.** `‖[H^(a)(u), H^(b)(u)]‖ ≤ 9×10⁻¹⁶` for all pairs and `u`. The map
   `a ↦ H^(a)(u)` is linear and its image is the **full 3-dim commutant** of `H(u)` (a regular
   3×3 matrix's commutant is `span{I,H,H²}`, dim 3). So the family *is* the maximal abelian
   algebra containing `H(u)`.
2. **Quadratic partner reduces.** Therefore `H(u)² = Σ_k c_k(u)·H^(a_k)(u)` with reconstruction
   residual `8×10⁻¹⁴` and `c_k(u)` **linear in `u`** (`u²`-coefficients ~10⁻⁴ fit-noise). The
   time-quadratic commuting partner is *not a new object* — it is the linear family with
   linear-in-`u` scalar coefficients, exactly as claimed.
3. **Shared eigenbasis.** Commuting ⇒ simultaneously diagonalizable: the entire family shares the
   **same eigenvectors `φ_i(u)`** at every `u` (overlap defect ≤3×10⁻¹⁶). Members differ *only* in
   eigenvalues `E_i^(a)(u)`, which are **linear in `a`**.

---

## 2. Bearing on viability — the negative half (the Abelian ceiling)

A commuting partner carries no information beyond the shared spectrum and eigenbasis. Since the
quadratic partner is reducible, **going quadratic adds nothing**: the commutant is exhausted by
`{I, H, H²}`, all sharing `φ_i(u)`. There is no extra independent conserved quantity to exploit.

Concretely, what the commuting structure delivers is precisely the **Abelian data**:
- the eigenvalues `E_i(u)` → WKB dynamical phases `∫E_i du` → Dykhne exponents at the complex
  turning points `E_i=E_j` (= our `Q₄` windows) → the **Brundobler–Elser extreme survivals** (exact);
- the Chernyak–Sinitsyn **`τ`-invariance** (P depends only on the invariant combinations).

It does **not** deliver the non-adiabatic *mixing* — the Stückelberg interference that our anchor
experiment showed **dominates** the middle survival (2.5×–104× over the incoherent product).

> **Verdict for S2 (commuting-partner / constraint-closure route):** it has an **Abelian ceiling**.
> It can re-derive BE + the adiabatic-limit asymptotics + the `τ`-invariance constraints — all of
> which we already have — but it **cannot reach the middle-survival prefactor**. The ring structure
> *proves* this is not a matter of working harder: the needed information is simply not in the
> commutant. This explains why Chernyak–Sinitsyn obtained only asymptotic results, why there is
> "no general Dykhne analog," and why our prefactor is large.

---

## 3. Bearing on viability — the constructive half (what the ring *does* buy)

The shared eigenbasis cleanly **separates** the problem. In the adiabatic frame the Schrödinger
generator is
```
G^(a)(u) = diag(E_1^(a),E_2^(a),E_3^(a))(u)  −  i·W(u),     W_ij(u) = ⟨φ_i(u)|dφ_j(u)/du⟩,
```
and because the eigenbasis is shared, **`W(u)` is `a`-independent** — a *fixed* geometric
(Stückelberg) connection. All slope dependence sits in the diagonal, **linearly in `a`**. Hence
```
P^(a) = | holonomy of [ fixed off-diagonal W ] twisted by [ a-linear diagonal phases ] |².
```
This is a genuine simplification and the right target structure:
- The hard, non-Abelian object `W(u)` is computed **once**, independent of the slopes; the slope
  dependence is "trivial" (linear-in-`a` phases). This is the precise reason the phase one-form is
  linear in `L_H ∝ a` while the transport data is `a`-independent (consistent with the assay).
- `W(u)` is built from the **Cauchy eigenvectors**, so it is an **algebraic connection on the
  spectral curve**. The holonomy of an algebraic connection twisted by linear phases is an
  **isomonodromy / KZ-type problem** — whose answer is generically a special function (the
  Kampé de Fériet object of BBGY), **not elementary**. This independently reinforces the v2 verdict
  that no fully-elementary generic closed form should be expected.

---

## 4. Net consequence for the strategy ranking

- **S2 (time-quadratic commuting partner) is demoted from "constraint closure that might pin the
  middle survival" to "delivers only Abelian data (BE + adiabatic + `τ`-invariance)."** The ring
  structure rules it out as a route to the prefactor. Bank what it gives; stop expecting more.
- **The home run is NOT the commuting partner.** Malikis–Cheianov's path-deformation uses a
  **zero-curvature operator `Ê` that does *not* commute with `H`** (`[Ê,H]≠0`): a genuinely
  **non-Abelian** object, *outside* the ring and therefore *not* subject to this reduction. So the
  ring insight does not kill S1 — it **relocates the crux** to a sharp, well-posed question:

  > **Does Type-1 admit a non-commuting zero-curvature operator `Ê` that separates the crossings?**

  Unlike the bow-tie (an `su(2)`/`su(3)` spin representation, where `Ê` comes from the Lie
  structure), Type-1 is **Cauchy/Gaudin** data with no manifest spin-rep structure, so the
  existence of a useful non-Abelian `Ê` is genuinely open. That is the next thing to settle — and
  it is the *only* remaining route among the four papers that could yield an exact, elementary `P`.

- **Refined picture of the answer.** The ring structure recasts `P` as the non-Abelian holonomy of
  a *fixed algebraic connection* `W` twisted by *`a`-linear phases* — an isomonodromy/KZ problem.
  The middle-survival prefactor is the non-Abelian holonomy of `W`; absent a non-Abelian `Ê` that
  trivializes it, the honest expectation is a Kampé-de-Fériet / hypergeometric closed form, with
  the elementary content confined to BE and the `ε_CS=0`-type slices.
