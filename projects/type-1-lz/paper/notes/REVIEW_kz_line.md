# Deep-verification review: the isomonodromy/KZ/Kampé de Fériet line

> **UPDATE (gate-test #1 run, `gate_test_genus.md`): this review's A4 verdict is itself corrected.**
> A4 below refuted "Kampé de Fériet" by arguing the curve is genus-1 elliptic. Gate-test #1 shows
> that argument used the *wrong* curve (the phase double-cover `μ²=Q₄`, genus 1) instead of the
> *spectral* curve `Σ: χ_H=0`, which is **genus-0 rational** (global Gaudin parametrization
> `(E,u)=(m/p,n/p)`, verified) with a **structural node = the exact crossing**. Net: the function
> class is **NOT elliptic/Painlevé**; it is **genus-0 confluent-Heun / Kampé de Fériet** — so the
> Kampé de Fériet *lineage* is reinstated (though genus-0 ≠ elementary; confluent-Heun connection
> coefficients are themselves generically hard). The other findings (A2 KZ-Fuchsian-reduction
> unproven; A3 product-form integral circular for generic N=3; A6 confluence uncontrolled) **stand**.


**Reviewer role:** tough-but-fair peer review, *before* committing to theory-building.
**Hypothesis under review (H-KZ):** "The Type-1 N=3 transition amplitudes are connection/Stokes
coefficients of a KZ/Gaudin flat connection; they admit an explicit Euler integral representation
`ψ(u) ≈ ∫_Γ exp(iuℓ(λ)) ∏_i(λ−ε_i)^{α_i} R(λ) dλ` with exponents `α_i ∝ Γ_ij` and saddles at the
`Q₄` windows; the N=3 value is a **Kampé de Fériet** function of two window-actions; `P=|S|²`, with
BE and special slices as elementary degenerations."

Method: decompose into physical/mathematical assumptions, decontextualize each, check
independently, and label every invalidating element **fundamental** (kills it) or **non-fundamental**
(fixable later).

---

## Assumption ledger

### A1 — `P` is the Stokes data of the linear ODE with a rank-2 irregular point at `u=∞`.
- A1a the ODE has an irregular singularity at `∞` (degree-1 coefficient ⇒ Poincaré rank 2): **TRUE,
  established** (N=2 is Weber's equation exactly).
- A1b `P=|connection/Stokes data|²`: **TRUE, standard** for LZ as ODE asymptotics.
- **Verdict: SOUND, not novel.** This framing is correct and is the solid backbone. No flaw.

### A2 — The relevant connection is a **Fuchsian KZ/Gaudin** connection (poles at `ε_i`, residues from `γ`).
- A2a our couplings have Gaudin/Cauchy form: **plausible, unverified** that the residues satisfy the
  KZ `Ω_ij` (split-Casimir) structure with the braid relations. The commuting-family check we did
  (`[H^a,H^b]=0`) is necessary but **not** the same as flatness of a `λ`-connection.
- A2b **CONFLATION (fundamental).** Two different connections are being merged: (i) the connection
  in *time* `u` whose Stokes data is `P` — **irregular** at `∞`; (ii) a KZ connection in the
  *spectral/parameter* variable — **Fuchsian** at the `ε_i`. Equating `P` with KZ connection
  coefficients needs a specific transform/duality. A naive Laplace transform of `iψ'=(H₀+uA)ψ`
  gives `Aφ'(λ)=(H₀−iλ)φ`, whose coefficient is **linear in `λ`** ⇒ *still irregular at `λ=∞`*, **not
  Fuchsian.** So the reduction to a Fuchsian KZ system is **not established for the linear model.**
- **Verdict: fundamental gap.** The KZ identification is borrowed from BBGY's *hyperbolic* model;
  it does not transfer automatically. Revise to "isomonodromy connection" (true) rather than "KZ"
  (unproven here).

### A3 — Explicit **product-form Euler integral** with `α_i ∝ Γ_ij`, saddles `=` `Q₄` windows.
- A3a existence of a closed contour representation for *generic* N=3: **WEAKEST LINK.** Demkov–Osherov,
  bow-tie, and Lin–Sinitsyn close *because of special structure* (one moving level / a single
  crossing point / a 2nd-order reduction). **Simulation-review failure mode:** Laplace-transforming
  the N=3 system leaves a **3rd-order** ODE in `λ`; its solution is generically *not*
  `∏(λ−ε_i)^{α_i}`. The product-form integrand is the *answer* to that 3rd-order equation — so
  writing it down is **circular** unless the equation degenerates (2nd-order/hypergeometric), which
  is exactly the Lin–Sinitsyn `ε₀=0,g₁₃=0` locus we already showed is *outside* Type-1.
- A3b `α_i ∝ Γ_ij`: **guess by analogy**, unverified.
- A3c saddles `=` `Q₄` windows: a **check, not a given**; the saddles of `iuℓ+Σα_i/(λ−ε_i)` need not
  coincide with the gap-closing turning points unless the representation is the correct one.
- **Verdict: fundamental gap for generic Type-1.** The "explicit" representation is only explicit on
  a degenerate (hypergeometric) sub-locus; for the generic 3rd-order case it is unproven and likely
  not product-form.

### A4 — The N=3 value is a **Kampé de Fériet** (two-variable hypergeometric) function.
- **REFUTED AS STATED (fundamental for this specific claim).** Kampé de Fériet / `₂F₁`-type
  functions live on **genus-0 (rational)** curves. Our phase curve `μ²=Q₄` is **genus-1 elliptic**
  (verified: `Q₄` has 4 distinct roots, two conjugate pairs, for every sample). The assay itself
  treats `I_X` as *periods of the elliptic curve `μ²=Q₄`*. An elliptic connection problem points to
  **Heun / confluent-Heun functions or elliptic/theta functions**, *not* Kampé de Fériet. BBGY's
  Kampé de Fériet arose from the *hyperbolic* model, whose relevant curve is rational (genus 0);
  the Cauchy/linear Type-1 model has a *genuinely elliptic* phase curve, so the transcendental class
  is different and higher.
- Worse case (also live): the connection coefficient of an **irregular** isomonodromy problem need
  not be a *classical* named function at all — it can be a **Painlevé-type isomonodromic
  τ-function**, i.e. a genuinely transcendental object beyond Heun.
- **Verdict: the specific "Kampé de Fériet" identification is unjustified and probably wrong.** The
  honest function-class ladder is: elementary (degenerate loci) → `₂F₁`/Kampé (genus-0 sub-loci) →
  **Heun/elliptic (genus-1, the generic case here)** → Painlevé τ (generic irregular isomonodromy).

### A5 — Elementary degenerations (BE, slices) recovered.
- BE recovery via residue-collapse of the closed window periods: **TRUE, established in the assay** —
  but this only validates the **diagonal/Abelian** part. **Non-fundamental:** does not test the hard
  off-diagonal object.

### A6 — Confluence from BBGY's hyperbolic `Â+B̂/t` to the linear `H₀+uA` is controlled.
- **Fundamental technical gap (already flagged).** Confluence merges singularities and can change
  the Stokes structure discontinuously; results do not transfer for free.

---

## Cross-checks

- **Observation review (supportive, with a twist):** an *elliptic* (genus-1) structure naturally
  produces the **multiple oscillations** Lin–Sinitsyn observed and the overlap-driven enhancement we
  measured (2.5–104×) — via theta/elliptic functions near colliding periods. This actually
  *strengthens* the Heun/elliptic reading over the rational Kampé de Fériet reading.
- **Self-consistency with the ring result (good):** "fixed connection `W` + a-linear phases" is
  intact; the review only changes *which* transcendental class the holonomy of `W` belongs to.

---

## Verdict

- **Solid core — KEEP (evidence ladder: analytically-derived/established):** (i) `P` is a
  connection/Stokes coefficient of the rank-2 irregular ODE; (ii) the C–S `τ`-deformations are
  *isomonodromic* (preserve `P`); (iii) the commuting partner yields only Abelian data (the proven
  ceiling); (iv) the answer is "closed form in special functions, not elementary" *at the level of
  existence/organization*. This is a correct and useful framework.
- **Specific computational claims — REVISE / partially REFUTE (evidence ladder: → conjecture, one
  sub-claim → refuted):** "KZ Fuchsian reduction" (A2), "explicit product-form Euler integral" (A3),
  and especially "**Kampé de Fériet**" (A4) are **not justified** for the linear/elliptic Type-1
  model. A4-as-stated is **refuted** by the genus-1 finding; the live class is **Heun/elliptic or
  Painlevé**.

**Do NOT proceed to theory-building around "Kampé de Fériet."** The framework (isomonodromy +
Abelian ceiling + connection-coefficient) is ready to build on; the *function-class and
representation* are not, and building on the wrong class would be the doomed-weeks scenario this
review exists to prevent.

## Decisive, cheap next tests (gate theory-building on these)

1. **Pin the genus/monodromy.** Determine the curve actually governing the **off-diagonal** Stokes
   data (the discriminant locus of `χ_H`, i.e. the sextic gap curve), compute its genus and the
   order/singularity-type of its Picard–Fuchs system. Genus-0 ⇒ hypergeometric/Kampé possible;
   genus-1 ⇒ Heun/elliptic; irregular+nonclassical ⇒ Painlevé τ. **This single computation fixes the
   target class.**
2. **Test the representation, don't assume it.** Laplace-transform the N=3 system symbolically; check
   whether the resulting 3rd-order `λ`-ODE factorizes / reduces to 2nd order on the benchmark strata
   (it must, for a product-form integral to exist). If it stays irreducibly 3rd-order, the explicit
   Euler representation does not exist and the object is Heun/Painlevé.
3. **Saddle ↔ window check.** If a candidate integrand is proposed, verify numerically that its
   saddles coincide with the `Q₄` turning points before trusting it.

## Meta-review note (recurring error pattern)

Same class of error as the earlier falsified product-ansatze: **committing to a specific
closed-form class by analogy to an adjacent solved model without first computing the
monodromy/genus of *our* curve.** Add to the standing checklist: *"Before naming the
special-function class, compute the genus and singularity structure of the governing curve for the
actual model — do not import it from a neighbor."*
