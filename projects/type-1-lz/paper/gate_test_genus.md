# Gate-test #1 result: the spectral curve is genus-0 (rational) — function class corrected

**Purpose (from `REVIEW_kz_line.md`):** determine the genus and singularity structure of the curve
governing the off-diagonal Stokes data, to fix the special-function class *before* theory-building.
Reproduce: the sympy/numpy snippets in the session log; all results to machine precision.

**This result corrects BOTH earlier takes**: the optimistic "Kampé de Fériet" *and* the pessimistic
"genus-1 → Heun/elliptic refutation." The truth: **the spectral curve is genus-0 rational**, with a
**structural node = the exact crossing**.

---

## 1. Findings (all verified to ≤10⁻¹³)

**(a) The spectral curve `Σ: χ_H(E,u)=0` is rational (genus 0) — universally for Type-1.**
The Gaudin/Bethe structure gives a *global rational parametrization*
```
E(λ) = m(λ)/p(λ),    u(λ) = n(λ)/p(λ),     p=∏(λ−ε_i),  n=quad(γ²),  m=quad(a·γ²).
```
Verified `|m(λ)/p(λ) − eigval(H(n/p))| ≤ 9×10⁻¹⁶`. A curve with a global rational parametrization
**is** genus 0. This is structural (holds for all Type-1 params), not sample-specific.

**(b) Riemann–Hurwitz, done correctly, agrees — once the node is identified.**
`D(u)=Disc_E χ_H` is a sextic with root multiplicities **[1,1,1,1,2]**. The double root is **not** a
branch point: the monodromy around it is the **identity** `[0,1,2]` (vs a transposition `[0,2,1]`
around each simple root). So it is a **node** — two sheets *cross* without braiding. True branch
points: the **4 simple roots** (transpositions). Normalization: smooth 3:1 cover of `P¹` with 4
simple branch points ⇒ `2g−2 = 3(−2)+4 = −2 ⇒ g = 0`. Consistent with (a).

**(c) The node is the "exact crossing."** It sits at a **real** `u* ≈ −0.249` where two eigenvalues
are **exactly degenerate** (gap `= 0`, eigenvalues `[−1.995, −0.119, −1.995]`). Present with the same
structure in every sample ⇒ **structural for Type-1**. This is precisely the user's flagged "exact
crossing"; the four simple (complex) branch points are the "avoided crossings at the real parts of
the `Q₄` conjugate pairs." A real eigenvalue crossing of a Hermitian family is codimension-2
(non-generic) — Type-1's structure *forces* the line `H₀+uA` through a diabolical point (the
crossing pair's effective coupling vanishes at `u*`).

**(d) Reconciling the genus-1 `μ²=Q₄`.** The curve `μ²=Q₄` (genus 1, 4 distinct roots) is the
**WKB phase/action double-cover** — the surface on which `√Q₄` (the action one-form) is
single-valued — **not** the spectral curve `Σ`. Its contribution to the *exponents* is residue-type
(third-kind), which is why the BE exponents collapse to elementary residues. The *sheet/connection*
structure lives on `Σ` (genus 0). My earlier refutation conflated these two curves — that step was
wrong.

---

## 2. Corrected function-class verdict

- **NOT elliptic / NOT Painlevé.** The spectral curve is rational; the earlier "genus-1 ⇒ Heun/
  elliptic" inference is **withdrawn** (it used the phase double-cover, not `Σ`).
- **Genus-0 special-function class is reinstated.** Exponents are elementary (genus-0 periods are
  residues — matches BE exactly). The remaining (prefactor) transcendence is a genus-0
  *connection-coefficient* problem.
- **The class is set by the marked points, not the genus.** On `P¹_λ` the connection has: 3 finite
  poles `{ε_i}`, an **irregular** point at `∞` (rank 2, the LZ growth), and the **node** (exact
  crossing) as a distinguished marked point. Three-plus-irregular/marked points ⇒ **confluent-Heun /
  Kampé de Fériet class** (genus-0, multi-singular), **not** plain `₂F₁` and **not** elliptic.

## 3. Honest residual (does this revive an elementary closed form?)

Partly, with a caveat. Genus-0 is good news — it rules out the hardest (elliptic/Painlevé) outcomes
and is consistent with the BBGY Kampé de Fériet lineage. **But genus-0 ≠ elementary:** the connection
coefficients of confluent-Heun equations are themselves generically *not* known in closed form. So
the working expectation becomes:
- **Exponents:** elementary (residues), already in hand (BE).
- **Prefactor (middle survival / off-diagonals):** a genus-0 **confluent-Heun / Kampé de Fériet**
  connection coefficient — closed form in named special functions, *probably not elementary*, but a
  far more tractable and better-characterized object than an elliptic/Painlevé transcendent.
- **The exact crossing (node) is a genuine reduction handle:** as a marked point with *elementary
  (logarithmic) local connection data*, it may lower the effective number of accessory parameters —
  the concrete thing for theory-building to exploit (this is the rigorous form of the earlier "H-B
  exact-crossing → P₂→₂" idea).

## 4. Net for theory-building

Gate-test #1 **passes in the constructive direction**: the curve is rational, the class is genus-0
confluent-Heun/Kampé de Fériet (not elliptic/Painlevé), and the exact crossing is a structural node
that both *characterizes* Type-1 and offers a reduction lever. Theory-building can now proceed on a
correct footing, targeting:
1. the connection problem on `P¹_λ` with poles `{ε_i}`, irregular point at `∞`, and the node `u*`;
2. exploiting the node's elementary local data to reduce accessory parameters;
3. with the standing expectation of a confluent-Heun/Kampé de Fériet closed form for the prefactor,
   benchmarked against the `10⁻⁹` harness and the measured 2.5–104× enhancement.

Still open (gate-test #2): whether the 3rd-order `λ`-connection reduces (factorizes) on any Type-1
sub-locus — which would drop confluent-Heun → hypergeometric → elementary there.
