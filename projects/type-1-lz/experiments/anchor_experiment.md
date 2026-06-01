# Day 1–2 anchor experiment: design, the Type-1 vs Chernyak–Sinitsyn mismatch, and first results

**Question that forced this note:** the strategy memo proposed reproducing the Chernyak–Sinitsyn
(C–S, arXiv:2006.15144) exact `ε₀=0` result as a validation anchor. But C–S Eq. 11 is the
**general** 3-state LZ Hamiltonian, not our Type-1 Cauchy family. *Do their results even apply to
our matrices?* This note resolves that before running anything, then runs the corrected experiment.

Reproduce with `python3 anchor_experiment.py` (needs numpy + scipy).

---

## 1. The dictionary: Type-1 → C–S frame (it exists, but lands outside their solvable locus)

C–S Eq. 11 is `H₃ = diag(b₁t, ε, −b₂t) + {g_ij}`: two **pure-slope outer** levels and one
**flat middle** level. Our Type-1 is `H(u) = H₀ + u·diag(a)` with the Cauchy coupling
`(H₀)_ij = γ_iγ_j(a_i−a_j)/(ε_i−ε_j)`.

**Mapping.** Order levels by slope `a_lo<a_mid<a_hi`. A time shift `u→u+t₀` plus a diagonal gauge
`(a_mid·u + c)·I` (neither changes transition probabilities) flattens the middle level and zeroes
the two outer constant offsets. Solving the two offset conditions:
```
t₀ = (H₀[hi,hi] − H₀[lo,lo]) / (a_lo − a_hi),   c = a_lo·t₀ + H₀[lo,lo]
ε_CS = a_mid·t₀ + H₀[mid,mid] − c           (the C–S middle energy)
b₁ = a_hi − a_mid,  b₂ = a_mid − a_lo        (outer slopes, both > 0)
```
The three couplings are the gauge-invariant off-diagonal entries of `H₀`, in particular the
**outer–outer coupling** `g_OUTER = (H₀)[lo,hi] = γ_loγ_hi(a_lo−a_hi)/(ε_lo−ε_hi)`.

**The crucial obstruction.** C–S's exactly-solvable point (confluent hypergeometric, their Eq. 49)
requires **both** `ε_CS = 0` **and** `g_OUTER = 0`. In Type-1:
- `ε_CS = 0` is reachable — one algebraic condition on `{γ,ε,a}` (numerically: scanning `a_mid`
  drives `ε_CS` through 0 near `a_mid≈1.997` for the canonical `ε,γ`).
- `g_OUTER = 0` is **never** reachable — it vanishes only if some `γ=0` (forbidden) or two slopes
  coincide (forbidden). Computed values: canonical `+0.72`, sampleB `+0.576`, sampleC `+0.936`.

> **Verdict.** The C–S → Type-1 dictionary is exact and useful, but the C–S *solvable* slice (the
> bow-tie point with zero outer–outer coupling) is **outside** the Type-1 family. So **C–S Eq. 49
> is not a valid anchor for Type-1.** Any anchor must be built from quantities intrinsic to Type-1.
> This is the corrected design below.

What *does* transfer from C–S: (i) the **method** — the Dykhne complex-turning-point semiclassics,
which is exactly our `Q₄`-window construction and applies to any spectrum; (ii) the **organizing
picture** — reduce the 3-state to combinations of two-level factors; (iii) the **warning** — for
`ε₀≠0` (and here always `g_OUTER≠0`) there is "no general analog of the Dykhne formula," i.e. the
prefactor is the hard, open part.

---

## 2. Corrected experiment (intrinsic to Type-1)

We never invoke an out-of-family exact result. Instead we self-calibrate on the one thing that
**is** exact for Type-1 — the Brundobler–Elser (BE) extreme-level survivals — and then measure the
open middle survival against the incoherent-product baseline.

**Ground truth.** Integrate `i U'(u) = H(u) U`, `U(−T)=I`, over `[−T,T]` (DOP853, rtol 1e-12).
`P[n←m] = |U_nm|²` (diabatic = adiabatic asymptotically). Rows sum to 1 to machine precision.

**Calibration (validates the integrator, no free parameters).** BE gives the two extreme-slope
survivals exactly:
```
P[lo survives] = exp(−2π(Γ_{lo,mid}+Γ_{lo,hi})),   P[hi survives] = exp(−2π(Γ_{hi,mid}+Γ_{hi,lo}))
Γ_ij = γ_i²γ_j²|a_i−a_j| / (ε_i−ε_j)²
```
The benchmark must reproduce these as `T→∞`.

**The open quantity.** The middle-slope survival `P_mid` is what BE/the windows cannot see. The
incoherent-grid baseline (independent sequential LZ) predicts
`P_mid^inc = exp(−2π(Γ_{mid,lo}+Γ_{mid,hi}))`. The deviation `P_mid / P_mid^inc` is precisely the
interference + multistate-Dykhne-prefactor content — the target of the whole programme.

**Convergence.** Survivals approach the asymptotic values with the characteristic oscillatory
`~1/T` Stückelberg tail (BE error oscillates ~1e-3 at T=40–160; `P_mid` stable to ~1%). The
project's interaction-picture harness removes the fast diabatic phases and reaches 1e-9 — that is
the oracle for the real runs; here T-convergence only needs to certify the qualitative result.

---

## 3. First results

| sample | P_lo (BE) | P_hi (BE) | **P_mid (benchmark)** | P_mid^inc (incoherent) | **ratio** |
|---|---|---|---|---|---|
| canonical | 0.0746 (0.0747) | 0.129 (0.129) | **0.212** | 0.0843 | **2.5×** |
| sampleB | 0.00074 (0.00040) | 0.056 (0.050) | **0.0169** | 0.000163 | **104×** |

**The headline.** The middle survival is **not** a small correction to the incoherent product —
it is **2.5× to 100× larger**, growing sharply as the two avoided crossings overlap (sampleB is the
strongly-overlapping regime). So the "prefactor/interference" we set out to *measure as a small
remainder* is in fact the **dominant** contribution to `P_mid`. The incoherent grid doesn't just
mispredict the middle survival by a little; it underestimates it by orders of magnitude.

**Physical reading.** Coherence strongly *protects* middle-level population relative to sequential
independent passages — constructive return of amplitude through the two crossings. This is the
Stückelberg interference between the two `Q₄` windows, now quantified.

---

## 4. Strategic implications (updates the v2 ranking)

1. **S2 (adiabatic Dykhne leading exponent) cannot reach `P_mid`.** Since the prefactor dominates
   (factor 2.5–100), the leading semiclassical exponent is the *wrong order of magnitude* for the
   middle survival in the overlapping regime. S2's value is confined to the BE extreme entries
   (exact) and the deep-adiabatic limit; it is **not** a route to the middle survival. Demote its
   claim accordingly.
2. **This sharpens what the off-diagonal/middle object must be.** A factor of 100 enhancement that
   grows with crossing overlap is a genuine two-window interference amplitude — consistent with the
   literature's "no Dykhne analog" / Kampé-de-Fériet transcendentality, and *inconsistent* with any
   product-of-real-tunneling-factors closed form (re-confirms the falsified ansatz class).
3. **It raises the stakes for S1 (zero-curvature factorization).** If a zero-curvature `Ê` exists
   and the path separates the crossings, the *exact* `S=∏S_ij` would have to reproduce this 2.5–100×
   enhancement through the *phases* of the two-level factors (Stückelberg), not their moduli. The
   decisive S1 test is therefore: does the path-deformed product, with its interference phases
   intact, match `P_mid` — or does an irreducible vertex survive exactly where the enhancement is
   largest (sampleB)?

**Next concrete step (S1 attempt):** construct `Ê` from the time-quadratic commuting partner and
test whether the deformed two-level product reproduces the measured `P_mid` (2.5× canonical, 104×
sampleB). sampleB is now the sharpest discriminator in the suite — any candidate closed form must
hit `P_mid ≈ 0.017`, two orders above the incoherent baseline.
