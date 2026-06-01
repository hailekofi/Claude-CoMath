# Nomenclature (authoritative) — Type-1 N=3 Landau–Zener

Canonical symbol conventions for the whole project. Established with the user 2026-06-01 to
remove a `Γ`/`S` collision. **`Γ` and `S` have exactly the meanings below and no others.**

## Spectral / Cauchy geometry
- `λ` : spectral parameter on the rational spectral curve; `u = n(λ)/p(λ) = Σ_k γ_k²/(λ−ε_k)`,
  `p(λ)=∏_k(λ−ε_k)`. The `λ_j(u)` are the three roots of `q_u = u·p − n`.
- **`Γ_j`** (single index, RESERVED): the canonical Cauchy form factor / sheet normalization,
  ```
  Γ_j = δ_j · ( Σ_k γ_k² / (λ_j − ε_k)² )^{-1/2} = δ_j · (−u'(λ_j))^{-1/2},   δ_j = ±1.
  ```
- **`S_{ij}`** (RESERVED): the canonical Cauchy object built from the form factors,
  ```
  S_{ij} = Γ_i Γ_j / (λ_i − λ_j).
  ```

## Couplings, widths- `s_{ij} := γ_iγ_j/(ε_i−ε_j)` (SIGNED Cauchy half-width) — deliberately lowercase, echoing
  `S_{ij}=Γ_iΓ_j/(λ_i−λ_j)` (with γ↔Γ, ε↔λ).
- `w_{ij} := |2 s_{ij}| = 2|γ_iγ_j|/|ε_i−ε_j|` — the (slope-free) avoided-crossing width.
- Off-diagonal Hamiltonian coupling `(H₀)_{ij} = s_{ij}(a_i−a_j) = γ_iγ_j(a_i−a_j)/(ε_i−ε_j)`.

## Dynamical quantities
- **Brundobler–Elser (BE) pairwise survival exponent** `= s_{ij}² |a_i−a_j|`
  ( `= γ_i²γ_j²|a_i−a_j|/(ε_i−ε_j)²` ). **NOT denoted `Γ`.** Written symbol-free as
  `s_{ij}²|a_i−a_j|`. Extreme-level survival `= ∏_{m} exp(−2π s_{nm}²|a_n−a_m|)`.
- **Canonical-frame Coulomb coefficient** `c_i = Σ_{j≠i} s_{ij}² (a_i−a_j)` (signed; `Σ_i c_i=0`).
- **`𝒮`** : the time-evolution **scattering matrix** (the transition amplitudes), `P_{x→j}=|𝒮_{xj}|²`.
  Regularized forms `𝒮_IP`, `𝒮_canon`. The MC-type factorization is `𝒮=∏_{i<j} 𝒮_{ij}`.
  In code (ASCII), the scattering matrix is named `Smat` (= `𝒮`); `s_ij`/`w_ij` as above.

## Channel ordering conventions (READ — three distinct orderings; do not conflate)
- **Basis / matrix index = ε-ordering (fixed).** States are labeled by index `i ↔ (γ_i,ε_i,a_i)` with
  `ε_0<ε_1<ε_2` (Params convention). The transition matrix `P[i,j]` (and `𝒮`, `H₀`, the spectral
  curve, the form factors `Γ_j` on `λ_j∈(ε_{j-1},ε_j)`) is **ε-indexed**. This is THE basis convention.
- **Extreme/middle = SLOPE-ordering.** Brundobler–Elser: the two **extreme-slope** levels
  (`argsort(a)[0]`,`argsort(a)[2]`) have exact survivals; the **middle-slope** level
  `m=argsort(a)[1]` is the OPEN one. The open middle survival is
  **`P₂→₂ ≡ P_mid := P[m,m]`, `m=argsort(a)[1]`** — the *middle-SLOPE* diagonal entry. It equals the
  literal adiabatic survival of the middle eigenstate (asymptotically adiabatic = diabatic).
  **It is NOT `P[1,1]` (ε-index-1) in general** — slope-middle = ε-middle *only* when the slopes are
  monotonic in ε (true for the `canonical` and `sampleB` benchmarks; false for generic samples). All
  middle-survival code identifies `m` via `argsort(a)` (verified consistent: anchor_experiment,
  num_S12, restart_probe, ws_ch).
- **Energy ordering = NOT a labeling convention.** The instantaneous energy rank *reverses* between
  `u=±∞` (it only swaps the two extremes; the middle is rank-2 at both ends). Do **not** use energy
  rank to label channels — use ε (basis) and slope (BE roles). [Retraction: an earlier framing of
  `P₂→₂` "by energy rank" was a red herring; the correct statement is slope-middle, ε-indexed.]

## Quick map (old → canonical)
| old / colliding | canonical |
|---|---|
| `Γ_{ij}` = BE exponent | `s_{ij}²|a_i−a_j|` (no `Γ`) |
| `Γ` (the Coulomb `c_i` build) | `c_i = Σ_j s_{ij}²(a_i−a_j)` |
| `S`, `S_IP`, `S_canon` (scattering) | `𝒮`, `𝒮_IP`, `𝒮_canon` (code: `Smat`) |
| `S=∏S_{ij}` (factorization) | `𝒮=∏𝒮_{ij}` |
| `Γ_j`, `S_{ij}=Γ_iΓ_j/(λ_i−λ_j)` | UNCHANGED — these are the canonical reserved symbols |

**Note:** files under `uploads/` are the user's original source materials and are left as-is;
this convention governs all session artifacts in `paper/`, `experiments/`, and `RESEARCH_LOG.md`.
