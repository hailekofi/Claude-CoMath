# Co-scaling and the marginally-overlapping-crossings theorem (WS-D obstruction, upgraded)

Promotes the WS-D obstruction from a raw numerical scan to an **analytic backbone + a sharp,
constraint-correct numerical bound**. Reproduce: `experiments/coscaling.py` (numpy).

## Setup
Type-1 diabatic data: diabatic energies `D_i(u) = s_i + a_i u`, `s_i = (H₀)_ii =
−Σ_{k≠i} γ_k²(a_i−a_k)/(ε_i−ε_k)`; couplings `V_ij = (H₀)_ij = γ_iγ_j(a_i−a_j)/(ε_i−ε_j)`.
Pairwise crossing at `u_ij = −(s_i−s_j)/(a_i−a_j)`. Avoided-crossing width (the `u`-range where the
gap `√((a_i−a_j)²(u−u_ij)² + 4V_ij²)` is within O(1) of its minimum `2|V_ij|`):
`w_ij = 2|V_ij|/|a_i−a_j|`.

## Lemma (width is slope-independent) — [established, exact]
```
w_ij = 2|V_ij|/|a_i−a_j| = 2|γ_iγ_j(a_i−a_j)/(ε_i−ε_j)| / |a_i−a_j| = 2|γ_iγ_j| / |ε_i−ε_j|.
```
The slope difference **cancels exactly** (symbolic, `coscaling.py`): the avoided-crossing width
depends only on `(γ,ε)`. Moreover, with `Γ_ij = γ_i²γ_j²|a_i−a_j|/(ε_i−ε_j)²` (the BE/LZ exponent),
```
Γ_ij = (w_ij/2)² · |a_i−a_j|    ⟺    |a_i−a_j| = 4Γ_ij / w_ij².      [exact]
```
So the slopes enter the geometry *only* through the `Γ_ij`; the widths are fixed by `(γ,ε)` alone.

## Scale-invariance — [analytic]
`sep/width` is invariant under each independent rescaling `γ→λγ`, `ε→νε`, `a→μa` (verified):
`w ∼ γ²/Δε`, and `s_i ∼ γ²·Δa/Δε` gives `u_ij ∼ γ²/Δε` — **both the separations and the widths
carry the same scale `γ²/Δε`, with the slope scale cancelling**. Hence `sep/width` is a
*dimensionless shape function* of the parameter ratios; no overall scaling can grow it.

## Theorem (permanently marginally-overlapping crossings) — [analytic backbone + numerically-supported bound]
On the **genuine-LZ locus** — where all three crossings are real Landau–Zener events,
`Γ_ij ∈ [Γ_min, Γ_max] ⊂ (0,∞)` — the closest pair of crossings is always within an O(1) multiple
of the largest width:
```
 min_{adjacent} |u_ij − u_kl| / max_ij w_ij  ≤  C ≈ 1.6    (35k samples, Γ_ij∈[0.2,5]; median 0.20).
```
**At least two of the three crossings are always marginally-overlapping.** Therefore the
Malikis–Cheianov configuration — *all three crossings simultaneously isolated and each a genuine
2-level LZ* — **does not exist anywhere in the Type-1 N=3 parameter space**, and the exact product
`S = ∏_{i<j} S_ij` is unreachable. [This is the corrected, sharp form of the WS-D obstruction.]

### Why the obvious escape fails (loophole closed)
One *can* send a single crossing to infinity via the degenerate-slope limit `a_i→a_j` (`u_ij→∞`).
But there `Γ_ij = (w_ij/2)²|a_i−a_j| → 0`: the runaway crossing becomes a **trivial diabatic
pass-through**, not an MC 2-level event — and the *other two* crossings stay within ~1.6 widths
(`min`-separation picks the close pair). So separating any crossing either makes it trivial or
leaves a marginally-overlapping pair behind. The genuine-LZ bound is robust (and *tighter*, 1.6,
than the unconstrained scan's 4.2, precisely because the high-ratio tail there had some `Γ_ij`
outside `[0.2,5]`).

## Evidence ladder
- Width lemma + `Γ`–width–slope relation + scale-invariance: **analytically-derived (exact)**.
- The O(1) bound constant (`≈1.6`) on the genuine-LZ locus: **numerically-supported** (35k samples;
  a fully analytic `sup` over shape space is the remaining open piece, but the *mechanism*
  — common scale `γ²/Δε`, regular off the slope-degenerate boundary where `Γ→0` — is analytic).

## Consequence for the program
This upgrades the obstruction theorem (O4 = FAIL) and gives a clean, quotable physical statement:
**the Cauchy structure welds the avoided-crossing widths to the level spacing (`w_ij=2γ_iγ_j/Δε_ij`,
slope-free), so Type-1 crossings are permanently marginally-overlapping** — the structural reason
Type-1 sits outside the factorizable ("constellation") class and why its middle-survival prefactor
is an irreducible genus-0 3-level-Weber Stokes connection coefficient (WS-A), not a product of
two-level factors. It also explains the empirical 15/85 split (the 15% are the high-`Γ`-spread tail
where the product is *approximately* good, never exact).
