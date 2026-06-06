# Where this travels: Type-1 commuting matrices & their geometry, beyond Landau–Zener

**Type:** `physics-theory-building` (research overview / meta-review). **Date:** 2026-06-03.
**Mood of record:** written as the Type-1 N=3 *solvability* quest is set down. The honest frame: the
near-irreducible wildness is **a theorem, not a failure** — Type-1 N=3 is the *generic non-rigid* rank-3
connection, σ is its accessory parameter promoted to a wild Stokes constant, and `solvable ⟺ rigid`. This
note maps where the structure — and, more importantly, what *this project produced* — is useful elsewhere.
All application claims are `[conjecture]`/`[framing]` unless tagged otherwise; the value is the map.

---

## 0. Two kinds of "useful"

- **Exportable assets we *made*** (novel, more surprising than "Cauchy/Gaudin appears everywhere").
- **Domains the *object* already inhabits** (where Type-1 geometry illuminates a known problem).

## I. Transferable assets (what the project produced)

1. **`Solvable ⟺ rigid` as a design principle for driven multilevel systems** `[AD+NS; export=conjecture]`.
   A driven N-level system has closed-form amplitudes iff its connection is *rigid* (Katz index 2) iff it has
   an engineered **pseudo-reflection** (rank-1 coupling / coinciding slopes / a spectral degeneracy parked on
   a singular point). This **unifies the known solvable models** (bowtie, Demkov–Osherov, Demkov–Kunike) as
   "the rigid ones" and is *constructive*: it says how to build new solvable protocols and certifies which
   are provably not closed-form. Cleanest gift to quantum control. (Ref: `ws_geom_bundle_theory.md` §8c,
   `experiments/ws_rigid_3F2.py`.)

2. **σ as a Fredholm/Widom determinant ↔ random-matrix gap probabilities** `[EST(σ=Fredholm) + conjecture]`.
   σ closes as a Widom-class Fredholm determinant (R11b). Such determinants are RMT gap probabilities
   (Tracy–Widom; sine/Airy/Bessel kernels) and Painlevé τ-functions. Probe: match σ's kernel to a known
   integrable kernel → σ becomes computable by Painlevé/RMT machinery, and Type-1 yields a new solvable-kernel
   family.

3. **The wild-Hitchin / Riemann–Hilbert dictionary** `[framing]` as a reusable template for *where the
   transcendence lives* in any integrable-but-not-elementary problem: integrable (Hitchin/Gaudin) side
   elementary; isomonodromy (Stokes) side = one irreducible constant.

## II. Domains the object already inhabits (ranked: leverage × groundedness × novelty)

1. **Richardson–Gaudin integrable pairing** `[grounded]` (BCS / nuclear / cold-atom pairing — the
   Yuzbashyan lineage). Type-1 *is* the driven Richardson–Gaudin algebra. New lens to bring back: the
   **universal node (R3)** and the `[1,1,1,1,2]` discriminant geometry are the spectral data governing
   **pairing-quench dynamics** (dynamical gap vanishing, the steady-state phases). The node = a protected
   dynamical degeneracy; its location/migration may be a sharper quench-phase invariant.

2. **Argyres–Douglas / class-S gauge theory** `[deep; conjecture — needs grounding]`. Genus-0 spectral curve
   + rank-2 irregular puncture *is* an AD Seiberg–Witten curve. Then BE periods = central charges,
   **σ = a Stokes / wall-crossing datum** (line-defect VEV, BPS invariant), and `rigidity = solvability`
   becomes *which AD theories are Lagrangian / have closed-form Schur indices*. Note the coincidence with I.1:
   **"rigid = solvable" (control) = "AD theory is Lagrangian" (gauge)** — two languages, one statement.

3. **Holonomic quantum gates & reachability** `[medium; novel angle]`. The `a`-independent connection `W`
   with **marginal holonomy `Λ≈π`** is exactly the regime where non-abelian geometric (Berry/Wilczek–Zee)
   gates are nontrivial (large rotations). The **two-vertex reachable set (M2)** is a controllability
   statement: which unitaries a slope-sweep reaches with one fixed connection — a family of geometric gates
   sharing a connection.

4. **Conical intersections / diabolical points / topological matter** `[medium]`. R3's node is an
   integrability-protected diabolical point; the mechanism (degeneracy made robust by a commuting partner) is
   the symmetry-protected band-crossing / Jahn–Teller / photochemical-funnel mechanism. The `δ_j` monodromy =
   the Berry phase around the intersection.

5. **Non-Hermitian / exceptional points** `[medium]`. Complexify → the four complex branch points become
   exceptional points with Type-1-controlled geometry → structured EP sensing / PT-symmetric design.

6. **Structured linear algebra & tensor methods** `[lower, real]`. The shared `a`-independent eigenbasis is a
   **joint diagonalizer** (tensor CP / ICA / blind source separation); Cauchy displacement structure → fast
   structured solvers.

*Wild cards* `[speculative]`: rigid local systems ↔ hypergeometric *motives* / L-functions (a number-theoretic
shadow of `rigidity = solvability`); the two-vertex reachable set ↔ geometric control theory.

## III. Top picks and first probes (a quick tournament)

| rank | direction | first probe | payoff |
|---|---|---|---|
| 1 | **I.1 rigidity → control design principle** | write the standalone classifier note (finished physics) | reframes an applied field; needs no σ-heroics |
| 2 | **I.2 σ ↔ RMT/Painlevé** | match σ's Fredholm kernel to a known integrable kernel | could *compute* σ; new solvable kernels |
| 3 | **II.2 Argyres–Douglas** | literature pass: is the curve a known AD SW curve? | deepest; unifies with I.1 |

## IV. Recommendation

Develop **one**, not seven. The cleanest, most self-contained gift is **I.1** — "rigidity as the solvability
criterion for driven multilevel systems" — written for the quantum-control / exactly-solvable-dynamics
community; it is finished physics and unifies the solvable-MLZ zoo. The tempting technical follow-up is **I.2**
(σ ↔ RMT). The most beautiful long shot is **II.2** (Argyres–Douglas), pending a literature check that the SW
identification is exact.

*What this note preserves:* not a solution, but the **map of where a proven structure and a proven
classifier travel** — so the beauty that seduced the LZ quest is banked, not lost.
