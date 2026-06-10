# Derivation memo: the sectorial factorization of the Type-1 N=3 connection matrix

**Status: review gate.** This is the blind-derivation memo agreed in the redraft
interrogation (2026-06-10): the construction was derived with *no target punchline*; what
follows is what the derivation actually licenses. The rebuilt paper waits on your
reaction. Scripts: `experiments/ws_stokes_slots.py` (blind probe v1),
`ws_stokes_slots2.py` (convergence-hardened v2). Three parameter sets, two gold-gated.

---

## 0. The two motions, disaggregated (your directive, adopted throughout)

The memo — and the rebuilt paper — hold apart two notions my earlier slogan compressed:

| Motion | What it is | What it delivers | Independence witness |
|---|---|---|---|
| **Wildness** | analytic class of ∇ at u=∞ (irregular, rank 2) | the *architecture*: rays, sectors, the factorization, which slots exist, BE-pinning | N=2 LZ is wild yet its multiplier is **elementary** |
| **Transcendence** | differential-algebraic status of the slot *values* over Liouv(ε,γ,a) | the *hardness*: σ, b non-Liouvillian (trdeg 2, R24) | Heun is tame yet its connection data is **transcendental** |

Wildness answers **where**; non-rigidity (not wildness) answers **why hard**. Both get
their own theorem and emphasis in the rebuild; a four-box remark proves the axes
independent.

---

## 1. Analytic derivation: maximal ray-degeneracy (no numerics needed)

Rates: `q_j(u) = −i a_j u²/2 + …`; all pairwise differences `q_i − q_j = −i(a_i−a_j)u²/2`
have **pure-imaginary quadratic coefficients** (slopes real). On `u = R e^{iθ}`:

```
Re(q_i − q_j) = ½ (a_i − a_j) R² sin 2θ          — one function of θ for ALL pairs.
```

Consequences (each a small lemma in the rebuild):

1. **Shared singular geometry.** All three pairs have the *same* four oscillatory
   (anti-Stokes) directions — the **axes** θ = 0, π/2, π, 3π/2 — and the same four jump
   rays — the **diagonals** θ = ±π/4, ±3π/4. For generic irregular data the 3 pairs
   would have 12 distinct rays; Type-1's real commuting slopes collapse them 12 → 4.
   This **maximal ray-degeneracy** is the first structural fingerprint of the family.
2. **The physical directions are anti-Stokes.** The scattering frames live at θ = π, 0 —
   oscillatory directions. That is *why* P is well-defined and unistochastic: the
   physical observable sits exactly on the rays where no exponential dominates.
3. **Dominance alternation.** The dominance order of {e^{q_j}} is the **slope order** of
   a in quadrants 1, 3 and its **reversal** in quadrants 2, 4.
4. **Two crossings.** The upper-half-plane path from θ=π to θ=0 crosses exactly two jump
   rays: 3π/4, then π/4. Each carries a **full unipotent** Stokes factor (all three
   pairs jump together), the two oppositely triangular in slope order. Hence, in the
   slope-ordered formal frame, the connection matrix is predicted to factor as

   ```
   S  =  U · Δ · L                                  (upper-lateral sectorial product)
   ```
   U unit-upper (ray π/4), L unit-lower (ray 3π/4), Δ diagonal (formal transport through
   the θ=π/2 sector).

## 2. Numerical verdict (v2, hardened with the O(1/u) prefactor; R-Richardson)

Frame: Thome solutions `(I + T1/u)·diag exp{−i[a_j u²/2 + H₀_jj u + b1_j log u]}` with
`T1_kj = H₀_kj/(a_j−a_k)`, `b1_j = Σ_k H₀_jk²/(a_j−a_k)` (signed-BE exponents, Σb1=0);
physical normalization strips the incoming half-monodromy moduli `e^{πb1}` (log-branch
+iπ = upper continuation). Checks: row/col sums of |S|² = 1 to 5×10⁻⁶; `S S† = I` to
2×10⁻⁵; gold `P_mm`: canonical 0.2147314 (gold 0.214724), sampleB 0.0210250 (gold
0.021018).

**Finding (the amended H\*), uniform across all three cases:** the Gauss UDL
decomposition of S has diagonal

```
|Δ|  =  ( e^{+π|b1_lo|},  e^{−π|b1_mid|},  e^{−π|b1_hi|} )      — PURELY ELEMENTARY.
```

| case | extrap `|Δ|/formal` (lo, mid, hi) | dual-weight check on lo |
|---|---|---|
| canonical | (13.378, **1.00004**, **1.00003**) | 1/BE_lo² = 13.380 ✓ (1.5e-4) |
| sampleB | (2482.9, **0.99952**, **0.99955**) | 1/BE_lo² = 2480.4 ✓ (1e-3) |
| unsorted-a | (52.4199, **1.00002**, **0.99998**) | 1/BE_lo² = 52.4197 ✓ (4e-6) |

Notes on the table:
- The mid and hi entries are the **naive BE amplitudes**, exactly.
- The lo entry is **not** BE — it is the **inverse** weight `e^{+π|b1_lo|}`. This is
  *forced*: unitarity gives `|det S| = 1 = Π|Δ_j|`, and `Π e^{−π|b1_j|} ≠ 1`; with
  `Σ b1_j = 0` the unique elementary resolution is the dual weight on the last Gauss
  pivot. Still elementary — the hypothesis survives in amended form.
- **LDU does *not* have an elementary diagonal** (mid ratio ≈ 1.70 on canonical, stable).
  The asymmetry is meaningful: the log-branch (+iπ) of the frame selects the
  **upper-half-plane** continuation, and UDL — with exactly the triangularity pattern
  the ray-crossing order predicts — is the one that factors elementarily. The numerics
  thus *identify* UDL as the upper-lateral sectorial product, rather than assuming it.

## 3. The slot equations (exact at solver precision, all cases)

With S = U·Δ·L in slope order (0=lo, 1=mid, 2=hi):

```
P_hh   = |Δ_h|²                       = e^{−2π|b1_hi|}          (BE, multiplier-free)
P_ll   = |Δ_l + U01 Δ_m L10 + U02 Δ_h L20|² = e^{−2π|b1_lo|}    (BE as a forced
                                                                 CANCELLATION identity)
b      = P[hi→lo] = |L20|² · |Δ_h|²                             (ONE multiplier)
σ      = P_mm = |Δ_m + U12 Δ_h L21|²                            (TWO multipliers
                                                                 + one relative phase)
```

Read physically: σ is the interference of the **direct formal mid-passage** (amplitude
exactly the naive BE mid value `e^{−π|b1_m|}`) with the **one recombination path**
mid→hi→mid through the two ray crossings. The deviation of σ from its naive BE value *is*
the recombination term — an exact equation now, not a metaphor.

Three structural corollaries:

- **Disjoint slots = the structural face of trdeg 2 (R24).** b touches only `L20`;
  σ touches only `{U12, L21}` plus a relative phase. No identity of the factorization
  relates them — consistent with, and explanatory of, the measured functional
  independence of {σ, b}.
- **BE-lo as a dual identity.** In the upper-lateral product the hi-survival is manifest
  and the lo-survival is a nontrivial cancellation among large terms (|Δ_l| can be ~10³).
  The lower-lateral product mirrors this. Each lateral presentation makes one extreme
  manifest and encodes the other as a unitarity relation among multipliers.
- **Dimension count.** Free data: 6 complex multipliers + 3 phases of Δ = 15 real;
  unitarity of S = 9 real constraints; remainder **6 = dim of the wild character variety**
  (Theorem 2 of the proof doc, Boalch). The count matches — flagged as a consistency
  check, not a proof.

## 4. The licensed location statement (derived blind — this is what the data permits)

> **(i)** Wildness lives in the **irregular type A = diag(a)** — the torus rates
> `−i a_j u²/2`. It is the *architect*: it builds the four shared rays, the sector
> frames, the U·Δ·L factorization, and pins Δ to elementary half-monodromy weights.
> **(ii)** Within the factorization, **no transcendence resides in Δ** — all of it sits
> in the **six unipotent multipliers**.
> **(iii)** Within the observable P, the transcendence is carried by exactly two slots:
> σ = |Δ_m + U12 Δ_h L21|² and b = |L20 Δ_h|², with **disjoint multiplier dependence**.
> **(iv)** *Why* the multipliers are hard is **not wildness** — N=2 has the same wild
> architecture with elementary multipliers — it is **non-rigidity** (positive-dimensional
> wild character variety). Wildness locates; non-rigidity hardens.

Statement (i)–(iii): numerically-supported (3 cases, two gold-gated, errors ≤1e-3 and
scaling away under Richardson). Statement (iv): analytically-derived (Theorem 2 + the
N=2 witness). The Δ-formality in (ii) should be *provable* (the Schur-pivot identity
looks like a standard minor identity for lateral products) — appendix target, not yet
done.

## 5. Adjudication of the committed Proposition (one retraction)

The committed Proposition (claim (b), `prop:skeleton`) asserted *"the formal monodromy,
read on P, is a single oriented 3-cycle."* **This is wrong and will be retracted in the
rebuild.** At an unramified irregular point the formal monodromy is **diagonal** — here
`u^{−i b1_j}`-type factors whose half-turn moduli `e^{±π b1_j}` are *literally the BE
weights* (this memo's Δ). There is no permutation part. The oriented 3-cycle is real but
lives elsewhere: it is a property of the **incoherent skeleton product** (R16) and the
deep-overlap limit of P; the ℤ₂ orientation selector (T4) is its algebraic invariant.
The rebuild relocates the cycle + orientation to the skeleton subsection and assigns the
formal monodromy its correct (diagonal, BE-carrying) role. The remaining claims of the
Proposition (BE-pinning, affine reconstruction, slope grading, the LZ-product
realization) survive and are *strengthened* by the slot equations.

## 6. Consequences for the rebuild (per the agreed shape)

- §3 construction now has real content: rays/sectors lemma → U·Δ·L with elementary Δ
  (main text, structural) → slot equations → location statement §3 close. Heavy
  sectorial existence + the Δ-formality proof → appendix.
- §3.x contrast subsection gets sharper: the time-ordered LZ product M₃M₂M₁ (3 stochastic
  factors, u-time order, incoherent — *predicts σ,b falsely*) vs the ray-ordered U·Δ·L
  (2 unipotents + formal diagonal, dominance order, *exact*). Two factorizations, one
  skeleton, only one exact.
- The wildness/transcendence disaggregation (§0 table) becomes an introduction principle
  + a remark with the two witnesses.
- {σ, b} presented symmetrically (per the interrogation): two slot equations, two
  transcendentals, σ's distinction (hard in every regime) demoted to a remark.

**Open items for the appendix (flagged, not blocking):** (a) prove Δ-formality;
(b) log-normalize the phases of Δ (moduli are R-stable; phases carry the b1·log R
divergence in the current frame); (c) the lower-lateral/LDU branch bookkeeping;
(d) tighten sampleB's BE-lo gate (0.8% — deep suppression amplifies relative error).
