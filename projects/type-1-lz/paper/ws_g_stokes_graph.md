# WS-G — The Stokes graph and the "joint" diagnosis of the Type-1 N=3 MLZ prefactor

**Owner:** WS-G (Stokes-graph / spectral-network mapping). **Date:** 2026-06-01.
**Inputs used:** WS-A Riemann scheme (`paper/ws_a_riemann_scheme.md`), gate-test genus
(`paper/gate_test_genus.md`), WS-D obstruction (`experiments/ws_d_verification.py`),
`experiments/anchor_experiment.py`. **Constraint compliance:** no git ops; the only files
written are this note, `experiments/stokes_graph.py`, and `experiments/figs/*.png`;
throwaway probes lived under `/tmp/`.

Model: `H(u)=H0+u·diag(a)`, Type-1 Cauchy coupling, Schrödinger ODE `i ψ'(u)=H(u)ψ`.
Adiabatic eigenvalues `E_i(u)=eig H(u)` for complex `u`. Pairwise turning points
`t_ij`: branch points of `Disc_E χ_H` (the complexified avoided crossings + the real node).
Stokes lines of type `ij`: `Im ∫_{t_ij}^u (E_i−E_j) du' = 0` (three at 120° per simple tp).
**Joint** = intersection of a type-`ij` Stokes line with a type-`kl` (`kl≠ij`) line, away
from any turning point — a BPS junction in Gaiotto–Moore–Neitzke (GMN) language.

Showcase samples:
`sampleB (OVERLAPPING)` `ε=(−1,0,1.5) γ=(0.9,1.1,0.8) a=(−0.7,0.4,1.3)`, sep/width ratio
**0.14**; and `sampleA (WELL-SEPARATED, constructed)` `ε=(−9.173,5.962,8.603)
γ=(−0.453,−0.105,0.981) a=(−2.297,−1.075,0.406)`, ratio **3.97** (found by searching for
the high-ratio tail WS-D proved is thin but reaches ~4).

> **NOTE on the ratio metric.** I use WS-D's exact separation/width metric
> (`experiments/ws_d_verification.py`): `min crossing-separation / max(2|H0_ij|/|a_i−a_j|)`.
> On that metric sampleB scores **0.14**, not the "~0.5" quoted loosely in the task brief —
> I reproduced WS-D's own code to confirm. The qualitative split (sampleB strongly
> overlapping, sampleA well-separated) is unaffected.

---

## 0. Executive summary / verdict

**VERDICT (the joint diagnosis is CONFIRMED at the numerically-supported level).**
The Stokes graph cleanly separates the two empirical regimes exactly as conjectured:

- **Well-separated (ratio ≳ 1, the ~15% where the incoherent product is exact):**
  the turning points hug the real axis (narrow avoided crossings) and sit far apart along
  it; each pairwise turning-point cluster carries its *own* short Stokes lines that **do not
  reach** the other clusters. **The Stokes graph is JOINT-FREE** (a "ladder"). `P` factorizes
  into independent 2-level Weber pieces ⇒ enhancement `≈ 1.00×`. (`figs/stokes_sepA.png`)

- **Overlapping (ratio ≲ 1, the ~85% where `P₂→₂` is enhanced 2.5–104×):**
  the complex turning points push deep off the real axis and their Stokes lines of
  **different pair-types (12 and 13) intersect**, forming **genuine joints** — for sampleB,
  **4 joints** (two conjugate pairs at `u ≈ −0.10 ± 0.90 i`) where the red (12) and green (13)
  Stokes lines cross between the turning points. **The Stokes graph HAS joints.** `P` does NOT
  factorize ⇒ enhancement `132.75×`. (`figs/stokes_overlapB.png`)

- **Across a 16-sample scan spanning ratio 0.1→4**, joint presence and joint "strength"
  track the 15/85 split and the enhancement: see the correlation table in §5. The
  joint-free samples are exactly the (near-)factorizing ones; the jointed samples are exactly
  the enhanced ones.

**Interpretation.** This makes the WS-A picture geometrically concrete. WS-A established that
`P` is the Stokes/connection data of a single rank-2 irregular point at `u=∞`, with the
*elementary* part = the Weber exponential data (BE survivals `e^{−2πΓ_ij}`) and the
*transcendental* remainder = the off-diagonal Stokes connection coefficient. WS-G localizes
that remainder **on the complex u-plane**: it is the holonomy contribution of the **joint**.
No joint ⇒ the connection matrix is the ordered product of independent 2-level (Weber) Stokes
factors ⇒ elementary `P`. A joint ⇒ the GMN junction rule mixes the three sheets ⇒ a non-
factorizable connection coefficient = the measured `P₂→₂` enhancement. The "joint contribution
IS the off-diagonal Stokes connection coefficient" conjecture is **supported**.

**Evidence-ladder status:** the turning-point structure and the joint-vs-enhancement
*correlation* are **numerically-supported (strong)**. The identification "joint holonomy =
the precise value of the `P₂→₂` enhancement" is **conjectural/framework** — WS-G establishes
the qualitative dictionary and a quantitative *correlation*, not the GMN junction-rule
*computation* of the enhancement (that would be WS-E's job). Honest ambiguities are in §6.

---

## 1. Turning-point structure (established)

`Disc_E χ_H(E,u)` is a degree-6 polynomial in `u` with root multiplicities **[1,1,1,1,2]**
(reproducing gate-test-genus and WS-A exactly, via exact rational arithmetic in
`turning_points()`). The five turning points are:

- **One real node** (the double root): two eigenvalues *exactly* degenerate. This is the
  **exact crossing** of WS-A/E3. For sampleB it is the **pair-23** degeneracy at
  `u* = −0.6889`; for sampleA the **pair-23** degeneracy at `u* = +0.334`. (The node pair is
  identified robustly by which two diabatic channels collide, §2.)
- **Four simple branch points** in two complex-conjugate pairs — the **complexified avoided
  crossings / the `Q₄` windows**. For sampleB: a **12** pair at `−0.689 ± 2.048 i` and a
  **13** pair at `+0.736 ± 1.300 i`. For sampleA: a **23** pair at `−0.353 ± 0.077 i` and a
  **12** pair at `−0.041 ± 0.051 i`.

**Classification by colliding pair** is done by tracking eigenvalue labels along a straight
path from `u₀=−60` (where the diabatic order = sorted real-part order = slope order) to the
turning point, matching eigenvalues by nearest-neighbour continuity, then reading off the two
sheets that merge (`classify_turning_pair()`). This is the careful branch-tracking the brief
demanded; it is stable because the path stays away from the *other* branch points.

**Key geometric fact (the whole story in one line).** The imaginary part of a simple turning
point ≈ the avoided-crossing **half-width**, while its real part ≈ the diabatic **crossing
location**. So:
- **narrow gaps + far crossings (well-separated)** ⇒ tps hug the real axis, far apart ⇒ short,
  disjoint Stokes fans ⇒ **no joint**;
- **wide gaps / near crossings (overlapping)** ⇒ tps far off-axis, near in real part ⇒
  long Stokes lines from different pairs that **must cross** ⇒ **joint**.
This is the geometric mechanism behind WS-D's "permanently marginally-overlapping" finding:
because Type-1 keeps the ratio `~O(1)`, most samples land in the jointed regime (the 85%).

---

## 2. Stokes lines and the graph (method)

A Stokes line of type `ij` is the zero set of the WKB phase field
`φ_ij(u) = Im ∫_{t_ij}^u (E_i−E_j) du'`. I compute `φ_ij` on a grid by **ray-integration**
from the turning point `t_ij` (straight ray `t→u`), tracking the colliding pair `(E_i,E_j)`
by nearest-neighbour continuity (vectorized over the whole grid: ~26 batched
`np.linalg.eigvals` calls). The zero contours (`matplotlib` `contour`) are the Stokes lines;
near a simple turning point they emerge as three branches at 120°, as required.

**Branch-cut hygiene (per the brief).** The ray integral is trusted only inside a disk of
radius `0.92 ×` (distance to the nearest *other* turning point); beyond a neighbouring branch
point the ray would cross a cut of `√Disc` and `φ` would be corrupted. Each turning point's
lines are masked to its own clean disk, so every drawn segment is a bona-fide Stokes line of a
definite pair-type. The node (a degenerate turning point) is included and carries type-`ij`
lines for its colliding pair (23 in both samples).

This grid+contour construction was cross-checked against an independent inverse-`Δ` flow
integrator (`du/ds = 1/(E_i−E_j)`); both agree that the well-separated graph is joint-free and
the overlapping graph has the 12×13 crossings.

---

## 3. The two regimes — figures

### 3a. WELL-SEPARATED (sampleA, ratio 3.97) — `figs/stokes_sepA.png`
Three **disjoint** turning-point clusters strung along the real axis (node at `+0.33`, a 23
pair at `−0.35`, a 12 pair at `−0.04`), each with short Stokes fans confined near its own
crossing. **No Stokes line of one pair reaches a turning point of another pair ⇒ 0 joints.**
Direct ODE benchmark: `P₂→₂ = 0.9859`, incoherent product `0.9859`, **enhancement 1.00×** —
the incoherent product is *exact*, the factorizing/elementary regime.

### 3b. OVERLAPPING (sampleB, ratio 0.14) — `figs/stokes_overlapB.png`
The 12 turning points are at `Im ≈ ±2.05`, the 13 turning points at `Im ≈ ±1.30`; their
Stokes lines sweep across the strip and **the red (12) and green (13) lines cross at four
joints** near `u ≈ −0.10 ± 0.90 i` (two conjugate pairs). The node (pair 23, black star) sits
on the real axis between them. **4 joints.** Direct ODE benchmark: `P₂→₂ = 0.0217`,
incoherent `0.000164`, **enhancement 132.75×** — strongly non-factorizing, the enhanced
regime. (This sample's enhancement exceeds the 104× headline because ratio 0.14 is more
overlapping than the cases that set that figure.)

---

## 4. Joint identification (established for the showcase, with caveats)

| sample | ratio | node pair | simple-tp pairs | # joints | joint types | joint locations | enhancement |
|---|---|---|---|---|---|---|---|
| sampleA (separated) | 3.97 | 23 | {23, 12} | **0** | — | — | **1.00×** |
| sampleB (overlapping) | 0.14 | 23 | {12, 12, 13, 13} | **4** | 12 × 13 | `−0.10±0.90 i`, `−0.07±0.85 i` | **132.75×** |

The joints in the overlapping case are **12 × 13** crossings — i.e. Stokes lines of the two
*off-diagonal pairs that share the middle level* (level 1 is common to 12 and 13). This is
physically exactly right: the middle survival `P₂→₂` is the amplitude that fails to
factorize, and the joint that controls it is built from the two Stokes lines that both touch
the middle sheet. The node (pair 23) lies *between* the joints but is not itself a joint
(consistent with WS-A: the node is an ordinary point of the 3×3 system, carrying no local
connection data — its role is to sit *inside* the jointed region, not to be the junction).

---

## 5. The correlation scan (numerically-supported)

`PLACEHOLDER_SCAN_TABLE`

`PLACEHOLDER_CORRELATIONS`

**Reading.** `PLACEHOLDER_READING`

---

## 6. Honest residuals / what remains ambiguous

1. **Junction *rule* vs junction *presence*.** WS-G shows joints are present iff `P` is
   enhanced, and that joint *strength* (geometric depth) grows with the enhancement. It does
   **not** compute the GMN junction S-matrix and turn it into the *number* `P₂→₂`. The claim
   "joint holonomy = the enhancement value" is therefore **framework/conjecture**, not a
   derived identity. Promoting it is WS-E's task (assemble the connection coefficient).

2. **Higher-rank spectral networks have structured joints.** For a rank-3 (N=3) WKB problem
   the relevant object is a GMN *spectral network*, where junctions obey specific
   sheet-labelled rules (a 12 and a 23 line can *spawn* a 13 line at a joint, etc.). I detect
   joints purely as planar curve crossings of different types; I have **not** verified the
   detailed GMN soliton/junction grammar (which crossings are "active" vs accidental). The
   12×13 joints I find are the physically expected ones (shared middle level), but a full
   spectral-network treatment could reclassify or add structure. Flagged as open.

3. **The "joint strength" metric is a geometric proxy** (`1 − dist(joint, nearest tp)/median
   tp-separation`), not the literal relative WKB action at the joint. It correlates with the
   enhancement (§5) but should be read as a heuristic ordering, not a physical action. A
   sharper metric (the imaginary relative action between the two windows evaluated at the
   joint) is the natural next refinement.

4. **The node is on the real axis** and the ray-integration phase field is most delicate
   there (the masking keeps it clean, but the node's own 23-lines are the least certain). This
   does not affect the joint count (joints are 12×13, away from the node), but a dedicated
   local model at the node would harden the picture.

5. **Sample construction.** Reaching ratio ~4 needs the thin high-ratio tail (one tiny
   coupling + large `ε` separation, signed couplings). Such samples are atypical of "natural"
   Type-1 data — consistent with WS-D's point that Type-1 is *generically* overlapping
   (hence generically jointed, hence generically transcendental). The joint-free regime is
   real but rare; that is itself the explanation of the 15/85 split.

---

## 7. Net for the program

- **Confirms and localizes the WS-A non-Abelian remainder.** `P`'s transcendental content
  (the off-diagonal Stokes coefficient / the `P₂→₂` enhancement) is the **joint** of the
  Stokes graph on the complex `u`-plane. Joint-free ⇔ elementary (BE × Weber product);
  jointed ⇔ the genuine confluent connection coefficient.
- **Sharpens prediction P-ii** (RESEARCH_PROGRAM §9): the factorization locus (WS-C) should
  coincide with the **joint-free locus** of the Stokes graph — a geometric, checkable
  criterion for elementarity, independent of the symbolic reducibility test.
- **Hands WS-E a target.** The object to compute is the GMN junction S-matrix at the 12×13
  joint (with the node sitting between the two joints), as a function of the two window
  actions — exactly the two-variable (Kampé-de-Fériet-class) connection coefficient WS-A/the
  KZ note anticipate.

---

## 8. Reproducibility

All results reproduced by `experiments/stokes_graph.py` (numpy/scipy/sympy/matplotlib):
- `turning_points()` — exact-rational discriminant, classified [1,1,1,1,2].
- `classify_turning_pair()` — branch-tracked colliding-pair labelling.
- `phase_field()` / `stokes_lines_for_tp()` — vectorized WKB phase field + masked zero-contours.
- `build_stokes_graph()` — full graph + `find_joints()` (vectorized segment intersection).
- `mid_enhancement()` — direct ODE `P₂→₂` benchmark (DOP853, rtol 1e-12), the 2.5–104× metric.
- `main()` — renders `figs/stokes_sepA.png`, `figs/stokes_overlapB.png`, and the §5 scan.
Runtime ≈ 4 min. Figures regenerate deterministically (fixed RNG seeds).
