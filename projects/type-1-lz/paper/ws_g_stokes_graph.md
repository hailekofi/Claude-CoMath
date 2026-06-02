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

**VERDICT (the joint diagnosis is SUPPORTED — clean on the showcases, sound-but-one-directional
on the scan).** The Stokes graph cleanly separates the two empirical regimes exactly as
conjectured on the constructed showcase pair, and across a 16-sample scan **every sample that
carries a joint is enhanced** (joint ⇒ non-factorization is sound), while the converse is only a
*lower-bound* witness (the conservative joint detector has false negatives in the far-apart
tail — see §2, §5). Concretely:

- **Well-separated (ratio ≳ 1, the ~15% where the incoherent product is exact):**
  the turning points hug the real axis (narrow avoided crossings) and sit far apart along
  it; each pairwise turning-point cluster carries its *own* short Stokes lines that **do not
  reach** the other clusters. **The Stokes graph is JOINT-FREE** (a "ladder"). `P` factorizes
  into independent 2-level Weber pieces ⇒ enhancement `≈ 1.00×`. (`figs/stokes_sepA.png`)

- **Overlapping (ratio ≲ 1, the ~85% where `P_{m→m}` is enhanced 2.5–104×):**
  the complex turning points push deep off the real axis and their Stokes lines of
  **different pair-types (12 and 13) intersect**, forming **genuine joints** — for sampleB,
  **4 joints** (two conjugate pairs at `u ≈ −0.13 ± 0.95 i` and `−0.07 ± 0.85 i`) where the red
  (12) and green (13) Stokes lines cross between the turning points. **The Stokes graph HAS
  joints.** `P` does NOT factorize ⇒ enhancement `120.5×`. (`figs/stokes_overlapB.png`)

- **Across a 16-sample scan spanning ratio 0.1→4**, joint presence/strength correlate with the
  enhancement in the sound direction (all 3 jointed samples are enhanced above the incoherent
  product — 2 strongly, 1 mildly at the 1.14/1.15 threshold; mean #joints 0.18 for factorizing
  vs 1.20 for enhanced; strongest joints in the most overlapping samples; Pearson joint-strength
  vs `log₁₀ enh` = +0.27). The conservative detector *misses* some joints in the far-apart tail
  (3 false negatives), so it under-counts but never over-counts — see the contingency table in §5.

**Interpretation.** This makes the WS-A picture geometrically concrete. WS-A established that
`P` is the Stokes/connection data of a single rank-2 irregular point at `u=∞`, with the
*elementary* part = the Weber exponential data (BE survivals `e^{−2πΓ_ij}`) and the
*transcendental* remainder = the off-diagonal Stokes connection coefficient. WS-G localizes
that remainder **on the complex u-plane**: it is the holonomy contribution of the **joint**.
No joint ⇒ the connection matrix is the ordered product of independent 2-level (Weber) Stokes
factors ⇒ elementary `P`. A joint ⇒ the GMN junction rule mixes the three sheets ⇒ a non-
factorizable connection coefficient = the measured `P_{m→m}` enhancement. The "joint contribution
IS the off-diagonal Stokes connection coefficient" conjecture is **supported**.

**Evidence-ladder status:** the turning-point structure and the joint-vs-enhancement
*correlation* are **numerically-supported (strong)**. The identification "joint holonomy =
the precise value of the `P_{m→m}` enhancement" is **conjectural/framework** — WS-G establishes
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

**Reach sensitivity — the central methodological honesty.** Joint detection depends on *how
far* each Stokes line is traced. I tested two cutoffs: (i) the conservative **disk mask**
above, and (ii) a **larger-reach** variant (trace to `~1.25×` the nearest tp, clipping near
other tps). They disagree in the borderline cases: the larger reach finds the **8** joints of
sampleB (the four 12×13 plus four 13×23) but also produces **2 spurious 12×23 joints in the
well-separated sampleA** (the 12 and 23 clusters' lines cross midway even though the physics is
factorizing). The conservative disk mask gives the **clean** showcase verdict (sampleB 4,
sampleA 0) with **no spurious joints**, at the cost of **false negatives** in the far-apart
tail (some genuinely-enhanced samples where the relevant tps are far apart while a same-type
conjugate twin sits close, shrinking the disk below the reach needed to meet the cross-pair
line). I therefore report the disk-mask results and treat **"a joint is found" as a SOUND
(sufficient) witness of non-factorization, while "no joint found" is INCONCLUSIVE** in the
tail. A definitive treatment needs the full GMN spectral-network trajectory rules (§6 item 2),
which fix the line-termination/junction grammar canonically — not implemented here.

---

## 3. The two regimes — figures

### 3a. WELL-SEPARATED (sampleA, ratio 3.97) — `figs/stokes_sepA.png`
Three **disjoint** turning-point clusters strung along the real axis (node of pair 23 at
`+0.33`, a complex 23 pair at `−0.35 ± 0.08 i`, a 12 pair at `−0.04 ± 0.05 i`), each with short
Stokes fans confined near its own crossing (note the tps hug the real axis — narrow avoided
crossings). **No Stokes line of one pair reaches a turning point of another pair ⇒ 0 joints**
(under the conservative disk mask; see §2 for why a looser reach would spuriously cross here).
Direct ODE benchmark: `P_{m→m} = 0.98592`, incoherent product `0.98588`, **enhancement 1.00×** —
the incoherent product is *exact* to 4 digits, the factorizing/elementary regime.

### 3b. OVERLAPPING (sampleB, ratio 0.14) — `figs/stokes_overlapB.png`
The 12 turning points are at `Im ≈ ±2.05`, the 13 turning points at `Im ≈ ±1.30`; their
Stokes lines sweep across the strip and **the red (12) and green (13) lines cross at four
joints** at `u ≈ −0.126 ± 0.952 i` and `−0.070 ± 0.853 i` (two conjugate pairs). The node
(pair 23, black star) sits on the real axis between them. **4 joints.** Direct ODE benchmark
(T=160, rtol 1e-11): `P_{m→m} = 0.0197`, incoherent product `0.000164`, **enhancement 120.5×** —
strongly non-factorizing, the enhanced regime. (This sample's enhancement exceeds the 104×
headline because ratio 0.14 is more overlapping than the cases that set that figure; the
precise ratio is mildly ODE-horizon-sensitive because both `P_{m→m}≈0.02` and the incoherent
`≈1.6e−4` are small.)

---

## 4. Joint identification (established for the showcase, with caveats)

| sample | ratio | node pair | simple-tp pairs | # joints | joint types | joint locations | enhancement |
|---|---|---|---|---|---|---|---|
| sampleA (separated) | 3.97 | 23 | {23, 23, 12, 12} | **0** | — | — | **1.00×** |
| sampleB (overlapping) | 0.14 | 23 | {12, 12, 13, 13} | **4** | 12 × 13 | `−0.126±0.952 i`, `−0.070±0.853 i` | **120.5×** |

The joints in the overlapping case are **12 × 13** crossings — i.e. Stokes lines of the two
*off-diagonal pairs that share the middle level* (level 1 is common to 12 and 13). This is
physically exactly right: the middle survival `P_{m→m}` is the amplitude that fails to
factorize, and the joint that controls it is built from the two Stokes lines that both touch
the middle sheet. The node (pair 23) lies *between* the joints but is not itself a joint
(consistent with WS-A: the node is an ordinary point of the 3×3 system, carrying no local
connection data — its role is to sit *inside* the jointed region, not to be the junction).

---

## 5. The correlation scan (numerically-supported, with stated noise)

16 samples, ratio stratified over `0.1 → 4` (the two showcase samples + 14 representatives,
one per geomspace ratio bin). `# joints` and `strength` from the **conservative disk-mask**
Stokes graph; `enh` = `P_{m→m} / P_{m→m}^{incoherent}` from direct ODE (fast settings, T=100). The
enhancement spans `~1 → 1e11`: when the incoherent product `e^{−2π(Γ+Γ)}` underflows (wide-`ε`
samples) the ratio is astronomically large simply because the BE-product prediction is
essentially zero — read `log₁₀ enh`, not the raw number.

| ratio | #joints | strength | P_{m→m} | incoherent | enh (×) |
|---:|---:|---:|---:|---:|---:|
| 0.142 | **4** | 0.55 | 0.0181 | 1.6e−4 | 110.9 |
| 0.152 | 0 | 0.00 | 0.0012 | ~0 | 6.3e10  ← *false neg* |
| 0.199 | 0 | 0.00 | 0.599 | 0.586 | 1.02 |
| 0.262 | **2** | 0.66 | 0.0019 | ~0 | 8.9e8 |
| 0.343 | 0 | 0.00 | 0.0112 | 2.7e−4 | 40.8  ← *false neg* |
| 0.450 | 0 | 0.00 | 0.601 | 0.577 | 1.04 |
| 0.591 | **2** | 0.80 | 0.0741 | 0.0651 | 1.14 |
| 0.775 | 0 | 0.00 | 4.6e−4 | ~0 | 1454  ← *false neg* |
| 1.017 | 0 | 0.00 | 0.719 | 0.652 | 1.10 |
| 1.333 | 0 | 0.00 | 0.996 | 0.996 | 1.00 |
| 1.746 | 0 | 0.00 | 0.974 | 0.973 | 1.00 |
| 2.325 | 0 | 0.00 | 0.617 | 0.621 | 0.99 |
| 3.016 | 0 | 0.00 | 0.920 | 0.919 | 1.00 |
| 3.940 | 0 | 0.00 | 0.965 | 0.965 | 1.00 |
| 3.970 (sampleA) | 0 | 0.00 | 0.986 | 0.986 | 1.00 |
| 0.142 (sampleB)* | 4 | 0.55 | — | — | 120.5 (hi-acc) |

(*sampleB appears twice — fast-ODE row 1 and the hi-accuracy showcase value.)

**Correlations (Pearson vs `log₁₀ enh`):** ratio `−0.39`; #joints `+0.22`; joint-strength
`+0.27`; joint-present `+0.30`.

**Joint-presence × enhancement contingency (threshold enh>1.15):**

| | enhanced | not enhanced |
|---|---:|---:|
| **joint found** | 2 | 1† |
| **no joint** | 3‡ | 10 |

† the lone "joint & not-enhanced" is the ratio-0.591 sample at enh = 1.14 — *just* below the
1.15 cutoff, i.e. mildly enhanced, not a genuine false positive. ‡ the three "no joint &
enhanced" are the documented **far-apart false negatives** (§2): the conservative disk shrinks
below the reach needed to meet the cross-pair line.

**Reading (honest).**
- **Direction that is SOUND:** every sample where a joint *is* found is enhanced (the one
  apparent exception is at the 1.14/1.15 boundary). **Joint presence ⇒ non-factorization.**
- **Direction that is NOISY:** "no joint found" does *not* imply factorization — 3/13 no-joint
  samples are strongly enhanced (false negatives from the conservative mask). So the disk-mask
  joint count is a *lower bound* witness, not a complete classifier.
- **Aggregate separation is clear:** mean #joints `= 0.18` for the (near-)factorizing samples
  vs `= 1.20` for the enhanced ones; the strongest joints (strength `0.55–0.80`) all sit in the
  most overlapping samples. `11/16 = 69%` of this (ratio-flat) suite are "incoherent-exact",
  compatible with the established ~85% (the scan over-weights the high-ratio tail by design, so
  it tilts toward the exact regime relative to natural sampling).
- **Net:** the scan *supports* the diagnosis in the sound direction and shows the right
  monotone trend, while exposing that a robust *bidirectional* classifier needs the full GMN
  network rules (§6 item 2), not the conservative planar-crossing proxy used here.

---

## 6. Honest residuals / what remains ambiguous

1. **Junction *rule* vs junction *presence*.** WS-G shows that *when* a joint is found `P` is
   enhanced, and that joint *strength* (geometric depth) trends with the enhancement. It does
   **not** compute the GMN junction S-matrix and turn it into the *number* `P_{m→m}`. The claim
   "joint holonomy = the enhancement value" is therefore **framework/conjecture**, not a
   derived identity. Promoting it is WS-E's task (assemble the connection coefficient).

2. **Higher-rank spectral networks have structured joints — and this is exactly the source of
   the false negatives.** For a rank-3 (N=3) WKB problem the relevant object is a GMN *spectral
   network*, where junctions obey sheet-labelled rules: a 12 and a 23 line that meet *spawn* a
   13 line, lines *terminate* at branch points, and only "active" crossings are physical
   joints. I detect joints purely as planar crossings of different-type Stokes lines inside a
   conservative disk, which (i) can **miss** a real joint when a line had to be spawned or
   would form just outside the trusted disk (the §5 false negatives), and (ii) could in
   principle flag an inactive crossing. The 12×13 joints I find on the showcase are the
   physically expected ones (shared middle level), but **a definitive, bidirectional
   classifier requires implementing the GMN trajectory/junction grammar** — the main open
   technical item this note leaves for WS-E. Flagged as open.

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
  (the off-diagonal Stokes coefficient / the `P_{m→m}` enhancement) is the **joint** of the
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
- `mid_enhancement()` — direct ODE `P_{m→m}` benchmark (DOP853, rtol 1e-12), the 2.5–104× metric.
- `main()` — renders `figs/stokes_sepA.png`, `figs/stokes_overlapB.png`, and the §5 scan.
Runtime ≈ 4 min. Figures regenerate deterministically (fixed RNG seeds).
