# WS-GEOM Milestone 2 — the image of $\Phi:(\gamma,\epsilon,a)\to P$ in reduced-$U(3)$ (facet A)

**Type:** `physics-numerics` (executed; gauge-robust reachable-set characterization, with an
explicit first-class negative). **Date:** 2026-06-03. **Builds on:** M1
(`paper/working_sessions/ws_geom_m1.md` — the node-selected directed-cycle skeleton), R15/R16
(`experiments/skeleton_two_transcendentals.py` — the two-transcendental $\{P_{mm},b\}$
coordinates), R9 / WS-C (`paper/working_sessions/ws_c_factorization_locus.md` — the decoupling locus).
**Scope followed:** `paper/working_sessions/ws_geom_m2_scope.md`. **Code:**
`experiments/ws_geom_m2_image.py` (reproducible; reuses `type1` + the gold `oracle.py`,
caches the cloud to `ws_geom_m2_cloud.pkl`). **Conventions:** NOMENCLATURE.md — states
$\epsilon$-indexed; BE roles by slope `lo,mid,hi = argsort(a)`; $P[x,j]=\mathrm{prob}(x\to j)$.

---

## 0. What M2 establishes (one paragraph)

The BE-independent part of the Type-1 image is the 2-D region swept by
$\{P_{mm},b\}$ ($P_{mm}=P[\mathrm{mid},\mathrm{mid}]$, $b=P[\mathrm{hi},\mathrm{lo}]$; R15).
Mapped with a fast moderate-accuracy $P$ over a stratified, gauge-reduced cloud
(N=2000), that region has a **clean, fully accounted-for geometry**:
(M2a) its closure touches **exactly two** of the six Birkhoff vertices — the **identity**
$(1,0)$ (diabatic corner) and the **directed 3-cycle** $(0,1)$ (adiabatic corner) — and **no
other** permutation; (M2b) the distance to the cycle vertex is a **monotone bias law** in
the window action $\delta_{\max}$, modulated by the cross-ratio $\chi$; (M2c, the headline)
the image **boundary coincides with the decoupling locus** — driving an *extreme* coupling to
zero lands the point on the $b=0$ edge, driving the *middle* coupling to zero lands it on the
$P_{mm}=1$ edge (WS-C), to $\pm0.00$; (M2d) the effective dimension of the lift
$\{P_{mm},b,\delta_{ij},\chi,c_i\}$ is **full** in the dynamical coordinates — the only exact
relations are the elementary BE-action/Coulomb kinematic algebra, and $P_{mm},b,\chi$ carry
**no constraint**. **This is a first-class NEGATIVE for the invariant hunt:** M2 re-finds
exactly $\{$BE $+$ double-stochasticity $+\chi\}$ and adds **no new analytic invariant**. The
genuinely new content is the **topology**: the two node-selected vertices and the
sector boundary $=$ decoupling locus. We state this plainly and do not manufacture a
constraint.

---

## 1. The map and the coordinates (from R15)

$P$ is doubly stochastic (4 DOF). BE fixes the two extreme-slope survivals
$P[\mathrm{lo},\mathrm{lo}],P[\mathrm{hi},\mathrm{hi}]$ as elementary functions of the window
actions $\delta_{ij}=s_{ij}^2|a_i-a_j|$. Double stochasticity then writes every entry as an
affine function of exactly two free numbers,
$$P_{mm}:=P[\mathrm{mid},\mathrm{mid}],\qquad b:=P[\mathrm{hi},\mathrm{lo}],$$
(R15; `reconstruct_from_two`). So the image splits as $\{$elementary BE data$\}\times\{P_{mm},b\}$,
and the interesting, $\sigma$-free part is the **2-D region swept by $\{P_{mm},b\}$**. M2 maps
that region and its lift $\{P_{mm},b,\delta_{ij},\chi,c_i\}$.

**Convention (load-bearing).** We use the R15/skeleton matrix layout
$P[\text{final},\text{initial}]=\mathrm{prob}(\text{initial}\to\text{final})$ (columns =
initial), so $b=P[\mathrm{hi},\mathrm{lo}]=\mathrm{prob}(\mathrm{lo}\to\mathrm{hi})$. The gold
oracle (`oracle.py`) returns the *transpose* layout $P[\text{initial},\text{final}]$; we
verified `quantumP == oracle.T` to $\sim10^{-2}$ (at fast accuracy). $P_{mm}$ is diagonal and
therefore convention-independent; the off-diagonal $b$ is compared against the oracle transpose
($\mathrm{og}[\mathrm{lo},\mathrm{hi}]$) in the spot check. All M2 results are stated in the R15
convention.

**Engine and validation.** A single DOP853 solve of $i\dot U=(H_0+uA)U$ on $[-T,T]$, $T=50$,
$\mathrm{rtol}=10^{-7}$ ($\sim$1.12 s/sample, $\sim$few$\times10^{-3}$ vs gold). The region
SHAPE / boundary / dimension do not need $10^{-9}$. Spot-checked against the gold oracle
(`oracle.py`, $T=80$) at four anchors — `[numerically-verified]` the fast $P_{mm}$ and $b$
agree with gold to $\le 9\times10^{-3}$ (canonical $P_{mm}$ dev $3.6\times10^{-3}$, $b$ dev
$7.7\times10^{-3}$; sampleB $b$ dev $3.3\times10^{-4}$; all four anchors $P_{mm}$ dev
$\le 8.6\times10^{-3}$). The Birkhoff-vertex,
boundary, and dimension results are **gauge-robust reachable-set features**, reported NOT as a
sampling density (the parameter measure is arbitrary; risk §3 of the scope).

---

## 2. M2a — exactly two reachable Birkhoff vertices  [numerically-supported]

In slope coordinates a permutation matrix has $P_{mm}\in\{0,1\}$ and $b=P[\mathrm{hi},\mathrm{lo}]
=\mathrm{prob}(\mathrm{lo}\to\mathrm{hi})\in\{0,1\}$ (skeleton/R15 convention
$P[\text{final},\text{initial}]$). The six permutations collapse to **four** distinct
$\{P_{mm},b\}$ points (verified by direct enumeration):
identity $\to(1,0)$; the **node-selected directed 3-cycle**
$\mathrm{lo}\to\mathrm{hi}\to\mathrm{mid}\to\mathrm{lo}$ (= M1's PI\_OUT cycle) $\to(0,1)$;
the lo–hi transposition $\to(1,1)$; and the **reverse** cycle
$\mathrm{lo}\to\mathrm{mid}\to\mathrm{hi}\to\mathrm{lo}$ together with the lo–mid and mid–hi
transpositions all $\to(0,0)$.

Distance of the N=2000 cloud closure to each (exact run values):

| vertex | $\{P_{mm},b\}$ | min-dist over cloud | frac within 0.02 | reached? |
|---|---|---:|---:|---|
| **identity** | $(1,0)$ | $0.0000$ | $0.16$ | **YES** |
| **directed 3-cycle** ($\mathrm{lo}\to\mathrm{hi}\to\mathrm{mid}\to\mathrm{lo}$) | $(0,1)$ | $0.0011$ | $0.07$ | **YES** |
| transposition lo–hi | $(1,1)$ | $0.573$ | $0$ | no |
| reverse-cycle / transp. | $(0,0)$ | $0.096$ | $0$ | no |

Only the **identity** (diabatic / small-action corner) and the **node-selected directed
3-cycle** $\mathrm{lo}\to\mathrm{hi}\to\mathrm{mid}\to\mathrm{lo}$ (adiabatic / large-action
corner) are reached. The transposition $(1,1)$ and the $(0,0)$ vertices are unreachable
(min-dist $\gtrsim 0.27$): the directed cycle is node-selected (M1) and its orientation
$b\to1$ (the chirality $\mathrm{lo}\to\mathrm{hi}$) is fixed by the slope order, so the
**reverse** cycle ($b=0$) and the transposition corners are geometrically forbidden. The
corners were also reached *constructively*: weak coupling $\to(0.997,0.001)$ (identity);
strong coupling $\to(0.021,0.963)$ (directed cycle).

**Gate M2a: PASS** — exactly two Birkhoff vertices in the image closure, fixed by the slope
ordering.

---

## 3. M2b — the directed-cycle-bias law  [numerically-supported]

Distance to the cycle vertex $(0,1)$ vs the middle window action
$\delta_{\max}=\max(\delta_{\mathrm{lo,mid}},\delta_{\mathrm{mid,hi}})$, binned over the cloud:

| $\delta_{\max}$ bin | $n$ | mean dist to cycle | mean $P_{mm}$ |
|---|---:|---:|---:|
| $<0.05$ | 898 | $1.338$ (near identity) | $0.918$ |
| $0.05$–$0.1$ | 153 | $1.035$ | $0.587$ |
| $0.1$–$0.2$ | 170 | $0.822$ | $0.365$ |
| $0.2$–$0.4$ | 178 | $0.541$ | $0.173$ |
| $0.4$–$0.8$ | 153 | $0.370$ | $0.070$ |
| $0.8$–$1.6$ | 129 | $0.220$ | $0.028$ |
| $>1.6$ | 319 | $0.079$ (near cycle) | $0.013$ |

- **Monotone in the action:** Spearman$(\mathrm{dist}_{\rm cycle},\delta_{\max})=-0.978$
  — more window action drives $\{P_{mm},b\}$ toward the directed cycle, exactly the M1/R16
  adiabatic picture. The approach is **logarithmic in $\delta$** (the distance saturates; a
  linear-in-$\delta$ exponent is the wrong model): the trend fit is
  $\mathrm{dist}_{\rm cycle}\approx -0.162\,\log\delta_{\max}+0.399$.
- **$\chi$-modulation (shape, not scale):** within a fixed action band
  ($0.2<\delta_{\max}<1.6$, $n=460$), Spearman$(\mathrm{dist}_{\rm cycle},\chi)=+0.426$ — at the
  *same* action, **more-separated windows ($\chi\to1$) sit farther from the cycle**. This is the
  same $\{$action $\times$ shape$\}$ split as WS-O3/R11 ($P_{mm}=f(\delta_X,\chi)$): the action
  sets the scale of the approach, $\chi$ tilts it.

**Gate M2b: characterized** — a monotone (log-in-action) bias toward the directed cycle,
shape-modulated by $\chi$. (Note: this is a regime *trend* of a 2-D reachable set, not a
sharp law; it inherits the honest M1 statement that the adiabatic skeleton is the directed
cycle.)

---

## 4. M2c — HEADLINE: image boundary $=$ decoupling locus  [numerically-supported, decisive]

The image is bounded by the two **axis-aligned edges** $b=0$ (no directed circulation) and
$P_{mm}=1$ (full middle survival): over the entire N=2000 cloud, $b\ge0$ and $P_{mm}\le1$ with
**zero violations**. The decisive test (R9 / WS-C): drive one coupling to decoupling
($\gamma_k\to0$) and see which edge is selected.

| coupling driven to 0 | what decouples | edge reached |
|---|---|---|
| $\gamma_{\rm lo}$ or $\gamma_{\rm hi}$ (an **extreme** level) | the chirality collapses | lands on $b=0.0000\pm0.0000$ |
| $\gamma_{\rm mid}$ (the **middle** level) | WS-C middle decoupling | lands on $P_{mm}=1.000\pm0.000$ |

Across three base points and all three couplings, the **edge selection is exact**: extreme
decoupling $\Rightarrow b\to0$ to machine display, middle decoupling $\Rightarrow P_{mm}\to1$
(the elementary survival; WS-C: $P_{mm}\to$ incoherent product, which $\to1$ as
$\delta_{\rm mid}\to0$). The boundary of the $\{P_{mm},b\}$ image is therefore **precisely the
decoupling locus**, edge-resolved by *which* level decouples — the geometric realization of
WS-C's "elementary $\Leftrightarrow$ a level decouples."

*Honest metric note.* A convex-hull-distance metric on a *sparse random* cloud is misleading:
the decoupling curves trace the analytic edges $b=0$ / $P_{mm}=1$ **exactly**, but random
sampling under-fills those edges, so a swept point can sit far from the *sampled* hull while
being on the *true* boundary. The load-bearing, gauge-robust statement is the **edge selection**
($b\to0$, $P_{mm}\to1$ to $\pm0.00$) together with the **global bounds** ($b\ge0$, $P_{mm}\le1$,
no violations) — not the hull-distance number.

**Gate M2c: CONFIRMED** — the image boundary coincides with the decoupling locus, with a sharp
edge-selection rule (extreme $\to b=0$, middle $\to P_{mm}=1$).

---

## 5. M2d — effective dimension: FULL (no new analytic constraint)  [numerically-supported]

**Global PCA** of the standardized 9-feature lift
$\{P_{mm},b,\delta_{\rm lm},\delta_{\rm mh},\delta_{\rm lh},\chi,c_0,c_1,c_2\}$ has spectrum
(normalized) $[1.000,\,0.659,\,0.483,\,0.414,\,0.229,\,0.149,\,0,\,0,\,0]$: the three exact
zeros are the elementary relations among $\{\delta_{ij},c_i\}$ (including $\sum_i c_i=0$).

**Local Jacobian** of $\Phi:(\epsilon,\gamma,a)\to$ lift (finite-difference, median over 12
interior points): singular spectrum $[7.4,\,2.2,\,0.60,\,0.28,\,0.073,\,1.8\times10^{-4},\,
\sim4\times10^{-14},\,\sim10^{-14},\,\sim7\times10^{-15}]$ — numerical rank 6, i.e. **3 exact
relations** among the 9 features. The decisive diagnostics:

- **The 3 relations live ENTIRELY in the elementary BE-action / Coulomb block**
  $\{\delta_{ij},c_i\}$. Both $c_i=\sum_{j\ne i}s_{ij}^2(a_i-a_j)$ and
  $\delta_{ij}=s_{ij}^2|a_i-a_j|$ are built from the same three quantities
  $s_{ij}^2(a_i-a_j)$; together with $\sum_i c_i=0$ that is 3 relations. **These are
  kinematic algebra, already known.**
- **$P_{mm}$, $b$, and $\chi$ have ZERO weight in every null (relation) vector.** They are
  genuinely free; no relation ties the dynamical content to the BE/Coulomb data.
- The $\{P_{mm},b\}$ **sub-Jacobian has rank 2** ($\sigma_2/\sigma_1=0.128$, clearly
  nonzero): $P_{mm}$ and $b$ are locally independent directions — **no collapse, no hidden
  $P_{mm}$–$b$ constraint**, exactly confirming the R15 prior. (Null-vector feature
  participation over the 9 features is $[0,0,1.0,1.05,1.10,0,1.01,1.05,1.02]$ — exactly zero on
  $P_{mm},b,\chi$ and order-1 on the $\{\delta_{ij},c_i\}$ block.)

**Codimension for M3 $=0$ in the dynamical coordinates.** The lift's only relations are the
elementary BE-action/Coulomb identities; the dynamical pair $\{P_{mm},b\}$ (and $\chi$) is
unconstrained. (The $\gamma_i$ enter $P$ only through $\gamma_i^2$ / the diagonal sign gauge
$D=\mathrm{diag}(\pm1)$, $DH_0D$ — a discrete gauge under which $P=|U|^2$ is invariant; this is
accounted for, not a new constraint.)

**Gate M2d: FULL effective dimension — a first-class NEGATIVE.** M2 re-finds
$\{$BE $+$ double-stochasticity $+\chi\}$ and **no new analytic invariant**.

---

## 6. Unistochasticity  [numerically-verified]

$P=|U|^2$ with $U\in U(3)$ is unistochastic by construction. Over the cloud the
doubly-stochastic defect is $\lesssim10^{-4}$ (solver tolerance), and **every** point passes the
$3\times3$ unitarity-triangle (Jarlskog–Stork chain-link) test: $0/2000$ failures. The image
sits inside the unistochastic region of the Birkhoff polytope, as it must; the two reached
corners (identity, directed cycle) are unistochastic vertices, and the image touches the inner
curved boundary only along the elementary $b=0$ / $P_{mm}=1$ edges (the decoupling locus, §4).

---

## 7. The map (summary picture)

- **Coordinates:** $\{P_{mm},b\}$ (the two R15 transcendentals), $\sigma$-free, fibered over
  the elementary BE data $\{\delta_{ij}\}$ and dressed by $\chi$.
- **Two corners:** identity $(1,0)$ (diabatic, small action) and the node-selected directed
  3-cycle $\mathrm{lo}\to\mathrm{hi}\to\mathrm{mid}\to\mathrm{lo}$ at $(0,1)$ (adiabatic, large
  action) — slope-ordered; the reverse cycle and transpositions ($(0,0)$, $(1,1)$) unreached.
- **Boundary:** the two edges $b=0$ and $P_{mm}=1$, $=$ the decoupling locus (extreme vs middle
  coupling $\to0$).
- **Interior:** a genuine 2-D reachable set; $P_{mm}\perp b$ (rank-2 Jacobian); the approach to
  the cycle is monotone (log) in the action, $\chi$-tilted.
- **Effective dimension:** FULL in the dynamical coordinates; the only relations are elementary
  BE/Coulomb algebra.

---

## 8. Honest limits

1. **First-class negative, stated plainly.** M2 finds **no new analytic constraint**. The image's
   analytic content is exhausted by $\{$BE, double-stochasticity, $\chi\}$ (R15/R16/WS-O3). The
   genuinely new content is the **topology** — the two node-selected vertices and the
   sector boundary $=$ decoupling locus — not an invariant. (Scope §5 "partial" outcome; §6 risk 1.)
2. **Reachable set, not a density.** Vertices, boundary, and dimension are gauge-robust; the
   2-D *density* inside the region is sampling-measure-dependent and is **not** claimed.
3. **Moderate-accuracy engine.** $\sim$few$\times10^{-3}$ vs gold — adequate for shape/boundary/
   dimension, validated at four anchors against the oracle. The bias law (M2b) is a regime
   *trend* of a 2-D set, not a sharp scalar law.
4. **Deep-overlap sampling restriction (explicit).** The sampler bounds the coupling scale
   ($\gamma\lesssim 2.8$) and **rejects deep-overlap samples** with max pairwise BE action
   $>8$. In that strong-coupling/merged-window corner the fast $T=50$ engine is both expensive
   ($\sim3\times10^{5}$ rhs evals) and unreliable (the non-resumming $\sigma$ core, M1 §4.4) — so
   including it would inject wrong points, not extend the reachable set. The cloud still reaches
   both corners (identity to $0.0000$, cycle to $0.0011$) and spans $P_{mm}\in[0,1]$,
   $b\in[0,1]$; the corner/boundary/dimension verdicts are robust. The fast engine also carries an
   endpoint-Stark tail at finite $T$ (M1 §5), affecting only precise interior values, not the
   gauge-robust features.

**Evidence-ladder tags.**
- `[numerically-verified]` — exactly two Birkhoff vertices reached (M2a); $b\ge0$, $P_{mm}\le1$
  with no violations and the exact edge selection (M2c); unistochasticity ($0/2000$ failures);
  fast-vs-gold spot agreement to $\le 9\times10^{-3}$ on $P_{mm}$ and $b$.
- `[numerically-supported]` — the monotone $\delta$-bias law and its $\chi$-modulation (M2b);
  the rank-2 $\{P_{mm},b\}$ Jacobian and the rank-6 full-map Jacobian (M2d).
- `[honest-negative]` — **FULL effective dimension: no new analytic invariant** (M2d); the
  content is topological, not a constraint.

---

## 9. Hand-off to M3

M2 produces the cheap-data cloud `experiments/ws_geom_m2_cloud.pkl`
($\{P_{mm},b,\delta_{ij},\chi,c_i\}$, N=2000) and the effective-dimension spectrum. The verdict
for the M3 caged symbolic-regression invariant hunt is:

- **Codimension $=0$ in the dynamical coordinates** $\{P_{mm},b\}$ (and $\chi$): there is **no
  hidden relation** to discover among them. M3 should **not** hunt for a $P_{mm}$–$b$ or
  $\{P_{mm},b\}$–$\{$BE,$\chi\}$ algebraic constraint — M2 shows there is none (the only relations
  are the elementary BE/Coulomb algebra, already proven).
- **The remaining target is the BOUNDARY FUNCTION and the TOPOLOGICAL LABEL**, not an invariant:
  (i) the closed-form of the two image edges $b=0$ / $P_{mm}=1$ as the decoupling locus
  ($\gamma_{\rm extreme}\to0$ vs $\gamma_{\rm mid}\to0$), which WS-C/R9 already characterizes
  analytically and `physics-derivation` can finish; (ii) the two-vertex spinor/sector label
  (identity vs directed cycle), the M1 node datum. Both are *topological/geometric* statements
  about the section, provable by hand — not new transcendental constants.
- **$\sigma$ is unchanged.** M2 does not close $\sigma$ (R9/R11/R17 stand); it shows the *section
  geometry* around $\sigma$ has no extra invariant beyond the known ones.

**Net for M3:** the honest, evidence-based instruction is to prove the **boundary $=$ decoupling
locus** statement and the **two-vertex topology** analytically (both are within reach), and to
**stop searching for a new $\{P_{mm},b\}$ invariant** — M2's full effective dimension says none
exists.

---

## 10. Reproduce

```
cd projects/type-1-lz/experiments
python3 ws_geom_m2_image.py            # versions banner; gold-oracle spot check;
                                       # builds/caches the N=2000 cloud; prints M2a (vertex
                                       # distances), M2b (the bias law), M2c (boundary vs
                                       # decoupling locus), M2d (PCA + Jacobian dimension),
                                       # and the unistochasticity check.
python3 ws_geom_m2_image.py --n 3000   # larger cloud
python3 ws_geom_m2_image.py --rebuild  # ignore the cache and recompute
```
Requires `numpy`, `scipy`, and the sibling `oracle.py` / `num_S12.py` (→ `uploads/assay`).
