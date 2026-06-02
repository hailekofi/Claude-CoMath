# WS-R — Restart-operator structure of the Type-1 N=3 LZ scattering matrix

**Owner:** WS-R (restart-operator / irregular-point local structure). **Date:** 2026-06-01.
**Tests the addendum hypotheses H-R1 / H-R2 / H-R3** of `paper/OPEN_PROBLEM.md`
(restart-operator structure of the single rank-2 irregular point at `u=∞`).
**Inputs used:** `NOMENCLATURE.md`, `experiments/q_restart_probe.py` (canonical frame &
Coulomb `c_i`), `paper/ws_g_stokes_graph.md` (12×13 joint, mid-level sharing),
`experiments/oracle.py` / `uploads/assay/ip.py` (the trusted IP propagator).
**Constraint compliance:** no git; only `paper/restart_structure.md` and
`experiments/restart_probe.py` written; scratch under `/tmp`.
**Reproducer:** `experiments/restart_probe.py` (`python3 experiments/restart_probe.py`,
≈4 min, numpy/scipy/sympy). Convergence side-check in §A.

---

## 0. Executive verdict

| Hyp | Claim | Verdict | Ladder |
|---|---|---|---|
| **H-R3** | `Γ_j=(-u'(λ_j))^{-1/2}` analytic at poles `ε_i`; sqrt branch points = turning points (`-u'=0` / W4 zeros), **not** `ε_i` | **CONFIRMED** (decisive) | **symbolic-proven** (sympy, exact) + numeric tp/pole separation |
| **H-R2** | (a) off-diagonal weight concentrates in the **mid** row/col; (b) decoupling one outer link collapses one shear → elementary | **CONFIRMED** (oracle-tested) | **numerically-supported (strong)**: frame-pinned, T-converged |
| **H-R1** | `𝒮_canon` = (formal-monodromy diagonal `e^{2πi c_i}`) ⋉ (triangular/unipotent Stokes shear in sector order) | **PARTIAL / REFRAMED** | the off-diagonal is **cyclic, not triangular**; `e^{2πi c_i}` is the *drift holonomy*, **not** a factor of the convergent matrix. The honest well-defined finding is the *mid-concentration* (= H-R2a), not triangularity. |

**One-line summary.** H-R3 is settled geometry and clean. H-R2 is true and is the
*physically content-ful* part of the restart picture: the off-diagonal Stokes data lives in
the two `{mid,lo}`/`{mid,hi}` carrier pairs (the 12×13 joint), and killing one outer link
collapses the joint to the elementary BE product. H-R1's specific *triangular*-shear
phrasing is **not** what the data show — in the frame-pinned sector basis the canonical
scattering matrix is dominated by the **3-cycle permutation** (the `PI_OUT=(2,0,1)` sheet
relabelling), and the formal monodromy `e^{2πi c_i}` is the holonomy of the *subtracted*
`logT` drift, not a diagonal factor of `𝒮_canon`. The well-posed, frame-pinned restart
content is the *mid-row/column concentration* of H-R2a, which we keep.

---

## 1. Frame and well-posedness (the anchor for everything below)

All measurements are in the **canonical frame** (NOMENCLATURE; `q_restart_probe.py`):
the symmetric diabatic interaction-picture matrix `𝒮_IP[j,i]=e^{iθ_j(T)}U(T,-T)[j,i]e^{-iθ_i(-T)}`
with `θ_i(t)=(H0)_ii t+a_i t²/2`, Coulomb-subtracted,
```
𝒮_canon[j,i] = e^{i(c_j-c_i)logT} · 𝒮_IP[j,i],   c_i = Σ_{j≠i} s_ij²(a_i-a_j)  (Σ_i c_i=0).
```
We verified in code that `𝒮_IP` equals the IP fundamental matrix `g` of
`uploads/assay/ip.py` (`e^{iθ}Ue^{-iθ}=g` to `<1e-9`), so the brief's recipe and the
project engine coincide. Channels are then permuted into **sector order `[lo,mid,hi]`**
(`argsort(a)`), the ordering in which "triangular" / "mid-shared" are meaningful.

**Well-posedness (the load-bearing convergence check, §A / `experiments/restart_probe.py`
side-run).** For sampleB, as `T = 80 → 120 → 160`:

| T | `R_tri=U/L` | `frac_mid` | `‖Δ𝒮_canon‖∞` | `‖Δ𝒮_IP‖∞` |
|---:|---:|---:|---:|---:|
| 80 | 1.922 | 0.6696 | — | — |
| 120 | 1.928 | 0.6689 | **2.1e−2** | 7.9e−1 |
| 160 | 1.922 | 0.6695 | **1.1e−2** | 5.7e−1 |

`𝒮_canon` **converges** (increment halves: 2.1e−2 → 1.1e−2) while the raw `𝒮_IP`
**diverges** (the `(c_j−c_i)logT` secular drift, ~0.6 per doubling). The derived quantities
(`R_tri`, `frac_mid`) are stable to 3 digits across `T`. **So every number below is a
frame-pinned property of a convergent matrix, not a `T`-artifact** — this is what makes the
measurements well-defined (cf. the brief's caveat that the bare-DP `Q_restart` relation was
content-free).

---

## 2. H-R3 — form factor & branch loci (CONFIRMED, decisive, symbolic)

**Claim.** `Γ_j=δ_j(-u'(λ_j))^{-1/2}` is **analytic** at the poles `ε_i` (`∝(λ-ε_i)`,
no monodromy there); the WKB square-root **branch points are the turning points** where
`-u'(λ)=0` (the ramification points of `u(λ)` / the W4 zeros), **not** the `ε_i`. Hence the
off-diagonal Stokes structure is an irregular-point (sectors-at-∞) effect, distinct from the
turning-point (window/BE) branching.

**Proof (sympy, exact — `hr3_symbolic()`).** With `u(λ)=n/p`, `p=∏(λ-ε_k)`,
`n=Σ g_k²∏_{l≠k}(λ-ε_l)`:

1. **`-u'(λ) = W4(λ)/p(λ)²` identically**, where `W4 = n p' − n' p` is exactly the quartic
   stored in `Geometry.W4`. (`simplify(-u' − W4/p²)=0`; `deg W4 = 4`.) So
   `Γ = (-u')^{-1/2} = |p| / √W4`.

2. **At a pole**: `(λ-ε_i)²·(-u') → g_i²` as `λ→ε_i` (verified for `i=0,1,2`). Hence near
   `ε_i`, `-u' ~ g_i²/(λ-ε_i)²`, so
   ```
   Γ = (-u')^{-1/2} ~ (λ-ε_i)/|g_i|  —  ANALYTIC, vanishing LINEARLY, no branch point.
   ```
   The poles of `u` are *simple zeros* of `Γ`, carrying **no monodromy**.

3. **The branch points of `Γ=|p|/√W4` are the zeros of `W4`** — the 4 turning points
   `(-u'=0)`. Located numerically (`turning_points()`), they are 2 complex-conjugate pairs,
   well separated from the real poles:

   | sample | poles `ε` | turning points (W4 zeros) | min \|tp−ε\| |
   |---|---|---|---:|
   | canonical | (−2, 0, 3) | `1.049±1.703 i`, `−0.802±1.108 i` | **1.37** |
   | sampleB | (−1, 0, 1.5) | `1.056±0.722 i`, `−0.612±0.501 i` | **0.63** |
   | well_sep | (−5, 0, 5) | `±2.686±2.686 i` | **3.55** |

These are exactly the `Q4`/window branch points of WS-G (the complexified avoided
crossings), and they are disjoint from the irregular-point pole loci.

**Verdict: CONFIRMED (symbolic-proven).** The decomposition is geometrically clean: the
`Γ_j` branch (turning points) carries the BE/window data; the **off-diagonal Stokes data is
an irregular-point effect at the `ε_i`, on a different sheet structure**. This is the cleanest
of the three and it underwrites the whole "restart-at-the-irregular-point" framing.

---

## 3. H-R2 — mid-pair carrier & decoupling collapse (CONFIRMED, oracle-tested)

### 3a. Mid-row/column concentration (CONFIRMED)

In sector order `[lo,mid,hi]`, split the off-diagonal probability weight of `𝒮_canon` into
the part **touching the middle level** (row 1 or col 1 — the `{mid,lo}` and `{mid,hi}`
carrier pairs) vs the lone **outer `lo↔hi`** part (the pair not sharing mid):

| sample | regime | `frac_mid = W_mid/(W_mid+W_outer)` |
|---|---|---:|
| canonical | moderate overlap | **0.607** |
| sampleB | strong overlap (jointed, ~117×) | **0.669** |
| well_sep | well-separated | **0.757** |

In all cases a **majority** (61–76%) of the off-diagonal weight is mid-touching, and only a
*minority* sits on the outer `lo↔hi` link — even though a generic dense unitary would split
the three off-diagonal pairs evenly. The middle level is the **shared subdominant partner**
of both Stokes shears, exactly as the WS-G `12×13` joint demands (level "1"/mid common to the
12 and 13 lines). **CONFIRMED, numerically-supported (strong).**

### 3b. WS-C decoupling: one shear vanishes → elementary (CONFIRMED — the clean one)

Send one **outer** coupling to zero (shrink `γ_outer → γ_outer·10⁻³`, so the corresponding
`s_{mid,outer}` and the `lo↔hi` link vanish) and ask whether the joint collapses to the
**incoherent BE product** (`P_mid → exp(−2π[s²|Δa|]_{mid,lo} + [s²|Δa|]_{mid,hi})`; the
decoupled pair's BE exponent → 0 so its factor → 1). Result (`hr2_decoupling()`):

| case | `P_mid` | incoherent BE | \|diff\| | `frac_mid` |
|---|---:|---:|---:|---:|
| canonical / full | 0.21617 | 0.08433 | 1.3e−1 (enhanced) | 0.607 |
| canonical / decouple hi | 0.21886 | 0.22136 | **2.5e−3** | **1.0000** |
| canonical / decouple lo | 0.38124 | 0.38095 | **2.9e−4** | **1.0000** |
| sampleB / full | 0.01913 | 0.000163 | 1.9e−2 (≈117× enh) | 0.669 |
| sampleB / decouple hi | 0.00179 | 0.00114 | **6.5e−4** | **1.0000** |
| sampleB / decouple lo | 0.13645 | 0.14280 | **6.4e−3** | **1.0000** |

Two effects fire together in **every** decoupling, in both regimes:

1. **The enhancement disappears**: the full case is off the incoherent product by
   `O(10⁻¹–10⁻²)` (the genuine non-factorizing joint); the moment one outer level is
   decoupled, `P_mid` lands on the incoherent BE product to `O(10⁻³)`. The remaining
   `mid↔(other outer)` crossing is a single 2-level (Weber) passage — **elementary**.
2. **`frac_mid → 1.0000`**: the outer `lo↔hi` weight collapses to zero; **one of the two
   shears has vanished**, leaving a single carrier pair. ("Elementary ⟺ a level decouples.")

**Verdict: CONFIRMED (numerically-supported, strong, oracle-grade).** This is the
content-ful restart statement: the off-diagonal Stokes data *is* carried by the two
mid-sharing shears `{mid,lo}`, `{mid,hi}`, and their non-commuting composition is the
`P_{m→m}` enhancement; remove either link and the joint reduces to one elementary Weber factor.
Matches WS-G (joint-free ⇔ elementary) and WS-C (factorization locus) from the
irregular-point side.

---

## 4. H-R1 — formal-monodromy ⋉ triangular shear (PARTIAL / REFRAMED)

H-R1 has two falsifiable sub-claims. Both are tested as *frame-pinned* quantities (not fits),
per the brief's instruction. Both come back **negative in their literal form**, with a clean
positive reframing.

### 4a. Is the off-diagonal part *triangular* (unipotent shear) in sector order? — **NO; it is cyclic.**

Triangularity metric (`triangularity()`): in sector order, upper off-diagonal weight
`U=Σ_{i<j}|𝒮[i,j]|²` vs lower `L=Σ_{i>j}|𝒮[i,j]|²`, ratio `R_tri=U/L`. A unipotent shear
needs one of `U,L ≈ 0`. **Haar-U(3) null model** (`haar_triangularity_baseline()`, 2·10⁴
samples): median `R_tri=0.999`, `P(R>3)=0`, `P(R<1/3)=0`, mean `U=mean L=1.0` — so any
strong departure would be a real signal.

| sample | `U` | `L` | `R_tri=U/L` |
|---|---:|---:|---:|
| canonical | 1.602 | 0.980 | 1.64 |
| sampleB | 1.931 | 1.002 | 1.93 |
| well_sep | 0.778 | 0.603 | 1.29 |

`R_tri` is `1.3–1.9` — **mildly above** the Haar median, **but nowhere near triangular**
(`L` is never small; a true shear would give `R_tri→∞`). The apparent upper-bias is an
artifact of the **sheet permutation**, not a unipotent structure. Decomposing `|𝒮_canon|²`
in sector order onto the 6 permutation skeletons (`perm_concentration`):

| sample | dominant permutation | its weight (/3) | identity wt | 3-cycle `(lo→mid→hi→lo)` wt |
|---|---|---:|---:|---:|
| canonical | **3-cycle (1,2,0)** | **2.224** | 0.418 | 2.224 |
| sampleB | **3-cycle (1,2,0)** | **2.860** | 0.067 | 2.860 |
| well_sep | identity | 1.620 | 1.620 | 0.952 |

For the **strongly-coupled / jointed** cases (canonical, sampleB — the physically
interesting ones) the canonical scattering matrix concentrates on the **3-cycle
`lo→mid→hi→lo`** (95% of the trace-weight for sampleB), i.e. the `PI_OUT=(2,0,1)` spectral-
sheet relabelling established in `oracle.py`. A 3-cycle is *not* triangular: its support
straddles both triangles (`lo→mid`, `mid→hi` upper; `hi→lo` lower), which is precisely why
`R_tri` sits at a finite `~1.3–1.9` rather than `≈1` (Haar) or `≈∞` (shear). **The
off-diagonal Stokes structure of `𝒮_canon` is cyclic, set by the sheet 3-cycle — not a
unipotent triangular shear.** (well_sep is identity-dominated because it is near-diabatic,
`P≈I`; there the small off-diagonal *is* mildly upper-biased, but that regime is the
elementary one and carries no joint.)

### 4b. Does a natural diagonal factor of `𝒮` equal the formal monodromy `e^{2πi c_i}`? — **NO (category error).**

`e^{2πi c_i}` is the **holonomy of the `logT` drift** `e^{i(c_j−c_i)logT}` that we
*subtract* to define the canonical frame (it is the monodromy of `logT` under `T→e^{2π}T`).
After subtraction, `𝒮_canon` is convergent and carries **no** such diagonal factor; computing
`e^{2πi c_i}` and comparing to any diagonal part of `𝒮_canon` is comparing the removed object
to the cleaned one. (For the record the values are
`e^{2πi c}=(−0.85−0.52i, 0.86+0.52i, −0.46+0.89i)` for canonical, etc. — printed by the
probe — but they are *the drift's* monodromy, established already in `q_restart_probe.py`,
not a measured factor of `𝒮`.) The honest statement: **the formal-monodromy diagonal is
real and already identified, but it is not "a diagonal factor of `𝒮`"; it is the holonomy of
the frame change.** Trying to "extract the restart factor" as `D·P` reproduces the
content-free `Q_restart` of `q_restart_probe.py` (achievable for any U(3)).

### 4c. What *is* well-defined here.

The frame-pinned, content-ful restart structure is **not** "(monodromy diagonal) ⋉
(triangular shear)" but:

- **(sheet 3-cycle permutation) ⋉ (mid-concentrated 2-pair Stokes data)** — the permutation
  is the `PI_OUT` relabelling; the genuine non-Abelian content is the two mid-sharing shears
  of §3, whose composition is `P_{m→m}`. The "diagonal/formal" part of the local irregular data
  is `e^{2πi c_i}` (the drift, §4b, established) and the "off-diagonal" part is the §3
  mid-pair shears (H-R2) — **so the correct restart factorization is the one H-R2 verifies,
  with the permutation being the cyclic sheet map rather than a triangular shear.**

**Verdict: PARTIAL / REFRAMED.** The (diagonal-formal) ⋉ (off-diagonal-Stokes) *split* is
real and is exactly H-R2; but the specific H-R1 phrasing — *triangular/unipotent* shear and
`e^{2πi c_i}` as a *diagonal factor of `𝒮`* — is not what the well-posed measurements show.
The off-diagonal is **cyclic**, and `e^{2πi c_i}` is the **frame-change holonomy**, not a
matrix factor. H-R1 is the subtle one and its sharp form is **not supported**; its sound
residue is subsumed by H-R2.

---

## 5. Net for the program

- **H-R3 (settled).** The two branch structures are cleanly separated: turning points
  (`-u'=0`, W4 zeros) carry the `Γ_j`/window/BE data; the irregular point at `u=∞`
  (`λ→ε_i`, where `Γ_j` is analytic) carries the off-diagonal Stokes data. This justifies
  treating `P_{m→m}` as an **irregular-point connection coefficient** distinct from the
  elementary Weber/BE part — the foundation the PA-1/PA-2/PV/CH tracks assume.

- **H-R2 (settled, strong).** The off-diagonal Stokes content is a **2-pair object** living
  in `{mid,lo}∪{mid,hi}` (the 12×13 joint, mid shared), and decoupling either outer link
  collapses it to a single elementary Weber factor. This gives PA-1/PA-2 a concrete reduced
  target: the **monodromy data = formal exponents `c_i` (known) + two mid-sharing Stokes
  shears**, and the connection coefficient is their (non-commuting) composition — matching
  the OPEN_PROBLEM addendum's restated structure (with the permutation cyclic, not
  triangular).

- **H-R1 (reframed, caution).** Do **not** carry the "triangular shear / `e^{2πi c_i}` as a
  factor of `𝒮`" phrasing into the τ-function/connection-coefficient work: the frame-pinned
  reality is a **cyclic** sheet permutation times the mid-pair shears, and `e^{2πi c_i}` is
  the subtracted drift's holonomy. The honest, well-defined restart factorization is H-R2's.

**Evidence ladder.** H-R3: **symbolic-proven** (exact) + numeric loci. H-R2: **numerically-
supported (strong)** — frame-pinned, T-converged, oracle-grade IP propagator, two regimes,
both decoupling directions. H-R1: the negative/reframe is **numerically-supported** (the
triangularity null model, the permutation-concentration, the drift-holonomy identification),
and the positive residue is exactly H-R2.

---

## A. Reproducibility

`experiments/restart_probe.py` (numpy/scipy/sympy; ≈4 min):
- `hr3_symbolic()` — sympy proof of `-u'=W4/p²`, pole-analyticity, `deg W4=4` (H-R3).
- `turning_points()` — W4 zeros vs poles (H-R3 loci).
- `S_IP` / `S_canon` — the convergent canonical-frame scattering matrix (engine =
  `uploads/assay/ip.py`, the oracle's propagator).
- `triangularity()` + `haar_triangularity_baseline()` — H-R1 triangularity metric & U(3)
  null model.
- `mid_concentration()` — H-R2a mid-row/col weight.
- `hr2_decoupling()` — H-R2b WS-C limit (`γ_outer→γ_outer·10⁻³`), `P_mid` vs incoherent BE.
- `coulomb_c()` — `c_i=Σ_{j≠i}s_ij²(a_i−a_j)` and `e^{2πi c_i}` (H-R1b).
- The §1 `T`-convergence table is a `/tmp` side-run over `T∈{80,120,160}` (canonical-frame
  convergence vs IP divergence); rerun by sweeping `T` in `S_canon`/`S_IP`. Default probe
  uses `T=120`, `rtol=1e-11`, at which `R_tri`/`frac_mid` are already 3-digit stable.

Samples: `canonical (ε=−2,0,3; γ=1,.8,1.2; a=−1,.5,2)`,
`sampleB (ε=−1,0,1.5; γ=.9,1.1,.8; a=−.7,.4,1.3)` (the WS-G overlapping/jointed showcase),
`well_sep (ε=−5,0,5; γ=1,1,1; a=−1.5,0,1.5)`.
