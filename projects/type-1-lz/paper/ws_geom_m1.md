# WS-GEOM Milestone 1 — the adiabatic-W Magnus skeleton of S in U(3), order by order

**Type:** `physics-derivation` + `physics-numerics` (validated result, with honest limits). **Date:**
2026-06-03. **Companion:** `paper/ws_geom_scope.md` §1, `paper/primer_magnus_feynman_geometry.md`.
**Code:** `experiments/ws_geom_magnus.py` (reproducible; reuses `type1` + the gold `oracle.py`).
**Conventions:** NOMENCLATURE.md — states ε-indexed (`ε_0<ε_1<ε_2`); BE roles by slope `argsort(a)`;
`P[x,j] = prob(x→j)`.

---

## 0. What M1 establishes (one paragraph)

The Type-1 N=3 scattering matrix `S ∈ U(3)` has an explicit, layered geometric skeleton:
**order 0 is the directed-3-cycle permutation** fixed by the adiabatic energy reordering between
`u=∓∞` — a *topological* datum, the `a`-independent eigenbundle monodromy — and the **higher Magnus
orders are the smooth `W`-holonomy dressing** of that permutation. We derive `Ω_1, Ω_2` explicitly in the
`a`-independent dressed coupling `W̃`, give the Feynman dictionary (sheets = lines, `W` = vertices,
dynamical phase = edges), and validate term by term against the gold oracle. **Order 0 = the directed
cycle is proven exactly** (a clean permutation matrix, matching the oracle's dominant structure to ~1e-4).
The term-by-term convergence holds in the **adiabatic regime**: at every adiabatic sample the error
decreases with Magnus order (`e0 > e1 > e2`), and the whole expansion converges to the oracle as the
system is made more adiabatic. The convergence is **modest, not fast** (order 1 already captures the
leading non-adiabatic correction; order 2 adds a smaller interference gain to a ~5×10⁻³ floor), and we
report honestly **where it does not work**: the small-coupling *diabatic* corner (order 0 = the cycle is
the wrong leading term there) and the *deep-overlap* core (the series does not resum — the same boundary
as WS-O3/R17, where `σ` lives).

---

## 1. Setup — the adiabatic frame and the dressed coupling

The commuting family (R1, `ring_structure.py`) shares one `a`-independent eigenbasis `Φ(u)=[φ_1 φ_2 φ_3]`,
`H(u)=H_0+uA` with `A=diag(a)`. In the adiabatic frame `ψ=Φ(u)χ`:

```
i χ' = [ D(u) − i W(u) ] χ,    D=diag(E_i(u)),    W_ij = ⟨φ_i|∂_u φ_j⟩   (anti-Hermitian, off-diagonal).
```

`W` is the **geometric seed** — `a`-independent (R1/F1). We use the **exact Hellmann–Feynman** form (no
finite differencing of eigenvectors):

```
W_ij(u) = ⟨φ_i| A |φ_j⟩ / (E_j(u) − E_i(u)),    i ≠ j,    W_ii = 0.
```

Going to the doubly-rotating (dressed) frame strips `D`; the dressed coupling is

```
W̃_ij(u) = W_ij(u) · exp( i ∫_{-T}^{u} (E_i − E_j) du' ).
```

The diabatic scattering matrix is

```
S = Φ(+T) · [ T-exp ∫ (−W̃) du ] · Φ(−T)^†          (diagonal Stark phases drop from P = |S|^2).
```

**Numerical care (the error-prone part).** The eigenframe is built by **overlap continuation**, NOT raw
`eigh` sorting: at each `u` the eigenvectors are matched to the previous `u` by maximum overlap (so labels
do not reorder at avoided crossings) and sign-continued. `W` uses the Hellmann–Feynman quotient.
The dynamical phase is the cumulative trapezoid of `E_i−E_j`. (`experiments/ws_geom_magnus.py:adiabatic_frame, dressed_coupling`.)

---

## 2. The Magnus expansion, order by order

`U_dress = exp(Ω)`, `Ω = Ω_0 + Ω_1 + Ω_2 + …` (Magnus = log of the ordered exponential; each `Ω_n`
anti-Hermitian ⇒ `S ∈ U(3)` at **every** truncation order — a built-in consistency check Dyson lacks).

**Order 0 (`W̃ = 0`).**
```
S_0 = Φ(+T) Φ(−T)^†.
```
`|S_0|^2` is a **permutation matrix** — the adiabatic-following map. Because the standard-basis eigenframes
at `∓T` are ordered by adiabatic energy and that ordering **reverses** between the two ends (it swaps the
two extreme-slope levels; the middle stays rank-2), the permutation is the **directed 3-cycle**. In the
oracle convention (`PI_IN=(0,1,2)`, `PI_OUT=(2,0,1)`) it is `x → PI_OUT[x]`, i.e. `(0→2, 1→0, 2→1)`.

**Order 1 (one dressed hop).**
```
Ω_1 = − ∫_{-T}^{+T} W̃(u) du.
```
A single inter-sheet hop dressed by the accumulated dynamical phase — the leading non-adiabatic
(LZ/Stückelberg) amplitude. `S^{(1)} = Φ(+T) exp(Ω_1) Φ(−T)^†`, `P^{(1)} = |S^{(1)}|^2`.

**Order 2 (two hops / first interference).**
```
Ω_2 = (1/2) ∫∫_{u_1 > u_2} [ W̃(u_1), W̃(u_2) ] du_1 du_2.
```
The commutator of two non-commuting hops — the genuine two-path interference (the leading `P_mm`-type
content). `S^{(2)} = Φ(+T) exp(Ω_1 + Ω_2) Φ(−T)^†`, `P^{(2)} = |S^{(2)}|^2`.

Verified anti-Hermiticity (machine precision): `||Ω_1+Ω_1^†|| ~ 1e-16`, `||Ω_2+Ω_2^†|| ~ 1e-16`, and
`||exp(Ω_1+Ω_2)^† exp(Ω_1+Ω_2) − I|| ~ 1e-15` — `S` stays in `U(3)` at every order.

---

## 3. The Feynman dictionary (primer §ii, made explicit)

The Dyson/Magnus series is a sum over hop-histories on three sheets (the three adiabatic levels):

| diagram element | MLZ object | `a`-dependence |
|---|---|---|
| **line** `i` | adiabatic sheet / level `E_i(u)` of the spectral curve | eigenvectors `φ_i` are `a`-independent |
| **vertex** `W̃_ij(u)` | a hop `j→i` at "time" `u`, amplitude `W_ij(u)=⟨φ_i|A|φ_j⟩/(E_j−E_i)` | **`W` is `a`-independent — geometric** |
| **edge (propagator)** | dynamical phase `exp(i∫(E_i−E_j))` between consecutive hops | **carries ALL the `a`-dependence** (via `E_i ⊃ u·a`) |
| order 0 (0 vertices) | adiabatic following → **directed-cycle permutation `S_0`** | topological; `a`-independent |
| order 1 (1 vertex) | one hop → leading LZ/Stückelberg amplitude → BE factors | phase edge is `a`-dependent |
| order 2 (the `Ω_2` commutator) | two interfering hops → leading `P_mm` interference | two phase edges |
| resummed all orders | the full holonomy = `σ` (no closed form, R9/R11/R17) | — |

**Key structural fact (Type-1).** The *vertices `W` are fixed geometric data* (`a`-independent, R1); the
*`a`-dependence lives entirely in the edge phases*. The diagram expansion therefore cleanly **separates
the geometry (vertices) from the elementary kinematics (phase edges)** — the cleanest organization of the
holonomy. We verified `W` is `a`-independent directly (`ring_structure.py`/`deformation_family_probe.py`, F1).

---

## 4. The decisive term-by-term validation

### 4.1 Order 0 = the directed cycle (DECISIVE — passes)

Canonical adiabatic anchor `eps=(-2,0,3)`, `gam=(1,0.8,1.2)`, `a=(-1,0.5,2)`:

```
|S_0|^2 (P[x,j]) =  [[0.0001 0.0001 0.9998]      order-0 permutation (x→j):     (2, 0, 1)
                     [0.9997 0.0001 0.0001]      directed 3-cycle (convention): (2, 0, 1)   MATCH
                     [0.0001 0.9998 0.0001]]      oracle dominant per row:       (2, 0, 1)   MATCH
```

`|S_0|^2` is a permutation matrix to ~1e-4, equals the directed 3-cycle, and equals the oracle's dominant
entry per row. **Order 0 = the directed cycle is established.**

### 4.2 The M1 gate — error vs Magnus order AND adiabaticity (adiabatic regime)

Hold `eps=(-2,0,3)`, `a=(-1,0.5,2)`; scale `gam` up → wider avoided crossings → more adiabatic (order-0
cycle is the correct skeleton, `e0` small and shrinking). Errors `e_n = max|P^{(n)} − P_oracle|`,
**T-averaged** over `T∈{70,80,90}` to cancel the endpoint Stark oscillation (see §5). Oracle at `T=120`.

| `gam`-scale | `δ_max` | `Λ=∫‖W̃‖` | `e0` | `e1` | `e2` | `e1/e0` | `e2/e1` | `e0>e1>e2` |
|---|---|---|---|---|---|---|---|---|
| 1.0 | 0.240 | 2.99 | 3.25e-1 | 7.29e-2 | 5.87e-2 | 0.225 | 0.805 | **True** |
| 1.3 | 0.685 | 2.98 | 1.29e-1 | 1.46e-2 | 1.04e-2 | 0.113 | 0.711 | **True** |
| 1.6 | 1.573 | 2.96 | 6.47e-2 | 1.11e-2 | 1.05e-2 | 0.171 | 0.951 | **True** |
| 2.0 | 3.840 | 2.94 | 2.69e-2 | 6.66e-3 | 6.63e-3 | 0.248 | 0.996 | **True** |

**Gate verdict (both required directions hold):**
- **(a) error decreases with Magnus ORDER** at every sample: `e0 > e1 > e2` in all four rows. ✓
- **(b) the whole expansion converges to the oracle as the system becomes more adiabatic**: `e0`
  (1/100 of which is the cycle-defect) falls `3.3e-1 → 2.7e-2`, `e1` falls `7.3e-2 → 6.7e-3`,
  `e2` falls `5.9e-2 → 6.6e-3` — all monotone toward 0. ✓

**Accuracy and honest qualifiers.** Order 1 captures the **leading** non-adiabatic correction (`e1/e0 ~
0.1–0.25`: a ~5–10× reduction over order 0). Order 2 (interference) adds a **further but smaller** gain
(`e2/e1 ~ 0.7` at moderate adiabaticity, → 1 at strong coupling): the absolute error floors at ~5×10⁻³,
because (i) the dressed action `Λ=∫‖W̃‖ ≈ 3` is **not small** across this sweep (it is nearly
scale-invariant — the bare‐coupling shrinkage is offset by the phase de-suppression), so the natural
Magnus small parameter never gets tiny here; and (ii) the finite-`T` endpoint Stark tail (§5) sets a
residual floor. So M1 demonstrates a **genuine, validated, order-improving geometric skeleton** with a
**clean order-0 = directed cycle and a clear order-1 leading correction**, but it is **not** a rapidly
convergent `δ²`-type series — the interference order is a modest, floor-limited correction. This is the
honest content of the M1 gate.

### 4.3 Honest limit — the diabatic corner (`δ→0` via small `gam`)

Scaling `gam` **down** at fixed `eps,a` also sends `δ→0`, but it **narrows** the crossings → the passage
is **diabatic**, the oracle `P → identity`, and **order 0 = the directed cycle is the WRONG leading term**
(`e0 → 1`). The series still converges term-by-term (`e2 < e1`) but it must spend its low orders rotating
the cycle back toward the identity. Representative oracle `P` at `gam`-scale 0.45:

```
P_oracle = [[0.899 0.051 0.050]      ->  I-distance 0.10  (near-diabatic, NOT the directed cycle).
            [0.065 0.905 0.031]
            [0.036 0.045 0.919]]
```

Full diabatic-corner table (`gam`-scale down, `T=70`):

| `gam`-scale | `δ_max` | `e0` | `e1` | `e2` | I-dist(`P_oracle`) |
|---|---|---|---|---|---|
| 0.45 | 0.0098 | 9.55e-1 | 4.21e-2 | 3.19e-2 | 0.101 |
| 0.30 | 0.0019 | 9.92e-1 | 2.09e-2 | 7.32e-3 | 0.021 |
| 0.20 | 0.0004 | 9.98e-1 | 2.11e-2 | 1.94e-3 | 0.004 |

`e0 → 1` (the directed cycle is the wrong leading term), `e2 < e1` (the series still converges), and
`P_oracle → I` (near-diabatic).

**Lesson (load-bearing).** `δ` small does **not** by itself mean "directed cycle dominates": there are two
opposite small-`δ` limits. The directed-cycle skeleton is the **adiabatic** one (wide crossings). The M1
gate is therefore stated in the adiabatic regime (scaling `gam` up), not the naive `δ→0` corner. The scope
note's phrasing "`δ≲0.25` ⟺ `W̃` small ⟺ directed cycle" conflated the two; the precise statement is:
**order 0 = directed cycle in the ADIABATIC regime, where the dressed-hop content (not `δ`) is small.**

### 4.4 Honest limit — deep overlap (no resummation)

Deep overlap (large `δ` via merged/strong crossings, `Λ` large) is where the series **does not resum** —
the same boundary as WS-O3/R17, the home of `σ`. The order-≥1 corrections **stop reducing the error and
even increase it**:

| sample | `δ_max` | `Λ` | `e0` | `e1` | `e2` | `e2<e1` |
|---|---|---|---|---|---|---|
| deep1 (merged eps) | 5.73 | 2.93 | 7.86e-3 | 1.75e-2 | 1.74e-2 | True but **e1 > e0** |
| deep2 (strong gam) | 14.6 | 2.90 | 7.14e-3 | 7.66e-3 | 7.55e-3 | True but **e1 > e0**, flat |
| sampleB (overlap) | 1.08 | 3.19 | 7.32e-2 | 1.34e-2 | 1.55e-2 | **False (e2 > e1)** |

The non-convergence signature is unambiguous: at deep1/deep2 **order 1 makes the error WORSE than order 0**
(`e1 > e0`) and order 2 does not recover; at sampleB **order 2 increases the error over order 1**
(`e2 > e1`). `P^{(1,2)}` are not controlled approximations there. This is the expected, honest deep-overlap
core: the geometric *skeleton* (order 0 + the low-order dressing) is the M1 target, and the resummation
boundary **is** the deep-overlap transcendental `σ`. (Note: small `e0` at deep1/deep2 is incidental — the
directed cycle happens to be near `P_oracle` for those parameters — not a sign of convergence; the
diagnostic is that the *corrections fail to systematically reduce* the error.)

---

## 5. The endpoint Stark oscillation (a bookkeeping caveat, resolved by T-averaging)

The dressed coupling has a Fresnel-type tail: `W_ij ~ 1/u²` decays but the dynamical phase
`∫(E_i−E_j) ~ (a_i−a_j)u²/2` keeps the off-diagonal entries from decoupling, so
`W̃_ij ~ e^{i c u²}/u²` and `Ω_1, Ω_2` pick up an `O(sin(cT²)/T²)` **endpoint oscillation** at finite `T`.
At a single `T` the per-order errors oscillate (e.g. `e1` ran 2.5e-3, 1.9e-2, 2.0e-3, 1.2e-2 over
`T=40,60,80,100`). **Averaging over a few incommensurate `T` cancels the oscillation** and exposes the
genuine truncation-free Magnus error (the §4.2 table is T-averaged). The physical fix beyond M1 would be
the canonical log-`T` Stark-phase subtraction (`num_S12.S12_canonical_phase`); for the M1 *probability*
gate, T-averaging suffices. This is exactly the phase/branch-bookkeeping risk flagged in `ws_geom_scope.md`
§5.2 — it bit, and it is resolved/quarantined, not faked.

---

## 6. What M1 delivers for the geometric theory

1. **The topological layer is explicit and exact.** `S_0 =` directed 3-cycle `(2,0,1)` — the
   `a`-independent eigenbundle reordering between `u=∓∞`. This is the sector of `U(3)` the dynamics lives
   in, fixed by geometry, not by any dynamical detail. (Proven to ~1e-4; a clean permutation.)
2. **The analytic dressing is order-resolved and validated.** `Ω_1` (one hop, BE-leading) and `Ω_2` (two
   hops, leading interference) are explicit in the `a`-independent seed `W`, anti-Hermitian (unitary at
   every order), and they reduce the error to the oracle term-by-term in the adiabatic regime; the whole
   expansion converges to the oracle as adiabaticity grows.
3. **The Feynman separation is verified.** Vertices `W` = geometry (`a`-independent); edges = dynamical
   phase (all the `a`-dependence). This is the bridge that turns "`S` is the holonomy of `W` on `Σ`" into a
   computable, layered object — topology at order 0, analytic dressing in the higher orders, `σ` quarantined
   as the all-orders resummation.
4. **The boundaries are honest.** Not the naive `δ→0` (that is the diabatic corner, wrong skeleton); the
   deep-overlap core does not resum (where `σ` lives). The convergence in between is order-improving but
   modest (Λ ≈ 3, floor ~5e-3), so M1 is a *validated skeleton*, not a fast series.

**Evidence-ladder tags.**
- `[analytically-derived]` — `Ω_0` (directed cycle from the reordering), `Ω_1`, `Ω_2` explicit; the
  vertex/edge `a`-(in)dependence; anti-Hermiticity ⇒ unitarity at every order.
- `[numerically-verified, machine-precision]` — order-0 = directed-cycle permutation = oracle dominant
  (~1e-4); `Ω_n` anti-Hermitian (~1e-16); `W` `a`-independence (R1/F1).
- `[numerically-supported]` — the M1 gate: term-by-term decrease `e0>e1>e2` and convergence to the oracle
  as adiabaticity grows (T-averaged, four samples). Modest convergence, ~5e-3 floor — reported as such.
- `[honest-negative / boundary]` — the diabatic corner (cycle wrong) and the deep-overlap non-resummation;
  the finite-`T` endpoint Stark oscillation (resolved by T-averaging).

**What M1 does NOT claim.** No closed form for `σ` (R9/R11/R17 stand). No fast/`δ²` convergence. No
resummation in deep overlap. The directed cycle is the *adiabatic* skeleton only. Order 2 is a modest,
floor-limited interference correction, not a precision term at these couplings.

---

## 7. Reproduce

```
cd projects/type-1-lz/experiments
python3 ws_geom_magnus.py     # versions banner; order-0 permutation vs directed cycle;
                              # the T-averaged M1 gate table; the diabatic corner; the
                              # deep-overlap non-convergence; the unitarity check.
```
Requires `numpy`, `scipy`, and the sibling `oracle.py` (→ `uploads/assay`).
