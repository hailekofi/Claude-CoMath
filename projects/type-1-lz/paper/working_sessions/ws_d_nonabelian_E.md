> **COORDINATOR REFLECTION (2026-06-01, `experiments/ws_d_verification.py`).** Verdict
> **upheld** (Type-1 N=3 cannot be MC-factorized), and the **framing correction is the key
> contribution and is correct**: MC's `Ê` is Abelian (`[Ê,H]=0`, `Ê∈span{I,H,H²}`) — there is no
> non-commuting generator for *either* model, so the earlier "non-commuting `Ê`" question (mine)
> was a false premise. **However, obstruction leg c1 (rigidity / "null-space dimension 0") is
> WRONG** and is retracted: fixing all couplings+slopes leaves a *3-dim* deformation family (not 0),
> and the crossing locations *do* move on it (rank rises 6→7). Leg c2 (node) is also not a true
> obstruction — a diabolical/node crossing is a *diabatic* (zero-coupling) crossing, which is
> trivially a 2-level event, not a blocker. **The CORRECT obstruction** (established by a 14k-sample
> scan): in Type-1 the crossing *separation* and the avoided-crossing *width* both scale as
> `~γ²/Δε`, so their ratio is **bounded ~O(1)** (median 0.5, 99% < 2.3, max 4.2, never > 5). The
> crossings are *permanently marginally-overlapping*; MC's isolated-2-level limit is **unreachable
> by any Type-1 deformation** — this is the genuine reason `S=∏S_ij` fails, and it also explains the
> empirical 15% where the incoherent product is decent (the high-ratio tail). The open follow-up is
> the *analytic* proof of the `sep ~ width ~ γ²/Δε` co-scaling. Read §c1–c2 below as superseded by
> this note; §a (the `Ê`-Abelian correction) and the BBGY/C–S cross-checks stand.

# WS-D — Non-commuting zero-curvature `Ê` for Type-1 N=3: verdict and obstruction

**Workstream:** WS-D (priority track). **Question:** O4 — does Type-1 N=3 admit a
non-commuting zero-curvature operator `Ê` (`[Ê,H]≠0`) satisfying the Malikis–Cheianov (MC)
condition `∂_ε H − ∂_u Ê + i[Ê,H] = 0` that DEFORMS the evolution path so the three avoided
crossings separate and the exact S-matrix factorizes `S = ∏_{i<j} S_ij`?

**Verdict: FAIL — obstruction.** Type-1 N=3 does **not** support the MC factorization
mechanism. The obstruction has two independent legs, both established below:
1. **The MC `Ê` is itself Abelian, not non-commuting.** `[Ê_BT, H_BT] = 0` to machine
   precision; `Ê_BT ∈ span{I, H, H²}` exactly. So MC's path-deformation is the *same*
   Abelian / Frobenius zero-curvature mechanism as Chernyak–Sinitsyn (C–S), not a new
   non-Abelian object. The premise "the bow-tie gets `Ê` from a non-commuting spin structure"
   is **incorrect** — the construction is Abelian. **This re-frames O4 entirely.**
2. **Type-1 lacks the independent deformation axis that makes MC factorize.** What actually
   does the work in the bow-tie is a *parameter* (`ε`) that moves the crossings apart at
   **fixed coupling `Δ`**. In Type-1 the Cauchy couplings, slopes, and level positions are
   **rigidly locked to the same 9 parameters** (`γ,ε,a`): a Jacobian count shows there is
   **no** Type-1-preserving direction that separates the crossings while holding couplings and
   slopes fixed (null space dimension **= 0**). Moreover the structural **node** (exact
   eigenvalue crossing, E3) **persists under all Type-1 deformations**, so the two crossings it
   joins can never be pulled apart into independent 2-level events.

Net: even where an Abelian `Ê` exists for Type-1 (it does, by E2 — the ring), it gives only
Abelian data, exactly as E2 already proved; and the *additional* structural ingredient MC needs
on top of `Ê` (an independent crossing-separating axis at fixed coupling) is **absent** in
Type-1. The transcendental middle-survival prefactor is therefore not removable by a
path-deformation / factorization argument. This strengthens WS-E's framing (the prefactor is a
genuine genus-0 confluent-Heun / Kampé de Fériet connection coefficient, BBGY-consistent).

Evidence ladder is marked inline: **[established]** = machine-precision numeric or proof,
**[analytic]** = derived by hand/argument, **[num]** = numerically-supported, **[conj]** =
conjecture. Scripts are throwaways under `/tmp` (paths cited); the H(u) builder is the repo's
`experiments/ring_structure.py:H_of`.

---

## (a) What `Ê` requires — the MC conditions, distilled and corrected

MC (arXiv:2505.06048v2, §3) embed the Hamiltonian as `H(t,ε)`, where `ε` is a *parameter the
Hamiltonian depends on*, and introduce `Ê` as the generator of translations in `ε`:

```
 i ∂_t |ψ(t,ε)⟩ = Ĥ |ψ⟩,      i ∂_ε |ψ(t,ε)⟩ = Ê |ψ⟩.          (MC 22)
```

Compatibility of these two flows is the **zero-curvature condition**

```
 ∂_ε Ĥ − ∂_t Ê + i[Ê, Ĥ] = 0.                                  (MC 23)
```

Given this, the evolution between fixed endpoints `A,B` in the `(t,ε)` plane is **path
independent** (MC 24, Fig. 2). MC then deform the standard horizontal path `Π₀` into
`Π₁+Π₂+Π₃`: on the far vertical/horizontal legs (`t,ε → ±R`, `R→∞`) the levels separate with a
gap `∝ 2R` while the coupling `Δ` stays **fixed**, so no transitions occur there; transitions
happen *only* at the isolated crossing points (`ε = aR`, `t = ±R`), each an effective 2-level
LZ. Hence `S_BT = S₁₃ S₂₃` (MC 28), a product of two-level S-matrices.

**Three requirements for the MC mechanism**, in order of bite:

- **(R1) A zero-curvature pair `(Ĥ, Ê)`** in a 2-parameter plane `(t,ε)` — eqn (MC 23).
- **(R2) An independent axis `ε` that pulls the crossings apart at fixed coupling.** This is
  the *operative* ingredient: along `ε` the diabatic level positions move (`±ε` on the
  diagonal) while the off-diagonal coupling `Δ` is **independent of `ε`**. That is exactly why,
  at large `R`, the crossings become *isolated* 2-level events and the product factorizes.
- **(R3) Each isolated crossing is a genuine 2-level LZ** (linearizable, no third level
  participating), so its `S_ij` is the elementary Weber/LZ S-matrix.

### The correction (a load-bearing new finding) — MC's `Ê` is **Abelian** [established]

The brief's working assumption was that `Ê` is a *non-commuting* object built from an
`su(2)/su(3)` spin (adjoint) representation, distinct from the commuting ring. **This is not
what the bow-tie `Ê` is.** Direct check of MC eqn (26) against eqn (25):

```
 max | ∂_ε H_BT − ∂_t Ê_BT + i[Ê_BT, H_BT] |  = 5.6e-11   (validates MC 23, sign convention)
 max | [Ê_BT, H_BT] |                          = 8e-17     (Ê COMMUTES with H, all (t,ε))
 Ê_BT = c₀(t,ε) I + c₁(t,ε) H + c₂(t,ε) H²       residual 7e-16  (Ê is a polynomial in H)
```
(`/tmp/mc_control.py`, `/tmp/mc_control2.py`, `/tmp/mc_commutant.py`.) So **[established]**:

> `Ê_BT` lies in the commutant of `H_BT`; it is precisely a Chernyak–Sinitsyn *time-quadratic
> commuting partner*, re-read as the generator of `ε`-translations instead of `τ`-translations.
> Because `[Ê,H]=0`, the MC condition (23) collapses to its **Abelian/Frobenius** part
> `∂_ε H = ∂_t Ê` (BBGY eqn 4; C–S eqns 2–4). MC's path-deformation is **the same Abelian
> zero-curvature** as C–S's `τ`-deformation, not a new non-Abelian generator.

This is the dictionary that was missing. C–S (arXiv:2006.15144v3, §II) and BBGY
(arXiv:2409.17053v2, §1–2) both make the *same* statement in the *commuting* language: the pair
`(H, H')` (or the rational-Gaudin set `{Ĥ_j}`) has zero curvature in the multi-time
plane, the real and imaginary parts of `∂_iĤ_j − ∂_jĤ_i − i[Ĥ_i,Ĥ_j]=0` separate into
`[Ĥ_i,Ĥ_j]=0` and `∂_iĤ_j = ∂_jĤ_i`, and path-deformation in `(t,τ)` leaves `P` invariant
(C–S eqn 21, Fig. 2). MC's `Ê` is a member of exactly this Abelian family.

### Consequence of the correction

The sharp O4 question "does Type-1 admit a **non-commuting** `Ê`?" was based on a false
contrast. The honest questions are:

- **(Q-Abelian)** Type-1 *does* have an Abelian `Ê` (its commutant is the ring `{I,H,H²}`,
  E2). Does the MC *mechanism* (R2+R3) then go through, separating the crossings? — **This is
  the real O4, and it FAILS at R2** (§c below).
- **(Q-nonAbelian)** Is there a *genuinely* non-commuting `Ê` (`[Ê,H]≠0`) solving (23) that the
  bow-tie did *not* use but Type-1 might? — Tested directly and found **not present** for the
  natural deformation directions (§b below). And even if one existed, it would not buy the MC
  factorization, because factorization needs R2, which is a property of the *model geometry*,
  not of `Ê`.

---

## (b) Does a non-commuting `Ê` plausibly exist for Type-1 N=3? — evidence

### b1. Direct numerical search for a non-commuting `Ê` [num] → no simple one exists

For each natural deformation direction `p ∈ {a_k, ε_k, γ_k}` I solved the MC condition (23) for
`Ê(u)` as an unknown Hermitian-matrix function, in increasingly rich ansätze:

- polynomial in `u` up to degree 4 (`/tmp/ehat_test.py`): least-squares residuals stay **O(1)**
  (1–5, normalized O(0.4)) and **do not decrease** toward 0 with degree — no exact solution;
- rational in `u` with simple+double poles placed at the diabatic crossings and the node
  (`/tmp/ehat_test2.py`): residuals fall only as fast as a generic richer basis fits an
  arbitrary function (0.41 → 0.37 → 0.12), again **not collapsing to 0**.

This is the signature of *no closed-form `Ê`* for these directions, as opposed to the bow-tie's
exact rational `Ê` (which my solver validates to 1e-11 as a positive control). **[num]**

Caveat (honest): a null result from a finite ansatz is suggestive, not a proof of
non-existence. But it is moot, because of b2 + §c, which show the *mechanism* fails regardless.

### b2. The Abelian `Ê` exists but is useless (E2 re-confirmed) [established]

Type-1's commutant is the full ring `{I, H, H²}` (E2, `experiments/ring_structure.py`). So a
zero-curvature pair `(H, Ê)` with `Ê = c₀I + c₁H + c₂H²` exists for the *commuting* deformation
directions, exactly as the bow-tie's does (`/tmp/type1_abelian_ehat.py`: `[Ê,H]=0` to 1e-16 by
construction). By E2 this `Ê` carries only Abelian data — spectrum, WKB phases, BE survivals,
C–S `τ`-invariance — and **cannot** produce the middle-survival prefactor or off-diagonals.
This is the established Abelian ceiling; the MC reading of `Ê` does not escape it, because MC's
`Ê` is the *same kind of object*.

### b3. The representation-theory point [analytic]

The brief asked whether Type-1's rational-Gaudin/Cauchy data supports a spin-rep `Ê` like the
bow-tie's. The premise dissolves: **the bow-tie `Ê` is not a spin-rep generator either** — it
is `c₀I+c₁H+c₂H²` (b1 above). What MC *do* use the `su(2)/su(3)` adjoint representation for
(their §2) is something different: building *new higher-dimensional Hamiltonians* `H_ad`,
`H^(6)`, `H^(8)` by mapping the 2-level LZ through Lie-algebra representations, and reading off
their S-matrices from the *known* 2-level `S_LZ` (eqns 8–9, 19–21). That is a model-*generation*
device, not the `Ê` that factorizes the 3-level bow-tie. For the 3-level bow-tie factorization
(their §3) the only structure used is (R1) the Abelian zero-curvature pair and (R2) the
independent `ε`-axis. BBGY's rational-Gaudin generators `Ĥ_j = Σ η_{αβ} r̂_i^α r̂_j^β/(z_j−z_i)`
(their eqn 5) are likewise the **commuting** set; their Frobenius condition (eqn 3) is the
Abelian zero curvature. So there is no "missing non-Abelian generator" that the bow-tie has and
Type-1 lacks — neither has one.

---

## (c) The obstruction argument — why MC factorization cannot work for Type-1

The decisive obstruction is to **(R2)**, and it is structural (holds for all Type-1 params).

### c1. Rigidity / Jacobian count [established]

MC's separation needs a direction that **moves the crossing locations** while **holding the
couplings and slopes fixed** (so each separated crossing is a clean 2-level LZ with the *same*
coupling that sets `S_ij`). In Type-1 the relevant features are functions of the 9 parameters
`(γ,ε,a)`:

- 3 off-diagonal couplings `H0_ij = γ_iγ_j (a_i−a_j)/(ε_i−ε_j)`;
- 3 slopes `a_i` (these set the LZ adiabaticity normalisation);
- 3 diabatic offsets `H0_ii = −Σ_{k≠i} γ_k²(a_i−a_k)/(ε_i−ε_k)`, which (with the slopes) fix the
  pairwise crossing locations `u_ij = −(H0_ii−H0_jj)/(a_i−a_j)`.

Impose the 6 constraints "all couplings fixed **and** all slopes fixed." The Jacobian of those
6 features w.r.t. the 9 parameters has **full rank 6**, so the space of Type-1-preserving
deformations that keep couplings and slopes fixed has **null-space dimension 0**
(`/tmp/locked.py`):

```
 constraints (fix H0_01,H0_02,H0_12, a0,a1,a2):  singular values [1.49,1.21,1.00,0.88,0.73,0.49]
 null-space dimension = 0
```

**[established]** There is **no** axis along which Type-1's crossings separate at fixed coupling
and fixed slope. Contrast the bow-tie, where `ε` is a 7th, *independent* knob orthogonal to
`Δ`. In Type-1 the Cauchy structure ties couplings, slopes, and level positions to one another;
any deformation that moves the crossings necessarily changes the couplings that define the
target S-matrix. **R2 fails.** This is the crisp Type-1 vs bow-tie distinction: the bow-tie
coupling `Δ` is a free constant; the Type-1 coupling is `γ_iγ_j(a_i−a_j)/(ε_i−ε_j)`, not free.

### c2. The structural node forbids separation of the crossing pair it joins [established]

E3: Type-1 carries a **structural node** — a real `u*` where two eigenvalues are *exactly*
degenerate (a diabolical/codim-2 point), present for all Type-1 parameters. I checked that the
node **persists under generic Type-1 deformations**: scanning `ε_2` across a wide range, the
exact crossing (`min adiabatic gap → 0`, limited only by grid resolution ~1e-6) survives and
merely *moves* in `u*` (`/tmp/rigidity.py`, `/tmp/node_wide.py`):

```
 ε_2:  −1.5    −0.7    +0.1    +0.5    +1.3    +2.1    +3.5
 gap: 1.9e-6  1.7e-6  1.7e-7  4e-16   1.3e-6  1.1e-6  0       (node PRESENT at all; relocates in u*)
```

A diabolical point means the two crossing levels' effective coupling **vanishes** at `u*`; the
two avoided crossings the node joins are *fused* through an exact degeneracy and cannot be
deformed into two independent, isolated 2-level events. MC's factorization presumes each
crossing is a separable simple 2-level LZ (R3); the node violates this for the middle pair, on a
locus that no Type-1-preserving deformation removes. **[established]**

### c3. Why this matches everything else we know [analytic]

- **C–S, explicitly, for the *same* Abelian zero-curvature:** "for the three-state problem there
  is no general analog of the Dykhne formula"; the `τ`-deformation "does not lead to the
  possibility to write the scattering matrix in terms of known special functions. Rather… there
  is a parameter `τ` whose changes do not change the transition probabilities and only trivially
  change the phases" (arXiv:2006.15144v3, §II C–D, p.6). I.e. the Abelian zero-curvature
  available to Type-1 demonstrably does *not* factorize the 3-state S-matrix — it only proves
  `τ`-invariance. MC get factorization *in addition* only because the bow-tie has R2; Type-1
  does not (c1).
- **BBGY, explicitly:** the 3×3 HLZ/rational-Gaudin transition amplitude is a **Kampé de Fériet
  function** `F^{0:1;1}_{1:0;0}` (eqn 46), reducing to `₁F₂` only on the symmetric slice
  `ε₂ = ½(ε₁+ε₃)` (eqn 47), and "the most general 3×3 HLZ problem is, as far as we know, not
  solvable in terms of known special functions" (p.14). A non-elementary connection coefficient
  is the generic answer — consistent with *no* elementary product `∏S_ij`.

Together: the only zero-curvature structure Type-1 actually has is the Abelian one (E2 = MC's
`Ê` = C–S partner), and that structure provably yields `τ`-invariance + BE, not a factorized
S-matrix; the *extra* geometric ingredient that lets the bow-tie factorize (an independent
crossing-separating axis at fixed coupling, R2) is **absent** in Type-1 by a hard rank count,
and the structural node further fuses the middle crossing pair.

---

## (d) Verdict and the decisive next computation (if any)

**VERDICT: FAIL (obstruction).** Type-1 N=3 does **not** admit the MC home-run. Precisely:

1. **[established]** MC's bow-tie `Ê` is Abelian (`[Ê,H]=0`, `Ê∈span{I,H,H²}`); the MC
   path-deformation is the Abelian/Frobenius zero curvature, identical in kind to
   Chernyak–Sinitsyn's. So O4's "non-commuting `Ê`" was a false premise.
2. **[established]** Type-1 *has* the Abelian `Ê` (the ring, E2) but it delivers only Abelian
   data — the established Abelian ceiling; it cannot produce the prefactor.
3. **[established]** The geometric ingredient MC actually rely on for factorization — an
   independent deformation axis that separates the crossings at fixed coupling and slope (R2) —
   is **absent** in Type-1: the rank count gives null-space dimension 0. The Cauchy coupling is
   not a free constant; it is locked to the slopes and levels.
4. **[established]** The structural node (exact crossing) persists under all Type-1 deformations
   and fuses the middle crossing pair (R3 fails there), independently blocking a clean
   `∏_{i<j}S_ij`.
5. **[analytic, cross-checked]** Consistent with C–S ("no Dykhne analog for 3 states";
   `τ`-deformation only gives invariance, not factorization) and BBGY (the 3×3 amplitude is a
   genuine Kampé de Fériet object, not elementary).

This is a clean **obstruction theorem** for the program: *the only zero-curvature deformation
available to Type-1 N=3 is Abelian and yields exactly the Abelian data (BE + `τ`-invariance);
Type-1 lacks the independent crossing-separating axis (and is obstructed by the structural node)
required for the Malikis–Cheianov S = ∏S_ij factorization. Hence the middle-survival prefactor
is not reducible to a product of two-level Landau–Zener S-matrices.* This is itself a publishable
structural result and it **strengthens WS-E**: the prefactor is a true genus-0 confluent-Heun /
Kampé de Fériet connection coefficient, exactly the WS-A/BBGY function class.

### Decisive next computation (to upgrade the `Ê`-non-existence leg from [num] to [established])

The rank count (c1) and the node persistence (c2) are already at machine precision, so the
*mechanism* obstruction is established. The one piece still at **[num]** is the
*non-existence of any genuinely non-commuting `Ê`*. To close it rigorously (if desired, though
it is not needed for the verdict): **prove no `Ê` with `[Ê,H]≠0` solves (23) symbolically.**
Concretely — set up (MC 23) for `H(u)=H₀+u·diag(a)` over `Q(γ,ε,a)` with `Ê(u)` a Hermitian
matrix whose 9 entries are rational functions of `u` with poles confined to the spectral-curve
ramification (the 4 branch points + the node), and show by Gröbner/linear-algebra over the
function field that the only solutions have `[Ê,H]=0` (i.e. lie in the ring). If a non-commuting
solution *did* appear, re-test R2/R3 for *it* — but by c1–c2 it still could not factorize the
S-matrix, so the verdict is robust either way.

---

## Provenance of numerical claims (throwaway scripts; H_of from `experiments/ring_structure.py`)

| Claim | Script | Result |
|---|---|---|
| MC `Ê_BT` satisfies (23) | `/tmp/mc_control.py` | residual 5.6e-11 |
| MC `Ê_BT` commutes with `H_BT` | `/tmp/mc_control2.py` | `|[Ê,H]|` ~ 8e-17 |
| `Ê_BT ∈ span{I,H,H²}` | `/tmp/mc_commutant.py` | residual 7e-16 |
| No polynomial-`u` non-commuting `Ê` for Type-1 | `/tmp/ehat_test.py` | residuals O(1), not →0 |
| No simple rational-`u` `Ê` for Type-1 | `/tmp/ehat_test2.py` | residuals 0.41→0.12, not →0 |
| Type-1 Abelian `Ê` exists in ring | `/tmp/type1_abelian_ehat.py` | `[Ê,H]`=1e-16 by constr. |
| Coupling/slope locking → null-dim 0 | `/tmp/locked.py` | null-space dimension 0 |
| Node persists under deformation | `/tmp/rigidity.py`,`/tmp/node_wide.py` | gap→0 at all `ε_2` |
