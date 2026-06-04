# WS-GEOM Milestone 3 — scope: FORMALIZE the topological sector (a derivation, not an invariant hunt)

**Type:** `physics-derivation` scoping (roadmap for review; not executed). **Date:** 2026-06-03.
**Builds on:** M1 (`ws_geom_m1.md` — node-selected directed-cycle skeleton), M2 (`ws_geom_m2.md` —
two-vertex image, edge-selection boundary, **full effective dimension**), R2/R3 (genus-0 curve + universal
node), R9 (classifier). **Status:** SCOPED.

---

## 0. Reframing (forced by M2d)

The original M3 was a *symbolic-regression invariant hunt*. **M2d kills that target**: the effective
dimension of the lift `{P_mm,b,δ_ij,χ,c_i}` is *full* (the only PCA zeros are the known elementary
relations, incl. `Σc_i=0`), so **there is no hidden analytic invariant to find** — the section's analytic
content is exhausted by BE + double-stochasticity + `χ`. M3 is therefore **not a search**; it is a
`physics-derivation` task: **prove, to physics-referee standard, the topological statements that M1+M2
established numerically.** All targets are `σ`-free.

## 1. Targets (each a topological/limit statement, to be derived)

### T1 [keystone] — the node$\to$directed-cycle theorem
*Claim:* for generic Type-1 N=3, the adiabatic-following permutation `S_0` (the order-0 Magnus term),
**continued through the unique real node** (R3), is the **directed 3-cycle** with orientation fixed by the
slope ordering — *not* the energy-sorted transposition.
*Strategy:* eigenbundle monodromy on the genus-0 spectral curve `Σ` (R2). Away from the node the adiabatic
levels do not cross (von Neumann–Wigner), so energy-rank is preserved and the bare reordering between
`u=∓∞` is the **extreme-swap transposition** `(lo\,hi)`. The **node** is a real transverse double point
where two specific sheets cross; continuation through it applies a **transposition of those two sheets**.
Show: (i) *which* two sheets the node swaps (from the node's sheet labels at `E_*`), and (ii) that the
composition (extreme-swap)$\circ$(node-swap) is the **3-cycle**. This is a finite algebraic/topological
statement; verify symbolically on the canonical + random samples (M1's overlap-continuation already
exhibits it).
*Evidence target:* `analytically-derived` (composition of two explicit transpositions), numerics in hand.

### T2 — two-vertex reachability
*Claim:* the image closure touches **exactly** the identity `(1,0)` and the directed cycle `(0,1)`; the
reverse cycle and the transpositions are unreachable.
*Strategy:* the only permutation endpoints are the two free limits — the **diabatic limit** (couplings
`→0`: `S→` identity in the diabatic basis) and the **adiabatic limit** (gaps `→∞`: `S→` the
node-selected cycle, T1). Argue no other permutation arises as a limit point (any finite coupling gives a
non-permutation interior point; the two permutation corners are the only boundary vertices). Cross-check
M2a (min-dist `~0.27` to the forbidden corners).
*Evidence target:* `analytically-derived` (two limits) + `numerically-supported` (M2a).

### T3 — the edge-selection boundary $=$ decoupling locus
*Claim:* the `{P_mm,b}` image boundary is the decoupling locus (R9), edge-resolved: an **extreme**
coupling `→0` ⇒ `b→0`; the **middle** coupling `→0` ⇒ `P_mm→1`.
*Strategy:* when coupling `γ_k→0`, level `k` decouples and `S` block-reduces. Derive: (a) extreme
decoupling removes the directed circulation through that level ⇒ chirality `b=P_{hi,lo}→0` (the cycle
needs all three links); (b) middle decoupling ⇒ the middle is a spectator ⇒ `P_mm→` its single-crossing
survival `→1` as `δ_mid→0`. Tie to R9 (elementary ⇔ a level decouples) and R15's affine map. Verify
(M2c, coordinator-checked: extreme`→b=0`, middle`→P_mm=1` to `±0.000`).
*Evidence target:* `analytically-derived` from the block structure + `numerically-supported` (M2c).

### T4 [stretch] — the spinor double-cover label
*Claim:* the eigenvector sign sector `δ_j=±1` (NOMENCLATURE) and the node's action on it constitute the
remaining topological datum; the eigenframe monodromy around the branch points of `Σ` is the double cover.
*Strategy:* track `δ_j` around the four complex branch points + the node of `Σ`; relate the global sign
sector to the order-0 permutation (T1). This is the deepest piece and **may remain partial** — scope it as
a stretch; a clean statement of the `δ_j` monodromy representation would complete the topological label,
but a partial/heuristic result is an acceptable M3 outcome.

## 2. Method & honesty
- This is **topology/limit analysis**, not a small-parameter expansion — no controlled-vs-uncontrolled
  ladder; instead each claim is a finite algebraic/topological statement with a symbolic check.
- **No `σ` content** anywhere: M3 formalizes the *rigid* sector; the transcendental dressing is untouched
  (R9/R11/R17 stand).
- Validation: every target has M1/M2 numerics already in hand; the derivations must reproduce them.
- Honest outcome ladder: **T1 closed** ⇒ the geometric theory has a *proven* topological core (the headline
  win). T1+T2+T3 ⇒ the full rigid skeleton is derived. T4 partial is acceptable.

## 3. Milestones
- **M3a:** T1 (node$\to$cycle) — derive + symbolic-verify. *Gate:* (extreme-swap)$\circ$(node-swap) `=`
  3-cycle as an identity over the sample family.
- **M3b:** T2 + T3 — the two-vertex limits and the edge-selection boundary, derived from the
  diabatic/adiabatic limits and the decoupling block structure.
- **M3c [stretch]:** T4 — the `δ_j` double-cover monodromy.
- **M3d:** fold the proven statements into the paper's geometric-theory section (R18--R19 → upgrade tags
  from NS to AD where T1--T3 close).

## 4. Deliverables (when executed)
`paper/working_sessions/ws_geom_m3.md` (the derivations T1--T4 with symbolic checks + evidence tags) and any supporting
`experiments/ws_geom_m3_*.py` (the permutation-composition / branch-monodromy checks). No `σ`, no closed
form claimed.

## 5. What M3 does NOT do
It does **not** search for an analytic invariant (M2d shows there is none), does **not** touch `σ`, and
does **not** claim a closed form. It upgrades the M1/M2 *numerical* topological observations to *derived*
statements — completing the rigid (knowable) half of the `U(3)`-selection theory and leaving `σ` as the
explicitly-quarantined irreducible remainder.
