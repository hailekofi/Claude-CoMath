# WS-GEOM M3 T1 — the node$\to$directed-cycle theorem  [analytically-derived]

**Type:** `physics-derivation` (with symbolic/numeric verification). **Date:** 2026-06-03.
**Companion:** `experiments/ws_geom_m3_t1.py`. **Builds on:** M1 (order-0 skeleton), R2 (genus-0 curve),
R3 (universal real node), NOMENCLATURE (the three orderings).

---

## Theorem (T1)
For generic Type-1 N=3, the **order-0 adiabatic-following permutation** — the leading term of the
adiabatic-`W` Magnus expansion (R18) — when the eigenframe is **continued through the unique real node**
(R3), is a **directed 3-cycle**, with orientation fixed by the spectral data. It is *not* the
energy-sorted extreme-swap transposition.

## The three orderings (state the theorem in the right basis)
Per NOMENCLATURE, three orderings must be kept distinct:
1. **ε-ordering** — the **fixed basis/diabatic labels** `i` with `ε_0<ε_1<ε_2`. The propagator `S`, the
   matrix `P`, and the standard basis vectors `e_k` are **ε-indexed**.
2. **Slope-ordering** — `lo,mid,hi = argsort(a)`, the BE/dynamical *roles*.
3. **Energy-rank** — the instantaneous rank of `E_k(u)`; **end-dependent** (it reverses between
   `u=∓∞`), and therefore **never a channel label**.
The permutations below act on the ε-indexed labels; their *content* is dictated by the slope ordering.

## Proof
**Step 1 — the bare reordering is the extreme-SLOPE swap.** As `u→±∞`, `H(u)=H_0+uA→uA`, so
`E_k(u)→a_k u` and the eigenvectors `→ e_k` (diabatic basis; `A=diag(a)`). Hence the instantaneous
energy ordering is **slope-ascending** at `u→+∞` and **slope-descending** at `u→−∞`: the energy rank of
each diabatic state *reverses* between the two ends. Energy-sorting connects rank-`r` at `−∞` to rank-`r`
at `+∞`, so in diabatic (ε-indexed) labels it realizes the **full order-reversal** of the slope ordering.
For `N=3` the order-reversal `(\mathrm{lo},\mathrm{mid},\mathrm{hi})\mapsto(\mathrm{hi},\mathrm{mid},\mathrm{lo})`
is exactly the **transposition of the two extreme-slope states**, `τ_ext=(\mathrm{lo}\ \mathrm{hi})`
(middle fixed). This depends only on `a_k` being distinct — hence it is **universal** (every sample).
*Verified:* energy-sort perm `=(2,1,0)` on canonical + 4 random draws (in the slope basis; `(2,1,0)` is
the extreme transposition).

**Step 2 — the node is an adjacent-rank transposition.** By Owusu–Wagh–Yuzbashyan (R3) there is a unique
real `u_*` where two instantaneous eigenvalues *cross*. A real level crossing is between **energy-adjacent**
levels (they coincide, hence are neighbors just off the crossing). Continuing the eigenframe *through*
`u_*` follows each eigenvector along its branch across the crossing, which **interchanges the two adjacent
energy ranks** relative to energy-sorting: a transposition `τ_node` of an adjacent pair. *Verified:* the
continued frame differs from the energy-sorted one by exactly one adjacent swap.

**Step 3 — composition is a 3-cycle.** The continued order-0 permutation is therefore
`π_{\rm cont} = τ_ext ∘ τ_node`. Two transpositions that **share exactly one index** compose to a
**3-cycle** (the elementary identity `(a\,b)(a\,c)=(a\,c\,b)`). Since `τ_ext=(\mathrm{lo}\ \mathrm{hi})`
and `τ_node` is an adjacent pair — `(\mathrm{lo}\ \mathrm{mid})` or `(\mathrm{mid}\ \mathrm{hi})`, each
sharing one index with `(\mathrm{lo}\ \mathrm{hi})` — the product is a directed 3-cycle. *Verified:*
`π_{\rm cont}=(1,2,0)` (a 3-cycle) on all samples, and `π_{\rm cont}=τ_ext∘τ_node` holds exactly. ∎

**Step 4 — orientation.** Which 3-cycle (the two orientations) is fixed by *which* adjacent pair the node
swaps, i.e. by the node's sheet labels at `(u_*,E_*)` — a function of the spectral data. On the sampled
family the node uniformly swaps the lower adjacent pair, giving the cycle
`\mathrm{lo}\to\mathrm{hi}\to\mathrm{mid}\to\mathrm{lo}`. *(The 3-cycle conclusion is general; the uniform
orientation is verified on the family, not separately proved — a minor owed item.)*

## Why this answers "energy-sort `=(2,1,0)` on every sample" and the ε-ordering question
The extreme-swap is the diabatic image of the **asymptotic energy-order reversal**: `E_k\sim a_k u`
flips its ordering between `∓∞`, and for three levels a full reversal swaps the ends and fixes the middle.
It is the **extreme-SLOPE** pair (`argsort(a)[0]\leftrightarrow argsort(a)[2]`) that swaps — *expressed in
the ε-indexed basis*. In the verification samples the slopes happened to be ε-monotone (`a_0<a_1<a_2`
alongside `ε_0<ε_1<ε_2`), so the slope-extremes coincided with ε-indices `0,2` and the permutation read
`(2,1,0)`; for a generic sample with non-ε-monotone slopes it swaps ε-indices `argsort(a)[0]` and
`argsort(a)[2]`, which need not be `0` and `2`. **This is precisely why the strict ε-ordering matters:**
the *basis* is ε-indexed (fixed), the *roles* (which states swap) are slope-ordered, and the
*energy* ordering is the derived, end-reversing quantity that produces the swap — and *because* it reverses
between the two ends, energy rank can never label a channel (NOMENCLATURE). T1 is the cleanest illustration
of the convention: a geometric statement that is only correct when the three orderings are not conflated.

## Evidence
- `[analytically-derived]` Steps 1–3 (asymptotics + crossing adjacency + the transposition-composition
  identity); the 3-cycle conclusion is general.
- `[numerically-verified]` energy-sort `=` extreme-swap, node `=` adjacent swap, continued `=` 3-cycle `=`
  `τ_ext∘τ_node`, on canonical + 4 random samples (`ws_geom_m3_t1.py`).
- `[numerically-supported / owed]` the uniform **orientation** across the family (verified, not proved).
- Regime caveat (M1): the order-0 skeleton equals the *full* oracle permutation only in the adiabatic
  regime; for diabatic/small-action samples the oracle sits near the identity. T1 concerns the order-0
  skeleton, not the full dynamics.

## T2 — two-vertex reachability  [analytically-derived + numerically-verified]
*Claim:* in the two free coordinates `{P_mm,b}` (R15), the image closure of `Φ` touches **exactly two** of
the six Birkhoff vertices — the identity `(1,0)` and the node-selected directed cycle `(0,1)` — and no
other. (The six permutations map to `{P_mm,b}` as: identity`→(1,0)`; directed cycle`→(0,1)`; extreme-swap
`(lo\,hi)→(1,1)`; reverse cycle, `(lo\,mid)`, `(mid\,hi)` `→(0,0)`.)
*Proof.* A vertex (a permutation matrix `P`) requires every transition probability to be `0` or `1`, which
occurs only in a deterministic limit.
- **Diabatic limit** `γ→0`: `H_0→0` (its entries `∝γ_iγ_j,γ_k^2`), so `H(u)→uA=\mathrm{diag}(a)u` is
  diagonal in the diabatic basis — no transitions, `S→\mathbb 1`, `{P_mm,b}→(1,0)` (identity).
- **Adiabatic limit** `γ→∞` (gaps `→∞`): evolution is adiabatic, the system follows the instantaneous
  eigenstates, and by **T1** the node-continued following permutation is the directed 3-cycle, so
  `{P_mm,b}→(0,1)`.
- **No other vertex is a limit point.** The extreme-swap `(lo\,hi)` `(1,1)` is the *energy-sorted* order-0,
  but the node (R3, always present) converts it to the 3-cycle, so it is never the adiabatic limit; it is
  not the diabatic limit (identity) either. The reverse cycle is excluded because the node fixes the
  orientation (T1, Step 4). The remaining transpositions are limit points of neither limit. Finite
  couplings give interior (non-permutation) points interpolating between `(1,0)` and `(0,1)`. ∎
*Verified* (`ws_geom_m3_t2t3.py`): `gscale=0.03→(1.000,0.000)`, `gscale=2.0→(0.024,1.000)`; the `(0,0)`
and `(1,1)` corners stay `>0.5` away throughout (and M2a's cloud: min-dist `~0.003,0.005` to the two
reached vertices, `~0.27` to the forbidden ones).

## T3 — edge-selection boundary $=$ decoupling locus  [analytically-derived + numerically-verified]
*Claim:* the `{P_mm,b}` image boundary is the decoupling locus (R9), edge-resolved: a **middle** coupling
`→0` ⇒ `P_mm→1`; an **extreme** coupling `→0` ⇒ `b→0`.
*Proof (decoupling block structure).* When `γ_k→0`, level `k` decouples and `H` block-reduces to a 2-level
problem on the other two `⊕` a trivial spectator `k`.
- **Middle decoupling** (`γ_mid→0` kills both `s_{mid,lo},s_{mid,hi}∝γ_mid`): the middle is a spectator,
  surviving with probability `1`, so `P_mm=P[mid,mid]→1` — the `{P_mm,b}` point lands on the `P_mm=1`
  edge (the residual lo–hi 2-level fixes `b∈[0,1]`).
- **Extreme decoupling** (`γ_lo→0` or `γ_hi→0`): the directed cycle `lo→hi→mid→lo` needs circulation
  through all three levels; removing an extreme breaks it, so the `lo→hi` amplitude `b=P[hi,lo]→0` — the
  point lands on the `b=0` edge.
Global bounds `b≥0`, `P_mm≤1` hold always (probabilities), so these are genuine boundary edges; the
decoupling locus traces them, edge-selected by *which* level decouples. This is the geometric realization
of R9 (**elementary ⇔ a level decouples**). ∎
*Verified* (`ws_geom_m3_t2t3.py`, M2c): `γ_mid→0 ⇒ P_mm=1.0000`; `γ_lo→0 ⇒ b=0.0000`; `γ_hi→0 ⇒ b=0.0000`.

## Status / hand-off
T1 [keystone], **T2, and T3 are closed** (analytically-derived + numerically-verified) ⇒ the geometric
theory (R18--R19) now has a *fully derived* topological skeleton: the leading `U(3)` permutation is a
node-selected directed 3-cycle (T1); the reachable set is the two-vertex region between identity and that
cycle (T2); its boundary is the decoupling locus, edge-selected (T3). The **only remaining M3 item is T4**
— the `δ_j` spinor double-cover label, which would also settle the *orientation* of the cycle (T1 Step 4)
in general — a stretch, possibly partial. With T1--T3 done, the rigid (σ-free) half of the
`U(3)`-selection theory is derived; `σ` (the analytic dressing) remains the sole irreducible remainder.
