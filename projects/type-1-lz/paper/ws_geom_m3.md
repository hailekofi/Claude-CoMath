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

## Status / hand-off
T1 [keystone] **closed** (analytically-derived, verified) ⇒ the geometric theory (R18) now has a *derived*
topological core: the leading `U(3)` permutation is a node-selected directed 3-cycle. Remaining M3:
T2 (two-vertex reachability), T3 (edge-selection boundary = decoupling locus), T4 (the `δ_j` double-cover
orientation/label — which would also settle Step 4's orientation in general).
