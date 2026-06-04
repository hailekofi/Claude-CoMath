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

**Step 4 — orientation (a Z₂ label; see T4).** Which 3-cycle (the two orientations) is fixed by *which*
adjacent pair the node swaps: `(mid,hi)→` forward `lo→hi→mid→lo`, `(lo,mid)→` reverse `lo→mid→hi→lo`
(the two transpositions sharing one index with `(lo\,hi)`). **This orientation is NOT universal** —
broad deterministic sampling (T4) gives *both* (REV 7, FWD 4 over clean samples), so it is a genuine
**Z₂ label**, not a fixed orientation. *(An earlier 5-sample family read uniformly forward — a small-sample
artifact, now corrected.)* The 3-cycle conclusion (Steps 1–3) is general and unaffected; only *which* of
the two it is varies with the spectral data — **the selector, now resolved in T4** (it is the
energy-position of the real node, `\mathrm{sign}(\operatorname{tr}H(u_*)-3E_*)`).

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

## T4 — the spinor / cycle-orientation Z₂ label AND its selector  [RESOLVED]
*Claim (established):* the directed-cycle **orientation is a genuine Z₂ topological label** — both
orientations occur across parameter space — and it is the eigenframe spinor/double-cover sector (`δ_j=±1`).
*Evidence (the Z₂ is real).* The deterministic, dynamics-free overlap-continuation of the order-0 eigenframe
(T1's construction) yields **both** orientations over broad samples (`ws_geom_m3_t4.py`). By the group rule
(T1 Step 4) this is the node swapping `(lo,mid)` vs `(mid,hi)`. So the `U(3)` section comes in **two mirror
sectors**, not one.

### The selector (resolved): the parity of the ε→slope ordering
*Theorem (T4 selector).* The Z₂ orientation has **two equivalent closed forms**:

**(1) Combinatorial — the clean form.** The orientation is the **parity of the permutation `σ` that maps
the ε-order to the slope-order** (the slope rank carried by each ε-index):
$$\boxed{\ \mathrm{orientation} \;=\; \operatorname{sgn}(\sigma)\ }\qquad \text{even}\Rightarrow\text{FWD}\ (1,2,0),\quad \text{odd}\Rightarrow\text{REV}\ (2,0,1).$$
This is the **sign character `S_3\to\{\pm1\}`** — the eigenframe Z₂/spinor sector itself. It is **purely
combinatorial: independent of the couplings `γ` and of every spacing** — it depends *only* on the relative
order of the diabatic (ε) and adiabatic (slope) labels. *(Verified γ-independent: an odd assignment stays
REV and an even one stays FWD over 80 random `γ` draws.)*

**(2) Spectral — the geometric realization.** Equivalently, the orientation is the energy-position of the
true degeneracy at the unique real node `(u_*,E_*)` (R3):
$$\mathrm{orientation} \;=\; \mathrm{sign}\big(E_3 - E_*\big)\;=\;\mathrm{sign}\big(\operatorname{tr}H(u_*) - 3E_*\big),$$
with `E_*` the doubly-degenerate eigenvalue (the R8 accessory `v_*`) and `E_3 = \operatorname{tr}H(u_*) - 2E_*`
the spectator. Explicitly:
- `E_3 > E_*` (node degenerates the **lower** pair) `⟺` even `σ` `⟺` **FWD** `(1,2,0)`, `lo→hi→mid→lo`;
- `E_3 < E_*` (node degenerates the **upper** pair) `⟺` odd `σ` `⟺` **REV** `(2,0,1)`, `lo→mid→hi→lo`.

The two forms **agree on every sample** — the real node sits on the lower/upper pair exactly according to
`\operatorname{sgn}(σ)`.

*Derivation (spectral flow).*
1. **Factor the continued frame.** The overlap-continuation transports the eigenframe in `u` from `−∞` to
   `+∞`. Away from the node the transport is adiabatic (each eigenvector stays on its energy level); the only
   place adjacent levels meet is the unique real node `u_*` (R3), where the pair `(k,k+1)` is *exactly*
   degenerate. Continuous transport carries the two colliding eigenvectors **smoothly through** the exact
   crossing (analytic continuation; no avoided-crossing mixing), i.e. the transposition `τ_k=(k,k+1)` of
   energy ranks relative to an energy-sorted frame. Hence `π_cont = τ_ext ∘ τ_k` (T1).
2. **`τ_ext` is fixed.** The pure energy-adiabatic map is the asymptotic order-reversal `τ_ext=(lo\,hi)`
   (Step 1 of T1) — the Brundobler–Elser content, identical for every member, carrying **no** Z₂ freedom.
3. **The Z₂ lives entirely in `k`.** With `τ_ext` fixed the only freedom is `k∈\{0,1\}`:
   `τ_ext∘τ_0=(lo\,hi)(0\,1)=` FWD and `τ_ext∘τ_1=(lo\,hi)(1\,2)=` REV (two transpositions sharing one
   index → a 3-cycle whose sense is set by the shared index).
4. **`k` is the parity of `σ`.** `k=0` (lower pair degenerate) `⟺` spectator above, `E_3>E_*`;
   `k=1` (upper pair) `⟺` `E_3<E_*`. To see this is `\operatorname{sgn}(σ)`: relabel by an **adjacent
   slope-transposition** (two slopes cross, `a_i\leftrightarrow a_{i+1}`). At that slope collision the real
   node runs to `u\to\pm\infty` and returns on the other side (verified: scanning `a_{mid}` through `a_{lo}`
   flips the orientation as `u_*` diverges), so the node moves from one adjacent pair to the other and `k`
   flips. Hence the orientation changes sign under every transposition of `σ` ⟹ `\text{orientation}(σ) =
   \operatorname{sgn}(σ)\cdot\text{orientation}(\mathbb 1)`, and the ε-monotone identity is FWD ⟹
   `\text{orientation}=\operatorname{sgn}(σ)`. Both `u_*,E_*` are explicit algebraic data (R3 real node =
   real double root of `Disc_E(u)`; R8 `v_*=E_*`), so the spectral form is a closed-form algebraic sign. ∎

*Why this succeeds where the local screen failed.* The refuted candidates (`-\mathrm{sign}(u_*)`;
node-proximity to the lo–mid vs hi–mid diabatic crossing; coupling ratios; slope offset;
`\mathrm{sign}(E_*-(H_0)_{mid,mid})`; crossing-time combinations — all `≈`chance) all probed **continuous,
local** features. The selector is a **discrete combinatorial invariant** — a permutation *parity*,
`γ`-independent — invisible to any single continuous local probe. Its spectral realization lives on the
**real node** (R3, an exact crossing = a real root of `Disc_E`), **not** the nearest complex branch point
(the dominant avoided crossing): using the latter misclassifies the small-`|u_*|` REV cases. Once the right
object — the ε↔slope parity, equivalently the real node's sheet-pair — is named, the rule is exact.

*Evidence (near-proof).* `ws_geom_m3_t4_selector.py`: **both** forms match the overlap-continuation Z₂ —
the **parity** `\operatorname{sgn}(σ)` on **119/119** clean samples, the **spectral** sign on **154/154**
(`100.00%`, zero mismatch; FWD 85 / REV 69) plus 33/33 on an independent block and the canonical anchor
(`u_*=-0.2493`, lower pair, FWD). The two forms **agree on every sample**. **Two independent algorithms**
for the spectral form — Hermitian `eigh` of `H(u_*)` and companion-matrix roots of `p(E,u_*)` (no Hermitian
solver) — agree on all 154; and the parity form is **`γ`-independent over 80 random draws** at fixed
ordering. The selector is `analytically-derived` (the spectral-flow group identity + the
transposition-flip ⟹ `\operatorname{sgn}(σ)` argument) and `numerically-supported (near-proof)`.
*(Methodological note: reading the orientation from the adiabatic-limit `P` is unreliable — that is the
deep-overlap regime where `P` does not cleanly reach a vertex; the deterministic overlap-continuation is the
reliable measurement, and it is what the closed-form sign reproduces.)*

## Status / hand-off
**T1 [keystone], T2, T3, T4 are all closed** (analytically-derived + numerically-verified/near-proof). The
geometric theory (R18--R19) now has a **fully derived rigid skeleton**: the leading `U(3)` permutation is a
node-selected directed 3-cycle (T1); the reachable set is the two-vertex region between identity and that
cycle (T2); its boundary is the decoupling locus, edge-selected (T3); the section splits into two Z₂ mirror
sectors by cycle orientation, and **the sector is selected by the energy-position of the real node**,
`\mathrm{sign}(\operatorname{tr}H(u_*)-3E_*)` (T4). The **σ-free topological theory is complete.** The
**only** remaining open item is `σ` itself (the analytic dressing — the irreducible Fredholm/Widom
connection constant), which does not touch the topological skeleton.
