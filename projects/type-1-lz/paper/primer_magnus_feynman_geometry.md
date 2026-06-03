# Primer — Magnus expansions, Feynman graphs, and a geometric theory of U(3) selection

**Audience:** a physicist who wants the conceptual scaffolding behind the WS-GEOM workstream
(`paper/ws_geom_scope.md`). **Goal:** explain (i) what a Magnus expansion is, (ii) how it becomes a
Feynman-graph structure for our multistate Landau–Zener (MLZ) problem, and (iii) how that structure
fleshes out a *geometric* theory of how the Type-1 Hamiltonian selects its slice of `U(3)`. Self-contained;
no new claims (a pedagogical companion).

---

## (i) Magnus expansions — the log of an ordered exponential

The propagator of a linear time-dependent system `ψ'(t)=A(t)ψ` is the **time-ordered exponential**
`U(t)=\mathcal T\exp\!\int_0^t A(s)\,ds`. Time-ordering is the nuisance: `A` at different times need not
commute, so you cannot just exponentiate `\int A`.

Two ways to expand it:

- **Dyson series** (ordinary perturbation theory): `U = 1 + \int A + \iint_{t_1>t_2} A(t_1)A(t_2) + …` —
  a sum of *ordered products*. Truncating it **breaks unitarity** (a finite sum of products is not
  unitary).
- **Magnus expansion:** write `U=\exp\Omega(t)` with a **single** exponent
  `Ω = Ω_1+Ω_2+Ω_3+…`, where
  ```
  Ω_1 = ∫ A(t_1) dt_1,
  Ω_2 = ½ ∫∫_{t_1>t_2} [A(t_1),A(t_2)] dt_1 dt_2,
  Ω_3 = (nested double commutators), …
  ```
  Each term is built from **nested commutators** of `A`. The magic: if `A` is anti-Hermitian
  (`A^†=−A`), then every `Ω_n` is anti-Hermitian, so `U=e^Ω` is **exactly unitary at every truncation
  order**. Magnus expands the *logarithm* of the propagator; Dyson expands the propagator itself.

**Convergence / regime.** The Magnus series converges when `∫‖A‖\,dt` is small — i.e. when the coupling
is weak or the evolution short. For us (next section) `A=−i\tilde W` is the *adiabatic* coupling, small
when the energy gaps are large (the well-separated / generic regime). Deep overlap = large `∫‖\tilde W‖`
= the series must be resummed (this is exactly where our `σ` lives and why it is hard).

**Why we want it here.** (1) Unitarity at every order keeps `S∈U(3)` honest — a built-in check. (2) The
*log* structure separates the "free" part (order 0) from genuine multi-event corrections cleanly. (3) The
commutators encode interference (two non-commuting hops) directly, which is the physics of `\Pmid`.

## (ii) From the series to Feynman graphs — the MLZ adiabatic-`W` expansion

Put the MLZ problem in the **adiabatic frame**: `ψ=Φ(u)χ`, `Φ` the (Type-1: `a`-independent) instantaneous
eigenbasis, so `iχ'=[D(u)−iW(u)]χ` with `D=\mathrm{diag}(E_i)` and `W_{ij}=⟨φ_i|φ_j'⟩` the off-diagonal
**derivative coupling** (the geometric seed). Strip the diagonal by going to the doubly-rotating frame:
the dressed coupling is
```
   \tilde W_{ij}(u) = W_{ij}(u)\; e^{\,i\int^u (E_i−E_j)\,du'} .
```
Now read the Dyson/Magnus series as **diagrams on three sheets** (the three adiabatic levels):

- **Lines** = the three adiabatic sheets `i=1,2,3` (the levels of the spectral curve).
- **Vertices** = insertions of `\tilde W_{ij}`: a **hop** from sheet `j` to sheet `i` at "time" `u`,
  weighted by the geometric amplitude `W_{ij}(u)` (the non-adiabatic coupling).
- **Propagators (edges)** = the accumulated **dynamical phase** `e^{i\int(E_i−E_j)}` between consecutive
  hops — a pure phase set by the energy gaps.
- **A term of order `n`** = a sum over all ordered sequences of `n` hops (paths that walk across the
  sheets), each weighted by `∏(vertices)×∏(edge phases)`, integrated over the ordered hop-times.

This is a genuine Feynman expansion: *sheets are particles, `W` is the interaction vertex, the dynamical
phase is the propagator, and the transition amplitude is the sum over histories.* The dictionary:

| diagram element | MLZ object |
|---|---|
| order 0 (no vertices) | adiabatic following → the **directed-cycle permutation** `S_0` |
| 1 vertex | a single non-adiabatic hop → leading LZ/Stückelberg amplitude (→ BE factors) |
| 2 vertices (the `Ω_2` commutator) | two interfering hops → the leading `\Pmid` interference |
| resummed all-orders | the full holonomy = `σ` (no closed form) |

**The key structural fact for Type-1.** `W` is `a`-**independent** (the family shares the eigenbasis,
R1). So the **vertices are fixed geometric data**, and the **`a`-dependence lives entirely in the edge
phases** (the dynamical actions). The diagram expansion therefore *separates the geometry (vertices) from
the elementary kinematics (phase edges)* — the cleanest possible organization of the holonomy.

## (iii) Toward a geometric theory of how H selects `U(3)`

The transition matrix `S∈U(3)` is the **holonomy** of the adiabatic connection (`W` is a non-abelian
Berry connection on the eigenbundle over the `u`-line, or over the spectral curve `Σ`). "How does `H`
select its slice of `U(3)`?" becomes "what is the holonomy representation of this bundle?" The Magnus/
Feynman structure makes the answer *layered*, and each layer is a geometric object:

1. **Order 0 = a topological datum.** `S_0` is the permutation by which the adiabatic energy ordering
   **reverses** between `u=∓∞`. That permutation (our directed 3-cycle) plus the **spinor double-cover
   signs** `δ_j=±1` (the sign of each eigenvector, which flips on encircling a branch point of `Σ`;
   NOMENCLATURE `Γ_j=δ_j(…)`) is a **discrete/topological label** — the sector of `U(3)` the dynamics
   lives in, fixed by the *geometry of the eigenbundle*, not by any dynamical detail.
2. **Higher orders = the smooth dressing inside that sector.** The `W`-vertices and phase-edges deform
   `S` within the topological sector. This is the *analytic* content — the BE factors, the interference,
   and (resummed) `σ`. The image region (facet A) is the set swept out as `(γ,ε,a)` vary; its **boundary**
   is where a hop amplitude `W_{ij}` vanishes — a level decouples — exactly the classifier locus (R9).

So the **geometric theory** reads:
> `U(3)` selection `=` (a **topological permutation/spinor sector** from the eigenbundle monodromy on the
> spectral curve) `×` (a **smooth holonomy dressing** from the `W`-vertices and dynamical-phase edges),
> with the **decoupling locus as the boundary** of the swept region.

The Feynman expansion is the *bridge*: it turns the abstract statement "S is the holonomy of `W` on `Σ`"
into a computable, layered object — topology explicit at order 0, analytic dressing in the higher orders,
and the hard transcendental `σ` quarantined as the all-orders resummation. A geometric theory needs
exactly this: a way to say *which* part is rigid/topological (and therefore knowable in closed form) and
*which* part is the irreducible transcendental — and to compute the boundary between them.

**What a candidate "invariant" would be, in this language.** Not a closed form for `σ`, but a
*characteristic-class*-type object: a winding/permutation label, the parity of the spinor sector, or a
relation among the BE actions `δ_{ij}` and the cross-ratio `χ` that is conserved along the family or that
cuts the boundary. The WS-GEOM numerics hunt for one; the Magnus/Feynman structure tells you *where to
look* (order 0 for the topology; the boundary `W_{ij}→0` for the analytic edge).

---

### One-paragraph summary
A Magnus expansion writes the propagator as `e^{Ω}` with `Ω` a sum of nested commutators of the coupling —
unitary at every order and convergent in the weak-coupling (generic) regime. For Type-1 MLZ the coupling
is the `a`-independent adiabatic seed `W`, and the expansion becomes a Feynman series on three sheets:
geometric `W`-vertices joined by dynamical-phase edges, summed over hop-histories. Order 0 is the
directed-cycle permutation set by the eigenbundle's energy-reordering and spinor signs (a **topological**
selection of the `U(3)` sector); the higher orders are the smooth `W`-holonomy dressing (the **analytic**
content, with `σ` as its non-closed resummation); the boundary of the swept region is the decoupling
locus. That layered decomposition — topology `×` holonomy dressing, with an explicit boundary — is the
skeleton of a geometric theory of how the Hamiltonian selects its slice of `U(3)`, and it is entirely
`σ`-free except for the one resummed datum we have proven is irreducible.
