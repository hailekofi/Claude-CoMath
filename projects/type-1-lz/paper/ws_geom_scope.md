# WS-GEOM — scope: the geometric selection of U(3) sections (adiabatic-W Magnus + image-region)

**Type:** `physics-derivation` + `physics-numerics` scoping roadmap (NOT a finished result). **Date:**
2026-06-02. **Companion primer:** `paper/primer_magnus_feynman_geometry.md`. **Status:** scoped, not
executed.

---

## 0. The question (refined intent)

Beyond the closed-form question (settled: no closed form; `σ` irreducibly transcendental, R9/R11/R17),
a *geometric* question is wide open and barely touched by the present work:

> **How does the Type-1 Hamiltonian select its slice of `U(3)`?** The map `Φ:(γ,ε,a)→S∈U(3)` (the
> propagator) lands the transition data in a particular, structured part of `U(3)`. What is that
> structure, and what selects it?

After end-rephasing the physical content is in reduced-`U(3)` (dim 4), and with ~4 essential physical
parameters `Φ` is generically `4→4` — so the "section" is **not a thin slice** but a full-dimensional
**region with structured boundary and a topological label**. The question splits cleanly:

- **Facet A — the image region.** *Which* part of the unistochastic body is occupied, and how? Known
  bias: toward the **directed 3-cycle** (R16) — a cluster near a permutation vertex of the Birkhoff
  polytope, biased by slope ordering, with boundary the *classifier locus* (a level decouples, R9).
- **Facet B — the holonomy / topological selection.** *How* the path `u↦U(u,−∞)` over the
  spectral-curve double cover lands there: the holonomy of the `a`-independent seed `W` on the
  eigenbundle, and the spinor/sign sector (the `δ_j=±1` form-factor signs, the eigenframe monodromy
  around the branch points and the node).

**Liberating point:** both facets are answerable **without closing `σ`**. `σ` is the *value* of one
holonomy datum; the *geometry* (image region, bundle, topology, perturbative skeleton) is real structure
regardless. This workstream targets that geometry.

## 1. Part I (`physics-derivation`): the adiabatic-`W` Magnus expansion — facet B

**Setup.** Adiabatic frame `ψ=Φ(u)χ`, `Φ=[φ_1\,φ_2\,φ_3]` the `a`-independent shared eigenbasis (R1).
`iχ'=[D^{(a)}(u)−iW(u)]χ`, `D^{(a)}=diag(E_i^{(a)})`, `W_{ij}=⟨φ_i|φ_j'⟩` (anti-Hermitian, off-diagonal,
`a`-independent, geometric). Go to the doubly-rotating frame (strip `D`): the dressed coupling
`\tilde W_{ij}(u)=W_{ij}(u)\,e^{i\int^u(E_i−E_j)}`. Then
`S = Φ(+∞)\,[\,\text{T-exp}\!\int(−\tilde W)\,]\,Φ(−∞)^†` (up to the diagonal Stark phases that drop from
`P`).

**The expansion (Magnus / Dyson in `\tilde W`).**
- **Order 0 (`W=0`, perfect adiabatic following).** `S_0 = Φ(+∞)·\text{diag(phases)}·Φ(−∞)^†`. Since
  `Φ(±∞)` are the standard-basis frames *ordered by adiabatic energy*, and that ordering **reverses**
  between `∓∞`, `|S_0|^2` is a **permutation** — the **directed 3-cycle** (R16's leading structure),
  fixed by the slope ordering. This is the *topological* skeleton of the section.
- **Order 1.** `Ω_1=−\int \tilde W\,du`: a single inter-sheet hop dressed by the dynamical phase. Gives
  the leading non-adiabatic (LZ/Stückelberg) amplitudes; `|·|^2` recovers the BE factors to leading
  order.
- **Order 2.** `Ω_2=\tfrac12\iint[\,\tilde W(u_1),\tilde W(u_2)\,]`: two hops / the first interference —
  the genuine `\Pmid`-type two-path content at leading order.

**Tasks.** (a) Derive `Ω_1,Ω_2` explicitly; identify the directed-cycle order-0 term and the first hop +
interference corrections; keep `W` (geometric, `a`-independent) separate from the phase edges
(elementary, `a`-dependent). (b) Validate **term by term** against the oracle in the **generic regime**
(`δ≲0.25`, where `\tilde W` is small): `S_0`, `S_0+Ω_1`, `S_0+Ω_1+Ω_2` should approach the oracle with
the expected order in the gap. (c) Map the result onto the Feynman-graph structure (primer §ii): sheets =
lines, `W` = vertices, dynamical phases = propagators.

**Honest limits.** The series does **not** resum in deep overlap (`δ≳0.5`) — same boundary as WS-O3/R17;
that is fine, the geometric *skeleton* (order 0 + low orders) is the target, and the resummation boundary
*is* the deep-overlap core. Unitarity: Magnus (exp of anti-Hermitian `Ω`) preserves `S∈U(3)` at each
order — a built-in consistency check the Dyson truncation lacks.

## 2. Part II (`physics-numerics`): map the image region + invariant hunt — facet A

**Tasks.** (a) Generate large `(γ,ε,a)→P` data (cheap: oracle ~s/sample; fast solver / WS-O3 instant).
(b) Characterize the **image region** in reduced-`U(3)`: the directed-cycle bias (distance to the nearest
permutation vertex vs `δ`), the effective dimension, and the **boundary** (test that it coincides with the
classifier/decoupling locus, R9). (c) **Invariant hunt** — caged symbolic regression (the *structure
detector*, NOT a surrogate, NOT a closed-form-finder): search for a combination of
`{Pmid,b,δ_{ij},χ,c_i,…}` that is **conserved across the family** or that **cuts the boundary**. Any
candidate → hand to `physics-derivation` to prove. Expected candidates are *topological/boundary* objects
(a winding/permutation label, the spinor sign sector, a relation among BE actions), not a codimension
constraint (the image is full-4-dim).

**Honest limits.** ML earns its keep ONLY as a conjecture-generator for facet A; do not let it become a
surrogate (we have the exact oracle) or a formula-finder (no closed form exists, R9/R11). The image may be
"as structured as we already know" — if the invariant hunt finds nothing beyond BE+DS+χ, that is a
first-class negative (the section's analytic structure is exhausted; only the topological label remains).

## 3. Milestones
1. **M1 (derivation):** `Ω_1,Ω_2` explicit; order-0 = directed cycle proven; Feynman-graph dictionary
   written. *Gate:* term-by-term convergence to the oracle in the generic regime.
2. **M2 (numerics):** image-region map; boundary = classifier locus confirmed; directed-cycle-bias law.
3. **M3 (invariant):** one candidate invariant from the caged regression, **proven** or precisely refuted.
4. **M4 (synthesis):** the geometric statement — *the section = (order-0 topological permutation from the
   eigenbundle/spinor data) dressed by the `W`-holonomy, with boundary the decoupling locus* — written as
   a paper section, `σ`-free.

## 4. What success / partial / failure mean (evidence-ladder honest)
- **Success:** an explicit order-0+1+2 geometric skeleton matching the oracle generically + a proven
  invariant/topological label ⇒ a genuine geometric theory of the selection (`analytically-derived` +
  `numerically-supported`).
- **Partial:** the skeleton closes but no new invariant ⇒ the selection's *analytic* content is BE+DS+χ,
  the *topological* content is the order-0 permutation/spinor sector (still a real result).
- **Failure:** the Magnus series is uninformative even at low order (unlikely — order 0 is exact and
  geometric) ⇒ record and stop.

## 5. Risks
1. The geometric "hard part" still bottoms out at `σ` (the resummed holonomy) — but the *skeleton* and the
   *topology* do not, so the deliverable is honest and `σ`-free.
2. Phase/branch bookkeeping (the `e^{i\int(E_i−E_j)}` edges, the spinor signs `δ_j`) is fiddly; track the
   double-cover signs carefully (NOMENCLATURE `Γ_j=δ_j(...)`).
3. Invariant hunt may return nothing new (first-class negative, see §4 partial).

## 6. Reusable assets
deformation_family_probe.py / ring_structure.py (the seed `W`, shared eigenbasis), skeleton_two_
transcendentals.py (the directed cycle, `{Pmid,b}`), ws_o3_uniform.py (the semiclassical leading terms,
the generic-regime validation), oracle.py / num_S12.py (gold + cheap data), ws_g_stokes_graph.md (the
Stokes/cover geometry), the demoted WKB/selector "common cover & spinor lift" material (now re-approached
with the seed-`W` understanding). Primer: `paper/primer_magnus_feynman_geometry.md`.
