# Where Type-M travels, v2: outstanding physics problems the structure can attack

**Type:** `physics-intuition` (Generation) + research overview. **Date:** 2026-07-09.
**Supersedes/extends:** `where_this_travels.md` (2026-06-03). That note was written as the Type-1
N=3 *solvability* quest was set down and is Type-1-N=3-centric. This v2 answers the broader
question — *the Type-M classification is completely parameterized, exactly solvable, and carries
nontrivial geometry; what outstanding physics problems can that structure help solve?* — and folds
in what was established since (R31–R42): the Möbius/arity-4 invariant theory, the nodal-cubic
genus dichotomy ("Type-M ⇔ genus M−1" [NS/conjecture at general M]), the wildness theorem
(`solvable ⟺ rigid`), the two-spheres ε↔a duality flag (R39), and the anchored factorization
S = W₊D₊·C·D₋W₋⁻¹. All export claims tagged; prior art checked where flagged.

---

## 0. What the structure actually exports (the asset list, updated)

An N×N family H(u) = H₀ + uA of type M (Owusu–Wagh–Yuzbashyan 2009 / Owusu–Yuzbashyan 2011,
J. Phys. A 44, 395302) gives, *simultaneously*:

1. **Complete parameterization** — the integrable matrices form an *explicitly parameterized
   algebraic variety* inside matrix space (coordinates {γ, ε, a}-type data). You can sample it,
   project onto it, and move exactly along or exactly transverse to it. No other nontrivial
   integrable class offers this.
2. **Exact spectra on an explicit curve** — eigenvalues are algebraic on the spectral curve;
   Type-1 is rational (nodal, genus 0 at N=3, R33); genus grows with type
   ["Type-M = genus M−1", NS at N=3, conjecture generally]. DOS, gap structure, discriminant
   geometry are exact.
3. **Theorem-level exact crossings** — integrability forces real parameter-crossings violating
   Wigner–von Neumann (OWY 2009; R3). The nodes are algebraically located and *protected by the
   commuting partner*, not by symmetry.
4. **Shared eigenbasis / exact adiabatic connection** — the commuting family shares one
   a-independent eigenframe; the derivative connection W(u) = ⟨φᵢ|φⱼ′⟩ is a single geometric
   object for the whole family (R1), rational on the curve in the λ-frame (R34/R36).
5. **The wildness theorem** — `closed-form driven dynamics ⟺ rigid connection` (Katz);
   Type-1 N=3 is the *generic non-rigid* rank-3 point, and the solvable MLZ zoo (bowtie,
   Demkov–Osherov, …) is exactly the rigid stratum. A classifier, constructive in both directions.
6. **Invariant coordinates** — P = F(be×3, χ) on a 4-dim moduli with full Möbius covariance
   (R31/R32, N=3): the natural "shape/scale" variables for any question you put to the family.

Two kinds of target below: **(I) outstanding problems in physics** where these assets give real
leverage, ranked; **(II) deep/structural long shots**. Fields where the family already lives
(Richardson–Gaudin pairing, central spin, sectors of 1D Hubbard/Heisenberg) make the exports
physical, not toy.

---

## I. Ranked problem pool

### P1 — A quantitative theory of integrability breaking (the "quantum KAM laboratory") ★flagship
**The outstanding problem.** There is no analytic theory of the Poisson→Wigner–Dyson crossover:
how weak a perturbation destroys integrability, how the threshold scales with system size
(claims of exponentially small thresholds are debated), and — more basically — there is no good
*definition of the distance from integrability* for a given Hamiltonian. "Glimmers of a quantum
KAM theorem" (Brandino–Caux–Konik) is the state of the art; it is a recognized open problem.
**Type-M leverage (three independent handles, no other class has any of them):**
 (a) *Distance-to-integrability order parameter.* Because the integrable set is an explicit
     variety, `dist(H, Type-M variety)` is a computable number (constrained least squares over
     {γ,ε,a}). Correlate it with spectral/AGP chaos measures → a candidate order parameter the
     field lacks. [conjecture; probe cheap]
 (b) *The type index M is a discrete integrability dial.* Prior art (CHECKED): Scaramazza–
     Shastry–Yuzbashyan, PRE 94, 032106 (2016) — typical integrable matrices are Poisson iff the
     number of commuting partners n ≳ log N; below that, level repulsion *within the exactly
     solvable class*. So the Type-M hierarchy already *crosses the Berry–Tabor boundary while
     staying exactly solvable* — a solvable model OF the transition itself, essentially
     unexploited. The curve geometry (genus growing with M) is the natural candidate for the
     analytic mechanism. [grounded + conjecture]
 (c) *Exact transverse perturbation theory.* KAM-type analysis needs exact knowledge of the
     tangent/normal decomposition at the integrable point; the parameterization provides it.
**First probes.** ⟨r⟩ vs M at fixed N (extend `ws_type_M_levelstats.py`, which found the sharp
Poisson→GOE crossover under +εGOE); dist-to-variety vs ⟨r⟩ on a breaking path; crossover scale
vs (N, M) against Rosenzweig–Porter scaling.
**Risk.** SSY 2016 owns the statics baseline — the *new* content must be the geometric mechanism,
the order parameter, and the crossover dynamics, not the Poisson statement.

### P2 — Exact adiabatic gauge potentials: benchmark for the AGP chaos program + exact counterdiabatic driving
**The outstanding problem.** The AGP norm is the most sensitive known probe of quantum chaos
(Pandey–Claeys–Campbell–Polkovnikov–Sels, PRX 10, 041017 (2020)) and the generator of
counterdiabatic (shortcut-to-adiabaticity) protocols — but it is essentially numerics-only:
there is no nontrivial family where the AGP is exact, so threshold claims (e.g. exponentially
small breaking scales) lack an analytic anchor, and CD driving for real many-level systems
(annealing bottlenecks, pairing sweeps) has no exactly solvable nontrivial benchmark.
**Type-M leverage.** For the commuting family the adiabatic connection *is* the a-independent
W (asset 4) — rational data on the curve; its norm, pole structure at the node, and behavior at
avoided crossings are exactly computable, member-by-member across the whole family. Breaking
integrability then deforms an *exactly known* AGP. Because Type-1 sectors are Richardson–Gaudin
(BCS pairing, central spin), the exact CD protocol is for a physical sweep, not a toy.
[AD for the family-direction AGP; export = conjecture]
**First probes.** Closed form for ‖W‖² from the curve (N=3 first, then N); compare its
divergence structure at the node/avoided crossings with the Pandey et al. scaling; write the
exact CD term for a Richardson sweep and measure the fidelity gain.
**Risk.** Must be careful *which* deformation direction (a-direction exact vs ε,γ-directions);
prior art on free-fermion AGP (Pozsgay et al., SciPost 17, 075 (2024)) shows the community is
moving here — timing favors speed.

### P3 — `Solvable ⟺ rigid` as the classification of driven multilevel protocols (finished physics)
Unchanged top pick from v1, now theorem-backed (wildness proof): a driven N-level protocol has
closed-form amplitudes iff its connection is rigid; the solvable zoo is the rigid stratum; the
criterion is *constructive* (how to engineer solvability: rank-1 coupling, slope coincidence,
degeneracy parked on a singular point) and *certifying* (proves the generic case has no closed
form). Outstanding problem solved: "why these solvable models and no others" — a 90-year zoo
unified. Cheapest deliverable: the standalone classifier note for the quantum-control /
exactly-solvable-dynamics community. [AD+NS; no σ-heroics needed]

### P4 — Statistics and lifting of protected degeneracies (diabolical points, avoided-crossing networks)
**Outstanding.** Distributions of diabolical points / conical intersections and their lifting
under perturbation matter in photochemistry (funnels), molecular magnets, adiabatic quantum
computing (bottleneck gaps) — and are treated phenomenologically; there is no solvable model of
diabolical-point statistics.
**Leverage.** OWY crossings are algebraically located nodes of an explicit discriminant
([1,1,1,1,2] real-root structure); breaking integrability lifts them in a computable way; the
δⱼ monodromy is the Berry phase. Type-M gives the first *ensemble* of guaranteed, located,
protected crossings. [conjecture; probe: node-lifting statistics under +εGOE, already scoped in
the 2026-06-03 log entry]

### P5 — Eigenvector statistics / non-ergodic phases in a solvable ensemble
**Outstanding.** The debates around Rosenzweig–Porter-type non-ergodic extended phases and
multifractality hinge on models where eigenvector statistics are only partially controllable.
**Leverage.** Type-M eigenvectors are explicit rational objects on the curve (adjugate columns,
R34) — an ensemble with *exact eigenvector* control and Poisson levels; its overlap/multifractal
spectra are computable, giving an analytic baseline the RP debates lack. [conjecture]

### P6 — Pairing-quench dynamics (Richardson–Gaudin, v1 II.1)
The node/discriminant geometry as sharper invariants for BCS quench phase diagrams (dynamical
gap vanishing, steady-state classification); node migration as a quench-phase invariant.
[grounded lineage; conjectural payoff]

---

## II. Deep / structural long shots (kept, upgraded)

- **σ ↔ Fredholm/Painlevé kernels** (v1 I.2): match σ's Widom-class determinant to an integrable
  kernel → σ computable + a new solvable-kernel family for RMT gap probabilities. [EST(σ=Fredholm)
  + conjecture]
- **Type-M genus hierarchy ↔ Argyres–Douglas / class-S** (v1 II.2, upgraded by R33): with
  "type = genus M−1", the *whole classification* maps to a hierarchy of SW-type curves, not one
  curve; σ = Stokes/wall-crossing datum; and `rigid = solvable` (control) = "the AD theory is
  Lagrangian" (gauge) — one statement, two languages. Needs the literature identification pass.
- **Bispectral ε↔a duality** (R39, new since v1): the Gaudin sphere (marked points at ε, weights
  γ²) vs the sweep sphere (punctures at slopes a) exchange smells of an MTV-type bispectral
  duality of the *driven* family — if exact on P, a genuinely new duality of driven quantum
  systems. Cheap decisive probe already flagged.
- **Wild cards** (v1): joint diagonalization / tensor decomposition (shared eigenbasis as
  algorithmic primitive); non-Hermitian/EP design (complexified branch points); rigidity ↔
  hypergeometric motives.

---

## III. Tournament (leverage × groundedness × outstandingness × cost)

| rank | direction | why it wins | first probe cost |
|---|---|---|---|
| 1 | **P1+P2: the integrability-breaking laboratory** (distance order parameter + M-dial + exact AGP) | attacks a *recognized open problem* (quantum KAM / ETH onset) with three handles nobody else has; home-turf lineage (SSY 2016 is the ancestor, not the competitor) | days (extend levelstats; AGP from curve) |
| 2 | **P3: rigid ⟺ solvable classifier note** | finished physics; unifies the solvable zoo; zero new computation | write-up only |
| 3 | **P4: diabolical-point statistics** | first solvable model of degeneracy statistics; feeds annealing + photochemistry | days |
| 4 | **P5: eigenvector-solvable ensemble** | analytic entry into the RP/multifractality debates | ~week |
| 5 | **II: AD/genus hierarchy + duality probe** | deepest unification; blocked on a lit pass | lit pass + 1 probe |

## IV. Recommendation

Develop the **integrability-breaking laboratory (P1+P2)** as the flagship: it is the one
direction where *every* Type-M asset (complete parameterization, exact curve, protected nodes,
exact connection) bears on a single, hot, genuinely open problem — and where the 2016 SSY level
statistics paper proves the family's statistics are nontrivial *within* the solvable class.
Publish **P3** alongside as the self-contained gift. Keep the ε↔a duality probe (II) as the
cheap curiosity with the highest surprise-per-cost.

*Discipline (META_REVIEW):* every export above is `[conjecture]` until its first probe runs;
M-dial statics must cite SSY 2016 as prior art; AGP claims must state which deformation
direction is exact before any norm is quoted.
