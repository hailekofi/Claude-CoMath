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

---

## V. The isospectral-torus asset (Hermitian Type-M, M ≥ 2) — added 2026-07-10

**The asset (user, 2026-07-10).** For ansatz Type-M with M ≥ 2 there are **M−1 additional
parameters** deforming a real-symmetric Type-M family into a *Hermitian* Type-M family with
**identical spectrum**. I.e. over every real-symmetric Type-M point sits an explicit
(M−1)-dimensional **isospectral torus** of Hermitian Hamiltonians: eigenvalues (as functions of
u) frozen; eigenvectors, and everything built from them, free to move. Type-1 (M=1) has no such
torus — consistent with everything the N=3 program found.

### V.0 The structural identification (the headline conjecture)

**[conjecture, dimension-matched]** The M−1 isospectral phases are the **angle variables on the
Jacobian of the genus-(M−1) spectral curve** — the fiber of the spectrum map. The two
conjectures now support each other:
- "Type-M ⇔ genus M−1" (R33-adjacent) predicts a Jacobian of complex dimension M−1; its real
  slice (the fixed torus of the anti-holomorphic involution picking out real-symmetric members)
  has real dimension M−1 — **exactly the user's parameter count**. Type-1: g=0, trivial
  Jacobian, no phases. ✓
- This is precisely the **finite-gap paradigm** (Dubrovin–Novikov: the isospectral manifold of a
  g-gap potential is the Jacobian of the genus-g curve; angles = norming-constant phases), and
  the Beauville–Mumford/Hitchin picture: the Type-M variety is an **algebraic completely
  integrable system** — spectral data = action/base coordinates, the θ-torus = the fiber. The
  θ-flows are Lax (isospectral) flows: dH/dθ = [H, A_θ].
So Type-M is the matrix pencil analog of finite-gap theory, with the type index = number of
"gaps" + 1. First test: on the minimal example, track the eigenvector (Mumford) divisor vs θ and
check *linear* motion on Jac. **CRITICAL PRE-CHECK (M1 discipline, gates everything below):**
verify the θ-deformation is **not pure gauge** — not a u-independent unitary (e.g. diagonal
phase) conjugation of the whole pencil. If it were, every claim below trivializes. Expected
genuine (the count M−1 is already modulo the N−1 diagonal-phase gauge), but must be computed
once, not asserted.

### V.1 What the torus is, physically: a spectrum-frozen dial

Everything **spectral** is exactly frozen along θ: char-poly coefficients, gaps, DOS, level
statistics, discriminant/node locations, Stückelberg actions, BE exponents, the formal monodromy
c_i. Everything **eigenvector-borne** may move: transition amplitudes' Stokes content,
AGP/quantum-geometric tensor, Berry curvature, eigenstate expectation values and entanglement,
currents. That makes θ a *scalpel that separates spectral from non-spectral physics* — with
theorem-grade cleanliness, because the freezing is exact, not approximate. New problems this
opens:

**Q1 — Gap-vs-dynamics separation theorems (driven systems / adiabatic computing).**
Outstanding: adiabatic theorems and annealing runtimes are gap-based; how much dynamical
information do gaps actually carry? Probe: does the LZ transition matrix P move along θ at
frozen spectrum? *Either answer is a result:* θ-invariant ⇒ P is a function of spectral data
alone (a massive new solvability lever: compute P from the curve, no dynamics) and the R31/R32
arity theorem generalizes verbatim; θ-varying ⇒ an exactly parameterized family of Hamiltonians
with identical gap structure and computably different transition probabilities — the spread over
the torus IS the irreducible non-spectral (eigenvector/Stokes) content, quantifying the slack in
all gap-based bounds. Note the N=3 Type-1 arity theorem P = F(be×3, χ) used *spectral*
invariants only — consistent with M=1 having no torus; Q1 is its M ≥ 2 continuation.
[conjecture; probe = minimal Type-2 example, days]

**Q2 — Time-reversal breaking that is invisible to the spectrum (RMT foil; T-violation bounds).**
The θ's deform real-symmetric (TRS, β=1) into genuinely Hermitian (β=2) **with zero spectral
signature** — impossible for generic ensembles, where GOE→GUE is *defined* spectrally.
Integrability screens TRS breaking from the spectrum. Consequences: (a) a precise counterexample
class for the statistical-spectroscopy program that bounds T-violation (nuclear data) from level
statistics — exposes the genericity assumption those bounds need; (b) under integrability
breaking +εV, does the endpoint universality class (GOE vs GUE) depend on θ? If yes, θ *steers
the universality class while spectrally invisible at ε=0* — the flagship laboratory (§P1/P2)
gains a second, orthogonal dial: (M, θ) = (how integrable, how time-reversal-broken), with the
AGP (exponentially TRS-sensitive per the PRX 2020 line) as the natural detector; (c) conceptual
hook: single-system spectra cannot distinguish real from complex quantum mechanics here — an
exactly solvable illustration of why the real-vs-complex-QM falsification program needed
entangled tests. [conjecture/framing]

**Q3 — Pure-holonomy control: geometric phase at exactly frozen dynamical phase.**
Along θ the instantaneous spectrum — hence every dynamical phase ∫E dt — is *identically*
constant; only the connection/holonomy (Berry, Wilczek–Zee) moves. Outstanding: separating
geometric from dynamical phase is the perennial obstacle in holonomic gates and geometric-phase
metrology. Type-M with M ≥ 2 is a completely parameterized control manifold on which that
separation is exact by construction. Bonus geometry: the protected node persists at the same
(u*, E*) for **all** θ (spectrum frozen), so the degeneracy sweeps out a manifold of
**anomalous codimension** in (u, θ)-space — Hermitian degeneracies generically need codimension
3, integrability + isospectrality beat that by construction — an exactly solvable factory of
Weyl-point/monopole/diabolical structures with algebraically known loci, sharpening §P4 and the
OWY crossing theorem simultaneously. Probe: Berry curvature and Chern numbers on (u, θ) around
the persistent node in the minimal Type-2 example. [conjecture]

**Q4 — Eigenvector chaos/ETH at frozen (Poisson) spectrum.**
ETH is an eigenvector statement; level statistics is spectral. The field routinely infers
eigenstate properties from spectral probes. The θ-torus moves eigenstate structure (matrix
elements of observables, eigenstate entanglement) at exactly fixed Poisson spectrum — a
counterexample factory for "spectrum ⇒ eigenstates" inferences and a clean testbed for which
diagnostics probe which sector. Combines with §P5 (rational eigenvectors ⇒ analytic overlap
statistics as functions on the torus). [conjecture]

**Q5 — The wild-character-variety orbit (math-phys anchor; upgrades the AD long shot).**
Frozen spectrum ⇒ frozen formal data (exponents, irregular types); the θ-action on the driven
problem can move **only the Stokes matrices**: the isospectral torus maps into the wild
character variety at fixed formal data — plausibly a (Lagrangian) torus orbit, i.e. the
angle-fibration of the wild Hitchin system made explicit at matrix-pencil level. This is the
concrete bridge the §II AD/class-S item was missing: base = curve moduli (frozen), fiber = θ.
[conjecture; needs the gauge pre-check first]

### V.2 Minimal probe (fully specified, cheap)

N=3, generic real-symmetric pencil = the smooth-cubic genus-1 samples already in hand (R33a
generic samples) ⇒ Type-2 ⇒ **one** phase θ; driver A=diag(a) as always. Steps: (1) construct
the Hermitian isospectral extension (OY 2011 Hermitian Type-2 ansatz, or solve isospectrality
directly); (2) **gauge check** (V.0); (3) verify spectrum-freezing to machine precision along θ;
(4) compute P(θ), ‖AGP‖(θ), Berry curvature at the (persistent) degeneracy vs θ. Outcomes gate
Q1–Q5. All machinery (propagators, oracle-grade solvers, curve builders) exists in
`experiments/`.

### V.3 Effect on the tournament

The isospectral torus does not displace the flagship — it **arms it**: §P1+P2 becomes a
two-dial laboratory (M = integrability dial, θ = spectrally-invisible TRS dial, both exactly
solvable), and Q1's either-way payoff makes the minimal Type-2 probe arguably the highest
information-per-cost experiment now on the board — cheaper than the M-dial sweep and decisive
for the structure (gauge or genuine; spectral or not) before anything else is built on M ≥ 2.
