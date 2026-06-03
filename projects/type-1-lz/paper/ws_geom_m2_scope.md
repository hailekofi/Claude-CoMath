# WS-GEOM Milestone 2 — scope: mapping the image region of $\Phi$ in reduced-$U(3)$ (facet A)

**Type:** `physics-numerics` scoping (a roadmap for review, NOT executed). **Date:** 2026-06-02.
**Builds on:** M1 (`paper/ws_geom_m1.md` — the node-selected directed-cycle skeleton), R15/R16 (the
two-transcendental skeleton), R9 (the classifier), R3 (the node). **Status:** SCOPED — pausing for review
before launch.

---

## 0. Goal (facet A of the geometric question)

M1 answered facet B (the holonomy *skeleton*: order 0 = the node-selected directed cycle, dressed by the
`W`-holonomy). M2 answers **facet A**: *which* part of `U(3)` does Type-1 occupy, and how? Concretely,
characterize the **image of the map** `Φ:(γ,ε,a)↦P` — its shape, its corners, its **boundary**, its
effective dimension, and its topological sectors — using large volumes of cheap data. The deliverable is a
quantitative map of the section, `σ`-free.

## 1. The right coordinates (from R15)

The transition matrix `P` is doubly stochastic (4 DOF). BE fixes the two extreme survivals
`P^{BE}_{lo}, P^{BE}_{hi}` (elementary functions of the actions `δ_{ij}`), and double-stochasticity then
writes everything in terms of just **two** numbers `{P_mm, b}` (R15). So the *interesting* image — the part
not already given by BE — is the **2-dimensional region swept by `{P_mm, b}`**, and the full image is that
region fibered over the elementary BE data. M2 works primarily in the `{P_mm, b}` plane (and its lift to
`{P_mm, b, δ_{ij}, χ, c_i}`).

## 2. The objects to characterize

1. **The two reachable vertices (topology, from M1+R16).** As the window action `δ` ranges, `P`
   interpolates between two permutation vertices of the Birkhoff polytope: the **identity** (diabatic
   corner, small `δ`, `{P_mm,b}→(1,0)`) and the **directed 3-cycle** (adiabatic corner, large `δ`,
   `{P_mm,b}→(0,1)`). M1 showed the cycle is *node-selected*; the transposition and the other three
   permutations are **not** reached. *Claim to test:* the Type-1 image touches exactly these two of the six
   vertices, fixed by the slope ordering.
2. **The directed-cycle-bias law.** Quantify the approach: distance of `{P_mm,b}` to the cycle vertex as a
   function of `δ_{max}` (and the cross-ratio `χ`). Expect monotone in `δ` (M1/R16), modulated by `χ`.
3. **The boundary = the decoupling locus (the M2 gate).** The image's boundary in the `{P_mm,b}` region
   should coincide with the **classifier/decoupling locus** (R9: a level decouples ⇔ a hop amplitude
   `W_{ij}→0` ⇔ the elementary boundary). *Decisive test:* drive one coupling `γ_k→0` (or one
   `s_{mid,·}→0`) and check `{P_mm,b}` lands on the boundary of the swept region, and that the boundary so
   traced matches the region's empirical edge.
4. **Effective dimension (sets up M3).** Is the map `(γ,ε,a)↦{P_mm,b,δ_{ij},χ}` full-rank, or does Type-1
   lie on a lower-dimensional subvariety (a hidden constraint)? Measure the local Jacobian rank / a PCA of
   the data cloud. *Codimension = number of invariants* → hands the target to M3. (Prior: `{P_mm,b}` are
   independent of each other, R15, so the `{P_mm,b}` plane is genuinely 2D; the question is whether the
   *lift* including the BE data has a relation.)
5. **The unistochasticity boundary.** Confirm `P` stays inside the unistochastic region (it must,
   `P=|S|^2`), and note where/whether the image touches the curved inner boundary of the Birkhoff polytope.

## 3. Method (cheap data + honest analysis)

- **Data.** Generate large stratified `(γ,ε,a)` samples (e.g. $10^4$–$10^5$), `γ,ε,a` over physical
  ranges, gauge-reduced. For each: `{P_mm,b}` (fast solver / WS-O3 in the generic regime, oracle in deep
  overlap), and `{δ_{ij}, χ, c_i, sep/width}`. Cache. (Cheap: seconds/sample for the oracle; the uniform
  formula is instant where valid.)
- **Region map.** 2D density of `{P_mm,b}`; convex hull / α-shape for the boundary; corner identification
  (distance to each Birkhoff vertex). Color by `δ_{max}`, `χ`, and the order-0 sector (identity vs cycle,
  from M1) to expose the topological structure.
- **Boundary test.** Constrained sampling along `γ_k→0` and along the genuine-LZ locus; overlay the
  decoupling locus; quantify boundary ↔ decoupling-locus agreement.
- **Dimension.** Local PCA on neighborhoods in parameter space mapped to `{P_mm,b,δ,χ}`; report the
  singular-value spectrum (a gap ⇒ a constraint ⇒ codimension for M3).

## 4. Milestones / gates
- **M2a:** the `{P_mm,b}` image map with the two-vertex (identity + directed cycle) structure confirmed,
  the other vertices shown unreachable. *Gate:* exactly two Birkhoff vertices in the closure of the image.
- **M2b:** the directed-cycle-bias law `dist(cycle)` vs `(δ_{max}, χ)` fitted/characterized.
- **M2c (the headline gate):** the image boundary **=** the decoupling locus (R9), quantified.
- **M2d:** the effective-dimension spectrum → the codimension / candidate constraint count for M3.

## 5. What success / partial / failure mean (honest)
- **Success:** a clean quantitative section map — two node-selected vertices, a `δ`–`χ` bias law, boundary
  = decoupling locus, and a definite effective dimension. This is a genuine geometric characterization of
  the `U(3)` section, `σ`-free, and it tells M3 exactly how many invariants to hunt for.
- **Partial:** the vertex/boundary structure is clean but the dimension is full (no constraint) ⇒ the
  section's analytic content is exhausted by `{BE, P_mm, b}` and only the **topological** label (the
  two-vertex/spinor sector) is new — still a real result.
- **Failure (unlikely):** the image is featureless / the boundary does not track the decoupling locus ⇒
  record and reassess (would contradict R9; a surprise worth chasing).

## 6. Risks
1. **The image may be "as structured as we already know."** If M2 only re-finds BE+DS+`χ`, that is a
   first-class negative (see §5 partial), not a failure — but set expectations: the *new* content is most
   likely the **topology** (two node-selected vertices, the spinor sector boundary), not a new analytic
   constraint.
2. **Sampling/gauge bias.** The image shape depends on the parameter measure; report it as a *reachable
   set* (boundary, vertices, dimension — gauge-robust features), not a measure-dependent density.
3. **Deep-overlap data cost.** `{P_mm,b}` there needs the oracle (seconds), not the uniform formula; budget
   the deep-overlap samples (they are the ones near the cycle vertex).

## 7. Hand-off to M3
M2 produces the cheap-data cloud `{P_mm,b,δ_{ij},χ,c_i}` and the effective-dimension spectrum. M3 (the
caged symbolic-regression invariant hunt) runs on exactly this cloud, targeting the codimension M2 finds
(a boundary function, a topological label, or a relation among the BE data and `χ`), with every candidate
handed to `physics-derivation` to prove.

## 8. Reusable assets
oracle.py / num_S12.py (cheap+gold data), ws_o3_uniform.py (instant generic-regime `{P_mm,b}`),
skeleton_two_transcendentals.py (the `{P_mm,b}` coordinates, the directed cycle, the incoherent baseline),
ws_geom_magnus.py (the order-0 sector label), ws_c_factorization_locus.md (the classifier/decoupling
locus). Companion: `paper/ws_geom_scope.md`, `paper/primer_magnus_feynman_geometry.{md,tex}`.
