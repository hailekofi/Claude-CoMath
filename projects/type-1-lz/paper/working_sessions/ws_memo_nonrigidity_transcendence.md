# Memo addendum: the non-rigidity → transcendence derivation, and the self-containment plan

**Status: review-gate addendum** (companion to `ws_memo_sectorial_factorization.md`).
Two user directives (2026-06-10):
(I) make the non-rigidity→transcendence link a **derivation in its own right**, not the
Heun analogy; (II) make the paper **self-contained** — Type-1 and every intermediate
quantity introduced in text or cited to a published source.

---

# PART I — The derivation: non-rigidity ⟹ σ, b non-Liouvillian

## I.0 What is to be derived (and in which sense)

**Claim.** For generic Type-1 N=3 data, the connection constants σ = P_mm and b = P[hi,lo]
are **non-Liouvillian functions** of (ε, γ, a): they lie in no field built from
ℚ(ε,γ,a) by finitely many algebraic extensions, exponentials, and integrals. Moreover
this hardness is a consequence of **non-rigidity**, and of non-rigidity *alone* — wildness
is neither sufficient nor the cause.

Three distinct "transcendence" notions, kept apart (this is where hand-waving usually
enters):

- **(A) No closed-form solution in u.** The ODE i ψ'=(H₀+uA)ψ has no Liouvillian
  fundamental system. This is a statement about the differential Galois group *in u*.
- **(B) Connection constant non-elementary in the parameters.** σ(ε,γ,a) is non-Liouvillian
  over Liouv(ε,γ,a). *This is the claim of interest* (and the trdeg-2 statement, R24).
- **(C) σ a transcendental number** at a fixed point. Hopeless and irrelevant (Lindemann
  makes even e^{−2πδ} transcendental).

The derivation targets **(B)**. The decisive point of §I.1 is that (A) ≠ (B).

## I.1 Why wildness cannot be the source (the N=2 separating witness)

N=2 Landau–Zener (Weber/Demkov) is **wild** (rank-1 irregular point), and its solutions —
parabolic cylinder functions — are **non-Liouvillian in u**: statement (A) holds. Its
differential Galois group in u is non-solvable (exponential torus + the two opposite
Stokes unipotents generate a Borel-exceeding subgroup of GL₂). **Yet its connection
constant is elementary:** P_{1→1} = e^{−2πδ}. So:

> (A) holds for N=2 but (B) fails for N=2. Therefore (A) ⇏ (B), and the
> Galois-group-in-u / wildness argument — however true — *cannot* deliver the
> transcendence of a connection constant. The two motions are logically independent, and
> N=2 is the witness that separates them.

Consequence: the transcendence of σ is **not** inherited from wildness. We must find what
N=2 has that Type-1 N=3 lacks. The answer is **rigidity**.

## I.2 Rigid side: rigidity *computes* the constant elementarily

**Lemma (rigid ⇒ closed-form constant).** If an irreducible connection is rigid (Katz
index χ = 2), it is constructed from rank-1 connections by a finite sequence of middle
convolutions and tensor twists (Katz [Katz1996]; algorithmic form Dettweiler–Reiter
[DR2000, DR2007]). Each operation has an explicit Euler/Laplace integral realization, so
the connection/Stokes data are **explicit Γ-ratios** in the local exponents (the rigid
connection formula; for the wild-rigid case, Jimbo's and Its–Lisovyy–Tykhyy-type explicit
Stokes data). N=2 LZ realizes this: it is rigid, and the algorithm returns precisely
e^{−2πδ}.

So an elementary connection constant is the **footprint of rigidity**. This already gives
the contrapositive the "empty cell" theorem uses — but the *positive* derivation of
hardness is §I.3.

## I.3 Non-rigid side: the constant becomes an isomonodromy transcendent

**Step 1 (deformation exists ⟺ non-rigidity).** Non-rigidity here is concrete: Katz index
χ = 0, one accessory parameter (½(2−χ)=1); equivalently the wild character variety has
dim 6 > 0 (Boalch [Boalch2014]). The connection is therefore a **non-isolated** point of a
positive-dimensional moduli M_wild, and σ, b are non-constant coordinate functions on it.

**Step 2 (the deformation is an isomonodromy flow).** Put the de-confluenced Fuchsian
system (Lemma "tame de-confluence") in normal form: 3 true singular points (fixable at
0,1,∞) plus the accessory datum realized as an **apparent singularity** whose position λ
is free. Isomonodromic motion of λ — the deformation preserving all monodromy/Stokes data
— is governed by the **Schlesinger / Jimbo–Miwa–Ueno equations** [JMU1981]; for the
rank-and-point count here this nonlinear flow is of **Garnier type** (the rank-3,
three-point analog of the Heun↔Painlevé-VI correspondence). In the wild picture the same
flows are Boalch's **wild isomonodromy** on M_wild, which realizes Painlevé/Garnier
systems geometrically [Boalch2009, Boalch2014]. Along this flow, the connection constant
σ(λ) **solves a second-order nonlinear (Painlevé/Garnier) ODE** — not a linear one.

**Step 3 (irreducibility ⇒ non-Liouvillian).** The general solution of a Painlevé/Garnier
isomonodromy equation is **not Liouvillian** over the field of coefficients: it lies in no
Picard–Vessiot extension. For Painlevé I–VI this is the Nishioka–Umemura irreducibility
theorem [Nishioka1988, Umemura1990] and its Galois-groupoid / Malgrange proofs (Casale
[Casale2008]); for Painlevé VI specifically, Cantat–Loray prove Malgrange-irreducibility
via the dynamics on the character variety [CantatLoray2009] — directly our setting, since
M_wild *is* the character variety and the flow *is* the dynamics on it. Hence σ(λ), as a
function of the deformation, is non-Liouvillian; pulling back along the (algebraic)
embedding of Type-1 data into M_wild, σ is non-Liouvillian in (ε,γ,a).

**Step 4 (genericity — the honest caveat).** Non-Liouvillian holds **off the
classical-solution locus** — the proper (codim ≥ 1) subvariety where the Painlevé/Garnier
flow admits algebraic or Riccati (one-parameter elementary) solutions. Type-1 generic data
is non-classical, witnessed two ways already in hand:
- **R24** (trdeg 2): σ and b are *functionally independent* across the family — incompatible
  with both lying on a one-parameter classical solution sheet;
- **Prop. (certification)**: the gold-gated falsification of every finite Barnes-G / Euler-Γ
  product rules out the rigid/classical closed forms at the canonical anchor.

A dedicated confirmation is proposed in §I.5.

## I.4 The disaggregation, now mechanized (the payoff)

The derivation *builds in* the two motions instead of asserting them:

| | provides | N=2 LZ | Type-1 N=3 |
|---|---|---|---|
| **Wildness** (irregular type A) | the Stokes slots; non-Liouvillian u-solutions (A) | ✔ wild | ✔ wild |
| **Rigidity status** | whether the slot *value* is pinned by local data | **rigid → elementary** | **non-rigid → Painlevé/Garnier transcendent** |

> **Location statement (boxed in the paper).** σ is hard *because it is a connection
> constant of a non-rigid isomonodromy (Painlevé/Garnier) family* — **not** because the
> connection is wild. Wildness positions the constant (the U·Δ·L slots, Δ elementary);
> non-rigidity, through isomonodromy and Painlevé irreducibility, makes the slot value
> non-Liouvillian. N=2 is wild and rigid → elementary; Heun is tame and non-rigid →
> transcendental; Type-1 N=3 is wild ∧ non-rigid. The transcendence rides on the second
> axis only.

This is strictly stronger than the Heun analogy: Heun (tame, non-rigid, transcendental
connection constant via PVI) is now the **base case** of the *same* mechanism, not a loose
parallel.

## I.5 Rigor ceiling and numerical confirmation (honest scope)

- **Rigorous:** §I.1 (N=2 separation), §I.2 (rigid⇒Γ-ratio), Step 1–2 of §I.3 (the
  deformation and its isomonodromy/Painlevé–Garnier nature).
- **Rigorous for the model reduction, frontier in full:** Step 3. The irreducibility
  theorem is complete for Painlevé I–VI (incl. the PVI character-variety-dynamics proof
  that matches our framing); for the genuine **rank-3 Garnier** flow the irreducibility is
  established in major cases but not with PVI-level completeness. This is the **one real
  gap** — a known-hard problem in differential Galois theory, flagged as such, not a
  hand-wave.
- **Numerically-supported:** §I.4 non-classicality (R24, certification), **plus the
  differential-algebraic order test below (R27)**.
- **Numerical confirmation — DONE (R27, `ws_painleve_da_test.py`).** A differential-algebraic
  *order* test: along a coupling-scale slice (γ→λγ; σ sweeps [0.016, 0.905]) compute σ(λ) to
  ~3×10⁻⁷–4×10⁻⁶, build a Chebyshev interpolant, and test for each (order k, degree d) whether
  an algebraic ODE P(t,σ,…,σ^(k))=0 exists (smallest singular value of the column-normalized
  monomial-jet matrix). Two controls calibrate the detector at the only diagnostic cell
  (order-1, d2): **exp** (Liouvillian) → 4×10⁻¹⁵ *hit*; **Bessel J0** (holonomic) → 3×10⁻³
  *miss* at order-1, 7×10⁻¹² *hit* at order-2.
  **σ result:** order-1 d2 = **1.5×10⁻⁴ (miss)** — ~37× its fit floor and the same order as J0's
  clean miss; order-2 d2 = 3×10⁻⁶ (at floor, ~hit). The order-1 residual is **stable under
  resolution** (1.50×10⁻⁴ at deg-20 → 1.48×10⁻⁴ at deg-26 while the floor improved 60×), so it
  is a real ODE-miss, not fit error. **Conclusion:** σ(parameter) tracks the *transcendent*
  pattern (order-1 NO / order-2 YES), not the *elementary* pattern (order-1 YES). This
  **confirms σ is not order-1 differentially algebraic = non-Liouvillian-classical** (robust,
  numerically-supported) and is **consistent with an order-2 (Painlevé/Garnier-type)
  transcendent** (suggestive: the order-2 residual is at-floor, not orders below).
- **Remaining numerical item.** Pinning the *exact* order-2 (Painlevé) structure — order-2
  residual orders below floor, as J0 achieves — needs σ to ~10⁻⁹ on a smooth slice
  (oracle/mpmath connection solver), beyond the ~10⁻⁶ interaction-picture floor. Likewise a
  *true* isomonodromic deformation (λ = apparent-singularity position, formal data fixed) would
  identify the specific Garnier equation; both flagged for the appendix.

---

# PART II — Self-containment plan

Target: a reader fluent in either MLZ physics **or** irregular-connection geometry can
follow end-to-end. Every intermediate quantity is defined in one line in text **or** cited
to a published source. Proposed §1 ("Model and objects", ~1.5 pp) introduces the physics;
a one-paragraph "Geometric dictionary" introduces the connection language; the glossary
below is the checklist (each row = an in-text definition + a citation).

| Quantity / notion | One-line in-text definition | Primary reference(s) |
|---|---|---|
| Multistate Landau–Zener (MLZ) | i ψ'(u) = (H₀ + uA)ψ, A=diag(a), swept crossing | [BE1993], [Sinitsyn-rev] |
| **Type-1 family / integrability** | the commuting family {H₀(a)+u·diag(a)} sharing an eigenbasis ∀u; the maximal abelian (Gaudin) algebra | [OWY2010] (classification), [Yuzbashyan-Gaudin] |
| Cauchy/Gaudin coupling H₀ | (H₀)_ij = γ_iγ_j(a_i−a_j)/(ε_i−ε_j), diag by row-sum rule | [OWY2010] |
| Diabatic transition matrix P | P_{x→j}=|S_xj|², doubly stochastic (unistochastic) | [BE1993] |
| Brundobler–Elser (BE) law + δ_ij | extreme-slope survivals = e^{−2π Σδ}, δ_ij=γ_i²γ_j²|a_i−a_j|/(ε_i−ε_j)² | [BE1993]; proofs [Dobrescu-Sinitsyn], [Volkov-Ostrovsky] |
| σ = P_mm, b = P[hi,lo] | middle survival; one recombination off-diagonal | this paper |
| Connection / scattering matrix S | transport between canonical asymptotic frames at u=∓∞ | [Wasow1965], [JMU1981] |
| Exponential torus, rates q_i | formal diagonal factor exp(q_i), q_i=½ a_i u² | [Wasow1965], [LR2016] |
| Formal monodromy | the u^L (diagonal, BE-weight) factor of the Thome solution | [Wasow1965], [LR2016] |
| Stokes matrices / multipliers | sectorial jumps gluing the formal solution; the U,L of U·Δ·L | [LR2016], [Boalch2014] |
| Wild character variety M_wild | moduli of the connection with fixed formal type; dim 6 here | [Boalch2014] |
| Poincaré rank / wild vs tame | pole order −1 of the connection 1-form; irregular vs regular | [Wasow1965] |
| Katz rigidity index χ | χ=(2−|S|)N²+Σ dim Z(g_s); rigid ⟺ χ=2 | [Katz1996] |
| Middle convolution | the rank-1→rigid construction (rigid ⇒ Γ-ratio constants) | [Katz1996], [DR2000] |
| De-confluence / apparent singularity / accessory parameter | bounded-ramp tame limit; Heun/Garnier normal form; the free modulus λ | [JMU1981]; Heun: [Ronveaux] |
| Isomonodromy / Schlesinger–JMU | deformation preserving monodromy+Stokes; λ-flow | [JMU1981], [Boalch2009] |
| Painlevé/Garnier transcendent | the nonlinear ODE the connection constant solves along the flow | [JMU1981], [Boalch2009] |
| Irreducibility (non-Liouvillian) | Painlevé solutions lie in no Picard–Vessiot extension | [Nishioka1988],[Umemura1990],[Casale2008],[CantatLoray2009] |
| Liouvillian / differential Galois | closed form = solvable Galois group (Kolchin) | [Kolchin], [vdPutSinger2003] |
| Ramis density (local Galois) | torus+formal monodromy+Stokes generate the local Galois group | [vdPutSinger2003], [Ramis] |

**New bibliography keys to add** (beyond the existing Wasow, JMU, Katz, Boalch, OWY, BE,
Sinitsyn): Dettweiler–Reiter (DR2000/DR2007); Loday-Richaud *Divergent Series, Summability
and Resurgence II* (LR2016); van der Put–Singer *Galois Theory of Linear Differential
Equations* (vdPutSinger2003); Nishioka (1988); Umemura (1990); Casale (2008);
Cantat–Loray (2009); Ronveaux (Heun); Dobrescu–Sinitsyn and Volkov–Ostrovsky (BE proofs);
Yuzbashyan (Gaudin/Type-1). All to be verified against the published record before they
enter the .tex (no key cited that I have not confirmed exists).

**Scope discipline.** Standard objects (exp torus, Stokes, Katz index, isomonodromy) get a
one-line definition + cite, not a tutorial; the *new* content (Type-1 specifics, the
U·Δ·L factorization, the slot equations, this derivation) is developed in full. The intro
states the physical model first (a sweep through a 3-level crossing, the survival σ) so a
physics reader is grounded before any connection language appears.
