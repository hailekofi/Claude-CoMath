# Attack lines for closed-form Type-1 N=3 P — commuting partner / exact crossing / Q4-pair structure

Generated 2026-06-01 (physics-intuition → reflection → tournament). Target: middle survival
P_{2→2} and one independent off-diagonal, in closed form over {γ,ε,a}, no numerical ODE, ≤ a
handful of geometric/Abelian integrals.

## Grounding fact (numerically-supported, exact to 1e-17)
Across the commuting family (fixed γ,ε; varying slopes a): **p, n, W₄, Q₄, κ_S are IDENTICAL**;
`a` enters ONLY through the linear factor **L_H**, which sits in the phase one-form
ω ∝ √Q₄·L_H·W₄/p³ and **cancels** in the transport one-form η ∝ p²/(Q₄√W₄). So:
- the transport curve ν²=W₄ and its periods are **family-invariant** (a-independent);
- ALL slope dependence of P flows through L_H, and **L_H is linear in a**.

This is the precise content of "exploit the commuting partners."

---

## H-A: Slope-flow / isomonodromy in a (commuting-family deformation)
**Statement.** Since the connection's geometry (Q₄,W₄,p) is fixed and a enters linearly via L_H,
the family {H(a)} is an isomonodromic-type deformation. P(a) should satisfy a closed first-order
flow ∂P/∂a_k = [P, M_k] (or an affine system) with M_k built from L_H-residues at ε_k — integrable
in closed form, anchored by the exactly-known BE extreme entries at the family's boundary slopes.
**Why practical.** A flow with known anchor + algebraic generator integrates to algebraic/exponential
closed form; no ODE-in-time.
**Cheapest test.** Numerically build P(a) on a 1-param slope ray a(s)=a₀+s·e_k via the benchmark;
check d/ds log-derivative structure matches a candidate algebraic M_k(s). (commuting partners)
**Risk.** "P(a) closed flow" may require the amplitude (gauge) not just |amplitude|²; flow may not
close on probabilities alone.

## H-B: Exact crossing = removable / computable local monodromy fixes the middle datum
**Statement.** Type-1 N=3 has a guaranteed exact (real) level crossing. At that degeneracy the local
monodromy is a computable rotation in the 2-dim crossing eigenspace; the middle-level survival is
fixed by composing this exact local datum with the two BE extreme periods via unitarity + double
stochasticity. P_{2→2} = (algebraic function of the crossing data) — elementary in {γ,ε,a}.
**Why practical.** Exact crossing ⇒ a marked regular point with explicit local data; no resummation.
**Cheapest test.** Locate the exact crossing u_c (real root structure of resultant of q_u); extract
the 2×2 local rotation; test P_{2→2}^pred vs benchmark on the 3 canonical samples. (exact crossing)
**Risk.** The crossing may be diabatic (no gap) → trivial local monodromy → no new datum; need to
confirm it is an *avoided/active* crossing that carries phase.

## H-C: Avoided crossings at Re(Q₄ conjugate-pair) points carry the missing off-diagonal
**Statement.** The off-diagonal (interference) parameter is set by the avoided crossings located at
the real parts of the Q₄ conjugate-pair points (the window "centers" τ in geometry.windows). The
missing off-diagonal = a finite sum of two-window contributions with a relative phase = Re of a
mixed period (the part the genus-1 phase curve cannot see, living where the two windows' real parts
coincide). 
**Why practical.** A finite sum over the 2 windows + one real phase integral — exactly "a handful
of Abelian integrals."
**Cheapest test.** Compute the two window centers u_center(τ); form candidate
P_off = |√p_{X1} ± √p_{X2} e^{iφ}|² with φ = Re(mixed period); fit φ to benchmark, then check φ is
geometric (matches a computed period, not a free fit). (avoided crossings at Re Q₄ pairs)
**Risk.** φ might be a genuine transcendental Stokes constant (frame-map test warned of this) → not
a period → fails practicality. The test (φ geometric vs free) is decisive.

## H-D: KZ/Gaudin monodromy specialization (from literature, top of prior memo)
**Statement.** P is the connection matrix of the rational KZ/Gaudin flat connection; BBGY 2409.17053
N=3 solution specializes to H₀+uA, giving P_{2→2} and off-diagonal as monodromy/period data.
**Why practical.** Output is connection data in {γ,ε,a}. **Cheapest test.** Specialize their N=3
formula; demand BE entries drop out; benchmark. (commuting partners, integrability)
**Risk.** Their model is Â+B̂/t (hyperbolic); specialization to linear uA may be singular.

## H-E: Constraint closure (count conserved quantities vs unknowns)
**Statement.** Commuting triple (H_x,Y,R) + quadratic commuting operators (Chernyak-Sinitsyn) impose
linear S-matrix constraints; with BE (2 params) + unitarity + double stochasticity, count whether
the remaining 2 are pinned. **Why practical.** If it closes, algebraic. **Cheapest test.** Write the
∞-limit constraints from [H,Y]=[H,R]=0; count independent eqns vs unknowns. (commuting partners)
**Risk.** May only close adiabatically (leading order), not exactly.

---

## Reflection (deep-verification highlights)
- H-A assumes a flow closes on **probabilities**; physically the natural closed object is the
  **amplitude/gauge** S-matrix. Demote unless reformulated on amplitudes. Fundamental-ish.
- H-B hinges on the crossing being **phase-active** (gapped/avoided), not a bare diabatic crossing.
  Must verify before investing — but if active, this is the cleanest route to P_{2→2}.
- H-C's make-or-break is **φ geometric vs transcendental** — a single decisive, cheap numerical
  test already expressible in the existing harness (compare fitted φ to a computed mixed period).
- H-D/H-E import external machinery; both gated on reading 2409.17053 / 2006.15144 directly.
- Falsified-ansatz guard: none of these is a product-of-two-level-factors (that class is dead).

## Tournament (Elo seed 1200; criteria novelty×correctness×testability×PRACTICALITY)
| Rank | Line | Why it leads | Decisive cheap test |
|------|------|--------------|---------------------|
| 1 | **H-B exact crossing → P_{2→2}** | directly targets the middle datum the phase curve can't see; exact local data, no resummation; meets bar | is the guaranteed crossing phase-active? extract 2×2 local rotation |
| 2 | **H-C Re(Q₄-pair) avoided crossings → off-diagonal** | targets the *other* missing param; finite period sum; one decisive test (φ geometric?) | fit φ, test if φ equals a computed mixed period |
| 3 | **H-A slope-flow/isomonodromy** | uses commuting structure rigorously; closed flow if it closes on the right object | log-derivative of P(a) along a slope ray vs algebraic M_k |
| 4 | **H-D KZ/Gaudin** | literature solves N=3; but specialization risk | specialize BBGY; BE must drop out |
| 5 | **H-E constraint closure** | cheap to attempt; may only close adiabatically | count eqns vs unknowns |

**Meta-review.** H-B and H-C are complementary — together they target *both* missing parameters,
both meet the practicality bar, and both have a single decisive cheap test runnable in the existing
harness. They are the recommended near-term pair. H-A is the principled backbone if a clean
amplitude flow can be written. H-D/H-E remain the literature-import hedges from the strategy memo.
