# WS-GEOM Milestone 7 — multi-path / transfer-matrix sum (H3 reopened)

**Type:** `physics-literature` + `physics-numerics` + `physics-derivation`. **Date:** 2026-06-03.
**Builds on:** M5 (`ws_geom_m5.md` — the graph calculator, marginal `Λ≈π`), M6 (`ws_geom_m6_tournament.md`
— where H3 was deprioritized), R4/R5 (permanent overlap), R12 (`ws_o3_uniform.py` — the uniform law).
**Code:** `experiments/ws_geom_m7_multipath.py`. **Conventions:** NOMENCLATURE.md.

---

## 0. What M7 establishes (one paragraph)

The multi-path / independent-crossing transfer-matrix (ICTM) sum — Sinitsyn/Kayanuma-style coherent
interference of semiclassical trajectories through the crossing network — is a **genuinely distinct object**
from the Magnus tower (it is organized by *trajectory*, not by perturbative order), so the M5/Ω₃ negative
does not bear on it; M6 deprioritized it without test. Reopened properly (literature + clean numerics +
derivation), the verdict for Type-1 N=3 is a **first-class negative**: the ICTM does **not** improve on the
2-path uniform law (it is much worse), and the structural reason is R4 — the Cauchy weld **forbids crossing
separation**, so the independent-crossing approximation can never reach its exact regime and its error *is*
σ. This **demystifies** the method: multi-path is exact when crossings separate (or when a model's structure
sums the interference exactly, as in Sinitsyn's solvable 4-/6-state cases); Type-1 is the structural worst
case.

## 1. The distinction M6 missed

The M6 tournament ranked H3 last as "opaque matrix-product proliferation; marginal series" — conflating it
with the **Magnus tower** (organized by number of hops; asymptotic with `Λ≈π`; Ω₃ shown not to refine, M5).
The **multi-path sum** is different: a coherent sum over distinct *trajectories* through the crossing
network, each an *exact* per-crossing connection matrix times an accumulated phase. The Ω₃ negative says
nothing about it. (Correction logged: the low ranking was a too-quick dismissal.)

## 2. Literature grounding

Sinitsyn and collaborators obtain **exact** transition probabilities via path interference, but for
**special solvable models** (a four-state two-qubit model; a six-state model), where trajectories connecting
the same endpoints interfere coherently — explicitly *not* the incoherent LZ product. The
independent-crossing approximation factorizes into a product of 2-state jumps **only "if path interference
and accidental crossings are absent."** Type-1 N=3 general is not such a solvable model (its `P_mm` is the
σ-transcendental). Refs: Sinitsyn, *Exact transition probabilities in a six-state LZ system with path
interference*, arXiv:1501.06083; *Solvable four-state LZ model with path interference* (OSTI 1335604);
*Multipath interference in a multistate LZ model*; Kayanuma (independent-crossing approximation).

## 3. The ICTM construction (convention-safe)

Diabatic energies `d_i(u)=(H_0)_{ii}+a_i u`; constant couplings `V_{ij}=(H_0)_{ij}`; the 3 pairwise crossings
`u_{ij}=((H_0)_{jj}-(H_0)_{ii})/(a_i-a_j)`, ordered along `u`; partition `[-T,T]` at the midpoints between
consecutive crossings (one crossing per interval). Per interval, propagate the **exact isolated 2-level**
Hamiltonian `[[d_i,V_{ij}],[V_{ij},d_j]]` (DOP853) — this carries the correct Stokes phase automatically —
and give the spectator its diagonal diabatic phase; embed in 3×3 and multiply in `u`-order. The **only**
approximation is the factorization (each crossing treated as isolated 2-level + spectator); coherent
multi-path interference is fully included (it is a 3×3 matrix product). No hand-coded Stückelberg phases.

## 4. Numerics

**(B) ICTM vs the uniform law for `P_mm` (Type-1, canonical ε,a; coupling-scale sweep):**

| scale | ref `P_mm` | ICTM | uniform | \|ICTM−ref\| | \|uni−ref\| |
|---|---|---|---|---|---|
| 0.5 | 0.859 | 0.961 | 0.860 | 0.101 | **0.001** |
| 0.8 | 0.418 | 0.759 | 0.436 | 0.342 | **0.018** |
| 1.2 | 0.151 | 0.200 | 0.078 | **0.049** | 0.073 |
| 1.8 | 0.039 | 0.116 | 0.000 | 0.078 | **0.039** |
| 2.5 | 0.015 | 0.301 | 0.000 | 0.286 | **0.015** |

The uniform law wins 4/5; the ICTM is catastrophic at strong coupling. The full coherent path sum is **not**
a refinement for Type-1.

**(A) Validation, separable generic 3-level toy** (`b_i` free, crossings spreadable): the ICTM error
**decreases with crossing separation** (`0.16 → 0.043` as crossings `±5 → ±40`) but only slowly — consistent
with the ICA being asymptotic; we did *not* obtain a clean exact-in-the-limit validation, so the ICTM numbers
are indicative, not definitive.

## 5. The structural result (the real contribution)

**Type-1's Cauchy weld (R4/R5) prevents crossing separation.** Since `H_0 ∝ γ²`, turning the coupling *down*
sends all `(H_0)_{ii} → 0`, **collapsing every crossing onto `u=0`** (maximal overlap); turning it up merely
widens them. There is **no limit in which Type-1's crossings are isolated**, so the
independent-crossing/multi-path factorization can never reach its exact regime; its residual is the overlap
content = **σ**. This is the same `Λ≈π` permanent overlap that makes σ irreducible and the Magnus chart
marginal (M5/`where_we_are.tex`): the obstruction is one fact wearing three hats.

## 6. Verdict / evidence ladder

- `[numerically-supported]` ICTM ≫-worse than the uniform law for Type-1 `P_mm` (4/5 samples; §4).
- `[analytically-derived / structural]` R4 forbids crossing separation ⇒ multi-path is σ-limited for Type-1
  (§5); the separable-toy trend confirms separation (not the method) is what Type-1 denies.
- `[literature]` multi-path is exact only for special solvable models / separated crossings (§2).
- **H3 closes as a NEGATIVE for Type-1 N=3 — but a demystifying one:** it gives the precise condition
  (crossing separability) under which multi-path interference is the right tool, and shows Type-1 violates it
  maximally. Fully consistent with M5/M6/Ω₃: the uniform law is the best cheap estimate, σ the irreducible
  remainder. Folded into the working paper as R23(b).
