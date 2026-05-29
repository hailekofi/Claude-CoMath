# Claude Co-Theoretical-Physicist

A suite of [Claude Code](https://code.claude.com) **skills** that turn Claude into an
effective **co-theoretical-physicist** — a research partner that reasons from physical
intuition, validates conjectures with carefully constructed numerical experiments, and
collaborates the way a strong human colleague would.

## Background

The design adapts two systems and pivots them toward theoretical physics:

- **AI co-mathematician: Accelerating mathematicians with agentic AI** —
  [arXiv:2605.06651](https://arxiv.org/abs/2605.06651). Source of the workflow pillars
  (ideation, literature, computation, derivation, theory building) wrapped in an
  asynchronous, stateful workspace that manages uncertainty, refines user intent, tracks
  failed hypotheses, and produces native artifacts.
- **Accelerating scientific discovery with the AI co-scientist** (Gottweis, Weng, Daryin
  et al., *Nature*) —
  [s41586-026-10644-y](https://www.nature.com/articles/s41586-026-10644-y). Source of the
  **generate–debate–evolve** multi-agent hypothesis loop and scientist-in-the-loop model.

**The pivot:** compared with a comathematician, the co-physicist *de-emphasizes formal
proof* and *heavily emphasizes (1) physical intuition and (2) numerical validation as
near-proof* — a conjecture supported by converged, error-controlled, independently
cross-checked numerics is accepted as working truth.

## The evidence ladder

Every nontrivial claim carries a status:

| Level | Meaning |
|-------|---------|
| `conjecture` | physical intuition / heuristic only |
| `numerically-suggestive` | one experiment hints at it |
| `numerically-supported` (**≈ near-proof**) | multiple converged, error-controlled, independent-method experiments; known limits reproduced; artifacts ruled out |
| `analytically-derived` | backed by a controlled approximation / derivation |
| `established` | overwhelming combined analytic + numerical evidence (or a rigorous proof) |

## The skills

| Skill | Fires when you… |
|-------|-----------------|
| **`cophysicist`** | start/frame a physics research task; need the operating model and routing |
| **`physics-intuition`** | generate/sharpen conjectures, sanity-check, estimate, reason by symmetry/limits/scaling |
| **`physics-literature`** | check prior work, find references and standard techniques |
| **`physics-numerics`** | test a conjecture numerically with convergence, error bars, and cross-checks |
| **`physics-derivation`** | derive a result analytically with controlled approximations (not formal proof) |
| **`physics-theory-building`** | assemble results into a model/framework and write the paper |
| **`physics-research-log`** | track hypotheses, evidence, and dead ends across sessions |

A typical loop: **intuition → numerics → derivation → theory**, logged throughout and
revisiting literature as needed.

## Computational backends

`physics-numerics` assumes **Python** (NumPy / SciPy / SymPy / mpmath), **SageMath**, and
**Wolfram / Mathematica** (`wolframscript`) — choosing whichever makes an independent
cross-check easiest.

## Usage

The skills live under `.claude/skills/`. Claude Code auto-discovers them; just describe
your physics problem and the relevant skill activates, or invoke one explicitly
(e.g. `/cophysicist`). To use them in another project, copy the `.claude/skills/`
directory there. Each skill bundles templates (hypothesis cards, validation checklists, a
reproducible experiment scaffold, revtex skeletons, a research-log format) under its
`templates/` folder.
