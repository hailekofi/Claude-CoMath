# Claude Co-Theoretical-Physicist

A suite of [Claude Code](https://code.claude.com) **skills** that turn Claude into an
effective **co-theoretical-physicist** — a research partner that reasons from physical
intuition, validates conjectures with carefully constructed numerical experiments, and
collaborates the way a strong human colleague would.

## Background

The design adapts two systems — read in full as primary sources (PDFs in this repo) — and
pivots them toward theoretical physics:

- **AI co-mathematician: Accelerating mathematicians with agentic AI** (Zheng et al.,
  Google DeepMind, 2026) — [arXiv:2605.06651](https://arxiv.org/abs/2605.06651). Source of
  the **seven design principles** (embrace work beyond proof; refine intent; native
  artifacts; asynchronous steering; progressive disclosure; track uncertainty; preserve
  failed explorations), the stateful workspace with a **project coordinator** delegating to
  **parallel workstreams**, **hard programmatic constraints** against hand-waving, and the
  living **working paper**.
- **Accelerating scientific discovery with Co-Scientist** (Gottweis et al., *Nature*,
  2026) — [s41586-026-10644-y](https://www.nature.com/articles/s41586-026-10644-y). Source
  of the coalition of specialized agents — **Generation, Reflection, Ranking, Proximity,
  Evolution, Meta-review** — orchestrated by a Supervisor, the **Elo tournament** of
  scientific debates, and meta-review feedback (learning without retraining).

**The pivot:** compared with a comathematician, the co-physicist *de-emphasizes formal
proof* and *heavily emphasizes (1) physical intuition and (2) numerical validation as
near-proof* — a conjecture supported by converged, error-controlled, independently
cross-checked numerics is accepted as working truth. Theorem-proving becomes controlled
**derivation**, and a dedicated **reflection** (peer-review) and **tournament**
(rank-and-evolve) stage make the generate→reflect→rank→evolve loop explicit.

The two source papers are included in this repo as PDFs for reference.

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

Nine composable skills. Each maps onto an agent role from the source papers:

| Skill | Fires when you… | Source role |
|-------|-----------------|-------------|
| **`cophysicist`** | start/frame a research task; coordinate; route | Project coordinator / Supervisor |
| **`physics-intuition`** | generate/sharpen conjectures, sanity-check, estimate, reason by symmetry/limits/scaling | Generation |
| **`physics-literature`** | check prior work, find references and standard techniques | grounding / lit-review |
| **`physics-reflection`** | stress-test a claim: verify assumptions, novelty, failure modes | Reflection |
| **`physics-tournament`** | rank competing hypotheses and evolve the leaders | Ranking + Proximity + Evolution |
| **`physics-numerics`** | test a conjecture numerically with convergence, error bars, cross-checks | simulation / verification |
| **`physics-derivation`** | derive a result analytically with controlled approximations (not formal proof) | proving (informal) |
| **`physics-theory-building`** | assemble results into a model + working paper; research overview | Meta-review synthesis |
| **`physics-research-log`** | track hypotheses, evidence, dead ends; meta-review feedback | Context memory + Meta-review |

A typical loop: **intuition → reflection → tournament → numerics/derivation → theory**,
logged throughout in `physics-research-log` and revisiting `physics-literature` as needed.

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
