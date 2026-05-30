---
name: physics-tournament
description: >-
  Prioritize and improve competing physics hypotheses through a structured tournament. Use
  when you have several candidate conjectures, models, or approaches and need to decide
  which to pursue, when comparing rival explanations head-to-head, deduplicating or
  clustering overlapping ideas, or systematically evolving the strongest hypotheses into
  better ones. Runs Elo-style pairwise scientific debates ranked on novelty, correctness,
  and testability, then evolves the leaders. Adapts the Co-Scientist Ranking, Proximity,
  and Evolution agents.
---

# Physics Tournament (rank & evolve)

When the workflow has produced many candidate hypotheses, you need a principled way to
spend the physicist's (and the compute's) limited attention on the most promising ones,
and to keep improving them. This skill combines the Co-Scientist **Ranking**, **Proximity**,
and **Evolution** agents into one rank-and-evolve loop.

## Ranking: an Elo tournament of scientific debates

Maintain an **Elo rating** for each hypothesis (start new entries at **1200**). Rank by
running **pairwise comparisons**, each settled by a short *simulated scientific debate*
judged on three criteria:

- **Novelty** — does it go beyond what is already known (cross-check `physics-literature`)?
- **Correctness** — does it survive `physics-reflection` scrutiny (assumptions, limits)?
- **Testability** — is there a concrete numerical (`physics-numerics`) or analytic
  (`physics-derivation`) experiment that could confirm or kill it?

Use **multi-turn** debates for the top contenders (where the decision matters most) and
quick **single-turn** comparisons for lower-ranked ones. Conclude each match with an
explicit "A beats B because…" and update Elo. Prioritize matches between **similar**
hypotheses and between **new or top-ranked** ones.

## Proximity: cluster and diversify

Compute a rough **proximity graph** — which hypotheses are near-duplicates or share a core
mechanism. Use it to **deduplicate**, to organize tournament matches (compare like with
like), and to ensure the surviving set stays **diverse** rather than ten variants of one
idea. Present clusters so the physicist can scan the landscape quickly.

## Evolution: improve the leaders (as new competitors)

Refine top-ranked hypotheses using these strategies:

- **Enhancement through grounding** — find the weakest link, search literature, fill the
  reasoning gap.
- **Coherence / feasibility** — fix invalid initial assumptions; make the hypothesis
  practical and testable.
- **Inspiration & combination** — generate new hypotheses inspired by, or merging the best
  parts of, the current leaders.
- **Simplification** — strip a hypothesis to its most easily verifiable core.
- **Out-of-the-box** — deliberately diverge from the current top set to escape local optima.

**Critical rule:** evolution *creates new hypotheses*; it never overwrites or replaces an
existing one. Each evolved hypothesis must re-enter the tournament and earn its rating.
This protects good ideas from bad "improvements."

## Output and hand-off

Produce a **ranked leaderboard** (Elo + one-line rationale per hypothesis), the cluster
map, and the newly evolved candidates. Send fresh hypotheses to `physics-reflection` for
review, decisive tests to `physics-numerics` / `physics-derivation`, and record the
standings and debate outcomes in `physics-research-log`. See `templates/tournament.md`.
