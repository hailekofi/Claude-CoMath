---
name: physics-research-log
description: >-
  The stateful research workspace — maintains a persistent log of open questions,
  hypotheses with evidence-ladder status, numerical/analytic evidence, dead ends and WHY
  they failed, decisions, and next steps. Use at the start of any session to recover
  context, and whenever a hypothesis is proposed, tested, derived, refuted, or a direction
  is abandoned. This is what lets work survive across asynchronous sessions and prevents
  silently re-trying failed approaches.
---

# Physics Research Log

Theoretical-physics research is long, branching, and full of dead ends. Without a
persistent, honest log, context evaporates between sessions and the same failed approach
gets retried. This skill is the stateful workspace that mirrors a human collaborator's
notebook.

## When to use

- **Session start:** read `RESEARCH_LOG.md` (and any per-problem logs) before doing
  anything else, so you don't repeat refuted ideas or lose threads.
- **On every meaningful event:** a hypothesis proposed, a numerical test run, a derivation
  completed, a conjecture refuted, a direction abandoned, or a decision made.

## What to record

Maintain `RESEARCH_LOG.md` at the repo (or project) root, plus optional per-problem files.
For each hypothesis track its **evidence-ladder status**:
`conjecture → numerically-suggestive → numerically-supported (near-proof) →
analytically-derived → established`, or `refuted`.

Record explicitly:
- **Open questions** — what we are actually trying to settle, and the regime of interest.
- **Hypotheses** — statement, current status, and links to the hypothesis card / numerics
  checklist / derivation.
- **Evidence** — key numbers (with error bars), what convergence/cross-checks were done,
  plot/script paths.
- **Dead ends** — what was tried, *why it failed*, and whether it could work in a
  different regime. This is the most valuable and most often-skipped entry.
- **Decisions & next steps** — what we chose to pursue and what to do next session.

## Method

1. At session start, summarize the current state from the log in 3–5 lines before working.
2. As work proceeds, append timestamped entries; don't overwrite history.
3. Keep statuses current — promote/demote hypotheses as evidence changes.
4. When abandoning a direction, write the dead-end entry *before* moving on.

Use `templates/research-log.md` as the starting structure.
