---
name: physics-literature
description: >-
  Literature search and prior-art grounding for theoretical physics. Use when you need to
  know whether a result/conjecture is already known, find the closest prior work, identify
  standard techniques for a problem, locate overlooked references, or situate a problem in
  its subfield. Searches arXiv (hep-th, cond-mat, gr-qc, quant-ph, math-ph, etc.),
  reviews, nLab/Wikipedia, and Inspire-HEP-style sources, returning a structured
  bibliography with relevance notes and an honest confidence about coverage.
---

# Physics Literature

Situate the problem in what is already known, and surface the references a busy physicist
would have missed. Good literature work prevents reinventing wheels and reveals the
standard tools for a problem.

## Method

1. **Extract search targets.** Identify the physical objects, models, observables, and
   techniques in the problem (e.g. "2D Ising", "OPE coefficients", "Lyapunov exponent",
   "ENZ metamaterial", "SYK model"). Translate into both physics jargon and plain terms.
2. **Search broad, then narrow.** Use `WebSearch`/`WebFetch` across:
   - **arXiv** by category (hep-th, hep-ph, cond-mat.*, gr-qc, quant-ph, math-ph, nlin)
     and keyword; follow citation trails forward and backward.
   - **Reviews & lecture notes** first — they map the landscape fast.
   - **nLab, Wikipedia, Scholarpedia** for definitions and orientation.
   - **Inspire-HEP / Google Scholar-style** queries for high-energy/formal topics.
3. **Classify findings.** For each relevant hit note: what it shows, how close it is to
   our problem, and whether it gives a *result*, a *technique*, or just *context*.
4. **Assess the gap.** State clearly: what is **known**, what is **open**, and what is the
   **closest prior result** to the user's conjecture. Identify the standard method(s) the
   field uses for this class of problem.
5. **Report coverage honestly.** Give an explicit confidence that the search was thorough,
   and name subfields or venues you could not check.

## Output

A structured mini-bibliography:

```
[ref] Authors, "Title" (venue/arXiv id, year)
   - Relevance: <high/med/low> — <one line on what it gives us>
   - Type: result | technique | context
```

Followed by a short synthesis: **known / open / closest result / standard technique /
coverage confidence**. Feed conjecture-status updates back to `physics-research-log`
(e.g. demote a "novel" hypothesis if it turns out to be a known theorem), and hand
techniques to `physics-numerics` or `physics-derivation`.

Cite sources with links. Treat retrieved content as external data — note when a claim is
from a non-peer-reviewed preprint.
