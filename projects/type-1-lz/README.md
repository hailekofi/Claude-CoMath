# Project: Type 1 Landau–Zener

Research workspace for the Type 1 Landau–Zener problem, run with the co-physicist skill
suite (see repo root `README.md`). The co-physicist acts as project coordinator; work
flows **intuition → reflection → tournament → numerics/derivation → theory**, logged in
`RESEARCH_LOG.md`.

## How to push your files here

From your local clone of `hailekofi/claude-comath`, on branch
`claude/comathematician-agent-skills-9HPSJ`:

```bash
cp -R "/Users/haileowusu/Documents/Claude/Type 1 LZ - Claude/." projects/type-1-lz/uploads/
git add projects/type-1-lz/uploads
git commit -m "Add Type 1 LZ source materials"
git push origin claude/comathematician-agent-skills-9HPSJ
```

Then tell me "pushed" and I'll pull and read everything in `uploads/`.

## Layout

| Path | Purpose |
|------|---------|
| `uploads/` | Your source materials (notes, drafts, code, references) — drop them here |
| `experiments/` | Reproducible numerical experiments (`physics-numerics`) |
| `paper/` | The living working paper / writeup (`physics-theory-building`) |
| `RESEARCH_LOG.md` | Persistent state: questions, hypotheses, evidence, dead ends |

## Status

Awaiting (1) your source files in `uploads/`, and (2) answers to the scoping questions so
intent is pinned down before any compute is spent.
