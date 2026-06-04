# uploads/ — provided sources, code, and finished documents

- **`writeups/`** — finished `.tex`/`.pdf` documents: the PRL main text + supplement,
  phase briefs, summaries, and the formal notes.
- **`assay/`** — importable Python package used by `../experiments/`
  (`from assay import Params, Geometry`, etc.). Experiments add it to `sys.path`
  via the relative path `../uploads`, so **it must stay at `uploads/assay/`**.
- **`type1_nontrivial_assay_and_graphify_artifacts/`** — graphify knowledge-graph
  artifacts (`.dot`, `.graphml`, `.json`, `.svg`, `.pdf` + README).
- **`lz_sweep_fast.py`, `lz_sweep_ip.py`** — standalone LZ sweep scripts.

LaTeX build intermediates (`*.aux`, `*.log`, `*.out`, `*.toc`) are gitignored and
regenerated on compile.
