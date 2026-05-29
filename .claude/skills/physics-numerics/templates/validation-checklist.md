# Numerical Validation Checklist: <hypothesis id / quantity>

**Hypothesis under test:** <statement + the quantitative prediction>
**Regime tested:** <parameters / scales>

| Check | Status | Evidence / notes |
|-------|--------|------------------|
| Convergence study (vary grid/step/cutoff/samples/precision) | ☐ | <how it converges; extrapolated value> |
| Error quantification (statistical + systematic) | ☐ | <error bars, bootstrap, finite-size estimate> |
| Independent methods (≥2 different algorithms/tools) | ☐ | <method A vs method B agreement> |
| Known-limit reproduction (solvable/published cases) | ☐ | <which limits, agreement> |
| Parameter sweep across regime of interest | ☐ | <range covered; where it breaks> |
| Artifact elimination (round-off, finite-size, RNG, NaN) | ☐ | <precision raised, size scaled, seed fixed> |

**Result:** value = <x> ± <σ>
**Verdict (evidence-ladder level):**
  numerically-suggestive | numerically-supported (near-proof) | refuted
**Regime of validity established:** <where the claim holds / fails>
**Reproducibility:** script = <path>, tool versions = <...>, seed = <...>
