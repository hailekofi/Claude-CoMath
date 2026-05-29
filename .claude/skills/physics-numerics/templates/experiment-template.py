#!/usr/bin/env python3
"""Numerical experiment scaffold for the co-theoretical-physicist suite.

Fill in the header, implement `run()`, and keep every knob explicit so the
experiment is reproducible and the validation checklist can be applied.

Hypothesis under test : <state the conjecture and its quantitative prediction>
Regime / parameters   : <couplings, dimensions, scales, limits>
Method                : <algorithm; note the independent cross-check method>
Expected (if true)    : <what we should see>
Backend / versions    : <python/numpy/scipy/sympy/mpmath/sage/wolfram + versions>
"""

from __future__ import annotations

import argparse
import json
import platform
from dataclasses import asdict, dataclass

import numpy as np

SEED = 12345  # fix the seed for reproducibility


@dataclass
class Params:
    """All experiment knobs in one place (vary these for convergence studies)."""
    resolution: int = 128      # grid size / basis size / sample count
    cutoff: float = 1.0e-12    # truncation / tolerance
    n_samples: int = 10_000    # Monte Carlo samples, if applicable
    # add problem-specific parameters here


def run(p: Params, rng: np.random.Generator) -> dict:
    """Run one instance of the experiment and return results + uncertainties.

    Return a dict with at least: {"value": float, "error": float, ...}.
    """
    # TODO: implement the experiment.
    raise NotImplementedError


def convergence_study(p: Params, rng: np.random.Generator) -> list[dict]:
    """Re-run while doubling resolution to demonstrate convergence."""
    results = []
    for k in range(4):
        pk = Params(**{**asdict(p), "resolution": p.resolution * (2 ** k)})
        res = run(pk, rng)
        res["resolution"] = pk.resolution
        results.append(res)
    return results


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default="results.json", help="where to save results")
    ap.add_argument("--converge", action="store_true",
                    help="run a convergence study instead of a single point")
    args = ap.parse_args()

    rng = np.random.default_rng(SEED)
    p = Params()

    payload = {
        "params": asdict(p),
        "seed": SEED,
        "env": {"python": platform.python_version(), "numpy": np.__version__},
    }
    payload["results"] = (
        convergence_study(p, rng) if args.converge else [run(p, rng)]
    )

    with open(args.out, "w") as f:
        json.dump(payload, f, indent=2)
    print(f"saved {args.out}")
    print(json.dumps(payload["results"], indent=2))


if __name__ == "__main__":
    main()
