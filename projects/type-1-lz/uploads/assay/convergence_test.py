"""
Pairwise convergence and empirical agreement harness — Phase III shared
infrastructure.

Three APIs:

  * :func:`compare_routes` — pairwise residual :math:`\\|P_X - P_Y\\|_\\infty`
    between two ``formula(record) -> 3x3`` closed-form predictors. Returns
    per-record + per-stratum statistics.

  * :func:`compare_inter_window` — U(3) geodesic distance between two
    ``route_D(record) -> 3x3`` inter-window matrix predictors. Used when
    the predictive routes (A/B/C extensions) expose a thin ``hn_D``,
    ``cluster_D``, ``ddp_D`` accessor.

  * :func:`agree_with_empirical` — U(3) geodesic distance between a
    predicted D and the empirical D extracted from the IP propagator. This
    is the **true gate**: matching the empirical D is what residual ≤ 1e-8
    on P requires.

The U(3) geodesic distance is the gauge-invariant comparison metric
(see :func:`assay.d_extract.d_U3`).
"""

from __future__ import annotations

from typing import Callable, Dict, List, Optional

import numpy as np

from .d_extract import d_U3, empirical_D
from .validation_matrix import load_dataset, summarize


__all__ = [
    "compare_routes",
    "compare_inter_window",
    "agree_with_empirical",
]


def _per_stratum_summary(per_record: List[dict]) -> Dict[str, dict]:
    """Aggregate per-record residuals into per-stratum statistics."""
    by = {}
    for r in per_record:
        by.setdefault(r["stratum"], []).append(r["residual"])
    out = {}
    for st, vals in by.items():
        arr = np.array(vals, dtype=float)
        out[st] = {
            "n": int(len(arr)),
            "mean": float(arr.mean()),
            "median": float(np.median(arr)),
            "p90": float(np.percentile(arr, 90)),
            "p99": float(np.percentile(arr, 99)),
            "max": float(arr.max()),
        }
    return out


# ---------------------------------------------------------------------------
#  Compare two P-predictors (formulas)
# ---------------------------------------------------------------------------
def compare_routes(
    formula_X: Callable,
    formula_Y: Callable,
    dataset: Optional[list] = None,
    label_X: str = "X",
    label_Y: str = "Y",
) -> dict:
    """
    Pairwise comparison of two closed-form predictors.

    Parameters
    ----------
    formula_X, formula_Y : callable
        Each maps a ``Record`` to a 3x3 transition probability matrix.
    dataset : list, optional
        Validation dataset. If None, loaded from the default pickle via
        :func:`assay.validation_matrix.load_dataset`.
    label_X, label_Y : str
        Short identifiers for printout / table.

    Returns
    -------
    dict
        ``per_record`` (list of dicts with ``stratum``, ``residual``),
        ``per_stratum`` (dict of summary stats keyed by stratum tag),
        ``label_X``, ``label_Y``.
    """
    if dataset is None:
        dataset = load_dataset()
    per_record = []
    for rec in dataset:
        try:
            P_X = formula_X(rec)
            P_Y = formula_Y(rec)
            r = float(np.max(np.abs(P_X - P_Y)))
        except Exception as e:
            r = float("nan")
        per_record.append({"stratum": rec.stratum, "residual": r})
    return {
        "per_record": per_record,
        "per_stratum": _per_stratum_summary(per_record),
        "label_X": label_X,
        "label_Y": label_Y,
    }


# ---------------------------------------------------------------------------
#  Compare two D-predictors (inter-window matrices)
# ---------------------------------------------------------------------------
def compare_inter_window(
    route_D_X: Callable,
    route_D_Y: Callable,
    dataset: Optional[list] = None,
    label_X: str = "X",
    label_Y: str = "Y",
) -> dict:
    """
    Pairwise U(3) geodesic distance between two inter-window-matrix
    predictors.

    Parameters
    ----------
    route_D_X, route_D_Y : callable
        Each maps a ``Record`` (or ``Geometry``) to a 3x3 unitary
        ``D``. Concretely: ``assay.hitchin_holonomy.hn_D``,
        ``assay.voros_cluster.cluster_D``, ``assay.ddp.ddp_D``.
    dataset : list, optional
        Validation dataset.
    label_X, label_Y : str
        Short identifiers.

    Returns
    -------
    dict
        ``per_record`` (list with ``stratum``, ``residual`` = U(3)
        geodesic), ``per_stratum`` summary, ``label_X``, ``label_Y``.
    """
    if dataset is None:
        dataset = load_dataset()
    per_record = []
    for rec in dataset:
        try:
            D_X = route_D_X(rec)
            D_Y = route_D_Y(rec)
            r = d_U3(D_X, D_Y)
        except Exception:
            r = float("nan")
        per_record.append({"stratum": rec.stratum, "residual": r})
    return {
        "per_record": per_record,
        "per_stratum": _per_stratum_summary(per_record),
        "label_X": label_X,
        "label_Y": label_Y,
    }


# ---------------------------------------------------------------------------
#  Compare a D-predictor to the empirical D
# ---------------------------------------------------------------------------
def agree_with_empirical(
    route_D: Callable,
    dataset: Optional[list] = None,
    label: str = "route_D",
    convention: str = "phase_free",
    sign_A: int = +1,
    sign_B: int = +1,
) -> dict:
    """
    U(3) geodesic distance from a route's predicted D to the empirical D.

    Parameters
    ----------
    route_D : callable
        Maps a ``Record`` to a 3x3 unitary D.
    dataset : list, optional
        Validation dataset.
    label : str
        Short identifier for the route.
    convention : {'phase_free', 'stokes_phased'}
        Which convention to use for the empirical D extraction.
    sign_A, sign_B : int
        Stokes-phase signs (only used in stokes_phased convention).

    Returns
    -------
    dict
        ``per_record``, ``per_stratum``, ``label``, ``convention``.
    """
    if dataset is None:
        dataset = load_dataset()
    per_record = []
    for rec in dataset:
        try:
            D_pred = route_D(rec)
            emp = empirical_D(
                rec,
                convention=convention,
                sign_A=sign_A,
                sign_B=sign_B,
            )
            r = d_U3(D_pred, emp["D"])
        except Exception:
            r = float("nan")
        per_record.append({"stratum": rec.stratum, "residual": r})
    return {
        "per_record": per_record,
        "per_stratum": _per_stratum_summary(per_record),
        "label": label,
        "convention": convention,
    }


# ---------------------------------------------------------------------------
#  Pretty-print
# ---------------------------------------------------------------------------
def print_summary(result: dict) -> None:
    """Render a comparison dict's per-stratum table."""
    title_X = result.get("label_X", result.get("label", "X"))
    title_Y = result.get("label_Y", "")
    title = f"{title_X} vs {title_Y}" if title_Y else title_X
    print(f"\n--- {title} ---")
    print(f"  {'stratum':<8} {'n':>3} {'mean':>11} {'p50':>11} {'p90':>11} "
          f"{'p99':>11} {'max':>11}")
    print(f"  {'-'*70}")
    for st in sorted(result["per_stratum"].keys()):
        s = result["per_stratum"][st]
        print(f"  {st:<8} {s['n']:>3} {s['mean']:>11.3e} {s['median']:>11.3e} "
              f"{s['p90']:>11.3e} {s['p99']:>11.3e} {s['max']:>11.3e}")
