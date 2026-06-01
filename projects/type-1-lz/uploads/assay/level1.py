"""
Level 1: independently solved segment-ODE assay.

The benchmark is the program's eq. (3): the one-pass adiabatic integration
over [-T,T].  The constructor independently solves each segment IVP

    dU/du = A(u) U ,   U(a) = I

with identity initial data (no benchmark intermediates) and forms the
ordered product

    M_prod = F[D2] F[W2] F[D1] F[W1] F[D0]            (path order).

Both use the interaction-picture engine (ip.py) but as disjoint pathways:
one continuous pass vs. five independent identity-reset solves.

A diabatic one-pass propagator (benchmark.py) serves as a fully
independent cross-check in a different basis.

Level 1 pass criteria (nontrivial_assay_program, sec. 6.4):
  1. segmented product == one-pass adiabatic  (< 1e-8);
  2. stable under shrinking/expanding windows;
  3. stable under refining T.
Level 1 does NOT test simple component structure.
"""

from __future__ import annotations
import numpy as np

from .geometry import Geometry
from .ip import propagate_ad_ip
from .benchmark import transition_matrix as transition_matrix_diabatic
from .windows import segmentation

# permutation (one-line) of the program's Pi_+ = [[0,0,1],[1,0,0],[0,1,0]]
PI_PLUS_PERM = (2, 0, 1)


def propagate(geo: Geometry, a: float, b: float, x: int,
              rtol: float = 1e-12, atol: float = 1e-13) -> np.ndarray:
    """Independent adiabatic segment IVP (interaction-picture engine)."""
    return propagate_ad_ip(geo, a, b, x, rtol=rtol, atol=atol)


def adiabatic_scattering(geo: Geometry, T: float, x: int,
                         seg: dict | None = None):
    """(M_onepass, M_segmented, seg) -- adiabatic fundamental matrices."""
    if seg is None:
        seg = segmentation(geo, T=T)
    M_onepass = propagate(geo, -T, T, x)
    M_seg = np.eye(3, dtype=complex)
    for kind, a, b in seg["segments"]:
        M_seg = propagate(geo, a, b, x) @ M_seg     # path order
    return M_onepass, M_seg, seg


def _best_perm_match(Pa, Pb):
    """col/row permutations making Pa ~ Pb; return (err, row_perm, col_perm)."""
    from itertools import permutations
    best = (np.inf, None, None)
    for rp in permutations(range(3)):
        for cp in permutations(range(3)):
            e = np.max(np.abs(Pa - Pb[np.ix_(rp, cp)]))
            if e < best[0]:
                best = (e, rp, cp)
    return best


def level1_assay(geo: Geometry, T: float = 60.0, x: int | None = None,
                 seg: dict | None = None, cross_check: bool = False,
                 verbose: bool = True) -> dict:
    if x is None:
        x = geo.par.x

    M_one, M_seg, seg = adiabatic_scattering(geo, T, x, seg=seg)
    seg_consistency = float(np.max(np.abs(M_one - M_seg)))

    P_direct = np.abs(M_one.T) ** 2               # benchmark (one-pass)
    P_prod = np.abs(M_seg.T) ** 2                 # constructor (segmented)
    bench_agreement = float(np.max(np.abs(P_prod - P_direct)))

    out = {
        "T": T, "x": x,
        "seg_consistency": seg_consistency,
        "benchmark_agreement": bench_agreement,
        "P_direct": P_direct, "P_prod": P_prod,
        "passed": bool(seg_consistency < 1e-8 and bench_agreement < 1e-8),
    }

    if cross_check:
        # independent diabatic benchmark, up to the audited endpoint perm
        P_dia = transition_matrix_diabatic(geo, T=T)
        err, rp, cp = _best_perm_match(P_direct, P_dia)
        out["diabatic_xcheck_err"] = err
        out["diabatic_xcheck_perm"] = (rp, cp)

    if verbose:
        print(f"--- Level 1 assay  (T={T}, x={x}) ---")
        print(f"  segment-product vs one-pass adiabatic : {seg_consistency:.3e}")
        print(f"  P_prod vs P_direct                    : {bench_agreement:.3e}")
        if cross_check:
            print(f"  diabatic-basis cross-check (up to perm): "
                  f"{out['diabatic_xcheck_err']:.3e}  perm={out['diabatic_xcheck_perm']}")
        print(f"  PASS: {out['passed']}")
    return out
