"""
Level 0 control: exact segmentation identity (smoke test).

For boundaries  -T = t_0 < t_1 < ... < t_m = T  build segment factors
from the benchmark one-pass fundamental matrix,

    F[t_b, t_a] = U(t_b,-T) U(t_a,-T)^{-1} ,

and check  U(T,-T) = F[t_m,t_{m-1}] ... F[t_1,t_0].

This identity is mathematically tautological (nontrivial_assay_program,
section 5).  It validates indexing / multiplication order ONLY.  No claim
of independent continuation construction may rest on it.
"""

from __future__ import annotations
import numpy as np

from .geometry import Geometry
from .benchmark import fundamental_matrix


def segmentation_identity(geo: Geometry, T: float = 60.0,
                          boundaries=None, rtol: float = 1e-12,
                          atol: float = 1e-13) -> dict:
    if boundaries is None:
        boundaries = np.linspace(-T, T, 7)
    boundaries = np.asarray(boundaries, float)
    if boundaries[0] != -T or boundaries[-1] != T:
        raise ValueError("boundaries must start at -T and end at +T")

    # one-pass dense benchmark
    res = fundamental_matrix(geo, -T, T, rtol=rtol, atol=atol, dense=True)
    U_full = res.U
    sol = res.sol

    def U_at(t):
        return sol(t).reshape(3, 3)

    # segment factors from the global solution
    factors = []
    for b in range(1, boundaries.size):
        Ua = U_at(boundaries[b-1])
        Ub = U_at(boundaries[b])
        factors.append(Ub @ np.linalg.inv(Ua))

    prod = np.eye(3, dtype=complex)
    for F in factors:                       # path order: latest on the left
        prod = F @ prod

    err = float(np.max(np.abs(prod - U_full)))
    return {
        "n_segments": len(factors),
        "max_abs_recomposition_error": err,
        "passed": err < 1e-9,
    }
