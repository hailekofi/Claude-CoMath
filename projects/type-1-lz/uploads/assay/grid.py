"""
Closed-form transition matrix: the integrable multistate Landau-Zener grid.

Type-1 N=3 is an integrable multistate-LZ model.  Its three diabatic
levels (energies H0_ii + a_i u) cross pairwise at three points

    u_ij = (H0_jj - H0_ii) / (a_i - a_j) .

Each crossing of pair (i,j) is an ordinary 2-level Landau-Zener event with
the elementary adiabaticity parameter

    Gamma_ij = gam_i^2 gam_j^2 |a_i - a_j| / (eps_i - eps_j)^2 ,
    q_ij     = exp(-2 pi Gamma_ij)        (diabatic stay probability).

Integrability => the transition *probabilities* carry no Stokes-phase
interference, so the diabatic scattering matrix is the ordered product of
real rotations

    S = R(last) ... R(first) ,    R_ij = 2x2 rotation, cos^2 = q_ij,

ordered by crossing position u_ij, and  P_diabatic = S .^ 2  (elementwise).

The adiabatic transition matrix of the gauged ODE is P_diabatic up to the
fixed endpoint channel permutation.
"""

from __future__ import annotations
import numpy as np
from itertools import permutations

from .geometry import Geometry


def lz_parameters(geo: Geometry) -> dict:
    """Elementary 2-level LZ data q_ij = exp(-2 pi Gamma_ij)."""
    g2, a, eps = geo.g2, geo.a, geo.eps
    q, Gam = {}, {}
    for i, j in [(0, 1), (0, 2), (1, 2)]:
        G = g2[i] * g2[j] * abs(a[i] - a[j]) / (eps[i] - eps[j]) ** 2
        Gam[(i, j)] = G
        q[(i, j)] = float(np.exp(-2.0 * np.pi * G))
    return {"q": q, "Gamma": Gam}


def crossing_order(geo: Geometry) -> list:
    """The three pairs ordered by crossing position u_ij along the real axis."""
    H0 = geo.H0
    a = geo.a
    items = []
    for i, j in [(0, 1), (0, 2), (1, 2)]:
        u_ij = (H0[j, j] - H0[i, i]) / (a[i] - a[j])
        items.append((u_ij, (i, j)))
    items.sort()
    return [pair for _, pair in items]


def _stochastic(pair: tuple, q: float) -> np.ndarray:
    """
    2(+)1 doubly-stochastic Landau-Zener crossing matrix: in the (i,j)
    block, stay with probability q and jump with probability 1-q; the
    spectator channel passes through with probability 1.
    """
    i, j = pair
    M = np.eye(3)
    M[i, i] = q
    M[j, j] = q
    M[i, j] = 1.0 - q
    M[j, i] = 1.0 - q
    return M


def grid_P(geo: Geometry) -> np.ndarray:
    """
    Integrable LZ-grid transition probabilities: the *incoherent* ordered
    product of the doubly-stochastic crossing matrices,

        P = M(last) ... M(first) ,

    ordered by crossing position.  Integrability removes Stokes-phase
    interference, so probabilities compose classically.
    P[x, j] is the x->j transition probability (rows = incoming).
    """
    data = lz_parameters(geo)
    order = crossing_order(geo)
    P = np.eye(3)
    for pair in order:                       # path order: first applied first
        P = _stochastic(pair, data["q"][pair]) @ P
    return P


def best_match(P_grid: np.ndarray, P_bench: np.ndarray):
    """Best row/col permutation aligning the grid P to a benchmark P."""
    best = (np.inf, None, None)
    for rp in permutations(range(3)):
        for cp in permutations(range(3)):
            e = np.max(np.abs(P_grid - P_bench[np.ix_(rp, cp)]))
            if e < best[0]:
                best = (e, rp, cp)
    return best
