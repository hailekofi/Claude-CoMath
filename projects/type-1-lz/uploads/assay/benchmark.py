"""
Benchmark one-pass Landau-Zener propagator (diabatic form).

The Type-1 N=3 evolution in the diabatic basis is

    i dpsi/du = (H0 + u A) psi ,        H = H0 + u diag(a)

(LZ_summary eqs. for H_ij, H_ii).  The fundamental matrix solves

    dU/du = -i (H0 + u A) U ,    U(-T) = I .

Transition matrix:  P[x->j] = |U(T,-T)[j,x]|^2 .

The diabatic transition matrix equals the adiabatic one of the program up
to the fixed endpoint permutations (Pi_-, Pi_+); see permutation.py.
This module is the *benchmark* pathway only -- constructors must never read
its intermediate values.
"""

from __future__ import annotations
import numpy as np
from scipy.integrate import solve_ivp

from .geometry import Geometry


def fundamental_matrix(geo: Geometry, t0: float, t1: float,
                       U0: np.ndarray | None = None,
                       rtol: float = 1e-12, atol: float = 1e-13,
                       dense: bool = False):
    """
    Integrate dU/du = -i (H0 + u A) U from t0 to t1.

    Returns the solve_ivp result-like object holding U(t1) as `.U` (3x3),
    and, if dense=True, a dense-output callable `.sol`.
    """
    if U0 is None:
        U0 = np.eye(3, dtype=complex)
    H0 = geo.H0.astype(complex)
    A = geo.A.astype(complex)

    def rhs(t, y):
        U = y.reshape(3, 3)
        dU = -1j * (H0 + t * A) @ U
        return dU.ravel()

    sol = solve_ivp(rhs, (t0, t1), U0.ravel(), method="DOP853",
                    rtol=rtol, atol=atol, dense_output=dense)
    if not sol.success:
        raise RuntimeError(f"benchmark integration failed: {sol.message}")
    out = _Result()
    out.U = sol.y[:, -1].reshape(3, 3)
    out.t0, out.t1 = t0, t1
    if dense:
        out.sol = sol.sol
    return out


def adiabatic_frame(geo: Geometry, u: float) -> np.ndarray:
    """
    Matrix V whose column i is the instantaneous eigenvector of H(u)
    belonging to diabatic channel i (matched by max overlap with e_i).

    At large |u| the eigenvectors -> diabatic basis, so the matching is
    unambiguous; projecting onto this frame removes the O(1/|u|) basis
    mismatch and leaves only the exponentially small genuine adiabatic
    tail transitions.
    """
    w, V = np.linalg.eigh(geo.H(u))
    # match each eigenvector to the diabatic index of largest |component|
    perm = [-1, -1, -1]
    cols = list(range(3))
    for i in range(3):
        # diabatic channel i: pick remaining column with largest |V[i,col]|
        best = max(cols, key=lambda c: abs(V[i, c]))
        perm[i] = best
        cols.remove(best)
    Vm = V[:, perm].astype(complex)
    # fix global phase of each column (largest component real positive)
    for i in range(3):
        k = np.argmax(np.abs(Vm[:, i]))
        Vm[:, i] *= np.conj(Vm[k, i]) / abs(Vm[k, i])
    return Vm


def transition_matrix(geo: Geometry, T: float = 60.0,
                       rtol: float = 1e-12, atol: float = 1e-13) -> np.ndarray:
    """
    Direct benchmark transition matrix in the adiabatic convention,

        P[x, j] = | (V_+^H U(T,-T) V_-)[j, x] |^2 ,

    with V_- , V_+ the instantaneous adiabatic frames at -T, +T.
    Convergence in T is exponential (genuine adiabatic tail transitions
    are gap-suppressed).
    """
    res = fundamental_matrix(geo, -T, T, rtol=rtol, atol=atol)
    Vm = adiabatic_frame(geo, -T)
    Vp = adiabatic_frame(geo, T)
    M = Vp.conj().T @ res.U @ Vm        # adiabatic scattering matrix
    P = np.abs(M.T) ** 2                # P[x, j] = |M[j, x]|^2
    return P


def unitarity_defect(geo: Geometry, T: float = 60.0, **kw) -> float:
    """Max |U^dag U - I| over the one-pass fundamental matrix."""
    res = fundamental_matrix(geo, -T, T, **kw)
    U = res.U
    return float(np.max(np.abs(U.conj().T @ U - np.eye(3))))


class _Result:
    """Lightweight propagator result holder."""
    U: np.ndarray
    t0: float
    t1: float
    sol = None
