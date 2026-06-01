"""
Gauged adiabatic coefficient matrix for the Type-1, N=3 evolution.

The gauged adiabatic amplitudes solve (LZ_summary / primer eq. 1.5)

    f_i'(u) = - sum_{j!=i} S_ij(u) f_j(u) + i (E_i(u) - E_x(u)) f_i(u),

i.e.  dF/du = A(u) F  with

    A[i,i] = i (E_i - E_x),     A[i,j] = - S_ij   (i != j).

This is the *constructor* coefficient matrix used by the segment-ODE and
lateral-contour assays.  It is built independently of the benchmark
propagator (benchmark.py), so the two pathways are disjoint.

Derived quantities (LZ_summary):
    Sigma2_i = sum_j gam_j^2 / (lam_i - eps_j)^2
    Gamma_i  = delta_i / sqrt(Sigma2_i)
    S_ij     = Gamma_i Gamma_j / (lam_i - lam_j)        (S_ii = 0)
    E_i      = sum_j a_j gam_j^2 / (lam_i - eps_j)  =  m(lam_i)/p(lam_i)
    delta_0  = (1 - 2*theta(u)) sign(gam_0),   delta_{1,2} = sign(gam_{1,2})

For real ordered Type-1 data the three roots of q_u are real and interlace
the eps_j, so they are labelled unambiguously by interval:
    lam_1 in (eps_0, eps_1),  lam_2 in (eps_1, eps_2),  lam_0 the exterior.
S_ij is smooth through u = 0 even though lam_0 -> +-infinity there.
"""

from __future__ import annotations
import numpy as np

from .geometry import Geometry


def labelled_lambdas(geo: Geometry, u: float, hint=None) -> np.ndarray:
    """
    Return [lam_0, lam_1, lam_2] with lam_0 exterior, lam_1, lam_2 interior.

    For real u the roots are real and labelled by interval.  If `hint` (a
    previous [lam_0,lam_1,lam_2]) is given, label by nearest continuation
    instead -- needed for complex-u lateral contours.
    """
    if u == 0.0:
        u = 1e-30          # q_u keeps degree 3; lam_0 -> finite-but-huge
    roots = geo.lambdas(u)
    if hint is not None:
        out = np.empty(3, complex)
        pool = list(roots)
        for k in range(3):
            j = int(np.argmin([abs(r - hint[k]) for r in pool]))
            out[k] = pool.pop(j)
        return out
    # real-u interval labelling
    rr = np.sort(np.real(roots))
    e0, e1, e2 = geo.eps
    lam1 = lam2 = lam0 = None
    for r in rr:
        if e0 < r < e1 and lam1 is None:
            lam1 = r
        elif e1 < r < e2 and lam2 is None:
            lam2 = r
        else:
            lam0 = r
    if lam1 is None or lam2 is None or lam0 is None:
        # fallback: degenerate / outside-bracket -- use sorted order
        return np.array([roots[np.argmax(np.abs(roots))]]
                        + [r for r in roots], dtype=complex)[:3]
    return np.array([lam0, lam1, lam2], dtype=complex)


def adiabatic_data(geo: Geometry, u: float, hint=None) -> dict:
    """All derived adiabatic quantities at base point u."""
    lam = labelled_lambdas(geo, u, hint)
    eps = geo.eps
    g2 = geo.g2
    a = geo.a

    # Sigma2_i = sum g^2/(lam_i-eps)^2
    Sig2 = np.array([np.sum(g2 / (lam[i] - eps) ** 2) for i in range(3)])
    # sign factors
    sgn = np.sign(geo.gam)
    delta = np.array([
        (1.0 - 2.0 * (1.0 if u > 0 else 0.0)) * sgn[0],
        sgn[1], sgn[2]], dtype=complex)
    Gam = delta / np.sqrt(Sig2 + 0j)
    # E_i
    E = np.array([np.sum(a * g2 / (lam[i] - eps)) for i in range(3)])
    # S_ij
    S = np.zeros((3, 3), complex)
    for i in range(3):
        for j in range(3):
            if i != j:
                S[i, j] = Gam[i] * Gam[j] / (lam[i] - lam[j])
    return {"lam": lam, "Sigma2": Sig2, "Gamma": Gam, "E": E, "S": S,
            "delta": delta}


def coeff_matrix(geo: Geometry, u: float, x: int | None = None,
                 hint=None) -> np.ndarray:
    """
    The 3x3 gauged adiabatic coefficient matrix A(u):
        A[i,i] = i (E_i - E_x),   A[i,j] = -S_ij.
    """
    if x is None:
        x = geo.par.x
    d = adiabatic_data(geo, u, hint)
    E, S = d["E"], d["S"]
    A = -S.astype(complex)
    for i in range(3):
        A[i, i] = 1j * (E[i] - E[x])
    return A


def coeff_from_tracked(geo: Geometry, u: complex, lam: np.ndarray,
                       Gamma: np.ndarray, x: int) -> np.ndarray:
    """
    Gauged adiabatic coefficient matrix from externally-tracked (lam, Gamma).

    Used by the complex-contour integrator (lateral.py), which evolves the
    eigenvalues lam_i and the transport normalisations Gamma_i = 1/sqrt(Sig2_i)
    as ODE state so that branch continuity is maintained off the real axis.
    """
    eps, g2, a = geo.eps, geo.g2, geo.a
    E = np.array([np.sum(a * g2 / (lam[i] - eps)) for i in range(3)])
    A = np.zeros((3, 3), complex)
    for i in range(3):
        for j in range(3):
            if i != j:
                A[i, j] = -Gamma[i] * Gamma[j] / (lam[i] - lam[j])
    for i in range(3):
        A[i, i] = 1j * (E[i] - E[x])
    return A


def lambda_dot(geo: Geometry, u: complex, lam: np.ndarray) -> np.ndarray:
    """d(lam_i)/du from implicit differentiation of q_u(lam_i)=0."""
    p, n = geo.p, geo.n
    pd, nd = np.polyder(p), np.polyder(n)
    out = np.empty(3, complex)
    for i in range(3):
        dq_du = np.polyval(p, lam[i])               # ∂q/∂u
        dq_dl = u * np.polyval(pd, lam[i]) - np.polyval(nd, lam[i])  # ∂q/∂λ
        out[i] = -dq_du / dq_dl
    return out


def sqrtSig_dot(geo: Geometry, lam: np.ndarray, lamdot: np.ndarray,
                sqrtSig: np.ndarray) -> np.ndarray:
    """
    d/du of a continuous branch sqrtSig_i with sqrtSig_i^2 = Sigma2_i.
    Sigma2_i = sum g^2/(lam-eps)^2 ;  (Sigma2_i)' = -2 Sigma3_i lam_i' .
    """
    eps, g2 = geo.eps, geo.g2
    out = np.empty(3, complex)
    for i in range(3):
        Sig3 = np.sum(g2 / (lam[i] - eps) ** 3)
        out[i] = -Sig3 * lamdot[i] / sqrtSig[i]
    return out


def cauchy_frame(geo: Geometry, u: float, hint=None) -> np.ndarray:
    """
    Cauchy / pole frame  V[i,alpha] = Gamma_i gam_alpha / (lam_i - eps_alpha).

    Rows index adiabatic sheet i, columns index diabatic channel alpha.
    V is the change of basis between the adiabatic amplitudes f_i and the
    diabatic components; it is the algebraic seed of the Level-3 frame maps.
    """
    d = adiabatic_data(geo, u, hint)
    lam, Gam = d["lam"], d["Gamma"]
    V = np.empty((3, 3), complex)
    for i in range(3):
        for al in range(3):
            V[i, al] = Gam[i] * geo.gam[al] / (lam[i] - geo.eps[al])
    return V
