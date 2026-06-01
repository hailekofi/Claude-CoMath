"""
Interaction-picture propagators (fast engine).

Lab-frame integration of the Type-1 evolution is O(T^2) in step count
because the amplitudes oscillate at frequency proportional to |u|.  The
interaction picture factors out the diagonal phase, leaving slowly
varying amplitudes g that an adaptive integrator can step through the
tails in O(1) steps.

Two propagators:

  benchmark_ip   -- diabatic IP.  H = H0 + u A,  diagonal phase
                    phi_i(u) = H0_ii u + a_i u^2/2  is *analytic*
                    (a polynomial), so no phase integration is needed.
                    g_i' = -i sum_j H0_ij e^{i(phi_i-phi_j)} g_j .

  propagate_ad_ip -- gauged adiabatic IP (the constructor pathway).
                    State carries the two interior eigenvalues
                    lam_1, lam_2 (lam_0 = n2/u - p2 - lam_1 - lam_2 is
                    recovered exactly, avoiding its u=0 singularity) and
                    the phases phi_i = int (E_i - E_x).
                    g_i' = - sum_j S_ij e^{i(phi_j-phi_i)} g_j .

For a single pass U_f(b,a) = Phi(b) U_g(b,a) Phi(a)^{-1}, Phi diagonal.
"""

from __future__ import annotations
import numpy as np
from scipy.integrate import solve_ivp

from .geometry import Geometry
from .adiabatic import labelled_lambdas


# ===========================================================================
#  Diabatic interaction-picture benchmark
# ===========================================================================
def _benchmark_ip_g(geo: Geometry, T: float, rtol: float, atol: float):
    """Integrate the IP fundamental matrix g over [-T, T]; return (g, phi_fn)."""
    H0 = geo.H0
    a = geo.a
    H0d = np.diag(H0).copy()
    H0off = H0 - np.diag(H0d)

    def phi(u):                                  # analytic diagonal phase
        return H0d * u + a * (u * u) / 2.0

    def rhs(u, y):
        g = y.reshape(3, 3)
        ph = phi(u)
        # Vtilde_ij = -i H0off_ij exp(i(phi_i - phi_j))
        E = np.exp(1j * (ph[:, None] - ph[None, :]))
        Vt = -1j * H0off * E
        return (Vt @ g).ravel()

    sol = solve_ivp(rhs, (-T, T), np.eye(3, dtype=complex).ravel(),
                    method="DOP853", rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(f"benchmark_ip failed: {sol.message}")
    return sol.y[:, -1].reshape(3, 3), phi


def benchmark_ip_U(geo: Geometry, T: float = 60.0,
                   rtol: float = 1e-12, atol: float = 1e-13) -> np.ndarray:
    """Lab-frame one-pass diabatic fundamental matrix U(T,-T)."""
    g, phi = _benchmark_ip_g(geo, T, rtol, atol)
    Phi_T = np.exp(-1j * phi(T))
    Phi_mT = np.exp(-1j * phi(-T))
    # U_f = Phi(T) g Phi(-T)^{-1}
    return (Phi_T[:, None]) * g * (np.exp(1j * phi(-T))[None, :])


# ===========================================================================
#  Adiabatic interaction-picture propagator (constructor pathway)
# ===========================================================================
def _ad_ip_rhs_factory(geo: Geometry, x: int):
    eps = geo.eps
    g2 = geo.g2
    a = geo.a
    n2 = geo.n2
    p2 = geo.p2
    sgn = np.sign(geo.gam)

    def rhs(u, y):
        g = y[0:9].reshape(3, 3)
        lam1, lam2 = y[9], y[10]
        phi = y[11:14].real
        lam0 = n2 / u - p2 - lam1 - lam2
        lam = np.array([lam0, lam1, lam2])

        d0 = lam[:, None] - eps[None, :]          # lam_i - eps_j
        Sig2 = np.sum(g2[None, :] / d0 ** 2, axis=1)
        E = np.sum(a[None, :] * g2[None, :] / d0, axis=1)
        delta = np.array([sgn[0] * (1.0 if u < 0 else -1.0), sgn[1], sgn[2]])
        Gam = delta / np.sqrt(Sig2)

        # S_ij and the IP coefficient Vtilde
        Vt = np.zeros((3, 3), complex)
        for i in range(3):
            for j in range(3):
                if i != j:
                    S = Gam[i] * Gam[j] / (lam[i] - lam[j])
                    Vt[i, j] = -S * np.exp(1j * (phi[j] - phi[i]))
        dg = (Vt @ g).ravel()
        dlam1 = -1.0 / Sig2[1]
        dlam2 = -1.0 / Sig2[2]
        dphi = (E - E[x]).astype(complex)
        out = np.empty(14, complex)
        out[0:9] = dg
        out[9] = dlam1
        out[10] = dlam2
        out[11:14] = dphi
        return out

    return rhs


def propagate_ad_ip(geo: Geometry, a: float, b: float, x: int,
                    rtol: float = 1e-12, atol: float = 1e-13) -> np.ndarray:
    """
    Independent adiabatic segment IVP in the interaction picture.
    Returns the lab-frame fundamental matrix U_f(b,a) = Phi(b) g(b),
    with the segment phase reference phi(a) = 0.
    """
    lam = labelled_lambdas(geo, a)               # one root solve at segment start
    y0 = np.zeros(14, complex)
    y0[0:9] = np.eye(3, dtype=complex).ravel()
    y0[9] = lam[1]
    y0[10] = lam[2]
    # phi(a) = 0
    rhs = _ad_ip_rhs_factory(geo, x)
    sol = solve_ivp(rhs, (a, b), y0, method="DOP853", rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(f"propagate_ad_ip [{a},{b}] failed: {sol.message}")
    yf = sol.y[:, -1]
    g = yf[0:9].reshape(3, 3)
    phi = yf[11:14].real
    return (np.exp(1j * phi)[:, None]) * g       # U_f = diag(e^{i phi(b)}) g(b)


# ===========================================================================
#  Adiabatic IP propagator along a complex-u contour (Level 2)
# ===========================================================================
def propagate_ad_ip_contour(geo: Geometry, contour, x: int,
                            rtol: float = 1e-11, atol: float = 1e-12):
    """
    Integrate the gauged adiabatic IP evolution along a complex-u contour.

    `contour(s)` for s in [0,1] must return (u, du/ds).  The contour starts
    and ends on the real axis.  Eigenvalues lam_1,lam_2 and the transport
    normalisations sqrtSig_i (sqrtSig_i^2 = Sigma2_i, continuous branch) are
    carried as ODE state so branch continuity is maintained off-axis;
    lam_0 = n2/u - p2 - lam_1 - lam_2 is recovered exactly.

    Returns U_f = diag(e^{i phi(end)}) g(end), the lab-frame contour
    propagator from contour(0) to contour(1).
    """
    eps = geo.eps
    g2 = geo.g2
    a = geo.a
    n2 = geo.n2
    p2 = geo.p2
    sgn = np.sign(geo.gam)

    u0, _ = contour(0.0)
    lam = labelled_lambdas(geo, complex(u0).real if abs(complex(u0).imag) < 1e-14
                           else u0)
    d0 = lam[:, None] - eps[None, :]
    Sig2_0 = np.sum(g2[None, :] / d0 ** 2, axis=1)
    delta0 = np.array([sgn[0] * (1.0 if u0.real < 0 else -1.0), sgn[1], sgn[2]])
    sqrtSig0 = delta0 * np.sqrt(Sig2_0)            # real-axis branch at start

    y0 = np.zeros(17, complex)
    y0[0:9] = np.eye(3, dtype=complex).ravel()
    y0[9] = lam[1]
    y0[10] = lam[2]
    y0[11:14] = sqrtSig0
    # phi(start) = 0  -> y0[14:17] = 0

    def rhs(s, y):
        u, duds = contour(s)
        g = y[0:9].reshape(3, 3)
        lam1, lam2 = y[9], y[10]
        ss = y[11:14]
        phi = y[14:17]
        lam0 = n2 / u - p2 - lam1 - lam2
        lam = np.array([lam0, lam1, lam2])

        d0 = lam[:, None] - eps[None, :]
        Sig2 = np.sum(g2[None, :] / d0 ** 2, axis=1)
        Sig3 = np.sum(g2[None, :] / d0 ** 3, axis=1)
        E = np.sum(a[None, :] * g2[None, :] / d0, axis=1)

        ld1 = -1.0 / Sig2[1]
        ld2 = -1.0 / Sig2[2]
        ld0 = -n2 / u ** 2 - ld1 - ld2
        ldot = np.array([ld0, ld1, ld2])
        # d(sqrtSig_i)/du = -Sig3_i lam_i' / sqrtSig_i
        ssdot = -Sig3 * ldot / ss

        Gam = 1.0 / ss
        Vt = np.zeros((3, 3), complex)
        for i in range(3):
            for j in range(3):
                if i != j:
                    S = Gam[i] * Gam[j] / (lam[i] - lam[j])
                    Vt[i, j] = -S * np.exp(1j * (phi[j] - phi[i]))
        dg = (Vt @ g).ravel()
        dphi = E - E[x]

        out = np.empty(17, complex)
        out[0:9] = dg * duds
        out[9] = ld1 * duds
        out[10] = ld2 * duds
        out[11:14] = ssdot * duds
        out[14:17] = dphi * duds
        return out

    sol = solve_ivp(rhs, (0.0, 1.0), y0, method="DOP853", rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(f"contour propagation failed: {sol.message}")
    yf = sol.y[:, -1]
    g = yf[0:9].reshape(3, 3)
    phi = yf[14:17]
    return (np.exp(1j * phi)[:, None]) * g
