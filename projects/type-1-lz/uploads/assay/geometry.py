"""
Type-1, N=3 Landau-Zener: geometry auditor.

Family-invariant and member-dependent algebraic data, exactly as in the
primer (type1_n3_gluing_selector_primer_v6) and LZ_summary.

Conventions
-----------
Polynomials are stored as numpy coefficient arrays, HIGHEST degree first
(numpy.polyval / numpy.roots convention).

Parameters
----------
eps : real, strictly ordered  eps0 < eps1 < eps2
gam : nonzero real couplings
a   : real slopes
x   : preferred incoming channel index in {0,1,2}
"""

from __future__ import annotations
from dataclasses import dataclass
import numpy as np


# ---------------------------------------------------------------------------
#  Parameter container
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Params:
    eps: tuple
    gam: tuple
    a: tuple
    x: int = 0

    def __post_init__(self):
        e = np.asarray(self.eps, float)
        g = np.asarray(self.gam, float)
        if e.shape != (3,) or g.shape != (3,) or len(self.a) != 3:
            raise ValueError("N=3 only: eps, gam, a must each have length 3")
        if not (e[0] < e[1] < e[2]):
            raise ValueError("eps must be strictly ordered")
        if np.any(g == 0.0):
            raise ValueError("gam entries must be nonzero")
        if self.x not in (0, 1, 2):
            raise ValueError("x must be 0, 1 or 2")


# ---------------------------------------------------------------------------
#  Geometry: all algebraic data derived from a Params
# ---------------------------------------------------------------------------
class Geometry:
    """All p,n,m,q_u,Q4,W4,L_H,kappa_S data plus diabatic Hamiltonian."""

    def __init__(self, par: Params):
        self.par = par
        eps = np.asarray(par.eps, float)
        gam = np.asarray(par.gam, float)
        a = np.asarray(par.a, float)
        g2 = gam ** 2
        self.eps, self.gam, self.a, self.g2 = eps, gam, a, g2

        # --- p(lambda) = prod(lambda - eps_i), degree 3 (monic) ---------
        p2 = -eps.sum()
        p1 = eps[0]*eps[1] + eps[0]*eps[2] + eps[1]*eps[2]
        p0 = -eps[0]*eps[1]*eps[2]
        self.p = np.array([1.0, p2, p1, p0])
        self.p2, self.p1, self.p0 = p2, p1, p0

        # --- n(lambda), m(lambda): degree 2 -----------------------------
        # n_2 = sum g^2 ; n_1 = -sum g^2 sum_{b!=a} eps_b ; n_0 = sum g^2 prod_{b!=a} eps_b
        def quad_from(weights):
            c2 = weights.sum()
            c1 = -sum(weights[k] * (eps.sum() - eps[k]) for k in range(3))
            c0 = sum(weights[k] * np.prod([eps[b] for b in range(3) if b != k])
                     for k in range(3))
            return np.array([c2, c1, c0])

        self.n = quad_from(g2)
        self.m = quad_from(a * g2)
        self.n2, self.n1, self.n0 = self.n
        self.m2, self.m1, self.m0 = self.m

        # --- family member H: interpolation quadratic f_H(eps_a)=a_a ----
        # f_H(l) = alphaH l^2 + betaH l + gammaH
        V = np.vander(eps, 3)              # columns: l^2, l^1, l^0
        self.alphaH, self.betaH, self.gammaH = np.linalg.solve(V, a)

        # --- L_H(lambda) = betaH + alphaH*(sum eps - lambda) : linear ---
        sume = eps.sum()
        self.LH = np.array([-self.alphaH, self.betaH + self.alphaH * sume])

        # --- W4 = n p' - n' p  (quartic) --------------------------------
        pder = np.polyder(self.p)
        nder = np.polyder(self.n)
        self.W4 = np.polysub(np.polymul(self.n, pder), np.polymul(nder, self.p))
        self.W4 = _trim_to_degree(self.W4, 4)

        # --- Q4 = Disc_z(g_xi), explicit primer coefficients ------------
        n0, n1, n2 = self.n0, self.n1, self.n2
        q4 = n1**2 - 4*n0*n2
        q3 = 2*(-n0*n1 - 2*n0*n2*p2 + n1**2*p2 - n1*n2*p1 + 2*n2**2*p0)
        q2 = (-3*n0**2 - 6*n0*n2*p1 + n1**2*p2**2 + 6*n1*n2*p0
              - 2*n1*n2*p1*p2 + n2**2*p1**2)
        q1 = 2*(-n0**2*p2 - 2*n0*n1*p1 + n0*n1*p2**2 + n0*n2*p0
                - n0*n2*p1*p2 + 2*n1**2*p0 - n1*n2*p0*p2 + n2**2*p0*p1)
        q0 = -4*n0**2*p1 + n0**2*p2**2 + 4*n0*n1*p0 - 2*n0*n2*p0*p2 + n2**2*p0**2
        self.Q4 = np.array([q4, q3, q2, q1, q0])

        # --- kappa_S = prod g^2 * prod_{a<b}(eps_a-eps_b)^2 -------------
        disc_eps = ((eps[0]-eps[1])*(eps[0]-eps[2])*(eps[1]-eps[2]))**2
        self.kappaS = float(np.prod(g2) * disc_eps)

        # --- diabatic Hamiltonian  H(u) = H0 + u * diag(a) --------------
        H0 = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                if i != j:
                    H0[i, j] = gam[i]*gam[j]*(a[i]-a[j])/(eps[i]-eps[j])
        for i in range(3):
            H0[i, i] = -sum(g2[j]*(a[i]-a[j])/(eps[i]-eps[j])
                            for j in range(3) if j != i)
        self.H0 = H0
        self.A = np.diag(a)

    # -- diabatic Hamiltonian at base point u ---------------------------
    def H(self, u: float) -> np.ndarray:
        return self.H0 + u * self.A

    # -- spectral polynomial q_u(lambda) = u p - n ----------------------
    def q_u(self, u: float) -> np.ndarray:
        return np.array([u, u*self.p2 - self.n2,
                         u*self.p1 - self.n1, u*self.p0 - self.n0])

    def lambdas(self, u: float) -> np.ndarray:
        """Three roots of q_u (unsorted)."""
        return np.roots(self.q_u(u))

    # -- member-dependent coordinates -----------------------------------
    def y_coord(self, lam):
        """Phase coordinate y = Q4 L_H^2 / p^2."""
        return (np.polyval(self.Q4, lam) * np.polyval(self.LH, lam)**2
                / np.polyval(self.p, lam)**2)

    def r_coord(self, lam):
        """Transport coordinate r = W4 / L_H^4."""
        return np.polyval(self.W4, lam) / np.polyval(self.LH, lam)**4

    # -- root sets -------------------------------------------------------
    def Q4_roots(self) -> np.ndarray:
        return np.roots(self.Q4)

    def W4_roots(self) -> np.ndarray:
        return np.roots(self.W4)

    def LH_root(self) -> float:
        return float(-self.LH[1] / self.LH[0])  # lambda where L_H = 0

    # -- projected conjugate-pair Q4 windows ----------------------------
    def windows(self, theta: float = 0.0, tol: float = 1e-7):
        """
        Identify the two projected conjugate-pair Q4 windows.

        Returns a list of dicts ordered by the projection coordinate
        tau(z;theta) = Re(e^{-i theta} z), each with:
          roots   : the conjugate pair (lambda values)
          u_vals  : u(lambda) for those roots
          tau     : projection coordinate of the pair (real)
          u_center: real-axis window center  Re(mean(u_vals))
        """
        roots = self.Q4_roots()
        e_ith = np.exp(-1j*theta)
        tau = np.real(e_ith*roots)
        used = [False]*4
        pairs = []
        order = np.argsort(np.imag(roots))
        for idx in order:
            if used[idx]:
                continue
            # find conjugate partner
            best, bestd = None, np.inf
            for jdx in range(4):
                if jdx == idx or used[jdx]:
                    continue
                d = abs(roots[idx] - np.conj(roots[jdx]))
                if d < bestd:
                    best, bestd = jdx, d
            if best is None:
                continue
            used[idx] = used[best] = True
            pr = np.array([roots[idx], roots[best]])
            uv = np.array([self._u_of(l) for l in pr])
            pairs.append({
                "roots": pr,
                "u_vals": uv,
                "tau": float(np.mean(np.real(e_ith*pr))),
                "u_center": float(np.mean(np.real(uv))),
                "conj_defect": float(bestd),
            })
        pairs.sort(key=lambda d: d["tau"])
        return pairs

    def _u_of(self, lam):
        return np.polyval(self.n, lam) / np.polyval(self.p, lam)


# ---------------------------------------------------------------------------
#  helpers
# ---------------------------------------------------------------------------
def _trim_to_degree(coeffs, deg):
    """Pad/trim a numpy coeff array (high-first) to exactly degree `deg`."""
    c = np.asarray(coeffs, float)
    c = np.trim_zeros(c, 'f')
    if c.size == 0:
        c = np.array([0.0])
    if c.size < deg + 1:
        c = np.concatenate([np.zeros(deg + 1 - c.size), c])
    return c


def cross_check(geo: Geometry) -> dict:
    """Independent verification of the explicit Q4/W4 formulas."""
    # W4 via direct definition n p' - n' p (already used) vs primer 2.14
    n0, n1, n2 = geo.n0, geo.n1, geo.n2
    p0, p1, p2 = geo.p0, geo.p1, geo.p2
    W4_primer = np.array([
        n2, 2*n1, 3*n0 + n1*p2 - n2*p1,
        2*(n0*p2 - n2*p0), n0*p1 - n1*p0])
    # Q4 via discriminant of g_xi at a few sample points
    xs = np.array([0.3, -1.1, 2.4, 5.0])
    q4_explicit = np.polyval(geo.Q4, xs)
    q4_disc = np.array([_Q4_via_disc(geo, x) for x in xs])
    return {
        "W4_formula_max_abs_diff": float(np.max(np.abs(geo.W4 - W4_primer))),
        "Q4_disc_max_abs_diff": float(np.max(np.abs(q4_explicit - q4_disc))),
    }


def _Q4_via_disc(geo: Geometry, xi: float) -> float:
    """Q4(xi) = Disc_z(g_xi) computed from g_xi = (n(z)p(xi)-n(xi)p(z))/(z-xi)."""
    pxi = np.polyval(geo.p, xi)
    nxi = np.polyval(geo.n, xi)
    # g_xi(z) numerator = n(z) p(xi) - n(xi) p(z); divide by (z - xi)
    num = np.polysub(pxi*np.concatenate([[0.0], geo.n]), nxi*geo.p)
    g, rem = np.polydiv(num, np.array([1.0, -xi]))
    g = _trim_to_degree(g, 2)
    A, B, C = g
    return float(B*B - 4*A*C)
