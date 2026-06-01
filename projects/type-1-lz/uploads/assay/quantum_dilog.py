"""
Faddeev quantum dilogarithm Phi_b(z) and its semiclassical limit.

Faddeev's non-compact quantum dilogarithm
-----------------------------------------
For Im b > 0 (or |b|=1, b not a root of unity) and Im z inside the strip
|Im z| < |Im((b + b^{-1})/2)|, Faddeev's quantum dilogarithm is

    Phi_b(z) = prod_{n=0}^inf  (1 + e^{ 2 pi b   ( z + (n + 1/2)(b + b^{-1})/2 )}) /
                               (1 + e^{ 2 pi b^-1 ( z + (n + 1/2)(b + b^{-1})/2 )}) .

[Equivalent forms in the literature differ by sign conventions on the
exponent; we use the convention compatible with Allegretti's "Voros
symbols as cluster coordinates" so that the cluster mutation acts as

    Y -> Y . Phi_b(Y) / Phi_b(Y^{-1})   etc.]

Semiclassical limit
-------------------
For b -> 0+ (Planck constant epsilon = 2 pi b^2 -> 0) one has the
universal exponential form

    Phi_b(z) ~ exp( Li_2(-e^{2 pi b z}) / (2 pi i b^2) )       (b -> 0)

so the *classical limit* of the cluster Y-mutation is encoded in the
classical dilogarithm Li_2(-Y) where Y = exp(Voros symbol).  In the
WKB-Stokes setting the Voros symbol attached to a 1-cycle gamma is
V_gamma = (1/epsilon) oint_gamma sqrt(y) du = (1/epsilon) I_X, so
Y_gamma = exp(-V_gamma) is exponentially small at epsilon -> 0; the
mutation jump becomes

    log Y_new = log Y_old + Li_2(-Y_other_thru_arrow) / (2 pi i) .

For our application we set epsilon = 1 (the WKB pencil is exact, not
asymptotic) and Y_X = exp(I_X) with I_X purely imaginary in the
canonical sample, so |Y_X| = 1 and the mutation is unitary.

Module API
----------
Li2          : the principal-branch dilogarithm Li_2(z) for complex z.
phi_classical: the semiclassical Phi_b limit  exp(Li_2(-e^{2 pi z}) / (2 pi i))
phi_faddeev  : Faddeev's Phi_b(z) by the rapidly convergent infinite product.
mutation_Y   : cluster Y-mutation Y_k' = Y_k^{-1}, Y_j' = Y_j (1 + Y_k^{eps_jk})^{-eps_jk}.
mutation_X   : conjugate cluster X-mutation (Fock-Goncharov form).
"""

from __future__ import annotations
from typing import Iterable
import numpy as np
from scipy.special import spence


# ---------------------------------------------------------------------------
#  Classical dilogarithm (Spence/Jonquiere)
# ---------------------------------------------------------------------------
def _Li2_series(z: complex, n_terms: int = 400, tol: float = 1e-18) -> complex:
    """Maclaurin series for |z| <= 1/2."""
    s = 0.0 + 0.0j
    zk = z
    for k in range(1, n_terms + 1):
        term = zk / (k * k)
        s += term
        zk *= z
        if abs(term) < tol:
            break
    return s


def Li2(z) -> complex:
    """
    Principal-branch dilogarithm  Li_2(z) = -int_0^z log(1 - t)/t dt ,
    for complex z, by recursion-free dispatch over four regions:

      (R1)  |z| <= 1/2                  : Maclaurin
      (R2)  |1 - z| <= 1/2              : reflection through (1-z)
      (R3)  |z| > 2                     : inversion through 1/z
      (R4)  otherwise (|z|, |1-z| > 1/2): Landen through z/(z-1)

    Each branch reduces to (R1) in one step (no recursion).
    """
    z = complex(z)
    if z == 0:
        return 0.0 + 0.0j
    if z == 1:
        return np.pi ** 2 / 6.0 + 0.0j
    # (R1) direct series
    if abs(z) <= 0.5:
        return _Li2_series(z)
    # (R2) reflection: Li_2(z) = pi^2/6 - log(z) log(1-z) - Li_2(1-z)
    if abs(1.0 - z) <= 0.5:
        w = 1.0 - z
        if w == 0:
            return np.pi ** 2 / 6.0 + 0.0j
        return (np.pi ** 2 / 6.0 - np.log(z) * np.log(w) - _Li2_series(w))
    # (R3) inversion: Li_2(z) = -Li_2(1/z) - pi^2/6 - (1/2) log^2(-z)
    if abs(z) >= 2.0:
        w = 1.0 / z
        if abs(w) <= 0.5:
            inner = _Li2_series(w)
        else:
            # second-level Landen
            u = w / (w - 1.0)
            inner = -_Li2_series(u) - 0.5 * np.log(1.0 - w) ** 2
        return -inner - np.pi ** 2 / 6.0 - 0.5 * np.log(-z) ** 2
    # (R4) Landen: Li_2(z) = -Li_2(z/(z-1)) - (1/2) log^2(1-z)
    w = z / (z - 1.0)
    if abs(w) <= 0.5:
        inner = _Li2_series(w)
    else:
        # last resort: reflection on w
        u = 1.0 - w
        inner = (np.pi ** 2 / 6.0 - np.log(w) * np.log(u)
                 - _Li2_series(u))
    return -inner - 0.5 * np.log(1.0 - z) ** 2


# Real-z fast path through scipy (validated against the complex routine).
def Li2_real(x: float) -> float:
    """Real Li_2(x) via scipy.special.spence.  spence(y) = Li_2(1 - y)."""
    return float(spence(1.0 - x))


# ---------------------------------------------------------------------------
#  Faddeev quantum dilogarithm by the convergent infinite product
# ---------------------------------------------------------------------------
def phi_faddeev(z, b: complex, n_terms: int = 200, tol: float = 1e-15
                ) -> complex:
    """
    Faddeev Phi_b(z) by the convergent infinite product.

    Phi_b(z) = prod_n (1 + q_+ Q^{2n+1}) / (1 + q_- Q^{2n+1}^{-1}) ,

    where Q = e^{i pi b^2},   q_+ = e^{2 pi b z}, q_- = e^{2 pi b^{-1} z}.
    For Im(b^2) > 0 the product converges absolutely.

    For real b (the semiclassical limit) it diverges and the classical
    form `phi_classical` should be used.
    """
    z = complex(z)
    b = complex(b)
    if abs(b.imag) < 1e-12 and abs((b.conjugate() - b)) < 1e-12 and b.real ** 2 > 0:
        # purely real b -- product does not converge
        return phi_classical(z, b)
    Q = np.exp(1j * np.pi * b ** 2)               # |Q|<1 if Im(b^2)>0
    bb = b + 1.0 / b
    # Faddeev's defining product is most numerically stable in this form
    qp_base = np.exp(2 * np.pi * b * z)
    qm_base = np.exp(2 * np.pi * (1.0 / b) * z)
    result = 1.0 + 0.0j
    for n in range(n_terms):
        num = 1.0 + qp_base * Q ** (2 * n + 1)
        den = 1.0 + qm_base * Q ** (-(2 * n + 1))
        result *= num / den
        if abs(Q ** (2 * n + 1)) < tol:
            break
    return result


def phi_classical(z, b: complex = 0.0) -> complex:
    """
    Semiclassical Phi_b(z) limit:

        Phi_b(z) ~ exp( Li_2(-e^{2 pi b z}) / (2 pi i b^2) )       (b -> 0)

    With b = 1 (a fixed normalisation; the WKB pencil supplies the
    exact phase already, so the Planck weighting is absorbed elsewhere)
    this collapses to

        phi_classical(z) := exp( Li_2(-e^{2 pi z}) / (2 pi i) ) .

    This is the form used by mutation_Y to compute the modulus and phase
    of the Y-variable jump on a single Stokes-line crossing.
    """
    z = complex(z)
    arg = -np.exp(2 * np.pi * z)
    return np.exp(Li2(arg) / (2j * np.pi))


# ---------------------------------------------------------------------------
#  Cluster mutation: Y-pattern (Fock-Goncharov)
# ---------------------------------------------------------------------------
def mutation_Y(Y: np.ndarray, B: np.ndarray, k: int) -> np.ndarray:
    """
    Apply cluster Y-mutation at node k:

        Y_k' = Y_k^{-1}
        Y_j' = Y_j * (1 + Y_k^{sgn(B_jk)})^{- B_jk}    (j != k)

    Y   : shape (n,) complex vector of Y-variables.
    B   : shape (n, n) skew-symmetric exchange matrix (integer-valued).
    k   : node index to mutate.
    """
    Y = np.asarray(Y, dtype=complex).copy()
    Yk = Y[k]
    new = Y.copy()
    new[k] = 1.0 / Yk
    for j in range(len(Y)):
        if j == k:
            continue
        bjk = int(B[j, k])
        if bjk == 0:
            continue
        sgn = 1 if bjk > 0 else -1
        # numerically stable:  Y_j (1 + Y_k^sgn)^{-bjk}
        base = 1.0 + Yk ** sgn
        new[j] = Y[j] * base ** (-bjk)
    return new


def mutation_X(X: np.ndarray, B: np.ndarray, k: int) -> np.ndarray:
    """
    Same combinatorics as mutation_Y but recorded as the standard
    X-pattern (used to track signed monomial transformations of the
    Stokes matrices).  X_k -> X_k^{-1},
        X_j' = X_j  X_k^{max(B_jk, 0)} (1 + X_k)^{-B_jk}     if B_jk > 0
        X_j' = X_j  X_k^{max(-B_jk, 0)} (1 + X_k^{-1})^{-B_jk}  if B_jk < 0
    """
    X = np.asarray(X, dtype=complex).copy()
    Xk = X[k]
    new = X.copy()
    new[k] = 1.0 / Xk
    for j in range(len(X)):
        if j == k:
            continue
        bjk = int(B[j, k])
        if bjk == 0:
            continue
        if bjk > 0:
            new[j] = X[j] * Xk ** bjk * (1.0 + Xk) ** (-bjk)
        else:
            new[j] = X[j] * (1.0 + 1.0 / Xk) ** (-bjk)
    return new


# ---------------------------------------------------------------------------
#  Quantum-dilogarithm operator on a single 2x2 mutation
# ---------------------------------------------------------------------------
def stokes_matrix_2x2(p: float, phase: complex = 0.0) -> np.ndarray:
    """
    Reduced 2x2 Stokes/connection matrix for a single saddle crossing.

    On the semiclassical level (real b -> 0) the cluster mutation realises
    the 2-level Landau-Zener connection in the (i, j) block,

        T_X(p) = [[ sqrt(1-p)  -sqrt(p) e^{-i phase}],
                  [ sqrt(p) e^{i phase}    sqrt(1-p)]]

    with `p = exp(-|Im I_X|)` and `phase` the Stokes-corrected Voros phase
    `arg Phi_b(V_X)` carried by the mutation.  When `phase = 0` (purely
    imaginary Voros symbol, our canonical regime) this is the elementary
    Landau-Zener rotation.
    """
    c = np.sqrt(max(0.0, 1.0 - p))
    s = np.sqrt(max(0.0, p))
    e = np.exp(1j * phase)
    return np.array([[c, -s * np.conj(e)],
                     [s * e,    c]], dtype=complex)


__all__ = ["Li2", "Li2_real", "phi_faddeev", "phi_classical",
           "mutation_Y", "mutation_X", "stokes_matrix_2x2"]
