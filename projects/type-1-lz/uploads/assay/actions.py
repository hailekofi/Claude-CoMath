"""
Window actions: the Abelian/elliptic periods that set the closed-form
Landau-Zener exponents.

The phase one-form is  omega = sqrt(y) du  with  y = Q4 L_H^2 / p^2,
so  sqrt(y) = sqrt(Q4) L_H / p  and, since  du = u' dlam = -(W4/p^2) dlam,

    omega = - sqrt(Q4) L_H W4 / p^3  dlam .

The action attached to a projected Q4 window X = {q, qbar} is the period

    I_X = oint_X omega ,

the contour encircling the conjugate pair of Q4 roots of that window.
sqrt(Q4) is tracked continuously along the contour; encircling two of the
four branch points returns it to its branch, so I_X is a genuine period
of the elliptic curve mu^2 = Q4.  The Landau-Zener exponent is
delta_X = |Im I_X| / (2 pi)  (Dykhne form); the window transition
probability is  p_X = exp(-2 pi delta_X) = exp(-|Im I_X|).
"""

from __future__ import annotations
import numpy as np

from .geometry import Geometry


def _omega_lambda(geo: Geometry, lam, sqrtQ4):
    """omega / dlam = - sqrt(Q4) L_H W4 / p^3  (sqrtQ4 supplied, branch-tracked)."""
    LH = np.polyval(geo.LH, lam)
    W4 = np.polyval(geo.W4, lam)
    p = np.polyval(geo.p, lam)
    return -sqrtQ4 * LH * W4 / p ** 3


def window_action(geo: Geometry, win: dict, n: int = 4000) -> dict:
    """
    Period  I_X = oint_X omega  for one projected Q4 window.

    The contour is a circle in the lambda-plane centred on the midpoint of
    the conjugate Q4 pair, with radius 1.6x half the pair separation, so it
    encloses exactly those two branch points.  sqrt(Q4) is continued
    continuously around the loop.
    """
    q = win["roots"]
    centre = 0.5 * (q[0] + q[1])
    radius = 0.8 * abs(q[0] - q[1]) * 2.0          # enclose both, clear of others
    th = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    lam = centre + radius * np.exp(1j * th)
    dlam = 1j * radius * np.exp(1j * th) * (2.0 * np.pi / n)

    # branch-tracked sqrt(Q4) along the loop
    Q4v = np.polyval(geo.Q4, lam)
    s = np.sqrt(Q4v[0])
    sq = np.empty(n, complex)
    for k in range(n):
        cand = np.sqrt(Q4v[k])
        s = cand if abs(cand - s) < abs(-cand - s) else -cand
        sq[k] = s

    integrand = _omega_lambda(geo, lam, sq)
    I_X = np.sum(integrand * dlam)
    # closure check: sqrt(Q4) must return to its starting branch
    closure = abs(sq[-1] * np.exp(1j * (th[-1] + 2 * np.pi / n - th[0])) - sq[0])
    delta = abs(I_X.imag) / (2.0 * np.pi)
    return {
        "I_X": I_X,
        "delta": delta,
        "p_transition": float(np.exp(-abs(I_X.imag))),
        "centre": centre, "radius": radius,
        "branch_closure_defect": float(abs(np.sqrt(Q4v[-1]) -
                                           (sq[-1] if abs(sq[-1]-np.sqrt(Q4v[-1]))
                                            < abs(sq[-1]+np.sqrt(Q4v[-1])) else -sq[-1]))),
    }


def all_window_actions(geo: Geometry) -> list:
    """Window actions for both projected Q4 windows (theta=0)."""
    return [window_action(geo, w) for w in geo.windows()]
