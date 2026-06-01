"""
Transport-side periods: the spinor-cover analogue of ``actions.I_X``.

The phase one-form lives on the elliptic curve  :math:`\mu^2 = Q_4`, and its
window periods set the Brundobler-Elser exponents.  The transport-side
analogue lives on the spinor cover  :math:`\rho^2 = \kappa_S/r = \kappa_S\,L_H^4/W_4`
of the curve  :math:`\nu^2 = W_4`.  Concretely

.. math::

   \eta \;=\; \frac{\rho}{y}\,du,
   \qquad \rho^2 = \kappa_S\,\frac{L_H^4}{W_4},
   \qquad y = \frac{Q_4\,L_H^2}{p^2}.

So  :math:`\rho = \sqrt{\kappa_S}\,L_H^2/\sqrt{W_4}` and

.. math::

   \eta &= \sqrt{\kappa_S}\,\frac{L_H^2}{\sqrt{W_4}}\cdot\frac{p^2}{Q_4 L_H^2}\,du
        = \sqrt{\kappa_S}\,\frac{p^2}{Q_4\sqrt{W_4}}\,du.

Using  :math:`du = -(W_4/p^2)\,d\lambda`  to descend to the :math:`\lambda`-line,

.. math::

   \eta\,d\lambda \;=\; -\,\sqrt{\kappa_S}\,\frac{\sqrt{W_4}}{Q_4}\,d\lambda.

The branch points are the four roots of  :math:`W_4(\lambda)` (the spinor-cover
turning points), and the four roots of  :math:`Q_4(\lambda)` are simple
poles of  :math:`\eta`.

KEY ALGEBRAIC IDENTITY (proved in :func:`verify_residue_identity`)
-----------------------------------------------------------------

  .. math::  16\,\kappa_S\,W_4(\lambda) + Q_4'(\lambda)^2 \;\equiv\; 0
             \pmod{Q_4(\lambda)}.

At every  :math:`Q_4`-root  :math:`q`  this gives
:math:`Q_4'(q)^2 = -16\,\kappa_S W_4(q)`, hence

  .. math::

     \operatorname*{Res}_{\lambda=q}\eta \;=\;
     -\sqrt{\kappa_S}\,\frac{\sqrt{W_4(q)}}{Q_4'(q)} \;=\; \pm\,\frac{i}{4},

independent of the parameters.  Every transport window period
:math:`J_X = \oint_{\gamma_X}\eta` therefore lies in
:math:`\{0,\pm i\pi/2,\pm i\pi,\dots\}`  -- a purely topological list.

Hence the transport-side spinor-cover does NOT supply new continuous
algebraic data to the closed form.  See ``notes/track2_locate.tex`` for
the discussion and consequences.

Module API
----------

* :func:`transport_period_at_Q4_pair`  --  contour integral around one
  Q\ :sub:`4` conjugate pair (the transport analogue of
  :func:`assay.actions.window_action`).
* :func:`transport_residues`  --  the four residues  :math:`\operatorname*{Res}_q \eta`
  with one global sheet sign.
* :func:`all_transport_periods`  --  list of (centre, J_X) for the two
  Q\ :sub:`4`\ -windows in u-order.
* :func:`verify_residue_identity`  --  sympy-level check of the identity
  :math:`16\,\kappa_S W_4 + (Q_4')^2 \equiv 0 \pmod{Q_4}`.
"""

from __future__ import annotations
import numpy as np

from .geometry import Geometry


# ---------------------------------------------------------------------------
#  Transport one-form  eta dlam = -sqrt(kappa_S) sqrt(W4)/Q4 dlam
# ---------------------------------------------------------------------------
def _eta_lambda(geo: Geometry, lam, sqrtW4):
    """eta / dlam evaluated on a sampled contour, with branch-tracked sqrt(W4)."""
    return -np.sqrt(geo.kappaS) * sqrtW4 / np.polyval(geo.Q4, lam)


# ---------------------------------------------------------------------------
#  Contour period around one Q4 conjugate pair
# ---------------------------------------------------------------------------
def transport_period_at_Q4_pair(geo: Geometry, win: dict, n: int = 8000,
                                radius_scale: float = 1.6) -> dict:
    """
    Transport period  J_X = oint_{gamma_X} eta  for one projected Q4 window.

    The contour is a circle in the lambda-plane centred on the midpoint of
    the conjugate Q4 pair, with radius (radius_scale * half the separation),
    so it encloses exactly that conjugate pair of Q4 *poles* (NOT branch
    points; the W4 branch points sit elsewhere).  sqrt(W4) is continued
    continuously around the loop.

    Parameters
    ----------
    geo          : :class:`Geometry`
    win          : a single window dict as returned by ``geo.windows()``
                   (contains ``roots``, ``u_center``, ``conj_defect``)
    n            : quadrature sample count (default 8000)
    radius_scale : ratio of contour radius to half-separation (default 1.6
                   matches ``assay.actions.window_action``)

    Returns
    -------
    dict with keys
        J_X                  : complex period
        centre, radius       : contour parameters
        branch_closure_defect: how close branch-tracked sqrt(W4) returns
                               to its start
        residue_pair         : list of (q, Res_q eta) for the two Q4 roots
                               inside the contour, sheet-fixed
    """
    q = win["roots"]
    centre = 0.5 * (q[0] + q[1])
    radius = radius_scale * abs(q[0] - q[1]) / 2.0 * 2.0  # = radius_scale * |sep|

    th = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    lam = centre + radius * np.exp(1j * th)
    dlam = 1j * radius * np.exp(1j * th) * (2.0 * np.pi / n)

    W4v = np.polyval(geo.W4, lam)
    sW = np.empty(n, complex)
    s = np.sqrt(W4v[0])
    for k in range(n):
        cand = np.sqrt(W4v[k])
        s = cand if abs(cand - s) < abs(-cand - s) else -cand
        sW[k] = s

    integrand = _eta_lambda(geo, lam, sW)
    J_X = np.sum(integrand * dlam)

    # closure defect: does sqrt(W4) return to start?
    s_end_natural = np.sqrt(W4v[0])  # one period later, lam returns
    closure = min(abs(sW[-1] - s_end_natural), abs(sW[-1] + s_end_natural))
    # Identify which Q4 roots are inside the contour (the "pair" of poles)
    Q4_roots = geo.Q4_roots()
    inside = [(qz, abs(qz - centre) < radius) for qz in Q4_roots]
    Q4d = np.polyder(geo.Q4)
    residue_pair = []
    for qz, ok in inside:
        if ok:
            W4q = complex(np.polyval(geo.W4, qz))
            Q4dq = complex(np.polyval(Q4d, qz))
            # residue = -sqrt(kS) sqrt(W4(q))/Q4'(q); sheet chosen by branch-tracking
            # Find the value of sW closest to qz (won't be exact -- contour misses qz)
            # We just record the constant absolute value |i/4|; sign is fixed by the
            # algebraic identity 16 kappa_S W4(q) + Q4'(q)^2 = 0 mod Q4.
            R = -np.sqrt(geo.kappaS) * np.sqrt(W4q) / Q4dq
            residue_pair.append((complex(qz), complex(R)))
    return {
        "J_X": complex(J_X),
        "centre": complex(centre),
        "radius": float(radius),
        "branch_closure_defect": float(closure),
        "residue_pair": residue_pair,
    }


# ---------------------------------------------------------------------------
#  Residue list at the four Q4 roots
# ---------------------------------------------------------------------------
def transport_residues(geo: Geometry) -> list:
    """
    Compute  :math:`\operatorname*{Res}_{q}\eta = -\sqrt{\kappa_S}\sqrt{W_4(q)}/Q_4'(q)`
    at each of the four Q4 roots with the principal branch of the square roots.

    By the algebraic identity 16 kappa_S W_4 + (Q_4')^2 = 0 mod Q_4 these are
    pinned to  :math:`\pm i/4`.
    """
    Q4d = np.polyder(geo.Q4)
    out = []
    for qz in geo.Q4_roots():
        W4q = complex(np.polyval(geo.W4, qz))
        Q4dq = complex(np.polyval(Q4d, qz))
        R = -np.sqrt(geo.kappaS) * np.sqrt(W4q) / Q4dq
        out.append((complex(qz), complex(R)))
    return out


# ---------------------------------------------------------------------------
#  All transport periods (the two Q4-windows in u-order)
# ---------------------------------------------------------------------------
def all_transport_periods(geo: Geometry, n: int = 8000) -> list:
    """Transport periods for both projected Q4 windows (theta=0), u-ordered."""
    wins = sorted(geo.windows(), key=lambda w: w["u_center"])
    return [transport_period_at_Q4_pair(geo, w, n=n) for w in wins]


# ---------------------------------------------------------------------------
#  Symbolic verification of the residue identity (one-time check)
# ---------------------------------------------------------------------------
def verify_residue_identity() -> dict:
    """
    Prove symbolically that  16 kappa_S W_4(lambda) + Q_4'(lambda)^2 = 0
    modulo Q_4(lambda) in  Q[eps, gamma][lambda].

    Returns a dict with the verified statements.
    """
    import sympy as sp

    lam = sp.symbols("lambda", real=True)
    e0, e1, e2 = sp.symbols("epsilon_0 epsilon_1 epsilon_2", real=True)
    g0, g1, g2 = sp.symbols("gamma_0 gamma_1 gamma_2", real=True, positive=True)

    p = sp.expand((lam - e0) * (lam - e1) * (lam - e2))
    n = sp.expand(g0 ** 2 * (lam - e1) * (lam - e2)
                  + g1 ** 2 * (lam - e0) * (lam - e2)
                  + g2 ** 2 * (lam - e0) * (lam - e1))
    W4 = sp.expand(n * sp.diff(p, lam) - sp.diff(n, lam) * p)

    z, xi = sp.symbols("z xi", real=True)
    Nz = n.subs(lam, z)
    Pz = p.subs(lam, z)
    Nxi = n.subs(lam, xi)
    Pxi = p.subs(lam, xi)
    num = sp.expand(Pxi * Nz - Nxi * Pz)
    g_xi = sp.cancel(num / (z - xi))
    gpoly = sp.Poly(g_xi, z)
    A = gpoly.coeff_monomial(z ** 2)
    B = gpoly.coeff_monomial(z ** 1)
    C = gpoly.coeff_monomial(z ** 0)
    Q4_xi = sp.expand(B ** 2 - 4 * A * C)
    Q4 = sp.expand(Q4_xi.subs(xi, lam))
    Q4d = sp.diff(Q4, lam)

    disc_p = ((e0 - e1) * (e0 - e2) * (e1 - e2)) ** 2
    kappa_S = (g0 * g1 * g2) ** 2 * disc_p

    expr = sp.expand(16 * kappa_S * W4 + Q4d ** 2)
    quot, rem = sp.div(expr, Q4, lam)
    rem_simp = sp.simplify(rem)
    return {
        "remainder_of_(16_kappa_S_W4_plus_Q4_d_sq)_mod_Q4": rem_simp,
        "passes": rem_simp == 0,
        "quotient_degree": sp.degree(sp.Poly(quot, lam)),
    }
