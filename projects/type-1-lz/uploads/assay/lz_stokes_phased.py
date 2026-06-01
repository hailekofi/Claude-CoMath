"""
Stokes-phased Landau-Zener amplitude block — Phase III Route D.

The phase-free LZ block in :mod:`assay.closed_form` (``lz_rotation(p)``) is
the real-orthogonal :math:`R(\\theta)` with :math:`\\sin^2\\theta = p`.
The proper two-level Landau-Zener S-matrix in the labelled-adiabatic basis
carries Stokes phases :math:`\\varphi_S(\\delta)` on its diagonal entries:

.. math::

   N_{\\text{proper}}(p,\\delta) \\;=\\;
   \\begin{pmatrix}
     \\sqrt{1-p}\\,\\mathrm e^{-\\mathrm i\\varphi_S(\\delta)} & \\sqrt p \\\\
     -\\sqrt p & \\sqrt{1-p}\\,\\mathrm e^{+\\mathrm i\\varphi_S(\\delta)}
   \\end{pmatrix},

with

.. math::

   p \\;=\\; \\mathrm e^{-2\\pi\\delta}, \\qquad
   \\varphi_S(\\delta) \\;=\\; \\delta(\\log\\delta - 1) + \\frac{\\pi}{4}
                                + \\arg\\Gamma(1 - \\mathrm i\\delta).

This is the standard Joye-Kunz-Pfister convention (1991) used in the
multistate Landau-Zener literature. The Stokes phase
:math:`\\varphi_S(\\delta)` is the unique transcendental ingredient of the
elementary 2-level LZ amplitude, carrying everything beyond the
real-orthogonal rotation magnitude :math:`\\sqrt{1-p}, \\sqrt p`.

When the phase-free convention is used to factor
:math:`U_{\\text{bench}} = N_B \\cdot D \\cdot N_A`, the extracted D
absorbs the missing Stokes phases. The conjecture under test in Phase III
Route D is that :math:`D_{\\text{corrected}}` — extracted with the proper
:math:`N_X` containing :math:`\\varphi_S` — has a simpler regime-independent
structure than the phase-free :math:`D`.

The ``sign_A``, ``sign_B`` flags expose the convention's two sign
ambiguities (sign of :math:`\\varphi_S`, which off-diagonal carries the
minus sign) for the single-window calibration documented in the Phase III
plan.
"""

from __future__ import annotations

import numpy as np
from scipy.special import loggamma


__all__ = [
    "stokes_phase",
    "lz_block_proper",
    "embed_complex",
    "factor_D_proper",
]


def stokes_phase(delta: float) -> float:
    """
    LZ Stokes phase :math:`\\varphi_S(\\delta)`.

    Parameters
    ----------
    delta : float
        Landau-Zener adiabaticity parameter, :math:`\\delta = |\\im I_X|/(2\\pi)`.

    Returns
    -------
    float
        The Stokes phase (radians), defined for all :math:`\\delta \\ge 0`.

    Notes
    -----
    The formula has three terms:
      * :math:`\\delta(\\log\\delta - 1)` — the leading WKB phase
        (well-defined for :math:`\\delta > 0`; limit 0 as
        :math:`\\delta \\to 0^+`).
      * :math:`\\pi/4` — the standard turning-point Stokes contribution.
      * :math:`\\arg\\Gamma(1 - \\mathrm i\\delta)` — the exact non-perturbative
        correction (Berry's "matching" data).

    For :math:`\\delta \\to 0`: :math:`\\varphi_S \\to \\pi/4`.
    For :math:`\\delta \\to \\infty`: :math:`\\varphi_S \\sim \\delta\\log\\delta`
    (grows; the phase-free convention is recovered up to a globally
    accumulating phase that can be gauged away in :math:`P = |U|^2`).
    """
    d = float(delta)
    if d < 0:
        raise ValueError(f"delta must be >= 0, got {d!r}")
    # delta * (log delta - 1) handled robustly at delta = 0
    if d < 1e-15:
        wkb = 0.0
    else:
        wkb = d * (np.log(d) - 1.0)
    # arg Gamma(1 - i delta) via loggamma (complex)
    lg = loggamma(1.0 - 1j * d)
    arg_gamma = float(np.imag(lg))
    return wkb + 0.25 * np.pi + arg_gamma


def lz_block_proper(p: float, delta: float, sign: int = +1) -> np.ndarray:
    """
    Proper 2x2 LZ amplitude block with Stokes phase.

    Parameters
    ----------
    p : float
        Elementary jump probability, :math:`p = \\exp(-2\\pi\\delta)`.
    delta : float
        Landau-Zener adiabaticity (consistent with ``p`` via
        ``delta = -log(p) / (2 pi)``).
    sign : int
        Convention flag ``+1`` (default) or ``-1`` — selects the sign of
        :math:`\\varphi_S` in the diagonal entries. The "true" convention
        is pinned by the single-window calibration; both are equally valid
        up to an overall basis choice. See module docstring.

    Returns
    -------
    np.ndarray, shape (2, 2), complex
        ::

           N = [[ sqrt(1-p) * exp(-i sign * phi_S),   sqrt(p)             ],
                [-sqrt(p),                            sqrt(1-p) * exp(+i sign * phi_S) ]]

    Notes
    -----
    The phase-free limit (:math:`\\varphi_S \\to 0`) recovers
    :func:`assay.closed_form.lz_rotation` modulo a sign convention on the
    off-diagonal (here ``[0,1] = +sqrt(p)``, ``[1,0] = -sqrt(p)``; the
    phase-free version uses the opposite signs).
    """
    p = max(0.0, min(1.0, float(p)))
    c = float(np.sqrt(1.0 - p))
    s = float(np.sqrt(p))
    phi = sign * stokes_phase(float(delta))
    return np.array(
        [
            [c * np.exp(-1j * phi),  s + 0.0j],
            [-s + 0.0j,              c * np.exp(+1j * phi)],
        ],
        dtype=complex,
    )


def embed_complex(M2: np.ndarray, pair: tuple) -> np.ndarray:
    """
    Embed a complex 2x2 block onto 3x3 acting on ``pair``, identity on the
    spectator level.

    Parameters
    ----------
    M2 : np.ndarray, shape (2, 2), complex
        The 2x2 amplitude block to embed.
    pair : tuple of int, length 2
        Indices ``(i, j)`` with :math:`i, j \\in \\{0, 1, 2\\}`,
        :math:`i \\ne j`.

    Returns
    -------
    np.ndarray, shape (3, 3), complex
        Embedded 3x3 matrix with ``M2`` on the ``pair`` block and 1.0 on
        the spectator diagonal.

    Notes
    -----
    This is the complex-valued analogue of :func:`assay.closed_form.embed`,
    which silently casts to real. Routes that use ``lz_block_proper`` must
    use ``embed_complex`` to preserve the Stokes phase.
    """
    i, j = pair
    sp = ({0, 1, 2} - {i, j}).pop()
    M = np.zeros((3, 3), dtype=complex)
    M[i, i], M[i, j] = M2[0, 0], M2[0, 1]
    M[j, i], M[j, j] = M2[1, 0], M2[1, 1]
    M[sp, sp] = 1.0 + 0.0j
    return M


def factor_D_proper(
    U: np.ndarray,
    p_A: float,
    p_B: float,
    delta_A: float,
    delta_B: float,
    pair_A: tuple,
    pair_B: tuple,
    sign_A: int = +1,
    sign_B: int = +1,
) -> np.ndarray:
    """
    Solve D from :math:`U = N_B^{\\text{proper}} \\, D \\, N_A^{\\text{proper}}`.

    Parameters
    ----------
    U : np.ndarray, shape (3, 3), complex
        Empirical lab-frame fundamental matrix from
        :func:`assay.ip.propagate_ad_ip` on the validation sample.
    p_A, p_B : float
        Elementary jump probabilities at windows A and B.
    delta_A, delta_B : float
        Adiabaticities at windows A and B (consistent with ``p_X``).
    pair_A, pair_B : tuple of int
        Active level pairs (length 2) at the two windows.
    sign_A, sign_B : int
        Stokes-phase sign conventions; either ``+1`` (default) or ``-1``.
        Set by the single-window calibration in the Phase III plan.

    Returns
    -------
    np.ndarray, shape (3, 3), complex
        ``D = N_B^(-1) @ U @ N_A^(-1)``, a unitary 3x3 matrix.

    Notes
    -----
    Unlike :func:`assay.closed_form.lz_rotation` (which is real-orthogonal
    so :math:`N^{-1} = N^T`), ``N_proper`` is complex unitary, so we use
    the conjugate-transpose inverse :math:`N^{-1} = N^\\dagger`.
    """
    NA = embed_complex(lz_block_proper(p_A, delta_A, sign=sign_A), pair_A)
    NB = embed_complex(lz_block_proper(p_B, delta_B, sign=sign_B), pair_B)
    return NB.conj().T @ U @ NA.conj().T
