"""
Empirical extraction of the inter-window 3x3 matrix D — Phase III shared
infrastructure.

The factorisation
:math:`U_{\\text{bench}} = N_B \\, D \\, N_A` is the central decomposition
of the Phase II / III closed-form push: :math:`N_X` is the elementary
2x2 Landau-Zener amplitude block on the active pair of window :math:`X`,
:math:`D \\in \\mathrm U(3)` is the inter-window connection matrix that
carries the genuinely transcendental content.

This module promotes the diagnostic code from ``notes/phase2_*.py`` to
first-class assay primitives. Two convention paths are exposed:

  * :func:`factor_D` uses the phase-free :func:`assay.closed_form.lz_rotation`
    blocks (matches the existing Phase II diagnostic).
  * :func:`factor_D_proper` (re-export from :mod:`assay.lz_stokes_phased`)
    uses Stokes-phased blocks (Phase III Route D).

The U(3) geodesic distance :func:`d_U3` is the gauge-invariant metric for
comparing two D matrices, used throughout the convergence harness.

Functions
---------
factor_D
    Phase-free factor :math:`D = N_B^T U N_A^T` with real-orthogonal :math:`N_X`.
extract_rotation
    Decompose D as :math:`\\mathrm{diag}(\\mathrm e^{\\mathrm i\\chi}) \\cdot R_{(sp_A, sp_B)}(\\theta)`.
empirical_D
    Run the IP propagator on a validation Record and return the empirical
    D plus its pair/spectator metadata. Caches per-record to
    ``graphify-out/empirical_D_cache.pkl`` (independent of graphify, just a
    convenient location).
d_U3
    U(3) geodesic distance, gauge-invariant.
"""

from __future__ import annotations

import os
import pickle
from typing import Optional

import numpy as np
import scipy.linalg as spla

from .closed_form import lz_rotation, embed, window_pairs
from .actions import window_action
from .geometry import Params, Geometry
from .ip import propagate_ad_ip
from .lz_stokes_phased import factor_D_proper  # noqa: F401  (re-export)


__all__ = [
    "factor_D",
    "factor_D_proper",
    "extract_rotation",
    "empirical_D",
    "d_U3",
]


# ---------------------------------------------------------------------------
#  Phase-free factor (matches notes/phase2_D_fast.py)
# ---------------------------------------------------------------------------
def factor_D(
    U: np.ndarray,
    p_A: float,
    p_B: float,
    pair_A: tuple,
    pair_B: tuple,
) -> np.ndarray:
    """
    Phase-free factor :math:`D = N_B^T U N_A^T` with real-orthogonal LZ blocks.

    Identical to ``notes/phase2_D_fast.py::factor_D`` — promoted here as the
    canonical phase-free entry point.
    """
    NA = embed(lz_rotation(p_A), pair_A)
    NB = embed(lz_rotation(p_B), pair_B)
    # N is real orthogonal: N^{-1} = N^T
    return NB.T @ U @ NA.T


# ---------------------------------------------------------------------------
#  Rotation extraction (matches notes/phase2_D_fast.py::extract_rotation)
# ---------------------------------------------------------------------------
def extract_rotation(D: np.ndarray, sp_A: int, sp_B: int) -> dict:
    """
    Decompose :math:`D = \\mathrm{diag}(\\mathrm e^{\\mathrm i\\chi})
    \\cdot R_{(sp_A, sp_B)}(\\theta)` and report the residual.

    Returns
    -------
    dict
        ``theta``, ``chi`` (length-3 array of phases), ``resid``
        (max-norm of :math:`|D - \\mathrm{diag}\\cdot R|`), plus the
        off-diagonal phase arguments for diagnostic use.
    """
    i_c = list({0, 1, 2} - {sp_A, sp_B})[0]
    chi_c = float(np.angle(D[i_c, i_c]))
    s_mag = 0.5 * (abs(D[sp_A, sp_B]) + abs(D[sp_B, sp_A]))
    c_mag = 0.5 * (abs(D[sp_A, sp_A]) + abs(D[sp_B, sp_B]))
    theta = float(np.arctan2(s_mag, c_mag))
    chi_sp_A = float(np.angle(D[sp_A, sp_A]))
    chi_sp_B = float(np.angle(D[sp_B, sp_B]))
    chi = np.zeros(3)
    chi[sp_A] = chi_sp_A
    chi[sp_B] = chi_sp_B
    chi[i_c] = chi_c

    R = np.eye(3, dtype=complex)
    R[sp_A, sp_A] = np.cos(theta)
    R[sp_B, sp_B] = np.cos(theta)
    R[sp_A, sp_B] = -np.sin(theta)
    R[sp_B, sp_A] = np.sin(theta)
    D_model = np.diag(np.exp(1j * chi)) @ R
    resid = float(np.max(np.abs(D - D_model)))

    arg_off_ab = float(np.angle(D[sp_A, sp_B])) if abs(D[sp_A, sp_B]) > 1e-12 else 0.0
    arg_off_ba = float(np.angle(D[sp_B, sp_A])) if abs(D[sp_B, sp_A]) > 1e-12 else 0.0
    return {
        "theta": theta,
        "chi": chi,
        "resid": resid,
        "arg_off_AB": arg_off_ab,
        "arg_off_BA": arg_off_ba,
    }


# ---------------------------------------------------------------------------
#  Per-record empirical D (cached)
# ---------------------------------------------------------------------------
_CACHE_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "assay",
    "empirical_D_cache.pkl",
)


def _load_cache() -> dict:
    if os.path.exists(_CACHE_PATH):
        try:
            with open(_CACHE_PATH, "rb") as f:
                return pickle.load(f)
        except Exception:
            return {}
    return {}


def _save_cache(cache: dict) -> None:
    os.makedirs(os.path.dirname(_CACHE_PATH), exist_ok=True)
    with open(_CACHE_PATH, "wb") as f:
        pickle.dump(cache, f)


def empirical_D(
    record,
    T: float = 120.0,
    rtol: float = 1e-12,
    atol: float = 1e-13,
    use_cache: bool = True,
    convention: str = "phase_free",
    sign_A: int = +1,
    sign_B: int = +1,
) -> dict:
    """
    Compute the empirical 3x3 inter-window matrix D from a validation Record.

    Parameters
    ----------
    record : assay.validation_matrix.Record
        Holds ``eps``, ``gam``, ``a``, ``x``.
    T : float
        Half-width of the integration interval (the IP propagator runs
        over ``[-T, +T]``). Default 120, matching the Phase II benchmark.
    rtol, atol : float
        IP propagator tolerances. Default tight (1e-12 / 1e-13).
    use_cache : bool
        If True, lookup/store in :file:`assay/empirical_D_cache.pkl`. The
        cache key includes the convention and tolerances so different
        runs don't collide.
    convention : {'phase_free', 'stokes_phased'}
        Which LZ-block convention to use. Phase-free matches Phase II;
        Stokes-phased is Phase III Route D.
    sign_A, sign_B : int
        Stokes-phase sign conventions (only used if
        ``convention == 'stokes_phased'``).

    Returns
    -------
    dict
        ``D`` (3x3 complex), ``pair_A``, ``pair_B``, ``sp_A``, ``sp_B``,
        ``p_A``, ``p_B``, ``delta_A``, ``delta_B``, ``I_A``, ``I_B``,
        ``U`` (the underlying IP fundamental matrix).
    """
    par = Params(
        eps=tuple(record.eps),
        gam=tuple(record.gam),
        a=tuple(record.a),
        x=int(record.x),
    )
    geo = Geometry(par)

    wp = window_pairs(geo)
    wp_sorted = sorted(wp, key=lambda d: d["u_center"])
    pair_A = wp_sorted[0]["pair"]
    pair_B = wp_sorted[1]["pair"]
    sp_A = wp_sorted[0]["spectator"]
    sp_B = wp_sorted[1]["spectator"]
    I_A = window_action(geo, wp_sorted[0]["window"])["I_X"]
    I_B = window_action(geo, wp_sorted[1]["window"])["I_X"]
    p_A = float(np.exp(-abs(I_A.imag)))
    p_B = float(np.exp(-abs(I_B.imag)))
    delta_A = float(abs(I_A.imag) / (2.0 * np.pi))
    delta_B = float(abs(I_B.imag) / (2.0 * np.pi))

    # cache key from (eps, gam, a, x, T, rtol, atol, convention, signs)
    cache = _load_cache() if use_cache else {}
    key = (
        tuple(par.eps), tuple(par.gam), tuple(par.a), par.x,
        float(T), float(rtol), float(atol),
        str(convention), int(sign_A), int(sign_B),
    )
    if use_cache and key in cache:
        out = dict(cache[key])
        # restore numpy arrays from possibly-pickled state
        out["D"] = np.asarray(out["D"])
        out["U"] = np.asarray(out["U"])
        return out

    U = propagate_ad_ip(geo, -T, T, par.x, rtol=rtol, atol=atol)

    if convention == "phase_free":
        D = factor_D(U, p_A, p_B, pair_A, pair_B)
    elif convention == "stokes_phased":
        D = factor_D_proper(
            U, p_A, p_B, delta_A, delta_B, pair_A, pair_B,
            sign_A=sign_A, sign_B=sign_B,
        )
    else:
        raise ValueError(f"unknown convention {convention!r}")

    out = {
        "D": D,
        "U": U,
        "pair_A": pair_A,
        "pair_B": pair_B,
        "sp_A": sp_A,
        "sp_B": sp_B,
        "p_A": p_A,
        "p_B": p_B,
        "delta_A": delta_A,
        "delta_B": delta_B,
        "I_A": complex(I_A),
        "I_B": complex(I_B),
    }
    if use_cache:
        cache[key] = out
        _save_cache(cache)
    return out


# ---------------------------------------------------------------------------
#  U(3) geodesic distance (gauge-invariant)
# ---------------------------------------------------------------------------
def d_U3(D1: np.ndarray, D2: np.ndarray) -> float:
    """
    Geodesic distance on :math:`\\mathrm U(3)`:

    .. math::

       d(D_1, D_2) \\;=\\; \\|\\,\\log(D_1^\\dagger \\, D_2)\\,\\|_F.

    This is gauge-invariant under right-multiplication by a global U(1)
    phase: :math:`d(D_1, \\mathrm e^{\\mathrm i\\alpha} D_2) = d(D_1, D_2)`
    is *false* in general — the U(1) phase factor enters as
    :math:`3\\alpha^2` under the trace. To make it strictly gauge-invariant
    in U(1) overall phase, the convention here projects out the overall
    phase by setting :math:`\\det D_1^\\dagger D_2` to unit modulus before
    taking ``logm``.
    """
    M = D1.conj().T @ D2
    # remove overall U(1) phase: divide by det(M)^(1/3) so det = 1
    # (det is a complex number with |det| = 1 since both are unitary).
    detM = np.linalg.det(M)
    M = M / (detM ** (1.0 / 3.0))
    L = spla.logm(M)
    return float(np.linalg.norm(L, "fro"))
