"""
Mixed-cover periods sensitive to both phase and transport sectors.

The phase one-form  :math:`\omega = \sqrt{Q_4}\,L_H\,du/p`  lives on
:math:`\mu^2 = Q_4`, and the transport one-form  :math:`\eta` constructed in
:mod:`assay.transport_actions` lives on :math:`\nu^2 = W_4`.  Their natural
common cover is  :math:`\Sigma_{\Delta,S}`, the fibre product of the two
elliptic curves; it has genus 5.

Several first-order candidates that sample BOTH branches simultaneously:

  M1  =  omega * eta / du
       =  - sqrt(kappa_S) * L_H * sqrt(W_4) / ( p * sqrt(Q_4) ) * dlam
       = - sqrt(kappa_S) * L_H / sqrt(Q_4) * (sqrt(W_4)/p) * dlam.

  M2  =  omega + eta
  M3  =  omega - eta
  M4  =  (omega^2 + eta^2)        -- on  Sigma_{Delta,S}, double cover
  M5  =  omega * eta
  Mlog =  d log(rho)              -- transport gauge differential

Periods are computed around each of the two projected Q4 conjugate-pair
windows.  Because W_4 branch points sit elsewhere, both sheets of the
spinor cover are picked up cleanly inside a Q4-window contour by branch
tracking sqrt(W_4) AND sqrt(Q_4) together.

This module also exposes :func:`compute_period_battery`, which evaluates
every candidate period on every dataset record (parallelised), and
:func:`make_formula_candidate`, which packages a period combination into
a 3x3 closed-form callable for use with :mod:`assay.validation_matrix`.
"""

from __future__ import annotations
import numpy as np
from typing import Callable, Dict, List, Optional, Tuple

from .geometry import Geometry, Params
from .actions import window_action
from .transport_actions import transport_period_at_Q4_pair


# ---------------------------------------------------------------------------
#  Contour quadrature with branch-tracked sqrt(Q4) AND sqrt(W4)
# ---------------------------------------------------------------------------
def _contour(win: dict, n: int, radius_scale: float = 1.6):
    q = win["roots"]
    centre = 0.5 * (q[0] + q[1])
    radius = radius_scale * abs(q[0] - q[1])
    th = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    lam = centre + radius * np.exp(1j * th)
    dlam = 1j * radius * np.exp(1j * th) * (2.0 * np.pi / n)
    return lam, dlam, centre, radius


def _branch_track(values):
    """Continuously track sqrt of a complex array around a closed contour."""
    n = len(values)
    out = np.empty(n, complex)
    s = np.sqrt(values[0])
    for k in range(n):
        cand = np.sqrt(values[k])
        s = cand if abs(cand - s) < abs(-cand - s) else -cand
        out[k] = s
    return out


def mixed_periods_one_window(geo: Geometry, win: dict, n: int = 8000,
                             radius_scale: float = 1.6) -> dict:
    """Compute every candidate period around a single Q4 window."""
    lam, dlam, centre, radius = _contour(win, n, radius_scale)

    Q4v = np.polyval(geo.Q4, lam)
    W4v = np.polyval(geo.W4, lam)
    sQ = _branch_track(Q4v)
    sW = _branch_track(W4v)

    LHv = np.polyval(geo.LH, lam)
    Pv = np.polyval(geo.p, lam)
    nuS = np.sqrt(geo.kappaS)

    # omega / dlam = -sqrt(Q4) L_H W4 / p^3   (from actions.py)
    omega = -sQ * LHv * W4v / Pv ** 3
    # eta / dlam = -sqrt(kS) sqrt(W4) / Q4    (from transport_actions.py)
    eta = -nuS * sW / Q4v
    # du / dlam = -W4 / p^2,  inv = -p^2 / W4
    inv_du = -Pv ** 2 / W4v

    # M1 = omega * eta / du = -sqrt(kS) * L_H * sqrt(W4) / (p * sqrt(Q4))
    M1 = -nuS * LHv * sW / (Pv * sQ)
    # M2, M3: omega +- eta
    M2 = omega + eta
    M3 = omega - eta
    # M4: omega^2 + eta^2  (second-order)
    M4 = omega ** 2 + eta ** 2
    # M5: omega * eta (single contour winding around eta-curve + omega-curve)
    M5 = omega * eta
    # Mlog: d log rho/d lam = d log(sqrt(kS) L_H^2 / sqrt(W4))/d lam
    #     = 2 L_H' / L_H - (1/2) W4' / W4
    LH_d = np.array([geo.LH[0]])  # derivative of linear L_H is alpha_H = -LH[0]
    W4_d = np.polyder(geo.W4)
    # Use derivatives of polynomials directly
    LH_der_val = geo.LH[0]
    Mlog = 2.0 * LH_der_val / LHv - 0.5 * np.polyval(W4_d, lam) / W4v
    # M6: (omega + eta)^2 / du
    M6 = (omega + eta) ** 2 * inv_du
    # M7: omega * eta * du (volume element-like)
    M7 = omega * eta * (-W4v / Pv ** 2)  # = omega * eta * du/dlam, but du/dlam already
    # Wait: integrand in lambda is (form)/dlam * dlam. omega is omega/dlam, so omega*dlam is the form.
    # omega*eta*du means multiplied by (du/dlam) -- actually the natural product is omega*eta as a
    # 2-form, but on the cover that's a 1-form once you fix variable. Skip M7 in main test.

    out = {
        "centre": complex(centre),
        "radius": float(radius),
        "I_omega": complex(np.sum(omega * dlam)),
        "I_eta": complex(np.sum(eta * dlam)),
        "I_M1": complex(np.sum(M1 * dlam)),
        "I_M2": complex(np.sum(M2 * dlam)),
        "I_M3": complex(np.sum(M3 * dlam)),
        "I_M4": complex(np.sum(M4 * dlam)),
        "I_M5": complex(np.sum(M5 * dlam)),
        "I_M6": complex(np.sum(M6 * dlam)),
        "I_Mlog": complex(np.sum(Mlog * dlam)),
    }
    return out


def all_mixed_periods(geo: Geometry, n: int = 8000,
                      radius_scale: float = 1.6) -> Tuple[dict, dict]:
    """Mixed-cover periods on both Q4-windows in u-order."""
    wins = sorted(geo.windows(), key=lambda w: w["u_center"])
    return tuple(mixed_periods_one_window(geo, w, n=n, radius_scale=radius_scale)
                 for w in wins)


# ---------------------------------------------------------------------------
#  Candidate magic-number library
# ---------------------------------------------------------------------------
def candidate_magic_numbers(periods_pair) -> Dict[str, float]:
    """
    Given the pair of (window-A, window-B) period dicts, build a battery of
    candidate "magic numbers" in [0, 1] that COULD encode a transition
    probability or survival probability.

    Returns a dict {name -> magic number in [0,1] (or possibly outside)}.
    """
    A, B = periods_pair
    out = {}
    # |Im| of each first-order period, both windows
    keys = ["I_omega", "I_eta", "I_M1", "I_M2", "I_M3", "I_M4", "I_M5",
            "I_M6", "I_Mlog"]
    for k in keys:
        out[f"exp_-|Im_{k}_A|"] = float(np.exp(-abs(A[k].imag)))
        out[f"exp_-|Im_{k}_B|"] = float(np.exp(-abs(B[k].imag)))
        out[f"exp_-|Re_{k}_A|"] = float(np.exp(-abs(A[k].real)))
        out[f"exp_-|Re_{k}_B|"] = float(np.exp(-abs(B[k].real)))
    # combinations
    for k in keys:
        sA = abs(A[k].imag)
        sB = abs(B[k].imag)
        out[f"exp_-(|Im_{k}_A|+|Im_{k}_B|)"] = float(np.exp(-(sA + sB)))
        out[f"exp_-|Im_{k}_A|-|Im_{k}_B|"] = out[f"exp_-(|Im_{k}_A|+|Im_{k}_B|)"]
        out[f"exp_-||Im_{k}_A|-|Im_{k}_B||"] = float(np.exp(-abs(sA - sB)))
        if sA > 0 or sB > 0:
            out[f"sin2(|Im_{k}_A|-|Im_{k}_B|/2)"] = float(np.sin(abs(sA - sB) / 2.0) ** 2)
    return out


# ---------------------------------------------------------------------------
#  Higher-level formula factory for use with assay.validation_matrix.validate
# ---------------------------------------------------------------------------
def be_with_mixed_middle(rec) -> np.ndarray:
    """
    Trial formula: BE for the two extremes; predicted middle survival from
    a transport-side or mixed-period candidate.

    Default: use survival = 1 - p_A - p_B + p_A p_B  (incoherent middle,
    classical reflection probability of the middle level off both
    crossings).  This is identical to assay.grid (more or less) so it's
    a baseline.

    Replaced by candidate combinations in validate_candidates().
    """
    a = np.array(rec.a)
    pi_in = np.argsort(a)
    pi_out = np.argsort(-a)
    n_max = int(np.argmax(a))
    n_min = int(np.argmin(a))
    n_mid = ({0, 1, 2} - {n_max, n_min}).pop()
    N = 3
    P = np.zeros((N, N))
    # BE extremes
    for n in (n_max, n_min):
        q_prod = 1.0
        for m in range(N):
            if m == n:
                continue
            pair = (n, m) if n < m else (m, n)
            q_prod *= rec.q_ij[pair]
        i_in = int(np.where(pi_in == n)[0][0])
        i_out = int(np.where(pi_out == n)[0][0])
        P[i_in, i_out] = q_prod
    # Middle survival from product (1 - pA)(1 - pB)
    par = Params(eps=rec.eps, gam=rec.gam, a=rec.a, x=rec.x)
    geo = Geometry(par)
    wins = sorted(geo.windows(), key=lambda w: w["u_center"])
    p_X = []
    for w in wins:
        I_X = window_action(geo, w)["I_X"]
        p_X.append(float(np.exp(-abs(I_X.imag))))
    pA, pB = p_X
    P_mid_survival = (1.0 - pA) * (1.0 - pB)
    i_in_mid = int(np.where(pi_in == n_mid)[0][0])
    i_out_mid = int(np.where(pi_out == n_mid)[0][0])
    P[i_in_mid, i_out_mid] = P_mid_survival
    return P


# ---------------------------------------------------------------------------
#  Stokes-data fallback: ordinary LZ Stokes phase at each Q4 turning point
# ---------------------------------------------------------------------------
def lz_stokes_phase(delta: float) -> float:
    """
    Standard LZ Stokes phase.

      phi_S(delta) = delta * (log(delta) - 1) + (pi/4) + arg Gamma(1 - i*delta).

    Here ``delta`` is the local adiabaticity parameter (Massey's 'delta',
    real).  When delta -> 0 the phase tends to pi/4.
    """
    from scipy.special import loggamma
    d = float(delta)
    if abs(d) < 1e-14:
        return float(np.pi / 4)
    val = d * (np.log(d) - 1.0) + np.pi / 4
    val += float(np.imag(loggamma(complex(1.0, -d))))
    return val


def adiabaticity_at_window(geo: Geometry, win: dict) -> float:
    """
    Local LZ adiabaticity parameter delta = |I_X| / (2 pi) at a Q4 window.
    (This is the Dykhne exponent; the Stokes phase uses the same delta.)
    """
    I_X = window_action(geo, win)["I_X"]
    return float(abs(I_X.imag) / (2.0 * np.pi))


# ---------------------------------------------------------------------------
#  Battery: compute candidate periods on every record
# ---------------------------------------------------------------------------
def compute_period_battery(records) -> List[Dict[str, complex]]:
    """For each record, return a flattened dict of all mixed periods (both windows)."""
    out = []
    for rec in records:
        par = Params(eps=rec.eps, gam=rec.gam, a=rec.a, x=rec.x)
        geo = Geometry(par)
        try:
            A, B = all_mixed_periods(geo)
        except Exception:
            out.append({})
            continue
        flat = {}
        for k, v in A.items():
            flat[f"{k}_A"] = v
        for k, v in B.items():
            flat[f"{k}_B"] = v
        # Per-window absolute imaginary parts for convenience
        flat["abs_Im_omega_A"] = abs(complex(A["I_omega"]).imag)
        flat["abs_Im_omega_B"] = abs(complex(B["I_omega"]).imag)
        out.append(flat)
    return out
