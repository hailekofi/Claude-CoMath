"""
oracle.py  --  WS-F gold-standard P(gamma, eps, a) oracle for the Type-1, N=3 MLZ problem.

This module is the *trusted* high-accuracy oracle for the full 3x3 doubly-stochastic
transition matrix P of the Type-1 Cauchy multistate Landau-Zener model

      H(u) = H0 + u * diag(a),     i dpsi/du = H(u) psi,   u in (-inf, +inf)
      (H0)_ij = gamma_i gamma_j (a_i - a_j)/(eps_i - eps_j)        (i != j)
      (H0)_ii = - sum_{k!=i} gamma_k^2 (a_i - a_k)/(eps_i - eps_k)

It is a clean wrapper over the project's interaction-picture harness
``uploads/assay/ip.py`` (the adiabatic-IP propagator ``propagate_ad_ip``), with the
two fixes that the harness needed to reach gold-standard accuracy faithfully:

  FIX 1 (convention).  The diabatic <-> adiabatic endpoint permutation used throughout
  the assay (``pi_in = argsort(a)``, ``pi_out = argsort(-a)``) is WRONG.  The IP
  propagator labels its three channels by the *spectral sheet* (lam_0 exterior to the
  poles, lam_1 in (eps0,eps1), lam_2 in (eps1,eps2)).  Because the spectral roots
  interlace the real poles for all real u, the sheet <-> diabatic-channel map is a
  FIXED permutation independent of a:
        incoming (u -> -inf):  sheet i  ==  diabatic channel i      PI_IN  = (0,1,2)
        outgoing (u -> +inf):  sheet 0->dia 2, 1->dia 0, 2->dia 1   PI_OUT = (2,0,1)
  (Verified against the BE extreme survivals landing exactly on the diabatic diagonal,
  and against the asymptotic energy slopes E_i(u)/u -> a_{map(i)} for >190 random
  samples; see oracle_report.md.)

  FIX 2 (extrapolation order).  The probability truncation error of the adiabatic-IP
  propagator decays as C / T^4 (the oscillatory 1/T amplitude tails cancel in the
  doubly-stochastic probability).  ``validation_matrix.benchmark_P`` Richardson-
  extrapolates with the 8:1 (1/T^3) weight (8 P(2T) - P(T))/7, which is mismatched and
  *under*-corrects.  The correct combination is the 16:1 (1/T^4) weight
        P_extrap = (16 P(2T) - P(T)) / 15 .
  (Verified: successive halving of the truncation error by ~16x; see report.)

The oracle returns the full 3x3 P in the **diabatic** basis (rows = incoming diabatic
channel x, columns = outgoing diabatic channel j, P[x,j] = prob x -> j), together with a
conservative per-entry error estimate.

CONVENTIONS (see oracle_report.md for the full pinned-down statement)
---------------------------------------------------------------------
  * Levels are indexed 0,1,2 in the SAME order as eps (strictly increasing eps).
  * "Slope order" lo/mid/hi = argsort(a): lo = arg min a, hi = arg max a, mid = middle.
  * The two EXACT Brundobler-Elser survivals are the extreme-SLOPE diagonal entries
        P[lo,lo] = prod_{m!=lo} exp(-2 pi Gamma_{lo,m}),
        P[hi,hi] = prod_{m!=hi} exp(-2 pi Gamma_{hi,m}),
        Gamma_ij = gamma_i^2 gamma_j^2 |a_i - a_j| / (eps_i - eps_j)^2 .
  * The OPEN middle survival is  P[mid,mid]  (mid = the middle-slope diabatic level).

Requires: numpy, scipy, and the assay package importable from ``../uploads``.
"""
from __future__ import annotations

import os
import sys
import numpy as np

# --- make the assay package importable --------------------------------------
_HERE = os.path.dirname(os.path.abspath(__file__))
_UPLOADS = os.path.join(_HERE, "..", "uploads")
if _UPLOADS not in sys.path:
    sys.path.insert(0, _UPLOADS)

from assay import Params, Geometry                      # noqa: E402
from assay.ip import propagate_ad_ip                    # noqa: E402
from assay.adiabatic import adiabatic_data              # noqa: E402


# ===========================================================================
#  Fixed spectral-sheet <-> diabatic-channel permutations (FIX 1)
# ===========================================================================
PI_IN = (0, 1, 2)     # diabatic channel reached by spectral sheet i at u -> -inf
PI_OUT = (2, 0, 1)    # diabatic channel reached by spectral sheet i at u -> +inf


def slope_permutation(geo: Geometry, u: float) -> tuple:
    """
    Empirical sheet -> diabatic map at base point u, from the asymptotic energy
    slope E_i(u)/u matched to the nearest a_k.  Used only to VERIFY PI_IN/PI_OUT
    (and to detect a node where interlacing fails).  Returns a length-3 tuple.
    """
    a = np.asarray(geo.a, float)
    E = adiabatic_data(geo, u)["E"].real
    return tuple(int(np.argmin(np.abs(a - E[i] / u))) for i in range(3))


# ===========================================================================
#  Type-1 algebra helpers
# ===========================================================================
def Gamma_ij(eps, gam, a, i, j) -> float:
    """Pairwise Brundobler-Elser exponent Gamma_ij = g_i^2 g_j^2 |a_i-a_j|/(eps_i-eps_j)^2."""
    return gam[i] ** 2 * gam[j] ** 2 * abs(a[i] - a[j]) / (eps[i] - eps[j]) ** 2


def be_survivals(eps, gam, a) -> dict:
    """
    Exact BE extreme-slope survivals and the incoherent middle-survival baseline.
    Returns dict with keys lo, mid, hi (slope indices), P_lo, P_hi (exact), P_mid_inc.
    """
    eps = np.asarray(eps, float); gam = np.asarray(gam, float); a = np.asarray(a, float)
    lo, mid, hi = (int(k) for k in np.argsort(a))
    P_lo = np.exp(-2 * np.pi * (Gamma_ij(eps, gam, a, lo, mid) + Gamma_ij(eps, gam, a, lo, hi)))
    P_hi = np.exp(-2 * np.pi * (Gamma_ij(eps, gam, a, hi, mid) + Gamma_ij(eps, gam, a, hi, lo)))
    P_mid_inc = np.exp(-2 * np.pi * (Gamma_ij(eps, gam, a, mid, lo) + Gamma_ij(eps, gam, a, mid, hi)))
    return dict(lo=lo, mid=mid, hi=hi, P_lo=float(P_lo), P_hi=float(P_hi),
                P_mid_inc=float(P_mid_inc))


# ===========================================================================
#  Core: diabatic P at a single truncation T (adiabatic-IP engine)
# ===========================================================================
def _P_diabatic_at_T(geo: Geometry, T: float, rtol: float, atol: float) -> np.ndarray:
    """
    Diabatic transition matrix P[x,j] = prob(x -> j) at truncation half-window T,
    from the adiabatic-IP propagator, reindexed by the fixed PI_IN/PI_OUT maps.

    propagate_ad_ip(geo,-T,T,x) returns the lab-frame fundamental column for incoming
    spectral sheet x; |U[j,x]|^2 over j is the adiabatic transition x->j.  We then
    place adiabatic (i,j) into diabatic (PI_IN[i], PI_OUT[j]).
    """
    Pa = np.empty((3, 3))
    for x in range(3):
        U = propagate_ad_ip(geo, -T, T, x, rtol=rtol, atol=atol)
        Pa[x, :] = np.abs(U[:, x]) ** 2          # Pa[x,j] = adiabatic x->j
    Pd = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            Pd[PI_IN[i], PI_OUT[j]] = Pa[i, j]
    return Pd


# ===========================================================================
#  Public oracle
# ===========================================================================
def oracle_P(eps, gam, a, T: float = 80.0,
             rtol: float = 1e-13, atol: float = 1e-14,
             richardson: bool = True,
             verify_convention: bool = True) -> dict:
    """
    Gold-standard diabatic transition matrix P(gamma, eps, a).

    Parameters
    ----------
    eps, gam, a : length-3 model parameters (eps strictly increasing, gam nonzero).
    T           : base truncation half-window.  The propagator is run at T and 2T and
                  16:1 Richardson-extrapolated in the cutoff (tail ~ 1/T^4).
    rtol, atol  : DOP853 tolerances for the IP IVP.
    richardson  : if False, return the raw P(2T) without extrapolation.
    verify_convention : if True, check the fixed PI_IN/PI_OUT against the asymptotic
                  energy-slope map and flag any mismatch (e.g. near a node).

    Returns
    -------
    dict with:
      P        : (3,3) diabatic transition matrix, P[x,j] = prob x->j  (gold standard).
      err      : (3,3) conservative per-entry error estimate.
      P_T,P_2T : the raw matrices feeding the Richardson combination.
      slope_order : dict(lo,mid,hi,...) from be_survivals (slope indices + exact BE).
      be_delta : dict of |P_oracle - P_exact| for the two BE extreme survivals.
      doubly_stochastic_defect : max |rowsum-1|, |colsum-1|.
      convention_ok : bool (PI_IN/PI_OUT verified) ; convention_detail on mismatch.
      meta     : dict of run parameters.
    """
    eps = tuple(float(v) for v in eps)
    gam = tuple(float(v) for v in gam)
    a = tuple(float(v) for v in a)
    geo = Geometry(Params(eps=eps, gam=gam, a=a, x=0))

    conv_ok, conv_detail = True, {}
    if verify_convention:
        pin = slope_permutation(geo, -1.0e5)
        pout = slope_permutation(geo, 1.0e5)
        conv_ok = (pin == PI_IN and pout == PI_OUT)
        conv_detail = dict(pin_measured=pin, pout_measured=pout,
                           PI_IN=PI_IN, PI_OUT=PI_OUT)

    P_T = _P_diabatic_at_T(geo, T, rtol, atol)
    P_2T = _P_diabatic_at_T(geo, 2.0 * T, rtol, atol)

    if richardson:
        P = (16.0 * P_2T - P_T) / 15.0
        # conservative per-entry error: max of (Richardson correction residual) and
        # (distance from the finer raw grid).  The true error is ~ |P-P_2T|/15-ish;
        # we report the larger, safe surrogate.
        err = np.maximum(np.abs(P - P_2T), np.abs(P_2T - P_T) / 15.0)
    else:
        P = P_2T.copy()
        err = np.abs(P_2T - P_T)

    bes = be_survivals(eps, gam, a)
    lo, hi = bes["lo"], bes["hi"]
    be_delta = {"lo": float(abs(P[lo, lo] - bes["P_lo"])),
                "hi": float(abs(P[hi, hi] - bes["P_hi"]))}

    ds_defect = float(max(np.max(np.abs(P.sum(axis=1) - 1.0)),
                          np.max(np.abs(P.sum(axis=0) - 1.0))))

    return dict(P=P, err=err, P_T=P_T, P_2T=P_2T,
                slope_order=bes, be_delta=be_delta,
                doubly_stochastic_defect=ds_defect,
                convention_ok=conv_ok, convention_detail=conv_detail,
                meta=dict(eps=eps, gam=gam, a=a, T=T, rtol=rtol, atol=atol,
                          richardson=richardson))


# ===========================================================================
#  Stratified benchmark suite
# ===========================================================================
#  Each stratum gives (eps, gam, a) and a short description.  The near-node and
#  near-degenerate-slope entries are documented in oracle_report.md.
SUITE = {
    # canonical reference (research program)
    "canonical": dict(
        eps=(-2.0, 0.0, 3.0), gam=(1.0, 0.8, 1.2), a=(-1.0, 0.5, 2.0),
        desc="canonical reference (well-separated, moderate overlap)"),
    # strong-interference overlapping case (the 100x middle-survival enhancement)
    "sampleB": dict(
        eps=(-1.0, 0.0, 1.5), gam=(0.9, 1.1, 0.8), a=(-0.7, 0.4, 1.3),
        desc="overlapping crossings; strong interference, ~100x P_mid enhancement"),
    # well-separated crossings (large eps gaps, modest slopes)
    "well_separated": dict(
        eps=(-5.0, 0.0, 5.0), gam=(1.0, 1.0, 1.0), a=(-1.5, 0.0, 1.5),
        desc="well-separated crossings (large eps gaps)"),
    # weak coupling (small gam): nearly-diabatic, P ~ I
    "weak_coupling": dict(
        eps=(-2.0, 0.0, 3.0), gam=(0.35, 0.30, 0.40), a=(-1.0, 0.5, 2.0),
        desc="weak coupling (small gamma): near-diabatic passage"),
    # strong coupling (large gam): nearly-adiabatic
    "strong_coupling": dict(
        eps=(-2.0, 0.0, 3.0), gam=(2.2, 2.0, 2.4), a=(-1.0, 0.5, 2.0),
        desc="strong coupling (large gamma): near-adiabatic passage"),
    # near-degenerate slope (two slopes close => one crossing very slow/adiabatic)
    "near_deg_slope": dict(
        eps=(-2.0, 0.0, 3.0), gam=(1.0, 0.8, 1.2), a=(-1.0, 0.45, 0.55),
        desc="near-degenerate slopes (a_mid ~ a_hi): one near-adiabatic crossing"),
    # near-node: smallest achievable real-axis eigenvalue gap (avoided crossing
    # closest to the spectral-curve node).  Parameters set in build_suite().
    # (filled at import time below)
}


def _near_node_params():
    """
    Near-node stratum: a bounded parameter set whose minimum real-axis eigenvalue
    gap is the smallest we could find (the avoided crossing closest to the node of
    the spectral curve).  For generic real Type-1 data the three spectral roots
    strictly interlace the real poles, so an EXACT real-u degeneracy does not occur;
    this is the closest physical approach.  See oracle_report.md.
    """
    return dict(eps=(-2.0, -1.4, 3.0), gam=(1.5, 1.5, 0.6), a=(-2.0, 1.9, 2.1),
                desc="near-node: small real-axis avoided-crossing gap (node approach)")


SUITE["near_node"] = _near_node_params()


def run_suite(T: float = 80.0, rtol: float = 1e-13, atol: float = 1e-14,
              names=None, verbose: bool = True) -> dict:
    """
    Run the gold oracle on every stratum (or the subset ``names``).  Returns
    {name: oracle_P(...) result dict}.  Prints a compact table if verbose.
    """
    out = {}
    names = names or list(SUITE.keys())
    for name in names:
        spec = SUITE[name]
        res = oracle_P(spec["eps"], spec["gam"], spec["a"], T=T, rtol=rtol, atol=atol)
        res["desc"] = spec["desc"]
        out[name] = res
        if verbose:
            so = res["slope_order"]
            mid = so["mid"]
            print(f"\n[{name}]  {spec['desc']}")
            print(f"  eps={spec['eps']} gam={spec['gam']} a={spec['a']}")
            with np.printoptions(precision=10, suppress=True, linewidth=120):
                print("  P (diabatic, P[x,j]=x->j) =")
                print("   ", str(res["P"]).replace("\n", "\n    "))
            print(f"  max per-entry err est = {np.max(res['err']):.2e}")
            print(f"  doubly-stochastic defect = {res['doubly_stochastic_defect']:.2e}")
            print(f"  BE deltas: lo={res['be_delta']['lo']:.2e}  hi={res['be_delta']['hi']:.2e}")
            print(f"  P_mid (open) = {res['P'][mid,mid]:.10f}   "
                  f"incoherent = {so['P_mid_inc']:.10f}   "
                  f"ratio = {res['P'][mid,mid]/so['P_mid_inc']:.3f}")
            if not res["convention_ok"]:
                print(f"  !! CONVENTION MISMATCH: {res['convention_detail']}")
    return out


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Type-1 N=3 LZ gold-standard oracle")
    ap.add_argument("--T", type=float, default=80.0)
    ap.add_argument("--rtol", type=float, default=1e-13)
    ap.add_argument("--atol", type=float, default=1e-14)
    ap.add_argument("--only", type=str, default=None,
                    help="comma-separated stratum names (default: all)")
    args = ap.parse_args()
    names = args.only.split(",") if args.only else None
    run_suite(T=args.T, rtol=args.rtol, atol=args.atol, names=names)
