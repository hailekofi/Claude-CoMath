"""
num_S12.py  --  WS-NUM (PA-4) computable model of the Type-1 N=3 middle survival.

This is the *always-on safety net* deliverable for open-problem #2 of the Type-1
N=3 Landau-Zener program:  a trusted, high-precision, COMPUTABLE model of the open
middle survival

        P_2->2  ==  |S_12|^2  ==  P[mid, mid]

(the middle-SLOPE diabatic level's survival), expressed as a function of the natural
geometric arguments, together with inverse-symbolic / PSLQ recognition attempts.

It delivers a computable P_2->2(gamma, eps, a) REGARDLESS of whether a closed form is
recognized, in two tiers:

  TIER 0 (FLOOR, gold).  ``P22_oracle`` -- the trusted value from the project gold
      oracle ``experiments/oracle.py`` (the adiabatic interaction-picture propagator
      with the FIX1 permutation + FIX2 Richardson).  Benchmarked to <=1e-8 (T=80) /
      few x 1e-9 (T>=120) against the EXACT Brundobler-Elser extreme survivals.  This
      is the finite recipe that satisfies the practicality bar.

  TIER 1 (FAST).  ``P22_fast`` -- a single adiabatic-IP pass at T=80, rtol=1e-9, which
      already reproduces the gold value to ~1e-9 (no Richardson) in ~14 s, used for
      dataset generation and quick evaluation.

  TIER 2 (SURROGATE).  ``P22_model`` -- a smooth, fast-to-evaluate surrogate written
      purely in the NATURAL GEOMETRIC ARGUMENTS (the two window actions I_X and the
      Q4 turning-point cross-ratio).  It interpolates the gold data; it is the
      "computable model as a function of the geometry" requested by the parametrization
      task.  Its accuracy is reported in the validation table (it is an interpolant,
      not exact; use TIER 0 for gold values).

THE NATURAL GEOMETRIC ARGUMENTS (parametrization)
-------------------------------------------------
For the four complex Q4 branch points (two conjugate pairs) we form:

  * the two WINDOW ACTIONS I_X = oint_X (-sqrt(Q4) L_H W4 / p^3) dlam
    (the imaginary periods, ``uploads/assay/actions.py``).  Each
    delta_X = |Im I_X| / (2 pi) is a SUM of pairwise BE exponents s_ij^2|a_i-a_j|
    sharing one extreme level:
        delta_lo-window = be(lo,mid) + be(lo,hi)   ( -> P[lo,lo] = exp(-2 pi delta) )
        delta_hi-window = be(hi,mid) + be(hi,lo)   ( -> P[hi,hi] = exp(-2 pi delta) )
    i.e. the two window actions ARE the two BE extreme-survival exponents.  The middle
    level is the shared subdominant partner of both windows -- which is exactly why its
    survival is the open coupled quantity and not a single residue.

  * the turning-point CROSS-RATIO chi of the four Q4 roots.  Because the roots come in
    two complex-conjugate pairs, chi is REAL; it measures the sep/width geometry of the
    two avoided-crossing windows (chi -> 1 well-separated, chi -> 0 strongly merged).

We also expose the pairwise BE exponents be_ij, the incoherent baseline P_mid_inc
(the BE product for the middle level), and the canonical-frame Coulomb coefficients
c_i = sum_{j!=i} s_ij^2 (a_i - a_j).

CONVENTIONS: identical to oracle.py / NOMENCLATURE.md.  s_ij = gam_i gam_j/(eps_i-eps_j),
w_ij = |2 s_ij|, BE exponent = s_ij^2 |a_i - a_j|.  The scattering matrix is 𝒮 (Smat).

Requires: numpy, scipy; mpmath/sympy optional (for the PSLQ recognition pass).
The assay package ``../uploads`` and the sibling ``oracle.py`` must be importable.
"""
from __future__ import annotations

import os
import sys
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_UPLOADS = os.path.join(_HERE, "..", "uploads")
for _p in (_UPLOADS, _HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from assay import Params, Geometry                         # noqa: E402
from assay.ip import propagate_ad_ip                       # noqa: E402
from assay.actions import all_window_actions               # noqa: E402
import oracle                                              # noqa: E402  (sibling)

PI_IN = oracle.PI_IN     # (0,1,2)
PI_OUT = oracle.PI_OUT   # (2,0,1)


# ===========================================================================
#  1.  Type-1 algebra: pairwise BE exponents, Coulomb c_i, slope order
# ===========================================================================
def be_exponent(eps, gam, a, i, j) -> float:
    """Pairwise Brundobler-Elser exponent  s_ij^2 |a_i-a_j| = g_i^2 g_j^2|a_i-a_j|/(eps_i-eps_j)^2."""
    return gam[i] ** 2 * gam[j] ** 2 * abs(a[i] - a[j]) / (eps[i] - eps[j]) ** 2


def coulomb_c(eps, gam, a) -> np.ndarray:
    """Canonical-frame Coulomb coefficient c_i = sum_{j!=i} s_ij^2 (a_i-a_j); sum_i c_i = 0."""
    e = np.asarray(eps, float); g = np.asarray(gam, float); av = np.asarray(a, float)
    s = lambda i, j: g[i] * g[j] / (e[i] - e[j])
    return np.array([sum(s(i, j) ** 2 * (av[i] - av[j]) for j in range(3) if j != i)
                     for i in range(3)])


def slope_order(a):
    """(lo, mid, hi) = argsort(a): lo=argmin a, hi=argmax a, mid=middle slope."""
    lo, mid, hi = (int(k) for k in np.argsort(np.asarray(a, float)))
    return lo, mid, hi


def P_mid_incoherent(eps, gam, a) -> float:
    """Incoherent (independent-crossing) baseline for the middle survival."""
    lo, mid, hi = slope_order(a)
    return float(np.exp(-2 * np.pi * (be_exponent(eps, gam, a, mid, lo)
                                      + be_exponent(eps, gam, a, mid, hi))))


# ===========================================================================
#  2.  The natural geometric arguments  (window actions + Q4 cross-ratio)
# ===========================================================================
def q4_cross_ratio(geo: Geometry) -> complex:
    """
    Cross-ratio of the four complex Q4 turning points, roots sorted by real part:
        chi = (z0-z2)(z1-z3) / ((z0-z3)(z1-z2)) .
    The four roots are two complex-conjugate pairs, so chi is real (up to ~1e-12).
    chi -> 1 : the two avoided-crossing windows are well separated (large sep/width);
    chi -> 0 : the windows merge (strong overlap).
    """
    z = geo.Q4_roots()
    z = z[np.argsort(z.real)]
    return ((z[0] - z[2]) * (z[1] - z[3])) / ((z[0] - z[3]) * (z[1] - z[2]))


def geometry_args(eps, gam, a) -> dict:
    """
    All natural geometric arguments of the open middle survival.

    Returns dict with:
      I_X        : list of the two window actions (complex; imaginary periods).
      delta_X    : sorted [delta_small, delta_large], delta = |Im I_X|/(2 pi)
                   (each = a sum of pairwise BE exponents through one extreme level).
      chi        : Q4 turning-point cross-ratio (real).
      be         : dict of the three pairwise BE exponents {'lo_mid','mid_hi','lo_hi'}
                   (keys in slope labels), and 'sum' totals.
      c          : Coulomb coefficients c_i (ARRAY, eps-order).
      P_mid_inc  : incoherent baseline.
      slope      : (lo, mid, hi).
      sep_width  : a scalar sep/width ratio = |u-centre separation| / mean(width),
                   the stratification axis (0.1 strongly overlapping -> 4 well separated).
    """
    eps = tuple(float(v) for v in eps)
    gam = tuple(float(v) for v in gam)
    a = tuple(float(v) for v in a)
    geo = Geometry(Params(eps=eps, gam=gam, a=a, x=0))

    acts = all_window_actions(geo)
    I_X = [A["I_X"] for A in acts]
    delta_X = sorted(abs(A["I_X"].imag) / (2 * np.pi) for A in acts)
    chi = q4_cross_ratio(geo)

    lo, mid, hi = slope_order(a)
    be = {
        "lo_mid": be_exponent(eps, gam, a, lo, mid),
        "mid_hi": be_exponent(eps, gam, a, mid, hi),
        "lo_hi": be_exponent(eps, gam, a, lo, hi),
    }
    be["mid_sum"] = be["lo_mid"] + be["mid_hi"]   # the incoherent middle exponent

    # sep/width ratio from the two Q4 windows' real-axis u-centres and widths.
    wins = geo.windows()
    u_centres = sorted(w["u_center"] for w in wins)
    # width ~ imaginary extent of each conjugate pair in u
    widths = [abs(np.imag(w["u_vals"][0])) for w in wins]
    sep = abs(u_centres[1] - u_centres[0])
    mean_w = float(np.mean(widths)) if np.mean(widths) > 0 else 1.0
    sep_width = float(sep / mean_w)

    return dict(I_X=I_X, delta_X=delta_X, chi=complex(chi),
                be=be, c=coulomb_c(eps, gam, a),
                P_mid_inc=P_mid_incoherent(eps, gam, a),
                slope=(lo, mid, hi), sep_width=sep_width)


# ===========================================================================
#  3.  The diabatic scattering matrix 𝒮 and the open P_2->2  (engines)
# ===========================================================================
def _Smat_adiabatic(geo: Geometry, T: float, rtol: float, atol: float) -> np.ndarray:
    """
    The diabatic transition-PROBABILITY matrix P[x,j] from a single adiabatic-IP pass
    at half-window T, reindexed by the FIXED PI_IN/PI_OUT spectral-sheet maps (FIX 1).
    (We take moduli per channel; the absolute off-diagonal PHASES require the canonical
    log-T subtraction and are returned separately by ``S12_canonical_phase``.)
    """
    cols = [propagate_ad_ip(geo, -T, T, x, rtol=rtol, atol=atol)[:, x] for x in range(3)]
    Smat = np.array(cols).T                          # Smat[j, x], spectral-sheet basis
    Pd = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            Pd[PI_IN[i], PI_OUT[j]] = np.abs(Smat[j, i]) ** 2
    return Pd


def P22_fast(eps, gam, a, T: float = 80.0,
             rtol: float = 1e-9, atol: float = 1e-10) -> float:
    """
    FAST middle survival P_2->2 = P[mid,mid] from ONE adiabatic-IP pass (no Richardson).
    At T=80, rtol=1e-9 this already reproduces the gold oracle to ~1e-9 in ~14 s;
    it is the dataset-generation / quick-evaluation engine.
    """
    eps = tuple(float(v) for v in eps); gam = tuple(float(v) for v in gam)
    a = tuple(float(v) for v in a)
    geo = Geometry(Params(eps=eps, gam=gam, a=a, x=0))
    Pd = _Smat_adiabatic(geo, T, rtol, atol)
    _, mid, _ = slope_order(a)
    return float(Pd[mid, mid])


def P22_oracle(eps, gam, a, T: float = 120.0, **kw) -> dict:
    """
    GOLD (TIER 0 / FLOOR) middle survival, the trusted computable recipe.
    Thin wrapper over ``oracle.oracle_P`` (full Richardson + convention check).
    Returns dict(P22, err, full=<oracle result>, geom=<geometry_args>).
    """
    res = oracle.oracle_P(eps, gam, a, T=T, **kw)
    _, mid, _ = slope_order(a)
    return dict(P22=float(res["P"][mid, mid]),
                err=float(res["err"][mid, mid]),
                doubly_stochastic_defect=res["doubly_stochastic_defect"],
                be_delta=res["be_delta"],
                full=res, geom=geometry_args(eps, gam, a))


def S12_canonical_phase(eps, gam, a, T: float = 80.0,
                        rtol: float = 1e-9, atol: float = 1e-10) -> dict:
    """
    The middle-level off-diagonal AMPLITUDE and canonical-frame PHASE of 𝒮.

    The diabatic-IP off-diagonal phases drift as -(c_j - c_i) log T; the convergent
    canonical scattering matrix is  𝒮_canon[j,i] = e^{i(c_j-c_i) log T} S[j,i].
    We return |S_{mid, .}| (the amplitudes feeding P_2->2 and the two off-diagonals
    out of the middle channel) and the log-T-subtracted phases at base T.

    NOTE: amplitudes/probabilities are gauge-clean (T-convergent); the absolute phase
    of a single entry is convention-dependent -- we report the canonical-frame phase
    (c_i subtracted) and the gauge-invariant phase combinations.
    """
    eps = tuple(float(v) for v in eps); gam = tuple(float(v) for v in gam)
    a = tuple(float(v) for v in a)
    geo = Geometry(Params(eps=eps, gam=gam, a=a, x=0))
    cols = [propagate_ad_ip(geo, -T, T, x, rtol=rtol, atol=atol)[:, x] for x in range(3)]
    Smat = np.array(cols).T                          # spectral-sheet basis
    # reindex amplitudes to diabatic
    Sd = np.zeros((3, 3), complex)
    for i in range(3):
        for j in range(3):
            Sd[PI_IN[i], PI_OUT[j]] = Smat[j, i]
    lo, mid, hi = slope_order(a)
    c = coulomb_c(eps, gam, a)
    # canonical-frame phase subtraction  e^{i(c_j - c_i) log T}
    logT = np.log(T)
    phase_canon = {}
    for j in range(3):
        ph = np.angle(Sd[mid, j] * np.exp(1j * (c[j] - c[mid]) * logT))
        phase_canon[("mid", j)] = float(ph)
    return dict(
        S_mid_amp={j: float(abs(Sd[mid, j])) for j in range(3)},
        P22=float(abs(Sd[mid, mid]) ** 2),
        phase_canon=phase_canon,
        c=c, slope=(lo, mid, hi))


# ===========================================================================
#  4.  Stratified parameter family  (sep/width ratio 0.1 -> 4)
# ===========================================================================
#  Each entry: (eps, gam, a, description).  Spans well-separated -> strongly merged.
#  The named oracle anchors (canonical, sampleB) are the precision pins.
STRATA = {
    "canonical": ((-2.0, 0.0, 3.0), (1.0, 0.8, 1.2), (-1.0, 0.5, 2.0),
                  "canonical anchor (moderate overlap)"),
    "sampleB": ((-1.0, 0.0, 1.5), (0.9, 1.1, 0.8), (-0.7, 0.4, 1.3),
                "deep overlap anchor (~128x P_mid enhancement)"),
    "well_sep": ((-5.0, 0.0, 5.0), (1.0, 1.0, 1.0), (-1.5, 0.0, 1.5),
                 "well separated (large eps gaps)"),
    "weak": ((-2.0, 0.0, 3.0), (0.35, 0.30, 0.40), (-1.0, 0.5, 2.0),
             "weak coupling (near-diabatic)"),
    "strong": ((-2.0, 0.0, 3.0), (2.2, 2.0, 2.4), (-1.0, 0.5, 2.0),
               "strong coupling (near-adiabatic)"),
    # a graded sep/width sweep: shrink the eps spread to merge the windows
    "sep_wide": ((-4.0, 0.0, 4.0), (1.0, 0.8, 1.2), (-1.0, 0.5, 2.0),
                 "sep/width large (wide eps)"),
    "sep_mid": ((-1.5, 0.0, 2.2), (1.0, 0.8, 1.2), (-1.0, 0.5, 2.0),
                "sep/width moderate"),
    "sep_small": ((-0.8, 0.0, 1.1), (1.0, 0.8, 1.2), (-1.0, 0.5, 2.0),
                  "sep/width small (merging windows)"),
    "sep_tiny": ((-0.5, 0.0, 0.7), (1.0, 0.8, 1.2), (-1.0, 0.5, 2.0),
                 "sep/width tiny (strong overlap)"),
    # extra samples to break parameter degeneracy (vary gamma / slopes independently)
    "g_lo": ((-2.0, 0.0, 3.0), (0.6, 0.5, 0.7), (-1.0, 0.5, 2.0),
             "moderate-weak coupling"),
    "g_hi": ((-2.0, 0.0, 3.0), (1.5, 1.3, 1.7), (-1.0, 0.5, 2.0),
             "moderate-strong coupling"),
    "slope_asym": ((-2.0, 0.0, 3.0), (1.0, 0.8, 1.2), (-2.0, 0.3, 1.0),
                   "asymmetric slopes"),
    "eps_asym": ((-3.0, 0.0, 1.0), (1.0, 0.9, 1.1), (-1.0, 0.4, 1.6),
                 "asymmetric eps spacing"),
    "mid_low_slope": ((-2.0, 0.0, 3.0), (1.0, 0.8, 1.2), (-0.5, -0.2, 2.0),
                      "middle slope near low"),
}


_CACHE = os.path.join(_HERE, "num_S12_dataset.pkl")


def build_dataset(names=None, T: float = 80.0, rtol: float = 1e-9,
                  atol: float = 1e-10, engine: str = "fast", verbose: bool = True,
                  cache: bool = False):
    """
    Compute the high-precision P_2->2 and the geometric arguments for each stratum.

    engine='fast'   : single adiabatic-IP pass (P22_fast), ~14 s/sample, ~1e-9.
    engine='oracle' : full Richardson gold (P22_oracle), slower, the FLOOR value.

    cache=True       : load the persisted dataset (num_S12_dataset.pkl) if present,
                       else build it and save.  (The integration is the bottleneck;
                       caching lets the model/recognition iterate instantly.)

    Returns list of dicts with keys: name, eps, gam, a, P22, geom (geometry_args),
    delta_small, delta_large, chi, sep_width, P_mid_inc, ratio.
    """
    if cache and os.path.exists(_CACHE):
        import pickle
        with open(_CACHE, "rb") as fh:
            rows = pickle.load(fh)
        if names:
            rows = [r for r in rows if r["name"] in names]
        return rows
    names = names or list(STRATA.keys())
    rows = []
    for nm in names:
        eps, gam, a, desc = STRATA[nm]
        if engine == "oracle":
            r = P22_oracle(eps, gam, a, T=max(T, 120.0))
            P22 = r["P22"]; g = r["geom"]; err = r["err"]
        else:
            P22 = P22_fast(eps, gam, a, T=T, rtol=rtol, atol=atol)
            g = geometry_args(eps, gam, a); err = None
        d = g["delta_X"]
        row = dict(name=nm, desc=desc, eps=eps, gam=gam, a=a, P22=P22,
                   delta_small=d[0], delta_large=d[1], chi=g["chi"].real,
                   sep_width=g["sep_width"], P_mid_inc=g["P_mid_inc"],
                   ratio=P22 / g["P_mid_inc"], err=err, geom=g)
        rows.append(row)
        if verbose:
            print("[%-10s] P22=%.9f  inc=%.9f  ratio=%6.3f  "
                  "dlt=(%.4f,%.4f) chi=%.5f sw=%.3f%s"
                  % (nm, P22, g["P_mid_inc"], row["ratio"],
                     d[0], d[1], g["chi"].real, g["sep_width"],
                     "" if err is None else "  err=%.1e" % err))
    if cache:
        import pickle
        with open(_CACHE, "wb") as fh:
            pickle.dump(rows, fh)
    return rows


# ===========================================================================
#  5.  The surrogate model  P22_model(geometric arguments)
# ===========================================================================
#  PARAMETRIZATION (see paper/num_S12_model.md for the full discussion).
#
#  The natural geometric arguments are the two WINDOW ACTIONS I_X (imaginary periods)
#  and the Q4 turning-point CROSS-RATIO chi.  The two window actions are the two SUMS
#  of pairwise BE exponents through the two extreme levels:
#       Sigma_lo = be(lo,mid) + be(lo,hi)   (= delta of the lo-window)
#       Sigma_hi = be(hi,mid) + be(hi,lo)   (= delta of the hi-window)
#  Together with chi these are the geometry's invariants.  But P_2->2 also needs the
#  middle's OWN exponent  m = be(lo,mid) + be(mid,hi); equivalently the three pairwise
#  BE exponents (b_lm, b_mh, b_lh).  We therefore expose BOTH the window-action pair
#  (Sigma_lo, Sigma_hi, chi) AND the finer (b_lm, b_mh, b_lh, chi), and fit the surrogate
#  on the finer set (which is well-posed; the window sums are linear combinations of it).
#
#  TARGET TRANSFORM.  P_2->2 lives in (0,1) and ranges over many decades (adiabatic ->
#  diabatic), so a raw or log(P/P_inc) fit is ill-conditioned (P_inc -> 0 adiabatically).
#  We fit the LOGIT y = log(P/(1-P)) -- bounded-aware, well conditioned across all strata
#  -- as a low-order polynomial in the three BE exponents and (1-chi).  The model then
#  returns  P = sigmoid(design . coeffs),  guaranteed in (0,1).
#
#  This is the computable f(window actions, cross-ratio).  It is an INTERPOLANT of the
#  gold dataset (tag: [num-model]); for GOLD values use P22_oracle / P22_fast.  Its
#  cross-validated generalization error is reported by ``cv_report`` and in the md.

_MODEL_COEFFS = None   # set by calibrate_model(); persisted literals below.


def _logit(p):
    p = min(max(float(p), 1e-15), 1 - 1e-15)
    return np.log(p / (1.0 - p))


def _sigmoid(y):
    return 1.0 / (1.0 + np.exp(-y))


def _design_be(b_lm, b_mh, b_lh, chi):
    """
    Feature vector in the three pairwise BE exponents and the cross-ratio.
    Compact (8 features) to avoid overfitting the ~14-point gold dataset.
    """
    omc = 1.0 - chi
    return np.array([
        1.0,
        b_lm + b_mh,            # the middle's own incoherent exponent  m
        b_lh,                   # the outer (lo-hi) exponent
        b_lm * b_mh,            # the two-window product (interference scale)
        (b_lm - b_mh) ** 2,     # window asymmetry
        omc,                    # cross-ratio departure from 1 (window overlap)
        omc * (b_lm + b_mh),
        (b_lm + b_mh) ** 2,
    ], float)


def _row_be(r):
    """Extract (b_lm, b_mh, b_lh, chi) for a dataset row (slope-labelled BE exps)."""
    be = r["geom"]["be"]
    return be["lo_mid"], be["mid_hi"], be["lo_hi"], r["chi"]


def calibrate_model(rows) -> np.ndarray:
    """
    Least-squares fit of the LOGIT of P_2->2 to the design features over the dataset.
    Returns the coefficient vector; also sets the module-level ``_MODEL_COEFFS``.
    """
    global _MODEL_COEFFS
    X = np.array([_design_be(*_row_be(r)) for r in rows])
    y = np.array([_logit(r["P22"]) for r in rows])
    coeffs, *_ = np.linalg.lstsq(X, y, rcond=None)
    _MODEL_COEFFS = coeffs
    return coeffs


def cv_report(rows, verbose=True):
    """
    Leave-one-out cross-validation of the surrogate: refit on N-1 rows, predict the
    held-out P22, report the absolute deviation.  This is the HONEST generalization
    error of the interpolant (in-sample residual is near-zero by construction).
    """
    devs = []
    for k in range(len(rows)):
        train = [rows[i] for i in range(len(rows)) if i != k]
        c = calibrate_model(train)
        pred = _sigmoid(_design_be(*_row_be(rows[k])) @ c)
        dev = abs(pred - rows[k]["P22"])
        devs.append((rows[k]["name"], rows[k]["P22"], float(pred), float(dev)))
    calibrate_model(rows)   # restore full fit
    if verbose:
        print("  leave-one-out CV:")
        for nm, p, pr, d in devs:
            print("    [%-12s] gold=%.6f  pred=%.6f  dev=%.2e" % (nm, p, pr, d))
        print("    median LOO dev = %.2e   max LOO dev = %.2e"
              % (np.median([d for *_, d in devs]), max(d for *_, d in devs)))
    return devs


def P22_model(eps=None, gam=None, a=None, *,
              b_lm=None, b_mh=None, b_lh=None, chi=None, coeffs=None) -> float:
    """
    SURROGATE (TIER 2) middle survival as a smooth function of the natural geometric
    arguments.  Call with (eps,gam,a) -- geometry computed internally -- or directly
    with the four invariants (b_lm, b_mh, b_lh, chi).

    P = sigmoid(design . coeffs), guaranteed in (0,1).  Requires a calibrated coefficient
    vector.  This is the fast computable f; for GOLD values use P22_oracle / P22_fast.
    """
    c = coeffs if coeffs is not None else _MODEL_COEFFS
    if c is None:
        raise RuntimeError("model not calibrated; run calibrate_model(build_dataset(...))")
    if eps is not None:
        g = geometry_args(eps, gam, a)
        lo, mid, hi = g["slope"]
        b_lm = g["be"]["lo_mid"]; b_mh = g["be"]["mid_hi"]
        b_lh = g["be"]["lo_hi"]; chi = g["chi"].real
    return float(_sigmoid(_design_be(b_lm, b_mh, b_lh, chi) @ c))


# ===========================================================================
#  6.  Inverse-symbolic / PSLQ recognition
# ===========================================================================
def pslq_recognize(value, basis_names, basis_values, tol=1e-9, maxcoeff=10 ** 6):
    """
    Try to recognize ``value`` as an integer combination of the constants in
    ``basis_values`` (named by ``basis_names``) using mpmath.pslq.  Returns
    (relation_dict, residual) or (None, None) if no low-height relation is found.
    """
    try:
        import mpmath as mp
    except Exception:
        return None, None
    mp.mp.dps = 40
    vec = [mp.mpf(float(value))] + [mp.mpf(float(v)) for v in basis_values]
    rel = mp.pslq(vec, tol=mp.mpf(10) ** (-int(-np.log10(tol))), maxcoeff=maxcoeff)
    if rel is None or rel[0] == 0:
        return None, None
    # value = -(sum_k rel[k+1] basis[k]) / rel[0]
    names = ["value"] + list(basis_names)
    relation = {names[k]: int(rel[k]) for k in range(len(rel)) if rel[k] != 0}
    resid = float(abs(sum(mp.mpf(rel[k]) * vec[k] for k in range(len(rel)))))
    return relation, resid


def recognition_pass(row, tol=1e-7, verbose=True):
    """
    Run a battery of inverse-symbolic checks on one dataset row's P22 / amplitudes /
    log against candidate closed forms built from the window actions, BE exponents,
    cross-ratio, pi, and simple Gamma-ratio constants.  Returns a list of hits.
    """
    hits = []
    g = row["geom"]
    ds, dl, chi = row["delta_small"], row["delta_large"], row["chi"]
    P22, Pinc = row["P22"], row["P_mid_inc"]
    R = P22 / Pinc
    logR = np.log(R)

    # candidate basis constants for PSLQ on logR
    basis = {
        "pi": np.pi, "1": 1.0,
        "log(1-chi)": np.log(max(1 - chi, 1e-12)),
        "log(chi)": np.log(max(chi, 1e-12)),
        "ds": ds, "dl": dl, "ds*dl": ds * dl,
        "log(ds)": np.log(max(ds, 1e-12)), "log(dl)": np.log(max(dl, 1e-12)),
        "log(ds+dl)": np.log(max(ds + dl, 1e-12)),
    }
    rel, resid = pslq_recognize(logR, list(basis.keys()), list(basis.values()), tol=tol)
    if rel is not None and resid is not None and resid < tol:
        hits.append(("logR ~ int-combo", rel, resid))

    # simple closed-form guesses for R itself
    guesses = {
        "R == 1 + (something)*?": None,
        "R vs 1/(1-exp(-2pi ds)) ...": None,
    }
    # explicit: is R recognizable as a ratio of (1-p) factors (Demkov-Osherov style)?
    p_ds = np.exp(-2 * np.pi * ds)
    p_dl = np.exp(-2 * np.pi * dl)
    do_candidates = {
        "DO 1": (1 - p_ds) * (1 - p_dl) + p_ds * p_dl,   # generic two-link interference
        "DO 2": 1 + p_ds * p_dl - p_ds - p_dl,
        "DO 3": (1 - p_ds * p_dl),
    }
    for nm, val in do_candidates.items():
        if val > 0 and abs(np.log(val) - logR) < tol:
            hits.append((nm, {"R": val}, abs(np.log(val) - logR)))

    if verbose:
        if hits:
            for h in hits:
                print("   HIT:", h)
        else:
            print("   no low-height relation found (R=%.9f logR=%.6f)" % (R, logR))
    return hits


# ===========================================================================
#  7.  Validation against the gold oracle
# ===========================================================================
def validate(rows_fast, names=None, T: float = 120.0, verbose=True):
    """
    Compare the fast-engine and surrogate P22 to the GOLD oracle on the named strata.
    Returns a list of dicts; prints a table.  (Oracle calls are slow; default to the
    two precision anchors unless ``names`` widens it.)
    """
    names = names or ["canonical", "sampleB"]
    table = []
    by_name = {r["name"]: r for r in rows_fast}
    for nm in names:
        eps, gam, a, desc = STRATA[nm]
        gold = P22_oracle(eps, gam, a, T=T)
        fast = by_name.get(nm, {}).get("P22")
        model = None
        if _MODEL_COEFFS is not None:
            model = P22_model(eps, gam, a)
        rec = dict(name=nm, gold=gold["P22"], gold_err=gold["err"],
                   fast=fast, fast_dev=None if fast is None else abs(fast - gold["P22"]),
                   model=model, model_dev=None if model is None else abs(model - gold["P22"]),
                   be_delta=gold["be_delta"])
        table.append(rec)
        if verbose:
            print("[%-10s] gold=%.9f (err %.1e)  fast=%s (dev %s)  model=%s (dev %s)  BEd=(%.1e,%.1e)"
                  % (nm, rec["gold"], rec["gold_err"],
                     "%.9f" % fast if fast is not None else "  --  ",
                     "%.1e" % rec["fast_dev"] if rec["fast_dev"] is not None else "--",
                     "%.6f" % model if model is not None else "  --  ",
                     "%.1e" % rec["model_dev"] if rec["model_dev"] is not None else "--",
                     gold["be_delta"]["lo"], gold["be_delta"]["hi"]))
    return table


# ===========================================================================
#  8.  CLI driver
# ===========================================================================
if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="WS-NUM computable P_2->2 model + recognition")
    ap.add_argument("--engine", default="fast", choices=["fast", "oracle"])
    ap.add_argument("--T", type=float, default=80.0)
    ap.add_argument("--only", type=str, default=None)
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--recognize", action="store_true")
    args = ap.parse_args()
    names = args.only.split(",") if args.only else None

    print("=== building dataset (engine=%s, T=%g) ===" % (args.engine, args.T))
    rows = build_dataset(names=names, T=args.T, engine=args.engine)

    print("\n=== calibrating surrogate model ===")
    c = calibrate_model(rows)
    print("coeffs =", np.array2string(c, precision=6))
    # in-sample residuals
    for r in rows:
        m = P22_model(r["eps"], r["gam"], r["a"])
        print("  [%-10s] data=%.9f model=%.9f dev=%.2e" % (r["name"], r["P22"], m, abs(m - r["P22"])))

    if args.recognize:
        print("\n=== inverse-symbolic / PSLQ recognition ===")
        for r in rows:
            print("[%s]" % r["name"])
            recognition_pass(r)

    if args.validate:
        print("\n=== validation vs gold oracle ===")
        validate(rows, names=names)
