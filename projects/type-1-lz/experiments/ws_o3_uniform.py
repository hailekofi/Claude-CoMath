"""
ws_o3_uniform.py  --  WS-O3 uniform semiclassical (Dykhne + Stueckelberg) model of the
Type-1 N=3 middle survival P_mm = P[mid,mid].

This is the deliverable for open-question O3 (SESSION_SYNTHESIS sec.6): a *uniform*
semiclassical, computable approximation to the open middle survival built from the
single-crossing Dykhne (Landau-Zener-Weber) building blocks composed with a Stueckelberg
interference phase that is DRESSED by the spectral-geometry shape (the cross-ratio chi),
so it does not blow up as Type-1's two avoided crossings marginally overlap (R4).

The construction (full derivation + evidence-ladder tags in
paper/ws_o3_uniform_asymptotics.md):

  The middle-slope level passes two avoided crossings, mid-lo and mid-hi, with the
  intermediate lo-hi crossing sitting between them in u (verified ordering). The two
  per-crossing adiabaticities (= BE/window actions = LZ adiabaticity parameters) are

        delta_lo = s_{mid,lo}^2 |a_mid - a_lo|   (= b_lm)
        delta_hi = s_{mid,hi}^2 |a_mid - a_hi|   (= b_mh)
        delta_lh = s_{lo ,hi }^2 |a_lo - a_hi |  (= b_lh)   (the outer link)

  with single-crossing diabatic STAY probabilities  p_X = exp(-2 pi delta_X).

  Two semiclassical Feynman paths return the middle to itself:
    (S) STAY diabatic at both its crossings        amplitude^2 = p_lo p_hi
    (R) jump out at one crossing, traverse the lo-hi link, jump back at the other
                                                    amplitude^2 = (1-p_lo)(1-p_hi) p_lh
  Uniform two-path survival:

      P_mm  ~  p_lo p_hi  +  (1-p_lo)(1-p_hi) p_lh
                 +  2 sqrt( p_lo p_hi (1-p_lo)(1-p_hi) p_lh ) * cosPhi      (*)

  The leading (parameter-free) UNIFORM term is the incoherent stay+return sum
  A0 + Aret; it is finite and well-behaved through the overlap (it does NOT use a
  divergent inter-crossing dynamical phase -- Type-1's crossings never separate enough
  to develop one; R4). The interference cosPhi is the OVERLAP-DRESSED correction:
  numerically it is NOT a rapidly oscillating Stueckelberg phase but a smooth, monotone
  function of the two window actions and the shape chi:

      cosPhi  ~  -(delta_lo + delta_hi) * ( 1 - K (1 - chi) )                (**)

  i.e. the two paths sit near quadrature (Phi ~ pi/2, cosPhi ~ 0 as delta -> 0) and tilt
  by an amount set by the middle's TOTAL window action (delta_lo+delta_hi) modulated by
  the SHAPE chi.  K is a single O(1) constant calibrated on the gold dataset.

This cleanly exhibits the R11 structure  P_mm = f({two window actions}, chi):
  - the two BE window actions enter through p_lo, p_hi (scale);
  - the cross-ratio chi enters ONLY through the shape factor (1 - K(1-chi)) (shape).

REGIME OF VALIDITY (honest):
  near-exact (<~1-2%) in the separated / weak-to-moderate regime (delta_lo,delta_hi <~ 0.25);
  ~3-6% at moderate overlap (delta <~ 0.45);
  FAILS (two-path truncation breaks; >~10%) in the deep adiabatic-overlap regime
  (delta >~ 0.5, e.g. sampleB, strong) where R4's permanent extreme-extreme coupling
  dominates and no per-crossing factorization is controlled.

Reuses the project builder `type1` and the gold oracle (experiments/oracle.py via
num_S12.P_mm_fast / oracle.oracle_P). Run:  python ws_o3_uniform.py
"""
from __future__ import annotations
import os, sys
import numpy as np
from scipy.special import loggamma

_HERE = os.path.dirname(os.path.abspath(__file__))
_UPLOADS = os.path.join(_HERE, "..", "uploads")
for _p in (_UPLOADS, _HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import oracle          # noqa: E402  gold oracle (sibling)
import num_S12         # noqa: E402  fast engine + geometry args (chi, deltas)


# ---------------------------------------------------------------------------
#  Project builder (use exactly as specified)
# ---------------------------------------------------------------------------
def type1(eps, gam, a):
    eps = np.array(eps, float); gam = np.array(gam, float); a = np.array(a, float)
    g2 = gam ** 2
    H0 = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            if i != j:
                H0[i, j] = gam[i] * gam[j] * (a[i] - a[j]) / (eps[i] - eps[j])
        H0[i, i] = -sum(g2[k] * (a[i] - a[k]) / (eps[i] - eps[k])
                        for k in range(3) if k != i)
    return H0, np.diag(a)


# ---------------------------------------------------------------------------
#  Geometric building blocks
# ---------------------------------------------------------------------------
def be_exponent(eps, gam, a, i, j) -> float:
    """Pairwise BE / window action / LZ adiabaticity  delta_ij = s_ij^2 |a_i-a_j|."""
    return gam[i] ** 2 * gam[j] ** 2 * abs(a[i] - a[j]) / (eps[i] - eps[j]) ** 2


def stokes_phase(delta: float) -> float:
    """Single-crossing Stueckelberg/Stokes phase  phi = d(ln d -1)+pi/4+arg Gamma_E(1-i d)."""
    if delta <= 0:
        return np.pi / 4
    return (delta * (np.log(delta) - 1.0) + np.pi / 4.0
            + float(np.angle(np.exp(loggamma(1.0 - 1j * delta)))))


# ---------------------------------------------------------------------------
#  The uniform formula
# ---------------------------------------------------------------------------
#  Calibrated overlap-dressing constant K in eq (**) -- the SINGLE fitted O(1) number.
#  Fit on the gold dataset (scale x gamma x shape sweep); see __main__ / the md.
K_OVERLAP = 3.0


def uniform_P_mm(eps, gam, a, K: float = K_OVERLAP, return_parts: bool = False):
    """
    Uniform semiclassical middle survival P_mm in {delta_lo, delta_hi, chi}.

    eps, gam, a : length-3 model params (eps strictly increasing).
    K           : overlap-dressing constant (eq **).
    Returns float P_mm, or (P_mm, dict-of-parts) if return_parts.
    """
    eps = np.array(eps, float); gam = np.array(gam, float); a = np.array(a, float)
    lo, mid, hi = (int(k) for k in np.argsort(a))            # slope order

    d_lo = be_exponent(eps, gam, a, mid, lo)
    d_hi = be_exponent(eps, gam, a, mid, hi)
    d_lh = be_exponent(eps, gam, a, lo, hi)

    p_lo = np.exp(-2 * np.pi * d_lo)
    p_hi = np.exp(-2 * np.pi * d_hi)
    p_lh = np.exp(-2 * np.pi * d_lh)

    A0 = p_lo * p_hi                                  # STAY path
    Aret = (1 - p_lo) * (1 - p_hi) * p_lh             # RETURN-through-link path
    Ac = 2.0 * np.sqrt(max(A0 * Aret, 0.0))           # coherent amplitude

    chi = float(num_S12.q4_cross_ratio(
        num_S12.Geometry(num_S12.Params(eps=tuple(eps), gam=tuple(gam),
                                        a=tuple(a), x=0))).real)

    # uniform overlap-dressed interference (eq **): near-quadrature + shape tilt
    cosPhi = -(d_lo + d_hi) * (1.0 - K * (1.0 - chi))
    cosPhi = float(np.clip(cosPhi, -1.0, 1.0))        # uniform clip (keeps it physical)

    P = A0 + Aret + Ac * cosPhi
    P = float(np.clip(P, 0.0, 1.0))
    if return_parts:
        return P, dict(d_lo=d_lo, d_hi=d_hi, d_lh=d_lh, chi=chi,
                       A0=A0, Aret=Aret, Ac=Ac, cosPhi=cosPhi,
                       phi_stokes=stokes_phase(d_lo) + stokes_phase(d_hi))
    return P


def _incoh(eps, gam, a) -> float:
    """Parameter-free uniform leading term  A0 + Aret  (no fitted constant, no chi)."""
    eps = np.array(eps, float); gam = np.array(gam, float); a = np.array(a, float)
    lo, mid, hi = (int(k) for k in np.argsort(a))
    p_lo = np.exp(-2 * np.pi * be_exponent(eps, gam, a, mid, lo))
    p_hi = np.exp(-2 * np.pi * be_exponent(eps, gam, a, mid, hi))
    p_lh = np.exp(-2 * np.pi * be_exponent(eps, gam, a, lo, hi))
    return float(p_lo * p_hi + (1 - p_lo) * (1 - p_hi) * p_lh)


# ---------------------------------------------------------------------------
#  Validation sweep against the oracle
# ---------------------------------------------------------------------------
# Gold anchors (oracle_report.md / num_S12.GOLD_ANCHORS), exact BE-verified.
GOLD_ANCHORS = dict(num_S12.GOLD_ANCHORS)


def _build_sweep():
    """
    A reproducible sweep across sep/width and a chi sweep, spanning the regime
    boundary.  Returns list of (name, eps, gam, a).  Uses the project's STRATA plus a
    scale x gamma family (fixed shape) and a shape family (moves chi).
    """
    rows = []
    for nm, (eps, gam, a, desc) in num_S12.STRATA.items():
        rows.append((nm, eps, gam, a))
    base_eps = np.array((-2.0, 0.0, 3.0)); base_gam = np.array((1.0, 0.8, 1.2))
    base_a = (-1.0, 0.5, 2.0)
    for s in (0.8, 1.0, 1.3, 1.7, 2.2, 3.0):
        for gs in (0.7, 0.85, 1.0, 1.2):
            rows.append(("scl_s%.1f_g%.2f" % (s, gs),
                         tuple(s * base_eps), tuple(gs * base_gam), base_a))
    shapes = [
        ((-3.0, 0.0, 1.0), (1.0, 0.9, 1.1), (-1.0, 0.4, 1.6)),
        ((-1.0, 0.0, 4.0), (1.0, 0.8, 1.2), (-1.0, 0.5, 2.0)),
        ((-2.0, 0.0, 3.0), (1.0, 0.8, 1.2), (-2.0, 0.3, 1.0)),
        ((-2.5, 0.0, 2.5), (1.0, 0.7, 1.3), (-1.2, 0.6, 1.8)),
        ((-1.5, 0.0, 3.5), (0.9, 0.8, 1.1), (-1.3, 0.2, 1.7)),
    ]
    for k, (eps, gam, a) in enumerate(shapes):
        for s in (1.5, 2.0, 2.8):
            rows.append(("shp%d_s%.1f" % (k, s), tuple(s * np.array(eps)), gam, a))
    return rows


_GOLD_CACHE = os.path.join(_HERE, "ws_o3_gold_cache.pkl")


def _gold_values(rows, engine="fast", T=70.0, use_cache=True):
    """Compute (and cache) gold P_mm for each sweep row. Caching keyed by (name,eps,gam,a)."""
    import pickle
    cache = {}
    if use_cache and os.path.exists(_GOLD_CACHE):
        with open(_GOLD_CACHE, "rb") as fh:
            cache = pickle.load(fh)
    out = {}
    dirty = False
    for (nm, eps, gam, a) in rows:
        key = (nm, tuple(np.round(eps, 9)), tuple(np.round(gam, 9)), tuple(np.round(a, 9)),
               engine)
        if key in cache:
            out[nm] = cache[key]
            continue
        if engine == "oracle":
            r = oracle.oracle_P(eps, gam, a, T=max(T, 120.0))
            lo, mid, hi = (int(k) for k in np.argsort(np.array(a, float)))
            val = float(r["P"][mid, mid])
        else:
            val = num_S12.P_mm_fast(eps, gam, a, T=T, rtol=3e-9)
        cache[key] = val
        out[nm] = val
        dirty = True
    if use_cache and dirty:
        with open(_GOLD_CACHE, "wb") as fh:
            pickle.dump(cache, fh)
    return out


def run_validation(engine: str = "fast", T: float = 70.0, verbose: bool = True):
    """
    Evaluate the uniform formula vs the oracle across the sweep; print an error table
    stratified by max(delta) (the regime axis) and report max/RMS errors.

    engine='fast'   : num_S12.P_mm_fast (single IP pass, ~1e-9 vs gold; fast).
    engine='oracle' : full Richardson gold (slow).
    """
    rows = _build_sweep()
    golds = _gold_values(rows, engine=engine, T=T)
    recs = []
    for (nm, eps, gam, a) in rows:
        P_gold = golds[nm]
        P_u, parts = uniform_P_mm(eps, gam, a, return_parts=True)
        P_inc = _incoh(eps, gam, a)
        dmax = max(parts["d_lo"], parts["d_hi"])
        recs.append(dict(name=nm, gold=P_gold, uniform=P_u, incoherent=P_inc,
                         dmax=dmax, chi=parts["chi"], d_lo=parts["d_lo"],
                         d_hi=parts["d_hi"], err_u=abs(P_u - P_gold),
                         err_inc=abs(P_inc - P_gold)))

    if verbose:
        print("\n%-16s %8s %8s %8s | %8s %8s  d_lo  d_hi   chi" %
              ("name", "gold", "uniform", "incoh", "err_uni", "err_inc"))
        for r in sorted(recs, key=lambda r: r["dmax"]):
            print("%-16s %8.5f %8.5f %8.5f | %8.4f %8.4f  %.3f %.3f  %.3f" %
                  (r["name"], r["gold"], r["uniform"], r["incoherent"],
                   r["err_u"], r["err_inc"], r["d_lo"], r["d_hi"], r["chi"]))

    print("\n=== accuracy by regime (regime axis = max(delta_lo, delta_hi)) ===")
    print("%-22s %4s  %-22s  %-22s" % ("regime", "n", "UNIFORM (max / rms)",
                                       "incoherent (max / rms)"))
    bins = [(0.0, 0.10, "separated   d<0.10"),
            (0.0, 0.25, "weak-mod    d<0.25"),
            (0.0, 0.45, "moderate    d<0.45"),
            (0.45, 1e9, "deep overlap d>0.45"),
            (0.0, 1e9, "ALL")]
    for lo_, hi_, label in bins:
        sub = [r for r in recs if lo_ <= r["dmax"] < hi_] if lo_ > 0 else \
              [r for r in recs if r["dmax"] < hi_]
        if not sub:
            continue
        eu = np.array([r["err_u"] for r in sub])
        ei = np.array([r["err_inc"] for r in sub])
        print("%-22s %4d  %7.4f / %7.4f      %7.4f / %7.4f" %
              (label, len(sub), eu.max(), np.sqrt(np.mean(eu ** 2)),
               ei.max(), np.sqrt(np.mean(ei ** 2))))
    return recs


def check_anchors(verbose: bool = True):
    """Check the uniform formula directly on the gold anchors (slope-middle survival)."""
    out = {}
    for nm, gold in GOLD_ANCHORS.items():
        eps, gam, a, desc = num_S12.STRATA[nm]
        P_u, parts = uniform_P_mm(eps, gam, a, return_parts=True)
        P_inc = _incoh(eps, gam, a)
        out[nm] = dict(gold=gold, uniform=P_u, incoherent=P_inc,
                       err=abs(P_u - gold), dmax=max(parts["d_lo"], parts["d_hi"]),
                       chi=parts["chi"])
        if verbose:
            print("[%-12s] gold=%.6f uniform=%.6f (err %.4f) incoh=%.6f  "
                  "dmax=%.3f chi=%.3f" %
                  (nm, gold, P_u, abs(P_u - gold), P_inc, out[nm]["dmax"], out[nm]["chi"]))
    return out


def recalibrate_K(engine="fast", T=70.0):
    """
    Re-fit the single overlap constant K (eq **) on the separated+moderate sweep
    (delta<0.45) by least squares in P-space.  Returns K_opt and the achieved errors.
    """
    from scipy.optimize import minimize_scalar
    rows = _build_sweep()
    golds = _gold_values(rows, engine=engine, T=T)
    data = []
    for (nm, eps, gam, a) in rows:
        _, parts = uniform_P_mm(eps, gam, a, return_parts=True)
        if max(parts["d_lo"], parts["d_hi"]) < 0.45:
            data.append((eps, gam, a, golds[nm]))

    def loss(K):
        return float(np.mean([(uniform_P_mm(e, g, a_, K=K) - P) ** 2
                              for (e, g, a_, P) in data]))
    res = minimize_scalar(loss, bounds=(-2, 8), method="bounded")
    return res.x, np.sqrt(res.fun)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="WS-O3 uniform semiclassical P_mm")
    ap.add_argument("--engine", default="fast", choices=["fast", "oracle"])
    ap.add_argument("--T", type=float, default=70.0)
    ap.add_argument("--recalibrate", action="store_true")
    args = ap.parse_args()

    print("=== gold anchors (slope-middle survival) ===")
    check_anchors()

    if args.recalibrate:
        Kopt, rmsfit = recalibrate_K(engine=args.engine, T=args.T)
        print("\nrecalibrated K = %.4f  (fit rms in P over delta<0.45 = %.4f); "
              "module uses K=%.2f" % (Kopt, rmsfit, K_OVERLAP))

    print("\n=== full validation sweep vs oracle (engine=%s) ===" % args.engine)
    run_validation(engine=args.engine, T=args.T)
