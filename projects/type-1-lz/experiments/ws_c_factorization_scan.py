"""
ws_c_factorization_scan.py  --  WS-C: locate the FACTORIZATION locus of the Type-1 N=3 MLZ
transition matrix P, and test P-ii (factorizable == joint-free == trivial-middle-coupling).

Objective (program prediction P-ii): determine the exact locus in {gamma,eps,a} where the
middle survival P_mm reduces to the incoherent product of two-level LZ factors
   P_mid^{inc} = exp(-2 pi (Gamma_mid,lo + Gamma_mid,hi)),
and correlate that against (a) the sep/width separation ratio, (b) the minimum middle
coupling strength min(Gamma_mid,lo, Gamma_mid,hi), and (c) WS-G joint count.

Metrics per sample:
  enh        = P_mid(exact oracle) / P_mid(incoherent)         [the headline; >1 = non-factorizing]
  logmiss    = |log10 P_mid(exact) - log10 P_mid(inc)|         [order-robust factorization defect]
  absdef     = |P_mid(exact) - P_mid(inc)|                      [absolute defect, the *honest* one]
  Gmid_lo, Gmid_hi   = the two middle-level BE exponents
  Gmid_prod  = Gmid_lo * Gmid_hi                                [trivial-coupling diagnostic; ->0]
  Gmid_min   = min(Gmid_lo, Gmid_hi)                            [weak-link diagnostic]
  ratio      = sep/width (coscaling)
  Gmin,Gmax  = min/max over all three Gamma_ij (genuine-LZ band check)

Reproduce: python ws_c_factorization_scan.py  (writes /tmp/ws_c_scan.npz + prints summary).
numpy/scipy + oracle.py (imports the WS-F gold oracle).
"""
from __future__ import annotations
import numpy as np
from oracle import oracle_P, Gamma_ij, be_survivals
from coscaling import geometry as cos_geometry


# ---------------------------------------------------------------------------
def sample_metrics(eps, gam, a, T=60.0, full=True, rtol=1e-12, atol=1e-13):
    """Return a dict of factorization metrics for one parameter point.
    full=True calls the gold oracle (slow); full=False returns only the cheap
    geometric/BE diagnostics (no ODE)."""
    eps = np.asarray(eps, float); gam = np.asarray(gam, float); a = np.asarray(a, float)
    bes = be_survivals(eps, gam, a)
    lo, mid, hi = bes["lo"], bes["mid"], bes["hi"]
    Gmid_lo = Gamma_ij(eps, gam, a, mid, lo)
    Gmid_hi = Gamma_ij(eps, gam, a, mid, hi)
    # all three Gamma_ij
    Gall = [Gamma_ij(eps, gam, a, i, j) for (i, j) in [(0, 1), (0, 2), (1, 2)]]
    # sep/width ratio (coscaling geometry)
    u, w, G = cos_geometry(gam, eps, a)
    uv = np.sort(list(u.values()))
    ratio = min(uv[1] - uv[0], uv[2] - uv[1]) / max(w.values())
    P_inc = bes["P_mid_inc"]
    out = dict(lo=lo, mid=mid, hi=hi,
               Gmid_lo=Gmid_lo, Gmid_hi=Gmid_hi,
               Gmid_prod=Gmid_lo * Gmid_hi, Gmid_min=min(Gmid_lo, Gmid_hi),
               Gmin=min(Gall), Gmax=max(Gall), ratio=ratio, P_inc=P_inc)
    if full:
        r = oracle_P(tuple(eps), tuple(gam), tuple(a), T=T, rtol=rtol, atol=atol)
        P_mid = float(r["P"][mid, mid])
        out["P_mid"] = P_mid
        out["enh"] = P_mid / P_inc if P_inc > 0 else np.inf
        out["absdef"] = abs(P_mid - P_inc)
        # order-robust log defect (guard underflow)
        lpm = np.log10(max(P_mid, 1e-300)); lpi = np.log10(max(P_inc, 1e-300))
        out["logmiss"] = abs(lpm - lpi)
        out["ds_defect"] = r["doubly_stochastic_defect"]
        out["be_delta"] = max(r["be_delta"]["lo"], r["be_delta"]["hi"])
        out["err"] = float(np.max(r["err"]))
    return out


# ---------------------------------------------------------------------------
def named_strata():
    """The WS-F/WS-G named strata + the WS-G sampleA (factorizing, separated)."""
    return {
        "canonical":       dict(eps=(-2, 0, 3),       gam=(1, 0.8, 1.2),         a=(-1, 0.5, 2)),
        "sampleB_overlap": dict(eps=(-1, 0, 1.5),     gam=(0.9, 1.1, 0.8),       a=(-0.7, 0.4, 1.3)),
        "wsg_sampleA_sep": dict(eps=(-9.173, 5.962, 8.603), gam=(-0.453, -0.105, 0.981), a=(-2.297, -1.075, 0.406)),
        "well_separated":  dict(eps=(-5, 0, 5),       gam=(1, 1, 1),             a=(-1.5, 0, 1.5)),
        "weak_coupling":   dict(eps=(-2, 0, 3),       gam=(0.35, 0.30, 0.40),    a=(-1, 0.5, 2)),
        "near_deg_slope":  dict(eps=(-2, 0, 3),       gam=(1, 0.8, 1.2),         a=(-1, 0.45, 0.55)),
    }


if __name__ == "__main__":
    print("=== Named strata: enhancement vs trivial-middle-coupling diagnostics ===")
    hdr = f"{'name':16s} {'enh':>9s} {'absdef':>9s} {'Gmid_lo':>9s} {'Gmid_hi':>9s} {'Gmid_min':>9s} {'ratio':>7s}"
    print(hdr)
    for name, s in named_strata().items():
        # wsg_sampleA has eps spanning ~18 -> needs relaxed rtol / smaller T to be fast;
        # absdef-level factorization classification is robust to this.
        wide = (max(s["eps"]) - min(s["eps"])) > 8.0
        T = 25.0 if wide else 50.0
        rtol = 1e-9 if wide else 1e-11
        m = sample_metrics(s["eps"], s["gam"], s["a"], T=T, rtol=rtol, atol=rtol * 1e-1)
        print(f"{name:16s} {m['enh']:9.3f} {m['absdef']:9.2e} {m['Gmid_lo']:9.4f} "
              f"{m['Gmid_hi']:9.4f} {m['Gmid_min']:9.4f} {m['ratio']:7.3f}", flush=True)
