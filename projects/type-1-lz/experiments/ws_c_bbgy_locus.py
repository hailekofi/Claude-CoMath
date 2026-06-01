"""
ws_c_bbgy_locus.py  --  WS-C step 2(b): test whether the BBGY symmetric reduction
eps_2 = (eps_1 + eps_3)/2 is a FACTORIZING / elementary locus for the LINEAR Type-1 model.

BBGY (Bychkov-...-Yuzbashyan-class hyperbolic 3-state model) reduces their 3x3 Kampe de
Feriet to a 1F2 on the symmetric line.  We ask: does the *linear* Type-1 model inherit an
elementary (incoherent-product) middle survival when the middle pole sits at the midpoint
eps_2 = (eps_1+eps_3)/2 ?

We measure absdef = |P_mid(oracle) - P_mid(incoherent)| on the symmetric slice, scanning the
*other* free parameters (gam, a), and compare to off-symmetric controls.  If the symmetric
line were a genuine elementary locus, absdef -> 0 there independent of gam,a.

Reproduce: python ws_c_bbgy_locus.py.  numpy + oracle.py.
"""
from __future__ import annotations
import numpy as np
from oracle import oracle_P, be_survivals, Gamma_ij
from coscaling import geometry as cos_geometry


def absdef(eps, gam, a, T=40.0):
    # absdef-level factorization classification only needs ~1e-4 accuracy; use a
    # relaxed tolerance / shorter window so the scan is fast.
    bes = be_survivals(eps, gam, a)
    mid = bes["mid"]
    r = oracle_P(tuple(eps), tuple(gam), tuple(a), T=T, rtol=1e-10, atol=1e-11)
    Pm = float(r["P"][mid, mid])
    return Pm, bes["P_mid_inc"], abs(Pm - bes["P_mid_inc"])


def ratio_of(eps, gam, a):
    u, w, G = cos_geometry(gam, eps, a)
    uv = np.sort(list(u.values()))
    return min(uv[1] - uv[0], uv[2] - uv[1]) / max(w.values())


if __name__ == "__main__":
    rng = np.random.default_rng(7)
    print("=== BBGY symmetric slice eps_2 = (eps_1+eps_3)/2  vs  off-symmetric controls ===")
    print(f"{'case':22s} {'eps2':>8s} {'P_mid':>9s} {'P_inc':>9s} {'absdef':>10s} {'ratio':>7s}")

    e1, e3 = -2.0, 3.0
    e2_sym = 0.5 * (e1 + e3)   # = 0.5
    # Pick several (gam,a) at moderate overlap so the prefactor is genuinely active.
    cases = [
        ("gA", (1.0, 0.8, 1.2), (-1.0, 0.5, 2.0)),
        ("gB", (1.1, 0.9, 1.0), (-0.8, 0.3, 1.5)),
    ]
    for tag, gam, a in cases:
        # symmetric
        eps_s = (e1, e2_sym, e3)
        Pm, Pi, d = absdef(eps_s, gam, a)
        print(f"{tag+'  SYMMETRIC':22s} {e2_sym:8.3f} {Pm:9.4f} {Pi:9.4f} {d:10.3e} {ratio_of(eps_s,gam,a):7.3f}", flush=True)
        # off-symmetric controls (same e1,e3; shift e2)
        for e2 in (e2_sym - 0.9, e2_sym + 0.9):
            eps_o = (e1, e2, e3)
            Pm, Pi, d = absdef(eps_o, gam, a)
            print(f"{tag+'  off(e2=%.2f)'%e2:22s} {e2:8.3f} {Pm:9.4f} {Pi:9.4f} {d:10.3e} {ratio_of(eps_o,gam,a):7.3f}", flush=True)
        print()
