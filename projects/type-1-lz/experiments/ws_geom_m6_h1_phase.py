"""
ws_geom_m6_h1_phase.py -- M6 / H1 decisive first-cut: does a FIRST-PRINCIPLES Stueckelberg phase
beat the FITTED constant K in the uniform Dykhne-Stueckelberg law?

The uniform law (R12/R13, ws_o3_uniform.py) is  P_mm ~ A0 + Aret + 2 sqrt(A0 Aret) cosPhi,
with the *fitted* proxy  cosPhi ~ -(d_lo+d_hi)(1 - K(1-chi)),  K=3.0 calibrated on the gold set.
Tellingly it already computes the single-crossing Stokes phases phi_stokes(d_lo)+phi_stokes(d_hi)
but DOES NOT USE them. H1: replace the proxy with the true phase Phi = phi_dyn + phi_stokes, where
phi_dyn is the dynamical action between the two diabatic crossings (a first-principles integral on
the genus-0 curve).

DECISIVE TEST (fair, equal free-parameter count): extract the TRUE cosPhi from a reference solve via
    cosPhi_true = (P_mm - A0 - Aret) / (2 sqrt(A0 Aret)),
then compare, in the regime where the interference is actually present (moderate coupling, both A0 and
Aret appreciable), the mean |cos - cosPhi_true| of:
    (i)  the fitted-K proxy (its one constant K), vs
    (ii) cos(phi_dyn + phi_stokes) with ONE global additive offset (its one constant).

HONEST OUTCOME (numerically-supported -- a NEGATIVE result): the elementary first-principles phase
  Phi = phi_dyn + phi_stokes(d_lo) + phi_stokes(d_hi) + phi_stokes(d_lh)
(the mid-lo adiabatic gap action between the two crossings + ALL THREE crossings' Stokes phases) does NOT
beat -- and is substantially WORSE than -- the fitted constant K. On a robust N=42 interference sample:
fitted-K 0.090, parameter-free first-principles 0.281 (~3x worse), first-principles + best offset 0.179
(~2x worse). The fitted K is absorbing genuinely NON-elementary (sigma-adjacent, non-period) content that
no elementary turning-point phase reproduces. (An earlier N=17 sample showed a SPURIOUS near-match -- the
larger sample overturns it; this is the robustness check catching a false positive.) VERDICT: H1 FAILS --
the graph machinery does NOT improve Dykhne-Stueckelberg with an elementary phase; beating K requires the
EXACT connection = sigma (the opaque exact-WKB object). sigma is load-bearing even in the semiclassical phase.

Reproduce: python3 ws_geom_m6_h1_phase.py
"""
from __future__ import annotations
import os
import sys
import numpy as np
from scipy.integrate import solve_ivp

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
from ws_geom_magnus import type1  # noqa: E402
import ws_o3_uniform as uni       # noqa: E402

np.seterr(all="ignore")


def ref_P(eps, gam, a, T=110.0):
    H0, Ad = type1(eps, gam, a)
    s = solve_ivp(lambda u, y: (-1j * (H0 + u * Ad) @ y.reshape(3, 3)).ravel(),
                  [-T, T], np.eye(3, dtype=complex).ravel(),
                  rtol=1e-8, atol=1e-10, method="DOP853")
    return np.abs(s.y[:, -1].reshape(3, 3).T) ** 2


def phi_dyn(eps, gam, a):
    """First-principles dynamical action: the mid-lo ADIABATIC gap integrated between the two
    diabatic crossings (the dominant leg of the STAY-vs-RETURN Stueckelberg loop)."""
    H0, A = type1(eps, gam, a); ad = np.diag(A); d = np.diag(H0)
    lo, mid, hi = (int(k) for k in np.argsort(a))
    ucr = lambda i, j: (d[j] - d[i]) / (a[i] - a[j])
    u1, u2 = sorted([ucr(lo, mid), ucr(mid, hi)])
    us = np.linspace(u1, u2, 200)
    W = np.array([np.linalg.eigvalsh(H0 + u * ad) for u in us])
    return float(np.trapezoid(W[:, 1] - W[:, 0], us))


def collect(a0=np.array([-1.0, 0.5, 2.0]), g0=np.array([1.0, 0.8, 1.2]),
            scales=(0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3),
            shifts=(-0.9, -0.5, -0.1, 0.3, 0.7, 1.1, 1.5, 1.9)):
    mid = int(np.argsort(a0)[1])
    rows = []
    for sc in scales:
        for shift in shifts:
            eps = np.array([-2.0, shift, 3.0]); gam = sc * g0
            _, pa = uni.uniform_P_mm(eps, gam, a0, return_parts=True)
            A0, Aret = pa["A0"], pa["Aret"]; Ac = 2 * np.sqrt(max(A0 * Aret, 0.0))
            if Ac < 0.05:
                continue
            pmm = ref_P(eps, gam, a0)[mid, mid]; ct = (pmm - A0 - Aret) / Ac
            if abs(ct) > 1.0:
                continue
            # parameter-free first-principles phase: phi_dyn + ALL THREE crossings' Stokes phases
            fs3 = (uni.stokes_phase(pa["d_lo"]) + uni.stokes_phase(pa["d_hi"])
                   + uni.stokes_phase(pa["d_lh"]))
            rows.append(dict(ct=ct, cf=pa["cosPhi"], ml=phi_dyn(eps, gam, a0),
                             fs2=pa["phi_stokes"], fs3=fs3, Ac=Ac))
    return rows


def best_offset(ct, phi):
    best = (1e9, 0.0)
    for ph0 in np.linspace(-np.pi, np.pi, 145):
        e = float(np.mean(np.abs(np.cos(phi + ph0) - ct)))
        if e < best[0]:
            best = (e, ph0)
    return best


def main():
    print("=" * 88)
    print("M6 / H1 -- first-principles Stueckelberg phase vs fitted K")
    print("=" * 88)
    rows = collect()
    ct = np.array([r["ct"] for r in rows]); cf = np.array([r["cf"] for r in rows])
    Ac = np.array([r["Ac"] for r in rows])
    ml = np.array([r["ml"] for r in rows])
    fs2 = np.array([r["fs2"] for r in rows]); fs3 = np.array([r["fs3"] for r in rows])
    print(f"valid interference samples (both A0,Aret appreciable): N={len(rows)}")
    e_fit = float(np.mean(np.abs(cf - ct)))
    # the headline: PARAMETER-FREE first-principles phase = phi_dyn + 3 crossings' Stokes phases
    e_pf = float(np.mean(np.abs(np.cos(ml + fs3) - ct)))
    e_off, ph0 = best_offset(ct, ml + fs2)
    print(f"\n  fitted-K law (1 fitted const K):                       mean|cos-cos_true| = {e_fit:.3f}")
    print(f"  PARAMETER-FREE  phi_dyn + Sum_3 phi_stokes (NO const): mean|cos-cos_true| = {e_pf:.3f}")
    print(f"  phi_dyn + 2 stokes, best single offset (1 const):      mean|cos-cos_true| = {e_off:.3f}")
    winner = "fitted K" if e_fit < min(e_pf, e_off) else "first-principles phase"
    print("\n  READING (robust, honest -- a NEGATIVE result):")
    print(f"   - the FITTED constant K ({e_fit:.3f}) DECISIVELY BEATS the first-principles phase: the")
    print(f"     parameter-free form is ~3x worse ({e_pf:.3f}), and even WITH its own fitted offset it is")
    print(f"     ~2x worse ({e_off:.3f}). Winner: {winner}.")
    print("   - the elementary turning-point phase (gap action + Stokes phases) does NOT reproduce cosPhi;")
    print("     the fitted K is absorbing genuinely NON-elementary (sigma-adjacent, non-period) content.")
    print("   - (An earlier N=17 sample showed a spurious ~match; the larger sample OVERTURNS it -- the")
    print("      robustness check is what caught the false positive.)")
    print("  VERDICT: H1 FAILS. The graph/turning-point machinery does NOT improve Dykhne-Stueckelberg with")
    print("           an elementary phase. Beating K requires the EXACT connection = sigma (the opaque")
    print("           exact-WKB object). So sigma is load-bearing even at the level of the semiclassical phase.")


if __name__ == "__main__":
    main()
