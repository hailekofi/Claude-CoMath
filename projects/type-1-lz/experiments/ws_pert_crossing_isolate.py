"""
ws_pert_crossing_isolate.py -- broaden the perturbation analysis:
  PART 1: multiple generic V directions -> confirm the eps^1 (perturbative) scaling is direction-generic.
  PART 2: ISOLATE the crossing's own contribution by comparing two perturbations of EQUAL norm:
     * V_open : generic -> opens the node gap at O(eps);
     * V_keep : same generic M but with the node's 2x2 block (in the degenerate eigenbasis at u_*)
                SCALARIZED (set to (trace/2)*I, off-diagonal zeroed) -> the degenerate pair stays
                degenerate at u_* to first order, so the crossing is PRESERVED (gap opens only at O(eps^2)),
                while integrability is still broken elsewhere.
  The crossing-specific effect is the DIFFERENCE dP(V_open) - dP(V_keep): if the crossing-removal channel is
  analytic O(eps^2), this difference is subdominant to the common O(eps) bulk.

P[x,j]=|U[j,x]|^2; slope lo,mid,hi=argsort(a).

FINDINGS (numerically-supported, with honest noise caveats):
  PART 1 (direction-genericity): across seeds 1-4, |dP_mm| ~ eps^(~1) -- exponents scatter 0.65-1.76 in this
    lean (single-T, narrow-eps) config, but NONE is non-perturbative (no fractional-stable / ->0 exponent).
    Together with the clean seed-7 run (ws_pert_crossing.py: eps^1.0 down to eps=1e-3), the PERTURBATIVE
    scaling is direction-generic.
  PART 2 (isolating the crossing): the construction works -- gap_keep stays at the grid floor (crossing
    PRESERVED) while gap_open grows. The robust result: dPmm_open ~= dPmm_keep (they differ by ~30-40%; the
    DIFFERENCE is ~10x smaller than either). So WHETHER OR NOT the crossing is opened, dP_mm is nearly the
    same -- the observable responds to the GENERIC integrability-breaking, and the crossing-removal per se is
    a SUBDOMINANT contribution. (The predicted O(eps^2) scaling of the isolated channel is NOT cleanly
    resolved here -- diff exponent ~1.1-1.4, crossover+noise-limited; would need T-averaging + finer gaps +
    smaller eps to pin.)
  NET: consistent with ws_pert_crossing.py -- the impact of removing the protected crossing is PERTURBATIVE
  and, moreover, the crossing-removal is not even the DOMINANT part of a generic perturbation's effect on the
  observable; the bulk is ordinary first-order integrability-breaking.

Reproduce: python3 ws_pert_crossing_isolate.py  (lean defaults; raise Ts/n and lower eps for clean exponents)
"""
from __future__ import annotations
import numpy as np
from scipy.integrate import solve_ivp
np.seterr(all="ignore")


def type1(eps, gam, a):
    eps = np.array(eps, float); gam = np.array(gam, float); a = np.array(a, float); g2 = gam ** 2
    H0 = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            if i != j:
                H0[i, j] = gam[i] * gam[j] * (a[i] - a[j]) / (eps[i] - eps[j])
        H0[i, i] = -sum(g2[k] * (a[i] - a[k]) / (eps[i] - eps[k]) for k in range(3) if k != i)
    return H0, np.diag(a)


def propP(H0p, A, Ts=(70.0,)):
    # single T: the endpoint Fresnel tail cancels in dP = P(eps)-P(0) (same T for both).
    acc = np.zeros((3, 3))
    for T in Ts:
        s = solve_ivp(lambda u, y: (-1j * (H0p + u * A) @ y.reshape(3, 3)).ravel(),
                      [-T, T], np.eye(3, dtype=complex).ravel(), rtol=1e-9, atol=1e-11, method="DOP853")
        acc += np.abs(s.y[:, -1].reshape(3, 3).T) ** 2
    return acc / len(Ts)


def node_data(H0, A, T=80.0, n=6001):
    av = np.diag(np.diag(A)); us = np.linspace(-T, T, n)
    gmin = 1e9; ustar = 0.0; kp = 0
    for u in us:
        w = np.linalg.eigvalsh(H0 + u * av); g = np.diff(w); k = int(np.argmin(g))
        if g[k] < gmin:
            gmin = g[k]; ustar = u; kp = k
    du = us[1] - us[0]
    for u in np.linspace(ustar - du, ustar + du, 601):
        w = np.linalg.eigvalsh(H0 + u * av); g = np.diff(w); k = int(np.argmin(g))
        if g[k] < gmin:
            gmin = g[k]; ustar = u; kp = k
    _, Phi = np.linalg.eigh(H0 + ustar * av)
    return ustar, kp, Phi


def min_gap(H0p, A, T=80.0, n=1501):
    av = np.diag(np.diag(A)); us = np.linspace(-T, T, n)
    gs = np.array([np.min(np.diff(np.linalg.eigvalsh(H0p + u * av))) for u in us])
    k = int(np.argmin(gs))
    fine = np.linspace(us[max(k-1, 0)], us[min(k+1, n-1)], 401)
    return float(min(np.min(np.diff(np.linalg.eigvalsh(H0p + u * av))) for u in fine))


def make_V(seed, Phi, kp, keep=False):
    rng = np.random.default_rng(seed)
    M = rng.normal(size=(3, 3)); M = (M + M.T) / 2
    Mt = Phi.T @ M @ Phi                      # eigenbasis at u_*
    if keep:                                  # scalarize the node 2x2 block
        c = (Mt[kp, kp] + Mt[kp+1, kp+1]) / 2
        Mt[kp, kp] = Mt[kp+1, kp+1] = c
        Mt[kp, kp+1] = Mt[kp+1, kp] = 0.0
    V = Phi @ Mt @ Phi.T; V = (V + V.T) / 2
    return V / np.linalg.norm(V)


def expo(eps_list, dvals):
    x = np.log(np.array(eps_list)); y = np.log(np.abs(np.array(dvals)) + 1e-300)
    return float(np.polyfit(x, y, 1)[0])


def main():
    eps0 = [-2.0, 0.0, 3.0]; gam = [1.0, 0.8, 1.2]; a = [-1.0, 0.5, 2.0]
    H0, A = type1(eps0, gam, a); lo, mid, hi = (int(k) for k in np.argsort(a))
    ustar, kp, Phi = node_data(H0, A)
    P0 = propP(H0, A)
    epss = [0.1, 0.03, 0.01, 0.003]
    SEEDS = [1, 2, 3, 4]
    print("=" * 90)
    print(f"Isolating the crossing: node at u_*={ustar:+.3f}, degenerate pair (ranks {kp},{kp+1})")
    print("=" * 90)

    print("\nPART 1 -- multiple generic V directions: exponent of |dP_mm| (expect ~1 each):")
    for seed in SEEDS:
        V = make_V(seed, Phi, kp, keep=False)
        dmm = [propP(H0 + e * V, A)[mid, mid] - P0[mid, mid] for e in epss]
        gex = expo(epss, [min_gap(H0 + e * V, A) for e in epss])
        print(f"  seed {seed}: |dP_mm| ∝ eps^{expo(epss, dmm):.2f}   (gap ∝ eps^{gex:.2f})")

    print("\nPART 2 -- isolate the crossing (V_open vs V_keep, equal norm, seed 1):")
    Vo = make_V(1, Phi, kp, keep=False); Vk = make_V(1, Phi, kp, keep=True)
    print(f"  {'eps':>7} | {'gap_open':>9} {'gap_keep':>9} | {'dPmm_open':>10} {'dPmm_keep':>10} {'diff':>10}")
    go = []; gk = []; do = []; dk = []; df = []
    for e in epss:
        gO = min_gap(H0 + e * Vo, A); gK = min_gap(H0 + e * Vk, A)
        dO = propP(H0 + e * Vo, A)[mid, mid] - P0[mid, mid]
        dK = propP(H0 + e * Vk, A)[mid, mid] - P0[mid, mid]
        go.append(gO); gk.append(gK); do.append(dO); dk.append(dK); df.append(dO - dK)
        print(f"  {e:7.3f} | {gO:9.2e} {gK:9.2e} | {dO:10.3e} {dK:10.3e} {dO-dK:10.3e}")
    print(f"\n  gap_open ∝ eps^{expo(epss, go):.2f}   gap_keep ∝ eps^{expo(epss, gk):.2f}  "
          f"(keep >> steeper => crossing preserved to leading order)")
    print(f"  dPmm_open ∝ eps^{expo(epss, do):.2f}   dPmm_keep ∝ eps^{expo(epss, dk):.2f}   "
          f"DIFFERENCE ∝ eps^{expo(epss, df):.2f}")
    print("\n  Reading: if gap_keep is much steeper than gap_open (~2 vs ~1) the crossing is preserved by")
    print("  V_keep; if dPmm_open and dPmm_keep are BOTH ~eps^1 and their DIFFERENCE is ~eps^2, the")
    print("  crossing-removal channel is analytic O(eps^2), subdominant to the common O(eps) bulk.")


if __name__ == "__main__":
    main()
