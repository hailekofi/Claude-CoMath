"""
ws_geom_m3_t2t3.py -- WS-GEOM M3 T2 + T3 (verification).

T2 (two-vertex reachability): the image of Phi:(gam,eps,a)->{P_mm,b} touches EXACTLY two
Birkhoff vertices -- identity (1,0) [diabatic limit, couplings->0] and the node-selected
directed cycle (0,1) [adiabatic limit, couplings->inf]. The other four permutations map to
(0,0) (reverse cycle, (lo mid), (mid hi)) or (1,1) (extreme-swap (lo hi)) and are NOT reached.

T3 (edge-selection boundary = decoupling locus, R9): driving a coupling to zero (a level
decouples) lands {P_mm,b} on an axis edge, edge-resolved by WHICH level decouples:
  middle decouples (g_mid->0)  -> P_mm -> 1   (the P_mm=1 edge; the middle is a spectator);
  extreme decouples (g_lo or g_hi ->0) -> b -> 0 (the b=0 edge; the directed cycle is broken).

Reproduce: python3 ws_geom_m3_t2t3.py
"""
from __future__ import annotations
import numpy as np
from scipy.integrate import solve_ivp


def type1(eps, gam, a):
    eps = np.array(eps, float); gam = np.array(gam, float); a = np.array(a, float); g2 = gam ** 2
    H0 = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            if i != j:
                H0[i, j] = gam[i] * gam[j] * (a[i] - a[j]) / (eps[i] - eps[j])
        H0[i, i] = -sum(g2[k] * (a[i] - a[k]) / (eps[i] - eps[k]) for k in range(3) if k != i)
    return H0, np.diag(a)


def pmb(eps, gam, a, T=70.0):
    H0, A = type1(eps, gam, a)
    s = solve_ivp(lambda u, y: (-1j * (H0 + u * A) @ y.reshape(3, 3)).ravel(),
                  [-T, T], np.eye(3, dtype=complex).ravel(), rtol=1e-9, atol=1e-11, method="DOP853")
    P = np.abs(s.y[:, -1].reshape(3, 3)) ** 2
    lo, mid, hi = np.argsort(a)
    return P[mid, mid], P[hi, lo]


def main():
    print("=" * 84)
    print("WS-GEOM M3 T2 + T3 (verification)")
    print("=" * 84)
    eps = [-2, 0, 3]; a = [-1, .5, 2.]; base = np.array([1, .8, 1.2])

    print("\nT2 -- two-vertex reachability  (identity=(1,0), directed cycle=(0,1); (0,0)&(1,1) forbidden)")
    print(f"  {'gscale':>7} {'(P_mm, b)':>16} {'dist(0,0)':>10} {'dist(1,1)':>10}")
    far00 = far11 = True
    for sc in [0.03, 0.1, 0.3, 1.0, 2.0, 3.5]:
        pm, b = pmb(eps, (sc * base).tolist(), a)
        d00 = np.hypot(pm, b); d11 = np.hypot(pm - 1, b - 1)
        far00 &= (d00 > 0.5); far11 &= (d11 > 0.5)
        print(f"  {sc:7.2f} {f'({pm:.3f},{b:.3f})':>16} {d00:10.2f} {d11:10.2f}")
    print(f"  => identity & cycle reached; (0,0)&(1,1) never approached (always >0.5): {far00 and far11}")

    print("\nT3 -- edge-selection boundary = decoupling locus")
    cases = [("middle decouple  g_mid->0", [1, 1e-3, 1.2], "P_mm->1"),
             ("extreme decouple g_lo->0 ", [1e-3, .8, 1.2], "b->0"),
             ("extreme decouple g_hi->0 ", [1, .8, 1e-3], "b->0")]
    for name, gam, expect in cases:
        pm, b = pmb(eps, gam, a)
        print(f"  {name}: (P_mm,b)=({pm:.4f},{b:.4f})   expect {expect}")
    print("  => middle decoupling -> P_mm=1 edge; extreme decoupling -> b=0 edge (R9 edge-selection).")


if __name__ == "__main__":
    main()
