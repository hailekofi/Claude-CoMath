"""
ws_geom_m3_t4.py -- WS-GEOM M3 T4 (the spinor / cycle-orientation Z2 label).

OUTCOME (T4 RESOLVED): the directed-cycle ORIENTATION is a genuine Z2 topological label
-- BOTH orientations occur across parameter space -- and its SELECTOR is now derived
(see ws_geom_m3_t4_selector.py).

  * Z2 ESTABLISHED: the deterministic, dynamics-free overlap-continuation of the order-0 eigenframe
    (T1's geometric construction) yields BOTH orientations over broad samples:
        forward  (1,2,0) lo->hi->mid->lo   <=>  node swaps slope-pair (mid,hi)
        reverse  (2,0,1) lo->mid->hi->lo   <=>  node swaps slope-pair (lo,mid)
    [(lo hi) o (adjacent pair) = the two opposite 3-cycles -- group-theory half, derived in T1.]
    This corrects T1 Step-4's "uniform orientation" (a 5-sample artifact). The Z2 is the eigenframe
    spinor/double-cover sector (delta_j = +-1).

  * SELECTOR RESOLVED (ws_geom_m3_t4_selector.py): orientation = sign(tr H(u_*) - 3 E_*), the
    energy-position of the REAL node (R3) -- LOWER pair degenerate => FWD, UPPER pair => REV. Matches this
    Z2 on 154/154 clean samples (two independent algorithms agree). The earlier "orientation = -sign(u_*)"
    was REFUTED (~92% on one seed, 33% on another); the local predictor screen failed because the selector
    is the GLOBAL energy-position of the protected real crossing, NOT a local sign and NOT the nearest
    complex branch point.

Reproduce: python3 ws_geom_m3_t4.py   (deterministic continuation; no oracle calls)
"""
from __future__ import annotations
import numpy as np
from collections import Counter


def type1(eps, gam, a):
    eps = np.array(eps, float); gam = np.array(gam, float); a = np.array(a, float); g2 = gam ** 2
    H0 = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            if i != j:
                H0[i, j] = gam[i] * gam[j] * (a[i] - a[j]) / (eps[i] - eps[j])
        H0[i, i] = -sum(g2[k] * (a[i] - a[k]) / (eps[i] - eps[k]) for k in range(3) if k != i)
    return H0, np.diag(a)


def cont_orient(eps, gam, a, T=80.0, n=120001):
    """Deterministic overlap-continuation of the order-0 eigenframe through u (incl. the node);
    returns the permutation in slope-rank labels."""
    H0, A = type1(eps, gam, a)
    us = np.linspace(-T, T, n)
    _, V = np.linalg.eigh(H0 + us[0] * A); Vs = V.copy(); Vp = V.copy()
    for u in us[1:]:
        _, V = np.linalg.eigh(H0 + u * A)
        M = np.abs(Vp.conj().T @ V); V = V[:, np.argsort(np.argmax(M, axis=0))]
        for k in range(3):
            if np.vdot(Vp[:, k], V[:, k]).real < 0:
                V[:, k] *= -1
        Vp = V
    sig = [int(np.argmax(np.abs(Vs[:, k]))) for k in range(3)]
    pi = [int(np.argmax(np.abs(Vp[:, k]))) for k in range(3)]
    P = np.zeros((3, 3))
    for k in range(3):
        P[pi[k], sig[k]] = 1
    order = np.argsort(a); Ps = P[np.ix_(order, order)]
    return tuple(int(np.argmax(Ps[r])) for r in range(3))


def label(p):
    return "FWD lo->hi->mid->lo" if p == (1, 2, 0) else ("REV lo->mid->hi->lo" if p == (2, 0, 1) else f"other{p}")


def main():
    print("=" * 80)
    print("WS-GEOM M3 T4 -- the cycle-orientation Z2 (spinor/double-cover) label")
    print("=" * 80)
    tally = Counter()
    for seed in range(40):
        rng = np.random.default_rng(100 + seed)
        eps = sorted(rng.uniform(-2, 2, 3))
        if min(np.diff(eps)) < 0.3:
            continue
        gam = list(rng.uniform(-1.5, 1.5, 3)); a = list(rng.uniform(-2, 2, 3))
        if min(np.abs(np.diff(sorted(a)))) < 0.3 or min(abs(x) for x in a) < 0.2:
            continue
        tally[cont_orient(eps, gam, a)] += 1
    print("Deterministic overlap-continuation orientation over broad samples (slope labels):")
    for p, c in tally.most_common():
        print(f"  {p}: {c:2d}   {label(p)}")
    both = (1, 2, 0) in tally and (2, 0, 1) in tally
    print(f"\nBOTH orientations occur => Z2 label is REAL: {both}")
    print("Selector RESOLVED in ws_geom_m3_t4_selector.py: orientation = sign(tr H(u_*) - 3 E_*),")
    print("the energy-position of the real node (LOWER pair => FWD, UPPER pair => REV); 154/154 clean.")


if __name__ == "__main__":
    main()
