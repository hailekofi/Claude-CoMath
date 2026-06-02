"""
skeleton_two_transcendentals.py  --  the "elementary BE skeleton + two transcendentals" structure.

Establishes, for the Type-1 N=3 transition matrix P (doubly stochastic, P[j,k]=prob k->j,
slope order lo,mid,hi = argsort(a)):

  R15  BE + double-stochasticity determine P up to EXACTLY TWO real numbers:
       P_mm := P[mid,mid] (middle survival) and b := P[hi,lo] (lo->hi, a chirality).
       Every entry is an explicit affine function of {P_ll^BE, P_hh^BE, P_mm, b}.
       P is strongly NON-symmetric (||P-P^T|| ~ 0.6-1.0), so b is independent of P_mm
       (no symmetry collapse to a single unknown).

  R16  Regime decomposition vs the incoherent (BE + crossing geometry, NO interference)
       model: P_mm is the UNIVERSAL hard number (incoherent off by up to ~0.99 in deep
       overlap), while b is a SOFT transcendental -- elementary (cyclic) in BOTH the
       separated and the deep-overlap limits, only mildly interference-dressed in between.
       Transcendence localizes on P_mm because the middle level is the only one that
       RECOMBINES (it passes both the mid-lo and mid-hi crossings).

Reproducible: python3 skeleton_two_transcendentals.py
Gold context: the quantum P here uses a single adaptive high-order solve (DOP853, T=90);
values agree with oracle.py / num_S12 to ~1e-3 (enough for the structural statements;
the gold anchors are P_mm canonical 0.214724, sampleB 0.021018).
"""
from __future__ import annotations
import numpy as np
from scipy.integrate import solve_ivp


def type1(eps, gam, a):
    eps = np.array(eps, float); gam = np.array(gam, float); a = np.array(a, float); g2 = gam**2
    H0 = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            if i != j:
                H0[i, j] = gam[i]*gam[j]*(a[i]-a[j])/(eps[i]-eps[j])
        H0[i, i] = -sum(g2[k]*(a[i]-a[k])/(eps[i]-eps[k]) for k in range(3) if k != i)
    return H0, np.diag(a)


def quantumP(eps, gam, a, T=90.0):
    """P[j,k] = |U_jk|^2, U = U(+T,-T) in the diabatic (standard) basis; columns = initial."""
    H0, A = type1(eps, gam, a)
    sol = solve_ivp(lambda u, y: (-1j*(H0 + u*A) @ y.reshape(3, 3)).ravel(),
                    [-T, T], np.eye(3, dtype=complex).ravel(),
                    rtol=1e-12, atol=1e-13, method="DOP853")
    return np.abs(sol.y[:, -1].reshape(3, 3))**2


def delta(eps, gam, a, i, j):
    s = gam[i]*gam[j]/(eps[i]-eps[j])
    return s*s*abs(a[i]-a[j])                       # BE window action s_ij^2 |a_i-a_j|


def incoherentP(eps, gam, a):
    """BE + crossing geometry, NO interference: independent sequential 2-level crossings,
    composed in time order.  The 'classical skeleton' baseline."""
    eps = np.array(eps, float); gam = np.array(gam, float); a = np.array(a, float)
    H0, _ = type1(eps, gam, a)
    cross = []
    for i in range(3):
        for j in range(i+1, 3):
            uij = (H0[j, j]-H0[i, i])/(a[i]-a[j])      # diabatic crossing time
            p = np.exp(-2*np.pi*delta(eps, gam, a, i, j))
            cross.append((uij, i, j, p))
    cross.sort(key=lambda c: c[0])                     # time order
    T = np.eye(3)
    for _, i, j, p in cross:
        M = np.eye(3); M[i, i] = p; M[j, j] = p; M[i, j] = 1-p; M[j, i] = 1-p
        T = M @ T
    return T


def reconstruct_from_two(Pll, Phh, Pmm, b, lo, mid, hi):
    """The affine reconstruction: all 9 entries from {P_ll^BE, P_hh^BE, P_mm, b}, via
    double-stochasticity.  Returns the slope-ordered 3x3 (rows=final, cols=initial)."""
    P = np.zeros((3, 3))
    P[lo, lo] = Pll; P[hi, hi] = Phh; P[mid, mid] = Pmm; P[hi, lo] = b
    P[mid, lo] = 1 - Pll - b
    P[hi, mid] = 1 - b - Phh
    P[mid, hi] = Pll + b - Pmm
    P[lo, hi] = 1 - (Pll + b - Pmm) - Phh
    P[lo, mid] = 1 - Pll - P[lo, hi]
    return P


def main():
    print("="*84)
    print("SKELETON + TWO TRANSCENDENTALS  (Type-1 N=3 transition matrix structure)")
    print("="*84)

    samples = {
        "canonical": ([-2, 0, 3], [1, 0.8, 1.2], [-1, 0.5, 2]),
        "generic1": ([-1.4, 0.2, 1.7], [1.1, -0.7, 0.9], [-1.2, 0.6, 1.5]),
        "generic2(non-mono a)": ([-1, 0, 1], [1, 1, 1], [0.5, -1.0, 1.3]),
    }

    print("\n--- R15: BE fixes 2 entries; matrix = affine in {P_mm, b}; P is NOT symmetric ---")
    for name, (eps, gam, a) in samples.items():
        eps = np.array(eps, float); gam = np.array(gam, float); a = np.array(a, float)
        lo, mid, hi = np.argsort(a)
        P = quantumP(eps, gam, a)
        Pll_be = np.exp(-2*np.pi*(delta(eps, gam, a, lo, mid)+delta(eps, gam, a, lo, hi)))
        Phh_be = np.exp(-2*np.pi*(delta(eps, gam, a, hi, lo)+delta(eps, gam, a, hi, mid)))
        be_err = max(abs(P[lo, lo]-Pll_be), abs(P[hi, hi]-Phh_be))
        Prec = reconstruct_from_two(P[lo, lo], P[hi, hi], P[mid, mid], P[hi, lo], lo, mid, hi)
        rec_err = np.max(np.abs(P - Prec))
        asym = np.max(np.abs(P - P.T))
        print(f"  [{name:20s}] BE fixes P_ll,P_hh (err {be_err:.1e}); "
              f"reconstruct-from-2 err {rec_err:.1e}; ||P-P^T||={asym:.2f} (NOT symmetric)")

    print("\n--- R16: P_mm is the hard number; b is a SOFT transcendental (incoherent vs quantum) ---")
    regimes = {
        "separated  d~0.006": ([-3, 0, 3], [0.45, 0.4, 0.5], [-1.3, 0.2, 1.6]),
        "moderate   d~0.24 (canonical)": ([-2, 0, 3], [1, 0.8, 1.2], [-1, 0.5, 2]),
        "deep overlap d~5.7": ([-1, 0, 1], [1.5, 1.4, 1.5], [-1, 0.3, 1.4]),
    }
    print(f"  {'regime':32s} {'P_mm:Q':>8} {'P_mm:incoh':>10} {'|dP_mm|':>8} | "
          f"{'b:Q':>7} {'b:incoh':>8} {'|db|':>7}")
    for name, (eps, gam, a) in regimes.items():
        eps = np.array(eps, float); gam = np.array(gam, float); a = np.array(a, float)
        lo, mid, hi = np.argsort(a)
        Q = quantumP(eps, gam, a); I = incoherentP(eps, gam, a)
        pm_q, pm_i = Q[mid, mid], I[mid, mid]
        b_q, b_i = Q[hi, lo], I[hi, lo]
        print(f"  {name:32s} {pm_q:8.4f} {pm_i:10.4f} {abs(pm_q-pm_i):8.4f} | "
              f"{b_q:7.4f} {b_i:8.4f} {abs(b_q-b_i):7.4f}")
    print("\n  => b collapses to its incoherent (cyclic) value in BOTH limits (|db|~1e-3);")
    print("     P_mm is off by up to ~0.99 incoherently in deep overlap -- the universal hard one.")


if __name__ == "__main__":
    main()
