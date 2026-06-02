"""
deformation_family_probe.py  --  the commuting family as an isomonodromic deformation.

Tests the user's idea: map the 2-parameter commuting family (excluding identity) to the
two transcendentals {P_mm, b} of the transition matrix.

FINDINGS (coordinator-run 2026-06-02):
  F1 [established, machine prec]  Different slope vectors a (same gamma,eps) COMMUTE and SHARE
     one a-INDEPENDENT eigenbasis phi_i(u): ||[H^a,H^a']||~1e-16, shared-eigvec defect ~1e-15.
     => the derivative-coupling seed  W_ij = <phi_i|phi_j'>  is a-INDEPENDENT: the ENTIRE
     transcendental content of the WHOLE 2-parameter family is ONE geometric object W.
  F2 [numerically-supported]  Along a path a(t) in the family, the two transcendentals
     {P_mm(a), b(a)} vary SMOOTHLY (and so do the elementary BE actions delta_ij).  Each
     member's {P_mm,b} = holonomy of the FIXED W against that member's a-linear phases.

INTERPRETATION (honest):
  - The a-INVARIANT is the SEED W (the connection), NOT the two transcendentals (which vary).
    This is the isomonodromy structure: W = monodromy/Stokes seed; a = deformation times;
    {P_mm(a),b(a)} = tau-data flowing along the deformation (a candidate Schlesinger/Garnier flow).
  - So the family ORGANIZES the transcendentals (one shared rigid seed) but does NOT COMPUTE
    them (Abelian ceiling, R1): no a-flow produces W; the flow only transports it.
  - No free deep-overlap shortcut: the holonomy of W with deep-overlap phases is still the hard
    connection problem.  The genuine lever is the (frontier) a-flow / Schlesinger derivation with
    the separated-regime elementary limit as initial condition.

Reproduce: python3 deformation_family_probe.py
"""
from __future__ import annotations
import numpy as np
from itertools import permutations
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


def H_of(eps, gam, a, u):
    H0, A = type1(eps, gam, a)
    return H0 + u*A


def quantumP(eps, gam, a, T=90.0):
    H0, A = type1(eps, gam, a)
    sol = solve_ivp(lambda u, y: (-1j*(H0 + u*A) @ y.reshape(3, 3)).ravel(),
                    [-T, T], np.eye(3, dtype=complex).ravel(),
                    rtol=1e-11, atol=1e-12, method="DOP853")
    return np.abs(sol.y[:, -1].reshape(3, 3))**2


def delta(eps, gam, a, i, j):
    s = gam[i]*gam[j]/(eps[i]-eps[j]); return s*s*abs(a[i]-a[j])


def foundation(eps, gam, a_list):
    print("F1: family members (same gamma,eps; different slope a) -- commute & share eigenbasis?")
    for u in [-2.0, -0.3, 1.0, 3.0]:
        Hs = [H_of(eps, gam, a, u) for a in a_list]
        c = max(np.max(np.abs(Hs[0] @ Hk - Hk @ Hs[0])) for Hk in Hs[1:])
        _, V0 = np.linalg.eigh(Hs[0]); _, V1 = np.linalg.eigh(Hs[1])
        ov = np.abs(V0.conj().T @ V1)
        defect = min(np.max(np.abs(ov[list(p), :] - np.eye(3))) for p in permutations(range(3)))
        print(f"   u={u:+.1f}  max||[H0,Hk]||={c:.1e}  shared-eigvec defect={defect:.1e}")


def deformation(eps, gam, a_sep, a_deep, n=7):
    print("\nF2: deformation path a(t)=(1-t)a_sep+t a_deep in the family -- transcendentals vary smoothly")
    print(f"   {'t':>4} {'d_lo,mid':>9} {'d_mid,hi':>9} | {'P_mm':>8} {'b=lo->hi':>9}")
    a_sep = np.array(a_sep, float); a_deep = np.array(a_deep, float)
    for t in np.linspace(0, 1, n):
        a = (1-t)*a_sep + t*a_deep
        lo, mid, hi = np.argsort(a); P = quantumP(eps, gam, a)
        print(f"   {t:4.2f} {delta(eps,gam,a,lo,mid):9.3f} {delta(eps,gam,a,mid,hi):9.3f} | "
              f"{P[mid,mid]:8.4f} {P[hi,lo]:9.4f}")


def main():
    print("="*78)
    print("COMMUTING FAMILY AS ISOMONODROMIC DEFORMATION  (Type-1 N=3)")
    print("="*78)
    eps, gam = [-2, 0, 3], [1, 0.8, 1.2]
    rng = np.random.default_rng(0)
    foundation(eps, gam, [[-1, 0.5, 2], [0.3, -0.8, 1.4], list(rng.uniform(-1.5, 1.5, 3))])
    deformation(eps, gam, [-0.30, 0.10, 0.42], [-1.5, 0.45, 2.1])
    print("\n=> ONE a-independent seed W powers the whole family; {P_mm,b} are its a-varying")
    print("   holonomy (the isomonodromy/Schlesinger deformation data). Seed invariant, NOT the")
    print("   two numbers; no free deep-overlap shortcut (Abelian ceiling).")


if __name__ == "__main__":
    main()
