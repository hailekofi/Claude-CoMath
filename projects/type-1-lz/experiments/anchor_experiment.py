"""
Day 1-2 anchor experiment for Type-1 N=3 Landau-Zener (standalone, reproducible).

Two parts:
  (A) DICTIONARY   Type-1 {gam,eps,a}  ->  Chernyak-Sinitsyn (arXiv:2006.15144) Eq.(11)
                   frame, showing the C-S solvable bow-tie point is OUTSIDE Type-1.
  (B) EXPERIMENT   benchmark P by direct integration, self-calibrated on the EXACT
                   Brundobler-Elser extreme survivals, vs the incoherent-product
                   prediction for the (open) middle survival.

Depends only on numpy + scipy.  This is a cross-check replica of the project's
geometry.py formulas; the project's own interaction-picture harness is the 1e-9 oracle.
"""
import numpy as np
from scipy.integrate import solve_ivp


def type1(eps, gam, a):
    """Type-1 Cauchy Hamiltonian data:  H(u) = H0 + u*diag(a)."""
    eps = np.asarray(eps, float); gam = np.asarray(gam, float); a = np.asarray(a, float)
    g2 = gam ** 2
    H0 = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            if i != j:
                H0[i, j] = gam[i] * gam[j] * (a[i] - a[j]) / (eps[i] - eps[j])
        H0[i, i] = -sum(g2[k] * (a[i] - a[k]) / (eps[i] - eps[k]) for k in range(3) if k != i)
    return H0, np.diag(a)


def be_exp(eps, gam, a, i, j):
    """Pairwise Brundobler-Elser exponent be_exp_ij."""
    return gam[i] ** 2 * gam[j] ** 2 * abs(a[i] - a[j]) / (eps[i] - eps[j]) ** 2


# ---------------------------------------------------------------- (A) dictionary
def to_CS_frame(H0, a):
    """
    Map Type-1 H(u)=H0+u diag(a) to Chernyak-Sinitsyn Eq.(11): two pure-slope outer
    levels + one FLAT middle level at energy eps_CS, via a time-shift u->u+t0 and a
    diagonal gauge (a_mid*u + c)*I.  Returns the C-S parameters; couplings are the
    (gauge-invariant) off-diagonal entries of H0.
    """
    lo, mid, hi = np.argsort(a)
    t0 = (H0[hi, hi] - H0[lo, lo]) / (a[lo] - a[hi])
    c = a[lo] * t0 + H0[lo, lo]
    eps_CS = a[mid] * t0 + H0[mid, mid] - c
    return dict(order=(lo, mid, hi), eps_CS=eps_CS, b1=a[hi] - a[mid], b2=a[mid] - a[lo],
                g_lomid=H0[lo, mid], g_midhi=H0[mid, hi], g_OUTER=H0[lo, hi])


# --------------------------------------------------------------- (B) benchmark P
def benchmark_P(eps, gam, a, T=160.0, rtol=1e-12, atol=1e-13):
    """P[n<-m] = |U_nm|^2 from i U' = H(u) U, U(-T)=I.  (diabatic = adiabatic asympt.)"""
    H0, A = type1(eps, gam, a)
    rhs = lambda u, y: (-1j * ((H0 + u * A) @ y.reshape(3, 3))).ravel()
    sol = solve_ivp(rhs, [-T, T], np.eye(3, dtype=complex).ravel(),
                    rtol=rtol, atol=atol, method='DOP853')
    return np.abs(sol.y[:, -1].reshape(3, 3)) ** 2


if __name__ == "__main__":
    samples = {
        "canonical": ((-2., 0., 3.), (1., .8, 1.2), (-1., .5, 2.)),
        "sampleB":   ((-1., 0., 1.5), (.9, 1.1, .8), (-.7, .4, 1.3)),
    }
    for lbl, (eps, gam, a) in samples.items():
        eps = np.array(eps); gam = np.array(gam); a = np.array(a)
        H0, _ = type1(eps, gam, a)
        d = to_CS_frame(H0, a); lo, mid, hi = d["order"]
        P = benchmark_P(eps, gam, a)
        BE_lo = np.exp(-2 * np.pi * (be_exp(eps, gam, a, lo, mid) + be_exp(eps, gam, a, lo, hi)))
        BE_hi = np.exp(-2 * np.pi * (be_exp(eps, gam, a, hi, mid) + be_exp(eps, gam, a, hi, lo)))
        inc   = np.exp(-2 * np.pi * (be_exp(eps, gam, a, mid, lo) + be_exp(eps, gam, a, mid, hi)))
        print(f"\n[{lbl}]  C-S frame: eps_CS={d['eps_CS']:+.4f}, g_OUTER={d['g_OUTER']:+.4f} (never 0 in Type-1)")
        print(f"  calib  P_lo={P[lo,lo]:.6f} (BE {BE_lo:.6f})   P_hi={P[hi,hi]:.6f} (BE {BE_hi:.6f})")
        print(f"  OPEN   P_mid={P[mid,mid]:.6f}  incoherent={inc:.6f}  ratio={P[mid,mid]/inc:.2f}")
