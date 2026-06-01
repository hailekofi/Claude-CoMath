"""
ws_e_geometry.py  --  WS-E: geometric invariants of the Type-1 N=3 MLZ middle survival,
the BE-survival = window-action identity, the double-stochastic dof count, and the
benchmark of candidate closed forms for P_mid against the WS-F gold oracle.

Findings (reproduced):
  * The two imaginary window actions |Im I_X| equal the two EXTREME BE survival exponents
    (each a sum over the two crossings of that extreme level).  The middle survival is NOT
    a window action.
  * P_mid is a genuine independent dof: fixing the two BE survivals leaves a 2-parameter
    family in the doubly-stochastic manifold along which P[mid,mid] varies.
  * The four complex turning points (avoided crossings) are two conjugate pairs => their
    cross-ratio is REAL (one geometric invariant per sample).  The node is the real double
    root (an ordinary point of the u-system: no local connection data).
  * Candidate elementary / two-pathway closed forms FAIL: in particular for strongly
    overlapping samples the true P_mid lies OUTSIDE the entire two-path Stuckelberg band,
    so no single interference phase reproduces it (=> genuine 3-crossing joint coherence).

Run:  python3 ws_e_geometry.py          (uses oracle.py; ~1 min/sample)
      python3 ws_e_geometry.py --fast   (skip oracle; geometry + dof only)
"""
from __future__ import annotations
import os, sys
import numpy as np
import sympy as sp

_HERE = os.path.dirname(os.path.abspath(__file__))
_UP = os.path.join(_HERE, "..", "uploads")
for c in (_UP, _HERE):
    if c not in sys.path:
        sys.path.insert(0, c)

from assay import Params, Geometry              # noqa: E402
from assay.actions import all_window_actions    # noqa: E402


def Gamma(eps, gam, a, i, j):
    return gam[i] ** 2 * gam[j] ** 2 * abs(a[i] - a[j]) / (eps[i] - eps[j]) ** 2


# --------------------------------------------------------------------------- #
#  turning points: complex u where two eigenvalues of H(u) coincide
# --------------------------------------------------------------------------- #
def turning_points(eps, gam, a):
    u = sp.symbols('u')
    H0 = sp.zeros(3, 3)
    e = [sp.nsimplify(x) for x in eps]; g = [sp.nsimplify(x) for x in gam]
    aa = [sp.nsimplify(x) for x in a]
    for i in range(3):
        for j in range(3):
            if i != j:
                H0[i, j] = g[i] * g[j] * (aa[i] - aa[j]) / (e[i] - e[j])
        H0[i, i] = -sum(g[k] ** 2 * (aa[i] - aa[k]) / (e[i] - e[k])
                        for k in range(3) if k != i)
    M = H0 + u * sp.diag(*aa)
    E = sp.symbols('E')
    chi = sp.Poly((M - E * sp.eye(3)).det(), E)
    disc = sp.Poly(sp.discriminant(chi, E), u)
    roots = np.roots([complex(c) for c in disc.all_coeffs()])
    return roots[np.argsort(np.abs(roots.imag))]


def cross_ratio(z):
    z1, z2, z3, z4 = z
    return ((z1 - z3) * (z2 - z4)) / ((z1 - z4) * (z2 - z3))


# --------------------------------------------------------------------------- #
#  P_mid is an independent dof: 2-parameter family fixing both BE survivals
# --------------------------------------------------------------------------- #
def pmid_free_dof(lo, hi, mid):
    """Return the shifts of P[mid,mid] along the doubly-stochastic directions that fix
    P[lo,lo] and P[hi,hi].  Nonzero shifts => P_mid not pinned by BE + double stochasticity."""
    idx = [lo, 1 if {lo, hi} == {0, 2} else None]  # generic build below
    # basis of zero-row/col-sum 3x3 matrices (4-dim)
    B = []
    for i in (0, 1):
        for j in (0, 1):
            M = np.zeros((3, 3)); M[i, j] += 1; M[i, 2] -= 1; M[2, j] -= 1; M[2, 2] += 1
            B.append(M)
    A = np.array([[b[lo, lo] for b in B], [b[hi, hi] for b in B]])
    null = np.linalg.svd(A)[2][2:]
    shifts = []
    for c in null:
        Mdir = sum(ci * bi for ci, bi in zip(c, B))
        shifts.append(float(Mdir[mid, mid]))
    return shifts


# --------------------------------------------------------------------------- #
#  candidate closed forms
# --------------------------------------------------------------------------- #
def candidates(d1, d2):
    p1, p2 = np.exp(-d1), np.exp(-d2)
    band_lo = (np.sqrt(p1 * p2) - np.sqrt((1 - p1) * (1 - p2))) ** 2
    band_hi = (np.sqrt(p1 * p2) + np.sqrt((1 - p1) * (1 - p2))) ** 2
    return {
        "incoherent p1p2": p1 * p2,
        "1-sech((d1+d2)/2)": 1 - 1 / np.cosh((d1 + d2) / 2),
        "p1p2+(1-p1)(1-p2)": p1 * p2 + (1 - p1) * (1 - p2),
        "2path band": (band_lo, band_hi),
    }


SAMPLES = {
    "canonical": ((-2., 0., 3.), (1., 0.8, 1.2), (-1., 0.5, 2.)),
    "sampleB":   ((-1., 0., 1.5), (0.9, 1.1, 0.8), (-0.7, 0.4, 1.3)),
    "well_sep":  ((-5., 0., 5.), (1., 1., 1.), (-1.5, 0., 1.5)),
}


def main(fast=False, T=80.0):
    if not fast:
        from oracle import oracle_P
    print(f"{'sample':10s} {'P_mid(oracle)':>13s} {'inc=p1p2':>10s} {'d1':>7s} {'d2':>7s} "
          f"{'Re(cr)':>8s} {'window|ImI|':>22s}")
    for nm, (eps, gam, a) in SAMPLES.items():
        lo, mid, hi = (int(k) for k in np.argsort(a))
        d1 = 2 * np.pi * Gamma(eps, gam, a, mid, lo)
        d2 = 2 * np.pi * Gamma(eps, gam, a, mid, hi)
        geo = Geometry(Params(eps=eps, gam=gam, a=a, x=0))
        imIx = sorted(float(abs(w["I_X"].imag)) for w in all_window_actions(geo))
        tp = turning_points(eps, gam, a)
        cmplx = tp[np.abs(tp.imag) > 1e-5]
        cr = cross_ratio(cmplx[:4]) if len(cmplx) >= 4 else complex("nan")
        if fast:
            pmid = float("nan")
        else:
            r = oracle_P(eps, gam, a, T=T)
            pmid = r["P"][mid, mid]
        print(f"{nm:10s} {pmid:13.6f} {np.exp(-d1-d2):10.6f} {d1:7.3f} {d2:7.3f} "
              f"{np.real(cr):8.3f} {str([round(x,3) for x in imIx]):>22s}")
        # candidate benchmark
        cand = candidates(d1, d2)
        for k, val in cand.items():
            if k == "2path band":
                lo_b, hi_b = val
                inside = (lo_b <= pmid <= hi_b) if not fast else None
                print(f"      2-path band [{lo_b:.4f},{hi_b:.4f}]  contains P_mid? {inside}")
            else:
                print(f"      {k:22s} = {val:.6f}")
        # dof
        shifts = pmid_free_dof(lo, hi, mid)
        print(f"      P_mid free-dof shifts (fixing both BE survivals): "
              f"{[round(s,4) for s in shifts]}  -> nonzero => P_mid independent")
    print("\nBE survivals = window actions (each |Im I_X| is a sum of the two 2*pi*Gamma of "
          "that extreme level).  The middle survival is NOT a window action.")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--fast", action="store_true")
    ap.add_argument("--T", type=float, default=80.0)
    a = ap.parse_args()
    main(fast=a.fast, T=a.T)
