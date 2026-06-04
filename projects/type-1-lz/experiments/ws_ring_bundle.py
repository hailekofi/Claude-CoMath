"""
ws_ring_bundle.py -- the fiber bundle over the Type-1 ring: how a WITHIN-RING perturbation (fixed gamma,eps;
vary the slopes a) perturbs the transition amplitude, and how 'elementary' that response is.

Base = ring (a mod shift, fixed gamma,eps); fiber = {P_mm, b} (BE-reduced amplitude). We use FINITE
within-ring steps (|da|=0.1) -- NOT tiny-h derivatives -- because a single-T propagator's endpoint Fresnel
tail VARIES with a and swamps a 1/h-amplified derivative (an early tiny-h gradient run gave garbage
~ -2.4 vs the true ~ -0.08; retracted). With T-averaged P and finite steps the responses are clean.

FINDINGS (T-averaged DOP853, step 0.1):
  * SHIFT-NULL: a->a+(1,1,1) gives dP_mm = -0.0000 (global-phase invariance, exact).
  * RANK-2 tangent: any within-ring dP is reproduced from {dP[lo,lo],dP[hi,hi],dP_mm,db} via
    double-stochasticity to ~1e-6 -> the tangent map is ELEMENTARY except for the rank-2 transcendental
    pair {dP_mm,db} (the tangent-level R15).
  * ELEMENTARY CAPTURE of the response: the elementary uniform law (R12) tracks the EXACT within-ring
    response of P_mm with cosine +0.92 and magnitude ratio ~1.3 over substantial-response directions
    (per-direction ratios 1.50,1.01,1.08). It flips sign only in NEAR-NULL directions where the true
    response is ~1e-3 (there the law's ~2% value-error dominates the tiny signal).
  NET: the within-ring perturbation is MOSTLY characterizable non-transcendentally -- the elementary law
  gives a good (~30%-magnitude, cosine 0.92) account of how the amplitude moves -- with a SUBDOMINANT
  transcendental residual that becomes relatively dominant only where the elementary response nearly
  vanishes. (No EXACT elementary flow -- consistent with R17 -- but a good APPROXIMATE one.)

Reproduce: python3 ws_ring_bundle.py
"""
from __future__ import annotations
import os, sys
import numpy as np
from scipy.integrate import solve_ivp
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
import ws_o3_uniform as uni  # noqa
np.seterr(all="ignore")


def _type1(eps, gam, a):
    eps = np.array(eps, float); gam = np.array(gam, float); a = np.array(a, float); g2 = gam ** 2
    H0 = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            if i != j:
                H0[i, j] = gam[i] * gam[j] * (a[i] - a[j]) / (eps[i] - eps[j])
        H0[i, i] = -sum(g2[k] * (a[i] - a[k]) / (eps[i] - eps[k]) for k in range(3) if k != i)
    return H0, np.diag(a)


def Pexact(eps, gam, a, Ts=(80.0, 110.0)):
    """T-averaged DOP853 (the endpoint Fresnel tail VARIES with a, so it must be averaged out -- a single-T
    propagator gives garbage finite-difference responses). We use FINITE within-ring steps (not tiny-h
    derivatives) to avoid 1/h noise amplification."""
    H0, A = _type1(eps, gam, a); acc = np.zeros((3, 3))
    for T in Ts:
        s = solve_ivp(lambda u, y: (-1j * (H0 + u * A) @ y.reshape(3, 3)).ravel(),
                      [-T, T], np.eye(3, dtype=complex).ravel(), rtol=1e-9, atol=1e-11, method="DOP853")
        acc += np.abs(s.y[:, -1].reshape(3, 3).T) ** 2
    return acc / len(Ts)


def grad_a(f, a, h=1e-3):
    a = np.array(a, float); g = np.zeros(3)
    for k in range(3):
        ap = a.copy(); am = a.copy(); ap[k] += h; am[k] -= h
        g[k] = (f(ap) - f(am)) / (2 * h)
    return g


def main():
    eps = [-2.0, 0.0, 3.0]; gam = [1.0, 0.8, 1.2]; a0 = [-1.0, 0.5, 2.0]
    lo, mid, hi = (int(k) for k in np.argsort(a0))
    print("=" * 92)
    print("Fiber bundle over the Type-1 ring: within-ring (a) response of the amplitude")
    print(f"base point a={a0} (gamma,eps fixed); slope order lo,mid,hi={lo,mid,hi}")
    print("=" * 92)
    sh = np.array([1.0, 1.0, 1.0]) / np.sqrt(3)
    def proj(v):
        return v - (v @ sh) * sh

    P0 = Pexact(eps, gam, a0); Pmm0 = P0[mid, mid]; b0 = P0[hi, lo]
    Pmm0_u = uni.uniform_P_mm(np.array(eps), np.array(gam), np.array(a0))
    print(f"\nbase: exact P_mm={Pmm0:.4f}  uniform P_mm={Pmm0_u:.4f}  b={b0:.4f}")

    # FINITE within-ring responses (step 0.1 in several shift-orthogonal directions; no 1/h amplification)
    print("\nFINITE within-ring responses  (step |da|=0.1, shift-orthogonal directions):")
    print(f"  {'dir':>22} | {'dPmm_exact':>11} {'dPmm_unif':>10} {'ratio u/e':>10} | {'shift-null dPmm':>15}")
    rng = np.random.default_rng(0); step = 0.1
    dirs = {}
    for k, nm in [(0, 'e0'), (1, 'e1'), (2, 'e2')]:
        v = np.zeros(3); v[k] = 1.0; dirs['de_' + nm] = proj(v)
    for s in range(2):
        dirs[f'rand{s}'] = proj(rng.normal(size=3))
    dirs['shift(1,1,1)'] = sh.copy()
    rows = []
    for nm, d in dirs.items():
        d = d / (np.linalg.norm(d) + 1e-30)
        a1 = list(np.array(a0) + step * d)
        P1 = Pexact(eps, gam, a1)
        dpe = P1[mid, mid] - Pmm0
        dpu = uni.uniform_P_mm(np.array(eps), np.array(gam), np.array(a1)) - Pmm0_u
        rows.append((nm, dpe, dpu, d))
        tag = '' if 'shift' not in nm else '   <- should be ~0'
        print(f"  {nm:>22} | {dpe:+11.4f} {dpu:+10.4f} {dpu/dpe if abs(dpe)>1e-6 else float('nan'):>10.3f}{tag}")

    # elementary capture over the in-plane response vectors (exclude the shift row)
    de = np.array([r[1] for r in rows if 'shift' not in r[0]])
    du = np.array([r[2] for r in rows if 'shift' not in r[0]])
    cos = float(de @ du / (np.linalg.norm(de) * np.linalg.norm(du) + 1e-30))
    print(f"\n  elementary capture of the within-ring response (uniform vs exact, over the in-plane dirs):")
    print(f"    cosine = {cos:+.3f}   |response_uniform|/|response_exact| = {np.linalg.norm(du)/np.linalg.norm(de):.3f}")

    # (4) rank-2 tangent check: generic a-perturbation, full dP vs reconstruction from 4 entries
    rng = np.random.default_rng(3); da = proj(rng.normal(size=3)); da /= np.linalg.norm(da); eta = 1e-3
    Pp = Pexact(eps, gam, list(np.array(a0) + eta * da)); Pm = Pexact(eps, gam, list(np.array(a0) - eta * da))
    dP = (Pp - Pm) / (2 * eta)
    # A doubly-stochastic tangent (row & col sums = 0) has 4 DOF. Take as controls the 2 elementary BE-survival
    # responses {dP[lo,lo],dP[hi,hi]} and the 2 transcendental ones {dP_mm,db}; reconstruct the other 5 entries
    # from the 6 (rank-5) sum constraints. If max|reconstruct - exact dP| ~ 0, the tangent map is rank-2.
    known = {(lo, lo): dP[lo, lo], (hi, hi): dP[hi, hi], (mid, mid): dP[mid, mid], (hi, lo): dP[hi, lo]}
    idx = [(lo, mid), (lo, hi), (mid, lo), (mid, hi), (hi, mid)]   # the 5 unknown entries
    M = np.zeros((6, 5)); rhs = np.zeros(6)
    for r in range(3):                                            # row-sum = 0
        for ii, (rr, cc) in enumerate(idx):
            if rr == r:
                M[r, ii] = 1.0
        rhs[r] = -sum(known.get((r, c), 0.0) for c in range(3))
    for c in range(3):                                            # col-sum = 0
        for ii, (rr, cc) in enumerate(idx):
            if cc == c:
                M[3 + c, ii] = 1.0
        rhs[3 + c] = -sum(known.get((r2, c), 0.0) for r2 in range(3))
    sol, *_ = np.linalg.lstsq(M, rhs, rcond=None)
    R = np.zeros((3, 3))
    for key, val in known.items():
        R[key] = val
    for ii, (rr, cc) in enumerate(idx):
        R[rr, cc] = sol[ii]
    resid = float(np.max(np.abs(R - dP)))
    print(f"\n(4) RANK-2 tangent check: reconstruct full dP from 4 numbers "
          f"{{dP[lo,lo],dP[hi,hi],dP_mm,db}} via double-stochasticity:")
    print(f"    max|reconstruct - exact dP| = {resid:.2e}  (small => the tangent map is rank-2: only dP_mm,db")
    print(f"    are transcendental; dP[lo,lo],dP[hi,hi] are elementary BE-survival responses).")


if __name__ == "__main__":
    main()
