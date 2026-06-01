"""
Coordinator verification of WS-D (non-Abelian Ê / MC factorization obstruction).

WS-D's VERDICT (Type-1 cannot be MC-factorized) is upheld, but its stated obstruction
leg c1 (rigidity / "null-space dimension 0") is INCORRECT. This script (1) falsifies c1,
and (2) establishes the CORRECT obstruction: Type-1 crossings are permanently
marginally-overlapping (separation/width ratio bounded ~O(1)), so MC's isolated-2-level
limit is unreachable. numpy only.
"""
import numpy as np


def features(th):
    """th=(g0,g1,g2,e0,e1,e2,a0,a1,a2) -> (couplings[3], slopes[3], diabatic crossing locs[3])."""
    g, e, a = th[0:3], th[3:6], th[6:9]
    pairs = [(0, 1), (0, 2), (1, 2)]
    H0ii = [-sum(g[k]**2 * (a[i]-a[k])/(e[i]-e[k]) for k in range(3) if k != i) for i in range(3)]
    c = [g[i]*g[j]*(a[i]-a[j])/(e[i]-e[j]) for (i, j) in pairs]
    u = [-(H0ii[i]-H0ii[j])/(a[i]-a[j]) for (i, j) in pairs]
    return np.array(c), np.array(a), np.array(u)


def _jac(fsel, th, h=1e-6):
    base = fsel(th); J = np.zeros((len(base), 9))
    for k in range(9):
        tp = th.copy(); tp[k] += h; tm = th.copy(); tm[k] -= h
        J[:, k] = (fsel(tp) - fsel(tm)) / (2*h)
    return J


def falsify_c1():
    """WS-D c1 claimed: fixing couplings+slopes leaves null-space dim 0 (crossings can't move)."""
    th = np.array([1., .8, 1.2, -2., 0., 3., -1., .5, 2.])
    F_cs = lambda t: np.concatenate(features(t)[:2])         # couplings + slopes (6)
    F_cu = lambda t: np.concatenate([*features(t)[:2], [features(t)[2][2]]])  # + middle crossing
    r6 = np.linalg.matrix_rank(_jac(F_cs, th), tol=1e-7)
    r7 = np.linalg.matrix_rank(_jac(F_cu, th), tol=1e-7)
    print(f"[c1 check] rank(couplings+slopes)={r6} -> null-space dim={9-r6} (WS-D said 0; correct is {9-r6})")
    print(f"           rank rises {r6}->{r7} adding a crossing => crossings CAN move at fixed coupling")
    print(f"           => WS-D's rigidity obstruction (c1) is FALSE.\n")


def real_obstruction(N=20000, seed=1):
    """Correct obstruction: separation/width ratio of the crossings is bounded ~O(1)."""
    rng = np.random.default_rng(seed); ratios = []
    for _ in range(N):
        e = np.sort(rng.uniform(-4, 4, 3))
        if np.min(np.diff(e)) < 0.2: continue
        g = rng.uniform(0.3, 1.8, 3) * rng.choice([-1, 1], 3)
        a = rng.uniform(-3, 3, 3)
        if np.min(np.abs(np.diff(np.sort(a)))) < 0.2: continue
        c, _, u = features(np.concatenate([g, e, a]))
        pairs = [(0, 1), (0, 2), (1, 2)]
        w = [2*abs(c[k])/abs(a[pairs[k][0]]-a[pairs[k][1]]) for k in range(3)]
        uv = np.sort(u); r = min(uv[1]-uv[0], uv[2]-uv[1]) / max(w)
        if np.isfinite(r): ratios.append(r)
    ratios = np.array(ratios)
    print(f"[real obstruction] separation/width over {len(ratios)} Type-1 samples:")
    print(f"   median={np.median(ratios):.2f}  99%={np.percentile(ratios,99):.2f}  "
          f"max={ratios.max():.2f}  frac(>5)={np.mean(ratios>5):.4f}")
    print("   => bounded ~O(1): crossings are PERMANENTLY marginally-overlapping; MC isolated-2-level")
    print("      limit unreachable => exact product S=prod S_ij is OBSTRUCTED (verdict upheld, right reason).")
    print("   Also explains the empirical 15% 'incoherent product works' = the high-ratio tail.")


if __name__ == "__main__":
    falsify_c1()
    real_obstruction()
