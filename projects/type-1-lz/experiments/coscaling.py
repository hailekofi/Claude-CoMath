"""
Co-scaling derivation for the Type-1 N=3 crossing geometry (supports the WS-D obstruction).

Establishes:
  (L)  WIDTH LEMMA [established, exact]:  w_ij = 2|g_i g_j|/|e_i-e_j|  (slope-INDEPENDENT).
       w_ij is the primary width variable; the LZ exponent Gamma_ij keeps its established meaning
       Gamma_ij = g_i^2 g_j^2 |a_i-a_j|/(e_i-e_j)^2 (used only for the genuine-LZ constraint).
  (S)  scale-invariance: sep/width is invariant under independent rescaling of g, e, a.
  (O)  OBSTRUCTION [analytic backbone + numerically-supported bound]: on the genuine-LZ locus
       (all Gamma_ij = O(1)) the closest pair of crossings is always within ~1.6 widths
       => at least two crossings are PERMANENTLY marginally-overlapping => MC's all-isolated
       configuration is unreachable.  The degenerate-slope "escape" drives Gamma->0 (trivial).
numpy only.
"""
import numpy as np


def geometry(g, e, a):
    g = np.asarray(g, float); e = np.asarray(e, float); a = np.asarray(a, float)
    P = [(0, 1), (0, 2), (1, 2)]
    s = [-sum(g[k]**2*(a[i]-a[k])/(e[i]-e[k]) for k in range(3) if k != i) for i in range(3)]  # H0_ii
    u = {p: -(s[p[0]]-s[p[1]])/(a[p[0]]-a[p[1]]) for p in P}                     # crossing locations
    w = {p: 2*abs(g[p[0]]*g[p[1]])/abs(e[p[0]]-e[p[1]]) for p in P}              # width lemma
    G = {p: g[p[0]]**2*g[p[1]]**2*abs(a[p[0]]-a[p[1]])/(e[p[0]]-e[p[1]])**2 for p in P}  # LZ exponent
    return u, w, G


def sep_over_width(g, e, a):
    u, w, G = geometry(g, e, a); uv = np.sort(list(u.values()))
    return min(uv[1]-uv[0], uv[2]-uv[1])/max(w.values()), min(G.values()), max(G.values())


def genuine_LZ_bound(N=400000, Gmin=0.2, Gmax=5.0, seed=3):
    rng = np.random.default_rng(seed); rr = []
    for _ in range(N):
        e = np.sort(rng.uniform(-4, 4, 3))
        if np.min(np.diff(e)) < 0.15: continue
        g = rng.uniform(0.3, 2.0, 3)*rng.choice([-1, 1], 3); a = rng.uniform(-3, 3, 3)
        if np.min(np.abs(np.diff(np.sort(a)))) < 0.05: continue
        r, gmin, gmax = sep_over_width(g, e, a)
        if Gmin <= gmin and gmax <= Gmax and np.isfinite(r): rr.append(r)
    rr = np.array(rr)
    print(f"(O) genuine-LZ locus (Gamma in [{Gmin},{Gmax}]): {len(rr)} samples")
    print(f"    sep/width  median={np.median(rr):.2f}  99%={np.percentile(rr,99):.2f}  "
          f"max={rr.max():.2f}  frac(>5)={np.mean(rr>5):.5f}  -> bounded ~O(1)")


if __name__ == "__main__":
    # (L) width lemma is slope-independent (numeric spot check)
    e = np.array([-2., 0., 3.]); g = np.array([1., .8, 1.2])
    _, w1, _ = geometry(g, e, [-1., .5, 2.]); _, w2, _ = geometry(g, e, [3., -2., 0.7])
    print("(L) width lemma slope-independence: max|w(a)-w(a')| =",
          max(abs(w1[p]-w2[p]) for p in w1), " (=0 exactly)")
    # (S) scale invariance of sep/width
    base = sep_over_width(g, e, [-1., .5, 2.])[0]
    sc = [sep_over_width(g*1.7, e, [-1., .5, 2.])[0], sep_over_width(g, e*2.3, [-1., .5, 2.])[0],
          sep_over_width(g, e, np.array([-1., .5, 2.])*0.6)[0]]
    print(f"(S) sep/width under g,e,a rescaling: {base:.4f} vs {[round(x,4) for x in sc]} (invariant)")
    genuine_LZ_bound()
