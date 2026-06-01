"""
restart_probe.py  --  WS-R: tests of the RESTART-OPERATOR structure of the Type-1 N=3
Landau-Zener scattering matrix at the single rank-2 irregular point u=inf.

Hypotheses under test (see paper/restart_structure.md for verdicts):

  H-R1 [lead, subtle]
      In the CANONICAL frame (symmetric diabatic IP, Coulomb-subtracted) with channels
      ordered by slope (lo<mid<hi = the sector ordering), does 𝒮_canon decompose as
          (diagonal)  ⋉  (triangular / unipotent OFF-diagonal Stokes shear)
      rather than a dense block?  We measure a WELL-DEFINED, frame-pinned quantity:
      the ratio of upper- vs lower-triangular off-diagonal weight in sector order, and
      compare to a Haar-random U(3) baseline.  We ALSO test whether any natural diagonal
      factor of 𝒮 equals the formal monodromy e^{2πi c_i}.

  H-R2 [oracle-testable]
      (a) Does the off-diagonal weight of 𝒮_canon concentrate in the MID row/column (the
          middle level is the shared subdominant partner of both {mid,lo} and {mid,hi}
          carrier pairs = WS-G's 12x13 joint)?
      (b) Under the WS-C trivial-coupling limit (send ONE outer-middle coupling s_{mid,outer}
          -> 0 by shrinking a gamma), does ONE shear vanish and the joint collapse to
          elementary (the BE-factorized / incoherent product)?

  H-R3 [decisive geometry, mostly symbolic]
      The form factor Γ_j=(-u'(λ_j))^{-1/2} is ANALYTIC at the poles ε_i (vanishes ∝(λ-ε_i),
      no monodromy there); the WKB square-root branch points are the TURNING POINTS where
      -u'(λ)=0 (the W4 zeros / ramification points of u(λ)), NOT the ε_i.  This says the
      off-diagonal Stokes structure is an irregular-point effect, distinct from the
      turning-point (window/BE) branching.

Conventions strict to NOMENCLATURE.md:
    s_ij = γ_iγ_j/(ε_i-ε_j);  w_ij=|2 s_ij|;  BE exponent = s_ij^2 |a_i-a_j|;
    c_i  = Σ_{j≠i} s_ij^2 (a_i-a_j)  (signed Coulomb coeff, Σ_i c_i=0);
    Γ_j  = δ_j (-u'(λ_j))^{-1/2}  (canonical Cauchy form factor).

Engine: the project's diabatic interaction-picture propagator (uploads/assay/ip.py),
the same one the gold oracle wraps.  numpy/scipy/sympy.
"""
from __future__ import annotations
import os, sys
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_UPLOADS = os.path.join(_HERE, "..", "uploads")
if _UPLOADS not in sys.path:
    sys.path.insert(0, _UPLOADS)

from assay import Params, Geometry           # noqa: E402
from assay.ip import _benchmark_ip_g         # noqa: E402


# ===========================================================================
#  Cauchy / Coulomb data
# ===========================================================================
def s_ij(eps, gam, i, j):
    return gam[i] * gam[j] / (eps[i] - eps[j])


def coulomb_c(eps, gam, a):
    """c_i = Σ_{j≠i} s_ij^2 (a_i-a_j)  (signed-BE Coulomb coefficient; Σ_i c_i = 0)."""
    eps = np.asarray(eps, float); gam = np.asarray(gam, float); a = np.asarray(a, float)
    return np.array([sum(s_ij(eps, gam, i, j) ** 2 * (a[i] - a[j])
                         for j in range(3) if j != i) for i in range(3)])


def be_exponent(eps, gam, a, i, j):
    """Brundobler-Elser pairwise survival exponent s_ij^2 |a_i-a_j|."""
    return s_ij(eps, gam, i, j) ** 2 * abs(a[i] - a[j])


# ===========================================================================
#  The interaction-picture / canonical scattering matrix
# ===========================================================================
def S_IP(geo: Geometry, T: float, rtol=1e-12, atol=1e-13) -> np.ndarray:
    """
    The symmetric diabatic interaction-picture scattering matrix 𝒮_IP[j,i] (row=out j,
    col=in i).  By the brief's recipe S_IP[j,i]=e^{iθ_j(T)}U(T,-T)[j,i]e^{-iθ_i(-T)} with
    θ_i(t)=(H0)_ii t + a_i t^2/2 ; this is EXACTLY the IP fundamental matrix g returned by
    _benchmark_ip_g (verified: e^{iθ}Ue^{-iθ}=g to <1e-9).
    """
    g, _ = _benchmark_ip_g(geo, T, rtol, atol)
    return g


def S_canon(geo: Geometry, T: float, c: np.ndarray, rtol=1e-12, atol=1e-13) -> np.ndarray:
    """𝒮_canon[j,i] = e^{i(c_j-c_i)logT} 𝒮_IP[j,i]  (Coulomb-subtracted, log-T convergent)."""
    g = S_IP(geo, T, rtol, atol)
    return np.exp(1j * (c[:, None] - c[None, :]) * np.log(T)) * g


# ===========================================================================
#  H-R1 : triangularity of the off-diagonal Stokes part in SECTOR (slope) order
# ===========================================================================
def slope_order(a):
    """Return (lo,mid,hi) = the sector ordering by increasing slope a."""
    lo, mid, hi = (int(k) for k in np.argsort(np.asarray(a, float)))
    return lo, mid, hi


def reorder_slope(M, a):
    """Permute a 3x3 matrix into sector (slope) order [lo,mid,hi]."""
    perm = list(slope_order(a))
    return M[np.ix_(perm, perm)]


def triangularity(Mp):
    """
    Given a matrix ALREADY in sector order [lo,mid,hi], return:
      U = Σ_{i<j}|M[i,j]|^2   (upper off-diag weight),
      L = Σ_{i>j}|M[i,j]|^2   (lower off-diag weight),
      R_tri = U/L  (>1 upper-triangular / shear-up, <1 lower, ~1 dense).
    'Triangular/unipotent shear' <=> one of U,L ~ 0.
    """
    iu = np.triu_indices(3, 1)
    il = np.tril_indices(3, -1)
    U = float(np.sum(np.abs(Mp[iu]) ** 2))
    L = float(np.sum(np.abs(Mp[il]) ** 2))
    return U, L, (U / L if L > 0 else np.inf)


def haar_u3(rng):
    z = (rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3))) / np.sqrt(2)
    q, r = np.linalg.qr(z)
    return q * (np.diag(r) / np.abs(np.diag(r)))


def haar_triangularity_baseline(nsamp=20000, seed=0):
    """Distribution of R_tri (and U,L) for Haar-random U(3): the null model."""
    rng = np.random.default_rng(seed)
    Rs, Us, Ls = [], [], []
    for _ in range(nsamp):
        U_, L_, R = triangularity(haar_u3(rng))
        Rs.append(R); Us.append(U_); Ls.append(L_)
    Rs = np.array(Rs)
    return dict(median_R=float(np.median(Rs)),
                frac_R_gt_3=float(np.mean(Rs > 3)),
                frac_R_lt_third=float(np.mean(Rs < 1 / 3)),
                mean_U=float(np.mean(Us)), mean_L=float(np.mean(Ls)))


# ===========================================================================
#  H-R2 : mid-row/column concentration of the off-diagonal weight
# ===========================================================================
def mid_concentration(Mp):
    """
    Mp in sector order [lo,mid,hi].  Off-diagonal weight that TOUCHES the middle level
    (row 1 or col 1) vs the lone lo<->hi off-diagonal weight (the pair NOT sharing mid).
    Returns (W_mid, W_outer, frac_mid).
    """
    offmask = ~np.eye(3, dtype=bool)
    P = np.abs(Mp) ** 2
    mid_mask = np.zeros((3, 3), bool)
    mid_mask[1, :] = True; mid_mask[:, 1] = True
    mid_mask &= offmask
    outer_mask = offmask & ~mid_mask        # exactly {(lo,hi),(hi,lo)}
    W_mid = float(P[mid_mask].sum())
    W_outer = float(P[outer_mask].sum())
    return W_mid, W_outer, W_mid / (W_mid + W_outer)


# ===========================================================================
#  H-R3 : form-factor analyticity and turning points (geometry)
# ===========================================================================
def turning_points(geo: Geometry):
    """The 4 zeros of W4 = n p' - n' p  (the points where -u'(λ)=W4/p^2 = 0)."""
    return np.roots(geo.W4)


def hr3_symbolic():
    """
    Symbolic confirmation (sympy):
      (1) -u'(λ) = W4(λ)/p(λ)^2  exactly  (W4 = n p' - n' p, degree 4);
      (2) near a pole ε_i,  -u' ~ g_i^2/(λ-ε_i)^2,  so Γ=( -u')^{-1/2} ~ (λ-ε_i)/|g_i|
          -> analytic, vanishes linearly (no branch point / monodromy at the poles);
      (3) the sqrt branch points of Γ are the W4 zeros (turning points), distinct from ε_i.
    Returns a dict of booleans / leading coefficients.
    """
    import sympy as sp
    l = sp.symbols('lambda')
    e = sp.symbols('e0 e1 e2', real=True)
    g = sp.symbols('g0 g1 g2', real=True)
    p = (l - e[0]) * (l - e[1]) * (l - e[2])
    n = (g[0] ** 2 * (l - e[1]) * (l - e[2]) + g[1] ** 2 * (l - e[0]) * (l - e[2])
         + g[2] ** 2 * (l - e[0]) * (l - e[1]))
    minus_up = sp.simplify(-sp.diff(n / p, l))
    W4 = sp.expand(n * sp.diff(p, l) - sp.diff(n, l) * p)
    id_ok = sp.simplify(minus_up - W4 / p ** 2) == 0
    degW4 = sp.degree(sp.Poly(W4, l))
    leadcoef = {}
    for i in range(3):
        lead = sp.simplify(sp.series(minus_up * (l - e[i]) ** 2, l, e[i], 1).removeO())
        leadcoef[i] = (sp.simplify(lead - g[i] ** 2) == 0)   # == g_i^2 ?
    return dict(minus_up_eq_W4_over_p2=bool(id_ok), degW4=int(degW4),
                pole_leading_is_gi2=leadcoef)


# ===========================================================================
#  Drivers
# ===========================================================================
SAMPLES = {
    "canonical": ((-2.0, 0.0, 3.0), (1.0, 0.8, 1.2), (-1.0, 0.5, 2.0)),
    "sampleB":   ((-1.0, 0.0, 1.5), (0.9, 1.1, 0.8), (-0.7, 0.4, 1.3)),
    "well_sep":  ((-5.0, 0.0, 5.0), (1.0, 1.0, 1.0), (-1.5, 0.0, 1.5)),
}


def run_one(name, eps, gam, a, T=120.0, rtol=1e-11, atol=1e-12, verbose=True):
    geo = Geometry(Params(eps=eps, gam=gam, a=a))
    c = coulomb_c(eps, gam, a)
    Sc = S_canon(geo, T, c, rtol, atol)
    Scp = reorder_slope(Sc, a)                       # sector order [lo,mid,hi]
    U, L, Rtri = triangularity(Scp)
    Wmid, Wouter, fmid = mid_concentration(Scp)
    P = np.abs(Scp) ** 2
    res = dict(name=name, c=c, S_canon_sector=Scp, P_sector=P,
               U=U, L=L, Rtri=Rtri, Wmid=Wmid, Wouter=Wouter, frac_mid=fmid,
               turning_points=turning_points(geo), eps=eps)
    if verbose:
        lo, mid, hi = slope_order(a)
        print(f"\n=== {name}  (slope order lo,mid,hi = {lo},{mid},{hi}) ===")
        print(f"  c_i = {np.round(c,4)}  (Σ={c.sum():.1e})   formal monodromy e^(2πi c_i):")
        print(f"        {np.round(np.exp(2j*np.pi*c),4)}")
        with np.printoptions(precision=4, suppress=True, linewidth=120):
            print("  |𝒮_canon|^2  in sector order [lo,mid,hi]:\n   ",
                  str(np.round(P, 4)).replace("\n", "\n    "))
        print(f"  H-R1  upper off-diag wt U={U:.4f}  lower L={L:.4f}  R_tri=U/L={Rtri:.3f}")
        print(f"  H-R2  mid-touching off-diag wt={Wmid:.4f}  outer(lo-hi)={Wouter:.4f}  "
              f"frac_mid={fmid:.3f}")
    return res


def hr2_decoupling(name, eps, gam, a, link_scale=1e-3, T=120.0,
                   rtol=1e-11, atol=1e-12, verbose=True):
    """
    H-R2(b): WS-C trivial-coupling limit.  Decouple ONE outer level (lo or hi) from the
    middle by shrinking its gamma, which sends the corresponding outer-middle coupling
    s_{mid,outer} -> 0 (and the lo<->hi link too).  The remaining mid<->(other outer) link
    is a single 2-level (Weber) crossing -> the joint should collapse to ELEMENTARY, i.e.
    the middle survival P_mid -> the incoherent BE product (here a single surviving factor,
    since the decoupled pair's BE exponent -> 0 so its factor -> 1).  We report the full
    case and both decoupled variants; one shear is expected to vanish in each.
    """
    lo, mid, hi = slope_order(a)
    out = {}
    geo0 = Geometry(Params(eps=eps, gam=gam, a=a))
    c0 = coulomb_c(eps, gam, a)
    Scp0 = reorder_slope(S_canon(geo0, T, c0, rtol, atol), a)
    P0 = np.abs(Scp0) ** 2
    inc0 = np.exp(-2 * np.pi * (be_exponent(eps, gam, a, mid, lo)
                               + be_exponent(eps, gam, a, mid, hi)))
    Wmid0, Wouter0, fmid0 = mid_concentration(Scp0)
    out["full"] = dict(P_mid=float(P0[1, 1]), inc=float(inc0),
                       Wmid=Wmid0, Wouter=Wouter0, frac_mid=fmid0)
    if verbose:
        print(f"  [{name}/full]        P_mid={P0[1,1]:.6f}  incoherent BE={inc0:.6f}  "
              f"|diff|={abs(P0[1,1]-inc0):.2e}  Wmid={Wmid0:.4f}")
    for tag, idx in [("decouple_hi", hi), ("decouple_lo", lo)]:
        g2 = list(gam); g2[idx] = gam[idx] * link_scale
        geo = Geometry(Params(eps=eps, gam=tuple(g2), a=a))
        c = coulomb_c(eps, tuple(g2), a)
        Scp = reorder_slope(S_canon(geo, T, c, rtol, atol), a)
        P = np.abs(Scp) ** 2
        inc = np.exp(-2 * np.pi * (be_exponent(eps, tuple(g2), a, mid, lo)
                                  + be_exponent(eps, tuple(g2), a, mid, hi)))
        Wmid, Wouter, fmid = mid_concentration(Scp)
        out[tag] = dict(P_mid=float(P[1, 1]), inc=float(inc),
                        Wmid=Wmid, Wouter=Wouter, frac_mid=fmid,
                        P_mid_minus_inc=float(abs(P[1, 1] - inc)))
        if verbose:
            print(f"  [{name}/{tag}] P_mid={P[1,1]:.6f}  incoherent BE={inc:.6f}  "
                  f"|diff|={abs(P[1,1]-inc):.2e}  Wmid={Wmid:.4f}  frac_mid={fmid:.4f}")
    return out


def main():
    print("#" * 78)
    print("# H-R3 (symbolic geometry): form factor analyticity & turning points")
    print("#" * 78)
    sym = hr3_symbolic()
    print(f"  -u'(λ) = W4/p^2 (exact, sympy): {sym['minus_up_eq_W4_over_p2']}")
    print(f"  deg W4 = {sym['degW4']} (=> 4 turning points)")
    print(f"  near pole ε_i: (λ-ε_i)^2·(-u') -> g_i^2 for i=0,1,2: "
          f"{[sym['pole_leading_is_gi2'][i] for i in range(3)]}")
    print("  => Γ_j=(-u')^{-1/2} ~ (λ-ε_i)/|g_i| : ANALYTIC at poles (no monodromy).")
    print("  => branch points of Γ = zeros of W4 (turning points, -u'=0) != ε_i.")

    print("\n" + "#" * 78)
    print("# H-R1 baseline: Haar-random U(3) triangularity null model")
    print("#" * 78)
    base = haar_triangularity_baseline()
    print(f"  Haar U(3): median R_tri={base['median_R']:.3f}, "
          f"P(R>3)={base['frac_R_gt_3']:.3f}, P(R<1/3)={base['frac_R_lt_third']:.3f}, "
          f"mean U={base['mean_U']:.3f}, mean L={base['mean_L']:.3f}")

    print("\n" + "#" * 78)
    print("# H-R1 & H-R2(a): 𝒮_canon triangularity + mid-concentration on samples")
    print("#" * 78)
    results = {}
    for name, (eps, gam, a) in SAMPLES.items():
        results[name] = run_one(name, eps, gam, a)
        print("  turning points (W4 zeros):", np.round(results[name]["turning_points"], 4))
        print("  poles eps:", eps, " min|tp-eps|:",
              min(abs(t - e) for t in results[name]["turning_points"] for e in eps))

    print("\n" + "#" * 78)
    print("# H-R2(b): WS-C decoupling limit -> does one shear vanish (elementary)?")
    print("#" * 78)
    for name in ("canonical", "sampleB"):       # the two jointed (non-elementary) cases
        eps, gam, a = SAMPLES[name]
        print(f"\n--- {name} ---")
        hr2_decoupling(name, eps, gam, a)

    return results, base, sym


if __name__ == "__main__":
    main()
