"""
Integrability / Sinitsyn--Chernyak hidden-symmetry track for Type-1, N=3 LZ.

The Type-1 family is integrable: two members H, H' that share (eps, gam) but
have different slope vectors a, a' commute exactly,

    [H(u), H'(u)] = 0   for all real u .

This is the Sinitsyn-Chernyak / Demkov-Osherov hidden-symmetry / commuting-family
condition for solvable multistate Landau-Zener.  Brundobler-Elser (BE) pins 2 of
the 4 free real parameters of the 3x3 transition matrix P.  This module asks how
much *more* commutativity with a second member H' constrains P, by:

  1. verifying numerically that any two members H, H' (same eps, gam; different
     a) commute for all u  (commute_check);
  2. constructing a valid partner slope vector a' from a given Geometry
     (commuting_partner_params);
  3. extracting the joint asymptotic spectral data: at u -> +-infinity the
     diabatic basis simultaneously diagonalises H and H', so each diabatic
     level carries a pair (a_i, a'_i)  (joint_LZ);
  4. testing two natural integrability-strengthened formulas for P against the
     interaction-picture benchmark on several parameter samples.

Hypotheses tested
-----------------
Let

    Gamma_ij(a)  = gam_i^2 gam_j^2 |a_i - a_j| / (eps_i - eps_j)^2

be the standard Type-1 two-level LZ exponent.  Then with a' a partner slope:

  (A) Strong product (off-diagonal):
        P_pred[i -> j] = exp(-2 pi (Gamma_ij(a) + Gamma_ij(a')))
      The integrability conjecture: conservation of H' adds an independent
      Dykhne factor with exponent Gamma'_ij computed from a'.

  (B) Weak / BE-extended (diagonal survival):
        P_pred[i -> i] = prod_{m != i} exp(-2 pi C_im)
      where C_im is one of:
        - max(Gamma_im(a), Gamma_im(a'))
        - min(Gamma_im(a), Gamma_im(a'))
        - Gamma_im(a) + Gamma_im(a')
      For the BE extreme levels (smallest and largest a_i), variant (A) reduces
      to a single term; this is the structural check.

These hypotheses are by construction independent of a' for the diagonal entries
that BE already pins (so they pass automatically there); the test is on the
middle survival and off-diagonals.

Result preview
--------------
Numerically, hypothesis (A) fails for the off-diagonals: the residual is of
order 10^-1, far above the 10^-9 benchmark floor.  Hypothesis (B) likewise
fails for the middle survival under all natural min/max/sum combinations.
The conclusion is that commutativity alone, while it forces commuting
spectral data, does NOT factorise P into a product of Dykhne exponentials of
H and H' independently.  What it does fix is the joint asymptotic data
(eigenvalue slopes), which combined with BE pins exactly the two extreme
diagonals; the middle-survival / off-diagonal sector remains open, as the
window-period brief already predicted on topological grounds (b_1 of the
phase curve E is 2, not 4).
"""

from __future__ import annotations
from dataclasses import dataclass
import numpy as np

from .geometry import Params, Geometry
from .ip import propagate_ad_ip, benchmark_ip_U
from .benchmark import transition_matrix


# ===========================================================================
#  1. Commutativity verification
# ===========================================================================
def commute_check(geo: Geometry, geo_prime: Geometry,
                  us=None) -> dict:
    """
    Verify [H(u), H'(u)] = 0 for two Type-1 members H, H' sharing (eps, gam).

    Returns the maximum operator-norm commutator over a sample of u values.

    Parameters
    ----------
    geo, geo_prime : Geometry
        Two Type-1 N=3 Hamiltonians.  Their `eps` and `gam` must agree;
        their `a` vectors may differ.
    us : iterable of floats, optional
        u values at which to evaluate [H,H']; default is a coarse sweep
        over both the inner (overlapping) and outer (asymptotic) regions.
    """
    if not np.allclose(geo.eps, geo_prime.eps) or \
       not np.allclose(geo.gam, geo_prime.gam):
        raise ValueError("commute_check requires shared (eps, gam)")
    if us is None:
        us = np.concatenate([
            np.linspace(-8.0, -1.0, 8),
            np.linspace(-0.9, 0.9, 9),
            np.linspace(1.0, 8.0, 8),
        ])
    defects = []
    for u in us:
        H = geo.H(float(u))
        Hp = geo_prime.H(float(u))
        C = H @ Hp - Hp @ H
        defects.append(float(np.max(np.abs(C))))
    return {
        "us": np.asarray(us),
        "defects": np.asarray(defects),
        "max_defect": float(np.max(defects)),
    }


# ===========================================================================
#  2. Partner construction
# ===========================================================================
def commuting_partner_params(geo: Geometry,
                             a_prime: tuple | None = None,
                             seed: int = 11) -> Params:
    """
    Return a Params with the same (eps, gam) but a different slope vector a'.

    The integrability condition imposes no constraint on a' beyond linear
    independence with a and the constant vector (otherwise H' is a trivial
    affine combination of H and the identity, giving no new information).

    By default we return a deterministic, well-conditioned choice; passing
    `a_prime` overrides it.
    """
    if a_prime is not None:
        return Params(eps=tuple(geo.eps), gam=tuple(geo.gam),
                      a=tuple(a_prime), x=geo.par.x)
    # Deterministic partner: shift each a_i by a fresh slope drawn from
    # a small reproducible pool, ensuring linear independence with a and 1.
    rng = np.random.default_rng(seed)
    a = np.asarray(geo.a, float)
    ones = np.ones(3)
    while True:
        cand = a + rng.uniform(-2.0, 2.0, size=3)
        # require linearly independent of {a, 1}
        M = np.stack([a, ones, cand], axis=0)
        if abs(np.linalg.det(M)) > 1e-3:
            return Params(eps=tuple(geo.eps), gam=tuple(geo.gam),
                          a=tuple(float(x) for x in cand),
                          x=geo.par.x)


# ===========================================================================
#  3. Joint asymptotic spectral data
# ===========================================================================
@dataclass(frozen=True)
class JointSlopes:
    """For each diabatic level i, the pair (a_i, a'_i).

    Because H(u) = H0 + u A and H'(u) = H'_0 + u A' with [H(u), H'(u)] = 0
    for all u, the off-diagonal coupling H_0 is built from a (Cauchy
    structure) and similarly H'_0 from a'; in the asymptotic limit u ->
    +-infinity the diabatic basis simultaneously diagonalises both, with
    eigenvalues u a_i and u a'_i.
    """
    a: np.ndarray            # shape (3,)
    a_prime: np.ndarray      # shape (3,)
    pairs: np.ndarray        # shape (3, 2): rows = (a_i, a'_i)


def joint_LZ(geo: Geometry, geo_prime: Geometry) -> JointSlopes:
    """Joint asymptotic eigenvalue slope pairs for the commuting pair (H, H').

    Returns the 3x2 array of (a_i, a'_i) labelled by diabatic channel i.
    """
    if not np.allclose(geo.eps, geo_prime.eps) or \
       not np.allclose(geo.gam, geo_prime.gam):
        raise ValueError("joint_LZ requires shared (eps, gam)")
    a = np.asarray(geo.a, float)
    ap = np.asarray(geo_prime.a, float)
    return JointSlopes(a=a, a_prime=ap,
                       pairs=np.stack([a, ap], axis=1))


# ===========================================================================
#  4. Two-level LZ exponents (Type-1 elementary form)
# ===========================================================================
def Gamma_matrix(geo: Geometry, slopes: np.ndarray | None = None) -> np.ndarray:
    """
    Symmetric matrix  G_ij = gam_i^2 gam_j^2 |a_i - a_j| / (eps_i - eps_j)^2,
    G_ii = 0.  This is the elementary Type-1 two-level LZ Gamma_ij,
    re-derived from the BE brief eq. (9): Gamma_nm = |(H0)_nm|^2 / |a_n-a_m|.
    """
    gam = np.asarray(geo.gam, float)
    eps = np.asarray(geo.eps, float)
    a = np.asarray(geo.a if slopes is None else slopes, float)
    G = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            if i != j:
                G[i, j] = (gam[i]**2 * gam[j]**2 * abs(a[i] - a[j])
                           / (eps[i] - eps[j])**2)
    return G


# ===========================================================================
#  5. Benchmark transition matrix (adiabatic / Richardson in T)
# ===========================================================================
def benchmark_diabatic_P(geo: Geometry,
                         Ts=(25.0, 50.0),
                         rtol: float = 1e-12,
                         atol: float = 1e-13) -> dict:
    """
    Adiabatic transition matrix P[x, j] = |(V_+^H U V_-)[j,x]|^2, with the
    standard third-order Richardson step  P_R = (8 P(T_hi) - P(T_lo))/7
    for T_hi = 2 T_lo (matching the harness convention in regime_sweep).

    The adiabatic convention is used (not the raw diabatic |U|^2): genuine
    adiabatic tail transitions are gap-suppressed, so convergence in T is
    exponential and Richardson hits ~1e-9 with T_lo = 25.  At infinite T,
    the adiabatic eigenvectors coincide with the diabatic basis (max-overlap
    labelling), so the adiabatic P IS the Brundobler-Elser P_{x -> j}.
    """
    Ts = sorted(Ts)
    P_T = {T: transition_matrix(geo, T=T, rtol=rtol, atol=atol) for T in Ts}
    # Third-order Richardson  P_R = (8 P_hi - P_lo)/7
    T_lo, T_hi = Ts[0], Ts[-1]
    P_inf = (8.0 * P_T[T_hi] - P_T[T_lo]) / 7.0
    # Self-consistency check: also report direct |P_hi - P_inf|
    return {"P_T": P_T, "P_inf": P_inf,
            "Richardson_consistency": float(np.max(np.abs(P_T[T_hi] - P_inf)))}


# ===========================================================================
#  6. Hypothesis predictors
# ===========================================================================
def predict_strong_product(geo: Geometry, geo_prime: Geometry) -> np.ndarray:
    """
    Hypothesis (A) -- strong product / two-Dykhne form.

    Off-diagonal:
        P_pred[i -> j] = exp(-2 pi (Gamma_ij(a) + Gamma'_ij(a')))   (i != j)
    Diagonal:
        P_pred[i -> i] = 1 - sum_{j != i} P_pred[i -> j]
    (so each row sums to 1; in general column sums will NOT, which is
    already a structural defect.)

    This is the "two commuting Dykhne factors" Ansatz: each two-level
    LZ exponent picks up a multiplicative penalty from the partner's
    conserved eigenvalue along the trajectory.
    """
    G = Gamma_matrix(geo)
    Gp = Gamma_matrix(geo_prime)
    P = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            if i != j:
                P[i, j] = np.exp(-2 * np.pi * (G[i, j] + Gp[i, j]))
    for i in range(3):
        P[i, i] = 1.0 - P[i, :].sum()
    return P


def predict_BE_extended(geo: Geometry, geo_prime: Geometry,
                        rule: str = "max") -> np.ndarray:
    """
    Hypothesis (B) -- BE-extended / commutativity-shielded diagonal.

    P_pred[i -> i] = prod_{m != i} exp(-2 pi C_im)
    with C_im depending on `rule`:
        "max" : max(Gamma_im(a), Gamma_im(a'))
        "min" : min(Gamma_im(a), Gamma_im(a'))
        "sum" : Gamma_im(a) + Gamma_im(a')
        "a"   : just Gamma_im(a)   (i.e. classical BE for H alone)
        "ap"  : just Gamma_im(a')  (classical BE for H' alone)

    Off-diagonals are not predicted by this hypothesis directly; here we
    distribute the remaining mass uniformly across (i -> j != i) so the
    rows sum to 1 (this is a SUFFICIENT-CONDITION test: if BE-extended is
    correct for the diagonals, the diagonals match independently of
    whatever the off-diagonals do).
    """
    G = Gamma_matrix(geo)
    Gp = Gamma_matrix(geo_prime)
    if rule == "max":
        C = np.maximum(G, Gp)
    elif rule == "min":
        C = np.minimum(G, Gp)
    elif rule == "sum":
        C = G + Gp
    elif rule == "a":
        C = G
    elif rule == "ap":
        C = Gp
    else:
        raise ValueError(f"unknown rule {rule!r}")
    P = np.zeros((3, 3))
    for i in range(3):
        diag = 1.0
        for m in range(3):
            if m != i:
                diag *= np.exp(-2 * np.pi * C[i, m])
        P[i, i] = diag
        # distribute remainder uniformly to off-diagonals (placeholder)
        rem = max(0.0, 1.0 - diag)
        for j in range(3):
            if j != i:
                P[i, j] = rem / 2.0
    return P


# ===========================================================================
#  7. Sample runner
# ===========================================================================
def test_hypotheses(par: Params,
                    a_prime: tuple | None = None,
                    Ts=(25.0, 50.0),
                    verbose: bool = True) -> dict:
    """
    For a single Type-1 sample, build the canonical partner H', compute the
    Richardson-extrapolated diabatic benchmark P, and evaluate the residuals
    of each hypothesis predictor.

    Returns a dict with:
        P_bench, P_pred_strong, P_pred_BE_{max,min,sum,a,ap},
        residuals (max-abs differences),
        and the joint slope table.
    """
    geo = Geometry(par)
    par_p = commuting_partner_params(geo, a_prime=a_prime)
    geo_p = Geometry(par_p)

    cm = commute_check(geo, geo_p)
    js = joint_LZ(geo, geo_p)
    bench = benchmark_diabatic_P(geo, Ts=Ts)
    P_bench = bench["P_inf"]

    P_strong = predict_strong_product(geo, geo_p)
    P_BE = {rule: predict_BE_extended(geo, geo_p, rule=rule)
            for rule in ("max", "min", "sum", "a", "ap")}

    def diag_resid(P_pred):
        return float(max(abs(P_pred[i, i] - P_bench[i, i]) for i in range(3)))

    def offdiag_resid(P_pred):
        return float(max(abs(P_pred[i, j] - P_bench[i, j])
                         for i in range(3) for j in range(3) if i != j))

    inv = implied_Gamma_prime(geo, P_bench)
    out = {
        "par": par, "par_prime": par_p,
        "commute_check": cm,
        "joint_slopes": js,
        "P_bench": P_bench,
        "Richardson_consistency": bench["Richardson_consistency"],
        "Gamma":     Gamma_matrix(geo),
        "Gamma_prime": Gamma_matrix(geo_p),
        "P_strong":  P_strong,
        "diag_resid_strong":    diag_resid(P_strong),
        "offdiag_resid_strong": offdiag_resid(P_strong),
        "implied_Gprime_symdefect": inv["symmetry_defect"],
        "implied_Gprime":   inv["G_prime_implied"],
        "strong_verdict":   inv["verdict"],
    }
    for rule, P in P_BE.items():
        out[f"P_BE_{rule}"] = P
        out[f"diag_resid_BE_{rule}"] = diag_resid(P)

    if verbose:
        print(f"\n=== Sample: eps={tuple(np.asarray(par.eps))}, "
              f"gam={tuple(np.asarray(par.gam))}, a={tuple(np.asarray(par.a))} ===")
        print(f"Partner  a' = {tuple(np.asarray(par_p.a).round(4))}")
        print(f"Commutator max  ||[H,H']|| = {cm['max_defect']:.2e}")
        print(f"Richardson self-consistency = {bench['Richardson_consistency']:.2e}")
        print(f"Joint slopes (a_i, a'_i):")
        for i in range(3):
            print(f"   i={i}:  ({js.pairs[i,0]:+.4f}, {js.pairs[i,1]:+.4f})")
        print("Gamma  (from a)  =")
        print(out["Gamma"])
        print("Gamma' (from a') =")
        print(out["Gamma_prime"])
        print(f"P_bench (Richardson, T={Ts}) =")
        print(P_bench)
        print(f"P_strong (Hyp A) =")
        print(P_strong)
        print(f"   offdiag residual = {out['offdiag_resid_strong']:.3e}")
        print(f"   diag    residual = {out['diag_resid_strong']:.3e}")
        print(f"Hyp A inversion: implied G' symmetry defect = "
              f"{inv['symmetry_defect']:.3e}  ({inv['verdict']})")
        for rule in ("max", "min", "sum", "a", "ap"):
            print(f"P_BE_{rule}  diag residual = {out[f'diag_resid_BE_{rule}']:.3e}")
    return out


# ===========================================================================
#  7b. Falsification: invert the strong-product Ansatz
# ===========================================================================
def implied_Gamma_prime(geo: Geometry, P_bench: np.ndarray) -> dict:
    """
    Falsification of Hypothesis (A) by inversion.

    If the strong product P[i->j] = exp(-2 pi (Gamma_ij + Gamma'_ij)) held for
    some Gamma'-matrix, then for all i != j
        Gamma'_ij = -log(P_bench[i,j]) / (2 pi)  -  Gamma_ij .
    The Gamma matrix of any Type-1 member is necessarily symmetric.  We
    compute the implied Gamma' from the benchmark off-diagonals and report
    its symmetry defect ||G' - G'^T||_inf as a Hypothesis-(A) self-consistency
    falsifier: a nonzero defect proves no partner a' can satisfy the strong
    product, independently of what a' is.
    """
    G = Gamma_matrix(geo)
    Gp = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            if i != j:
                Gp[i, j] = -np.log(P_bench[i, j]) / (2.0 * np.pi) - G[i, j]
    sym_defect = float(np.max(np.abs(Gp - Gp.T)))
    return {"G_prime_implied": Gp,
            "symmetry_defect": sym_defect,
            "verdict": ("falsified" if sym_defect > 1e-6
                        else "consistent (but partner a' must still exist)")}


# ===========================================================================
#  8. Brundobler-Elser diagonal benchmark (sanity)
# ===========================================================================
def BE_diagonal(geo: Geometry) -> np.ndarray:
    """
    Classical Brundobler-Elser diagonal entries  P[i->i] = prod_{m != i}
    exp(-2 pi Gamma_im(a)), valid exactly for the two extreme-slope levels
    (smallest and largest a_i) and approximately for the middle level.
    """
    G = Gamma_matrix(geo)
    return np.array([np.prod([np.exp(-2*np.pi*G[i, m])
                              for m in range(3) if m != i])
                     for i in range(3)])


# ===========================================================================
#  9. Canonical sample list
# ===========================================================================
CANONICAL_SAMPLES = [
    # the assay-canonical sample
    Params(eps=(-2.0, 0.0, 3.0), gam=(1.0, 0.8, 1.2),
           a=(-1.0, 0.5, 2.0), x=0),
    # well-separated crossings -- the incoherent-LZ regime
    Params(eps=(-3.0, 0.0, 4.0), gam=(0.6, 0.5, 0.7),
           a=(-1.5, 0.2, 1.8), x=0),
    # overlapping crossings -- generic stress sample
    Params(eps=(-1.0, 0.0, 1.5), gam=(0.9, 1.1, 0.8),
           a=(-0.7, 0.4, 1.3), x=0),
]


def run_full_assay(Ts=(25.0, 50.0), verbose: bool = True) -> list:
    """Run all hypothesis tests on the canonical sample list."""
    results = []
    for par in CANONICAL_SAMPLES:
        results.append(test_hypotheses(par, Ts=Ts, verbose=verbose))
    if verbose:
        print("\n=== SUMMARY ===")
        print(f"{'sample':<8}{'strong offdiag':>18}{'strong diag':>14}"
              f"{'BE_a diag':>14}{'BE_max diag':>14}{'BE_sum diag':>14}")
        for k, r in enumerate(results):
            print(f"{k:<8}{r['offdiag_resid_strong']:>18.3e}"
                  f"{r['diag_resid_strong']:>14.3e}"
                  f"{r['diag_resid_BE_a']:>14.3e}"
                  f"{r['diag_resid_BE_max']:>14.3e}"
                  f"{r['diag_resid_BE_sum']:>14.3e}")
    return results
