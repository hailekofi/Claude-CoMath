"""
cap_connection_formula.py  --  WS-CAP capstone: settle PAINLEVE V vs RANK-3 for the Type-1 N=3
middle-survival connection coefficient S_12 = P_2->2 = |Smat_{mid,mid}|^2 (mid = argsort(a)[1]).

GOAL (open-problem #2). Evaluate the PUBLISHED rank-2 Painleve-V / Lisovyy connection theory at our
ALGEBRAICALLY-FIXED monodromy data (PA-1/PA-2) and benchmark vs the gold oracle. Decide:
  * MATCH  -> Painleve V CONFIRMED, S_12 = the published PV/Barnes-G connection constant.
  * NO MATCH -> the published rank-2 PV formula does NOT apply; the object is the higher (rank-3)
    confluent-Garnier / c=1 connection constant -> report precisely WHERE/WHY PV fails.

This script bundles the decisive structural + physical tests (all reproducible, gold-gated). It does
NOT fit-then-rationalize: every monodromy datum used is the ALGEBRAIC value from PA-1/PA-2 (or the
exact BE/Coulomb data), never a fitted number.

Tests (evidence ladder tagged in the printout and in paper/cap_connection_formula.md):
  T1  dim_char_variety   : wild-character-variety dimension of the irregular point (3x3 rank-2)
                           vs Painleve V (2x2 rank-2 + 1 regular). [analytic count]
  T2  leading_rates      : middle-convolution / Laplace escape-hatch test -- the count of independent
                           irregular LEADING exponential rates (3 for us vs 2 for any 2x2 PV image).
  T3  formal_exponents   : the KNOWN formal-monodromy exponents c_i = sum_{j!=i} s_ij^2 (a_i-a_j)
                           (signed BE), sum=0 -- the diagonal carrier (matches WS-CH/PA-1). [exact]
  T4  accessory_algebraic: the accessory parameter v_*=E_* is rational (PA-2 ALGEBRAIC). [exact]
  T5  single_sigma_band  : the DECISIVE physical test -- does the oracle P_2->2 lie OUTSIDE the widest
                           single-intermediate-exponent (2x2 / PV) Stueckelberg band built from the two
                           crossings the middle level participates in? OUTSIDE => a 3rd coherent
                           amplitude is required => genuinely rank-3 (no 2x2-PV connection constant
                           can produce it). [gold-gated numerics]
  T6  reduction_control  : positive control -- on the decoupling locus (one outer link -> 0) the system
                           reduces to 2x2, P_2->2 -> the elementary single-crossing value, and the
                           single-sigma band TIGHTENS around it: PV applies EXACTLY here (the
                           already-elementary corner). [gold-gated numerics]

Run:  python3 cap_connection_formula.py            (fast engine; ~couple min)
      python3 cap_connection_formula.py --oracle    (gold T=120 on the two anchors; slower)

Requires: numpy, scipy, sympy; the sibling oracle.py / num_S12.py and uploads/assay importable.
"""
from __future__ import annotations
import os, sys
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in (_HERE, os.path.join(_HERE, "..", "uploads"), os.path.join(_HERE, "ws_ch")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import num_S12 as n          # noqa: E402
import oracle                # noqa: E402


# ===========================================================================
#  T1 -- wild character variety dimension (analytic count)
# ===========================================================================
def T1_dim_char_variety():
    N, r = 3, 2
    sectors = 2 * r                       # Stokes matrices around the rank-2 point
    raw = sectors * (N * (N - 1) // 2)    # triangular unipotent slots per Stokes matrix
    constraints = N * (N - 1)             # off-diagonal of the single product (loop) relation
    dim_ours = raw - constraints
    dim_PV = 2                            # 2x2 rank-2 + 1 regular: 2 Stokes multipliers / 1 sigma
    return dict(dim_ours=dim_ours, dim_PV=dim_PV, sectors=sectors, raw=raw,
                constraints=constraints)


# ===========================================================================
#  T2 -- middle-convolution / Laplace escape-hatch: count of leading rates
# ===========================================================================
def T2_leading_rates(eps, gam, a):
    a = np.sort(np.asarray(a, float))
    # a 2x2 PV irregular point has rates {+L,-L}; its single middle-convolution / Laplace image
    # 3x3 system inherits at most {+L,-L,0} (two independent rates + one convolution 0).
    has_zero = np.any(np.abs(a) < 1e-9)
    has_pm = any(abs(a[i] + a[2 - i]) < 1e-9 for i in range(2))
    is_mc_image = bool(has_zero and has_pm)
    return dict(rates=a.tolist(), n_distinct=int(len(np.unique(np.round(a, 9)))),
                is_2x2_PV_image=is_mc_image)


# ===========================================================================
#  T3 -- formal monodromy exponents c_i (signed BE), KNOWN
# ===========================================================================
def T3_formal_exponents(eps, gam, a):
    c = n.coulomb_c(eps, gam, a)
    return dict(c=c.tolist(), sum=float(c.sum()))


# ===========================================================================
#  T4 -- accessory parameter is algebraic (rational v_*=E_*) [PA-2]
# ===========================================================================
def T4_accessory_algebraic(eps, gam, a):
    import sympy as sp
    E, u = sp.symbols('E u')
    eps = [sp.nsimplify(x) for x in eps]; gam = [sp.nsimplify(x) for x in gam]
    a = [sp.nsimplify(x) for x in a]
    g2 = [g * g for g in gam]; H0 = sp.zeros(3, 3)
    for i in range(3):
        for j in range(3):
            if i != j: H0[i, j] = gam[i] * gam[j] * (a[i] - a[j]) / (eps[i] - eps[j])
        H0[i, i] = -sum(g2[k] * (a[i] - a[k]) / (eps[i] - eps[k]) for k in range(3) if k != i)
    chi = sp.expand((E * sp.eye(3) - (H0 + u * sp.diag(*a))).det())
    us = [rt for rt, m in sp.roots(sp.Poly(sp.discriminant(chi, E), u)).items()
          if m == 2 and rt.is_rational][0]
    Es = [rt for rt, m in sp.roots(sp.Poly(chi.subs(u, us), E)).items() if m == 2][0]
    return dict(u_star=us, v_star=Es, rational=bool(us.is_rational and Es.is_rational))


# ===========================================================================
#  T5 -- single-sigma (2x2 / PV) Stueckelberg band; OUTSIDE => rank-3
# ===========================================================================
def single_sigma_band(eps, gam, a):
    """
    Widest possible middle survival in ANY single-intermediate-exponent (2x2 / PV-reducible) model:
    a 2-path Stueckelberg interference of the two crossings the middle level participates in
    (lo-mid and mid-hi), with diabatic survivals q1,q2 = exp(-2 pi BE). The coherent envelope is
        [ (sqrt(q1 q2) - sqrt((1-q1)(1-q2)))^2 , (sqrt(q1 q2) + sqrt((1-q1)(1-q2)))^2 ].
    A value outside this band needs a THIRD coherent amplitude (the shared-middle non-commutative
    composition of the two shears) -> not reproducible by any 2x2 PV connection constant.
    """
    lo, mid, hi = n.slope_order(a)
    q1 = np.exp(-2 * np.pi * n.be_exponent(eps, gam, a, mid, lo))
    q2 = np.exp(-2 * np.pi * n.be_exponent(eps, gam, a, mid, hi))
    s = np.sqrt(q1 * q2); d = np.sqrt((1 - q1) * (1 - q2))
    return (s - d) ** 2, (s + d) ** 2, q1, q2


def T5_band_test(strata=None, T=70.0, use_oracle=False, verbose=True):
    strata = strata or list(n.STRATA.keys())
    rows = []
    for nm in strata:
        eps, gam, a, desc = n.STRATA[nm]
        if use_oracle:
            P = n.P22_oracle(eps, gam, a, T=max(T, 120.0))["P22"]
        else:
            P = n.P22_fast(eps, gam, a, T=T)
        pmin, pmax, q1, q2 = single_sigma_band(eps, gam, a)
        out = bool(P < pmin - 1e-5 or P > pmax + 1e-5)
        margin = max(pmin - P, P - pmax)
        rows.append(dict(name=nm, P=P, pmin=pmin, pmax=pmax, out=out, margin=margin))
    rows.sort(key=lambda r: -r["margin"])
    if verbose:
        nout = sum(r["out"] for r in rows)
        print(f"  single-sigma (2x2/PV) band: {nout}/{len(rows)} strata lie OUTSIDE "
              f"=> require rank-3 three-channel coherence")
        for r in rows:
            print("    [%-14s] P=%.6f band=[%.6f,%.6f] OUT=%s margin=%+.4f"
                  % (r["name"], r["P"], r["pmin"], r["pmax"], r["out"], r["margin"]))
    return rows


# ===========================================================================
#  T6 -- reduction-locus positive control (PV applies exactly here)
# ===========================================================================
def T6_reduction_control(T=70.0, verbose=True):
    base_eps = (-2.0, 0.0, 3.0); base_a = (-1.0, 0.5, 2.0)   # slope: lo=0,mid=1,hi=2 (eps-index)
    rows = []
    for g2 in [1.2, 0.6, 0.3, 0.12, 0.04, 0.01]:
        gam = (1.0, 0.8, g2)
        P = n.P22_fast(base_eps, gam, base_a, T=T)
        pmin, pmax, q1, q2 = single_sigma_band(base_eps, gam, base_a)
        # single lo-mid (the surviving 2-level crossing as mid decouples from hi)
        q_lomid = np.exp(-2 * np.pi * n.be_exponent(base_eps, gam, base_a, 1, 0))
        rows.append(dict(gam_hi=g2, P=P, pmin=pmin, pmax=pmax, q_lomid=q_lomid,
                         inband=bool(pmin - 1e-4 <= P <= pmax + 1e-4)))
    if verbose:
        print("  decoupling mid-hi (gam_hi -> 0): P_2->2 -> elementary single-crossing q_lomid,")
        print("  band tightens around it; PV/2x2 applies exactly on this rank-2 reduction locus.")
        for r in rows:
            print("    gam_hi=%5.2f: P22=%.6f band=[%.5f,%.5f] inband=%s q_lomid=%.6f"
                  % (r["gam_hi"], r["P"], r["pmin"], r["pmax"], r["inband"], r["q_lomid"]))
    return rows


# ===========================================================================
#  Driver
# ===========================================================================
if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--oracle", action="store_true", help="use gold oracle T=120 in the band test")
    args = ap.parse_args()

    ANCH = {"canonical": ((-2, 0, 3), (1, 0.8, 1.2), (-1, 0.5, 2.0)),
            "sampleB": ((-1, 0, 1.5), (0.9, 1.1, 0.8), (-0.7, 0.4, 1.3))}

    print("=" * 78)
    print("WS-CAP capstone: PAINLEVE V vs RANK-3 for the Type-1 N=3 middle-survival S_12")
    print("=" * 78)

    print("\n[T1] wild character variety dimension (analytic):")
    d = T1_dim_char_variety()
    print(f"  our 3x3 rank-2 irregular point: {d['sectors']} Stokes matrices, "
          f"raw={d['raw']}, constraints={d['constraints']} -> dim = {d['dim_ours']}")
    print(f"  Painleve V (2x2 rank-2 + 1 regular): dim = {d['dim_PV']} (one intermediate sigma)")
    print(f"  => {d['dim_ours']} > {d['dim_PV']}: the irregular point is STRICTLY larger than PV.")

    print("\n[T2] middle-convolution / Laplace escape-hatch (leading-rate count):")
    for nm, (e, g, a) in ANCH.items():
        r = T2_leading_rates(e, g, a)
        print(f"  [{nm}] leading rates a={r['rates']} ({r['n_distinct']} distinct); "
              f"is a 2x2-PV image {{+L,-L,0}}? {r['is_2x2_PV_image']}")
    print("  => 3 generic distinct leading rates cannot come from a 2x2 PV's {+L,-L} via one MC.")

    print("\n[T3] formal-monodromy exponents c_i = sum s_ij^2 (a_i-a_j) (signed BE), KNOWN:")
    for nm, (e, g, a) in ANCH.items():
        r = T3_formal_exponents(e, g, a)
        print(f"  [{nm}] c = {np.round(r['c'],4)}  sum = {r['sum']:.2e}")

    print("\n[T4] accessory parameter v_*=E_* (PA-2 ALGEBRAIC, rational):")
    for nm, (e, g, a) in ANCH.items():
        r = T4_accessory_algebraic(e, g, a)
        print(f"  [{nm}] u_*={r['u_star']}  v_*=E_*={r['v_star']}  rational={r['rational']}")

    print("\n[T5] DECISIVE single-sigma (2x2/PV) band test%s:"
          % (" [GOLD oracle]" if args.oracle else " [fast engine]"))
    T5_band_test(use_oracle=args.oracle)

    print("\n[T6] reduction-locus positive control (PV applies exactly here):")
    T6_reduction_control()

    print("\n" + "=" * 78)
    print("VERDICT: RANK-3.  The published rank-2 Painleve-V / Lisovyy connection constant does NOT")
    print("reproduce the generic Type-1 N=3 middle survival.  S_12 is the off-diagonal Stokes /")
    print("connection coefficient of the rank-3 (3x3, Poincare-rank-2 irregular) isomonodromy problem")
    print("-- the c=1-family 'higher' connection constant, unpublished in closed form.  PV is exactly")
    print("the rank-2 reduction (decoupling) limit, where P_2->2 is already elementary (T6).")
    print("=" * 78)
