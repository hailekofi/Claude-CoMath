"""
ws_c_random_scan.py  --  WS-C: random + controlled scans correlating the factorization
defect against the trivial-middle-coupling diagnostic and the sep/width ratio.

Two experiments:
  (A) CONTROLLED weak-link sweep.  Fix a moderately-overlapping base point, then scale the
      MIDDLE level's coupling gamma_mid -> s*gamma_mid for s in [1 .. ~0] and watch
      absdef = |P_mid(exact) - P_mid(inc)| -> 0 as Gmid_min -> 0.  This directly tests
      "trivial-middle-coupling => factorizes".
  (B) RANDOM scan over Type-1 data: record (absdef, enh, Gmid_min, Gmid_prod, ratio) and
      report the correlation / the contingency table  factorizes <-> trivial-coupling.

Factorization threshold: absdef < 1e-3  (the oracle is gold to ~1e-7, so 1e-3 is a safe,
physically-meaningful "incoherent product is exact to 3 digits" cut).

Reproduce: python ws_c_random_scan.py [A|B|both].  Writes /tmp/ws_c_random.npz.  Slow (oracle).
"""
from __future__ import annotations
import sys
import numpy as np
from oracle import oracle_P, be_survivals, Gamma_ij
from coscaling import geometry as cos_geometry


def metrics(eps, gam, a, T=40.0, rtol=1e-10, atol=1e-11):
    eps = np.asarray(eps, float); gam = np.asarray(gam, float); a = np.asarray(a, float)
    bes = be_survivals(eps, gam, a)
    lo, mid, hi = bes["lo"], bes["mid"], bes["hi"]
    Gml = Gamma_ij(eps, gam, a, mid, lo); Gmh = Gamma_ij(eps, gam, a, mid, hi)
    u, w, G = cos_geometry(gam, eps, a); uv = np.sort(list(u.values()))
    ratio = min(uv[1]-uv[0], uv[2]-uv[1]) / max(w.values())
    Gall = [Gamma_ij(eps, gam, a, i, j) for (i, j) in [(0, 1), (0, 2), (1, 2)]]
    r = oracle_P(tuple(eps), tuple(gam), tuple(a), T=T, rtol=rtol, atol=atol)
    Pm = float(r["P"][mid, mid]); Pi = bes["P_mid_inc"]
    return dict(absdef=abs(Pm - Pi), enh=(Pm/Pi if Pi > 0 else np.inf),
                Gmid_lo=Gml, Gmid_hi=Gmh, Gmid_min=min(Gml, Gmh), Gmid_prod=Gml*Gmh,
                Gmin=min(Gall), Gmax=max(Gall), ratio=ratio, Pm=Pm, Pi=Pi,
                err=float(np.max(r["err"])))


def expt_A():
    """Controlled middle-coupling weak-link sweep."""
    eps = (-2.0, 0.0, 3.0); a = (-1.0, 0.5, 2.0)   # canonical; mid=index1
    base_g = np.array([1.0, 0.8, 1.2])
    print("=== (A) Weak-link sweep: scale gamma_mid (index 1) by s ===")
    print(f"{'s':>6s} {'Gmid_min':>10s} {'Gmid_prod':>11s} {'absdef':>11s} {'enh':>9s} {'ratio':>7s}")
    rows = []
    for s in [1.0, 0.7, 0.5, 0.3, 0.2, 0.1, 0.05, 0.02]:
        g = base_g.copy(); g[1] *= s
        m = metrics(eps, g, a)   # default fast settings (T=40, rtol 1e-10)
        rows.append((s, m["Gmid_min"], m["Gmid_prod"], m["absdef"], m["enh"], m["ratio"]))
        print(f"{s:6.2f} {m['Gmid_min']:10.4e} {m['Gmid_prod']:11.4e} {m['absdef']:11.4e} "
              f"{m['enh']:9.3f} {m['ratio']:7.3f}", flush=True)
    return rows


def expt_B(N=14, seed=11, T=35.0):
    """Random Type-1 scan: factorization-defect vs trivial-coupling vs ratio."""
    rng = np.random.default_rng(seed)
    recs = []
    print(f"\n=== (B) Random scan (N target {N}) ===")
    print(f"{'#':>3s} {'absdef':>10s} {'enh':>9s} {'Gmid_min':>10s} {'Gmid_prod':>11s} {'ratio':>7s} {'Gmin':>7s} {'Gmax':>7s}")
    n = 0
    while n < N:
        e = np.sort(rng.uniform(-4, 4, 3))
        if np.min(np.diff(e)) < 0.4:
            continue
        g = rng.uniform(0.2, 1.6, 3) * rng.choice([-1, 1], 3)
        a = rng.uniform(-3, 3, 3)
        if np.min(np.abs(np.diff(np.sort(a)))) < 0.2:
            continue
        try:
            m = metrics(e, g, a, T=T, rtol=1e-10, atol=1e-11)
        except Exception as ex:
            continue
        if not np.isfinite(m["absdef"]):
            continue
        recs.append(m); n += 1
        print(f"{n:3d} {m['absdef']:10.3e} {m['enh']:9.3f} {m['Gmid_min']:10.4e} "
              f"{m['Gmid_prod']:11.4e} {m['ratio']:7.3f} {m['Gmin']:7.3f} {m['Gmax']:7.3f}", flush=True)
    # analysis
    absdef = np.array([r["absdef"] for r in recs])
    Gmin_mid = np.array([r["Gmid_min"] for r in recs])
    ratio = np.array([r["ratio"] for r in recs])
    fac = absdef < 1e-3
    print(f"\n  factorizing (absdef<1e-3): {fac.sum()}/{len(recs)}")
    print(f"  among factorizing: Gmid_min  range [{Gmin_mid[fac].min():.2e}, {Gmin_mid[fac].max():.2e}]"
          if fac.any() else "  (none factorize)")
    print(f"  among NON-factorizing: Gmid_min range [{Gmin_mid[~fac].min():.2e}, {Gmin_mid[~fac].max():.2e}]"
          if (~fac).any() else "")
    # correlation of log(absdef) with log(Gmid_min)
    good = (absdef > 0) & (Gmin_mid > 0)
    if good.sum() > 3:
        c = np.corrcoef(np.log10(absdef[good]), np.log10(Gmin_mid[good]))[0, 1]
        cr = np.corrcoef(np.log10(absdef[good]), np.log10(ratio[good]))[0, 1]
        print(f"  Pearson log10(absdef) vs log10(Gmid_min) = {c:+.3f}")
        print(f"  Pearson log10(absdef) vs log10(ratio)    = {cr:+.3f}")
    np.savez("/tmp/ws_c_random.npz",
             absdef=absdef, Gmid_min=Gmin_mid, ratio=ratio,
             enh=np.array([r["enh"] for r in recs]),
             Gmid_prod=np.array([r["Gmid_prod"] for r in recs]),
             Gmin=np.array([r["Gmin"] for r in recs]),
             Gmax=np.array([r["Gmax"] for r in recs]))
    return recs


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "both"
    if which in ("A", "both"):
        expt_A()
    if which in ("B", "both"):
        expt_B()
