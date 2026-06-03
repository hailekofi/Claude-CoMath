"""
ws_geom_m3_t4_selector.py -- BIG SWING at the Z2 orientation selector.

CORE IDEA (why the local predictor screen failed): the orientation is NOT a local
parameter sign. It is the GLOBAL spectral-flow datum "which adjacent energy-rank pair
does the real node / dominant avoided crossing connect" -- the ENERGY-POSITION of the
node in the spectrum. We compute this from the EXPLICIT spectral curve:

  char poly of H(u)=H0+uA:   p(E,u) = E^3 - c2(u) E^2 + c1(u) E - c0(u),
    c2(u) = tr H(u)         (deg 1 in u),
    c1(u) = sum 2x2 minors  (deg 2 in u),
    c0(u) = det H(u)        (deg 3 in u).

  Branch points = roots of the E-discriminant Disc_E(u) (a polynomial in u, deg <= 6).
  The branch point nearest the real axis is the dominant avoided crossing (the node).
  At u_b = Re(root), the two colliding eigenvalues form the "close pair" with value E_*;
  the spectator level is E_3 = c2(u_b) - 2 E_*.

  PREDICTOR:  orientation = sign( E_3 - E_* )   [ lower-pair node vs upper-pair node ].

We test this DETERMINISTIC, CLOSED-FORM sign against the cont_orient ground-truth Z2.

Reproduce: python3 ws_geom_m3_t4_selector.py
"""
from __future__ import annotations
import numpy as np
from collections import Counter

np.seterr(all="ignore")


def type1(eps, gam, a):
    eps = np.array(eps, float); gam = np.array(gam, float); a = np.array(a, float); g2 = gam ** 2
    H0 = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            if i != j:
                H0[i, j] = gam[i] * gam[j] * (a[i] - a[j]) / (eps[i] - eps[j])
        H0[i, i] = -sum(g2[k] * (a[i] - a[k]) / (eps[i] - eps[k]) for k in range(3) if k != i)
    return H0, np.diag(a)


def cont_orient(eps, gam, a, T=80.0, n=120001):
    """Ground-truth Z2 (T4): deterministic overlap-continuation of the order-0 eigenframe."""
    H0, A = type1(eps, gam, a)
    us = np.linspace(-T, T, n)
    _, V = np.linalg.eigh(H0 + us[0] * A); Vs = V.copy(); Vp = V.copy()
    for u in us[1:]:
        _, V = np.linalg.eigh(H0 + u * A)
        M = np.abs(Vp.conj().T @ V); V = V[:, np.argsort(np.argmax(M, axis=0))]
        for k in range(3):
            if np.vdot(Vp[:, k], V[:, k]).real < 0:
                V[:, k] *= -1
        Vp = V
    sig = [int(np.argmax(np.abs(Vs[:, k]))) for k in range(3)]
    pi = [int(np.argmax(np.abs(Vp[:, k]))) for k in range(3)]
    P = np.zeros((3, 3))
    for k in range(3):
        P[pi[k], sig[k]] = 1
    order = np.argsort(a); Ps = P[np.ix_(order, order)]
    return tuple(int(np.argmax(Ps[r])) for r in range(3))


def char_polys(H0, A):
    """Return c2(u), c1(u), c0(u) as numpy poly1d (exact: degrees 1,2,3)."""
    a = np.diag(A)
    us = np.linspace(-4, 4, 9)
    C2 = []; C1 = []; C0 = []
    for u in us:
        H = H0 + u * np.diag(a)
        c2 = np.trace(H)
        # sum of principal 2x2 minors
        c1 = (H[0,0]*H[1,1]-H[0,1]*H[1,0] + H[0,0]*H[2,2]-H[0,2]*H[2,0] + H[1,1]*H[2,2]-H[1,2]*H[2,1])
        c0 = np.linalg.det(H)
        C2.append(c2); C1.append(c1); C0.append(c0)
    c2 = np.poly1d(np.polyfit(us, C2, 1))
    c1 = np.poly1d(np.polyfit(us, C1, 2))
    c0 = np.poly1d(np.polyfit(us, C0, 3))
    return c2, c1, c0


def disc_poly(c2, c1, c0):
    """E-discriminant of x^3 + B x^2 + C x + D with B=-c2, C=c1, D=-c0, as poly1d in u."""
    B = -c2; C = c1; D = -c0
    return 18*B*C*D - 4*B*B*B*D + B*B*C*C - 4*C*C*C - 27*D*D


def real_node(H0, A):
    """The R3 universal real node: the REAL double root of the E-discriminant (a TRUE
    degeneracy, not an avoided crossing). Returns u_* (and the char polys)."""
    c2, c1, c0 = char_polys(H0, A)
    roots = disc_poly(c2, c1, c0).roots
    real = roots[np.abs(roots.imag) < 1e-6].real            # the genuine real node(s)
    if len(real) == 0:                                       # fall back to least-imag
        z = roots[np.argmin(np.abs(roots.imag))]; real = np.array([z.real])
    # pick the root with the truly vanishing eigenvalue gap (the protected node)
    a = np.diag(A)
    best = None; bestgap = np.inf
    for ur in np.unique(np.round(real, 9)):
        w = np.linalg.eigvalsh(H0 + ur * np.diag(a)); g = np.min(np.diff(w))
        if g < bestgap:
            bestgap = g; best = ur
    return float(best), c2, c1, c0


def predict_orient(eps, gam, a):
    """CLOSED-FORM predictor: at the R3 real node, sign(E_3 - E_*) -- which adjacent
    energy-rank pair is the TRUE degeneracy (lower pair k=0 vs upper pair k=1)."""
    H0, A = type1(eps, gam, a)
    ustar, c2, c1, c0 = real_node(H0, A)
    w = np.linalg.eigvalsh(H0 + ustar * np.diag(np.diag(A)))  # 3 real eigenvalues, sorted
    k = int(np.argmin(np.diff(w)))                            # degenerate pair (k,k+1)
    Estar = 0.5 * (w[k] + w[k+1])
    E3 = w[2] if k == 0 else w[0]                             # spectator level
    return int(np.sign(E3 - Estar)), k, ustar


def predict_orient_companion(eps, gam, a):
    """INDEPENDENT method: same selector from the characteristic-polynomial coefficients
    ONLY (companion-matrix roots of p(E,u)); never calls a Hermitian eigensolver on H.
    Cross-checks predict_orient via a wholly different algorithm."""
    H0, A = type1(eps, gam, a)
    c2, c1, c0 = char_polys(H0, A)
    r = disc_poly(c2, c1, c0).roots
    rr = r[np.abs(r.imag) < 1e-4].real
    if len(rr) == 0:
        rr = np.array([r[np.argmin(np.abs(r.imag))].real])
    best = None; bestsplit = np.inf
    for ur in np.unique(np.round(rr, 8)):
        e = np.sort(np.roots([1, -c2(ur), c1(ur), -c0(ur)]).real)
        split = np.min(np.diff(e))
        if split < bestsplit:
            bestsplit = split; best = (ur, e)
    ur, e = best
    k = int(np.argmin(np.diff(e)))
    Estar = 0.5 * (e[k] + e[k+1]); E3 = e[2] if k == 0 else e[0]
    return int(np.sign(E3 - Estar))


def main():
    print("=" * 92)
    print("WS-GEOM M3 T4 selector -- closed-form energy-position predictor vs cont_orient Z2")
    print("=" * 92)
    # canonical anchor
    cs, ck, cu = predict_orient([-2, 0, 3], [1, .8, 1.2], [-1, .5, 2.])
    print(f"canonical: real node u_*={cu:+.4f}, degenerate pair k={ck} "
          f"({'lower' if ck==0 else 'upper'}), sign(E3-E*)={cs:+d}  "
          f"=> {'FWD' if cs>0 else 'REV'};  cont_orient={cont_orient([-2,0,3],[1,.8,1.2],[-1,.5,2.])}")
    rows = []
    confusion = Counter()
    for seed in range(40):
        rng = np.random.default_rng(100 + seed)
        eps = sorted(rng.uniform(-2, 2, 3))
        if min(np.diff(eps)) < 0.3:
            continue
        gam = list(rng.uniform(-1.5, 1.5, 3)); a = list(rng.uniform(-2, 2, 3))
        if min(np.abs(np.diff(sorted(a)))) < 0.3 or min(abs(x) for x in a) < 0.2:
            continue
        truth = cont_orient(eps, gam, a)
        if truth not in [(1, 2, 0), (2, 0, 1)]:
            continue
        s, k, ubr = predict_orient(eps, gam, a)
        s2 = predict_orient_companion(eps, gam, a)
        rows.append((truth, s, k, ubr, s2))
        confusion[(truth, s)] += 1

    print(f"\n{'truth':>10} {'sign(E3-E*)':>12} {'pair':>9} {'u_*':>9} {'companion':>10}")
    for truth, s, k, ubr, s2 in rows:
        lab = "FWD" if truth == (1, 2, 0) else "REV"
        print(f"  {lab:>8} {s:>12d} {('lower' if k==0 else 'upper'):>9} {ubr:>9.3f} {s2:>10d}")

    print("\nConfusion (truth, predictor-sign):")
    for key, c in sorted(confusion.items()):
        lab = "FWD" if key[0] == (1, 2, 0) else "REV"
        print(f"  ({lab}, sign={key[1]:+d}): {c}")

    # is it a deterministic 2-to-2 map?
    fwd_signs = set(s for (t, s, _, _, _) in rows if t == (1, 2, 0))
    rev_signs = set(s for (t, s, _, _, _) in rows if t == (2, 0, 1))
    clean = fwd_signs.isdisjoint(rev_signs) and len(fwd_signs) <= 1 and len(rev_signs) <= 1
    agree = all(s == s2 for (_, s, _, _, s2) in rows)
    print(f"\nFWD predictor-signs: {fwd_signs}   REV predictor-signs: {rev_signs}")
    print(f"DETERMINISTIC closed-form selector (clean 2-to-2 map): {clean}")
    print(f"Hermitian and companion-root methods agree on every sample: {agree}")
    if not clean:
        print("=> energy-position sign is NOT yet the clean selector; refine (see log).")


if __name__ == "__main__":
    main()
