"""
Frame-map test for the Level-3 window factorization.

Question: can the exact gauged-adiabatic window propagator
    W_real(X) = propagate_ad_ip(geo, l, r, x)
be written as  F_out . (2+1 block) . F_in^{-1}  with F_in, F_out built
ALGEBRAICALLY from the Cauchy/pole frame V[i,a] = Gam_i gam_a/(lam_i-eps_a)
(or related algebraic frames) evaluated at the window endpoints l, r?

A frame map M "block-diagonalizes" W into a (2+1) structure for crossing
pair (i,j) with spectator s if  M^{-1} W M  (or F_in^{-1} W F_out etc.)
has negligible entries in the spectator's off-diagonal row/column:
    rows/cols s coupling to {i,j} all ~ 0.

We try, for each window, every candidate (F_in, F_out) drawn from:
  - Cauchy frame V at l and r,
  - adiabatic eigenvector matrix of H(u) at l and r,
  - identity,
and report the residual block-off-diagonal magnitude.

No existing assay module is modified.
"""

from __future__ import annotations
import numpy as np

from .geometry import Geometry, Params
from .adiabatic import cauchy_frame, labelled_lambdas
from .ip import propagate_ad_ip
from .closed_form import window_pairs
from .windows import segmentation


# ---------------------------------------------------------------------------
#  algebraic frames at an endpoint
# ---------------------------------------------------------------------------
def adiabatic_eigvecs(geo: Geometry, u: float) -> np.ndarray:
    """
    Eigenvector matrix of H(u), columns = eigenvectors, ordered so that
    column k corresponds to labelled lam_k = labelled_lambdas[k].
    Returns a 3x3 matrix Vec with Vec[:,k] the k-th adiabatic eigenvector.
    """
    H = geo.H(u)
    w, vec = np.linalg.eigh(H)
    lam = labelled_lambdas(geo, u)
    # match eigenvalues w to labelled lam
    order = np.empty(3, int)
    pool = list(range(3))
    for k in range(3):
        j = int(np.argmin([abs(w[p] - lam[k].real) for p in pool]))
        order[k] = pool.pop(j)
    return vec[:, order].astype(complex)


def block_offdiag_norm(M: np.ndarray, pair: tuple) -> float:
    """
    Max |entry| of M that couples the spectator to the crossing pair.
    M is in a basis whose 3 indices are adiabatic sheets {0,1,2}; the
    spectator s should decouple => row s and col s off-diagonal ~ 0.
    """
    s = ({0, 1, 2} - set(pair)).pop()
    bad = []
    for k in range(3):
        if k != s:
            bad.append(abs(M[s, k]))
            bad.append(abs(M[k, s]))
    return float(max(bad))


def block_2x2(M: np.ndarray, pair: tuple) -> np.ndarray:
    i, j = sorted(pair)
    return M[np.ix_([i, j], [i, j])]


# ---------------------------------------------------------------------------
#  main test
# ---------------------------------------------------------------------------
def run(verbose: bool = True, half_width: float | None = None):
    par = Params(eps=(-2.0, 0.0, 3.0), gam=(1.0, 0.8, 1.2),
                 a=(-1.0, 0.5, 2.0), x=0)
    geo = Geometry(par)
    x = par.x

    # window intervals from the standard segmentation
    seg = segmentation(geo, T=60.0, half_width=half_width)
    wins = seg["windows"]                     # u-ordered, augmented with interval
    wp = window_pairs(geo)                    # u-ordered: crossed pair / spectator

    # theoretical jump probability per crossing pair:
    # Gamma_ij = gam_i^2 gam_j^2 |a_i-a_j| / (eps_i-eps_j)^2
    eps, gam, a = np.array(par.eps), np.array(par.gam), np.array(par.a)

    def Gamma_ij(i, j):
        return (gam[i]**2 * gam[j]**2 * abs(a[i]-a[j])
                / (eps[i]-eps[j])**2)

    results = []
    for w_idx in range(2):
        l, r = wins[w_idx]["interval"]
        pair = wp[w_idx]["pair"]
        spectator = wp[w_idx]["spectator"]

        W = propagate_ad_ip(geo, l, r, x, rtol=1e-13, atol=1e-14)

        # candidate algebraic frames at the two endpoints
        Vl = cauchy_frame(geo, l)
        Vr = cauchy_frame(geo, r)
        El = adiabatic_eigvecs(geo, l)
        Er = adiabatic_eigvecs(geo, r)
        I3 = np.eye(3, dtype=complex)

        frames = {
            "V(l)": Vl, "V(r)": Vr,
            "V(l)^-1": np.linalg.inv(Vl), "V(r)^-1": np.linalg.inv(Vr),
            "Eig(l)": El, "Eig(r)": Er,
            "I": I3,
        }

        # test  F_out^{-1} . W . F_in   for all combos (this conjugation
        # sends the diabatic/sheet structure of W into the F-basis).
        # Also the program form F_in^{-1} . W . F_out.
        trials = []
        for fin_name, Fin in frames.items():
            for fout_name, Fout in frames.items():
                # form A:  Fout^{-1} W Fin
                try:
                    MA = np.linalg.solve(Fout, W @ Fin)
                    trials.append((f"Fout^-1 . W . Fin  "
                                   f"[Fin={fin_name}, Fout={fout_name}]", MA))
                except np.linalg.LinAlgError:
                    pass
                # form B:  Fin^{-1} W Fout  (program's stated convention)
                try:
                    MB = np.linalg.solve(Fin, W @ Fout)
                    trials.append((f"Fin^-1 . W . Fout  "
                                   f"[Fin={fin_name}, Fout={fout_name}]", MB))
                except np.linalg.LinAlgError:
                    pass

        # raw W itself
        trials.append(("W (raw)", W))

        scored = []
        for label, M in trials:
            scored.append((block_offdiag_norm(M, pair), label, M))
        scored.sort(key=lambda t: t[0])

        Gam = Gamma_ij(*pair)
        p_theory = np.exp(-2*np.pi*Gam)

        results.append({
            "w_idx": w_idx, "interval": (l, r), "pair": pair,
            "spectator": spectator, "W": W,
            "scored": scored, "Gamma": Gam, "p_theory": p_theory,
        })

        if verbose:
            print(f"\n=== Window {w_idx}  u in [{l:.4f},{r:.4f}]  "
                  f"crossing pair {pair}  spectator {spectator} ===")
            print(f"  Gamma_ij = {Gam:.6e}   exp(-2 pi Gamma) = {p_theory:.6e}")
            print(f"  |W| =\n{np.abs(W)}")
            print(f"  best 6 frame maps by block-off-diagonal residual:")
            for resid, label, M in scored[:6]:
                print(f"    {resid:.3e}   {label}")
            # examine the best
            resid, label, M = scored[0]
            print(f"  --- BEST: {label}  residual={resid:.3e} ---")
            blk = block_2x2(M, pair)
            print(f"  2x2 block on pair {pair}:\n{blk}")
            print(f"  |2x2 block| =\n{np.abs(blk)}")
            # is it a rotation?  check |det|, columns orthonormal
            detb = np.linalg.det(blk)
            gram = blk.conj().T @ blk
            print(f"  det(block) = {detb:.6f}   |det| = {abs(detb):.6f}")
            print(f"  block^H block =\n{gram}")
            # extract jump probability from off-diagonal magnitude
            p_block = abs(blk[1, 0])**2 if blk.shape == (2, 2) else np.nan
            p_block2 = abs(blk[0, 1])**2
            print(f"  |block[1,0]|^2 = {p_block:.6e}   "
                  f"|block[0,1]|^2 = {p_block2:.6e}")
            print(f"  vs exp(-2 pi Gamma) = {p_theory:.6e}")

    # also: try conjugating the FULL window factor consistently --
    # does V(r)^{-1} W V(l) (transport from l-frame to r-frame) decouple?
    if verbose:
        print("\n\n=== summary: transport-consistent conjugation "
              "V(r)^-1 . W . V(l) ===")
        for res in results:
            w_idx = res["w_idx"]
            l, r = res["interval"]
            pair = res["pair"]
            W = res["W"]
            Vl = cauchy_frame(geo, l)
            Vr = cauchy_frame(geo, r)
            M = np.linalg.solve(Vr, W @ Vl)
            print(f"  window {w_idx}: block-off-diag residual "
                  f"= {block_offdiag_norm(M, pair):.3e}")

    return geo, results


def width_scan(verbose: bool = True):
    """
    Does the best frame-map block-off-diagonal residual VANISH as the
    window widens to contain the whole crossing?  If a genuine algebraic
    frame factorization exists, the residual -> 0 for wide windows.
    """
    par = Params(eps=(-2.0, 0.0, 3.0), gam=(1.0, 0.8, 1.2),
                 a=(-1.0, 0.5, 2.0), x=0)
    geo = Geometry(par)
    x = par.x
    wp = window_pairs(geo)
    centers = sorted(w["u_center"] for w in wp)
    gap = centers[1] - centers[0]

    print("\n\n############ WIDTH SCAN ############")
    print(f"window centers = {centers}   gap = {gap:.4f}")
    for hw_frac in (0.1, 0.2, 0.3, 0.4, 0.48):
        hw = hw_frac * gap
        seg = segmentation(geo, T=60.0, half_width=hw)
        wins = seg["windows"]
        print(f"\n--- half_width = {hw:.4f}  ({hw_frac} x gap) ---")
        for w_idx in range(2):
            l, r = wins[w_idx]["interval"]
            pair = wp[w_idx]["pair"]
            W = propagate_ad_ip(geo, l, r, x, rtol=1e-13, atol=1e-14)
            Vl = cauchy_frame(geo, l)
            Vr = cauchy_frame(geo, r)
            El = adiabatic_eigvecs(geo, l)
            Er = adiabatic_eigvecs(geo, r)
            cands = {
                "raw W": W,
                "Vr^-1 W Vl": np.linalg.solve(Vr, W @ Vl),
                "Vl^-1 W Vr": np.linalg.solve(Vl, W @ Vr),
                "Er^-1 W El": np.linalg.solve(Er, W @ El),
                "El^-1 W Er": np.linalg.solve(El, W @ Er),
                "Vr^H W Vl": Vr.conj().T @ W @ Vl,
            }
            best = min((block_offdiag_norm(M, pair), name)
                       for name, M in cands.items())
            allres = "  ".join(f"{n}:{block_offdiag_norm(M,pair):.2e}"
                               for n, M in cands.items())
            print(f"  win{w_idx} pair{pair}: best={best[1]} "
                  f"({best[0]:.3e})")
            print(f"      {allres}")


def wide_window_test(verbose: bool = True):
    """
    The decisive test.  Standard `segmentation` forces NARROW disjoint
    windows; the avoided crossings physically overlap, so the true window
    spans a wide u-range.  Split the real axis at the midpoint between the
    two crossing centres and ask whether ANY algebraic frame map decouples
    the spectator sheet of each (now full) window propagator.
    """
    par = Params(eps=(-2.0, 0.0, 3.0), gam=(1.0, 0.8, 1.2),
                 a=(-1.0, 0.5, 2.0), x=0)
    geo = Geometry(par)
    x = par.x
    wp = window_pairs(geo)
    centers = sorted(w["u_center"] for w in wp)
    mid = 0.5 * (centers[0] + centers[1])
    T = 60.0
    intervals = [(-T, mid), (mid, T)]
    eps, gam, a = np.array(par.eps), np.array(par.gam), np.array(par.a)

    def Gamma_ij(i, j):
        return gam[i]**2 * gam[j]**2 * abs(a[i]-a[j]) / (eps[i]-eps[j])**2

    print("\n\n############ WIDE-WINDOW DECISIVE TEST ############")
    for w_idx, (l, r) in enumerate(intervals):
        pair = wp[w_idx]["pair"]
        spec = wp[w_idx]["spectator"]
        W = propagate_ad_ip(geo, l, r, x, rtol=1e-13, atol=1e-14)
        G = Gamma_ij(*pair)
        p_th = np.exp(-2*np.pi*G)
        # algebraic frames just inside the (otherwise singular) endpoints
        Le = l if abs(l) < T else -40.0
        Re = r if abs(r) < T else 40.0
        Vl, Vr = cauchy_frame(geo, Le), cauchy_frame(geo, Re)
        El, Er = adiabatic_eigvecs(geo, Le), adiabatic_eigvecs(geo, Re)
        cands = {
            "raw W": W,
            "Vr^-1 W Vl": np.linalg.solve(Vr, W @ Vl),
            "Vl^-1 W Vr": np.linalg.solve(Vl, W @ Vr),
            "Vr^-1 W Vr": np.linalg.solve(Vr, W @ Vr),
            "Vl^-1 W Vl": np.linalg.solve(Vl, W @ Vl),
            "Er^-1 W El": np.linalg.solve(Er, W @ El),
            "El^-1 W Er": np.linalg.solve(El, W @ Er),
        }
        print(f"\n== window {w_idx}  pair {pair}  spectator {spec} ==")
        print(f"   Gamma_ij={G:.4f}  exp(-2 pi Gamma)={p_th:.4f}")
        print(f"   |W| =\n{np.round(np.abs(W),5)}")
        for n, M in cands.items():
            print(f"   {n:14s} block-off-diag = "
                  f"{block_offdiag_norm(M, pair):.4e}")
        # best blend of the two Cauchy frames (algebraic-frame hypothesis)
        best = np.inf
        for a1 in np.linspace(-1, 2, 31):
            for a2 in np.linspace(-1, 2, 31):
                Fin = a1*Vl + (1-a1)*Vr
                Fout = a2*Vl + (1-a2)*Vr
                try:
                    M = np.linalg.solve(Fin, W @ Fout)
                except np.linalg.LinAlgError:
                    continue
                best = min(best, block_offdiag_norm(M, pair))
        print(f"   best blend of Cauchy frames: residual = {best:.4e}")


if __name__ == "__main__":
    run()
    width_scan()
    wide_window_test()
