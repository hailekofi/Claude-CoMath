"""
Level 2: lateral / median full-window assay.

Each projected Q4 conjugate-pair window X is reconstructed as the median
of two independent lateral continuations:

    W_X = Med(U_X+, U_X-) = U_X+ ( (U_X-)^{-1} U_X+ )^{-1/2} ,

where U_X+- are the gauged adiabatic propagators from l_X to r_X along
contours gamma_X+- that arc through the upper / lower half u-plane,
passing on the far side of the upper / lower Q4 turning point of the
window.  Because the turning point sits between the arc and the real
segment, gamma_X+- are not homotopic to the real axis: U_X+ and U_X-
genuinely differ by the window Stokes automorphism, and the median
symmetrises them.

The edge factors D_e are independent real-axis adiabatic solves.
The product Pi_+ D2 W2 D1 W1 D0 is compared to the one-pass benchmark.
"""

from __future__ import annotations
import numpy as np
from scipy.linalg import sqrtm

from .geometry import Geometry
from .ip import propagate_ad_ip, propagate_ad_ip_contour
from .windows import segmentation


# ---------------------------------------------------------------------------
#  contour
# ---------------------------------------------------------------------------
def make_arc(l: float, r: float, H: float):
    """
    Smooth half-plane arc from real l to real r peaking at imaginary part H:
        u(s)  = l + s (r-l) + i H sin(pi s)
        du/ds = (r-l)        + i H pi cos(pi s)
    H > 0 -> upper arc, H < 0 -> lower arc.
    """
    span = r - l

    def contour(s):
        u = l + s * span + 1j * H * np.sin(np.pi * s)
        duds = span + 1j * H * np.pi * np.cos(np.pi * s)
        return u, duds
    return contour


def singularity_u_values(geo: Geometry) -> np.ndarray:
    """u-images of the Q4 and W4 roots -- the branch points of A(u)."""
    sing = []
    for lam in np.concatenate([geo.Q4_roots(), geo.W4_roots()]):
        sing.append(np.polyval(geo.n, lam) / np.polyval(geo.p, lam))
    return np.array(sing)


def arc_clearance(contour, singents, n=200) -> float:
    """Minimum distance from the sampled arc to any singularity."""
    ss = np.linspace(0, 1, n)
    pts = np.array([contour(s)[0] for s in ss])
    d = np.min(np.abs(pts[:, None] - singents[None, :]))
    return float(d)


# ---------------------------------------------------------------------------
#  median
# ---------------------------------------------------------------------------
def median(Uplus: np.ndarray, Uminus: np.ndarray):
    """
    W = Med(U+,U-) = U+ ((U-)^{-1} U+)^{-1/2}, fixed (principal) branch.
    Returns (W, stokes_automorphism, sqrt_branch_ok).
    """
    Sigma = np.linalg.solve(Uminus, Uplus)         # (U-)^{-1} U+
    Sig_sqrt = np.asarray(sqrtm(Sigma), dtype=complex)   # sqrtm -> complex256
    branch_ok = bool(np.max(np.abs(Sig_sqrt @ Sig_sqrt - Sigma)) < 1e-8)
    W = Uplus @ np.linalg.inv(Sig_sqrt)
    return W, Sigma, branch_ok


# ---------------------------------------------------------------------------
#  per-window lateral construction
# ---------------------------------------------------------------------------
def window_factor(geo: Geometry, win: dict, x: int,
                   height_factor: float = 1.3) -> dict:
    """
    Build U_X+, U_X-, and the median window factor W_X for one window.
    `win` is a window dict from windows.segmentation (has 'interval',
    'u_vals').
    """
    l, r = win["interval"]
    h = float(np.max(np.abs(np.imag(win["u_vals"]))))   # turning-point height
    H = height_factor * h

    sing = singularity_u_values(geo)
    # own turning points are *inside* the arc by construction; clearance is
    # measured against the other six branch points
    own = win["u_vals"]
    others = np.array([s for s in sing
                       if min(abs(s - own[0]), abs(s - own[1])) > 1e-6])

    arc_up = make_arc(l, r, +H)
    arc_dn = make_arc(l, r, -H)
    clear_up = arc_clearance(arc_up, others)
    clear_dn = arc_clearance(arc_dn, others)

    Uplus = propagate_ad_ip_contour(geo, arc_up, x)
    Uminus = propagate_ad_ip_contour(geo, arc_dn, x)
    W, Sigma, branch_ok = median(Uplus, Uminus)

    # real-axis window propagator, for reference
    W_real = propagate_ad_ip(geo, l, r, x)

    return {
        "interval": (l, r), "H": H, "turn_height": h,
        "clearance": (clear_up, clear_dn),
        "U_plus": Uplus, "U_minus": Uminus,
        "W": W, "W_real": W_real,
        "stokes_automorphism": Sigma, "branch_ok": branch_ok,
    }


# ---------------------------------------------------------------------------
#  simplicity diagnostics
# ---------------------------------------------------------------------------
def simplicity_diagnostics(wf: dict) -> dict:
    Up, Um, W = wf["U_plus"], wf["U_minus"], wf["W"]
    Sigma = wf["stokes_automorphism"]
    # conjugation symmetry  U_- ~ conj(U_+)
    conj_sym = float(np.max(np.abs(Um - np.conj(Up))))
    # Stokes automorphism conditioning
    cond_Sigma = float(np.linalg.cond(Sigma))
    eig_Sigma = np.linalg.eigvals(Sigma)
    # structural zeros in W
    min_abs_entry = float(np.min(np.abs(W)))
    # median vs real-axis window propagator
    med_vs_real = float(np.max(np.abs(W - wf["W_real"])))
    # "dimensionality" of log W: singular values of logm
    from scipy.linalg import logm
    LW = np.asarray(logm(W), dtype=complex)
    sv = np.sort(np.abs(np.linalg.eigvals(LW)))[::-1]
    return {
        "conj_symmetry": conj_sym,
        "cond_stokes_automorphism": cond_Sigma,
        "stokes_eigvals": eig_Sigma,
        "min_abs_W_entry": min_abs_entry,
        "median_vs_realaxis": med_vs_real,
        "logW_eig_magnitudes": sv,
        "branch_ok": wf["branch_ok"],
    }


# ---------------------------------------------------------------------------
#  full Level 2 assay
# ---------------------------------------------------------------------------
def level2_assay(geo: Geometry, T: float = 50.0, x: int | None = None,
                 verbose: bool = True) -> dict:
    if x is None:
        x = geo.par.x
    seg = segmentation(geo, T=T)
    wins = seg["windows"]

    # window factors
    wfs = [window_factor(geo, w, x) for w in wins]
    diags = [simplicity_diagnostics(wf) for wf in wfs]

    # edge factors: independent real-axis adiabatic solves over the D segments
    edges = {}
    for kind, a, b in seg["segments"]:
        if kind == "D":
            edges[(a, b)] = propagate_ad_ip(geo, a, b, x)

    # assemble product in path order  D2 W2 D1 W1 D0
    M = np.eye(3, dtype=complex)
    wi = 0
    for kind, a, b in seg["segments"]:
        F = edges[(a, b)] if kind == "D" else wfs[wi]["W"]
        if kind == "W":
            wi += 1
        M = F @ M
    P_prod = np.abs(M.T) ** 2

    # benchmark: one-pass adiabatic
    M_bench = propagate_ad_ip(geo, -T, T, x)
    P_direct = np.abs(M_bench.T) ** 2
    bench_agreement = float(np.max(np.abs(P_prod - P_direct)))

    passed = bench_agreement < 1e-8 and all(d["branch_ok"] for d in diags)

    out = {"T": T, "x": x, "benchmark_agreement": bench_agreement,
           "passed": bool(passed), "window_factors": wfs,
           "diagnostics": diags, "P_prod": P_prod, "P_direct": P_direct}

    if verbose:
        print(f"--- Level 2 assay  (T={T}, x={x}) ---")
        for k, (wf, d) in enumerate(zip(wfs, diags)):
            print(f"  window X{k+1}  interval={wf['interval']}  "
                  f"H={wf['H']:.3f}  clearance={wf['clearance']}")
            print(f"    median vs real-axis window : {d['median_vs_realaxis']:.3e}")
            print(f"    conjugation symmetry       : {d['conj_symmetry']:.3e}")
            print(f"    Stokes automorphism cond   : {d['cond_stokes_automorphism']:.3e}")
            print(f"    min |W entry|              : {d['min_abs_W_entry']:.3e}")
            print(f"    |eig(log W)|               : {d['logW_eig_magnitudes']}")
            print(f"    sqrt branch ok             : {d['branch_ok']}")
        print(f"  product vs benchmark P : {bench_agreement:.3e}")
        print(f"  PASS: {passed}")
    return out
