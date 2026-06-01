"""
Closed-form transition matrix: the multistate Landau-Zener amplitude product.

The no-selector chamber is a Demkov-Osherov structure -- level 0 crosses
levels 1 and 2 at the two projected Q4 windows.  Dropping the semi-infinite
edge phases (diagonal -> drop from P = |S|^2), the scattering amplitude is

    S ~ Perm . N_B . D_mid . N_A

with
  N_A = 2-level LZ rotation of the window-A level pair, embedded 2 (+) 1,
  N_B = likewise for window B,
  D_mid = diag(e^{i chi_0}, e^{i chi_1}, e^{i chi_2}) the adiabatic phase
          between the two crossings (plus absorbed Stokes phases).

Each LZ rotation is a real rotation R(theta), sin^2(theta_X) = p_X, with
p_X = exp(-|Im I_X|) the elliptic-period jump probability (actions.py).
All non-trivial phases collect into D_mid; one of its three entries is a
gauge, leaving two physical edge phases.

P[x->j] = |(Perm N_B D_mid N_A)[j,x]|^2 .
"""

from __future__ import annotations
import numpy as np
from itertools import permutations
from scipy.optimize import least_squares

from .geometry import Geometry
from .actions import all_window_actions
from .ip import propagate_ad_ip


# ---------------------------------------------------------------------------
#  building blocks
# ---------------------------------------------------------------------------
def lz_rotation(p: float) -> np.ndarray:
    """2x2 real Landau-Zener rotation, jump probability p."""
    c, s = np.sqrt(1.0 - p), np.sqrt(p)
    return np.array([[c, -s], [s, c]])


def embed(M2: np.ndarray, pair: tuple) -> np.ndarray:
    """Embed a 2x2 block into 3x3 acting on `pair`, identity on the spectator."""
    i, j = pair
    sp = ({0, 1, 2} - {i, j}).pop()
    M = np.zeros((3, 3), complex)
    M[i, i], M[i, j] = M2[0, 0], M2[0, 1]
    M[j, i], M[j, j] = M2[1, 0], M2[1, 1]
    M[sp, sp] = 1.0
    return M


def scattering_model(pA: float, pB: float, pairA: tuple, pairB: tuple,
                     chi: np.ndarray, perm: tuple) -> np.ndarray:
    """S_model = Perm . N_B . D_mid . N_A   (3x3 amplitude matrix)."""
    NA = embed(lz_rotation(pA), pairA)
    NB = embed(lz_rotation(pB), pairB)
    Dmid = np.diag(np.exp(1j * np.asarray(chi)))
    M = NB @ Dmid @ NA
    P = np.eye(3)[list(perm)]                     # permutation matrix
    return P @ M


def model_P(pA, pB, pairA, pairB, chi, perm) -> np.ndarray:
    S = scattering_model(pA, pB, pairA, pairB, chi, perm)
    return np.abs(S.T) ** 2                        # P[x,j] = |S[j,x]|^2


# ---------------------------------------------------------------------------
#  identify which level pair each window crosses
# ---------------------------------------------------------------------------
def window_pairs(geo: Geometry) -> list:
    """
    For each window, the crossed level pair.  A Q4 root q lies on sheet i
    when Re(q) is in the interval that hosts lam_i; the spectator is i and
    the crossing pair is the complementary two.  Returned in u-order
    (crossing A = smaller u_center first).
    """
    e0, e1, e2 = geo.eps
    wins = sorted(geo.windows(), key=lambda w: w["u_center"])
    out = []
    for w in wins:
        re = float(np.real(w["roots"][0]))
        if e0 < re < e1:
            spectator = 1
        elif e1 < re < e2:
            spectator = 2
        else:
            spectator = 0
        pair = tuple(sorted({0, 1, 2} - {spectator}))
        out.append({"window": w, "pair": pair, "spectator": spectator,
                    "u_center": w["u_center"]})
    return out


# ---------------------------------------------------------------------------
#  fit the two edge phases (structure test)
# ---------------------------------------------------------------------------
def fit_closed_form(geo: Geometry, x: int = 0, T: float = 240.0,
                    verbose: bool = True) -> dict:
    # benchmark, T-converged
    P_bench = np.abs(propagate_ad_ip(geo, -T, T, x, rtol=1e-13, atol=1e-14).T) ** 2

    # window data: jump probabilities and crossed pairs, in u-order
    wp = window_pairs(geo)
    acts = {id(w["window"]): None for w in wp}
    # recompute actions in matching order
    wins_uorder = [w["window"] for w in wp]
    from .actions import window_action
    p = [float(np.exp(-abs(window_action(geo, w)["I_X"].imag)))
         for w in wins_uorder]
    pA, pB = p[0], p[1]
    pairA, pairB = wp[0]["pair"], wp[1]["pair"]

    best = None
    for perm in permutations(range(3)):
        def resid(chi2):
            chi = np.array([0.0, chi2[0], chi2[1]])
            return (model_P(pA, pB, pairA, pairB, chi, perm)
                    - P_bench).ravel()
        # multi-start over the two phases
        for c1 in np.linspace(0, 2*np.pi, 5, endpoint=False):
            for c2 in np.linspace(0, 2*np.pi, 5, endpoint=False):
                sol = least_squares(resid, [c1, c2], method="lm")
                err = float(np.max(np.abs(resid(sol.x))))
                if best is None or err < best["err"]:
                    best = {"err": err, "perm": perm,
                            "chi": (0.0, float(sol.x[0]), float(sol.x[1]))}

    out = {"pA": pA, "pB": pB, "pairA": pairA, "pairB": pairB,
           "P_bench": P_bench, **best}
    out["P_model"] = model_P(pA, pB, pairA, pairB,
                             np.array(best["chi"]), best["perm"])
    if verbose:
        print(f"--- closed-form structure fit  (x={x}) ---")
        print(f"  crossing A: pair {pairA}  p_A={pA:.6f}   "
              f"crossing B: pair {pairB}  p_B={pB:.6f}")
        print(f"  best perm={best['perm']}  edge phases chi={out['chi']}")
        print(f"  P_bench=\n{P_bench}")
        print(f"  P_model=\n{out['P_model']}")
        print(f"  max|P_model - P_bench| = {best['err']:.3e}")
    return out
