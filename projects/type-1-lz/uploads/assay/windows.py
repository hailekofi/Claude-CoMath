"""
Window / edge bookkeeping for the no-active-selector product

    S = Pi_+ . D2 . W_X2 . D1 . W_X1 . D0 . Pi_-          (eq. 1)

with the segments, in real-axis path order from u = -infinity,

    D0 : -T      -> l_X1
    W1 : l_X1    -> r_X1     (projected Q4 window X1)
    D1 : r_X1    -> l_X2
    W2 : l_X2    -> r_X2     (projected Q4 window X2)
    D2 : r_X2    -> +T

The two projected Q4 windows are ordered along the real u-axis by the
real part of u(q) for their conjugate-pair Q4 roots q (this is the
projection that the propagation -- which runs in u -- actually sees).

Endpoint permutations (nontrivial_assay_program, appendix B):
    Pi_- = I_3 ,   Pi_+ = [[0,0,1],[1,0,0],[0,1,0]] .
"""

from __future__ import annotations
import numpy as np

from .geometry import Geometry

PI_MINUS = np.eye(3)
PI_PLUS = np.array([[0, 0, 1],
                    [1, 0, 0],
                    [0, 1, 0]], dtype=float)


def segmentation(geo: Geometry, T: float = 60.0,
                 half_width: float | None = None,
                 theta: float = 0.0) -> dict:
    """
    Build the ordered segment list for the no-selector product.

    Returns dict with:
      segments : list of (kind, a, b)   in path order, kind in {'D','W'}
      windows  : the two window dicts (u-ordered), each augmented with
                 interval [l,r]
      centers  : real-axis window centers (sorted ascending)
    """
    wins = geo.windows(theta=theta)
    # order the two windows along the real u-axis
    wins = sorted(wins, key=lambda w: w["u_center"])
    centers = [w["u_center"] for w in wins]

    gap = centers[1] - centers[0]
    if half_width is None:
        # widest disjoint windows that still leave room for edge D1
        half_width = 0.3 * gap
    if 2 * half_width >= gap:
        raise ValueError("windows overlap: shrink half_width")

    l1, r1 = centers[0] - half_width, centers[0] + half_width
    l2, r2 = centers[1] - half_width, centers[1] + half_width
    if not (-T < l1 and r2 < T):
        raise ValueError("windows fall outside [-T, T]")

    wins[0]["interval"] = (l1, r1)
    wins[1]["interval"] = (l2, r2)

    segments = [
        ("D", -T, l1),
        ("W", l1, r1),
        ("D", r1, l2),
        ("W", l2, r2),
        ("D", r2, T),
    ]
    return {"segments": segments, "windows": wins, "centers": centers,
            "half_width": half_width, "T": T}
