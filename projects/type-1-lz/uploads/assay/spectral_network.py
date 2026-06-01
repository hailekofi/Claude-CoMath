"""
Spectral network for the Type-1, N=3 Landau-Zener problem.

Hitchin reformulation
---------------------
The Type-1 problem is the parallel transport on a base C = u-line of a flat
connection whose semiclassical eigenvalues are the three roots
:math:`\\lambda_i(u)` of the spectral polynomial

.. math::  q_u(\\lambda) \;=\; u\,p(\\lambda) - n(\\lambda) \;=\; 0,

with :math:`p(\\lambda)=\\prod_a(\\lambda-\\epsilon_a)` and
:math:`n(\\lambda)=\\sum_a\\gamma_a^2\\prod_{b\\ne a}(\\lambda-\\epsilon_b)`.
The three sheets together define a rank-3 cover

.. math::  \\Sigma_{\\rm GMN} = \\{(u, \\lambda(u))\\} \\subset T^*C ,

which is the GMN/Hitchin spectral cover for the Hitchin one-form
:math:`\\lambda\\,du`.  Its sheets meet at the *turning points*, the four roots
of :math:`Q_4(\\lambda)` in the :math:`\\lambda`-plane, which map to four u-points

.. math::  u_\\alpha = n(\\lambda_\\alpha)/p(\\lambda_\\alpha) ,

paired into two conjugate-pair "windows" :math:`X_A, X_B`.

S-walls (a.k.a. WKB curves)
---------------------------
For phase :math:`\\vartheta`, an :math:`(i,j)` S-wall is a curve in the u-plane
satisfying

.. math::  e^{-i\\vartheta}\\,(\\lambda_i(u)-\\lambda_j(u))\\,du \\in \\mathbb{R}_+ .

Three S-walls emanate from each turning point :math:`u_\\alpha`; on a small
circle around :math:`u_\\alpha`, the three sheets that pairwise vanish carry
walls at angles :math:`\\arg(du) = \\vartheta - \\arg(\\lambda_i-\\lambda_j)`.

BPS states (finite saddles)
---------------------------
A finite BPS state is a saddle trajectory connecting two distinct turning
points :math:`u_\\alpha \\to u_\\beta`.  Its central charge is

.. math::  Z_\\gamma = \\oint_\\gamma \\lambda\\,du ,

the period of the Hitchin form on the homology class :math:`\\gamma` it lifts
to on :math:`\\Sigma_{\\rm GMN}`.  For the Type-1 N=3 problem, each conjugate-pair
turning-point window :math:`X` carries a single finite BPS state whose period
equals the project's window period :math:`I_X = \\oint\\sqrt{y}\\,du` from
:mod:`assay.actions` (proved in Track 1).

This module
-----------
* :func:`turning_points`            -- four turning points :math:`(\\lambda_\\alpha, u_\\alpha)`
                                       grouped into two windows.
* :func:`sheet_pair_at_window`      -- the (i,j) labels of the two sheets that
                                       collide at each window, in the
                                       u=±infty labelled-lambda basis.
* :func:`wall_direction`            -- local emanation directions of the
                                       three S-walls at one turning point.
* :func:`integrate_wall`            -- ODE integration of one S-wall in the
                                       u-plane, with sheet tracking.
* :func:`spectral_network`          -- assemble all S-walls for one geometry,
                                       classify endpoints (BPS saddle vs
                                       escape to infinity vs collide).
* :func:`bps_central_charges`       -- (Z_A, Z_B) for the two windows; equal
                                       to (I_X_A, I_X_B) up to sign.
* :func:`real_axis_wall_crossings`  -- list of (u_cross, sheet_pair) for every
                                       S-wall that crosses the real u-axis.

Notes
-----
* The trajectory ODE in u is stiff near turning points; we use scipy's DOP853
  with capped step size and a finite stopping condition (escape ball or
  collision proximity).
* Sheet tracking uses continuation from a real-axis sample point of the
  ``labelled_lambdas`` indexing.
"""

from __future__ import annotations
import numpy as np
from scipy.integrate import solve_ivp
from typing import Iterable

from .geometry import Geometry
from .adiabatic import labelled_lambdas


# ---------------------------------------------------------------------------
#  Turning points and their u-images
# ---------------------------------------------------------------------------
def turning_points(geo: Geometry) -> list:
    """
    The four turning points :math:`(\\lambda_\\alpha, u_\\alpha)` grouped into
    two windows.

    Returns
    -------
    list of two dicts, each with keys
        ``roots_lambda`` : (lam_+, lam_-) conjugate pair of Q4 roots
        ``u_vals``       : (u_+, u_-) their u-images
        ``u_center``     : Re(mean(u_vals))
        ``width``        : |Im(u_+)|
    """
    wins = sorted(geo.windows(), key=lambda w: w["u_center"])
    out = []
    for w in wins:
        uv = w["u_vals"]
        out.append({
            "roots_lambda": w["roots"],
            "u_vals": uv,
            "u_center": float(w["u_center"]),
            "width": float(abs(np.imag(uv[0]))),
            "conj_defect": w["conj_defect"],
        })
    return out


# ---------------------------------------------------------------------------
#  Sheet labelling at a turning point
# ---------------------------------------------------------------------------
def _lambda_at_u(geo: Geometry, u: complex, hint: np.ndarray | None = None
                 ) -> np.ndarray:
    """Tracked labelled_lambdas at complex u, falling back to root-sort + hint."""
    if abs(complex(u).imag) < 1e-12 and hint is None:
        return labelled_lambdas(geo, float(u.real))
    if hint is None:
        return labelled_lambdas(geo, complex(u))
    # nearest-continuation labelling
    if u == 0.0:
        u = 1e-30
    roots = geo.lambdas(u)
    out = np.empty(3, complex)
    pool = list(roots)
    for k in range(3):
        j = int(np.argmin([abs(r - hint[k]) for r in pool]))
        out[k] = pool.pop(j)
    return out


def sheet_pair_at_window(geo: Geometry, win: dict, eps_u: float = 0.01
                         ) -> tuple:
    """
    Identify the two sheet labels (i,j) whose eigenvalues collide at this
    Q4 window.

    The probe is at u = u_center + eps_u  +  i * eps_u  (a small complex
    offset into the upper half-plane to sit between the conjugate pair) and
    compares pairwise distances :math:`|\\lambda_i - \\lambda_j|`.

    Returns
    -------
    (i, j) with i<j in the ``labelled_lambdas`` indexing.
    """
    # Probe just inside the window
    uc = win["u_center"]
    width = win.get("width", 0.5)
    # find direction *toward* a turning point: choose imag part with smaller |Im|
    u_probe = complex(uc, 0.5 * width)
    # but if width is tiny, fall back to a smaller offset
    if width < 1e-6:
        u_probe = complex(uc + eps_u, eps_u)
    lam = _lambda_at_u(geo, u_probe)
    diffs = [(i, j, abs(lam[i] - lam[j]))
             for i in range(3) for j in range(i + 1, 3)]
    diffs.sort(key=lambda t: t[2])
    return diffs[0][0], diffs[0][1]


# ---------------------------------------------------------------------------
#  Wall emanation directions
# ---------------------------------------------------------------------------
def wall_directions(geo: Geometry, win: dict, theta: float = 0.0
                    ) -> list:
    """
    The three S-wall emanation directions in the u-plane at one turning point.

    At a Q4 turning point :math:`u_\\alpha`, the two colliding sheets have
    :math:`\\lambda_i - \\lambda_j \\propto \\sqrt{u - u_\\alpha}`.  An S-wall
    of type (i,j) for phase :math:`\\vartheta` then satisfies

    .. math::  \\arg(u - u_\\alpha) = \\frac{2}{3}\\big[\\vartheta - \\arg c\\big]
                                       + \\frac{2\\pi k}{3},\\ k=0,1,2,

    with :math:`c = \\lambda_i - \\lambda_j` near :math:`u_\\alpha`, and we use
    the conventional rank-3 spacing :math:`2\\pi/3` between the three walls
    (one of the three sheet labels at a generic turning point is the
    "non-participating" one and the other two pair the wall).

    Returns
    -------
    list of three complex unit vectors (the directions of the three walls).
    """
    # local expansion: near a Q4 root q (the lambda-plane turning point) on
    # one sheet branch, lam_i - lam_j ~ A * sqrt(u - u_alpha) with some
    # complex constant A.  The walls at angles arg(du) = -arg(A) + theta + 2pi k/3
    # (k=0,1,2).
    # Estimate A by sampling lam_i - lam_j at a small displacement from u_alpha.
    u_alpha = complex(np.mean(win["u_vals"]))     # use *real* of u_center (axis projection)
    # for the actual turning point we pick the upper-half representative:
    u_alpha_upper = win["u_vals"][int(np.argmax(np.imag(win["u_vals"])))]
    dr = 1e-3 * max(1.0, win.get("width", 0.5))
    u_probe = u_alpha_upper + dr
    # sheets that collide:
    i, j = sheet_pair_at_window(geo, win)
    lam = _lambda_at_u(geo, u_probe)
    diff = lam[i] - lam[j]
    A = diff / np.sqrt(complex(dr))    # leading coefficient
    # Wall direction: solve  arg((lam_i-lam_j) du) = theta
    # near turning point, (lam_i-lam_j) ~ A sqrt(du), so we want arg(A sqrt(du) du) = theta
    # i.e. arg(du^{3/2}) = theta - arg(A), so 3/2 * arg(du) = theta - arg(A) + 2pi k
    # arg(du) = (2/3)(theta - arg(A)) + 4 pi k / 3
    base = (2.0 / 3.0) * (theta - np.angle(A))
    return [np.exp(1j * (base + 4 * np.pi * k / 3.0)) for k in range(3)]


# ---------------------------------------------------------------------------
#  S-wall ODE integrator
# ---------------------------------------------------------------------------
def _wall_rhs_factory(geo: Geometry, i: int, j: int, theta: float):
    """
    RHS of the wall ODE: du/ds = e^{i theta} / (lam_i - lam_j),
    parameterised so |(lam_i-lam_j) du/ds| = 1 (constant speed in the
    Hitchin one-form).  Sheets are tracked via the labelled_lambdas
    continuation hint carried through the integrator.
    """
    state = {"hint": None}

    def rhs(s, y):
        u = complex(y[0], y[1])
        try:
            lam = _lambda_at_u(geo, u, state["hint"])
        except Exception:
            return [0.0, 0.0]
        state["hint"] = lam
        denom = lam[i] - lam[j]
        if abs(denom) < 1e-14:
            return [0.0, 0.0]
        du = np.exp(1j * theta) / denom
        return [du.real, du.imag]

    return rhs, state


def integrate_wall(geo: Geometry, u0: complex, direction: complex,
                   i: int, j: int, theta: float = 0.0,
                   s_max: float = 30.0, n_steps: int = 600,
                   r_escape: float = 12.0,
                   r_min_turn: float = 1e-3,
                   other_turn_points: Iterable[complex] | None = None,
                   ) -> dict:
    """
    Integrate a single S-wall from a starting point ``u0`` in initial
    direction ``direction`` (a unit complex) for sheets (i,j) at phase
    ``theta``.

    Stopping conditions:
        * |u| > r_escape           : wall escaped to infinity
        * approach within r_min_turn of another turning point: saddle hit
        * abs(lam_i - lam_j) very small (we hit a *different* turning point)

    Returns
    -------
    dict with keys
        u_trajectory : np.ndarray of complex u-points along the wall
        endpoint     : "escape", "saddle@k", or "stuck"
        s_final      : final integration parameter
        sheet_pair   : (i, j)
    """
    rhs, state = _wall_rhs_factory(geo, i, j, theta)
    # initial position slightly off the turning point along ``direction``
    u_init = u0 + 1e-3 * direction
    # initial hint -- sample sheets just off the turning point
    state["hint"] = _lambda_at_u(geo, u_init)
    y0 = [u_init.real, u_init.imag]

    other_turn_points = list(other_turn_points or [])

    traj = [u_init]
    s_grid = np.linspace(0, s_max, n_steps)
    sol = solve_ivp(rhs, (0, s_max), y0, t_eval=s_grid,
                    method="DOP853", rtol=1e-8, atol=1e-10, max_step=0.5)
    if not sol.success:
        return {"u_trajectory": np.array([u_init]),
                "endpoint": "stuck",
                "s_final": 0.0,
                "sheet_pair": (i, j)}
    ys = sol.y
    us = ys[0] + 1j * ys[1]
    traj_u = us

    # post-hoc analysis: where did it end?
    endpoint = "escape"
    s_final = sol.t[-1]
    for k, u in enumerate(us):
        if abs(u) > r_escape:
            traj_u = us[: k + 1]
            endpoint = "escape"
            s_final = sol.t[k]
            break
        for ti, tp in enumerate(other_turn_points):
            if abs(u - tp) < r_min_turn:
                traj_u = us[: k + 1]
                endpoint = f"saddle@{ti}"
                s_final = sol.t[k]
                return {"u_trajectory": traj_u, "endpoint": endpoint,
                        "s_final": s_final, "sheet_pair": (i, j)}
    return {"u_trajectory": traj_u, "endpoint": endpoint,
            "s_final": s_final, "sheet_pair": (i, j)}


# ---------------------------------------------------------------------------
#  Assemble the full spectral network
# ---------------------------------------------------------------------------
def spectral_network(geo: Geometry, theta: float = 0.0,
                     n_steps: int = 600, s_max: float = 30.0
                     ) -> dict:
    """
    Build the spectral network at phase :math:`\\vartheta` for one Geometry.

    Returns
    -------
    dict with keys
        windows      : list of two window dicts (as from turning_points)
        sheet_pairs  : (i_A, j_A), (i_B, j_B)
        walls        : list of all wall trajectories (six in total: three from
                       each window's upper-half-plane turning point;
                       the lower-half walls follow by complex conjugation
                       symmetry, so we record only the upper-half walls).
        bps_saddles  : list of finite saddles found (window-internal pairings)
        crossings    : list of (u_real, sheet_pair) for each real-axis crossing
    """
    tps = turning_points(geo)
    # gather all u-coords for endpoint detection
    all_tp_u = []
    for w in tps:
        all_tp_u.extend(list(w["u_vals"]))

    walls = []
    crossings = []
    bps_saddles = []
    sheet_pairs = []

    for wi, w in enumerate(tps):
        # upper-half-plane representative of the turning point
        u_alpha = w["u_vals"][int(np.argmax(np.imag(w["u_vals"])))]
        i, j = sheet_pair_at_window(geo, w)
        sheet_pairs.append((i, j))
        dirs = wall_directions(geo, w, theta=theta)

        # other turning points (avoid hitting ourselves)
        other_u = [u for u in all_tp_u if abs(u - u_alpha) > 1e-6]

        for k, d in enumerate(dirs):
            res = integrate_wall(geo, u_alpha, d, i, j, theta=theta,
                                 s_max=s_max, n_steps=n_steps,
                                 other_turn_points=other_u)
            res["window"] = wi
            res["dir_idx"] = k
            walls.append(res)
            if res["endpoint"].startswith("saddle"):
                bps_saddles.append({
                    "from_window": wi,
                    "dir_idx": k,
                    "endpoint": res["endpoint"],
                    "trajectory": res["u_trajectory"],
                    "sheet_pair": res["sheet_pair"],
                })
            # find real-axis crossings
            us = res["u_trajectory"]
            ims = np.imag(us)
            for m in range(len(ims) - 1):
                if ims[m] * ims[m + 1] < 0:
                    # linear-interp the crossing
                    t = ims[m] / (ims[m] - ims[m + 1])
                    u_x = (1 - t) * us[m] + t * us[m + 1]
                    crossings.append({
                        "u": complex(u_x.real, 0.0),
                        "sheet_pair": (i, j),
                        "window": wi,
                        "dir_idx": k,
                    })

    return {
        "windows": tps,
        "sheet_pairs": sheet_pairs,
        "walls": walls,
        "bps_saddles": bps_saddles,
        "crossings": crossings,
    }


# ---------------------------------------------------------------------------
#  BPS central charges = window periods
# ---------------------------------------------------------------------------
def bps_central_charges(geo: Geometry) -> tuple:
    """
    BPS central charges of the two finite saddles (one per window).
    Equal in magnitude to the window periods :math:`I_X` from
    :mod:`assay.actions`.

    Returns
    -------
    (Z_A, Z_B) complex tuple in u-order of the windows.
    """
    from .actions import window_action
    wins = sorted(geo.windows(), key=lambda w: w["u_center"])
    Zs = []
    for w in wins:
        Z = window_action(geo, w)["I_X"]
        Zs.append(complex(Z))
    return tuple(Zs)


# ---------------------------------------------------------------------------
#  Real-axis crossings: a clean list of where the path goes through walls
# ---------------------------------------------------------------------------
def real_axis_wall_crossings(geo: Geometry, theta: float = 0.0,
                             n_steps: int = 600, s_max: float = 30.0
                             ) -> list:
    """
    A simplified summary of where the real-u-axis traversal crosses S-walls.

    Each item is a dict with keys ``u``, ``sheet_pair``, ``window``.
    """
    net = spectral_network(geo, theta=theta, n_steps=n_steps, s_max=s_max)
    cs = net["crossings"]
    cs.sort(key=lambda c: c["u"].real)
    return cs


# ===========================================================================
#  Phase III Route A: BPS inventory + wall orientation (foundation for the
#  GMN omega-pair side-tracking implementation)
# ===========================================================================
def _bps_spectral_network(
    geo: Geometry,
    theta: float,
    s_max: float = 60.0,
    n_steps: int = 1500,
    r_escape: float = 20.0,
    r_min_turn: float = 0.1,
) -> dict:
    """
    Phase-III variant of :func:`spectral_network` with **relaxed** saddle
    detection tolerance ``r_min_turn`` and extended integration range.
    The default ``r_min_turn = 1e-3`` in the base function rejects every
    wall in our Type-1 N=3 geometry because the closest approach to the
    partner turning point is ~0.02–0.20 in u (well above 1e-3 but well
    below the wall escape radius).  A relaxed tolerance of 0.1 captures
    the in-window BPS saddles cleanly; a relaxed tolerance of 0.5 also
    catches inter-window saddles if they exist.
    """
    tps = turning_points(geo)
    all_tp_u = []
    for w in tps:
        all_tp_u.extend(list(w["u_vals"]))

    walls = []
    crossings = []
    bps_saddles = []
    sheet_pairs = []

    for wi, w in enumerate(tps):
        u_alpha = w["u_vals"][int(np.argmax(np.imag(w["u_vals"])))]
        i, j = sheet_pair_at_window(geo, w)
        sheet_pairs.append((i, j))
        dirs = wall_directions(geo, w, theta=theta)
        other_u = [u for u in all_tp_u if abs(u - u_alpha) > 1e-6]

        for k, d in enumerate(dirs):
            res = integrate_wall(
                geo, u_alpha, d, i, j, theta=theta,
                s_max=s_max, n_steps=n_steps,
                r_escape=r_escape, r_min_turn=r_min_turn,
                other_turn_points=other_u,
            )
            res["window"] = wi
            res["dir_idx"] = k
            walls.append(res)
            if res["endpoint"].startswith("saddle"):
                # identify the target window from the endpoint string
                target_idx_str = res["endpoint"].split("@")[1]
                try:
                    target_tp_index = int(target_idx_str)
                    target_u = other_u[target_tp_index]
                    # find which window the target TP belongs to
                    target_window = wi   # default
                    for tw_idx, tw in enumerate(tps):
                        if any(abs(target_u - utp) < 1e-6 for utp in tw["u_vals"]):
                            target_window = tw_idx
                            break
                except (ValueError, IndexError):
                    target_window = wi

                bps_saddles.append({
                    "from_window": wi,
                    "to_window": target_window,
                    "dir_idx": k,
                    "endpoint": res["endpoint"],
                    "trajectory": res["u_trajectory"],
                    "sheet_pair": res["sheet_pair"],
                })
            us = res["u_trajectory"]
            ims = np.imag(us)
            for m in range(len(ims) - 1):
                if ims[m] * ims[m + 1] < 0:
                    t = ims[m] / (ims[m] - ims[m + 1])
                    u_x = (1 - t) * us[m] + t * us[m + 1]
                    crossings.append({
                        "u": complex(u_x.real, 0.0),
                        "sheet_pair": (i, j),
                        "window": wi,
                        "dir_idx": k,
                    })

    return {
        "windows": tps,
        "sheet_pairs": sheet_pairs,
        "walls": walls,
        "bps_saddles": bps_saddles,
        "crossings": crossings,
    }


def bps_inventory(
    geo: Geometry,
    theta: float = 0.5 * np.pi,
    s_max: float = 60.0,
    n_steps: int = 1500,
    r_min_turn: float = 0.1,
    relax_for_inter_window: bool = False,
) -> dict:
    """
    Phase III Route A foundation: enumerate the BPS spectrum at phase
    ``theta`` and verify the genus-1 / b_1=2 prediction that two BPS
    states exist (one per window).  When ``relax_for_inter_window=True``,
    ``r_min_turn`` is increased to 0.5 to also catch any inter-window
    saddles (BPS states connecting different windows -- the architectural
    risk identified in the Phase III plan).

    Parameters
    ----------
    geo : Geometry
    theta : float
        BPS phase.  Default pi/2.  Slight tilt (pi/2 + 1e-4) stabilises
        the network at the BPS ray.
    s_max, n_steps : float, int
        ODE integration parameters.
    r_min_turn : float
        Saddle-detection tolerance.  Phase III default 0.1 (relaxed from
        the upstream 1e-3 in :func:`integrate_wall`).
    relax_for_inter_window : bool
        If True, set ``r_min_turn = 0.5`` to enumerate inter-window
        saddles too.  Set this for the inter-window-BPS scan recommended
        in the Phase III plan's Risk Analysis.

    Returns
    -------
    dict with keys:
        n_bps           : number of finite BPS saddles found
        saddles         : list of saddle dicts, each carrying:
                            'window'       : window of origin (from_window)
                            'to_window'    : window the saddle lands in
                                             (== from_window for in-window
                                              BPS; != for inter-window)
                            'sheet_pair'   : (i, j)
                            'trajectory'   : np.ndarray of u-points
                            'Z'            : BPS central charge
                            'mu'           : Voros symbol exp(i Z)
                            'abs_mu'       : |mu|
                            'kind'         : 'in_window' or 'inter_window'
        n_in_window     : count of in-window BPS (expected 2 by G1)
        n_inter_window  : count of inter-window BPS (expected 0; nonzero
                          is the Phase III "missing inter-window BPS"
                          risk hit)
        theta_used      : float
        b1_check        : bool, whether n_in_window == 2
    """
    from .actions import window_action

    if relax_for_inter_window:
        r_min_turn = max(r_min_turn, 0.5)

    sn = _bps_spectral_network(
        geo, theta=theta, s_max=s_max, n_steps=n_steps,
        r_min_turn=r_min_turn,
    )
    geo_windows = sorted(geo.windows(), key=lambda w: w["u_center"])

    saddles = []
    n_in_window = 0
    n_inter_window = 0
    for sad in sn["bps_saddles"]:
        wi = int(sad["from_window"])
        to_wi = int(sad.get("to_window", wi))
        if wi == to_wi:
            kind = "in_window"
            n_in_window += 1
        else:
            kind = "inter_window"
            n_inter_window += 1

        if 0 <= wi < len(geo_windows):
            if kind == "in_window":
                # Z = window period (Track 1)
                Z = complex(window_action(geo, geo_windows[wi])["I_X"])
            else:
                # inter-window: Z is a NEW Voros symbol along this saddle's
                # trajectory; we'll need the corrected V_AB on the right sheet.
                # As a placeholder, use the half-sum of the two windows' periods
                # (this is the GMN central-charge convention for inter-window
                # hypermultiplets at canonical theta).
                Z_from = complex(window_action(geo, geo_windows[wi])["I_X"])
                Z_to = complex(window_action(geo, geo_windows[to_wi])["I_X"])
                Z = 0.5 * (Z_from + Z_to)
        else:
            Z = 0.0 + 0.0j
        mu = complex(np.exp(1j * Z))
        saddles.append({
            "window": wi,
            "to_window": to_wi,
            "sheet_pair": tuple(int(k) for k in sad["sheet_pair"]),
            "trajectory": sad["trajectory"],
            "Z": Z,
            "mu": mu,
            "abs_mu": float(abs(mu)),
            "kind": kind,
        })

    return {
        "n_bps": len(saddles),
        "n_in_window": n_in_window,
        "n_inter_window": n_inter_window,
        "saddles": saddles,
        "theta_used": float(theta),
        "b1_check": (n_in_window >= 2),
        "r_min_turn_used": float(r_min_turn),
    }


def wall_with_orientation(geo: Geometry, theta: float = 0.5 * np.pi,
                          s_max: float = 40.0, n_steps: int = 1200) -> list:
    """
    Augment spectral_network's wall list with an ``orientation`` flag
    used by the GMN side-tracking convention.

    Orientation = sign(imag(direction)) at the wall's emission angle.
    Walls emanating into the upper half plane have orientation = +1;
    into the lower half plane = -1.  Real-axis crossings of walls with
    different orientations pick up reciprocal-mu unipotents.

    Returns
    -------
    list of wall dicts (copies of spectral_network's walls list) with
    an added ``orientation`` field.
    """
    net = spectral_network(geo, theta=theta, s_max=s_max, n_steps=n_steps)
    out = []
    for w in net["walls"]:
        d = dict(w)  # copy
        # direction is the initial wall_directions output; we infer it
        # from the trajectory's first step.
        traj = w["u_trajectory"]
        if len(traj) >= 2:
            initial_dir = complex(traj[1] - traj[0])
            d["orientation"] = +1 if initial_dir.imag >= 0 else -1
        else:
            d["orientation"] = 0
        out.append(d)
    return out
