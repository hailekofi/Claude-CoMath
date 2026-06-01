"""
Track 5 cross-link assay: phase-side (Q_4 window periods I_X) combined with
transport-side (W_4 root residuals S_{+-}) invariants, exposed for use by the
validation harness.

The phase one-form is

    omega = - sqrt(Q_4) * L_H * W_4 / p^3  d lambda

with the corresponding period over a Q_4 conjugate-pair window giving the
elliptic period I_X. The Brundobler-Elser exponents 2 pi Gamma_{ij} are
imaginary parts of these phase-side periods (assay/actions.py, Track 1).

The transport curve Sigma_R is defined by  nu^2 = W_4(lambda), with the
*transport* discriminant locus C_R the u-image of the W_4 zero set. Following
the primer / gluing-selector / exact-insertion-rule notes, the selector
criterion at parameter theta is

    Im( e^{-i theta} * S_{+-}(u_*; H) )  = 0,

with S_{+-}(u_*) the phase mismatch around the transport ramification at
u_* = u(v_*), v_* a W_4 root. Concretely, we integrate omega from a real
reference u_0 (set to 0) up to u_* along a chord in the lambda plane, and
report S_{+-} = int_{u_0}^{u_*} omega.

This module provides:

    cross_link_invariants(geo) -> dict
        Phase-side: window periods I_X (in u-order), |Im I_X|/(2 pi).
        Transport-side: per W_4 root residuals S_{+-}(v_alpha) (one per
        independent conjugate pair, picking Im(v) > 0 representatives),
        their |Re|, |Im| parts, and the corresponding contour-style residuals
        from integrating omega around the W_4 conjugate-pair loops
        (genuine transport-period candidates on Sigma_R).

    Hypothesis formulas (predict a 3x3 P from a Record using the cross-link
    invariants):
        formula_H_A   --  product q_ij * |Im S_{+-}(v_alpha)| factor
        formula_H_B   --  exp(- |Re S_{+-}|) * q_ij ansatz
        formula_H_C   --  Im of W_4 conjugate-pair residual itself
        formula_H_D   --  mixed-cover linear combination phase + transport
        formula_H_E   --  fallback: identity + best linear combo fit
"""

from __future__ import annotations
import numpy as np
from typing import Dict, List, Tuple
from itertools import permutations

from .geometry import Geometry, Params
from .actions import all_window_actions, _omega_lambda


# ---------------------------------------------------------------------------
#  primitive: integrate omega along a chord in lambda, with branch tracking
# ---------------------------------------------------------------------------
def _omega_chord_integral(geo: Geometry,
                          lam_start: complex,
                          lam_end: complex,
                          n_steps: int = 2000,
                          sign0: complex = None) -> complex:
    """
    Integrate omega = -sqrt(Q4) L_H W_4 / p^3 d lambda along the straight
    chord [lam_start, lam_end] in the lambda plane.

    sqrt(Q_4) is initialised at lam_start by sign0 (or by np.sqrt) and
    continued by nearest-value tracking thereafter. Used for transport-side
    residual S_{+-} = int_{lam_real_ref}^{lam_v} omega along a chord, where
    the chord avoids polynomial zeros if the start is real and the end is
    a W_4 root with sufficiently small |Im|.

    Returns the (complex) line integral.
    """
    ts = np.linspace(0.0, 1.0, n_steps)
    lam = lam_start + (lam_end - lam_start) * ts
    dl = (lam_end - lam_start) / (n_steps - 1)
    Q4v = np.polyval(geo.Q4, lam)
    s = sign0 if sign0 is not None else np.sqrt(Q4v[0])
    sq = np.empty(n_steps, complex)
    for k in range(n_steps):
        cand = np.sqrt(Q4v[k])
        s = cand if abs(cand - s) < abs(-cand - s) else -cand
        sq[k] = s
    integrand = _omega_lambda(geo, lam, sq)
    return float(0.0) + np.trapezoid(integrand, dx=dl)


def _w4_loop_period(geo: Geometry,
                    centre: complex,
                    radius: float,
                    n: int = 4000) -> complex:
    """
    Period of omega around a circle in lambda enclosing two W_4 branch points
    (a single conjugate pair). The integrand is meromorphic on the lambda
    plane (no branch from sqrt(Q4) on the loop -- it picks up sign +1 around
    the W4 pair only if the loop does not cross Q4 cuts; if the loop is small
    enough it does not enclose any Q4 branch point).

    Returns the complex period (could be 2 pi i * sum of W_4-residues for
    poles inside if any; for the integrand of omega the W_4-zeros are not
    poles, so this is a contour-period candidate of the "mixed" / transport
    style).
    """
    th = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    lam = centre + radius * np.exp(1j * th)
    dlam = 1j * radius * np.exp(1j * th) * (2.0 * np.pi / n)
    Q4v = np.polyval(geo.Q4, lam)
    s = np.sqrt(Q4v[0])
    sq = np.empty(n, complex)
    for k in range(n):
        cand = np.sqrt(Q4v[k])
        s = cand if abs(cand - s) < abs(-cand - s) else -cand
        sq[k] = s
    integrand = _omega_lambda(geo, lam, sq)
    return np.sum(integrand * dlam)


# ---------------------------------------------------------------------------
#  Main computation: cross-link invariants for one Geometry
# ---------------------------------------------------------------------------
def cross_link_invariants(geo: Geometry) -> dict:
    """
    Compute both phase-side (Q4 window periods) and transport-side (W4 root
    residuals) invariants for a single Geometry.

    Returns a dict with keys
      I_X        : tuple[complex, complex]  phase-side window periods (u-order)
      Gamma_X    : tuple[float, float]      |Im I_X| / (2 pi) for each window
      v_W4       : list[complex]            two independent W_4 roots (Im > 0)
      S_pm       : list[complex]            S_{+-}(v) = int_0^{v} omega along
                                            chord, branch-tracked
      |Re S_pm|  : list[float]
      |Im S_pm|  : list[float]
      I_W4       : list[complex]            W_4 loop periods (one per pair)
      Im_I_W4    : list[float]              |Im I_W4|
      Re_I_W4    : list[float]              |Re I_W4|
      u_W4       : list[complex]            u-images of the W_4 roots picked

    The two conjugate pairs of W_4 each contribute one representative root
    (Im > 0) for S_{+-} and one circular period for I_W4.
    """
    # phase-side
    wins = sorted(geo.windows(), key=lambda w: w["u_center"])
    acts = []
    for w in wins:
        from .actions import window_action
        acts.append(window_action(geo, w))
    IX = tuple(a["I_X"] for a in acts)
    GammaX = tuple(abs(I.imag) / (2.0 * np.pi) for I in IX)

    # transport-side: W_4 roots and pairs
    w4_roots = np.array(geo.W4_roots())
    # pair the four roots into two conjugate pairs by minimum |z - conj(w)|
    used = [False]*4
    pairs = []
    order = np.argsort(np.imag(w4_roots))  # negative-imag first
    for idx in order:
        if used[idx]:
            continue
        best, bestd = None, np.inf
        for j in range(4):
            if j == idx or used[j]:
                continue
            d = abs(w4_roots[idx] - np.conj(w4_roots[j]))
            if d < bestd:
                best, bestd = j, d
        if best is None:
            continue
        used[idx] = used[best] = True
        # pick the Im > 0 element as representative
        if w4_roots[idx].imag > 0:
            rep = w4_roots[idx]
        else:
            rep = w4_roots[best]
        pairs.append({"rep": rep, "pair": (w4_roots[idx], w4_roots[best])})

    # sort pairs by Re(rep) for stable identity
    pairs.sort(key=lambda d: float(np.real(d["rep"])))

    # S_{+-}(v) chord integrals
    # Reference lambda: pick a real lambda not too close to any eps, taken
    # as the midpoint between the two outer eps. (A small real shift avoids
    # the LH zero exactly.)
    lam_ref = 0.5 * (geo.eps[0] + geo.eps[-1])
    # nudge if too close to LH root or any eps
    while any(abs(lam_ref - e) < 0.05 for e in geo.eps) or abs(lam_ref - geo.LH_root()) < 0.05:
        lam_ref += 0.07

    S_pm = []
    u_W4 = []
    I_W4 = []
    for d in pairs:
        v = d["rep"]
        # avoid hitting eps poles by going around them: split the chord into
        # two arcs that detour above the real axis (Im part already in v)
        # The integrand has order-3 poles at eps_i (real). The chord from
        # lam_ref (real) to v (Im > 0) generically misses the real eps_i,
        # so a single chord is fine in the upper half plane.
        S = _omega_chord_integral(geo, lam_ref, v, n_steps=2000)
        S_pm.append(complex(S))
        # u-image
        p_val = np.polyval(geo.p, v)
        n_val = np.polyval(geo.n, v)
        u_W4.append(complex(n_val / p_val))
        # W_4 loop period: radius ~ 1.5 x half-separation of the W_4 pair,
        # but careful not to enclose any Q_4 root (would re-introduce branch)
        v_bar = np.conj(v)
        r_w4 = 1.5 * abs(v - v_bar) / 2.0
        # ensure radius does not reach a Q_4 root
        Q4_roots = geo.Q4_roots()
        cw4 = 0.5 * (v + v_bar)
        for q in Q4_roots:
            if abs(q - cw4) < r_w4 * 1.1:
                r_w4 = min(r_w4, 0.6 * abs(q - cw4))
        # also avoid the eps poles
        for e in geo.eps:
            if abs(e - cw4) < r_w4 * 1.2:
                r_w4 = min(r_w4, 0.5 * abs(e - cw4))
        Iw4 = _w4_loop_period(geo, cw4, r_w4, n=4000)
        I_W4.append(complex(Iw4))

    return {
        "I_X": IX,
        "Gamma_X": GammaX,
        "v_W4": [d["rep"] for d in pairs],
        "S_pm": S_pm,
        "|Re S_pm|": [abs(z.real) for z in S_pm],
        "|Im S_pm|": [abs(z.imag) for z in S_pm],
        "I_W4": I_W4,
        "|Im I_W4|": [abs(z.imag) for z in I_W4],
        "|Re I_W4|": [abs(z.real) for z in I_W4],
        "u_W4": u_W4,
    }


# ---------------------------------------------------------------------------
#  Endpoint permutations: the IP adiabatic basis
# ---------------------------------------------------------------------------
# The IP harness uses *exterior-first interval* labelling for the adiabatic
# eigenvalues lam_i, with the exterior root lam_0 carrying the asymptotically
# extreme eigenvalue at each end. With the diabatic energies E_i(u) ~ a_i u,
#   adia 0 at u -> -inf  = diabatic argmin(a)  (highest energy at u<<0)
#   adia 0 at u -> +inf  = diabatic argmax(a)
# matching the existing assay.validation_matrix.be_only_formula convention.
# So the per-record permutations are
#   pi_in[i] = diabatic level connected to adia i at u -> -inf  = argsort(a)[i]
#   pi_out[i] = diabatic level connected to adia i at u -> +inf  = argsort(-a)[i]
# This is the convention validated to ~1e-9 against P_bench by the BE
# diagonal entries (see validation_matrix.be_only_formula).


def _endpoint_perms(a):
    """(pi_in, pi_out) for the IP basis given slope vector a."""
    a = np.asarray(a, float)
    return np.argsort(a), np.argsort(-a)


def _be_extreme_positions(rec):
    """Adia (in, out) coords of the BE-extreme + middle entries.
    pi_in/pi_out depend on the slope ordering."""
    a = np.array(rec.a, float)
    pi_in, pi_out = _endpoint_perms(a)
    n_max = int(np.argmax(a))
    n_min = int(np.argmin(a))
    n_mid = (set(range(3)) - {n_max, n_min}).pop()
    return {
        "n_max": n_max, "n_mid": n_mid, "n_min": n_min,
        "i_in_max": int(np.where(pi_in == n_max)[0][0]),
        "i_out_max": int(np.where(pi_out == n_max)[0][0]),
        "i_in_mid": int(np.where(pi_in == n_mid)[0][0]),
        "i_out_mid": int(np.where(pi_out == n_mid)[0][0]),
        "i_in_min": int(np.where(pi_in == n_min)[0][0]),
        "i_out_min": int(np.where(pi_out == n_min)[0][0]),
    }


def _BE_diagonal(rec) -> Tuple[Tuple[int, int], float, Tuple[int, int], float]:
    """
    Brundobler-Elser diagonal entries in the IP adiabatic endpoint basis.
    """
    pos = _be_extreme_positions(rec)

    def q_prod(n):
        p = 1.0
        for m in range(3):
            if m == n:
                continue
            pair = (n, m) if n < m else (m, n)
            p *= rec.q_ij[pair]
        return p

    return ((pos["i_in_max"], pos["i_out_max"]), q_prod(pos["n_max"]),
            (pos["i_in_min"], pos["i_out_min"]), q_prod(pos["n_min"]))


# ---------------------------------------------------------------------------
#  Hypothesis formulas
# ---------------------------------------------------------------------------
def _geo_of_record(rec):
    par = Params(eps=rec.eps, gam=rec.gam, a=rec.a, x=rec.x)
    return Geometry(par)


def _be_skeleton_P(rec, mid_survival: float) -> np.ndarray:
    """
    Build a 3x3 BE-skeleton with the two BE-extreme entries set exactly and
    the middle-survival entry set to `mid_survival`. The remaining off-
    diagonals are filled by minimum-norm doubly-stochastic completion within
    the row/column residuals.

    Returns the 3x3 prediction in the IP adia basis (rows = incoming
    adia, cols = outgoing adia).
    """
    pos = _be_extreme_positions(rec)
    # BE q-products
    def q_prod(n):
        p = 1.0
        for m in range(3):
            if m != n:
                pair = (n, m) if n < m else (m, n)
                p *= rec.q_ij[pair]
        return p
    qmax = q_prod(pos["n_max"])
    qmin = q_prod(pos["n_min"])
    qmid = max(0.0, min(1.0, mid_survival))

    # Place the three diagonal-style survivals at their (i_in, i_out)
    P = np.zeros((3, 3))
    P[pos["i_in_max"], pos["i_out_max"]] = qmax
    P[pos["i_in_min"], pos["i_out_min"]] = qmin
    P[pos["i_in_mid"], pos["i_out_mid"]] = qmid
    # Row and column sums determine the remaining four off-diagonal entries
    # uniquely (each row/col has one deficit slot in 3x3 with three diagonal
    # entries pinned -- but in our IP basis, the three "survivors" sit at
    # (i_in_n, i_out_n) which is *not* a diagonal of the matrix in general).
    # We need to allocate row deficits to the off-diagonal slots, subject to
    # doubly-stochastic constraints.
    # A simple closed-form for the deficit pattern in the IP basis:
    # Each row i has one "main" slot at (i, INV_PI_OUT_of_diabatic_PI_IN_i).
    # The remaining two cells in row i are off-diagonals.
    # We fill by Sinkhorn from a uniform off-diagonal seed.
    row_def = 1.0 - P.sum(axis=1)
    col_def = 1.0 - P.sum(axis=0)
    # Off-diagonal mask = cells where P is currently 0
    mask = (P == 0.0)
    if mask.sum() == 0:
        return P
    # Initialise off-diagonal cells uniformly so each row's deficit splits
    # equally between its two zero cells; same for columns. Use a small
    # Sinkhorn fixup to hit DS exactly.
    init = np.zeros_like(P)
    for i in range(3):
        nz = mask[i].sum()
        if nz > 0:
            init[i, mask[i]] = max(row_def[i], 0.0) / nz
    P_full = P + init
    # Sinkhorn-Knopp normalisation just on the off-diagonal cells, preserving
    # the three survivor entries.
    survivors = ~mask
    for _ in range(50):
        # row pass: scale off-diagonal cells in row i to absorb row deficit
        for i in range(3):
            rs = P_full[i].sum()
            if rs == 0: continue
            scale = (1.0 - P_full[i][survivors[i]].sum()) / (
                P_full[i][mask[i]].sum() + 1e-300)
            P_full[i][mask[i]] *= max(scale, 0.0)
        # column pass: same for columns
        for j in range(3):
            cs = P_full[:, j].sum()
            if cs == 0: continue
            scale = (1.0 - P_full[:, j][survivors[:, j]].sum()) / (
                P_full[:, j][mask[:, j]].sum() + 1e-300)
            P_full[:, j][mask[:, j]] *= max(scale, 0.0)
        # convergence check
        if (abs(P_full.sum(axis=1) - 1).max() < 1e-12 and
                abs(P_full.sum(axis=0) - 1).max() < 1e-12):
            break
    return P_full


def _be_only_predicted(rec) -> np.ndarray:
    """BE-only skeleton with middle survival = q_prod(n_mid) (i.e. naive)."""
    def q_prod(n):
        p = 1.0
        for m in range(3):
            if m != n:
                pair = (n, m) if n < m else (m, n)
                p *= rec.q_ij[pair]
        return p
    pos = _be_extreme_positions(rec)
    return _be_skeleton_P(rec, q_prod(pos["n_mid"]))


def _middle_q_pair(rec) -> float:
    """The single q_ij of the pair NOT involving the steepest/shallowest
    level. For 3 levels with sorted indices, this is q_{i_low, i_high} where
    i_low = min(non_extreme indices), but actually all three q_ij pairs
    are relevant. We use q_{n_min, n_max} (the BE-extreme pair) as the
    one which the middle level does NOT directly cross.

    Actually: the three crossings involve the three q_ij. The middle level
    survival probability under the grid is q_{mid, min} * q_{mid, max}.
    """
    pos = _be_extreme_positions(rec)
    others = [pos["n_min"], pos["n_max"]]
    pair = tuple(sorted(others))
    return rec.q_ij[pair]


def _Gamma_middle_pair(rec) -> float:
    """Gamma_{n_min, n_max} = the BE-extreme pair coupling."""
    pos = _be_extreme_positions(rec)
    i_, j_ = sorted([pos["n_min"], pos["n_max"]])
    return (rec.gam[i_]**2 * rec.gam[j_]**2 * abs(rec.a[i_] - rec.a[j_])
            / (rec.eps[i_] - rec.eps[j_])**2)


def formula_H_A(rec) -> np.ndarray:
    """
    H_A: middle-level survival = q_{n_min,n_max} * exp(-|Im S_{+-}|).
    Uses the larger of the two transport-side |Im S| residuals (the more
    decisive one).
    """
    cli = cross_link_invariants(_geo_of_record(rec))
    im_S = max(cli["|Im S_pm|"])
    q_pair = _middle_q_pair(rec)
    mid_survival = q_pair * float(np.exp(-im_S))
    return _be_skeleton_P(rec, mid_survival)


def formula_H_B(rec) -> np.ndarray:
    """
    H_B: middle survival = exp(-|Re S_{+-}|) * q_{n_min,n_max}, using the
    smaller |Re S_{+-}|.
    """
    cli = cross_link_invariants(_geo_of_record(rec))
    re_S = min(cli["|Re S_pm|"])
    q_pair = _middle_q_pair(rec)
    mid_survival = q_pair * float(np.exp(-re_S))
    return _be_skeleton_P(rec, mid_survival)


def formula_H_C(rec) -> np.ndarray:
    """
    H_C: middle survival exponent equals |Im I_{W4}| of the smaller-modulus
    W_4 loop period.  (Numerical-zero baseline; will fail by construction
    since the omega-integrand has no W_4 residue, but recorded as falsification.)
    """
    cli = cross_link_invariants(_geo_of_record(rec))
    im_w4 = min(cli["|Im I_W4|"])
    return _be_skeleton_P(rec, float(np.exp(-im_w4)))


def formula_H_D(rec) -> np.ndarray:
    """
    H_D: mixed-cover linear combination, middle survival:
       exp( - 2 pi * Gamma_{n_min, n_max}  -  |Im S_{+-}_smaller| )
    """
    cli = cross_link_invariants(_geo_of_record(rec))
    G_pair = _Gamma_middle_pair(rec)
    im_S = min(cli["|Im S_pm|"])
    mid_survival = float(np.exp(-2.0 * np.pi * G_pair - im_S))
    return _be_skeleton_P(rec, mid_survival)


def formula_H_E_grid_with_transport_correction(rec) -> np.ndarray:
    """
    H_E: incoherent LZ grid reindexed to IP adia basis using the per-record
    permutations pi_in = argsort(a), pi_out = argsort(-a).
    """
    from .grid import grid_P
    geo = _geo_of_record(rec)
    P_dia = grid_P(geo)
    pi_in, pi_out = _endpoint_perms(rec.a)
    return P_dia[np.ix_(pi_in, pi_out)]


# ---------------------------------------------------------------------------
#  Naive q_mid + DS completion ("middle-pair survival") -- best non-transport
#  baseline; we measure how much improvement we get over this with transport
# ---------------------------------------------------------------------------
def formula_baseline_q_mid_pair(rec) -> np.ndarray:
    """Mid survival = q_{n_min, n_max} (the BE-extreme pair), no transport."""
    return _be_skeleton_P(rec, _middle_q_pair(rec))


def formula_baseline_be_only(rec) -> np.ndarray:
    """BE-only with naive mid survival = q_prod(n_mid)."""
    return _be_only_predicted(rec)


def formula_baseline_BE_extreme_zero_middle(rec) -> np.ndarray:
    """BE diagonals; middle survival fixed to 0; DS completion."""
    return _be_skeleton_P(rec, 0.0)


def formula_baseline_BE_extreme_one_middle(rec) -> np.ndarray:
    """BE diagonals; middle survival fixed to 1; DS completion."""
    return _be_skeleton_P(rec, 1.0)


# ---------------------------------------------------------------------------
#  Doubly-stochastic projection / minimum-norm completion
# ---------------------------------------------------------------------------
def _project_doubly_stochastic(P: np.ndarray, n_iter: int = 200,
                                tol: float = 1e-12) -> np.ndarray:
    """
    Sinkhorn-Knopp projection onto doubly stochastic 3x3 matrices,
    starting from |P|. Used for hypothesis formula completions.
    """
    M = np.maximum(np.abs(P), 1e-300)
    for _ in range(n_iter):
        row_sums = M.sum(axis=1)
        M /= row_sums[:, None]
        col_sums = M.sum(axis=0)
        M /= col_sums[None, :]
        # check
        if (abs(row_sums - 1).max() < tol
                and abs(col_sums - 1).max() < tol):
            break
    return M


# ---------------------------------------------------------------------------
#  Best-of: enumerate permutations & sign assignments to compare against
#  P_bench. Useful as a *diagnostic* (to see if there is ANY assignment
#  reaching machine zero residual under a single transport+phase formula).
# ---------------------------------------------------------------------------
def best_match_among(formulas: List, rec) -> Tuple[str, float, np.ndarray]:
    bests = []
    for name, fml in formulas:
        try:
            P = np.asarray(fml(rec), float)
            err = float(np.max(np.abs(P - rec.P_bench)))
            bests.append((err, name, P))
        except Exception as e:
            bests.append((float("inf"), name + f" (err: {e!r})", None))
    bests.sort()
    e, n, p = bests[0]
    return n, e, p


# ---------------------------------------------------------------------------
#  Selector residual: probe specific theta-active locus
# ---------------------------------------------------------------------------
def selector_residual_theta_scan(geo: Geometry,
                                 theta_grid: np.ndarray) -> Dict[str, np.ndarray]:
    """
    For each theta, compute Im(e^{-i theta} S_{+-}(v_alpha)) for both
    W_4 roots. The selector becomes active where this vanishes
    (codim 1 in (theta, params)).
    """
    cli = cross_link_invariants(geo)
    S = cli["S_pm"]
    out = {}
    for k, S_k in enumerate(S):
        rsd = np.imag(np.exp(-1j * theta_grid) * S_k)
        out[f"v_{k}"] = rsd
    return out


__all__ = [
    "cross_link_invariants",
    "formula_H_A", "formula_H_B", "formula_H_C", "formula_H_D",
    "formula_H_E_grid_with_transport_correction",
    "formula_baseline_q_mid_pair",
    "formula_baseline_be_only",
    "formula_baseline_BE_extreme_zero_middle",
    "formula_baseline_BE_extreme_one_middle",
    "selector_residual_theta_scan",
    "best_match_among",
]
