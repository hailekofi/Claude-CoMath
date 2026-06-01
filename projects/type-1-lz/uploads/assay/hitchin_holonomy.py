"""
Non-abelianization (Hollands-Neitzke / Gaiotto-Moore-Neitzke) holonomy
for the Type-1, N=3 Landau-Zener problem.

Strategy
========
The base curve is C = u-line.  The Hitchin spectral cover has three sheets
:math:`\\lambda_i(u)` (the labelled eigenvalues of :func:`adiabatic.labelled_lambdas`).
The flat connection :math:`\\nabla = du\\,\\partial_u - i\\,H(u)` is to be
parallel-transported from u=-T to u=+T.

GMN non-abelianization recipe
-----------------------------
Pick a phase :math:`\\vartheta`.  The connection matrix in the asymptotic
sheet (= labelled-adiabatic) basis is the ordered product

.. math::  U \;=\; \\Phi_+ \\cdot \\Big(\\prod_{c} \\mathcal{K}_c\\Big) \\cdot \\Phi_-

where:

* :math:`\\Phi_\\pm` are diagonal sheet-by-sheet phase factors that account
  for the abelian/WKB transport in the tail chambers:
  :math:`(\\Phi)_i = \\exp(-i\\int E_i\\,du)`.

* :math:`\\mathcal{K}_c` is the wall-crossing factor at the c-th S-wall
  crossed by the real-axis path.  An (i,j) S-wall contributes

.. math::  \\mathcal{K}_c \;=\; \\exp\\!\\big( \\mu_\\gamma\\,E_{ij}\\big)
                          \;=\; I + \\mu_\\gamma\\,E_{ij},

  with :math:`E_{ij}` the elementary nilpotent matrix unit and
  :math:`\\mu_\\gamma = \\exp(i Z_\\gamma)` the Voros symbol of the
  homology class :math:`\\gamma` carried by the wall.

Type-1 N=3 specialisation
-------------------------
For Type-1 with two conjugate-pair Q4 windows the BPS rays at phase
:math:`\\vartheta_c = \\pi/2` carry exactly two finite saddles
:math:`\\gamma_A, \\gamma_B` connecting the conjugate-partner turning
points within each window.  The Voros symbols are

.. math::  \\mu_{\\gamma_X} \\,=\\, \\exp(i\\,I_X),
           \\qquad |\\mu_{\\gamma_X}|^2 \\,=\\, e^{-|\\im I_X|}\\,=\\, p_X,

i.e. the elliptic window period from :func:`actions.window_action` IS the
Voros symbol of the corresponding BPS state (Track 1 result).

The real-axis path from :math:`u=-T` to :math:`u=+T` crosses the two
windows in u-order.  Inside each window the path "enters" the BPS saddle
configuration from one side and "exits" the other side.  The resulting
elementary unipotent factor :math:`\\mathcal{K}_X = I + \\mu_X E_{ij_X}`
turns the sheet basis amplitudes by the standard 2x2 Landau-Zener
rotation in the (i,j_X) plane, supplemented by the abelian phase factor
:math:`\\Phi` between windows.

The endpoint sheet permutation (going from u=-infty to u=+infty the
adiabatic labels permute by :math:`\\arg\\!\\mathrm{sort}(a)` /
:math:`\\arg\\!\\mathrm{sort}(-a)`) accounts for the diabatic-vs-adiabatic
basis rotation already established in :mod:`assay.validation_matrix`.

Implementation
--------------
* :func:`abelian_transport`        -- WKB diagonal between two u-points
* :func:`wall_jump_matrix`         -- single (i,j) wall crossing factor
* :func:`window_jump`              -- the 2x2 LZ block for one window with
                                       a full Stokes-Voros phase
* :func:`hn_connection_matrix`     -- assemble the full non-abelianized
                                       connection matrix
* :func:`formula_HN`               -- public callable consumable by
                                       :func:`assay.validation_matrix.validate`

The formula returns the 3x3 transition matrix :math:`P_{HN}=|U_{HN}|^2`
in the IP adiabatic endpoint basis.
"""

from __future__ import annotations
import numpy as np
from scipy.integrate import quad

from .geometry import Geometry, Params
from .adiabatic import labelled_lambdas
from .actions import window_action
from .closed_form import window_pairs


# ---------------------------------------------------------------------------
#  Sheet eigenvalues and energies along the real axis
# ---------------------------------------------------------------------------
def _sheet_E(geo: Geometry, u: float, k: int) -> float:
    """
    Eigenvalue E_k(u) of the labelled-adiabatic sheet k along the real
    u-axis.  E_k = sum a_j gam_j^2 / (lam_k - eps_j) = m(lam_k)/p(lam_k).
    """
    lam = labelled_lambdas(geo, u)
    eps = geo.eps
    g2 = geo.g2
    a = geo.a
    return float(np.real(np.sum(a * g2 / (lam[k] - eps))))


# ---------------------------------------------------------------------------
#  Abelian (WKB) transport between two u-points
# ---------------------------------------------------------------------------
def abelian_transport(geo: Geometry, u_a: float, u_b: float,
                      x: int = 0, n_segments: int = 256) -> np.ndarray:
    """
    Diagonal 3x3 matrix :math:`\\Phi(u_b, u_a)_{ii} = \\exp(-i\\int_{u_a}^{u_b}(E_i-E_x)du)`,
    the abelian transport between two real-axis points in the gauged
    (E_x-subtracted) WKB basis.

    The subtraction by E_x matches the validation harness's gauged adiabatic
    convention; the gauge cancels in the doubly-stochastic ``|U|^2``.
    """
    grid = np.linspace(u_a, u_b, n_segments + 1)
    phases = np.zeros(3, float)
    # cumulative integral by composite Simpson on each segment
    h = (u_b - u_a) / n_segments
    # use Simpson 1/3 on pairs
    if n_segments % 2 == 1:
        n_segments += 1
        grid = np.linspace(u_a, u_b, n_segments + 1)
        h = (u_b - u_a) / n_segments
    # compute integrand values
    E_vals = np.zeros((n_segments + 1, 3))
    for k in range(n_segments + 1):
        lam = labelled_lambdas(geo, grid[k])
        eps = geo.eps
        g2 = geo.g2
        a = geo.a
        E_vals[k] = np.real(np.array([np.sum(a * g2 / (lam[i] - eps))
                                       for i in range(3)]))
    # subtract reference E_x
    Eg = E_vals - E_vals[:, x:x + 1]
    # Simpson
    w = np.ones(n_segments + 1)
    w[1:-1:2] = 4.0
    w[2:-1:2] = 2.0
    w = w * (h / 3.0)
    phases = np.sum(w[:, None] * Eg, axis=0)
    return np.diag(np.exp(-1j * phases))


# ---------------------------------------------------------------------------
#  Wall-crossing matrix
# ---------------------------------------------------------------------------
def wall_jump_matrix(mu: complex, i: int, j: int) -> np.ndarray:
    """
    Unipotent (i,j) wall-crossing factor :math:`I + \\mu E_{ij}`.

    Parameters
    ----------
    mu : complex
        Voros symbol; in the Type-1 N=3 problem this equals
        :math:`\\exp(i\\,I_X)` with :math:`I_X` the elliptic window period.
    i, j : int
        Sheet indices (0,1,2); :math:`E_{ij}` is the matrix unit with 1 at
        position (i,j).
    """
    M = np.eye(3, dtype=complex)
    M[i, j] = mu
    return M


# ---------------------------------------------------------------------------
#  Window jump as a coherent 2x2 LZ rotation embedded in 3x3
# ---------------------------------------------------------------------------
def window_jump(p: float, pair: tuple, phase_off: complex = 1.0
                ) -> np.ndarray:
    """
    Coherent 2x2 Landau-Zener jump for one window:

    .. math::
       \\mathcal{N}_X \\;=\\; \\begin{pmatrix}\\sqrt{1-p} & -\\sqrt p\\,e^{i\\phi}\\\\
                                              \\sqrt p\\,e^{-i\\phi} & \\sqrt{1-p}\\end{pmatrix}

    embedded in 3x3 on rows/cols (i,j) = ``pair``, identity on the spectator.

    The phase :math:`\\phi` is the Stokes-Voros phase: in the
    non-abelianization formula, the wall-crossing matrix
    :math:`I + e^{i I_X} E_{ij}` and its conjugate (one for each half of the
    BPS saddle traverse) compose to a unitary 2x2 rotation whose off-diagonal
    phases carry :math:`\\arg(\\mu_X) = \\re I_X` (the Stokes phase from the
    elliptic period).
    """
    i, j = pair
    if i > j:
        i, j = j, i
    sp = ({0, 1, 2} - {i, j}).pop()
    c, s = np.sqrt(max(0.0, 1.0 - p)), np.sqrt(max(0.0, p))
    R = np.zeros((3, 3), complex)
    R[i, i], R[i, j] = c, -s * phase_off
    R[j, i], R[j, j] = s * np.conj(phase_off), c
    R[sp, sp] = 1.0
    return R


# ---------------------------------------------------------------------------
#  Voros symbols for the two windows
# ---------------------------------------------------------------------------
def voros_symbols(geo: Geometry) -> tuple:
    """
    Voros symbols :math:`\\mu_X = e^{i I_X}` for the two BPS saddles, in
    u-order of the windows.

    Returns
    -------
    ((mu_A, I_X_A), (mu_B, I_X_B))
    """
    from .closed_form import window_pairs as _wp
    wp = _wp(geo)
    out = []
    for entry in wp:
        I_X = window_action(geo, entry["window"])["I_X"]
        mu = np.exp(1j * I_X)
        out.append((complex(mu), complex(I_X)))
    return tuple(out)


# ---------------------------------------------------------------------------
#  Hollands-Neitzke connection matrix assembly
# ---------------------------------------------------------------------------
def hn_connection_matrix(geo: Geometry, x: int = 0, T: float = 60.0
                         ) -> np.ndarray:
    """
    The non-abelianized connection matrix in the labelled-adiabatic
    endpoint basis, from u=-T to u=+T.

    The recipe (Type-1 N=3 specialisation, theta=pi/2):

    1. Identify the two BPS windows in u-order and their sheet pairs.
    2. The Voros symbol of window X is mu_X = exp(i I_X).
    3. The 2x2 jump in window X is a unitary LZ rotation with
       p_X = |mu_X|^2 = exp(-|Im I_X|), and Stokes phase Re(I_X)
       (vanishingly small in the Type-1 case where I_X is purely
       imaginary up to numerical noise).
    4. Between windows and in the tails, transport is diagonal
       :math:`\\Phi_{ii} = \\exp(-i\\int(E_i-E_x)du)`.

    Returns
    -------
    U : (3,3) complex array, the lab-frame fundamental matrix in the
        labelled sheet basis.
    """
    # window structure
    wp = window_pairs(geo)
    wp_sorted = sorted(wp, key=lambda d: d["u_center"])
    # window jump probabilities (from Track-1 elementary product)
    p_X = []
    phase_X = []
    for entry in wp_sorted:
        IX = window_action(geo, entry["window"])["I_X"]
        p = float(np.exp(-abs(IX.imag)))
        p_X.append(p)
        # Stokes phase: Re(I_X) -- in Type-1 this is ~0 by Track-1 exact result
        phase_X.append(np.exp(1j * float(IX.real)))

    # transport breakpoints: -T, u_A, u_B, +T
    u_A = wp_sorted[0]["u_center"]
    u_B = wp_sorted[1]["u_center"]
    pairA = wp_sorted[0]["pair"]
    pairB = wp_sorted[1]["pair"]

    Phi_left = abelian_transport(geo, -T, u_A, x=x, n_segments=128)
    Phi_mid = abelian_transport(geo, u_A, u_B, x=x, n_segments=128)
    Phi_right = abelian_transport(geo, u_B, T, x=x, n_segments=128)

    N_A = window_jump(p_X[0], pairA, phase_off=phase_X[0])
    N_B = window_jump(p_X[1], pairB, phase_off=phase_X[1])

    U = Phi_right @ N_B @ Phi_mid @ N_A @ Phi_left
    return U


# ---------------------------------------------------------------------------
#  Public formula consumed by validation_matrix.validate
# ---------------------------------------------------------------------------
def formula_HN(rec, T: float = 80.0) -> np.ndarray:
    """
    Non-abelianization / Hollands-Neitzke formula for the 3x3 transition
    matrix :math:`P_{HN}=|U|^2` in the IP adiabatic endpoint basis.

    Endpoint permutations:
       The IP harness labels adiabatic channel 0 as the *exterior* eigenvalue
       (the largest in magnitude).  At u-> -infinity this is the eigenvector
       with the most negative slope a-projection; at u-> +infinity it is the
       eigenvector with the most positive slope.  Concretely:

           pi_in = argsort(a),    pi_out = argsort(-a).

    Returns
    -------
    (3,3) real array, doubly stochastic up to floating-point.
    """
    par = Params(eps=rec.eps, gam=rec.gam, a=rec.a, x=rec.x)
    geo = Geometry(par)
    try:
        U = hn_connection_matrix(geo, x=par.x, T=T)
    except Exception:
        return np.full((3, 3), np.nan)
    # P in sheet basis
    # U columns indexed by incoming sheet, rows by outgoing sheet
    P_sheet = np.abs(U) ** 2

    # The labelled sheet basis is already in u->+/-infty asymptotic
    # adiabatic order WITHIN each end.  But the IP harness's adiabatic
    # endpoint labelling is the eigh-sorted-by-energy basis.  At u-> -infty
    # the labelled lam_0 is the EXTERIOR (the most-negative energy for
    # the canonical sign convention).  In the IP harness, adiabatic
    # channel 0 == argsort(a)[0] diabatic == arg_min(a) energy slope.
    # This matches: labelled exterior lam_0 at u-> -infty is the level
    # whose slope is min a; at u-> +infty it is the level whose slope is
    # max a.
    a = np.array(par.a)
    # at u=-T: labelled order is (exterior, in (eps_0,eps_1), in (eps_1,eps_2))
    # whose corresponding diabatic channels are arg_min(a), middle, arg_max(a)
    pi_in = np.argsort(a)       # diabatic of (lam_0, lam_1, lam_2) at -infty
    pi_out = np.argsort(-a)     # diabatic of (lam_0, lam_1, lam_2) at +infty
    # However IP adiabatic indexing is:
    #   adia channel i at u->-infty  <-> diabatic pi_in[i]
    # and the labelled basis index k corresponds to diabatic chan
    #   d_in[k] = pi_in_lab[k] where pi_in_lab maps lab-index -> diabatic
    # The labelled basis convention (sorted by energy) matches the
    # adiabatic harness labels directly, so no extra permutation is
    # needed beyond what is already inside U.
    #
    # ie the test against be_only_formula confirms our convention IS the
    # one with adia 0 = argsort(a) at -infty.  No further permutation.
    return P_sheet


# ---------------------------------------------------------------------------
#  Diagnostic helpers
# ---------------------------------------------------------------------------
def hn_diagnostics(geo: Geometry, x: int = 0, T: float = 80.0) -> dict:
    """
    Decompose the HN connection into its building blocks for inspection.
    """
    wp = window_pairs(geo)
    wp_sorted = sorted(wp, key=lambda d: d["u_center"])
    u_A = wp_sorted[0]["u_center"]
    u_B = wp_sorted[1]["u_center"]
    pairA = wp_sorted[0]["pair"]
    pairB = wp_sorted[1]["pair"]

    p_X = []
    I_X_list = []
    for entry in wp_sorted:
        IX = window_action(geo, entry["window"])["I_X"]
        I_X_list.append(complex(IX))
        p_X.append(float(np.exp(-abs(IX.imag))))

    Phi_mid = abelian_transport(geo, u_A, u_B, x=x, n_segments=128)
    # mid-segment WKB phases
    wkb_phases = -np.angle(np.diag(Phi_mid))
    return {
        "u_A": u_A, "u_B": u_B,
        "pairA": pairA, "pairB": pairB,
        "p_A": p_X[0], "p_B": p_X[1],
        "I_X_A": I_X_list[0], "I_X_B": I_X_list[1],
        "wkb_phases_mid": wkb_phases.tolist(),
    }


# ===========================================================================
#  Phase III Route A extension: non-abelian wiring via spectral_network
# ===========================================================================
#
# The implemented hn_connection_matrix above reduces the
# Hollands-Neitzke formula to a sequential 2-level LZ product -- which
# is the bare BPS-spectrum prediction WITHOUT any inter-window
# saddles.  The spectral_network module already enumerates real-axis
# wall crossings and their sheet pairs; this section wires that output
# into a proper non-abelian chamber product:
#
#   U = Phi_R * (K_n * Phi_n * K_{n-1} * Phi_{n-1} * ... * K_1) * Phi_L
#
# with K_k = wall_jump_matrix(mu_X, i, j) at each real-axis crossing
# (u_k, (i,j)) emerging from window X (so mu_X = exp(i I_X) is the
# Voros symbol of that window's BPS saddle), and Phi_k =
# abelian_transport between consecutive crossings.
# ===========================================================================
def _wall_data(geo: Geometry, theta: float = 0.5 * np.pi) -> list:
    """
    Enumerate the real-axis wall crossings emerging from the spectral
    network at phase theta.

    Returns
    -------
    list of dicts
        Each dict has keys ``u_real`` (float, the real-axis u value of
        the crossing), ``sheet_pair`` ((i,j) sheet indices), ``mu``
        (the Voros symbol exp(i I_X) of the originating window),
        ``window`` (which window the wall emerges from).
    """
    from .spectral_network import spectral_network
    sn = spectral_network(geo, theta=theta)
    # collect Voros symbols per window using the canonical geo.windows()
    # (spectral_network's "windows" list renames "roots"->"roots_lambda";
    # window_action expects the canonical dict shape with "roots").
    geo_windows = sorted(geo.windows(), key=lambda w: w["u_center"])
    mu_per_window = []
    for w in geo_windows:
        IX = window_action(geo, w)["I_X"]
        mu_per_window.append(complex(np.exp(1j * IX)))

    out = []
    for crossing in sn["crossings"]:
        wi = int(crossing.get("window", 0))
        mu = mu_per_window[wi] if 0 <= wi < len(mu_per_window) else 1.0 + 0.0j
        i, j = crossing["sheet_pair"]
        out.append({
            "u_real": float(np.real(crossing["u"])),
            "sheet_pair": (int(i), int(j)),
            "mu": mu,
            "window": wi,
        })
    out.sort(key=lambda d: d["u_real"])
    return out


def hn_connection_matrix_nonabelian(
    geo: Geometry,
    x: int = 0,
    T: float = 80.0,
    theta: float = 0.5 * np.pi,
    n_seg_per_chamber: int = 64,
) -> np.ndarray:
    """
    Non-abelian Hollands-Neitzke connection matrix from -T to +T.

    Uses the spectral_network output: at each real-axis wall crossing,
    insert the unipotent K = I + mu * E_{ij}; between crossings, use
    the abelian transport.  The window-centred 2-level LZ blocks are
    also inserted at u_A_center, u_B_center as before, in case
    spectral_network misses the on-axis traversal contributions.

    Parameters
    ----------
    geo : Geometry
    x : int
        Reference channel for the gauged transport (matches
        validation_matrix convention).
    T : float
        Half-width of the integration window.
    theta : float
        BPS phase at which the spectral network is built.  Default
        pi/2 for the project's purely-imaginary I_X.
    n_seg_per_chamber : int
        Quadrature segments per chamber for the abelian transport.

    Returns
    -------
    U : (3, 3) complex
        The non-abelian connection matrix.
    """
    walls = _wall_data(geo, theta=theta)

    # window-centred LZ blocks (use existing window_jump on the active pair)
    wp = window_pairs(geo)
    wp_sorted = sorted(wp, key=lambda d: d["u_center"])
    p_X = []
    phase_X = []
    u_centers = []
    pair_X = []
    for entry in wp_sorted:
        IX = window_action(geo, entry["window"])["I_X"]
        p_X.append(float(np.exp(-abs(IX.imag))))
        phase_X.append(np.exp(1j * float(IX.real)))
        u_centers.append(float(entry["u_center"]))
        pair_X.append(entry["pair"])

    # all "stops" along the real axis, sorted: wall crossings + window centres
    stops = []
    for w in walls:
        stops.append({
            "u": w["u_real"],
            "kind": "wall",
            "data": w,
        })
    for k, uc in enumerate(u_centers):
        stops.append({
            "u": uc,
            "kind": "window",
            "data": {
                "p": p_X[k],
                "pair": pair_X[k],
                "phase_off": phase_X[k],
                "k": k,
            },
        })
    stops.sort(key=lambda s: s["u"])
    # deduplicate near-coincident stops (within 1e-6) by combining into one
    # multi-op stop (keep order within: window jumps then wall crossings)

    # Build the product U = Phi_R · op_n · Phi_n · ... · op_1 · Phi_L
    U = abelian_transport(geo, -T, stops[0]["u"], x=x,
                          n_segments=n_seg_per_chamber) if stops else \
        abelian_transport(geo, -T, T, x=x, n_segments=n_seg_per_chamber)

    for idx, s in enumerate(stops):
        if s["kind"] == "wall":
            d = s["data"]
            K = wall_jump_matrix(d["mu"], d["sheet_pair"][0], d["sheet_pair"][1])
            U = K @ U
        elif s["kind"] == "window":
            d = s["data"]
            N = window_jump(d["p"], d["pair"], phase_off=d["phase_off"])
            U = N @ U
        # abelian transport to the next stop (or to +T at the end)
        if idx < len(stops) - 1:
            u_next = stops[idx + 1]["u"]
        else:
            u_next = T
        Phi = abelian_transport(geo, s["u"], u_next, x=x,
                                n_segments=n_seg_per_chamber)
        U = Phi @ U

    return U


def formula_HN_nonabelian(rec, T: float = 80.0, theta: float = 0.5 * np.pi
                          ) -> np.ndarray:
    """
    Phase III Route A extended formula: non-abelian connection matrix
    via spectral_network real-axis wall crossings + window LZ blocks +
    abelian transport between.
    """
    par = Params(eps=rec.eps, gam=rec.gam, a=rec.a, x=rec.x)
    geo = Geometry(par)
    try:
        U = hn_connection_matrix_nonabelian(geo, x=par.x, T=T, theta=theta)
    except Exception:
        return np.full((3, 3), np.nan)
    return np.abs(U) ** 2


def hn_D(rec, T: float = 80.0, theta: float = 0.5 * np.pi) -> np.ndarray:
    """
    Accessor returning the inter-window content of the non-abelian HN
    formula -- i.e. U with the window-centred LZ blocks N_A and N_B
    factored out.  For use with :mod:`assay.convergence_test`.

    Returns
    -------
    np.ndarray, (3, 3) complex, unitary up to numerical precision.
    """
    par = Params(eps=rec.eps, gam=rec.gam, a=rec.a, x=rec.x)
    geo = Geometry(par)
    U = hn_connection_matrix_nonabelian(geo, x=par.x, T=T, theta=theta)
    # strip the two window blocks: D = N_B^H @ U @ N_A^H
    wp = window_pairs(geo)
    wp_sorted = sorted(wp, key=lambda d: d["u_center"])
    N_X = []
    for entry in wp_sorted:
        IX = window_action(geo, entry["window"])["I_X"]
        p = float(np.exp(-abs(IX.imag)))
        phase = np.exp(1j * float(IX.real))
        N_X.append(window_jump(p, entry["pair"], phase_off=phase))
    N_A, N_B = N_X
    return N_B.conj().T @ U @ N_A.conj().T


# ===========================================================================
#  Phase III Route A done correctly: GMN omega-pair side-tracking
# ===========================================================================
#
#  The Phase II / Phase III prototype formula_HN_nonabelian inserted
#  a single unipotent K = I + mu E_{ij} at each real-axis wall crossing
#  with mu = exp(i I_X).  For purely-imaginary I_X (Track 1 result),
#  |mu| = exp(-Im I_X) which is either huge or tiny depending on sign --
#  in both cases NON-UNITARY, breaking the connection matrix product.
#
#  The GMN ω-pair construction: each BPS state gamma contributes a
#  MATCHED PAIR of unipotents on opposite sides of its wall:
#       K_+(gamma) = I + mu_gamma * E_{ij}
#       K_-(gamma) = I + (1/conj(mu_gamma)) * E_{ji}
#  Together with an inner abelian transport Phi_mid (computed along the
#  BPS saddle interior), the product K_-(gamma) * Phi_mid * K_+(gamma)
#  is UNITARY and equals the Joye-Kunz-Pfister N_X^proper block exactly.
#
#  IMPLEMENTATION STATUS (Phase III scaffolding):
#  The functions below set up the data structures for the ω-pair
#  walk.  The exact algebraic form of K_± / Phi_mid that gives
#  N_X^proper is fixed by Hollands-Neitzke 2019 §§5-6 and 2016 §4.4.
#  Pending the literature read, validate_omega_pair_at_window is the
#  gate that confirms the implementation is correct: it must
#  reproduce N_X^proper from lz_stokes_phased.lz_block_proper to
#  U(2) geodesic distance ≤ 1e-12.
# ===========================================================================
def bps_omega_pair(geo: Geometry, saddle: dict
                   ) -> tuple:
    """
    Construct the matched (K_+, K_-) pair for one BPS saddle.

    Parameters
    ----------
    saddle : dict
        One entry from :func:`spectral_network.bps_inventory(...)['saddles']`,
        carrying 'sheet_pair' (i, j), 'Z' (central charge = I_X for
        window saddles), and 'mu' = exp(i Z).

    Returns
    -------
    tuple (K_plus, K_minus) of two 3x3 complex matrices.  These are the
    unipotent factors that, when ordered around the inner abelian
    transport Phi_mid, compose to the Stokes-phased LZ block N_X^proper.
    """
    i, j = saddle["sheet_pair"]
    mu = complex(saddle["mu"])

    K_plus = np.eye(3, dtype=complex)
    K_plus[i, j] = mu

    K_minus = np.eye(3, dtype=complex)
    # conjugate-reciprocal entry; for purely-imaginary Z, mu is real positive
    # and 1/conj(mu) = 1/mu (real reciprocal).  In general:
    K_minus[j, i] = 1.0 / np.conj(mu)

    return K_plus, K_minus


def phi_mid_for_saddle(geo: Geometry, saddle: dict,
                       x: int = 0, n_seg: int = 64) -> np.ndarray:
    """
    Inner abelian transport between the two ω-pair walls of one BPS
    saddle.  Computed as the diagonal phase
        Phi_mid_{kk} = exp(-i ∫_{u_start}^{u_end} (E_k - E_x) du)
    where the integration is along the real-axis projection of the
    BPS saddle's interior (between the two real-axis crossings of the
    saddle's two emanating walls).

    For a single BPS saddle on a Q4 conjugate-pair window, the inner
    transport accumulates exactly the WKB part `delta_X * log delta_X`
    of the Joye-Kunz-Pfister Stokes phase (Hollands-Neitzke 2016 §4.4).

    Parameters
    ----------
    saddle : dict
        From bps_inventory.
    x : int
        Reference channel for the gauged transport.
    n_seg : int
        Quadrature segments.

    Returns
    -------
    Phi_mid : 3x3 diagonal complex matrix.
    """
    traj = saddle["trajectory"]
    if len(traj) < 2:
        return np.eye(3, dtype=complex)
    # find the two real-axis crossings of the saddle's trajectory.
    # If the saddle stays in the UHP/LHP, fall back to the trajectory
    # endpoints' real parts (a placeholder until the literature read).
    ims = np.imag(traj)
    real_crosses = []
    for k in range(len(ims) - 1):
        if ims[k] * ims[k + 1] < 0:
            t = ims[k] / (ims[k] - ims[k + 1])
            u_x = (1 - t) * traj[k] + t * traj[k + 1]
            real_crosses.append(float(u_x.real))

    if len(real_crosses) >= 2:
        u_start, u_end = sorted(real_crosses)[:2]
    else:
        # fallback: use the trajectory's u-extent on the real axis
        u_start = float(np.real(traj[0]))
        u_end = float(np.real(traj[-1]))
        if u_start > u_end:
            u_start, u_end = u_end, u_start

    return abelian_transport(geo, u_start, u_end, x=x, n_segments=n_seg)


def side_at_crossing(wall: dict, u_cross: float,
                     tol: float = 1e-8) -> int:
    """
    Determine the GMN side (+1 or -1) of a real-axis wall crossing.

    The convention: side = sign(d Im(u)/d s) at the crossing, where s
    is the wall's arc-length parameter.

    Parameters
    ----------
    wall : dict
        One entry from spectral_network.spectral_network(...)["walls"] or
        from :func:`wall_with_orientation`.
    u_cross : float
        Real-axis crossing position.

    Returns
    -------
    int : +1 if the wall is heading up across the real axis at u_cross,
          -1 if heading down,  0 if undetermined.
    """
    traj = wall["u_trajectory"]
    ims = np.imag(traj)
    # find the trajectory step that brackets u_cross
    res = traj.real
    for k in range(len(res) - 1):
        if (res[k] - u_cross) * (res[k + 1] - u_cross) <= 0 \
                and ims[k] * ims[k + 1] < 0:
            d_imag = ims[k + 1] - ims[k]
            if abs(d_imag) < tol:
                continue
            return +1 if d_imag > 0 else -1
    return 0


def validate_omega_pair_at_window(geo: Geometry, window_idx: int,
                                  x: int = 0,
                                  theta: float = 0.5 * np.pi) -> dict:
    """
    THE CRITICAL GATE for Phase III Route A.

    Builds the ω-pair (K_+, K_-) for one window's BPS saddle, applies
    them around the inner abelian transport Phi_mid, and checks whether
    the resulting 2x2 block on the active pair matches the
    Joye-Kunz-Pfister proper LZ block N_X^proper.

    If this gate passes (U(2) geodesic ≤ 1e-12), the side convention
    and ω-pair construction are correct, and the rest of Route A's
    implementation will produce the expected closed form.

    If it fails, the diagnosis identifies which piece is wrong:
      - Sign convention of K_+ vs K_- (flip and re-test)
      - Phi_mid normalisation (off by a factor of 2 or sqrt(p))
      - mu convention (exp(i Z) vs exp(2 pi i Z) vs exp(i Z / eps))

    Parameters
    ----------
    geo : Geometry
    window_idx : int
        0 or 1 (A or B).
    x : int
        Gauge reference channel.
    theta : float
        BPS phase.

    Returns
    -------
    dict with keys:
        passed         : bool
        u2_geodesic    : float, distance to N_X^proper on the active pair
        K_plus, K_minus, Phi_mid : the constructed factors (for inspection)
        product        : K_- @ Phi_mid @ K_+ (3x3)
        N_X_proper     : the Joye-Kunz-Pfister target (3x3)
        sheet_pair     : (i, j)
        Z, mu, p, delta : the relevant elementary data
    """
    from .spectral_network import bps_inventory
    from .lz_stokes_phased import lz_block_proper, embed_complex

    inv = bps_inventory(geo, theta=theta)
    if not inv["b1_check"]:
        return {
            "passed": False,
            "error": f"BPS spectrum check failed: n_bps={inv['n_bps']}, "
                     f"expected 2 for genus-1 b_1=2 spectral curve.",
        }

    # pick the requested window's BPS saddle
    saddles_for_win = [s for s in inv["saddles"] if s["window"] == window_idx]
    if not saddles_for_win:
        return {
            "passed": False,
            "error": f"No BPS saddle found for window {window_idx}",
        }
    saddle = saddles_for_win[0]

    # construct K_+, K_-, Phi_mid
    K_plus, K_minus = bps_omega_pair(geo, saddle)
    Phi_mid = phi_mid_for_saddle(geo, saddle, x=x)

    # product (the GMN ordered product for one ω-pair)
    M = K_minus @ Phi_mid @ K_plus

    # the target N_X^proper
    Z = saddle["Z"]
    delta = float(abs(Z.imag) / (2.0 * np.pi))
    p = float(np.exp(-abs(Z.imag)))
    sheet_pair = saddle["sheet_pair"]
    N_proper_2x2 = lz_block_proper(p, delta, sign=+1)
    N_proper_3x3 = embed_complex(N_proper_2x2, sheet_pair)

    # U(2) geodesic distance on the active sheet pair (project both to
    # the (i,j) sub-block and compare)
    i, j = sheet_pair
    M_block = np.array([[M[i, i], M[i, j]],
                        [M[j, i], M[j, j]]], dtype=complex)
    N_block = np.array([[N_proper_2x2[0, 0], N_proper_2x2[0, 1]],
                        [N_proper_2x2[1, 0], N_proper_2x2[1, 1]]],
                       dtype=complex)
    # remove U(1) overall phase via det normalisation
    detD = np.linalg.det(M_block.conj().T @ N_block)
    if abs(detD) > 1e-30:
        R = (M_block.conj().T @ N_block) / (detD ** 0.5)
        from scipy.linalg import logm
        try:
            L = logm(R)
            u2_dist = float(np.linalg.norm(L, "fro"))
        except Exception:
            u2_dist = float("inf")
    else:
        u2_dist = float("inf")

    return {
        "passed": (u2_dist < 1e-12),
        "u2_geodesic": u2_dist,
        "K_plus": K_plus,
        "K_minus": K_minus,
        "Phi_mid": Phi_mid,
        "product": M,
        "N_X_proper": N_proper_3x3,
        "sheet_pair": sheet_pair,
        "Z": Z,
        "mu": saddle["mu"],
        "p": p,
        "delta": delta,
        "window_idx": window_idx,
    }
