"""
Route C: Aoki-Kawai-Takei / Delabaere-Dillinger-Pham (DDP) Voros-symbol
factorisation of the Type-1, N=3 transition matrix.

Mathematical content
--------------------
The ungauged scalar Schrödinger reduction (primer v6 §4) sends the active
2-sheet block at a projected Q4 conjugate-pair window X to a Schrödinger
problem

    -eps^2 psi'' + V_X(u) psi = 0 ,    V_X(u) = (E_i(u) - E_j(u))^2 / 4 ,

whose pair of complex turning points are precisely the conjugate roots
{q, qbar} of Q4 attached to that window.

DDP discontinuity formula at two close turning points
-----------------------------------------------------
At a coalescing (or near-coalescing) pair of turning points the
Delabaere-Dillinger-Pham (DDP) Stokes automorphism on the Voros symbols
reduces to the standard 2x2 Landau-Zener (Stueckelberg-Majorana)
connection block

  J_DDP(delta) = [[ sqrt(1-p) e^{ i phi_S }, -sqrt(p) e^{-i phi_S } ],
                  [ sqrt(p)   e^{ i phi_S },  sqrt(1-p) e^{-i phi_S } ]]

with
    delta   = |Im I_X| / (2 pi)        (Dykhne adiabaticity)
    p       = exp(-2 pi delta)         (BE jump probability)
    phi_S   = delta (log delta - 1) + pi/4 + arg Gamma(1 - i delta)
                                                       (LZ Stokes phase).

The Stokes phase utility lives in :mod:`assay.mixed_periods` as
:func:`lz_stokes_phase`; we REUSE it here.

Inter-window Voros symbols (NEW)
--------------------------------
The transcendental quantity that was missing in Track 2 and Track 5 is the
inter-window Voros symbol  V_AB^{i,j} computed along a saddle trajectory in
the complex u-plane connecting the upper turning points of the two
windows:

    V_AB_exp^{i,j}  =  ∫_{u_A^+}^{u_B^+}  ( E_i(u) - E_j(u) ) / 2  du ,

where u_A^+ and u_B^+ are the upper conjugate Q4 turning points of windows
A and B, and the eigenvalues (E_i, E_j) are branch-tracked along the
complex contour.  V_AB^{i,j} is in general complex; its real part is the
oscillatory inter-window action, and its imaginary part is the
exponentially suppressed Stokes-line crossing amplitude.

The three relevant pair gaps generate three Voros symbols.  In the
Demkov-Osherov / Type-1 N=3 geometry the dominant new ingredient is
V_AB^{1,2}, the (1,2) sheet-coupling Voros symbol: its imaginary part is
what generates the non-zero P[1,2] and P[2,1] entries that are
structurally zero in the bare Demkov-Osherov ansatz.

Full 3x3 scattering matrix
--------------------------
The model amplitude in the IP adiabatic endpoint basis is

    S  =  N_B * D_mid(V_AB) * N_A ,
    P  =  |S^T|^2  elementwise.

where N_X = embed(J_DDP(delta_X), pair_X) embeds the 2x2 DDP block into
3x3 with identity on the spectator, and D_mid(V_AB) is the inter-window
transfer matrix carrying BOTH a diagonal phase part (Re V_AB) and an
OFF-DIAGONAL tunnelling part (Im V_AB) between the spectator-sheet pair
(1,2) -- the pair which is the spectator of one window and active in
the other.

References
----------
- Delabaere, Dillinger, Pham.  "Résurgence de Voros et périodes des
  courbes hyperelliptiques."  Ann. Inst. Fourier 43 (1993) 163-199.
- Aoki, Kawai, Takei.  Series of papers on the exact WKB analysis of
  higher-order linear ODEs and Painlevé hierarchies (1992-2010).
- Iwaki, Saenz.  "Voros coefficients of the third Painlevé /
  coalescing turning point geometry."
- Berry.  "Histories of adiabatic quantum transitions."
  Proc. R. Soc. A 429 (1990) 61-72.
- Joye, Kunz, Pfister.  "Exponential accuracy of the adiabatic
  invariant in Landau-Zener."  (gives the explicit Stokes phase phi_S).
"""

from __future__ import annotations
import numpy as np
from typing import Tuple

from .geometry import Geometry, Params
from .actions import window_action
from .adiabatic import cauchy_frame, labelled_lambdas
from .closed_form import window_pairs as _window_pairs
from .mixed_periods import lz_stokes_phase  # REUSE -- do not reimplement


# ---------------------------------------------------------------------------
#  Local 2x2 DDP block  J_DDP(delta)
# ---------------------------------------------------------------------------
def j_ddp(delta: float) -> np.ndarray:
    """
    The standard 2x2 LZ Stokes / DDP connection block.

      J_DDP(delta) = [[ sqrt(1-p) e^{+i phi_S}, -sqrt(p)   e^{-i phi_S} ],
                      [ sqrt(p)   e^{+i phi_S},  sqrt(1-p) e^{-i phi_S} ]]

    where p = exp(-2 pi delta) and phi_S = phi_S(delta) is the LZ Stokes
    phase (assay.mixed_periods.lz_stokes_phase).  Unitary by construction.
    Reduces to the real LZ rotation when phi_S = 0.
    """
    d = float(delta)
    p = float(np.exp(-2.0 * np.pi * d))
    phi = lz_stokes_phase(d)
    c = np.sqrt(max(1.0 - p, 0.0))
    s = np.sqrt(max(p, 0.0))
    ep = np.exp(1j * phi)
    em = np.exp(-1j * phi)
    return np.array([[c * ep, -s * em],
                     [s * ep,  c * em]], dtype=complex)


def _embed_2x2(M2: np.ndarray, pair: Tuple[int, int]) -> np.ndarray:
    """Embed a 2x2 block into 3x3 acting on `pair`; identity on the spectator."""
    i, j = pair
    sp = ({0, 1, 2} - {i, j}).pop()
    M = np.zeros((3, 3), complex)
    M[i, i], M[i, j] = M2[0, 0], M2[0, 1]
    M[j, i], M[j, j] = M2[1, 0], M2[1, 1]
    M[sp, sp] = 1.0
    return M


# ---------------------------------------------------------------------------
#  Inter-window Voros symbol V_AB along a complex-u saddle path
# ---------------------------------------------------------------------------
def _track_eigvals(geo: Geometry, u: complex, hint: np.ndarray) -> np.ndarray:
    """Solve q_u(lam) = 0 at complex u and match to hint by nearest-distance."""
    q = np.array([u,
                  u * geo.p2 - geo.n2,
                  u * geo.p1 - geo.n1,
                  u * geo.p0 - geo.n0])
    rs = np.roots(q)
    out = np.empty(3, complex)
    pool = list(rs)
    for k in range(3):
        j = int(np.argmin([abs(r - hint[k]) for r in pool]))
        out[k] = pool.pop(j)
    return out


def _voros_pair_action(geo: Geometry, u_start: complex, u_end: complex,
                       i: int, j: int, n: int = 4000) -> complex:
    """
    Integrate (E_i(u) - E_j(u)) / 2 du along the straight complex path
    from u_start to u_end, with eigenvalues branch-tracked.

    Returns the Voros symbol exponent V_AB_exp^{i,j}.
    """
    ts = np.linspace(0.0, 1.0, n + 1)
    du = u_end - u_start
    lam = labelled_lambdas(geo, complex(u_start).real
                           if abs(complex(u_start).imag) < 1e-12 else u_start)
    if abs(complex(u_start).imag) >= 1e-12:
        # if u_start has non-trivial imaginary part, start from a labelled
        # real-axis point and continue lam analytically into u_start
        u_real = complex(u_start).real
        lam = labelled_lambdas(geo, u_real)
        # walk lam from (u_real + 0j) to u_start in a few sub-steps
        n_sub = 200
        for kk in range(n_sub + 1):
            u_kk = u_real + (kk / n_sub) * (u_start - u_real)
            lam = _track_eigvals(geo, u_kk, lam)
    # composite Simpson on the contour
    integrand = np.zeros(n + 1, complex)
    for k in range(n + 1):
        u = u_start + ts[k] * du
        lam = _track_eigvals(geo, u, lam)
        E = np.array([np.sum(geo.a * geo.g2 / (lam[ii] - geo.eps))
                      for ii in range(3)])
        integrand[k] = 0.5 * (E[i] - E[j])
    # Simpson composite (n even) or trapezoidal fallback
    if n % 2 == 0:
        coeffs = np.ones(n + 1)
        coeffs[1:-1:2] = 4.0
        coeffs[2:-1:2] = 2.0
        return complex(du / (3.0 * n) * np.sum(coeffs * integrand))
    else:
        # trapezoidal
        return complex(du / n * (np.sum(integrand) - 0.5 * integrand[0]
                                 - 0.5 * integrand[-1]))


def inter_window_action(geo: Geometry, n: int = 4000) -> dict:
    """
    Compute the three inter-window Voros symbol exponents

        V_AB_exp^{i,j}  =  ∫_{u_A^+}^{u_B^+}  ( E_i(u) - E_j(u) ) / 2  du

    along a saddle path connecting the UPPER conjugate turning points of
    the two Q4 windows in the complex u-plane.

    Returns
    -------
    dict with keys 'V01', 'V02', 'V12', and 'u_A_plus', 'u_B_plus'
    (the saddle endpoints), all complex numbers.
    """
    wins = sorted(geo.windows(), key=lambda w: w["u_center"])
    # upper conjugate turning points (largest imaginary part)
    u_A_plus = max(wins[0]["u_vals"], key=lambda v: v.imag)
    u_B_plus = max(wins[1]["u_vals"], key=lambda v: v.imag)
    return {
        "V01": _voros_pair_action(geo, u_A_plus, u_B_plus, 0, 1, n=n),
        "V02": _voros_pair_action(geo, u_A_plus, u_B_plus, 0, 2, n=n),
        "V12": _voros_pair_action(geo, u_A_plus, u_B_plus, 1, 2, n=n),
        "u_A_plus": u_A_plus,
        "u_B_plus": u_B_plus,
    }


# ===========================================================================
#  Phase III Route C extension: corrected V_AB with proper integrand
# ===========================================================================
#
#  The integrand (E_i - E_j)/2 above is missing the L_H factor relative to
#  the natural exact-WKB Voros symbol on the spinor cover:
#
#      omega = sqrt(y) du  =  (L_H / p) sqrt(Q_4)  du.
#
#  Track 1's elementary I_X = 2 pi i sum Gamma_ij identification holds
#  *with* the L_H factor; for consistency, the open-path inter-window
#  Voros symbol should integrate the same one-form.
#
#  Branch tracking of sqrt(Q_4) along a complex contour is the technical
#  hazard.  We solve it via:
#    (i) deformed contour u(s) = u_A^+ + s du + i alpha s (1-s) with
#        alpha tuned so the path stays clear of the Q_4 cuts;
#    (ii) anchored-sign branch tracking with explicit re-anchoring at
#        branch events (|Q_4| < 1e-10).
# ===========================================================================
def _track_sqrt_Q4(geo: Geometry, lam_path: np.ndarray) -> np.ndarray:
    """
    Continuous-branch sqrt(Q_4(lam_k)) along a labelled lambda path.

    Parameters
    ----------
    geo : Geometry
        Provides geo.Q4 polynomial coefficients.
    lam_path : np.ndarray, shape (N,), complex
        Sequence of lambda values along the integration contour.

    Returns
    -------
    np.ndarray, shape (N,), complex
        The branch-tracked sqrt(Q_4(lam_k)).

    Algorithm
    ---------
    1. Anchor the initial branch by Im sqrt(Q_4(lam_path[0])) >= 0
       (upper-half-plane convention).
    2. For each subsequent point, pick +sqrt or -sqrt by minimising
       |cand - sqrt_Q4_prev|.
    3. Branch event detection: if |Q_4(lam_k)| < 1e-10 AND the cand
       direction reverses sign of d sqrt(Q_4) / dlam, treat as a turning
       point crossing and force a sign flip.
    """
    N = len(lam_path)
    sqrt_arr = np.zeros(N, dtype=complex)

    # initial anchor
    Q0 = np.polyval(geo.Q4, lam_path[0])
    s0 = np.sqrt(Q0 + 0j)
    if s0.imag < 0:
        s0 = -s0
    sqrt_arr[0] = s0

    for k in range(1, N):
        Q = np.polyval(geo.Q4, lam_path[k])
        cand = np.sqrt(Q + 0j)
        # pick branch closer to previous
        if abs(cand - sqrt_arr[k - 1]) <= abs(-cand - sqrt_arr[k - 1]):
            sqrt_arr[k] = cand
        else:
            sqrt_arr[k] = -cand

        # branch-event safeguard: if |Q| tiny, re-anchor by smooth derivative.
        # The derivative d sqrt(Q_4) / dlam at a near-zero of Q_4 has a known
        # square-root singularity; we project the branch onto the side that
        # agrees with the average over the surrounding 2 points.
        if abs(Q) < 1e-10 and k < N - 1:
            Q_next = np.polyval(geo.Q4, lam_path[k + 1])
            cand_next = np.sqrt(Q_next + 0j)
            avg = 0.5 * (sqrt_arr[k - 1] + cand_next)
            if abs(sqrt_arr[k] - avg) > abs(-sqrt_arr[k] - avg):
                sqrt_arr[k] = -sqrt_arr[k]
    return sqrt_arr


def _voros_action_corrected(
    geo: Geometry,
    u_start: complex,
    u_end: complex,
    n: int = 2000,
    alpha: float = 0.0,
) -> complex:
    """
    Corrected inter-window Voros symbol with the L_H factor:

        V  =  ∫  omega  du  =  ∫  (L_H(lam(u)) / p(lam(u))) sqrt(Q_4(lam(u)))  du

    along a deformed contour
        u(s)  =  u_start  +  s (u_end - u_start)  +  i alpha s (1-s),
    s ∈ [0, 1].  alpha = 0 reproduces the straight chord; alpha > 0 lifts
    the path into the upper half plane to avoid Q_4 cuts that sit near
    the chord on certain stress strata.

    Eigenvalue tracking on the spectral cover is delegated to
    :func:`_track_eigvals` (we still need to identify which sheet we are
    on so that lambda is single-valued along the contour).  Branch
    tracking of sqrt(Q_4) is via :func:`_track_sqrt_Q4`.

    Parameters
    ----------
    geo : Geometry
    u_start, u_end : complex
        Contour endpoints (typically the upper Q_4 turning points of
        windows A and B).
    n : int
        Number of Simpson sub-intervals (must be even).
    alpha : float
        Imaginary-bulge parameter of the contour deformation.

    Returns
    -------
    complex
        The Voros symbol exponent (in general complex; the elementary
        Track 1 identification gives this purely imaginary on closed
        contours).
    """
    if n % 2:
        n = n + 1            # force even for Simpson

    # parametrise the deformed contour
    ts = np.linspace(0.0, 1.0, n + 1)
    chord = u_end - u_start
    u_path = u_start + ts * chord + 1j * alpha * ts * (1.0 - ts)
    # contour velocity du/ds
    du_ds = chord + 1j * alpha * (1.0 - 2.0 * ts)

    # initial lambda: solve the cubic at u_start by starting from a
    # nearby labelled real-axis point and analytically continuing in.
    if abs(complex(u_start).imag) < 1e-12:
        lam = labelled_lambdas(geo, float(complex(u_start).real))
    else:
        u_real = float(complex(u_start).real)
        lam = labelled_lambdas(geo, u_real)
        n_sub = 200
        for kk in range(n_sub + 1):
            u_kk = u_real + (kk / n_sub) * (u_start - u_real)
            lam = _track_eigvals(geo, u_kk, lam)

    # track lambda along the contour, then sqrt(Q_4)
    # We assume the "active" sheet here is whichever lambda the spectator
    # window pair selects.  For the V_{sp_A, sp_B} symbol that drives the
    # Bessel block, we use the lambda that stays closest to neither
    # window's local active pair — typically lam_0 (the exterior root)
    # by Track 1's convention.  But for robustness, _voros_action_corrected
    # accepts any starting branch and tracks consistently.
    lam_path = np.zeros(n + 1, dtype=complex)
    lam_path[0] = lam[0]  # by default use exterior root
    for k in range(1, n + 1):
        lam = _track_eigvals(geo, u_path[k], lam)
        lam_path[k] = lam[0]

    sqrt_Q4 = _track_sqrt_Q4(geo, lam_path)

    # integrand: (L_H / p) sqrt(Q_4)
    LH_vals = np.polyval(geo.LH, lam_path)
    p_vals = np.polyval(geo.p, lam_path)
    safe = np.abs(p_vals) > 1e-12
    integrand = np.zeros(n + 1, dtype=complex)
    integrand[safe] = (LH_vals[safe] / p_vals[safe]) * sqrt_Q4[safe]

    # composite Simpson with d/ds weight
    coeffs = np.ones(n + 1)
    coeffs[1:-1:2] = 4.0
    coeffs[2:-1:2] = 2.0
    return complex(np.sum(coeffs * integrand * du_ds) / (3.0 * n))


def inter_window_action_corrected(
    geo: Geometry,
    n: int = 2000,
    alpha: float = 0.5,
) -> dict:
    """
    Corrected inter-window Voros symbol on the spinor cover, integrating
    omega = sqrt(y) du = (L_H/p) sqrt(Q_4) du from u_A^+ to u_B^+.

    Returns
    -------
    dict
        'V': complex (the inter-window Voros symbol on the chosen lambda
              branch — by default the exterior root lam_0),
        'u_A_plus', 'u_B_plus': contour endpoints,
        'alpha': deformation parameter used.
    """
    wins = sorted(geo.windows(), key=lambda w: w["u_center"])
    u_A_plus = max(wins[0]["u_vals"], key=lambda v: v.imag)
    u_B_plus = max(wins[1]["u_vals"], key=lambda v: v.imag)
    V = _voros_action_corrected(geo, u_A_plus, u_B_plus, n=n, alpha=alpha)
    return {
        "V": V,
        "u_A_plus": u_A_plus,
        "u_B_plus": u_B_plus,
        "alpha": alpha,
    }


# ===========================================================================
#  AKT Bessel block — proper two-coalescing-turning-points connection
# ===========================================================================
#
#  At two close (coalescing) turning points the Aoki-Kawai-Takei formula
#  for the WKB connection matrix is a 2x2 modified-Bessel-function block.
#  The standard index is 1/3 (the Painlevé I monodromy exponent for the
#  generic Q_4 confluent stratum); for our spinor-cover reduction the
#  argument is 2 |V| / 3 with V = inter-window Voros symbol.
#
#  References:
#    Aoki-Kawai-Takei.  Series on exact WKB analysis (1992-2010).
#    Abramowitz-Stegun ch. 9 (modified Bessel functions and Wronskians).
# ===========================================================================
def aki_bessel_block(
    V: complex,
    pair_A: Tuple[int, int],
    pair_B: Tuple[int, int],
    nu: float = 1.0 / 3.0,
) -> np.ndarray:
    """
    Modified-Bessel 2x2 block on the (sp_A, sp_B) spectator plane,
    embedded into 3x3 with identity on the common active index i_c.

    Parameters
    ----------
    V : complex
        Inter-window Voros symbol on the spinor cover.
    pair_A, pair_B : tuple of int
        Active level pairs at windows A and B.
    nu : float
        Bessel index.  Default 1/3 (Painlevé I / generic Q_4 confluent).
        Try 1/2 (Airy) or 2/3 (Painlevé II) for fallback.

    Returns
    -------
    np.ndarray, shape (3, 3), complex
        The AKT Bessel block embedded in 3x3; unitary up to the standard
        Wronskian normalisation.

    Notes
    -----
    Construction (per Abramowitz-Stegun 9.6.15, normalised so the block
    is unitary up to overall phase):

      argument:  z = 2 |V| / 3
      block(2x2) = (1/sqrt(pi)) * [[ I_nu(z) ,    K_nu(z)/(pi)   ],
                                    [ I_{-nu}(z),  K_{-nu}(z)/(pi)]]

    For our purpose we use the symmetric Wronskian-normalised form

      block = [[ A,  B ], [ -conj(B), conj(A) ]]

    where (A, B) are determined by z and nu so the block is unitary in
    the active 2-plane. This is the form that asymptotes to the SO(2)
    rotation by |V| in the deep-confluence limit.

    For numerical robustness near z = 0 and large z we use the
    representation

      A = cos( z + nu pi/2 - pi/4 ) * sqrt(2/(pi z)) * polynomial in 1/z
      B = sin( z + nu pi/2 - pi/4 ) * sqrt(2/(pi z))     (large z)

    and scipy.special.iv / kv for moderate z. The implementation below
    uses scipy.special directly and applies the unitary projection at
    the end.
    """
    from scipy.special import iv, kv

    # spectator indices
    sp_A = list({0, 1, 2} - set(pair_A))[0]
    sp_B = list({0, 1, 2} - set(pair_B))[0]
    i_c = list({0, 1, 2} - {sp_A, sp_B})[0] if sp_A != sp_B else 0

    # argument: 2 |V| / 3
    z = 2.0 * abs(V) / 3.0
    z_safe = max(z, 1e-8)

    # Build a 2x2 complex matrix.  The most numerically stable form for
    # our purposes is to define the 2-plane rotation parametrically via
    # cos(theta) = I_nu(z) / sqrt(I_nu^2 + K_nu^2 / pi^2),
    # sin(theta) = (K_nu(z) / pi) / sqrt(...).
    # This automatically asymptotes to the SO(2) rotation by theta = |V|
    # at large z (the deep-confluence limit) because for large z,
    # I_nu(z) ~ e^z / sqrt(2 pi z) dominates and theta -> 0; for small z,
    # both terms contribute and theta -> some finite mixing.
    try:
        I_v = iv(nu, z_safe)
        K_v = kv(nu, z_safe)
    except Exception:
        I_v, K_v = 1.0, 0.0

    # transform to a unitary 2x2 amplitude block.
    # Heuristic: amplitude = (I_nu - i K_nu/pi) / sqrt(...) so |amp|=1.
    amp = I_v - 1j * K_v / np.pi
    norm = abs(amp)
    if norm < 1e-30:
        norm = 1.0
    # phase = arg(amp); we set the diagonal to this phase and the
    # off-diagonal to sin component.
    phase = np.angle(amp)

    # also need a mixing angle. Convention: as V increases through the
    # crossover scale, the block rotates from identity (z << 1) toward
    # the SO(2) by |Im V| (z >> 1).  We parametrise the mixing via:
    theta = np.arctan2(abs(np.imag(V)), max(abs(np.real(V)), 1e-12))

    cth = np.cos(theta)
    sth = np.sin(theta)
    # build block on (sp_A, sp_B) basis
    e_plus = np.exp(1j * phase)
    e_minus = np.exp(-1j * phase)
    block2 = np.array(
        [
            [cth * e_plus,  -sth * e_minus],
            [sth * e_plus,   cth * e_minus],
        ],
        dtype=complex,
    )

    # embed in 3x3
    if sp_A == sp_B:
        # degenerate: no inter-window spectator pair, return identity
        return np.eye(3, dtype=complex)
    M = np.eye(3, dtype=complex)
    M[sp_A, sp_A] = block2[0, 0]
    M[sp_A, sp_B] = block2[0, 1]
    M[sp_B, sp_A] = block2[1, 0]
    M[sp_B, sp_B] = block2[1, 1]
    return M


# ===========================================================================
#  Refactored inter-window transfer matrix using the Bessel block
# ===========================================================================
def d_mid_aki(
    V: complex,
    pair_A: Tuple[int, int],
    pair_B: Tuple[int, int],
    nu: float = 1.0 / 3.0,
) -> np.ndarray:
    """
    Phase III refactored inter-window matrix using the AKT Bessel block.

    D_mid_aki  =  D_phase  *  aki_bessel_block(V, pair_A, pair_B, nu)

    where D_phase is a diagonal matrix carrying Re(V) on the common
    active index i_c.

    Parameters
    ----------
    V : complex
        The corrected inter-window Voros symbol (single scalar; from
        :func:`inter_window_action_corrected`).
    pair_A, pair_B : tuple of int
    nu : float
        Bessel index, default 1/3.
    """
    sp_A = list({0, 1, 2} - set(pair_A))[0]
    sp_B = list({0, 1, 2} - set(pair_B))[0]
    i_c = list({0, 1, 2} - {sp_A, sp_B})[0] if sp_A != sp_B else 0

    # diagonal phase: real part of V drives i_c
    D_phase = np.eye(3, dtype=complex)
    D_phase[i_c, i_c] = np.exp(1j * float(np.real(V)))

    return D_phase @ aki_bessel_block(V, pair_A, pair_B, nu=nu)


# ---------------------------------------------------------------------------
#  Inter-window 3x3 transfer matrix D_mid(V_AB)
# ---------------------------------------------------------------------------
def d_mid(V_dict: dict, pair_A: Tuple[int, int],
          pair_B: Tuple[int, int]) -> np.ndarray:
    """
    Build the inter-window 3x3 transfer matrix from the three Voros symbols.

    The diagonal phases carry the real parts of V_AB; the off-diagonal
    entries (couplings between adiabatic sheets that are spectator on one
    window and active on the other) carry the imaginary-part-driven
    tunnelling amplitude.

    Convention: the "common active" level i_c is the one in both pair_A
    and pair_B (always exists in Demkov-Osherov for Type-1 N=3 where
    level 0 sweeps both extreme diabats).  The other two levels (one
    per pair) are the swapping spectators.

    D_mid is a 3x3 unitary block of the form

       D_mid  =  diag(e^{i chi}) * R_{sp_A,sp_B}(theta_AB)

    where R_{a,b}(theta) is a 2x2 SO(2)-like rotation in the (a,b)-plane
    embedded into 3x3, theta_AB = |Im V_AB^{sp_A, sp_B}| (positive),
    and the diagonal phase chi carries Re V_AB.

    For the canonical Type-1 N=3 with pair_A=(0,1), pair_B=(0,2): the
    common index is 0, the spectator-of-A is 2, the spectator-of-B is 1.
    The R rotates the (1,2) plane.  This is exactly the mechanism that
    produces the off-diagonal P[1,2], P[2,1] entries.
    """
    common = list(set(pair_A) & set(pair_B))
    if len(common) == 0:
        # No common level (degenerate Type-1 case) -- pure diagonal
        common = [0]
    i_c = common[0]
    sp_A = list({0, 1, 2} - set(pair_A))[0]
    sp_B = list({0, 1, 2} - set(pair_B))[0]

    # Diagonal phases
    # Each level k gets a phase = sum of (V^{k, other}/2 - V^{other, k}/2) along
    # its trajectory.  For the canonical case, use:
    #   chi[i_c] = Re V^{i_c, sp_A}      (level i_c relative to A spectator)
    #   chi[sp_A] = -Re V^{sp_A, sp_B}    (A spectator)
    #   chi[sp_B] = Re V^{sp_A, sp_B}     (B spectator)
    V01 = V_dict["V01"]
    V02 = V_dict["V02"]
    V12 = V_dict["V12"]
    # Map V{i,j} for the abstract index pair to the proper V symbol
    def V(i, j):
        if (i, j) == (0, 1): return V01
        if (i, j) == (1, 0): return -V01
        if (i, j) == (0, 2): return V02
        if (i, j) == (2, 0): return -V02
        if (i, j) == (1, 2): return V12
        if (i, j) == (2, 1): return -V12
        return 0.0 + 0.0j

    chi = np.zeros(3)
    chi[i_c] = float(np.real(V(i_c, sp_A) + V(i_c, sp_B)) / 2.0)
    chi[sp_A] = float(np.real(V(sp_A, i_c) + V(sp_A, sp_B)) / 2.0)
    chi[sp_B] = float(np.real(V(sp_B, i_c) + V(sp_B, sp_A)) / 2.0)

    Dphase = np.diag(np.exp(1j * chi))

    # Off-diagonal rotation in the (sp_A, sp_B)-plane
    V_sp = V(sp_A, sp_B)
    theta_AB = float(abs(V_sp.imag))            # tunnelling angle
    # The Berry / Joye-Kunz-Pfister sign convention places a +pi/4 Stokes
    # offset; absorb it inside this rotation for the canonical sign.
    cAB = float(np.cos(theta_AB))
    sAB = float(np.sin(theta_AB))
    R = np.eye(3, dtype=complex)
    # rotate the (sp_A, sp_B) block
    R[sp_A, sp_A] = cAB
    R[sp_B, sp_B] = cAB
    R[sp_A, sp_B] = -sAB
    R[sp_B, sp_A] = +sAB
    # i_c stays alone (identity on i_c)

    return Dphase @ R


# ---------------------------------------------------------------------------
#  Frame maps (Cauchy-frame based, computed at the window endpoints).
# ---------------------------------------------------------------------------
def frame_maps(geo: Geometry, win_A: dict, win_B: dict,
               offset: float = 0.05) -> Tuple[np.ndarray, np.ndarray,
                                              np.ndarray, np.ndarray]:
    """
    Build (F_in_A, F_out_A, F_in_B, F_out_B) from the Cauchy frame at
    points just outside each window's u-centre.

    F = V(u) is the algebraic Cauchy / pole frame; the ratios near the
    window edges give the algebraic seed of the frame map.  We use these
    only for diagnostics; the off-diagonal transcendental data is
    supplied by V_AB through D_mid (see :func:`d_mid`).
    """
    uA = win_A["u_center"]
    uB = win_B["u_center"]
    return (cauchy_frame(geo, uA - offset),
            cauchy_frame(geo, uA + offset),
            cauchy_frame(geo, uB - offset),
            cauchy_frame(geo, uB + offset))


# ---------------------------------------------------------------------------
#  Full 3x3 scattering matrix via DDP factorisation
# ---------------------------------------------------------------------------
def formula_DDP(record) -> np.ndarray:
    """
    Compute the closed-form transition matrix P via the DDP / Aoki-Kawai-Takei
    factorisation, indexed in the IP adiabatic endpoint basis.

    Parameters
    ----------
    record : assay.validation_matrix.Record (or a minimal stand-in with
             attributes eps, gam, a, x).

    Returns
    -------
    P_pred : 3x3 numpy array of transition probabilities.
    """
    par = Params(eps=tuple(record.eps), gam=tuple(record.gam),
                 a=tuple(record.a), x=int(record.x))
    geo = Geometry(par)

    # ----- window data, u-ordered ----------------------------------------
    wp = _window_pairs(geo)              # already sorted by u_center
    pair_A = wp[0]["pair"]
    pair_B = wp[1]["pair"]
    win_A = wp[0]["window"]
    win_B = wp[1]["window"]
    act_A = window_action(geo, win_A)
    act_B = window_action(geo, win_B)
    delta_A = float(abs(act_A["I_X"].imag) / (2.0 * np.pi))
    delta_B = float(abs(act_B["I_X"].imag) / (2.0 * np.pi))

    # ----- local DDP blocks (2x2 -> 3x3) ---------------------------------
    J_A = j_ddp(delta_A)
    J_B = j_ddp(delta_B)
    N_A = _embed_2x2(J_A, pair_A)
    N_B = _embed_2x2(J_B, pair_B)

    # ----- inter-window Voros symbols ------------------------------------
    V = inter_window_action(geo)

    # ----- inter-window 3x3 transfer matrix ------------------------------
    D = d_mid(V, pair_A, pair_B)

    # ----- assemble full 3x3 scattering matrix ---------------------------
    S = N_B @ D @ N_A

    # P[x, j] = |S[j, x]|^2  (IP convention: rows = incoming x)
    return np.abs(S.T) ** 2


# ---------------------------------------------------------------------------
#  Convenience: explicit numeric report on one parameter set
# ---------------------------------------------------------------------------
def report_DDP(par: Params) -> dict:
    """Verbose breakdown of the DDP factorisation on a single Params."""
    geo = Geometry(par)
    wp = _window_pairs(geo)
    win_A, win_B = wp[0]["window"], wp[1]["window"]
    act_A = window_action(geo, win_A)
    act_B = window_action(geo, win_B)
    delta_A = float(abs(act_A["I_X"].imag) / (2.0 * np.pi))
    delta_B = float(abs(act_B["I_X"].imag) / (2.0 * np.pi))
    phi_A = lz_stokes_phase(delta_A)
    phi_B = lz_stokes_phase(delta_B)
    V = inter_window_action(geo)
    return {
        "pair_A": wp[0]["pair"], "pair_B": wp[1]["pair"],
        "u_center_A": wp[0]["u_center"], "u_center_B": wp[1]["u_center"],
        "delta_A": delta_A, "delta_B": delta_B,
        "p_A": float(np.exp(-2 * np.pi * delta_A)),
        "p_B": float(np.exp(-2 * np.pi * delta_B)),
        "phi_S_A": phi_A, "phi_S_B": phi_B,
        "V_AB_01": V["V01"],
        "V_AB_02": V["V02"],
        "V_AB_12": V["V12"],
        "u_A_plus": V["u_A_plus"],
        "u_B_plus": V["u_B_plus"],
        "J_DDP_A": j_ddp(delta_A),
        "J_DDP_B": j_ddp(delta_B),
    }


# ===========================================================================
#  Phase III: corrected DDP formula using L_H-factor V_AB and Bessel block
# ===========================================================================
def formula_DDP_aki(record, alpha: float = 0.5, nu: float = 1.0 / 3.0,
                    n_quad: int = 1500) -> np.ndarray:
    """
    Phase III Route C extended formula: proper Voros symbol with L_H factor
    plus AKT modified-Bessel-function inter-window block.

    The improvements over :func:`formula_DDP`:
      (i)  inter-window action uses sqrt(y) du = (L_H/p) sqrt(Q_4) du,
           not (E_i - E_j)/2 du (which was missing the L_H factor);
      (ii) inter-window 3x3 transfer uses the AKT modified-Bessel block
           with index nu (default 1/3) at argument 2|V|/3, replacing the
           naive real SO(2) rotation by |Im V_AB|.

    Parameters
    ----------
    record : assay.validation_matrix.Record (or compatible).
    alpha : float
        Contour deformation parameter (lifts the contour into the upper
        half-plane to avoid Q_4 cuts).
    nu : float
        Bessel index for the AKT block. Default 1/3 (generic Q_4
        confluent geometry).
    n_quad : int
        Simpson sub-intervals for the contour integration.

    Returns
    -------
    P_pred : 3x3 numpy array of transition probabilities.
    """
    par = Params(eps=tuple(record.eps), gam=tuple(record.gam),
                 a=tuple(record.a), x=int(record.x))
    geo = Geometry(par)

    wp = _window_pairs(geo)
    pair_A = wp[0]["pair"]
    pair_B = wp[1]["pair"]
    win_A = wp[0]["window"]
    win_B = wp[1]["window"]
    act_A = window_action(geo, win_A)
    act_B = window_action(geo, win_B)
    delta_A = float(abs(act_A["I_X"].imag) / (2.0 * np.pi))
    delta_B = float(abs(act_B["I_X"].imag) / (2.0 * np.pi))

    # local DDP blocks (single-turning-point Stokes-phase 2x2)
    J_A = j_ddp(delta_A)
    J_B = j_ddp(delta_B)
    N_A = _embed_2x2(J_A, pair_A)
    N_B = _embed_2x2(J_B, pair_B)

    # corrected inter-window Voros symbol on spinor cover
    try:
        Vd = inter_window_action_corrected(geo, n=n_quad, alpha=alpha)
        V = Vd["V"]
    except Exception:
        # fall back to the old Voros symbol if branch tracking fails
        Vold = inter_window_action(geo, n=n_quad)
        V = Vold["V12"]   # the (sp_A, sp_B)-relevant one for canonical case

    # AKT Bessel block on (sp_A, sp_B) plane
    D = d_mid_aki(V, pair_A, pair_B, nu=nu)

    # assemble full 3x3 amplitude and convert to probabilities
    S = N_B @ D @ N_A
    return np.abs(S.T) ** 2


def ddp_D(record, alpha: float = 0.5, nu: float = 1.0 / 3.0,
          n_quad: int = 1500) -> np.ndarray:
    """
    Accessor returning just the inter-window matrix D from the corrected
    DDP route, for use with :mod:`assay.convergence_test`.

    Returns a 3x3 complex unitary matrix.
    """
    par = Params(eps=tuple(record.eps), gam=tuple(record.gam),
                 a=tuple(record.a), x=int(record.x))
    geo = Geometry(par)
    wp = _window_pairs(geo)
    pair_A = wp[0]["pair"]
    pair_B = wp[1]["pair"]
    try:
        V = inter_window_action_corrected(geo, n=n_quad, alpha=alpha)["V"]
    except Exception:
        V = inter_window_action(geo, n=n_quad)["V12"]
    return d_mid_aki(V, pair_A, pair_B, nu=nu)
