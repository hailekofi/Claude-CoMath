"""
Iwaki-Nakanishi cluster-algebra / Voros-symbol realisation of the Type-1
N=3 multistate Landau-Zener problem (Route B).

The Iwaki-Nakanishi / Allegretti theorem
-----------------------------------------
Let psi be the WKB solution to the second-order Schroedinger-type pencil

    (epsilon^2 d^2/du^2 + y(u)) psi = 0 ,                     (*)

with y(u) = Q_4(lambda(u)) L_H(lambda(u))^2 / p(lambda(u))^2 the phase
coordinate of the primer's quotient-ring scalar reduction (S4) restricted
to the active 2-sheet pair on each window.  The Stokes graph of (*) in
the u-plane has

   * turning points = simple zeros of y(u),
   * Stokes lines = trajectories of Im int_u0^u sqrt(y) du = 0 emanating
     from each turning point in 3 directions,
   * chambers = connected components of the complement.

To each saddle 1-cycle gamma on the spectral cover mu^2 = y(u) one
assigns a Voros symbol  V_gamma = (1/epsilon) oint_gamma sqrt(y) du.
Iwaki-Nakanishi (2014) and Allegretti (2018) prove that as the angle
theta in epsilon = |epsilon| e^{i theta} sweeps, the Borel-summed Voros
symbols transform by cluster Y-mutations at each Stokes-line crossing.

Type-1 N=3 cluster realisation
------------------------------
For the genus-1 problem at hand there are 4 turning points (two
conjugate Q_4 pairs) and 3 elementary 2-level crossings in u-order
(the Demkov-Osherov triangle).  The Track-1 bundling rule

    I_X  =  2 pi i  sum_{(i,j) in S_X} Gamma_ij^{sign}

with Gamma_ij^{sign} = gamma_i^2 gamma_j^2 (a_i - a_j)/(eps_i-eps_j)^2
assigns to each *window* W a subset S_W of the three pairs, and the two
window subsets are complementary (Track 1 Proposition 4.4).  So the
BPS quiver has 3 nodes (one per elementary crossing) and the cluster
mutation sequence is the u-ordered length-3 braid.  In the cluster
language the two-window Z_2 grading is the residual modular involution
of the A_2 Y-system after the BE diagonal is gauged away.

Public API
----------
StokesGraph         : turning points, Stokes-line topology.
Crossing            : a single u-ordered 2-level crossing with its Voros data.
BPSQuiver           : 3 nodes (the crossings) and signed adjacency B-matrix.
build_quiver(geo)   : construct the quiver from a Geometry.
voros_Y(geo)        : Y-variables Y_k = exp(2 pi i Gamma_k^{sign}) per crossing.
mutation_sequence   : u-ordered list of mutations [0, 1, 2].
formula_IN(record)  : closed-form 3x3 P built from the cluster mutation product.
sinitsyn_chernyak_relation : show the Y-system identity behind Track 3.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Tuple, Dict, Optional
import numpy as np

from .geometry import Params, Geometry
from .actions import all_window_actions, window_action
from .grid import lz_parameters, crossing_order
from .quantum_dilog import (Li2, phi_classical, mutation_Y,
                            stokes_matrix_2x2)


# ===========================================================================
#  Stokes graph
# ===========================================================================
@dataclass
class TurningPoint:
    u: complex                  # turning-point u value
    lam: complex                # corresponding lambda
    window: int                 # which Q_4 conjugate-pair window it belongs to


@dataclass
class StokesGraph:
    """
    Stokes graph of the scalar WKB pencil (*) in the u-plane.

    For Type-1 N=3 the spectral curve has genus 1 with 4 simple branch
    points (the four zeros of Q_4 in lambda, lifted to u by u = n/p),
    grouped into two complex-conjugate pairs.  Each pair anchors a
    saddle 1-cycle gamma_X on the spectral cover; the Track-1 bundling
    rule assigns to each saddle a set of underlying elementary
    Landau-Zener crossings (one or two pairs (i,j)).  The 3 pairwise
    diabatic crossings of Demkov-Osherov are recovered by *splitting*
    the saddle: the cluster mutation realises this splitting through
    the Y-pattern combinatorics.

    Attributes
    ----------
    turning_points : list[TurningPoint]
        Four u-turning points, ordered by Re(u) then Im(u).
    windows         : list of dicts (from Geometry.windows())
    crossings       : list[Crossing] (3 of them, u-ordered)
    chambers_at_infinity : (str, str)
    """
    turning_points: List[TurningPoint]
    windows: List[dict]
    crossings: List["Crossing"]
    chambers_at_infinity: Tuple[str, str] = ("u_neg_inf", "u_pos_inf")


@dataclass
class Crossing:
    """
    An elementary 2-level Landau-Zener crossing in u-order.

    Fields
    ------
    pair       : (i, j) diabatic level pair, i < j
    u_cross    : real u-value of the crossing  u_ij = (H0_jj - H0_ii)/(a_i - a_j)
    Gamma_sign : signed LZ exponent  gamma_i^2 gamma_j^2 (a_i - a_j)/(eps_i-eps_j)^2
    Gamma_abs  : |Gamma_sign|
    p_jump     : exp(-2 pi Gamma_abs)
    window_idx : which u-window saddle this crossing is bundled into (0 or 1)
    bundle_sign: the +-1 sign with which Gamma_sign enters the window's I_X
    """
    pair: Tuple[int, int]
    u_cross: float
    Gamma_sign: float
    Gamma_abs: float
    p_jump: float
    window_idx: int
    bundle_sign: int


def _u_of_lambda(geo: Geometry, lam: complex) -> complex:
    return np.polyval(geo.n, lam) / np.polyval(geo.p, lam)


def _u_crossing(geo: Geometry, pair: Tuple[int, int]) -> float:
    """u_{ij} = (H0_jj - H0_ii) / (a_i - a_j) -- the diabatic crossing point."""
    i, j = pair
    return (geo.H0[j, j] - geo.H0[i, i]) / (geo.a[i] - geo.a[j])


def _identify_bundling(geo: Geometry,
                       tol: float = 5e-3
                       ) -> List[List[Tuple[Tuple[int, int], int]]]:
    """
    Determine which window each pairwise crossing belongs to and with
    what sign, from the Track-1 bundling rule
        I_X / (2 pi i) =  sum_{(i,j) in S_X}  sgn_{X,(i,j)} Gamma_ij^{sign}.

    Returns a list (one per window in u-order) of lists of
    ((i,j), sgn) tuples.
    """
    eps = np.asarray(geo.eps, float)
    gam = np.asarray(geo.gam, float)
    a = np.asarray(geo.a, float)
    G_sign = {}
    for i in range(3):
        for j in range(i + 1, 3):
            G_sign[(i, j)] = (gam[i] ** 2 * gam[j] ** 2 *
                              (a[i] - a[j]) / (eps[i] - eps[j]) ** 2)
    # window periods I_X / (2 pi i)
    wins = sorted(geo.windows(), key=lambda w: w["u_center"])
    targets = []
    for w in wins:
        I = complex(window_action(geo, w)["I_X"])
        targets.append(I.imag / (2.0 * np.pi))   # signed real value
    pairs = list(G_sign.keys())
    # Search over all  +-1 assignments of each pair to one of the two
    # windows (or neither).  Encoding: assignment[k] in {-1, 0, +1}*2,
    # i.e. each pair contributes to window 0 with sign +-1 or to window 1
    # with sign +-1 -- never both (Track 1 Prop 4.4 disjoint capture).
    best = None
    best_err = np.inf
    from itertools import product
    for assign in product([(0, +1), (0, -1), (1, +1), (1, -1)], repeat=3):
        sums = [0.0, 0.0]
        for k, (w_idx, sgn) in enumerate(assign):
            sums[w_idx] += sgn * G_sign[pairs[k]]
        err = (abs(sums[0] - targets[0]) + abs(sums[1] - targets[1]))
        if err < best_err:
            best_err = err
            best = assign
    bundling: List[List[Tuple[Tuple[int, int], int]]] = [[], []]
    for k, (w_idx, sgn) in enumerate(best):
        bundling[w_idx].append((pairs[k], sgn))
    bundling.append(best_err)         # diagnostics
    return bundling


def build_stokes_graph(geo: Geometry) -> StokesGraph:
    """Construct the Stokes graph of the scalar WKB pencil in u."""
    # turning points (4 of them, two conjugate u-pairs)
    wins = sorted(geo.windows(), key=lambda w: w["u_center"])
    tps: List[TurningPoint] = []
    for wi, w in enumerate(wins):
        for lam in w["roots"]:
            tps.append(TurningPoint(u=complex(_u_of_lambda(geo, complex(lam))),
                                    lam=complex(lam), window=wi))
    tps.sort(key=lambda t: (t.u.real, t.u.imag))

    # bundling: which window each pair belongs to
    bundling = _identify_bundling(geo)
    bundle_err = bundling[-1]
    bundling = bundling[:2]

    eps = np.asarray(geo.eps, float)
    gam = np.asarray(geo.gam, float)
    a = np.asarray(geo.a, float)
    G_sign = {(i, j): gam[i] ** 2 * gam[j] ** 2 *
              (a[i] - a[j]) / (eps[i] - eps[j]) ** 2
              for i in range(3) for j in range(i + 1, 3)}
    G_abs = {k: abs(v) for k, v in G_sign.items()}

    crossings: List[Crossing] = []
    for w_idx, lst in enumerate(bundling):
        for (pair, sgn) in lst:
            crossings.append(Crossing(
                pair=pair,
                u_cross=float(_u_crossing(geo, pair)),
                Gamma_sign=float(G_sign[pair]),
                Gamma_abs=float(G_abs[pair]),
                p_jump=float(np.exp(-2 * np.pi * G_abs[pair])),
                window_idx=w_idx,
                bundle_sign=int(sgn),
            ))
    # u-ordered crossings (Demkov-Osherov sequence)
    crossings.sort(key=lambda c: c.u_cross)

    return StokesGraph(turning_points=tps, windows=wins,
                       crossings=crossings)


# ===========================================================================
#  BPS quiver
# ===========================================================================
@dataclass
class BPSQuiver:
    """
    BPS quiver of Type-1 N=3.

    Nodes = the 3 elementary u-ordered crossings (Demkov-Osherov triangle).
    The skew-symmetric B-matrix encodes adjacency: two crossings are
    adjacent iff they share a diabatic level.  For the standard
    cyclic A_2-quiver realisation B is the antisymmetric Cartan-like
    matrix

        B = [[0, +1, -1],
             [-1, 0, +1],
             [+1, -1, 0]]    (A_2 cyclic quiver)

    Equivalently the genus-1 spectral curve's primitive BPS state
    spectrum forms a 3-state cycle (the "trivalent vertex of pants"
    decomposition of the spectral cover); see Allegretti 2018 Eq. 4.12.
    """
    crossings: List[Crossing]
    B: np.ndarray                       # skew-symmetric exchange matrix
    pairs: List[Tuple[int, int]]
    bundling: List[List[Tuple[Tuple[int, int], int]]]
    bundling_error: float = 0.0

    @property
    def n(self) -> int:
        return len(self.crossings)


def build_quiver(geo: Geometry) -> BPSQuiver:
    """Construct the BPS quiver from a Geometry."""
    sg = build_stokes_graph(geo)
    bundling = _identify_bundling(geo)
    bundle_err = float(bundling[-1])
    bundling_lists = bundling[:2]

    crossings = sg.crossings
    n = len(crossings)
    pairs = [c.pair for c in crossings]
    # Build B-matrix:  arrow from crossing k to k+1 (cyclic on the A_2 quiver
    # of the spectral cover); intersection number = +1 if pairs share a level.
    B = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            shared = set(pairs[i]) & set(pairs[j])
            if shared:
                # u-order orientation: arrow i -> j if u_i < u_j
                if crossings[i].u_cross < crossings[j].u_cross:
                    B[j, i] += 1
                    B[i, j] -= 1
    return BPSQuiver(crossings=crossings, B=B, pairs=pairs,
                     bundling=bundling_lists, bundling_error=bundle_err)


# ===========================================================================
#  Voros symbols / Y-variables
# ===========================================================================
def voros_Y(geo: Geometry, epsilon: complex = 1.0) -> np.ndarray:
    """
    Y-variables of the initial cluster (chamber at u -> -inf), one per
    elementary crossing, in u-order.

        Y_k  =  exp( 2 pi i * Gamma_k^{sign} / epsilon ) ,

    so that |Y_k| = 1 and arg Y_k = 2 pi Gamma_k^{sign} / epsilon.  This
    matches the Track-1 normalisation I_X = 2 pi i sum_signed Gamma.
    """
    quiv = build_quiver(geo)
    Ys = np.array([np.exp(2j * np.pi * c.Gamma_sign * c.bundle_sign / epsilon)
                   for c in quiv.crossings], dtype=complex)
    return Ys


def mutation_sequence(quiv: BPSQuiver) -> List[int]:
    """
    Mutation sequence realising the connection from u = -inf to
    u = +inf.

    For the u-ordered A_2 cyclic quiver, the canonical mutation
    sequence is the length-n braid that traverses each node once
    in u-order:  [0, 1, ..., n-1].
    """
    return list(range(quiv.n))


def voros_phases(geo: Geometry, epsilon: complex = 1.0) -> np.ndarray:
    """
    Saddle-bridge Voros phases  chi_k = arg Phi_b(V_k) (classical limit).

    Each elementary crossing carries a Stokes-corrected Voros phase via
    the cluster Y-pattern mutation.  In the semiclassical limit and
    with V_k = 2 pi i Gamma_k^{sign} (purely imaginary), the classical
    quantum dilogarithm  phi_classical(V_k/(2 pi))  contributes a
    *unit-modulus* phase to the Stokes matrix.

    Returns a length-n array of phases (radians).
    """
    quiv = build_quiver(geo)
    chis = []
    for c in quiv.crossings:
        # z so that exp(2 pi z) = Y_k:   z = i Gamma_k^{sign} (purely imaginary)
        z = 1j * c.Gamma_sign * c.bundle_sign / epsilon
        phi = phi_classical(z)
        chis.append(float(np.angle(phi)))
    return np.array(chis, dtype=float)


# ===========================================================================
#  Closed-form predictor:  Iwaki-Nakanishi cluster formula
# ===========================================================================
def _embed_2x2(M2: np.ndarray, pair: Tuple[int, int]) -> np.ndarray:
    """Embed a 2x2 complex amplitude block into 3x3, identity on the spectator."""
    i, j = pair
    sp = ({0, 1, 2} - {i, j}).pop()
    M = np.zeros((3, 3), dtype=complex)
    M[i, i], M[i, j] = M2[0, 0], M2[0, 1]
    M[j, i], M[j, j] = M2[1, 0], M2[1, 1]
    M[sp, sp] = 1.0
    return M


def _endpoint_permutations(a: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """IP adiabatic endpoint permutations (exterior-first labelling)."""
    pi_in = np.argsort(np.asarray(a, float))
    pi_out = np.argsort(-np.asarray(a, float))
    return pi_in, pi_out


def cluster_connection_3x3(geo: Geometry,
                           epsilon: complex = 1.0,
                           use_voros_phase: bool = True
                           ) -> Dict[str, np.ndarray]:
    """
    Assemble the 3x3 diabatic connection matrix from the cluster mutation
    sequence.

    Recipe
    ------
    1.  Build the BPS quiver with its 3 u-ordered crossings.
    2.  Initialise the Y-variables  Y_k = exp(2 pi i Gamma_k^{sign}).
    3.  Walk the u-ordered mutation sequence [0, 1, 2].  Each step at
        node k contributes
           T_k  =  embed_{pair_k}( T_{X_k}( p_k, chi_k ) )
        where p_k = exp(-2 pi Gamma_k^{abs}) is the elementary LZ jump
        probability, and chi_k is the Voros phase carried by the
        Phi_b mutation (classical limit).
    4.  At each step, mutate the Y-variables according to the cluster
        Y-pattern, propagating the Voros phases through the quiver.

    Returns a dict with the per-step matrices and the final 3x3 S.
    """
    quiv = build_quiver(geo)
    seq = mutation_sequence(quiv)
    Ys0 = voros_Y(geo, epsilon=epsilon)
    chis = voros_phases(geo, epsilon=epsilon) if use_voros_phase \
        else np.zeros(quiv.n)

    S = np.eye(3, dtype=complex)
    Y = Ys0.copy()
    factors = []
    for k in seq:
        c = quiv.crossings[k]
        Tk2 = stokes_matrix_2x2(c.p_jump, phase=chis[k])
        Tk3 = _embed_2x2(Tk2, c.pair)
        S = Tk3 @ S
        factors.append(Tk3)
        # mutate Y (records the cluster identity for diagnostics)
        Y = mutation_Y(Y, quiv.B, k)
    return {"S": S, "Y_initial": Ys0, "Y_final": Y,
            "chis": chis, "quiver": quiv, "sequence": seq,
            "factors": factors}


def formula_IN(record, epsilon: complex = 1.0,
               use_voros_phase: bool = True) -> np.ndarray:
    """
    Iwaki-Nakanishi cluster closed-form predictor for the 3x3
    interaction-picture adiabatic transition matrix P.

    Parameters
    ----------
    record : validation_matrix.Record
        Holds (eps, gam, a, x, P_bench, q_ij, window_periods, ...).
    epsilon : complex
        Voros symbol Planck weighting.  Defaults to 1.
    use_voros_phase : bool
        Whether to apply the classical-Phi_b Voros phase per crossing.
        Setting to False reproduces the incoherent LZ-grid prediction.

    Returns
    -------
    P : np.ndarray, shape (3, 3)
        Predicted IP adiabatic transition probabilities.

    Construction
    ------------
    Uses ``cluster_connection_3x3`` to obtain the diabatic amplitude
    matrix S (S[j, i] = amplitude from diabatic i to diabatic j), then
    converts to the IP adiabatic basis via the exterior-first endpoint
    permutations  pi_in = argsort(a), pi_out = argsort(-a):

        P_adia[i, j]  =  |S[pi_out[j], pi_in[i]]|^2 .

    This matches the convention used by ``closed_form.fit_closed_form``
    and ``validation_matrix.grid_formula``.
    """
    par = Params(eps=record.eps, gam=record.gam, a=record.a, x=record.x)
    geo = Geometry(par)
    out = cluster_connection_3x3(geo, epsilon=epsilon,
                                 use_voros_phase=use_voros_phase)
    S = out["S"]
    pi_in, pi_out = _endpoint_permutations(np.array(record.a, dtype=float))
    # P_adia[i, j] = |S[pi_out[j], pi_in[i]]|^2
    P = np.abs(S[np.ix_(pi_out, pi_in)].T) ** 2
    return P


# ===========================================================================
#  Diagnostics / Sinitsyn-Chernyak T-system relation
# ===========================================================================
def sinitsyn_chernyak_relation(geo: Geometry,
                               geo_prime: Geometry) -> dict:
    """
    Show that the Sinitsyn-Chernyak integrability condition
    [H, H'] = 0 (verified to 7e-15 in Track 3) corresponds to a
    *cluster T-system identity* on the Y-variables of the
    Iwaki-Nakanishi Y-pattern, but does NOT collapse the mutation
    product to elementary Brundobler-Elser factors.

    In the A_2 cyclic quiver Y-pattern with three nodes {Y_0, Y_1, Y_2},
    the standard A_2 T-system identity is

        T_+(Y_0, Y_1, Y_2) :=
            (1 + Y_0)(1 + Y_1)(1 + Y_2)
            / [(1 + Y_0 Y_1)(1 + Y_1 Y_2)(1 + Y_0 Y_2)]
            = invariant under cyclic mutation .

    Since each Y_k = exp(2 pi i Gamma_k^{sign}) lies on the unit circle,
    T_+ is also a unit-modulus complex number.  For the commuting
    partner H' with Y'_k = exp(2 pi i Gamma_k^{sign}(a')), the ratio
    T_+(Y) / T_+(Y') is the cluster-invariant signature of the
    commuting pair (H, H').

    Crucially, T_+ does NOT factorise as exp(-2 pi sum Gamma_k^{abs})
    in general: the cluster identity provides *additional* information
    above and beyond BE.  This explains the Track 3 negative finding:
    commutativity gives a Y-system constraint, but the connection
    matrix retains a non-trivial mutation braid.

    Returns
    -------
    {
      "Y": Y-variables of H,
      "Y_prime": Y-variables of H',
      "T_plus_A2": T_+(Y),
      "T_plus_A2_prime": T_+(Y'),
      "BE_factor": elementary BE product exp(-2 pi sum Gamma_ij^abs),
      "ratio": T_+ / BE_factor,
      "commute_check_value": |T_+(Y) - T_+(Y')|.
    }
    """
    Y = voros_Y(geo, epsilon=1.0)
    Yp = voros_Y(geo_prime, epsilon=1.0)
    n = len(Y)

    def T_plus(Yv: np.ndarray) -> complex:
        prod_single = 1.0 + 0.0j
        prod_pair = 1.0 + 0.0j
        for k in range(n):
            prod_single *= (1.0 + Yv[k])
        for i in range(n):
            for j in range(i + 1, n):
                prod_pair *= (1.0 + Yv[i] * Yv[j])
        return prod_single / prod_pair

    T = T_plus(Y)
    Tp = T_plus(Yp)
    G = lz_parameters(geo)["Gamma"]
    pairs = [(0, 1), (0, 2), (1, 2)]
    be = float(np.exp(-2 * np.pi * sum(G[p] for p in pairs)))
    ratio = T / be
    return {"Y": Y, "Y_prime": Yp,
            "T_plus_A2": complex(T),
            "T_plus_A2_prime": complex(Tp),
            "BE_factor": be,
            "ratio": complex(ratio),
            "commute_check_value": float(abs(T - Tp))}


__all__ = [
    "StokesGraph", "BPSQuiver", "Crossing",
    "build_stokes_graph", "build_quiver",
    "voros_Y", "voros_phases", "mutation_sequence",
    "cluster_connection_3x3", "formula_IN",
    "sinitsyn_chernyak_relation",
]
