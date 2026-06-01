"""
WS-G: Stokes graph (exact-WKB / spectral-network) for the Type-1 N=3 MLZ problem.

Tests the JOINT diagnosis:  P is the ordered product of local Stokes jump matrices
along the Stokes graph on the complex u-plane.  N=2 -> one turning-point pair, the
answer factorizes (= the LZ formula).  N=3 -> Stokes lines from DIFFERENT pairwise
turning points can INTERSECT, forming a JOINT (a BPS junction, Gaiotto-Moore-Neitzke),
where the local data does NOT factorize into independent 2-level pieces.

CONJECTURE under test:
  well-separated (the ~15%, incoherent product exact)  =  joint-free Stokes graph,
  overlapping     (the ~85%, P_2->2 enhanced 2.5-104x)  =  Stokes graph WITH a joint,
  and the joint contribution = the off-diagonal Stokes connection coefficient.

MODEL:  H(u) = H0 + u*diag(a),  Type-1 Cauchy coupling.
  Schrodinger ODE  i dpsi/du = H(u) psi.
  Adiabatic eigenvalues E_i(u) = eig of H(u) (complex u).
  Pairwise turning points t_ij : E_i(u)=E_j(u) (branch points of Disc_E chi_H).
  Stokes lines of type ij from t_ij : Im int_{t_ij}^u (E_i-E_j) du' = 0  (120 deg apart).
  Joint = intersection of an ij Stokes line with a kl (kl != ij) Stokes line
          away from any turning point.

Outputs:  experiments/figs/stokes_*.png  and a correlation table to stdout.

Run:  python3 experiments/stokes_graph.py
Deps: numpy, scipy, sympy, matplotlib.
"""
import os
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(HERE, "figs")
os.makedirs(FIGS, exist_ok=True)

PAIRS = [(0, 1), (0, 2), (1, 2)]
PAIRCOL = {(0, 1): "tab:red", (0, 2): "tab:green", (1, 2): "tab:blue"}
PAIRLAB = {(0, 1): "12", (0, 2): "13", (1, 2): "23"}


# --------------------------------------------------------------------------- model
def type1(eps, gam, a):
    """Type-1 Cauchy Hamiltonian data: H(u) = H0 + u*diag(a)."""
    eps = np.asarray(eps, float); gam = np.asarray(gam, float); a = np.asarray(a, float)
    g2 = gam ** 2
    H0 = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            if i != j:
                H0[i, j] = gam[i] * gam[j] * (a[i] - a[j]) / (eps[i] - eps[j])
        H0[i, i] = -sum(g2[k] * (a[i] - a[k]) / (eps[i] - eps[k]) for k in range(3) if k != i)
    return H0, np.diag(a)


def Gamma(eps, gam, a, i, j):
    return gam[i] ** 2 * gam[j] ** 2 * abs(a[i] - a[j]) / (eps[i] - eps[j]) ** 2


def crossing_u(eps, gam, a, i, j):
    """Diabatic crossing location u where the two diabatic energies a_i u + H0_ii cross."""
    H0, A = type1(eps, gam, a)
    return -(H0[i, i] - H0[j, j]) / (a[i] - a[j])


def width(eps, gam, a, i, j):
    """Avoided-crossing width 2|H0_ij|/|a_i-a_j| (the u-extent of the gap minimum region)."""
    H0, A = type1(eps, gam, a)
    return 2 * abs(H0[i, j]) / abs(a[i] - a[j])


def sep_width_ratio(eps, gam, a):
    """min crossing separation / max avoided-crossing width (WS-D's overlap metric)."""
    uv = np.sort([crossing_u(eps, gam, a, *p) for p in PAIRS])
    w = max(width(eps, gam, a, *p) for p in PAIRS)
    return min(uv[1] - uv[0], uv[2] - uv[1]) / w


# ----------------------------------------------------------------- turning points
def turning_points(eps, gam, a):
    """
    Roots in complex u of Disc_E(chi_H) -- the pairwise turning points.
    Uses exact rational arithmetic for the discriminant, then numeric roots.
    Returns: list of (u_complex, multiplicity_label) where label in {'simple','node'}.
    """
    u, E = sp.symbols('u E')
    epsS = [sp.nsimplify(x) for x in eps]
    gamS = [sp.nsimplify(x) for x in gam]
    aS = [sp.nsimplify(x) for x in a]
    H0 = sp.zeros(3, 3)
    for i in range(3):
        for j in range(3):
            if i != j:
                H0[i, j] = gamS[i] * gamS[j] * (aS[i] - aS[j]) / (epsS[i] - epsS[j])
        H0[i, i] = -sum(gamS[k] ** 2 * (aS[i] - aS[k]) / (epsS[i] - epsS[k])
                        for k in range(3) if k != i)
    M = H0 + u * sp.diag(*aS)
    chi = sp.Poly(sp.expand((E * sp.eye(3) - M).det()), E)
    disc = sp.Poly(sp.expand(sp.resultant(chi, chi.diff(E))), u)
    coeffs = [complex(c) for c in disc.all_coeffs()]
    roots = np.roots(coeffs)
    # classify: cluster near-equal roots (the real double root = the node)
    out = []
    used = [False] * len(roots)
    for i in range(len(roots)):
        if used[i]:
            continue
        cluster = [i]
        for j in range(i + 1, len(roots)):
            if not used[j] and abs(roots[i] - roots[j]) < 1e-5:
                cluster.append(j)
        for j in cluster:
            used[j] = True
        z = np.mean([roots[j] for j in cluster])
        out.append((z, 'node' if len(cluster) > 1 else 'simple'))
    return out


def eigs_sorted_real(u, H0, A):
    """Eigenvalues at real-ish u, returned in a deterministic order (by real part)."""
    w = np.linalg.eigvals(H0 + u * A)
    return w[np.argsort(w.real)]


def classify_turning_pair(t, eps, gam, a):
    """
    Determine which diabatic pair (i,j) collides at turning point t, by tracking
    eigenvalue labels from a real reference point where the diabatic order = sorted
    real order, along a straight path to t, matching eigenvalues by continuity.
    Returns (i,j) in diabatic-index convention (0,1,2 = low,mid,high slope ordering
    inherited from sorting eigenvalues by real part at large negative u).
    """
    H0, A = type1(eps, gam, a)
    # reference: large negative real u -> eigenvalues ~ a_i u, ordered by a_i.
    order = np.argsort(a)  # low,mid,high slope
    u0 = -50.0
    w = np.linalg.eigvals(H0 + u0 * A)
    # at u0 the eigenvalue ~ a_i u0 + H0_ii; sort eigs to match slope order
    ref = np.sort(w.real)  # ascending; since u0<0, smallest = largest a
    # build label array: idx in sorted-real -> diabatic channel
    # at u0<<0: E ~ a*u0, so ascending real <-> descending a.
    diab_of_sorted = order[::-1]  # sorted-ascending position k -> diabatic channel
    labels = list(diab_of_sorted)  # labels[k] = diabatic channel of k-th eigenvalue
    wv = np.sort(w.real) + 0j  # current eigenvalues (track these)
    wv = w[np.argsort(w.real)]
    # march along straight line u0 -> t, re-matching by nearest neighbour
    N = 1600
    path = np.linspace(u0, t, N)
    for n in range(1, N):
        wn = np.linalg.eigvals(H0 + path[n] * A)
        # match wn to wv by nearest neighbour (Hungarian-lite, 3x3)
        cost = np.abs(wv[:, None] - wn[None, :])
        match = _greedy_match(cost)
        wv = wn[match]
    # at t two eigenvalues collide: find the two closest
    d = [(abs(wv[p] - wv[q]), labels[p], labels[q]) for p in range(3) for q in range(p + 1, 3)]
    d.sort()
    i, j = sorted((d[0][1], d[0][2]))
    return (i, j)


def _greedy_match(cost):
    """Return permutation match[k]=col assigned to row k, greedy on a 3x3 cost."""
    n = cost.shape[0]
    match = [-1] * n
    taken = [False] * n
    order = np.dstack(np.unravel_index(np.argsort(cost, axis=None), cost.shape))[0]
    for r, c in order:
        if match[r] == -1 and not taken[c]:
            match[r] = c; taken[c] = True
    return match


# ------------------------------------------------ pairwise energy difference (branch)
def delta_pair(u, H0, A, pair, ref_eigs=None):
    """
    Signed energy difference E_i(u)-E_j(u) for the diabatic pair, with eigenvalues
    tracked by continuity from ref_eigs (a (labels, eigenvalues) tuple).  Returns the
    complex difference and the updated tracked eigenvalues.
    Diabatic labels are assigned at u0=-50 by slope order; we carry a persistent map.
    """
    raise NotImplementedError  # replaced by StokesTracer below


# -------------------------------------------------------------- Stokes line tracer
class EigenTracer:
    """
    Tracks the three eigenvalues of H(u) as labelled diabatic channels (0,1,2) along
    an arbitrary path, by nearest-neighbour continuity, seeded at large negative real u.
    Provides delta_ij(u) along a path = E_i - E_j with consistent branch.
    """
    def __init__(self, eps, gam, a):
        self.H0, self.A = type1(eps, gam, a)
        self.a = np.asarray(a, float)
        order = np.argsort(self.a)
        self.u_seed = -60.0
        w = np.linalg.eigvals(self.H0 + self.u_seed * self.A)
        idx = np.argsort(w.real)             # ascending real
        self.seed_eigs = w[idx]
        # at u_seed<<0, E~a*u_seed ; ascending real <-> descending a
        self.seed_labels = order[::-1]       # seed_labels[k] = diabatic channel of k-th

    def trace_path(self, path):
        """Return array (len(path),3) of eigenvalues in diabatic-channel order along path.
        path[0] need not be the seed; we first connect seed->path[0]."""
        # connect seed to path[0]
        pre = np.linspace(self.u_seed, path[0], 800)
        wv = self.seed_eigs.copy()
        lab = list(self.seed_labels)
        for n in range(1, len(pre)):
            wn = np.linalg.eigvals(self.H0 + pre[n] * self.A)
            cost = np.abs(wv[:, None] - wn[None, :])
            m = _greedy_match(cost)
            wv = wn[m]
        # now wv in same positional order as seed; reorder to diabatic channel
        chan = np.empty(3, complex)
        chan[lab] = wv  # lab[k]=channel -> chan[channel]=wv[k]
        out = np.empty((len(path), 3), complex)
        cur = wv.copy()
        # store first
        c0 = np.empty(3, complex); c0[lab] = cur; out[0] = c0
        for n in range(1, len(path)):
            wn = np.linalg.eigvals(self.H0 + path[n] * self.A)
            cost = np.abs(cur[:, None] - wn[None, :])
            m = _greedy_match(cost)
            cur = wn[m]
            cc = np.empty(3, complex); cc[lab] = cur; out[n] = cc
        return out


def phase_field(H0, A, t, GX, GY, nseg=26):
    """
    Robust Stokes-phase scalar field  phi(u) = Im int_t^u (E_i - E_j) du'  on the grid
    (GX,GY), integrated along the STRAIGHT RAY t->u, with the colliding pair (E_i,E_j)
    tracked by continuity from the turning point t.  Zero set of phi = the Stokes lines
    emanating from t (three branches at 120 deg near t).  This is the standard exact-WKB
    Stokes-line definition; the ray-integration with pair tracking handles the branch of
    sqrt(disc) cleanly in the neighbourhood of t (the only place we trust it -- far past
    another turning point the ray may cross a cut, so we MASK beyond the nearest other tp).

    VECTORIZED: all grid points are advanced TOGETHER over a common normalized parameter
    s in [0,1] (u(s)=t+s*(grid-t)); at each s we batch-diagonalize the N Hamiltonians and
    track the colliding pair by nearest-neighbour to the previous step.  ~nseg batched
    eigvals calls instead of N*nseg scalar ones.
    """
    ny, nx = GX.shape
    grid = (GX + 1j * GY).ravel()
    N = grid.size
    diff = grid - t                                  # (N,)
    a = np.diag(A).astype(float)                     # slopes
    # seed colliding pair just off t along each ray
    s0 = 1e-3
    Us = t + s0 * diff
    W = _batch_eigvals(H0, a, Us)                    # (N,3) eigenvalues
    # colliding pair at t = the two closest eigenvalues (per ray)
    gaps = np.stack([np.abs(W[:, p] - W[:, q]) for (p, q) in PAIRS], axis=1)  # (N,3)
    kmin = np.argmin(gaps, axis=1)
    pidx = np.array([PAIRS[k][0] for k in kmin])
    qidx = np.array([PAIRS[k][1] for k in kmin])
    ep = W[np.arange(N), pidx]
    eq = W[np.arange(N), qidx]
    integ = np.zeros(N, complex)
    uprev = Us.copy()
    ss = np.linspace(s0, 1.0, nseg)
    for s in ss[1:]:
        Uk = t + s * diff
        Wk = _batch_eigvals(H0, a, Uk)               # (N,3)
        # match ep,eq to the closest of the 3 new eigenvalues (independently; the two
        # tracked sheets are well separated from the third except exactly at a tp)
        dep = np.abs(Wk - ep[:, None]); newp = np.argmin(dep, axis=1)
        deq = np.abs(Wk - eq[:, None]); newq = np.argmin(deq, axis=1)
        # if both map to same sheet (near-collision), keep previous assignment for the
        # weaker match
        clash = newp == newq
        if np.any(clash):
            # for clashes, force q to the next-best distinct sheet
            for idx in np.nonzero(clash)[0]:
                order = np.argsort(np.abs(Wk[idx] - eq[idx]))
                newq[idx] = order[1] if order[0] == newp[idx] else order[0]
        ep = Wk[np.arange(N), newp]
        eq = Wk[np.arange(N), newq]
        integ += (ep - eq) * (Uk - uprev)
        uprev = Uk
    return integ.imag.reshape(ny, nx)


def _batch_eigvals(H0, a, Us):
    """Eigenvalues of H0 + u*diag(a) for an array Us of u-values. Returns (len(Us),3)."""
    Us = np.asarray(Us)
    M = H0[None, :, :].astype(complex) + Us[:, None, None] * np.diag(a)[None, :, :]
    return np.linalg.eigvals(M)


def contour_polylines(GX, GY, phi, level=0.0):
    """Extract zero-level contour polylines (lists of complex points) from a scalar field."""
    cs = plt.contour(GX, GY, phi, levels=[level])
    polylines = []
    # matplotlib >=3.8: use allsegs
    for seg in cs.allsegs[0]:
        if len(seg) >= 2:
            polylines.append(seg[:, 0] + 1j * seg[:, 1])
    plt.close('all')
    return polylines


def stokes_lines_for_tp(H0, A, t, all_tps, span, grid_n=110):
    """
    Stokes lines (zero-contours of phi) emanating from turning point t, masked to a disk
    of radius = 0.92 * (distance to the nearest OTHER turning point) so we never trust the
    ray integral past another branch point (where the cut would corrupt it).
    Returns list of complex polylines.
    """
    others = [abs(t - s) for s, _ in all_tps if abs(t - s) > 1e-6]
    rmask = 0.92 * min(others) if others else span
    rmask = min(rmask, span)
    pad = rmask * 1.15
    gx = np.linspace(t.real - pad, t.real + pad, grid_n)
    gy = np.linspace(t.imag - pad, t.imag + pad, grid_n)
    GX, GY = np.meshgrid(gx, gy)
    phi = phase_field(H0, A, t, GX, GY)
    # mask outside disk
    R = np.abs((GX - t.real) + 1j * (GY - t.imag))
    phi = np.where(R <= rmask, phi, np.nan)
    polylines = contour_polylines(GX, GY, phi)
    # keep only polyline pieces that actually touch the turning point neighbourhood
    keep = []
    for pl in polylines:
        if np.min(np.abs(pl - t)) < 0.25 * rmask:
            keep.append(pl)
    return keep, rmask


def _greedy_match_2x3(cost, ncol):
    """Match 2 rows to distinct columns among ncol, greedy."""
    n = cost.shape[0]
    match = [-1] * n
    taken = [False] * ncol
    order = np.dstack(np.unravel_index(np.argsort(cost, axis=None), cost.shape))[0]
    for r, c in order:
        if match[r] == -1 and not taken[c]:
            match[r] = c; taken[c] = True
    return match


# ----------------------------------------------------------------- joint detection
def segments_intersect(p1, p2, p3, p4):
    """Do segments p1p2 and p3p4 (complex points) intersect? Return intersection or None."""
    d1 = p2 - p1
    d2 = p4 - p3
    den = d1.real * d2.imag - d1.imag * d2.real
    if abs(den) < 1e-14:
        return None
    t = ((p3.real - p1.real) * d2.imag - (p3.imag - p1.imag) * d2.real) / den
    s = ((p3.real - p1.real) * d1.imag - (p3.imag - p1.imag) * d1.real) / den
    if -1e-9 <= t <= 1 + 1e-9 and -1e-9 <= s <= 1 + 1e-9:
        return p1 + t * d1
    return None


def _line_intersections(la, lb):
    """
    Vectorized: all intersection points of polyline la with polyline lb.
    Returns list of complex intersection points.  Uses bounding-box prefilter then
    exact segment-segment test on surviving candidate pairs.
    """
    if len(la) < 2 or len(lb) < 2:
        return []
    a0 = la[:-1]; a1 = la[1:]
    b0 = lb[:-1]; b1 = lb[1:]
    # bounding-box prefilter (broadcast over A x B segment grids)
    axmin = np.minimum(a0.real, a1.real)[:, None]; axmax = np.maximum(a0.real, a1.real)[:, None]
    aymin = np.minimum(a0.imag, a1.imag)[:, None]; aymax = np.maximum(a0.imag, a1.imag)[:, None]
    bxmin = np.minimum(b0.real, b1.real)[None, :]; bxmax = np.maximum(b0.real, b1.real)[None, :]
    bymin = np.minimum(b0.imag, b1.imag)[None, :]; bymax = np.maximum(b0.imag, b1.imag)[None, :]
    overlap = (axmin <= bxmax) & (axmax >= bxmin) & (aymin <= bymax) & (aymax >= bymin)
    ii, jj = np.nonzero(overlap)
    out = []
    for i, j in zip(ii, jj):
        x = segments_intersect(a0[i], a1[i], b0[j], b1[j])
        if x is not None:
            out.append(x)
    return out


def find_joints(lines_by_pair, turning_pts, min_dist_from_tp=0.15):
    """
    Find intersections of Stokes lines of DIFFERENT pair-types, away from any turning
    point.  Returns list of (u_intersection, pairA, pairB).
    """
    joints = []
    keys = list(lines_by_pair.keys())
    for ia in range(len(keys)):
        for ib in range(ia + 1, len(keys)):
            pa, pb = keys[ia], keys[ib]
            if pa == pb:
                continue
            for la in lines_by_pair[pa]:
                for lb in lines_by_pair[pb]:
                    for x in _line_intersections(la, lb):
                        if min(abs(x - t) for t, _ in turning_pts) > min_dist_from_tp:
                            joints.append((x, pa, pb))
    # dedupe joints that are essentially the same point
    ded = []
    for x, pa, pb in joints:
        if all(abs(x - y) > 0.05 or {pa, pb} != {qa, qb} for y, qa, qb in ded):
            ded.append((x, pa, pb))
    return ded


# -------------------------------------------------------- WKB action at a joint
def wkb_action_between(tracer, t, u_target, pair, n=600):
    """
    |Im int_t^{u_target} (E_i-E_j) du| along a straight path -- the relative WKB action.
    A joint that sits at small |Im action| relative to neighbouring turning points is
    'deep' inside the Stokes region (strong joint).
    """
    path = np.linspace(t, u_target, n)
    ev = tracer.trace_path(path)
    i, j = pair
    delta = ev[:, i] - ev[:, j]
    integ = np.trapz(delta, path) if hasattr(np, "trapz") is False else np.trapezoid(delta, path)
    return integ


# --------------------------------------------------------------- benchmark P_2->2
def benchmark_P(eps, gam, a, T=200.0, rtol=1e-12, atol=1e-13):
    H0, A = type1(eps, gam, a)
    rhs = lambda u, y: (-1j * ((H0 + u * A) @ y.reshape(3, 3))).ravel()
    sol = solve_ivp(rhs, [-T, T], np.eye(3, dtype=complex).ravel(),
                    rtol=rtol, atol=atol, method='DOP853')
    return np.abs(sol.y[:, -1].reshape(3, 3)) ** 2


def mid_enhancement(eps, gam, a):
    """P_mid / P_mid^incoherent (the 2.5-104x diagnostic), with mid = middle slope."""
    P = benchmark_P(eps, gam, a)
    lo, mid, hi = np.argsort(a)
    inc = np.exp(-2 * np.pi * (Gamma(eps, gam, a, mid, lo) + Gamma(eps, gam, a, mid, hi)))
    return P[mid, mid], inc, P[mid, mid] / inc


# ----------------------------------------------------------------------- plotting
def build_stokes_graph(eps, gam, a, span=6.0):
    """
    Build the full Stokes graph: for every turning point (simple AND the node), classify
    its colliding pair, extract its Stokes lines (zero-contours of the WKB phase), and
    detect joints (intersections of lines of DIFFERENT pair-type away from a tp).
    Returns (tps, lines_by_pair, tp_pairs, joints).
    """
    H0, A = type1(eps, gam, a)
    tps = turning_points(eps, gam, a)
    lines_by_pair = {p: [] for p in PAIRS}
    tp_pairs = []
    for t, kind in tps:
        pair = classify_turning_pair(t, eps, gam, a)
        tp_pairs.append((t, pair, kind))
        lines, _ = stokes_lines_for_tp(H0, A, t, tps, span)
        for ln in lines:
            lines_by_pair[pair].append(ln)
    joints = find_joints(lines_by_pair, tps)
    return tps, lines_by_pair, tp_pairs, joints


def plot_stokes_graph(eps, gam, a, fname, title, span=6.0):
    tps, lines_by_pair, tp_pairs, joints = build_stokes_graph(eps, gam, a, span=span)
    simple_tps = [(t, p) for t, p, k in tp_pairs if k != 'node']

    fig, ax = plt.subplots(figsize=(7.2, 6.4))
    # diabatic crossing real locations (vertical guides)
    for p in PAIRS:
        uc = crossing_u(eps, gam, a, *p)
        ax.axvline(uc, color=PAIRCOL[p], ls=':', lw=0.7, alpha=0.5)
    # Stokes lines
    drawn = set()
    for p in PAIRS:
        for ln in lines_by_pair[p]:
            lab = f"Stokes {PAIRLAB[p]}" if p not in drawn else None
            ax.plot(ln.real, ln.imag, color=PAIRCOL[p], lw=1.3, alpha=0.85, label=lab)
            drawn.add(p)
    # turning points
    for t, kind in tps:
        if kind == 'node':
            ax.plot(t.real, t.imag, 'k*', ms=14, zorder=5,
                    label='node (exact crossing)' if 'node' not in drawn else None)
            drawn.add('node')
        else:
            ax.plot(t.real, t.imag, 'ko', ms=7, mfc='white', mew=1.5, zorder=5,
                    label='turning point' if 'tp' not in drawn else None)
            drawn.add('tp')
    # joints
    for x, pa, pb in joints:
        ax.plot(x.real, x.imag, 'mD', ms=11, mfc='none', mew=2.2, zorder=6,
                label='JOINT' if 'joint' not in drawn else None)
        drawn.add('joint')
    ax.set_xlabel('Re u'); ax.set_ylabel('Im u')
    ax.set_title(title)
    ax.axhline(0, color='gray', lw=0.5, alpha=0.4)
    ax.legend(loc='upper right', fontsize=8, framealpha=0.9)
    ax.set_aspect('equal', adjustable='datalim')
    fig.tight_layout()
    fig.savefig(fname, dpi=140)
    plt.close(fig)
    return simple_tps, joints, tps


# ----------------------------------------------------------------- joint strength
def joint_strength(tps, joints):
    """
    Geometric 'depth' of a joint:  1 - dist(joint, nearest turning point)/median(tp-tp
    separation), clamped to [0,1].  A joint that sits well inside the turning-point cluster
    (small dist relative to the cluster scale) is 'strong' (deep in the Stokes region); a
    joint grazing the edge is 'weak'.  Returns (max_strength, n_joints).
    """
    tp = [t for t, _ in tps]
    if not joints:
        return 0.0, 0
    seps = [abs(tp[i] - tp[j]) for i in range(len(tp)) for j in range(i + 1, len(tp))]
    scale = np.median([s for s in seps if s > 1e-6]) if any(s > 1e-6 for s in seps) else 1.0
    strengths = []
    for x, pa, pb in joints:
        dmin = min(abs(x - t) for t in tp)
        strengths.append(max(0.0, 1.0 - dmin / scale))
    return float(np.max(strengths)), len(joints)


# ============================================================================ MAIN
def build_well_separated():
    """
    Construct a WELL-SEPARATED Type-1 sample (sep/width ratio ~3-4) by widening eps,
    using small SIGNED couplings, and modest slope gaps so the three diabatic crossings
    are far apart relative to the (narrow) avoided-crossing widths.  WS-D established the
    ratio is bounded ~O(1) with a thin tail to ~4 -- reaching ratio ~3-4 requires exactly
    this tuning (one tiny coupling, large eps separation).  Search for the largest ratio.
    """
    best = None
    rng = np.random.default_rng(3)
    for _ in range(200000):
        eps = np.sort(rng.uniform(-12, 12, 3))
        if np.min(np.diff(eps)) < 2.0:
            continue
        gam = rng.uniform(0.1, 1.0, 3) * rng.choice([-1, 1], 3)
        a = np.sort(rng.uniform(-4, 4, 3))
        if np.min(np.diff(a)) < 0.4:
            continue
        try:
            r = sep_width_ratio(eps, gam, a)
        except Exception:
            continue
        if np.isfinite(r) and (best is None or r > best[0]):
            best = (r, tuple(round(x, 4) for x in eps),
                    tuple(round(x, 4) for x in gam), tuple(round(x, 4) for x in a))
    return best


def random_sample(rng):
    """Random Type-1 sample over a space wide enough to span sep/width ratio ~0.1-4."""
    while True:
        eps = np.sort(rng.uniform(-12, 12, 3))
        if np.min(np.diff(eps)) < 0.5:
            continue
        gam = rng.uniform(0.1, 1.6, 3) * rng.choice([-1, 1], 3)
        a = np.sort(rng.uniform(-4, 4, 3))
        if np.min(np.diff(a)) < 0.4:
            continue
        try:
            r = sep_width_ratio(eps, gam, a)
        except Exception:
            continue
        if np.isfinite(r):
            return (tuple(float(x) for x in eps), tuple(float(x) for x in gam),
                    tuple(float(x) for x in a), r)


def main():
    print("=" * 78)
    print("WS-G  Stokes graph / joint diagnosis for Type-1 N=3 MLZ")
    print("=" * 78)

    # ---- the two showcase samples ----
    sampleB = ((-1., 0., 1.5), (.9, 1.1, .8), (-.7, .4, 1.3))   # OVERLAPPING (ratio ~0.5)
    rB = sep_width_ratio(*sampleB)
    print(f"\nOVERLAPPING sampleB: eps={sampleB[0]} gam={sampleB[1]} a={sampleB[2]}")
    print(f"  sep/width ratio = {rB:.3f}")

    ws = build_well_separated()
    sampleA = (tuple(float(x) for x in ws[1]),
               tuple(float(x) for x in ws[2]),
               tuple(float(x) for x in ws[3]))
    print(f"\nWELL-SEPARATED (constructed): ratio={ws[0]:.3f}")
    print(f"  eps={sampleA[0]} gam={sampleA[1]} a={sampleA[2]}")

    for lbl, smp in [("WELL-SEPARATED (sampleA)", sampleA),
                     ("OVERLAPPING (sampleB)", sampleB)]:
        r = sep_width_ratio(*smp)
        Pmid, inc, enh = mid_enhancement(*smp)
        fn = os.path.join(FIGS, f"stokes_{'sepA' if 'WELL' in lbl else 'overlapB'}.png")
        title = (f"{lbl}\nsep/width={r:.2f}  P_2->2={Pmid:.4f}  "
                 f"incoh={inc:.4f}  enh={enh:.2f}x")
        stps, joints, tps = plot_stokes_graph(*smp, fname=fn, title=title)
        strg, nj = joint_strength(tps, joints)
        print(f"\n[{lbl}]")
        print(f"  turning-point pairs: {[(round(t.real,3),round(t.imag,3),PAIRLAB[p]) for t,p in stps]}")
        print(f"  joints found: {nj}   max strength: {strg:.3f}")
        if joints:
            for x, pa, pb in joints:
                print(f"     joint @ {x.real:+.3f}{x.imag:+.3f}i  types {PAIRLAB[pa]}x{PAIRLAB[pb]}")
        print(f"  P_2->2={Pmid:.5f}  incoherent={inc:.5f}  enhancement={enh:.2f}x")
        print(f"  saved figure: {fn}")

    # --------------- correlation scan: ~16 samples, ratio 0.3..4 ----------------
    print("\n" + "=" * 78)
    print("CORRELATION SCAN: joint presence/strength vs P_2->2 enhancement")
    print("=" * 78)
    rng = np.random.default_rng(2024)
    rows = []
    # include the two showcase samples
    cand = [sampleB, sampleA]
    # gather a spread of ratios, then stratify into log-spaced ratio bins so the high-
    # ratio (well-separated) tail is represented, not just the dense ~O(0.5) bulk.
    pool = [random_sample(rng) for _ in range(3000)]
    bins = np.geomspace(0.1, 4.5, 15)
    picked = []
    for b0, b1 in zip(bins[:-1], bins[1:]):
        inbin = [s for s in pool if b0 <= s[3] < b1]
        if inbin:
            picked.append(min(inbin, key=lambda s: abs(s[3] - 0.5 * (b0 + b1))))
    for eps, gam, a, r in picked:
        cand.append((eps, gam, a))
    seen = set()
    print(f"\n{'ratio':>7} {'#joints':>7} {'strength':>8} {'P2->2':>9} {'incoh':>9} {'enh(x)':>8}")
    print("-" * 56)
    for smp in cand:
        key = tuple(round(x, 4) for t in smp for x in t)
        if key in seen:
            continue
        seen.add(key)
        try:
            r = sep_width_ratio(*smp)
            tps, lines_by_pair, tp_pairs, joints = build_stokes_graph(*smp, span=6.0)
            strg, nj = joint_strength(tps, joints)
            Pmid, inc, enh = mid_enhancement(*smp)
            rows.append((r, nj, strg, Pmid, inc, enh))
            print(f"{r:7.3f} {nj:7d} {strg:8.3f} {Pmid:9.5f} {inc:9.5f} {enh:8.2f}")
        except Exception as e:
            print(f"  (skip sample: {e})")
    rows.sort()

    # correlation summary
    rows_a = np.array([(r, nj, strg, enh) for r, nj, strg, _, _, enh in rows])
    if len(rows_a) > 3:
        rr, jj, ss, ee = rows_a.T
        print("\nCORRELATIONS (Pearson):")
        print(f"  ratio vs enhancement      : {np.corrcoef(rr, ee)[0,1]:+.3f}")
        print(f"  #joints vs enhancement    : {np.corrcoef(jj, ee)[0,1]:+.3f}")
        print(f"  joint-strength vs enhancement: {np.corrcoef(ss, ee)[0,1]:+.3f}")
        print(f"  ratio vs #joints          : {np.corrcoef(rr, jj)[0,1]:+.3f}")
        print(f"  ratio vs joint-strength   : {np.corrcoef(rr, ss)[0,1]:+.3f}")
        # the 15/85 split: incoherent exact <=> enh ~ 1
        exact = ee < 1.15
        print(f"\n  'incoherent exact' (enh<1.15): {np.sum(exact)}/{len(ee)} samples")
        if np.any(exact) and np.any(~exact):
            print(f"     mean #joints  exact={jj[exact].mean():.2f}  overlap={jj[~exact].mean():.2f}")
            print(f"     mean strength exact={ss[exact].mean():.3f}  overlap={ss[~exact].mean():.3f}")
    print("\nDone. Figures in", FIGS)


if __name__ == "__main__":
    main()
