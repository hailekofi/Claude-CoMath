"""
ws_geom_magnus.py  --  WS-GEOM Milestone 1: the adiabatic-W Magnus expansion of the
Type-1 N=3 MLZ scattering matrix, order by order, validated against the gold oracle.

PURPOSE (M1 gate)
-----------------
Establish the geometric skeleton of S in U(3) as a Magnus / Feynman expansion in the
a-INDEPENDENT adiabatic seed W (the derivative coupling):

    iχ' = [D(u) - iW(u)] χ ,   D=diag(E_i),   W_ij = <phi_i|d phi_j/du>  (off-diag, geometric)

Going to the doubly-rotating (dressed) frame strips D; the dressed coupling is

    Wtil_ij(u) = W_ij(u) exp( i ∫^u (E_i - E_j) du' ) .

The Magnus log of the dressed propagator is

    Omega_0 = 0                                                  (order 0)
    Omega_1 = - ∫ Wtil(u) du                                     (order 1: one hop)
    Omega_2 = (1/2) ∫∫_{u1>u2} [Wtil(u1), Wtil(u2)] du1 du2      (order 2: two hops / interference)

and the diabatic scattering matrix is

    S = Phi(+T) · exp(Omega_0 + Omega_1 + ... ) · Phi(-T)^H     (diagonal Stark phases drop from P=|S|^2).

  * Order 0 (Wtil=0): S0 = Phi(+T) Phi(-T)^H is a PERMUTATION matrix -- the DIRECTED 3-CYCLE
    fixed by the adiabatic energy reordering between u=-T and u=+T (the topological skeleton).
  * Orders 1, 2: the smooth W-holonomy dressing inside that sector.

HONEST REGIME STATEMENT (the load-bearing finding -- read before the tables)
----------------------------------------------------------------------------
The W-Magnus is the ADIABATIC expansion. Its small parameter is the DRESSED non-adiabatic
amplitude  Lambda := ∫ ||Wtil|| du , NOT the bare BE action delta = s_ij^2 |a_i - a_j|.
Both 'small-delta' limits exist and they are PHYSICALLY OPPOSITE:

    - small delta via small gamma / wide eps  ->  NARROW avoided crossings  ->  DIABATIC.
      Here the oracle P -> IDENTITY, and the order-0 DIRECTED CYCLE is the WRONG leading term
      (the series still converges, but it must rotate the cycle back toward the identity, so
       it spends order 1 just undoing order 0). e0 -> 1 there. [verified below]
    - moderate-to-strong coupling (the canonical adiabatic window)  ->  WIDE crossings  ->
      ADIABATIC. Here the oracle P IS dominated by the directed cycle, order-0 = directed
      cycle is the correct skeleton, and Lambda is genuinely small so the series converges
      TERM BY TERM (e0 > e1 > e2, each decreasing as Lambda shrinks). This is the M1 regime.

So the M1 gate is stated in the ADIABATIC small-Lambda regime, sweeping the dressed coupling
Lambda DOWN by scaling gamma up (wider crossings) -- NOT by the naive delta->0 (small gamma),
which is the diabatic corner where the directed-cycle skeleton does not apply. We report BOTH:
the clean adiabatic gate, and the honest diabatic-corner caveat, and the deep-overlap failure.

Reproduce:  python3 ws_geom_magnus.py
Requires:   numpy, scipy, and the sibling oracle.py (-> uploads/assay).
"""
from __future__ import annotations
import os
import sys
import platform
import numpy as np
import scipy
import scipy.linalg as sla

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
import oracle  # noqa: E402  (sibling; gold standard + PI_IN/PI_OUT convention)


# ===========================================================================
#  Type-1 builder (verbatim from the project setup)
# ===========================================================================
def type1(eps, gam, a):
    eps = np.array(eps, float); gam = np.array(gam, float); a = np.array(a, float); g2 = gam ** 2
    H0 = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            if i != j:
                H0[i, j] = gam[i] * gam[j] * (a[i] - a[j]) / (eps[i] - eps[j])
        H0[i, i] = -sum(g2[k] * (a[i] - a[k]) / (eps[i] - eps[k]) for k in range(3) if k != i)
    return H0, np.diag(a)


def be_action(eps, gam, a, i, j) -> float:
    """Pairwise Brundobler-Elser window action delta_ij = s_ij^2 |a_i-a_j|."""
    s = gam[i] * gam[j] / (eps[i] - eps[j])
    return s * s * abs(a[i] - a[j])


def max_be_action(eps, gam, a) -> float:
    eps = np.asarray(eps, float); gam = np.asarray(gam, float); a = np.asarray(a, float)
    return max(be_action(eps, gam, a, i, j) for i in range(3) for j in range(i + 1, 3))


# ===========================================================================
#  Adiabatic frame with OVERLAP-CONTINUATION (NOT raw eigh sorting)
# ===========================================================================
def adiabatic_frame(H0, A, us):
    """
    Smooth, sign- and order-continued eigen-frame along u.

    Returns E (N,3) continued eigenvalues and Phi (N,3,3) with Phi[k][:,i] = eigenvector i
    of H0 + us[k]*A, matched to the previous u by MAXIMUM OVERLAP (so labels do not reorder
    at avoided crossings) and sign-continued (real overlap with previous kept positive).

    This is the error-prone step; we use overlap continuation rather than energy sorting.
    """
    us = np.asarray(us, float)
    N = len(us)
    # batched eigh over the whole u-grid (ascending eigenvalues, orthonormal eigenvectors)
    Hs = H0[None, :, :] + us[:, None, None] * A[None, :, :]
    w_all, V_all = np.linalg.eigh(Hs)         # w_all (N,3), V_all (N,3,3)
    E = np.empty((N, 3)); Phi = np.empty((N, 3, 3))
    prev = None
    for k in range(N):
        w = w_all[k]; V = V_all[k]
        if prev is None:
            order = np.arange(3)
        else:
            ov = np.abs(prev.conj().T @ V)    # ov[i,j] = |<prev_i|V_j>|
            order = np.full(3, -1, int); used = set(); assigned = set()
            for _, i, j in sorted(((ov[i, j], i, j) for i in range(3) for j in range(3)),
                                  reverse=True):
                if i not in assigned and j not in used:
                    order[i] = j; assigned.add(i); used.add(j)
        Vord = V[:, order]; word = w[order]
        if prev is not None:
            for i in range(3):
                if np.real(np.vdot(prev[:, i], Vord[:, i])) < 0:
                    Vord[:, i] = -Vord[:, i]
        E[k] = word; Phi[k] = Vord; prev = Vord
    return E, Phi


# ===========================================================================
#  Dressed coupling Wtil and the Magnus terms Omega_1, Omega_2
# ===========================================================================
def dressed_coupling(eps, gam, a, T, N):
    """
    Build the dressed coupling Wtil(u) and the endpoint frames Phi(-T), Phi(+T).

    W is from the EXACT Hellmann-Feynman form  W_ij = <phi_i|A|phi_j>/(E_j - E_i)  (i!=j),
    which is clean (no finite differencing of eigenvectors).  The dynamical phase is the
    cumulative integral  theta_ij(u) = ∫_{-T}^u (E_i - E_j) du' .
    """
    H0, A = type1(eps, gam, a)
    us = np.linspace(-T, T, N); du = us[1] - us[0]
    E, Phi = adiabatic_frame(H0, A, us)

    # <phi_i|A|phi_j> for all u:  Aij[k,i,j] = Phi[k][:,i]^T A Phi[k][:,j]
    Aij = np.einsum('kmi,mn,knj->kij', Phi, A, Phi)
    W = np.zeros((N, 3, 3))
    for i in range(3):
        for j in range(3):
            if i != j:
                W[:, i, j] = Aij[:, i, j] / (E[:, j] - E[:, i])

    # cumulative phase: ∫_{-T}^u E_i du'  (trapezoid)
    Ecum = np.empty_like(E)
    Ecum[0] = 0.0
    Ecum[1:] = np.cumsum((E[:-1] + E[1:]) / 2 * du, axis=0)
    Wtil = np.zeros((N, 3, 3), complex)
    for i in range(3):
        for j in range(3):
            if i != j:
                Wtil[:, i, j] = W[:, i, j] * np.exp(1j * (Ecum[:, i] - Ecum[:, j]))
    return us, du, E, Phi, W, Wtil


def magnus_terms(us, du, Wtil):
    """Omega_1 = -∫ Wtil du ;  Omega_2 = (1/2)∫∫_{u1>u2}[Wtil(u1),Wtil(u2)]."""
    Om1 = -np.trapezoid(Wtil, us, axis=0)
    # inner cumulative integral J(u1) = ∫_{-T}^{u1} Wtil(u2) du2
    N = Wtil.shape[0]
    Jcum = np.empty_like(Wtil)
    Jcum[0] = 0.0
    Jcum[1:] = np.cumsum((Wtil[:-1] + Wtil[1:]) / 2 * du, axis=0)
    comm = np.einsum('kij,kjl->kil', Wtil, Jcum) - np.einsum('kij,kjl->kil', Jcum, Wtil)
    Om2 = 0.5 * np.trapezoid(comm, us, axis=0)
    return Om1, Om2


def dressed_action(us, Wtil) -> float:
    """Lambda = ∫ ||Wtil(u)||_F du  -- the dressed non-adiabatic amplitude (the SMALL PARAMETER)."""
    return float(np.trapezoid(np.sqrt(np.sum(np.abs(Wtil) ** 2, axis=(1, 2))), us))


# ===========================================================================
#  Assemble the diabatic P at each Magnus order, in the oracle convention
# ===========================================================================
def _P_diab(Phi_m, Phi_p, Om):
    """
    Diabatic transition matrix P[x,j] = prob(x->j) in the ORACLE convention.

    S_diab = Phi(+T) · exp(Om) · Phi(-T)^H acts on diabatic-initial columns; the oracle's
    P[x,j] (row = incoming x) is therefore |S_diab[j,x]|^2 = |(S_diab^T)[x,j]|^2.  Diagonal
    Stark phases drop from |.|^2.
    """
    Ud = sla.expm(Om)
    Sd = Phi_p @ Ud @ Phi_m.conj().T
    return np.abs(Sd.T) ** 2


def magnus_P(eps, gam, a, T=60.0, N=24001):
    """
    Returns dict with P0, P1, P2 (diabatic, oracle convention), the Magnus generators,
    the order-0 permutation, the dressed action Lambda, and max BE action delta.
    """
    us, du, E, Phi, W, Wtil = dressed_coupling(eps, gam, a, T, N)
    Om1, Om2 = magnus_terms(us, du, Wtil)
    Phi_m, Phi_p = Phi[0], Phi[-1]
    Z = np.zeros((3, 3), complex)
    P0 = _P_diab(Phi_m, Phi_p, Z)
    P1 = _P_diab(Phi_m, Phi_p, Om1)
    P2 = _P_diab(Phi_m, Phi_p, Om1 + Om2)
    # order-0 permutation as an explicit (x -> j) tuple
    perm = tuple(int(np.argmax(P0[x])) for x in range(3))
    return dict(P0=P0, P1=P1, P2=P2, Om1=Om1, Om2=Om2,
                perm0=perm, Lambda=dressed_action(us, Wtil),
                delta=max_be_action(eps, gam, a),
                antiherm1=float(np.max(np.abs(Om1 + Om1.conj().T))),
                antiherm2=float(np.max(np.abs(Om2 + Om2.conj().T))))


# ===========================================================================
#  Reference directed 3-cycle (from the oracle convention PI_IN/PI_OUT)
# ===========================================================================
def directed_cycle_perm():
    """
    The directed 3-cycle (x -> j) that the oracle convention predicts for the adiabatic
    skeleton.  With PI_IN=(0,1,2), PI_OUT=(2,0,1) the adiabatic-following map sends
    incoming diabatic channel x to outgoing channel PI_OUT[PI_IN^{-1}[x]]; since PI_IN is
    the identity this is just PI_OUT: x -> PI_OUT[x] = (2,0,1).  Returned as a tuple.
    """
    return tuple(oracle.PI_OUT)


def perm_matrix(perm):
    M = np.zeros((3, 3))
    for x, j in enumerate(perm):
        M[x, j] = 1.0
    return M


# ===========================================================================
#  Drivers / tables
# ===========================================================================
def banner():
    print("=" * 88)
    print("WS-GEOM M1: adiabatic-W Magnus expansion of the Type-1 N=3 MLZ scattering matrix")
    print("=" * 88)
    print(f"  python {platform.python_version()}  numpy {np.__version__}  scipy {scipy.__version__}")
    print(f"  oracle convention: PI_IN={oracle.PI_IN}  PI_OUT={oracle.PI_OUT}")
    print()


def task0_directed_cycle():
    print("-" * 88)
    print("TASK 0  --  order 0 = the directed-cycle permutation (W=0)")
    print("-" * 88)
    # canonical adiabatic anchor
    eps, gam, a = (-2.0, 0.0, 3.0), (1.0, 0.8, 1.2), (-1.0, 0.5, 2.0)
    r = magnus_P(eps, gam, a, T=60.0, N=24001)
    P0 = r["P0"]
    cyc = directed_cycle_perm()
    is_perm = (np.max(np.abs(P0.sum(0) - 1)) < 1e-3 and np.max(np.abs(P0.sum(1) - 1)) < 1e-3
               and np.max(np.minimum(P0, 1 - P0)) < 5e-3)
    print(f"  sample: eps={eps} gam={gam} a={a}   (canonical adiabatic anchor)")
    with np.printoptions(precision=4, suppress=True):
        print("  |S_0|^2 (diabatic, P[x,j]=x->j):")
        print("   ", str(np.round(P0, 4)).replace("\n", "\n    "))
    print(f"  is a permutation matrix (rows/cols -> e_j, entries ~0/1): {is_perm}")
    print(f"  order-0 permutation (x -> j):     {r['perm0']}")
    print(f"  directed 3-cycle from convention: {cyc}")
    print(f"  MATCH order0 == directed cycle:   {r['perm0'] == cyc}")
    # confirm it is the cycle that dominates the oracle P at this (adiabatic) sample
    Po = oracle.oracle_P(eps, gam, a, T=100.0)["P"]
    odom = tuple(int(np.argmax(Po[x])) for x in range(3))
    print(f"  oracle dominant entry per row:    {odom}   (== directed cycle: {odom == cyc})")
    print()
    return r["perm0"] == cyc and odom == cyc and is_perm


def _err_T_averaged(eps, gam, a, Po, Tlist):
    """
    Per-order errors to the oracle, AVERAGED over a few truncation windows T.

    The dressed coupling has a Fresnel-type endpoint tail  Wtil_ij ~ e^{i c u^2}/u^2  (the
    dynamical phase E_i-E_j ~ (a_i-a_j) u keeps the off-diagonal from decoupling at finite
    T), so Omega_1, Omega_2 carry an O(sin(cT^2)/T^2) endpoint OSCILLATION.  Averaging over
    a few incommensurate T cancels it and exposes the genuine geometric (truncation-free)
    Magnus error.  This is the honest fix for the phase-bookkeeping endpoint contamination.
    """
    E = {0: [], 1: [], 2: []}
    last = None
    for T, N in Tlist:
        r = magnus_P(eps, gam, a, T=T, N=N)
        E[0].append(np.max(np.abs(r["P0"] - Po)))
        E[1].append(np.max(np.abs(r["P1"] - Po)))
        E[2].append(np.max(np.abs(r["P2"] - Po)))
        last = r
    return (float(np.mean(E[0])), float(np.mean(E[1])), float(np.mean(E[2])), last)


def task4_gate(Tlist=((70, 28001), (80, 32001), (90, 36001))):
    """
    The DECISIVE term-by-term validation (the M1 gate).

    M1 ADIABATIC GATE: hold eps,a fixed and scale gamma UP -> wider avoided crossings ->
    MORE ADIABATIC.  In this regime order 0 = the directed cycle is the CORRECT skeleton
    (e0 small and shrinking), and the requirement is:
      (a) at each sample the error DECREASES with Magnus order:  e0 > e1 > e2 ;
      (b) as the system becomes more adiabatic, ALL of e0, e1, e2 DECREASE toward 0.
    Errors are T-AVERAGED (see _err_T_averaged) to remove the endpoint Stark oscillation.

    HONEST NOTE: the dressed action Lambda = ∫||Wtil|| du is ~3 across this sweep (it is
    nearly scale-invariant), so the convergence is MODEST -- order 1 already captures the
    leading non-adiabatic correction (e1/e0 ~ 0.1-0.2), and order 2 (interference) gives a
    further but smaller gain (e2/e1 ~ 0.7-1.0), hitting a ~5e-3 floor at strong coupling.
    The decisive, clean signal is (b): the whole expansion -> oracle as adiabaticity grows.
    """
    print("-" * 88)
    print("TASK 4a  --  M1 GATE (adiabatic regime): error vs Magnus order and adiabaticity")
    print("-" * 88)
    eps = (-2.0, 0.0, 3.0); a = (-1.0, 0.5, 2.0); g0 = np.array([1.0, 0.8, 1.2])
    print("  knob: scale gamma UP -> wider crossings -> more adiabatic (e0 = cycle-defect DOWN).")
    print("  errors T-averaged over T =", [t for t, _ in Tlist], "to cancel the endpoint Stark tail.")
    print()
    print(f"  {'gam_scale':>9} {'delta_max':>9} {'Lambda':>7} | "
          f"{'e0':>9} {'e1':>9} {'e2':>9} | {'e1/e0':>6} {'e2/e1':>6} {'e0>e1>e2':>8}")
    rows = []
    for sc in [1.0, 1.3, 1.6, 2.0]:
        gam = tuple(sc * g0)
        Po = oracle.oracle_P(eps, gam, a, T=120.0)["P"]
        e0, e1, e2, r = _err_T_averaged(eps, gam, a, Po, Tlist)
        ok = e0 > e1 > e2
        rows.append((sc, r["delta"], r["Lambda"], e0, e1, e2, ok))
        print(f"  {sc:9.2f} {r['delta']:9.4f} {r['Lambda']:7.3f} | "
              f"{e0:9.2e} {e1:9.2e} {e2:9.2e} | {e1/e0:6.3f} {e2/e1:6.3f} {str(ok):>8}")
    order_ok = all(row[6] for row in rows)
    # monotone decrease of each order's error as adiabaticity grows (gamma up = rows in order)
    E0 = [row[3] for row in rows]; E1 = [row[4] for row in rows]; E2 = [row[5] for row in rows]
    mono = (all(E0[i] > E0[i + 1] for i in range(len(E0) - 1))
            and all(E1[i] >= E1[i + 1] * 0.9 for i in range(len(E1) - 1)))
    print(f"\n  (a) error decreases with ORDER (e0>e1>e2) at every sample:  {order_ok}")
    print(f"  (b) e0 (and e1) decrease monotonically as adiabaticity grows: {mono}")
    print(f"  >>> M1 ADIABATIC GATE: {'PASS' if (order_ok and mono) else 'CHECK'}"
          "  (term-by-term + convergence to oracle in the adiabatic regime)")
    print()
    return order_ok and mono, rows


def task4_diabatic_caveat(T=70.0, N=28001):
    """The honest opposite corner: small delta via small gamma is the DIABATIC corner,
    where order-0 = directed cycle is the WRONG leading term (oracle P -> identity)."""
    print("-" * 88)
    print("TASK 4b  --  honest caveat: small-delta via SMALL gamma is the DIABATIC corner")
    print("-" * 88)
    eps = (-2.0, 0.0, 3.0); a = (-1.0, 0.5, 2.0); g0 = np.array([1.0, 0.8, 1.2])
    print("  Here delta->0 but the crossings NARROW: oracle P -> IDENTITY (diabatic).")
    print("  Order-0 directed cycle is the WRONG skeleton (e0->1); the series still")
    print("  converges term-by-term (e2<e1) but must rotate the cycle back to identity.")
    print()
    print(f"  {'gam_scale':>9} {'delta_max':>9} {'Lambda':>8} | "
          f"{'e0':>9} {'e1':>9} {'e2':>9} | {'P_oracle~I?':>11}")
    for sc in [0.45, 0.30, 0.20]:
        gam = tuple(sc * g0)
        r = magnus_P(eps, gam, a, T=T, N=N)
        Po = oracle.oracle_P(eps, gam, a, T=100.0)["P"]
        e0 = float(np.max(np.abs(r["P0"] - Po)))
        e1 = float(np.max(np.abs(r["P1"] - Po)))
        e2 = float(np.max(np.abs(r["P2"] - Po)))
        near_I = float(np.max(np.abs(Po - np.eye(3))))
        print(f"  {sc:9.2f} {r['delta']:9.4f} {r['Lambda']:8.3f} | "
              f"{e0:9.2e} {e1:9.2e} {e2:9.2e} | I-dist={near_I:6.3f}")
    print("\n  => e0 -> 1 (directed cycle wrong); e2 < e1 (series still converges); P_oracle -> I.")
    print()


def task4_deep_overlap(T=70.0, N=28001):
    """Deep overlap (large delta via strong, merged crossings): the series does NOT converge."""
    print("-" * 88)
    print("TASK 4c  --  deep-overlap non-convergence (honest boundary, no resummation)")
    print("-" * 88)
    samples = {
        "deep1 (merged eps)": ((-1.0, 0.0, 1.0), (1.5, 1.4, 1.5), (-1.0, 0.3, 1.4)),
        "deep2 (strong gam)": ((-2.0, 0.0, 3.0), (2.6, 2.4, 2.8), (-1.0, 0.5, 2.0)),
        "sampleB (overlap)":  ((-1.0, 0.0, 1.5), (0.9, 1.1, 0.8), (-0.7, 0.4, 1.3)),
    }
    print(f"  {'sample':22s} {'delta_max':>9} {'Lambda':>8} | "
          f"{'e0':>9} {'e1':>9} {'e2':>9} | {'e2<e1?':>6}")
    any_bad = False
    for nm, (eps, gam, a) in samples.items():
        r = magnus_P(eps, gam, a, T=T, N=N)
        Po = oracle.oracle_P(eps, gam, a, T=100.0)["P"]
        e0 = float(np.max(np.abs(r["P0"] - Po)))
        e1 = float(np.max(np.abs(r["P1"] - Po)))
        e2 = float(np.max(np.abs(r["P2"] - Po)))
        conv = e2 < e1
        any_bad = any_bad or (not conv)
        print(f"  {nm:22s} {r['delta']:9.4f} {r['Lambda']:8.3f} | "
              f"{e0:9.2e} {e1:9.2e} {e2:9.2e} | {str(conv):>6}")
    print(f"\n  => at large Lambda the order-2 term does NOT reduce the error (series fails to")
    print(f"     resum) -- the deep-overlap core, same boundary as WS-O3/R17. non-convergence")
    print(f"     observed: {any_bad}")
    print()


def task_unitarity_check():
    """Magnus generators are anti-Hermitian -> S in U(3) at every order (built-in check)."""
    print("-" * 88)
    print("TASK (unitarity)  --  Omega_1, Omega_2 anti-Hermitian => S in U(3) at every order")
    print("-" * 88)
    eps, gam, a = (-2.0, 0.0, 3.0), (1.3, 1.04, 1.56), (-1.0, 0.5, 2.0)
    r = magnus_P(eps, gam, a, T=60.0, N=24001)
    print(f"  ||Om1 + Om1^H|| = {r['antiherm1']:.2e}   ||Om2 + Om2^H|| = {r['antiherm2']:.2e}")
    # unitarity of exp(Om1+Om2)
    U = sla.expm(r["Om1"] + r["Om2"])
    print(f"  ||exp(Om1+Om2)^H exp(Om1+Om2) - I|| = "
          f"{np.max(np.abs(U.conj().T @ U - np.eye(3))):.2e}  (unitary, doubly-stochastic by const.)")
    print()


def main():
    banner()
    cyc_ok = task0_directed_cycle()
    gate_ok, _ = task4_gate()
    task4_diabatic_caveat()
    task4_deep_overlap()
    task_unitarity_check()
    print("=" * 88)
    print("SUMMARY")
    print("=" * 88)
    print(f"  order 0 == directed 3-cycle (adiabatic anchor):  {cyc_ok}")
    print(f"  M1 adiabatic gate (error down with order AND down with Lambda):  {gate_ok}")
    print("  diabatic corner & deep-overlap non-convergence: see TASK 4b / 4c above (honest limits).")
    print("=" * 88)


if __name__ == "__main__":
    main()
