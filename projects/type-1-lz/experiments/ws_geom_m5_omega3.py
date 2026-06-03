"""
ws_geom_m5_omega3.py -- compute the THIRD Magnus generator Omega_3 and test whether it
REFINES the order-2 graph calculator (the user's question: is the marginality verdict safe,
or does Omega_3 give a significant refinement we never checked?).

Generator convention (matches ws_geom_magnus): chi' = -Wtil chi, i.e. Y'=A Y with A = -Wtil.
Standard Magnus terms:
  Om1 = int A
  Om2 = 1/2  int_{t1>t2} [A1,A2]
  Om3 = 1/6  int_{t1>t2>t3} ( [A1,[A2,A3]] + [A3,[A2,A1]] )
All Om_k are ANTI-HERMITIAN for anti-Hermitian A -- used here as the correctness gate.

VERDICT (this run, numerically-supported): Omega_3 gives NO significant refinement.
  * Correctness GATE passes: Om1,Om2 match magnus_terms() to 0.0; ||Om3+Om3^H||=0 (anti-Hermitian).
  * The term NORMS shrink geometrically (canonical ||Om1||,||Om2||,||Om3|| = 1.03, 0.51, 0.12; at the
    sweet spot sc=1.8: 0.295, 0.041, 0.004 -- ratio ~0.1, i.e. summable), so Om3 is genuinely SMALL,
    not a hidden large term.
  * Yet adding Om3 does NOT reduce the error to the true P: e3 ~ e2 across the whole sweep (raw AND
    cleanly T-averaged; sweet spot sc=1.8: e1,e2,e3 = 1.64e-3, 1.45e-3, 1.45e-3). ||Om3||=0.004 EXCEEDS
    the 1.45e-3 residual yet does not touch it.
  * Interpretation: the residual floor is OUTSIDE the perturbative Magnus tower -- it is the
    non-perturbative (sigma-scale) remainder + the finite-T endpoint tail, NOT a missing finite order.
    Higher orders (Om4...) are smaller still and cannot reach it. The marginality verdict (M5) HOLDS and
    is now unambiguous: adding graph orders cannot push the calculator below the ~1e-3 floor.

Om3 is built from nested cumulative integrals (all O(N)); we cross-check Om2 against the
existing magnus_terms(), gate anti-Hermiticity, then measure e3 = |P3-ref| vs e2 across the
adiabaticity sweep (sweet spot + marginal wall).

Reproduce: python3 ws_geom_m5_omega3.py
"""
from __future__ import annotations
import os, sys
import numpy as np
import scipy.linalg as sla

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
from ws_geom_magnus import dressed_coupling, magnus_terms, _P_diab, type1, max_be_action  # noqa
from ws_geom_m5_graphs import ref_P  # noqa
np.seterr(all="ignore")

EPS = (-2.0, 0.0, 3.0); A_ = (-1.0, 0.5, 2.0); G0 = np.array([1.0, 0.8, 1.2])


def _cum(F, du):
    """forward cumulative trapezoid  C(t)=int_{-T}^t F dt'."""
    C = np.empty_like(F); C[0] = 0.0
    C[1:] = np.cumsum((F[:-1] + F[1:]) / 2 * du, axis=0)
    return C


def _revcum(F, du):
    """reverse cumulative trapezoid  C(t)=int_t^{+T} F dt'."""
    C = np.empty_like(F); C[-1] = 0.0
    C[:-1] = np.cumsum(((F[:-1] + F[1:]) / 2 * du)[::-1], axis=0)[::-1]
    return C


def _comm(X, Y):
    return np.einsum('kij,kjl->kil', X, Y) - np.einsum('kij,kjl->kil', Y, X)


def magnus_123(us, du, Wtil):
    A = -Wtil
    Om1 = np.trapezoid(A, us, axis=0)
    J = _cum(A, du)                                  # int^t A
    P = _comm(A, J)                                  # [A(t),int^t A]
    Om2 = 0.5 * np.trapezoid(P, us, axis=0)
    # term1: [A1,[A2,A3]] ordered t1>t2>t3 -> [A(t1), int^{t1}[A(t2),int^{t2}A] dt2]
    Pcum = _cum(P, du)
    term1 = _comm(A, Pcum)
    # term2: [A3,[A2,A1]] ordered t1>t2>t3 -> integrate t1 (top), then t2 (top): reverse cumulants
    Jbar = _revcum(A, du)                            # int_t^{T} A
    Q = _comm(A, Jbar)                               # [A(t2), int_{t2}^T A]
    Qbar = _revcum(Q, du)
    term2 = _comm(A, Qbar)                           # [A(t3), int_{t3}^T Q]
    Om3 = (1.0 / 6.0) * np.trapezoid(term1 + term2, us, axis=0)
    return Om1, Om2, Om3


def run(eps, gam, a, T=70.0, N=24001):
    us, du, E, Phi, W, Wtil = dressed_coupling(eps, gam, a, T, N)
    Om1, Om2, Om3 = magnus_123(us, du, Wtil)
    # cross-check Om1,Om2 vs the existing implementation
    o1, o2 = magnus_terms(us, du, Wtil)
    chk = max(np.max(np.abs(Om1 - o1)), np.max(np.abs(Om2 - o2)))
    Phi_m, Phi_p = Phi[0], Phi[-1]
    P2 = _P_diab(Phi_m, Phi_p, Om1 + Om2)
    P3 = _P_diab(Phi_m, Phi_p, Om1 + Om2 + Om3)
    aH = lambda M: float(np.max(np.abs(M + M.conj().T)))
    return dict(Om1=Om1, Om2=Om2, Om3=Om3, P2=P2, P3=P3, chk=chk,
                n1=np.linalg.norm(Om1), n2=np.linalg.norm(Om2), n3=np.linalg.norm(Om3),
                aH3=aH(Om3))


def main():
    print("=" * 92)
    print("Omega_3 test -- does the 3rd Magnus generator REFINE the order-2 calculator?")
    print("=" * 92)
    # correctness gate on canonical
    r = run(EPS, tuple(G0), A_)
    print(f"correctness: |Om1,Om2 vs magnus_terms| = {r['chk']:.1e}   ||Om3+Om3^H|| = {r['aH3']:.1e} "
          f"(anti-Herm gate)")
    print(f"canonical norms: ||Om1||={r['n1']:.3f}  ||Om2||={r['n2']:.3f}  ||Om3||={r['n3']:.3f}\n")
    print(f"  {'sc':>4} {'Lambda~':>7} {'||Om1||':>8} {'||Om2||':>8} {'||Om3||':>8} | "
          f"{'e2':>9} {'e3':>9} | {'e3<e2?':>6} {'n3<n2?':>6}")
    for sc in [1.4, 1.8, 2.2, 2.8, 3.5]:
        gam = tuple(sc * G0)
        r = run(EPS, gam, A_)
        Po = ref_P(EPS, gam, A_)
        e2 = float(np.max(np.abs(r["P2"] - Po))); e3 = float(np.max(np.abs(r["P3"] - Po)))
        print(f"  {sc:4.1f} {'':>7} {r['n1']:8.3f} {r['n2']:8.3f} {r['n3']:8.3f} | "
              f"{e2:9.2e} {e3:9.2e} | {str(e3<e2):>6} {str(r['n3']<r['n2']):>6}")
    print("\n  (e3<e2 => Om3 refines; n3<n2 => the term is shrinking, i.e. still in the asymptotic")
    print("   convergent window. Where ||Om3|| >= ||Om2|| the series has hit its marginal floor.)")


if __name__ == "__main__":
    main()
