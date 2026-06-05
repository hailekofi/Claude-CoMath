"""
ws_tanh_demkov_kunike.py -- the TANH (Demkov-Kunike) TAME analog of the LZ sweep.

Bounding the ramp turns the WILD (rank-2 irregular) LZ connection into a TAME (regular-singular,
hypergeometric-class) one -- demonstrating concretely that LZ's transcendental sigma is a consequence of
the UNBOUNDED linear ramp (the irregular singularity at u=inf), not of the level structure.

Models (same canonical gamma,eps,a; H0 the Type-1 Cauchy matrix):
   LINEAR (LZ, WILD):  H(u) = H0 + u*A          -- irregular rank-2 at u=inf; eigenvalues diverge ~ a_i u
   TANH  (DK, TAME):   H(u) = H0 + tanh(u)*A     -- regular singular at u=+-inf; eigenvalues SATURATE

Demonstrations:
 (A) SINGULARITY TYPE: tanh eigenvalues saturate (regular) vs linear diverge (irregular); the dynamical
     phase theta_ij = int(E_i-E_j)du scales ~ T^1 (regular: e^{i lambda u}) for tanh vs ~ T^2 (Fresnel:
     e^{i lambda u^2}) for linear -- the exact regular-vs-irregular signature.
 (B) CONDITIONING: tanh P converges EXPONENTIALLY in T (off-diagonal coupling decays ~ sech^2 u); the
     linear LZ P carries the slow oscillatory Fresnel endpoint tail (the Stokes/irregular wall).
 (C) CLOSED-FORM ANCHOR (N=2, rigorous): the sech pulse (Rosen-Zener; same hypergeometric family) has the
     EXACT elementary amplitude P = sin^2(Omega*pi/2). Verified to machine precision -> tame => closed form.
 (D) CONTRAST (N=3): P_mm^tanh is tame & well-conditioned; P_mm^linear is the transcendental sigma (~0.2147).

FINDINGS (numerically-supported, decisive):
  (A) theta_tanh ~ T^1.05 (regular, e^{i*lambda*u}); theta_linear ~ T^1.99 (Fresnel, e^{i*lambda*u^2},
      irregular). tanh eigenvalues saturate to eig(H0+-A) by u~5; linear diverge ~a_i*u.
  (B) P_mm^tanh converges EXPONENTIALLY to ~1e-12 (increments -2.7e-5, -1.1e-7, +2.2e-8, +8e-12 at
      T=6,8,12,16) -- the regular-singular hallmark (non-adiabatic coupling decays ~sech^2). P_mm^linear
      only oscillates with a slow ~1e-4 Fresnel tail (needs T-averaging) -- the wild/irregular hallmark.
  (C) sech/Rosen-Zener P = sin^2(Omega*pi/2) reproduced to 1e-12 -> tame => exact elementary closed form.
  (D) P_mm^tanh = 0.20383 (tame; same H0,A, just a BOUNDED ramp) vs P_mm^linear = 0.21472 (= sigma, wild).
  RIGIDITY CAVEAT (AD): tame != closed-form. In z=(1+tanh u)/2 the tanh model is a rank-N Fuchsian system
  with 3 regular singular points; Katz rigidity rig=(2-3)N^2+3N=N(3-N): N=2 => rig=2 (RIGID => 2F1, Gamma-
  ratios, elementary, the (C) anchor); N=3 => rig=0 (NON-rigid, 1 accessory parameter => HEUN class:
  holonomic/D-finite but NOT elementary). So for OUR N=3, Demkov-Kunike gives holonomic-but-non-rigid (a
  linear ODE replaces the wild isomonodromy), NOT closed form. sigma is the irregular/confluent limit of the
  Heun accessory parameter.
  NOTE (a real subtlety, fixed): the tanh channel amplitude must be read in the EIGENBASIS of H(+-inf)
  (which is constant-but-not-diagonal); the diabatic-basis P oscillates indefinitely -- itself the
  regular-singular e^{i*lambda*u} signature -- and does NOT converge.
  CONCLUSION: LZ's transcendental sigma is the price of the UNBOUNDED linear ramp (the rank-2 irregular
  point at u=inf); bound the ramp and the connection is tame (regular-singular, hypergeometric-class) with
  closed-form / exponentially-computable amplitudes.
"""
from __future__ import annotations
import numpy as np
from scipy.integrate import solve_ivp
np.seterr(all="ignore")


def type1(eps, gam, a):
    eps = np.array(eps, float); gam = np.array(gam, float); a = np.array(a, float); g2 = gam ** 2
    H0 = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            if i != j:
                H0[i, j] = gam[i] * gam[j] * (a[i] - a[j]) / (eps[i] - eps[j])
        H0[i, i] = -sum(g2[k] * (a[i] - a[k]) / (eps[i] - eps[k]) for k in range(3) if k != i)
    return H0, np.diag(a)


def ramp(kind):
    return (lambda u: u) if kind == "linear" else (lambda u: np.tanh(u))


def Pmatrix(H0, A, kind, T, rtol=1e-10):
    """Physical channel amplitude: project the propagator onto the EIGENBASES of H(-+T). For the tanh
    model H(+-inf)=H0+-A is constant-but-NOT-diagonal, so the asymptotic states Rabi-oscillate in the
    diabatic frame; the eigenbasis projection removes that and gives the true in/out channel amplitude.
    (For the linear LZ model the endpoint eigenbasis -> the A-diabatic basis, recovering the usual P.)
    Eigenvalue-sorted bases at both ends, so index 1 = the middle channel (P_mm = P[1,1])."""
    f = ramp(kind)
    _, Vm = np.linalg.eigh(H0 + f(-T) * A)
    _, Vp = np.linalg.eigh(H0 + f(T) * A)
    s = solve_ivp(lambda u, y: (-1j * (H0 + f(u) * A) @ y.reshape(3, 3)).ravel(),
                  [-T, T], np.eye(3, dtype=complex).ravel(), rtol=rtol, atol=1e-12, method="DOP853")
    U = s.y[:, -1].reshape(3, 3)
    S = Vp.conj().T @ U @ Vm
    return np.abs(S) ** 2


def theta_scaling(H0, A, kind):
    """phase theta_01 = int_0^T (E_1-E_0) du at growing T; fit slope of log|theta| vs log T."""
    f = ramp(kind); Ts = np.array([5.0, 10.0, 20.0, 40.0]); th = []
    for T in Ts:
        us = np.linspace(0, T, 4001)
        gap = np.array([np.diff(np.linalg.eigvalsh(H0 + f(u) * A))[0] for u in us])
        th.append(abs(np.trapezoid(gap, us)))
    slope = np.polyfit(np.log(Ts), np.log(np.array(th) + 1e-30), 1)[0]
    return slope, th[-1]


def main():
    eps = [-2.0, 0.0, 3.0]; gam = [1.0, 0.8, 1.2]; a = [-1.0, 0.5, 2.0]
    H0, A = type1(eps, gam, a); mid = 1
    print("=" * 92)
    print("TANH (Demkov-Kunike) TAME analog vs LINEAR (LZ) WILD sweep -- same canonical (gamma,eps,a)")
    print("=" * 92)

    print("\n(A) SINGULARITY TYPE -- eigenvalues at growing u, and the phase scaling theta ~ T^p:")
    for u in [1.0, 5.0, 20.0]:
        wl = np.linalg.eigvalsh(H0 + u * A); wt = np.linalg.eigvalsh(H0 + np.tanh(u) * A)
        print(f"  u={u:5.1f}:  linear E={np.round(wl,2)} (diverge ~a_i u)   tanh E={np.round(wt,3)} (saturate)")
    sl, _ = theta_scaling(H0, A, "linear"); st, _ = theta_scaling(H0, A, "tanh")
    print(f"  phase scaling: theta_linear ~ T^{sl:.2f}  (=>2: irregular, e^(i*lambda*u^2), Fresnel)")
    print(f"                 theta_tanh   ~ T^{st:.2f}  (=>1: regular,   e^(i*lambda*u),  oscillatory)")

    print("\n(B) CONDITIONING -- P_mm vs truncation T (tanh converges exponentially; linear has Fresnel tail):")
    print(f"  {'T':>5} | {'P_mm tanh':>10} {'d(tanh)':>10} | {'P_mm linear':>11} {'d(linear)':>10}")
    prevT = prevL = None
    for T in [4.0, 6.0, 8.0, 12.0, 16.0]:
        pt = Pmatrix(H0, A, "tanh", T)[mid, mid]; pl = Pmatrix(H0, A, "linear", T)[mid, mid]
        dt = "" if prevT is None else f"{pt-prevT:+.2e}"; dl = "" if prevL is None else f"{pl-prevL:+.2e}"
        print(f"  {T:5.1f} | {pt:10.5f} {dt:>10} | {pl:11.5f} {dl:>10}")
        prevT, prevL = pt, pl

    print("\n(C) CLOSED-FORM ANCHOR (N=2 sech/Rosen-Zener; same hypergeometric family):  P = sin^2(Omega*pi/2)")
    sx = np.array([[0, 1], [1, 0]], complex)
    for Om in [0.3, 0.7, 1.0, 1.5]:
        s = solve_ivp(lambda u, y: (-1j * (0.5 * Om / np.cosh(u)) * sx @ y.reshape(2, 2)).ravel(),
                      [-40, 40], np.eye(2, dtype=complex).ravel(), rtol=1e-11, atol=1e-13, method="DOP853")
        P12 = np.abs(s.y[:, -1].reshape(2, 2))[1, 0] ** 2
        exact = np.sin(Om * np.pi / 2) ** 2
        print(f"  Omega={Om:.2f}: numeric P={P12:.6f}  exact sin^2(Omega*pi/2)={exact:.6f}  "
              f"|diff|={abs(P12-exact):.1e}")

    print("\n(D) CONTRAST (N=3 middle survival):")
    pt = Pmatrix(H0, A, "tanh", 16.0)[mid, mid]
    print(f"  P_mm^tanh   = {pt:.5f}  (TAME but HEUN-class: regular-singular, holonomic, non-rigid, converged)")
    print(f"  P_mm^linear = 0.21472   (WILD: the transcendental sigma -- no closed form, the whole project)")
    print("\n  => bounding the ramp (tanh) removes the rank-2 irregular point: spectrum saturates, phase is")
    print("     regular (~u not ~u^2), scattering well-conditioned, amplitude HOLONOMIC (a linear ODE). But")
    print("     for N=3 it is NON-RIGID (Heun, 1 accessory parameter), NOT elementary closed form -- only")
    print("     N=2 (rigid, 2F1) is. sigma is the irregular limit of the Heun accessory parameter.")


if __name__ == "__main__":
    main()
