"""
Ring structure of the Type-1 commuting family and its bearing on the
time-quadratic commuting-partner approach.  Standalone (numpy only).

Establishes three facts (all to machine precision):
  (1) {H^(a)(u) = H0(a) + u diag(a)} (fixed gam,eps; varying slope a) is a
      mutually COMMUTING family  -> a commutative ring at each u.
  (2) It spans the full 3-dim commutant of H(u); hence the time-quadratic
      commuting partner H(u)^2 = sum_k c_k(u) H^(a_k)(u) with c_k(u) LINEAR in u.
  (3) Commuting => the whole family shares the SAME eigenvectors phi_i(u);
      members differ only in eigenvalues E_i^(a)(u), which are linear in a.

Consequence: in the adiabatic basis the Schrodinger generator is
      G^(a)(u) = diag(E^(a)(u)) - i W(u),   W_ij = <phi_i|d phi_j/du>,
with W(u) a-INDEPENDENT (fixed geometric/Stuckelberg connection) and the
diagonal a-LINEAR.  So P^(a) = |holonomy of (fixed W) twisted by (a-linear phases)|^2.
"""
import numpy as np


def H_of(eps, gam, a, u):
    eps = np.asarray(eps, float); gam = np.asarray(gam, float); a = np.asarray(a, float)
    g2 = gam ** 2
    H0 = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            if i != j:
                H0[i, j] = gam[i] * gam[j] * (a[i] - a[j]) / (eps[i] - eps[j])
        H0[i, i] = -sum(g2[k] * (a[i] - a[k]) / (eps[i] - eps[k]) for k in range(3) if k != i)
    return H0 + u * np.diag(a)


def report(eps, gam, A):
    # (1) mutual commutativity
    mc = max(np.max(np.abs(H_of(eps, gam, A[i], u) @ H_of(eps, gam, A[j], u)
                          - H_of(eps, gam, A[j], u) @ H_of(eps, gam, A[i], u)))
             for u in (-3.1, 0.4, 2.7) for i in range(3) for j in range(i + 1, 3))
    # (2) reducibility of H^2 with linear coefficients
    us = np.linspace(-4, 4, 17); C = []
    for u in us:
        B = np.stack([H_of(eps, gam, A[k], u).reshape(-1) for k in range(3)]).T
        Hsq = (H_of(eps, gam, A[0], u) @ H_of(eps, gam, A[0], u)).reshape(-1)
        c, *_ = np.linalg.lstsq(B, Hsq, rcond=None)
        C.append((np.linalg.norm(B @ c - Hsq), c))
    resid = max(r for r, _ in C)
    Cmat = np.array([c for _, c in C])
    u2 = max(abs(np.polyfit(us, Cmat[:, k], 3)[1]) for k in range(3))  # u^2 coeff magnitude
    # (3) shared eigenvectors
    sv = 0
    for u in (-2.0, 0.6, 3.3):
        _, V0 = np.linalg.eigh(H_of(eps, gam, A[0], u))
        for k in (1, 2):
            _, Vk = np.linalg.eigh(H_of(eps, gam, A[k], u))
            sv = max(sv, 1 - np.abs(V0.T @ Vk).max(axis=0).min())
    print(f"(1) max||[H^a,H^b]||           = {mc:.2e}   (commuting family)")
    print(f"(2) H^2 reconstruction residual = {resid:.2e}   (in linear-family span)")
    print(f"    max |u^2 coeff of c_k(u)|   = {u2:.2e}   (c_k LINEAR in u)")
    print(f"(3) 1 - min max eigvec overlap  = {sv:.2e}   (shared eigenbasis)")


if __name__ == "__main__":
    eps = (-2., 0., 3.); gam = (1., .8, 1.2)
    A = [np.array([-1., .5, 2.]), np.array([0.3, -0.7, 1.4]), np.array([2., 1., -0.5])]
    report(eps, gam, A)
