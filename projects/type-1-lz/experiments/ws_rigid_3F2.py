"""
ws_rigid_3F2.py -- a TAME + RIGID rank-3 connection: the generalized-hypergeometric (3F2) class,
and the demonstration that rigidity => closed-form amplitude for N=3.

Setup: tanh 3-level system i psi' = (H0 + tanh(u) A) psi. In z=(1+tanh u)/2 it is a rank-3 Fuchsian system
with residues  R0=(H0-A)/2i, R1=-(H0+A)/2i, Rinf=-iA  at z=0,1,inf. Katz rigidity:
    rig = (2-3)*9 + sum_i dimZ(R_i),  rigid <=> rig=2  =>  sum dimZ = 11 = 3+3+5
=> rigidity REQUIRES one residue with a repeated eigenvalue (a pseudo-reflection, dimZ=5).

(I) RIGIDITY INDEX (numerical):
     * TWO EQUAL SLOPES a=(a,a,c): Rinf=-i*diag(a,a,c) is a pseudo-reflection (dimZ=5) -> rig=2 (RIGID, 3F2).
     * DISTINCT slopes a=(a,b,c): all residues regular-semisimple (dimZ=3) -> rig=0 (NON-rigid, Heun) -- the
       Type-1 case.
(II) RIGID => CLOSED FORM (cleanest instance, the degenerate-pair / Demkov-Osherov):
     a degenerate pair {0,1} (slope a, scalar block) + level 2 (slope c) with a RANK-1 coupling. A DARK
     combination of the pair decouples exactly (P=1) and the bright-2 sector is a 2-level tanh-DK (rigid 2F1,
     closed form). The full 3-level amplitude = (2-level DK) (+) (dark = identity), verified to machine
     precision. So rig=2 => the N=3 amplitude is closed-form, in contrast to the rig=0 Type-1 sigma.
"""
from __future__ import annotations
import numpy as np
from scipy.integrate import solve_ivp
np.seterr(all="ignore")


def dimZ(M, tol=1e-6):
    """Centralizer dimension of a (semisimple) matrix = sum of (eigenvalue multiplicity)^2."""
    ev = np.linalg.eigvals(M); used = [False] * len(ev); d = 0
    for i in range(len(ev)):
        if used[i]:
            continue
        m = 0
        for j in range(len(ev)):
            if not used[j] and abs(ev[i] - ev[j]) < tol:
                used[j] = True; m += 1
        d += m * m
    return d


def rig_index(H0, A):
    R0 = (H0 - A) / (2j); R1 = -(H0 + A) / (2j); Rinf = -1j * A
    dz = (dimZ(R0), dimZ(R1), dimZ(Rinf))
    return (2 - 3) * 9 + sum(dz), dz


def U_prop(H0, A, T=18.0, rtol=1e-11):
    s = solve_ivp(lambda u, y: (-1j * (H0 + np.tanh(u) * A) @ y.reshape(-1, int(np.sqrt(y.size)))).ravel(),
                  [-T, T], np.eye(H0.shape[0], dtype=complex).ravel(), rtol=rtol, atol=1e-13, method="DOP853")
    n = H0.shape[0]
    return s.y[:, -1].reshape(n, n)


def main():
    print("=" * 92)
    print("A TAME + RIGID rank-3 connection (3F2) and the rigidity => closed-form demonstration")
    print("=" * 92)

    print("\n(I) RIGIDITY INDEX of the tanh 3-level system  rig=(2-3)*9 + sum dimZ(R0,R1,Rinf):")
    # generic H0 (Cauchy-like, just a generic symmetric) to expose only the slope structure
    rng = np.random.default_rng(2)
    M = rng.normal(size=(3, 3)); H0g = (M + M.T) / 2
    for label, a in [("TWO EQUAL SLOPES a=(0.6,0.6,2.0)", [0.6, 0.6, 2.0]),
                     ("DISTINCT SLOPES  a=(-1,0.5,2.0)  [Type-1]", [-1.0, 0.5, 2.0])]:
        A = np.diag(a); rg, dz = rig_index(H0g, A)
        verdict = "RIGID (3F2) => closed form" if rg == 2 else "NON-rigid (Heun) => accessory param / sigma"
        print(f"  {label}:  dimZ(R0,R1,Rinf)={dz}  sum={sum(dz)}  rig={rg}  -> {verdict}")

    print("\n(II) RIGID => CLOSED FORM (degenerate-pair / Demkov-Osherov instance):")
    a, c = 0.6, 2.0; g1, g2, d, d3 = 0.5, 0.7, 0.0, 0.3
    H0 = np.array([[d, 0, g1], [0, d, g2], [g1, g2, d3]], float)
    A = np.diag([a, a, c])
    rg, dz = rig_index(H0, A)
    print(f"  model: degenerate pair {{0,1}} (slope {a}), level 2 (slope {c}), rank-1 coupling g=({g1},{g2})")
    print(f"  rigidity: dimZ={dz} sum={sum(dz)} rig={rg}  (rigid)")
    G = np.hypot(g1, g2)
    bright = np.array([g1, g2, 0.0]) / G
    dark = np.array([g2, -g1, 0.0]) / G
    e3 = np.array([0.0, 0.0, 1.0])
    U = U_prop(H0, A)
    # dark channel decouples exactly:
    P_dark = abs(dark @ U @ dark) ** 2
    # bright-level2 reduced 2-level model and its propagator
    H0b = np.array([[d, G], [G, d3]], float); Ab = np.diag([a, c])
    Ub = U_prop(H0b, Ab)
    # embed: amplitude among {bright, 3} from the full propagator vs the 2-level
    full_b3 = np.array([[bright @ U @ bright, bright @ U @ e3],
                        [e3 @ U @ bright, e3 @ U @ e3]])
    print(f"  P(dark -> dark) = {P_dark:.12f}   (exact 1 => dark decouples)")
    print(f"  max| full {{bright,3}} block  -  2-level DK propagator | = {np.max(np.abs(full_b3 - Ub)):.2e}")
    print(f"  => the 3-level amplitude = (2-level tanh-DK, closed-form 2F1) (+) (dark = identity).")
    print(f"     RIGID (rig=2) => CLOSED FORM. Contrast: distinct-slope (rig=0) is Heun; its linear-ramp")
    print(f"     version is the transcendental sigma (P_mm, the whole project).")


if __name__ == "__main__":
    main()
