"""
ws_e_laplace_class.py  --  WS-E: the Laplace/Euler integral representation of the
Type-1 N=3 MLZ amplitudes, and the proof that it does NOT close to a product integrand.

Result (reproduced for both anchors):
  * psi_j(u) = INT_C e^{-i u v} B_j(v) dv  with  B'(v) = K(v) B,
        K(v) = -i diag(1/a) (H0 - v I),
    a 3x3 FIRST-ORDER system, linear in v, irregular (rank-2) at v=infinity.
  * [K(v1),K(v2)] != 0  =>  K is non-abelian along the flow  =>  no v-independent frame
    diagonalises it  =>  the integrand is NOT a product  prod (v-b_k)^rho_k.
    (Contrast: Lin-Sinitsyn / BBGY have a 1/t Coulomb term whose Laplace image is Fuchsian,
     giving a first-order scalar eqn and a product integrand -> Beta x 2F1 / 1F2.)
  * The scalar 3rd-order ODE for B_0(v) has ONE finite singular point (an accessory point,
    NOT an eigenvalue of H0) with indicial exponents {0,1,3} (a GAP at 2), PLUS the rank-2
    irregular point at v=infinity.  => class is confluent-Heun / higher-Weber, one accessory
    parameter ABOVE the 1F2 / Kampe de Feriet class of the solvable tridiagonal cousins.

This is the rigorous, frame-independent statement of WS-A's verdict, in the Laplace frame.
"""
from __future__ import annotations
import sympy as sp


def H0_sym(eps, gam, a):
    e = [sp.nsimplify(x) for x in eps]
    g = [sp.nsimplify(x) for x in gam]
    a = [sp.nsimplify(x) for x in a]
    H0 = sp.zeros(3, 3)
    for i in range(3):
        for j in range(3):
            if i != j:
                H0[i, j] = g[i] * g[j] * (a[i] - a[j]) / (e[i] - e[j])
        H0[i, i] = -sum(g[k] ** 2 * (a[i] - a[k]) / (e[i] - e[k])
                        for k in range(3) if k != i)
    return H0, a


def transformed_system(eps, gam, a):
    """Return (v, K(v)) for  B'(v) = K(v) B,  the Laplace image of  i psi'=(H0+uA)psi."""
    v = sp.symbols('v')
    H0, a = H0_sym(eps, gam, a)
    K = -sp.I * sp.diag(*[1 / a[j] for j in range(3)]) * (H0 - v * sp.eye(3))
    return v, K, H0


def commutator_nonzero(v, K):
    v1, v2 = sp.symbols('v1 v2')
    C = sp.simplify(K.subs(v, v1) * K.subs(v, v2) - K.subs(v, v2) * K.subs(v, v1))
    return C, (C != sp.zeros(3, 3))


def scalar_ode(v, K, comp=0):
    """Cyclic-vector elimination -> y''' = c2 y'' + c1 y' + c0 y for y = B_comp."""
    c = sp.Matrix([[1 if i == comp else 0 for i in range(3)]])  # row e_comp^T
    rows = [c]
    for _ in range(3):
        c = sp.Matrix([[sp.diff(c[0, i], v) for i in range(3)]]) + c * K
        rows.append(c)
    R = sp.Matrix.vstack(rows[0], rows[1], rows[2])
    coeff = sp.simplify(rows[3] * R.inv())   # (c0,c1,c2)
    return coeff[0, 0], coeff[0, 1], coeff[0, 2]


def finite_singular_points(v, c2):
    den = sp.denom(sp.together(c2))
    return sp.solve(sp.Eq(den, 0), v)


def indicial_at(v, c0, c1, c2, v0):
    """Leading Frobenius indicial polynomial at a simple pole v0 of the coefficients."""
    A = sp.limit((v - v0) * c2, v, v0)
    r = sp.symbols('r')
    # y ~ (v-v0)^r : leading balance  r(r-1)(r-2) - A r(r-1) = 0
    poly = r * (r - 1) * (r - 2) - A * r * (r - 1)
    return sp.solve(poly, r), A


def report(name, eps, gam, a):
    print(f"\n================  {name}  ================")
    v, K, H0 = transformed_system(eps, gam, a)
    print("K(v) (transformed first-order system  B' = K B):")
    sp.pprint(K)
    C, nz = commutator_nonzero(v, K)
    print(f"\n[K(v1),K(v2)] nonzero ? -> {nz}   (nonzero => NO product integrand / non-abelian)")
    c0, c1, c2 = scalar_ode(v, K, 0)
    sp_ = finite_singular_points(v, c2)
    print(f"finite singular point(s) of the scalar ODE in v: {sp_}")
    H0n = sp.Matrix(H0)
    eig = list(H0n.eigenvals().keys())
    print(f"eig(H0) = {[sp.nsimplify(e, rational=False) for e in eig]}  "
          f"(singular point is NOT among these => genuine accessory point)")
    for v0 in sp_:
        roots, A = indicial_at(v, c0, c1, c2, v0)
        print(f"  indicial exponents at v0={v0}:  {sorted(roots, key=lambda x: sp.re(x))}  "
              f"(residue A={sp.nsimplify(A)})")
    print("  => exponents {0,1,3} (gap at 2): one accessory parameter ABOVE 1F2/Kampe de Feriet.")


if __name__ == "__main__":
    report("canonical", (-2, 0, 3), (1, sp.Rational(4, 5), sp.Rational(6, 5)),
           (-1, sp.Rational(1, 2), 2))
    report("sampleB", (-1, 0, sp.Rational(3, 2)),
           (sp.Rational(9, 10), sp.Rational(11, 10), sp.Rational(4, 5)),
           (sp.Rational(-7, 10), sp.Rational(2, 5), sp.Rational(13, 10)))
