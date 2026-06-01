"""
PA-2 verification (coordinator): the accessory parameter of the Type-1 N=3 connection problem is
ALGEBRAIC. The apparent/accessory point v_* of the scalar Laplace ODE equals E_* = the doubly-
degenerate eigenvalue of H(u) at the universal real NODE u_* (the exact crossing) — a RATIONAL
number in {gamma,eps,a}. Since the whole ODE is built by rational operations from H0 and 1/a_j,
all local data at the rational point v_* is rational => the accessory parameter is NOT a free
transcendental modulus. Hence S_12 = a Painleve-V / confluent-Heun / c=1 connection coefficient with
algebraically-fixed monodromy data (a specific named constant; computable; not elementary).
sympy. Verified: canonical v_*=E_*=-748/375 (u_*=-187/750); sampleB -112147/44000 (-3031/4400).
"""
import sympy as sp
E,u=sp.symbols('E u')
def node_and_accessory(eps,gam,a):
    eps=[sp.Rational(x) for x in eps];gam=[sp.Rational(x) for x in gam];a=[sp.Rational(x) for x in a]
    g2=[g*g for g in gam];H0=sp.zeros(3,3)
    for i in range(3):
        for j in range(3):
            if i!=j:H0[i,j]=gam[i]*gam[j]*(a[i]-a[j])/(eps[i]-eps[j])
        H0[i,i]=-sum(g2[k]*(a[i]-a[k])/(eps[i]-eps[k]) for k in range(3) if k!=i)
    chi=sp.expand((E*sp.eye(3)-(H0+u*sp.diag(*a))).det())
    us=[r for r,m in sp.roots(sp.Poly(sp.discriminant(chi,E),u)).items() if m==2 and r.is_rational][0]
    Es=[r for r,m in sp.roots(sp.Poly(chi.subs(u,us),E)).items() if m==2][0]
    return us,Es   # node u_*, accessory v_*=E_* (both rational)
if __name__=="__main__":
    for nm,(e,g,a) in {"canonical":((-2,0,3),(1,sp.Rational(4,5),sp.Rational(6,5)),(-1,sp.Rational(1,2),2)),
                       "sampleB":((-1,0,sp.Rational(3,2)),(sp.Rational(9,10),sp.Rational(11,10),sp.Rational(4,5)),(sp.Rational(-7,10),sp.Rational(2,5),sp.Rational(13,10)))}.items():
        us,Es=node_and_accessory(e,g,a)
        print(f"{nm}: node u_*={us}  accessory v_*=E_*={Es}  (both rational => PA-2 ALGEBRAIC)")
