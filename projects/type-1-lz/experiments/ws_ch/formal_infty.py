"""
WS-CH step 1: formal Thome solutions at the rank-2 irregular point v=inf of the
Laplace-frame system  B'(v) = K(v) B,   K(v) = -i diag(1/a)(H0 - v I).

Goal: extract the diagonal formal data (the WKB exponentials + algebraic prefactor
exponents c_i) and check c_i = sum_{j!=i} s_ij^2 (a_i - a_j) (signed BE), matching
the formal monodromy e^{2 pi i c_i}.

Note: in the v-frame the leading term of K is  -i diag(1/a) * (-v I) = i v diag(1/a).
So the eigen-phases grow like  exp( i v^2 /(2 a_j) ) per channel -- a rank-2 (Poincare
rank 1 for the *system*, level-1 in the Borel sense) irregular point. The algebraic
prefactor exponent is the formal monodromy exponent c_j.
"""
import sympy as sp

def H0_sym(eps, gam, a):
    e=[sp.nsimplify(x) for x in eps]; g=[sp.nsimplify(x) for x in gam]; a=[sp.nsimplify(x) for x in a]
    H0=sp.zeros(3,3)
    for i in range(3):
        for j in range(3):
            if i!=j: H0[i,j]=g[i]*g[j]*(a[i]-a[j])/(e[i]-e[j])
        H0[i,i]=-sum(g[k]**2*(a[i]-a[k])/(e[i]-e[k]) for k in range(3) if k!=i)
    return H0,e,g,a

def s_ij(g,e,i,j):
    return g[i]*g[j]/(e[i]-e[j])

def c_i_target(e,g,a,i):
    return sum(s_ij(g,e,i,j)**2*(a[i]-a[j]) for j in range(3) if j!=i)

def formal_exponents(eps,gam,a):
    """
    Diagonalise the formal solution at v=inf for B' = K(v) B.
    K(v) = i*v*D + K0,  D=diag(1/a),  K0 = -i diag(1/a) H0.
    Substitute B = T(v) exp(Phi(v)) with T -> I, Phi diagonal.
    Leading: phase' ~ i v / a_j  => phase_j = i v^2/(2 a_j).
    Sub-leading algebraic exponent: from the 1/v term of the diagonalised connection.
    We compute it via the standard formal-diagonalisation recursion to O(1/v).
    """
    v=sp.symbols('v')
    H0,e,g,a=H0_sym(eps,gam,a)
    D=sp.diag(*[1/a[j] for j in range(3)])
    K=sp.I*v*D + (-sp.I*D*H0)     # = -i D (H0 - v I)
    # Formal diagonalisation: seek gauge G(v)=I + G1/v + ... s.t. G^{-1} K G - G^{-1}G'
    # is diagonal to the needed order. Leading diagonal = i v D (already diagonal),
    # off-diagonal of K0 must be removed at O(v^0) -> sets G1 (the 1/v gauge term),
    # and the residual O(1/v) diagonal gives the algebraic exponents.
    K0 = -sp.I*D*H0
    Dlead = sp.I*D   # coefficient of v in K
    # off-diagonal part of K0
    G1=sp.zeros(3,3)
    for i in range(3):
        for j in range(3):
            if i!=j:
                # [Dlead, G1]_{ij} = (Dlead_ii - Dlead_jj) G1_ij must cancel K0_{ij}
                # gauge eqn at O(v^0): (i/a_i - i/a_j) G1_ij = -K0_ij  (sign from G^{-1}K0 G - ...)
                denom = Dlead[i,i]-Dlead[j,j]
                G1[i,j]= K0[i,j]/denom
    # residual O(1/v) diagonal = diag of (K0 G1 - G1 K0 ... ) corrections; the algebraic
    # exponent c_j is the 1/v coefficient of the diagonalised connection's diagonal.
    # Standard result: c_j = sum_{k!=j} (K0_jk K0_kj)/(Dlead_jj - Dlead_kk)  (the 2nd-order
    # perturbative diagonal shift), interpreted as the v^{c_j} prefactor exponent after
    # integrating phase_j = i v^2/(2a_j) + c_j log v.
    c=[]
    for j in range(3):
        cj=sum(K0[j,k]*K0[k,j]/(Dlead[j,j]-Dlead[k,k]) for k in range(3) if k!=j)
        c.append(sp.simplify(cj))
    return c, [sp.simplify(c_i_target(e,g,a,i)) for i in range(3)], (e,g,a)

if __name__=="__main__":
    for name,eps,gam,a in [
        ("canonical",(-2,0,3),(1,sp.Rational(4,5),sp.Rational(6,5)),(-1,sp.Rational(1,2),2)),
        ("sampleB",(-1,0,sp.Rational(3,2)),(sp.Rational(9,10),sp.Rational(11,10),sp.Rational(4,5)),(sp.Rational(-7,10),sp.Rational(2,5),sp.Rational(13,10))),
    ]:
        c,ctar,_=formal_exponents(eps,gam,a)
        print(f"\n== {name} ==")
        print(" formal algebraic exponents c_j (from diagonalisation):", c)
        print(" signed-BE target  sum_j s_ij^2(a_i-a_j):           ", ctar)
        print(" match:", [sp.simplify(c[i]-ctar[i])==0 for i in range(3)])
        print(" sum c_j =", sp.simplify(sum(c)))
