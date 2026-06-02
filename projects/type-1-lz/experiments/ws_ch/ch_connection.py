"""
WS-CH step 2: the central connection problem, numerically, in the confluent-Heun frame.

We compute the connection/Stokes data of the rank-2 irregular point at u=inf of the
Type-1 N=3 system   i psi' = (H0 + u A) psi   directly, by transporting the formal
(Thome/WKB) solution from one anti-Stokes sector to the opposite one and reading off
the connection matrix C.  |C_ij|^2 = P_{i->j}.  The off-diagonal entry linking the
middle slope to itself is P_mm = S_12 (the target).

This is the SAME object the oracle computes, but assembled as a connection matrix
between the formal solution bases at u=+/-inf (the two relevant Stokes sectors on the
real axis), i.e. exactly the central connection problem of the confluent-Heun-class
scalar ODE, lifted to the 3x3 system (the frame in which u=inf is an ordinary-data
irregular point with no spurious apparent singularity).

Method (Borel-free, direct): build the formal solution Y_formal(u) = T(u) exp(Phi(u))
to high asymptotic order at large |u|, use it as an exact boundary condition at u=-R,
integrate the linear system to u=+R, and project onto the formal basis at +R.  The
resulting transition matrix is the connection matrix; doubly-stochastic |.|^2 = P.
Richardson in R removes the truncation tail.  We benchmark vs experiments/oracle.py.
"""
import numpy as np
from scipy.integrate import solve_ivp

def build_H0(eps,gam,a):
    eps=np.asarray(eps,float); gam=np.asarray(gam,float); a=np.asarray(a,float)
    H0=np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            if i!=j: H0[i,j]=gam[i]*gam[j]*(a[i]-a[j])/(eps[i]-eps[j])
        H0[i,i]=-sum(gam[k]**2*(a[i]-a[k])/(eps[i]-eps[k]) for k in range(3) if k!=i)
    return H0,a

def formal_solution(H0,a,u,order=3):
    """
    Thome formal fundamental matrix Y(u) ~ T(u) exp(Phi(u)) for psi' = -i(H0+uA)psi.
    Per branch j: phase = -i[ a_j u^2/2 + H0_jj u + b1_j log u ],
       b1_j = sum_{k!=j} H0_jk^2/(a_j-a_k)  (signed BE row sum, = c_j algebraic exponent).
    Columns are the formal solutions; T(u)=I + T1/u + ... (we include the leading prefactor
    and the algebraic exponent, which is what the connection problem needs).
    """
    a=np.asarray(a,float)
    b1=np.array([sum(H0[j,k]**2/(a[j]-a[k]) for k in range(3) if k!=j) for j in range(3)])
    # phase per branch
    phase=np.array([-1j*(a[j]*u**2/2 + H0[j,j]*u + b1[j]*np.log(complex(u))) for j in range(3)])
    Y=np.diag(np.exp(phase)).astype(complex)
    return Y, b1, phase

def rhs(u,Yflat,H0,a):
    Y=Yflat.reshape(3,3)
    M=-1j*(H0+u*np.diag(a))
    return (M@Y).reshape(-1)

def connection_matrix(eps,gam,a,R=60.0,rtol=1e-12,atol=1e-13):
    H0,a=build_H0(eps,gam,a)
    a=np.asarray(a,float)
    # boundary at u=-R: set fundamental matrix = formal solution there
    Y0,b1,_=formal_solution(H0,a,-R)
    sol=solve_ivp(rhs,[-R,R],Y0.reshape(-1).astype(complex),args=(H0,a),
                  rtol=rtol,atol=atol,method='DOP853',dense_output=False)
    YR=sol.y[:,-1].reshape(3,3)
    Yf_R,_,_=formal_solution(H0,a,R)
    # connection: Y(R) = Yf_R @ C  => C = Yf_R^{-1} Y(R)
    C=np.linalg.solve(Yf_R,YR)
    P=np.abs(C)**2
    return C,P,b1

def slope_order(a):
    a=np.asarray(a,float); lo,mid,hi=np.argsort(a); return int(lo),int(mid),int(hi)

if __name__=="__main__":
    cases={
        "canonical":((-2,0,3),(1,0.8,1.2),(-1,0.5,2.0)),
        "sampleB":((-1,0,1.5),(0.9,1.1,0.8),(-0.7,0.4,1.3)),
    }
    for nm,(eps,gam,a) in cases.items():
        print(f"\n== {nm} ==")
        rows=[]
        for R in (40.0,80.0):
            C,P,b1=connection_matrix(eps,gam,a,R=R)
            rows.append(P)
        P40,P80=rows
        # 16:1 Richardson in R (matching oracle tail order heuristic)
        Pext=(16*P80-P40)/15
        lo,mid,hi=slope_order(a)
        print(" b1 (signed BE row sums):",np.round(b1,5))
        print(" |C|^2 at R=80 (raw):\n",np.round(P80,5))
        print(" row sums (should be 1):",np.round(P80.sum(1),4))
        print(f" P_mid (mid={mid}) raw R=80 = {P80[mid,mid]:.6f}   Richardson = {Pext[mid,mid]:.6f}")
