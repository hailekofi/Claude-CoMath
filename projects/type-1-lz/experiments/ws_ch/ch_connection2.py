"""
WS-CH step 2 (corrected): central connection problem in the ADIABATIC (WKB) frame.

The Thome formal solution of  i psi' = (H0+uA)psi  at the rank-2 irregular point u=inf,
in flux-normalised form, is the WKB/adiabatic solution
    psi_j(u) ~ v_j(u) exp(-i \int^u E_j(u') du'),
with E_j(u), v_j(u) the instantaneous eigenpairs of H(u)=H0+uA.  This is the correct
normalisation in which the connection matrix C (psi^{out}_k expanded in psi^{in}_j) is
unitary, so |C_kj|^2 is doubly-stochastic = P_{j->k}.

We compute C as the central connection coefficient by transporting the adiabatic frame
from u=-R to u=+R (interaction picture), exactly the confluent-Heun central connection
problem of the scalar reduction lifted to the system.  Richardson in R.  Benchmark vs
experiments/oracle.py anchors.

This is the WS-CH route's *computable* deliverable: a connection-matrix value of the
genus-0 confluent-Heun-class irregular point, benchmarked.
"""
import numpy as np
from scipy.integrate import solve_ivp

def build_H(eps,gam,a):
    eps=np.asarray(eps,float); gam=np.asarray(gam,float); a=np.asarray(a,float)
    H0=np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            if i!=j: H0[i,j]=gam[i]*gam[j]*(a[i]-a[j])/(eps[i]-eps[j])
        H0[i,i]=-sum(gam[k]**2*(a[i]-a[k])/(eps[i]-eps[k]) for k in range(3) if k!=i)
    return H0,a

def adiabatic(H0,a,u):
    """instantaneous eigenpairs of H(u), sorted by eigenvalue (ascending)."""
    H=H0+u*np.diag(a)
    w,V=np.linalg.eigh(H)
    return w,V   # columns V[:,k] eigenvector for w[k]

def connection_adiabatic(eps,gam,a,R=200.0,rtol=1e-12,atol=1e-13,Nphase=4000):
    """
    Transport in the adiabatic interaction picture from -R to +R.
    State in diabatic basis psi; we propagate psi' = -i H(u) psi for a full fundamental
    matrix started in the adiabatic basis at -R, then project on the adiabatic basis at +R.
    Connection: C = Vout^{dagger-flux} ... actually C[k,j] = <out_k | propagated in_j>.
    """
    H0,a=build_H(eps,gam,a)
    a=np.asarray(a,float)
    win,Vin=adiabatic(H0,a,-R)
    # propagate the full fundamental matrix with columns = Vin (the incoming adiabatic basis)
    def rhs(u,Yf):
        Y=Yf.reshape(3,3); return (-1j*(H0+u*np.diag(a))@Y).reshape(-1)
    sol=solve_ivp(rhs,[-R,R],Vin.astype(complex).reshape(-1),
                  rtol=rtol,atol=atol,method='DOP853')
    YR=sol.y[:,-1].reshape(3,3)
    wout,Vout=adiabatic(H0,a,R)
    # C[k,j] = component of propagated incoming-j along outgoing-k = Vout^dagger YR
    C=Vout.conj().T@YR
    P=np.abs(C)**2
    return C,P,(win,wout,Vin,Vout)

def slope_map(eps,gam,a,R=200.0):
    """which adiabatic sheet -> which diabatic slope channel at +-R (by eigen-slope)."""
    H0,a=build_H(eps,gam,a); a=np.asarray(a,float)
    win,_=adiabatic(H0,a,-R); wout,_=adiabatic(H0,a,R)
    # incoming: sheet k slope ~ win/( -R ) -> nearest a; outgoing: wout/R -> nearest a
    sin=[int(np.argmin(np.abs(a-win[k]/(-R)))) for k in range(3)]
    sout=[int(np.argmin(np.abs(a-wout[k]/R))) for k in range(3)]
    return sin,sout

if __name__=="__main__":
    cases={"canonical":((-2,0,3),(1,0.8,1.2),(-1,0.5,2.0)),
           "sampleB":((-1,0,1.5),(0.9,1.1,0.8),(-0.7,0.4,1.3))}
    for nm,(eps,gam,a) in cases.items():
        print(f"\n== {nm} ==")
        Ps={}
        for R in (120.0,240.0):
            C,P,_=connection_adiabatic(eps,gam,a,R=R)
            Ps[R]=P
        Pext=(16*Ps[240.0]-Ps[120.0])/15
        sin,sout=slope_map(eps,gam,a)
        print(" incoming sheet->slope:",sin," outgoing sheet->slope:",sout)
        print(" |C|^2 (R=240, adiabatic-sheet basis):\n",np.round(Ps[240.0],6))
        print(" row sums:",np.round(Ps[240.0].sum(1),6)," col sums:",np.round(Ps[240.0].sum(0),6))
        # diabatic P[x_in_slope, j_out_slope] = P_adiab[k_in, k_out] with slope maps
        Pdia=np.zeros((3,3))
        for ki in range(3):
            for ko in range(3):
                Pdia[sin[ki],sout[ko]]=Pext[ki,ko]
        a_arr=np.asarray(a,float); lo,mid,hi=np.argsort(a_arr)
        print(" diabatic P (slope channels), Richardson:\n",np.round(Pdia,6))
        print(f" P_mid (mid slope={mid}) = {Pdia[mid,mid]:.6f}")
