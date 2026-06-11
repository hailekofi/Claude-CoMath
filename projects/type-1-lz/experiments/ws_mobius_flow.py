"""
ws_mobius_flow.py -- identify the hidden non-gauge P-preserving flow (R31) as the
SPECIAL (inversion) generator of the fractional-linear eps-covariance of the Type-1
family (OWY), with Mobius weight on gamma:

   eps_i -> eps_i/(1-tau eps_i),  gam_i -> gam_i/(1-tau eps_i),  a fixed.

Analytic facts (verified exactly here):
 (i)  s_ij = gam_i gam_j/(eps_i-eps_j) invariant EXACTLY (finite tau) -> off-diag H0,
      all pairwise be_ij invariant.
 (ii) H0(gam',eps',a) = H0 + tau[ G(tau) A - T(tau) I ],  G=sum gam_k^2/(1-tau eps_k),
      T=sum gam_k^2 a_k/(1-tau eps_k)  -> u-translation + global phase -> P INVARIANT.
 => proves P-invariance along the flow; with F1-F4 this DERIVES arity <= 4 and hence
    P = F(be_lm, be_mh, be_lh, chi)   analytically.
Checks: [1] cosine(predicted generator, measured hidden Y) after gauge-orthogonalization;
[2] closure identity (ii) to machine precision; [3] finite-flow P-invariance to solver
floor; [4] chi invariance along the flow.
"""
import os,sys,numpy as np
from scipy.integrate import solve_ivp
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
import num_S12 as NS
from assay import Params, Geometry

def H0_of(gam,eps,a):
    gam=np.asarray(gam,float);eps=np.asarray(eps,float);a=np.asarray(a,float)
    H0=np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            if i!=j: H0[i,j]=gam[i]*gam[j]*(a[i]-a[j])/(eps[i]-eps[j])
        H0[i,i]=-sum(gam[k]**2*(a[i]-a[k])/(eps[i]-eps[k]) for k in range(3) if k!=i)
    return H0

def P_of(gam,eps,a,Rs=(30.,60.),rtol=1e-11,atol=1e-12):
    H0=H0_of(gam,eps,a); a=np.asarray(a,float)
    def onerun(R):
        w,V=np.linalg.eigh(H0+(-R)*np.diag(a))
        sol=solve_ivp(lambda u,Y:(-1j*(H0+u*np.diag(a))@Y.reshape(3,3)).reshape(-1),
                      [-R,R],V.astype(complex).reshape(-1),rtol=rtol,atol=atol,method='DOP853')
        YR=sol.y[:,-1].reshape(3,3)
        w2,V2=np.linalg.eigh(H0+R*np.diag(a))
        C=V2.conj().T@YR; P=np.abs(C)**2
        sin=[int(np.argmin(np.abs(a-w[k]/(-R)))) for k in range(3)]
        sout=[int(np.argmin(np.abs(a-w2[k]/R))) for k in range(3)]
        Pd=np.zeros((3,3))
        for ki in range(3):
            for ko in range(3): Pd[sin[ki],sout[ko]]=P[ki,ko]
        return Pd
    return (16*onerun(Rs[1])-onerun(Rs[0]))/15

def chi_of(gam,eps,a):
    geo=Geometry(Params(eps=tuple(eps),gam=tuple(gam),a=tuple(a)))
    return float(np.real(NS.q4_cross_ratio(geo)))

CASES={
 "canonical":(np.array([1.,0.8,1.2]),np.array([-2.,0.,3.]),np.array([-1.,0.5,2.]),
   np.array([-0.4157,-0.0763,0.4622, 0.3865,-0.6312,0.2447, -0.0065,0.,0.0065])),
 "second":(np.array([0.9,1.1,0.8]),np.array([-1.,0.,1.5]),np.array([-0.7,0.4,1.3]),
   np.array([-0.5684,-0.133,0.5162, 0.2573,-0.5105,0.2532, 0.0315,-0.0021,-0.0295])),
}

for nm,(GAM,EPS,A,Ymeas) in CASES.items():
    print("="*80); print(f"CASE {nm}"); print("="*80)
    # [1] cosine between gauge-orthogonalized predicted generator and measured Y
    V=np.concatenate([GAM*EPS, EPS**2, np.zeros(3)])     # (gam eps, eps^2, 0)
    F=np.stack([np.concatenate([np.zeros(3),np.ones(3),np.zeros(3)]),
                np.concatenate([0.5*GAM,EPS,np.zeros(3)]),
                np.concatenate([np.zeros(3),np.zeros(3),np.ones(3)]),
                np.concatenate([np.zeros(3),EPS,2*A])],1)  # 9x4 gauge flows
    PF=F@np.linalg.pinv(F)
    Vp=(np.eye(9)-PF)@V; Vp/=np.linalg.norm(Vp)
    Ym=(np.eye(9)-PF)@Ymeas; Ym/=np.linalg.norm(Ym)
    print(f"[1] cos(angle) predicted-vs-measured hidden flow = {abs(Vp@Ym):.6f}   (1=identified)")
    # [2] closure identity at finite tau
    tau=0.05
    g2=GAM/(1-tau*EPS); e2=EPS/(1-tau*EPS)
    G=np.sum(GAM**2/(1-tau*EPS)); T=np.sum(GAM**2*A/(1-tau*EPS))
    D=H0_of(g2,e2,A)-(H0_of(GAM,EPS,A)+tau*(G*np.diag(A)-T*np.eye(3)))
    print(f"[2] closure |H0' - (H0 + tau[G A - T I])| = {np.max(np.abs(D)):.2e}   (0=exact)")
    # [3] finite-flow P-invariance
    P0=P_of(GAM,EPS,A); P1=P_of(g2,e2,A)
    print(f"[3] max|P' - P| along finite flow (tau=0.05) = {np.max(np.abs(P1-P0)):.2e}  (solver floor)")
    # [4] chi invariance
    c0=chi_of(GAM,EPS,A); c1=chi_of(g2,e2,A)
    print(f"[4] |chi' - chi| = {abs(c1-c0):.2e}")
