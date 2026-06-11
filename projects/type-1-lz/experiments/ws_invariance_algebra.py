"""
ws_invariance_algebra.py -- the invariance algebra of P on Type-1 parameter space
{gamma, eps, a}, and the exact functional arity of P (how many algebraic invariants).

GOAL (user): discover the algebraic structure of P anchored in {gamma,eps,a} beyond BE.
SEED TENSION: R11 says P_mm ~ f(Sigma_lo, Sigma_hi, chi) (3 invariants, but only ~2e-3
slice agreement); R24 says rank d(p_lo,p_hi,sigma,b)/d(params) = 4 (incompatible with
exact 3-invariant factorization). This experiment settles the arity.

METHOD
[1] Known P-preserving flows on the 9 parameters (analytic, verified here):
      F1 eps-shift            eps -> eps + c                  (H0 invariant)
      F2 (gamma,eps)-scale    eps -> l*eps, gam -> sqrt(l)gam (H0 invariant)
      F3 a-shift              a -> a + c                      (adds c*u*I: pure phase)
      F4 u-scale              eps -> m*eps, a -> m^2*a        (u -> u/m reparam)
[2] IN-FAMILY CLOSURE probes (the Type-1-specific miracles, if they hold):
      F5 u-shift   : exists (gam',eps') with H0(gam',eps',a) = H0(gam,eps,a) + tau*A ?
      F6 phase     : exists (gam',eps') with H0(gam',eps',a) = H0(gam,eps,a) + c*I ?
    Solved by least-squares over (gam',eps') with gauge fixed (eps'_0=eps_0, gam'_0=gam_0);
    4 unknowns vs 6 equations -> a vanishing residual is a structural fact, not generic.
    Each closure that holds adds one P-preserving flow (u-shift/global phase preserve P).
[3] FIBER-RANK (arity) test at the canonical anchor:
      I = (Sigma_lo, Sigma_hi, chi); K = ker dI (numerically, SVD).
      All flows lie in K (they preserve I).  dim K = 9 - 3 = 6.
      Compute d(sigma,b) restricted to K.  rank > 0 on the non-flow part of K
      => {Sigma_lo,Sigma_hi,chi} is INCOMPLETE => a 4th (or more) invariant exists;
      its gradient direction (printed) is the analytic construction target.

Reproduce: python3 ws_invariance_algebra.py    (imports num_S12 -> assay/oracle)
"""
import os, sys, numpy as np
from scipy.optimize import least_squares
from scipy.integrate import solve_ivp

HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE)
import num_S12 as NS                       # chi + BE machinery (assay-backed)
from assay import Params, Geometry

EPS0=np.array([-2.,0.,3.]); GAM0=np.array([1.,0.8,1.2]); A0=np.array([-1.,0.5,2.])

def H0_of(gam,eps,a):
    gam=np.asarray(gam,float);eps=np.asarray(eps,float);a=np.asarray(a,float)
    H0=np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            if i!=j: H0[i,j]=gam[i]*gam[j]*(a[i]-a[j])/(eps[i]-eps[j])
        H0[i,i]=-sum(gam[k]**2*(a[i]-a[k])/(eps[i]-eps[k]) for k in range(3) if k!=i)
    return H0

# ---------------- P solver (adiabatic projection, Richardson) ----------------
def P_of(gam,eps,a,Rs=(30.,60.),rtol=1e-10,atol=1e-11):
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
    P1,P2=onerun(Rs[0]),onerun(Rs[1])
    return (16*P2-P1)/15

def obs4(gam,eps,a):
    """(p_lo, p_hi, sigma, b) in slope order."""
    P=P_of(gam,eps,a); lo,mid,hi=np.argsort(np.asarray(a,float))
    return np.array([P[lo,lo],P[hi,hi],P[mid,mid],P[hi,lo]])

def invariants3(gam,eps,a):
    """(Sigma_lo, Sigma_hi, chi)."""
    lo,mid,hi=NS.slope_order(a)
    Slo=NS.be_exponent(eps,gam,a,lo,mid)+NS.be_exponent(eps,gam,a,lo,hi)
    Shi=NS.be_exponent(eps,gam,a,hi,mid)+NS.be_exponent(eps,gam,a,hi,lo)
    geo=Geometry(Params(eps=tuple(eps),gam=tuple(gam),a=tuple(a)))
    chi=float(np.real(NS.q4_cross_ratio(geo)))
    return np.array([Slo,Shi,chi])

# ---------------- [2] in-family closure probes ----------------
def closure_probe(target_shift,label,tau):
    """Find (gam',eps') with H0(gam',eps',A0) = H0(GAM0,EPS0,A0) + tau*target_shift.
       Gauge: eps'_0=EPS0[0], gam'_0=GAM0[0]. Unknowns: gam'_1,gam'_2,eps'_1,eps'_2."""
    H_target=H0_of(GAM0,EPS0,A0)+tau*target_shift
    def resid(x):
        g=np.array([GAM0[0],x[0],x[1]]); e=np.array([EPS0[0],x[2],x[3]])
        D=H0_of(g,e,A0)-H_target
        return np.array([D[0,1],D[0,2],D[1,2],D[0,0],D[1,1],D[2,2]])
    x0=np.array([GAM0[1],GAM0[2],EPS0[1],EPS0[2]])
    sol=least_squares(resid,x0,xtol=1e-15,ftol=1e-15,gtol=1e-15)
    res=np.linalg.norm(resid(sol.x))
    g=np.array([GAM0[0],sol.x[0],sol.x[1]]); e=np.array([EPS0[0],sol.x[2],sol.x[3]])
    print(f"  [{label}] tau={tau:+.3f}: residual {res:.2e}  "
          f"{'CLOSES IN-FAMILY' if res<1e-10 else 'does NOT close'}")
    if res<1e-10:
        print(f"        gam'={np.round(g,8)}  eps'={np.round(e,8)}")
        dP=np.max(np.abs(obs4(g,e,A0)-obs4(GAM0,EPS0,A0)))
        print(f"        P-invariance check: max|dP|={dP:.2e}")
    return res<1e-10,(g,e)

# ---------------- [3] fiber-rank ----------------
def fiber_rank():
    th0=np.concatenate([GAM0,EPS0,A0]); h=2e-3
    def unpack(t): return t[0:3],t[3:6],t[6:9]
    # gradients of invariants
    dI=np.zeros((3,9)); dO=np.zeros((4,9))
    for k in range(9):
        tp=th0.copy(); tp[k]+=h; tm=th0.copy(); tm[k]-=h
        dI[:,k]=(invariants3(*unpack(tp))-invariants3(*unpack(tm)))/(2*h)
        dO[:,k]=(obs4(*unpack(tp))-obs4(*unpack(tm)))/(2*h)
    # kernel of dI
    Uk,sk,Vk=np.linalg.svd(dI)
    print(f"  singular values of dI (3x9): {np.round(sk,5)}")
    K=Vk[3:,:].T                         # 9x6 basis of ker dI
    M=dO@K                               # 4x6: observables along the fiber
    Um,sm,Vm=np.linalg.svd(M)
    print(f"  singular values of d(p_lo,p_hi,sigma,b)|ker(dI) : {np.round(sm,6)}")
    print(f"  ==> rank on the fiber = {int(np.sum(sm>50*max(1e-12,sm[-1] if len(sm)>4 else 0)+1e-4))} "
          f"(tol 1e-4; flows account for rank 0)")
    # the direction in the fiber that moves sigma the most:
    sig_row=(dO[2,:]@K)                  # 1x6
    x=K@sig_row/np.linalg.norm(sig_row)  # 9-vector
    print(f"  sigma-moving fiber direction (gamma|eps|a components):")
    print(f"    dgam={np.round(x[0:3],4)} deps={np.round(x[3:6],4)} da={np.round(x[6:9],4)}")
    print(f"  |d sigma| along it: {abs(dO[2,:]@x):.5f}  (vs |grad sigma| = {np.linalg.norm(dO[2,:]):.5f})")
    return sm

if __name__=="__main__":
    print("="*84)
    print("INVARIANCE ALGEBRA OF P ON TYPE-1 PARAMETER SPACE  (canonical anchor)")
    print("="*84)
    print("\n[1] sanity: known flows preserve P (spot check, F2 and F3)")
    base=obs4(GAM0,EPS0,A0)
    l=1.37
    pF2=obs4(np.sqrt(l)*GAM0,l*EPS0,A0)
    pF3=obs4(GAM0,EPS0,A0+0.61)
    print(f"  base (p_lo,p_hi,sigma,b) = {np.round(base,6)}")
    print(f"  F2 scale  max|dP|={np.max(np.abs(pF2-base)):.2e}   F3 a-shift max|dP|={np.max(np.abs(pF3-base)):.2e}")

    print("\n[2] in-family closure probes (the Type-1 miracles?)")
    ok5,_=closure_probe(np.diag(A0),"F5: H0 -> H0 + tau*A  (u-shift)",0.20)
    ok6,_=closure_probe(np.eye(3),"F6: H0 -> H0 + c*I   (global phase)",0.20)
    nflows=4+int(ok5)+int(ok6)
    print(f"  flow count = 4 known + {int(ok5)}+{int(ok6)} closures = {nflows}  "
          f"=> generic arity bound: {9-nflows}")

    print("\n[3] fiber-rank test: is (Sigma_lo,Sigma_hi,chi) complete?")
    sm=fiber_rank()
    print("\nREADING: rank 0 on fiber => arity 3 (R11 exact; R24 rank-4 must be revisited).")
    print("rank>=1 => a 4th invariant exists beyond {Sigma_lo,Sigma_hi,chi}; its gradient")
    print("direction above is the analytic construction target.")
