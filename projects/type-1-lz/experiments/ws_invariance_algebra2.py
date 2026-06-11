"""
ws_invariance_algebra2.py -- decisive test of the arity-4 conjecture:

   P = F( be_lm, be_mh, be_lh, chi )      [pairwise BE exponents + Q4 cross-ratio]

Follow-up to ws_invariance_algebra.py, which found: fiber sv of d(p,p,sigma,b) on
ker d(Sigma_lo,Sigma_hi,chi) = (0.135, 1.7e-4, ...) -> a 4th invariant exists; the
natural algebraic candidate is the THIRD pairwise BE exponent be_lh (the windows are
Sigma_lo=be_lm+be_lh, Sigma_hi=be_mh+be_lh; R11 used only the sums).

Tests at TWO base points:
 [A] I3 = (Sigma_lo,Sigma_hi,chi)        : fiber rank (reproduce run 1)
 [B] I4 = (be_lm,be_mh,be_lh,chi)        : fiber rank -> AT FLOOR everywhere = conjecture
 [C] explained-direction check: the strong I3-fiber direction X has d(be_lh).X != 0
 If [B] collapses: ker dP is 5-dim = 4 gauge flows + ONE HIDDEN non-gauge P-preserving
 flow (the integrability fingerprint); its direction is printed for analytic recognition.

BE gradients are ANALYTIC (exact); chi and P by central FD (h=2e-3).
Reproduce: python3 ws_invariance_algebra2.py
"""
import os, sys, numpy as np
from scipy.integrate import solve_ivp

HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE)
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

def obs2(gam,eps,a):
    """(sigma, b): the only observables that can move on a BE-fixing fiber."""
    P=P_of(gam,eps,a); lo,mid,hi=np.argsort(np.asarray(a,float))
    return np.array([P[mid,mid],P[hi,lo]])

# -------- analytic BE-exponent gradients --------
def be_and_grad(gam,eps,a,i,j):
    """be_ij = gi^2 gj^2 |ai-aj| / (ei-ej)^2 and its gradient in (gam,eps,a) order (9,)."""
    g=np.asarray(gam,float);e=np.asarray(eps,float);av=np.asarray(a,float)
    d=av[i]-av[j]; s=np.sign(d); De=e[i]-e[j]
    be=g[i]**2*g[j]**2*abs(d)/De**2
    grad=np.zeros(9)
    grad[i]=2*be/g[i]; grad[j]=2*be/g[j]                       # d/dgam
    grad[3+i]=-2*be/De; grad[3+j]= 2*be/De                     # d/deps
    grad[6+i]= g[i]**2*g[j]**2*s/De**2                         # d/da_i
    grad[6+j]=-g[i]**2*g[j]**2*s/De**2
    return be,grad

def chi_of(gam,eps,a):
    geo=Geometry(Params(eps=tuple(eps),gam=tuple(gam),a=tuple(a)))
    return float(np.real(NS.q4_cross_ratio(geo)))

def analyse(name,GAM,EPS,A,h=2e-3):
    print("="*84); print(f"BASE {name}: eps={EPS} gam={GAM} a={A}"); print("="*84)
    th0=np.concatenate([GAM,EPS,A]).astype(float)
    lo,mid,hi=NS.slope_order(A)
    # analytic BE rows
    be_lm,g_lm=be_and_grad(GAM,EPS,A,lo,mid)
    be_mh,g_mh=be_and_grad(GAM,EPS,A,mid,hi)
    be_lh,g_lh=be_and_grad(GAM,EPS,A,lo,hi)
    print(f"  pairwise BE: be_lm={be_lm:.5f} be_mh={be_mh:.5f} be_lh={be_lh:.5f}")
    # FD rows for chi and obs2
    def unpack(t): return t[0:3],t[3:6],t[6:9]
    dchi=np.zeros(9); dO=np.zeros((2,9))
    for k in range(9):
        tp=th0.copy(); tp[k]+=h; tm=th0.copy(); tm[k]-=h
        dchi[k]=(chi_of(*unpack(tp))-chi_of(*unpack(tm)))/(2*h)
        dO[:,k]=(obs2(*unpack(tp))-obs2(*unpack(tm)))/(2*h)
    for tag,rows in (("I3=(Sig_lo,Sig_hi,chi)",[g_lm+g_lh,g_mh+g_lh,dchi]),
                     ("I4=(be_lm,be_mh,be_lh,chi)",[g_lm,g_mh,g_lh,dchi])):
        dI=np.vstack(rows)
        U,s,V=np.linalg.svd(dI)
        K=V[dI.shape[0]:,:].T
        M=dO@K
        sm=np.linalg.svd(M,compute_uv=False)
        print(f"  [{tag}]  fiber sv of d(sigma,b): {np.round(sm,7)}")
        if tag.startswith("I3"):
            # explained-direction check
            row=dO[0,:]@K; X=K@row/ (np.linalg.norm(row)+1e-300)
            print(f"      strong I3-direction: d(be_lh).X = {g_lh@X:+.5f}  "
                  f"(|grad be_lh|={np.linalg.norm(g_lh):.4f})  -> be_lh "
                  f"{'EXPLAINS it' if abs(g_lh@X)>1e-2 else 'does NOT explain it'}")
        else:
            # the hidden flow: direction in ker(dI4) maximizing nothing -- it's the
            # whole kernel minus 4 gauge flows; print an orthonormal basis vector of
            # ker(dI4) orthogonal to the gauge flows
            flows=[]
            g0,e0,a0=GAM.astype(float),EPS.astype(float),A.astype(float)
            flows.append(np.concatenate([np.zeros(3),np.ones(3),np.zeros(3)]))           # eps-shift
            flows.append(np.concatenate([0.5*g0,e0,np.zeros(3)]))                        # (gam,eps)-scale
            flows.append(np.concatenate([np.zeros(3),np.zeros(3),np.ones(3)]))           # a-shift
            flows.append(np.concatenate([np.zeros(3),e0,2*a0]))                          # u-scale
            F=np.vstack(flows).T                       # 9x4
            # project flows out of K's column space
            Q,_=np.linalg.qr(np.hstack([F,K]))
            # hidden directions: in span(K), orthogonal to span(F)
            PF=F@np.linalg.pinv(F)
            Kperp=(np.eye(9)-PF)@K
            Uh,sh,Vh=np.linalg.svd(Kperp,full_matrices=False)
            nh=int(np.sum(sh>1e-8))
            print(f"      ker(dI4) = 5-dim; non-gauge part dim = {nh-0 if nh<=5 else nh} "
                  f"(expect 1): sv {np.round(sh[:6],4)}")
            Xh=Uh[:,0]
            print(f"      HIDDEN FLOW direction (gam|eps|a): dgam={np.round(Xh[0:3],4)} "
                  f"deps={np.round(Xh[3:6],4)} da={np.round(Xh[6:9],4)}")
            print(f"      d(sigma,b) along it: {np.round(dO@Xh,6)}  (floor ~1e-4)")

if __name__=="__main__":
    analyse("canonical",np.array([1.,0.8,1.2]),np.array([-2.,0.,3.]),np.array([-1.,0.5,2.]))
    analyse("second",np.array([0.9,1.1,0.8]),np.array([-1.,0.,1.5]),np.array([-0.7,0.4,1.3]))
    print("\nVERDICT: if I4 fiber sv ~ floor at BOTH points => P = F(be_lm,be_mh,be_lh,chi)")
    print("(numerically-supported arity 4) and the printed hidden flow is a genuine non-gauge")
    print("P-preserving deformation = the Type-1 dynamical-symmetry fingerprint.")
