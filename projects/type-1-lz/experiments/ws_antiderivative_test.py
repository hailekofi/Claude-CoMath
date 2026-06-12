"""
ws_antiderivative_test.py -- is sigma an ANTIDERIVATIVE of a well-behaved function of a
parameter? (User challenge; corrects R27's over-strong reading.)

KEY FACT: erf = int e^{-t^2} is Liouvillian (quadrature of elementary) yet FAILS order-1
DA and passes order-2 -- the same signature sigma showed in R27. So R27 did NOT exclude
the antiderivative class. DISCRIMINATOR: the DERIVATIVE. If f = int(elementary exp-closure)
then f' IS order-1 DA (erf': g'=-2tg). Painleve-class f has f' ALSO order-1-failing.

Battery (coupling-scale slice, sigma to ~1e-6; ALL controls noise-matched at 1e-6):
  y_erf  = A*erf(...)+c      : order-1 MISS expected  (validates: antiderivatives evade R27)
  y_erf' = Gaussian          : order-1 HIT expected   (the discriminator works)
  J0(9t)                     : order-1 MISS           (holonomic benchmark)
  J0'(9t) = -9 J1(9t)        : order-1 MISS           (derivative of holonomic stays miss)
  sigma                      : order-1 MISS           (R27, reproduced)
  sigma'                     : ???   <- THE NEW DATUM
     HIT  => sigma' elementary on this slice: chase the closed form (bombshell)
     MISS => the antiderivative hypothesis fails one quadrature-depth deeper;
             supports Painleve-class (and analytic R26, where quadrature is a
             'classical' operation killed by irreducibility).
Reproduce: python3 ws_antiderivative_test.py
"""
import numpy as np
from scipy.integrate import solve_ivp
from numpy.polynomial import chebyshev as C
from itertools import combinations_with_replacement
from scipy.special import j0,j1,erf

rng=np.random.default_rng(2)

def build(eps,gam,a):
    eps=np.asarray(eps,float);gam=np.asarray(gam,float);a=np.asarray(a,float)
    H0=np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            if i!=j: H0[i,j]=gam[i]*gam[j]*(a[i]-a[j])/(eps[i]-eps[j])
        H0[i,i]=-sum(gam[k]**2*(a[i]-a[k])/(eps[i]-eps[k]) for k in range(3) if k!=i)
    return H0,a

def Pmid(eps,gam,a,Rs=(30.,60.),rtol=1e-10,atol=1e-11):
    H0,a=build(eps,gam,a)
    def onerun(R):
        w,V=np.linalg.eigh(H0+(-R)*np.diag(a))
        sol=solve_ivp(lambda u,Y:(-1j*(H0+u*np.diag(a))@Y.reshape(3,3)).reshape(-1),
                      [-R,R],V.astype(complex).reshape(-1),rtol=rtol,atol=atol,method='DOP853')
        YR=sol.y[:,-1].reshape(3,3)
        w2,V2=np.linalg.eigh(H0+R*np.diag(a))
        P=np.abs(V2.conj().T@YR)**2
        sin=[int(np.argmin(np.abs(a-w[k]/(-R)))) for k in range(3)]
        sout=[int(np.argmin(np.abs(a-w2[k]/R))) for k in range(3)]
        Pd=np.zeros((3,3))
        for ki in range(3):
            for ko in range(3): Pd[sin[ki],sout[ko]]=P[ki,ko]
        return Pd
    P1,P2=onerun(Rs[0]),onerun(Rs[1]); Pe=(16*P2-P1)/15
    lo,mid,hi=np.argsort(a); return Pe[mid,mid]

def monomials(X,deg):
    n,d=X.shape; cols=[]
    for o in range(0,deg+1):
        for idx in combinations_with_replacement(range(d),o):
            c=np.ones(n)
            for j in idx: c=c*X[:,j]
            cols.append(c)
    return np.array(cols).T

def smin(jet,deg):
    M=monomials(jet,deg); nm=np.linalg.norm(M,axis=0); nm[nm==0]=1
    s=np.linalg.svd(M/nm,compute_uv=False); return s[-1]/s[0]

def da_row(name,tn,fn,degfit,kshift=0):
    """order-1 test on the (kshift)-th derivative of f: jet (t, f^(k), f^(k+1))."""
    a,b=tn.min(),tn.max()
    xn=(2*tn-(a+b))/(b-a)
    coef=C.chebfit(xn,fn,degfit)
    te=np.linspace(a+0.06*(b-a),b-0.06*(b-a),120)
    xe=(2*te-(a+b))/(b-a); sc=2.0/(b-a)
    ders=[C.chebval(xe,coef)]
    ck=coef.copy()
    for _ in range(kshift+2):
        ck=C.chebder(ck)*sc; ders.append(C.chebval(xe,ck))
    jet=np.column_stack([te,ders[kshift],ders[kshift+1]])
    r=[smin(jet,d) for d in (2,3)]
    print(f"  {name:34s} order-1 on f^({kshift}):  d2 {r[0]:.2e}   d3 {r[1]:.2e}")
    return r

if __name__=="__main__":
    eps=(-2,0,3); gam0=np.array([1.,0.8,1.2]); a3=(-1.,0.5,2.)
    Nn=35; degfit=26
    xn=np.cos(np.pi*np.arange(Nn)/(Nn-1))
    t0,t1=0.45,2.3
    tn=0.5*(t1-t0)*(xn[::-1]+1)+t0
    print("building sigma(lambda-coupling slice), %d nodes ..."%Nn)
    sn=np.array([Pmid(eps,lam*gam0,a3) for lam in tn])
    print(f"  sigma range [{sn.min():.4f},{sn.max():.4f}]")
    noise=rng.standard_normal(Nn)*1e-6      # noise-match all controls to sigma's floor
    yerf=0.5*erf(1.2*(tn-1.2))+0.4+noise
    yj0=j0(9*tn)+noise
    print("\nDA order-1 battery (decisive cells d2/d3; calibrate by controls):")
    print("-"*78)
    da_row("CONTROL erf-type (antiderivative)",tn,yerf,degfit,kshift=0)
    da_row("CONTROL erf-type DERIVATIVE",tn,yerf,degfit,kshift=1)
    da_row("CONTROL J0 (holonomic)",tn,yj0,degfit,kshift=0)
    da_row("CONTROL J0 DERIVATIVE",tn,yj0,degfit,kshift=1)
    da_row("SIGMA",tn,sn,degfit,kshift=0)
    da_row("SIGMA DERIVATIVE  <-- the new datum",tn,sn,degfit,kshift=1)
    print("-"*78)
    print("READING: erf-control: f MISS / f' HIT = the discriminator. If sigma' sits at the")
    print("erf'-HIT level => sigma IS an antiderivative of an elementary function (chase it!).")
    print("If sigma' sits at the J0'-MISS level => antiderivative hypothesis excluded one")
    print("quadrature deep; consistent with Painleve class (R26).")
