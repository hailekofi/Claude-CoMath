"""
ws_painleve_da_test.py -- differential-algebraic ORDER test of sigma(parameter), the
numerically-accessible shadow of the "sigma is a Painleve/Garnier transcendent" claim.

LOGIC (the honest hierarchy of differential-algebraic complexity):
  elementary / Liouvillian  : satisfies a LOW-order algebraic ODE; in particular often an
                              ORDER-1 algebraic ODE  P(t, f, f') = 0  (e.g. f=e^{g}, g alg).
  holonomic / classical-spec: order-1 fails, but an ORDER-2 *linear* algebraic ODE holds
                              (e.g. Bessel J0).
  Painleve/Garnier transcdt : order-1 fails AND no order-1 algebraic ODE exists; an order-2
                              *nonlinear* algebraic ODE holds only along the isomonodromic
                              time; on a generic slice it is typically hyper-transcendental.
  hyper-transcendental      : satisfies NO algebraic ODE of any bounded order/degree.

TEST: along a 1-parameter slice t = a_mid (others fixed), compute sigma(t) to ~1e-6, build a
Chebyshev interpolant, get sigma', sigma'' by spectral differentiation, and for each
(order k, degree d) test whether a nonzero polynomial relation P(t, sigma, ..., sigma^(k))=0
exists -- via the smallest singular value of the column-normalized monomial-jet matrix.
  s_min ~ floor  => such an algebraic ODE EXISTS at that (k,d).
  s_min ~ O(1)   => none exists at that (k,d).

CONTROLS validate the detector at the SAME pipeline/accuracy:
  exp  (Liouvillian):  expect order-1, deg-2 hit (s_min -> floor).
  J0   (holonomic):    expect order-1 miss, order-2 hit.
The CLAIM the derivation needs: sigma FAILS order-1 at every tested degree (non-Liouvillian/
non-classical), distinguishing it from the exp control. Failing higher orders too is
consistent with a non-rigid isomonodromy transcendent on a non-isomonodromic slice.

HONEST CAVEATS: (i) failing at bounded (k,d) is strong evidence, not proof; (ii) the
accuracy floor caps the max order resolvable; (iii) "order-2 fails on this slice" does NOT
contradict the Painleve claim -- a generic parameter slice need not be the isomonodromic
time. The robust deliverable is the ORDER-1 FALSIFICATION = sigma non-Liouvillian-classical.

Reproduce: python3 ws_painleve_da_test.py   (numpy/scipy only)
"""
import numpy as np
from scipy.integrate import solve_ivp
from numpy.polynomial import chebyshev as C
from itertools import combinations_with_replacement
from scipy.special import j0

# ----------------------------------------------------------------- Type-1 solver
def build_H(eps,gam,a):
    eps=np.asarray(eps,float); gam=np.asarray(gam,float); a=np.asarray(a,float)
    H0=np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            if i!=j: H0[i,j]=gam[i]*gam[j]*(a[i]-a[j])/(eps[i]-eps[j])
        H0[i,i]=-sum(gam[k]**2*(a[i]-a[k])/(eps[i]-eps[k]) for k in range(3) if k!=i)
    return H0,a

def adiabatic_frame(H0,a,u):
    w,V=np.linalg.eigh(H0+u*np.diag(a)); return w,V

def connection_ip(eps,gam,a,R,rtol=1e-12,atol=1e-13):
    """P (diabatic, slope-ordered) via adiabatic projection at +/-R."""
    H0,a=build_H(eps,gam,a); a=np.asarray(a,float)
    win,Vin=adiabatic_frame(H0,a,-R)
    sol=solve_ivp(lambda u,Yf:(-1j*(H0+u*np.diag(a))@Yf.reshape(3,3)).reshape(-1),
                  [-R,R],Vin.astype(complex).reshape(-1),rtol=rtol,atol=atol,method='DOP853')
    YR=sol.y[:,-1].reshape(3,3)
    wout,Vout=adiabatic_frame(H0,a,R)
    Cmat=Vout.conj().T@YR
    P=np.abs(Cmat)**2
    sin=[int(np.argmin(np.abs(a-win[k]/(-R)))) for k in range(3)]
    sout=[int(np.argmin(np.abs(a-wout[k]/R))) for k in range(3)]
    Pdia=np.zeros((3,3))
    for ki in range(3):
        for ko in range(3): Pdia[sin[ki],sout[ko]]=P[ki,ko]
    return Pdia

def sigma_of(eps,gam,a,Rs=(40.,80.)):
    """Richardson(16:1) in R of P[mid,mid]."""
    a_arr=np.asarray(a,float); lo,mid,hi=np.argsort(a_arr)
    P1=connection_ip(eps,gam,a,Rs[0]); P2=connection_ip(eps,gam,a,Rs[1])
    Pext=(16*P2-P1)/15
    return Pext[mid,mid]

# ----------------------------------------------------------------- DA-order test
def monomials(X,deg):
    """All monomials of total degree<=deg in the columns of X (n_pts x n_vars)."""
    n,d=X.shape; cols=[]; names=[]
    for o in range(0,deg+1):
        for idx in combinations_with_replacement(range(d),o):
            c=np.ones(n)
            for j in idx: c=c*X[:,j]
            cols.append(c); names.append(idx)
    return np.array(cols).T, names

def da_residual(jet,deg):
    """Smallest singular value of the column-normalized monomial-jet matrix.
       jet: n_pts x n_vars (vars = t, f, f', ... f^(k)). Returns s_min/s_max."""
    M,_=monomials(jet,deg)
    norms=np.linalg.norm(M,axis=0); norms[norms==0]=1.0
    Mn=M/norms
    s=np.linalg.svd(Mn,compute_uv=False)
    return s[-1]/s[0], M.shape[1], jet.shape[0]

def cheb_jet(tnodes,fnodes,teval,kmax,degfit=None):
    """Chebyshev interpolant of f on [a,b]; returns jet [t,f,f',...,f^kmax] at teval."""
    a,b=tnodes.min(),tnodes.max()
    # map to [-1,1]
    xn=(2*tnodes-(a+b))/(b-a); xe=(2*teval-(a+b))/(b-a)
    deg=degfit if degfit else len(tnodes)-1
    coef=C.chebfit(xn,fnodes,deg)
    cols=[teval, C.chebval(xe,coef)]
    ck=coef.copy(); scale=2.0/(b-a)
    for _ in range(kmax):
        ck=C.chebder(ck)*scale
        cols.append(C.chebval(xe,ck))
    return np.column_stack(cols), coef

def run_block(name,tnodes,fnodes,degfit):
    a,b=tnodes.min(),tnodes.max()
    teval=np.linspace(a+0.06*(b-a),b-0.06*(b-a),120)
    print(f"\n--- {name} ---")
    jet2,coef=cheb_jet(tnodes,fnodes,teval,2,degfit=degfit)
    tail=np.abs(coef[-4:]).max()/ (np.abs(coef).max()+1e-300)
    print(f"  Chebyshev fit deg={degfit}: coeff tail/peak = {tail:.2e} (smoothness/noise floor)")
    for k in (1,2):
        jet=jet2[:,:k+2]   # t,f,...,f^k
        row=[]
        for d in (2,3,4):
            s,ncol,npt=da_residual(jet,d)
            row.append(f"d{d}:{s:.2e}")
        print(f"  order-{k} algebraic ODE residual  " + "  ".join(row))
    return

if __name__=="__main__":
    eps=(-2,0,3); gam=(1,0.8,1.2); alo,ahi=-1.0,2.0
    # accuracy gate at canonical a_mid=0.5
    s_can=sigma_of(eps,gam,(alo,0.5,ahi))
    s_can_hi=sigma_of(eps,gam,(alo,0.5,ahi),Rs=(80.,160.))
    print("="*78); print("ACCURACY GATE"); print("="*78)
    print(f"  sigma(a_mid=0.5) R(40,80)={s_can:.7f}  R(80,160)={s_can_hi:.7f}  "
          f"gold=0.214724  |diff|={abs(s_can-0.214724):.1e}  self={abs(s_can-s_can_hi):.1e}")

    # build the slice sigma(a_mid) on Chebyshev nodes
    Nn=33
    xn=np.cos(np.pi*np.arange(Nn)/(Nn-1))          # Chebyshev extrema in [-1,1]
    t0,t1=0.0,1.0
    tnodes=0.5*(t1-t0)*(xn[::-1]+1)+t0             # ascending in [0,1]
    print("\nbuilding sigma(a_mid) on %d Chebyshev nodes in [%.2f,%.2f] ..."%(Nn,t0,t1))
    snodes=np.array([sigma_of(eps,gam,(alo,t,ahi)) for t in tnodes])
    print("  sigma range [%.5f, %.5f]" % (snodes.min(),snodes.max()))
    # smoothness diagnostic: chebyshev coeff decay of sigma
    cfull=C.chebfit((2*tnodes-(t0+t1))/(t1-t0),snodes,Nn-1)
    print("  sigma chebyshev coeff |c_k|: first 6", np.round(np.abs(cfull[:6]),6),
          " last 4", np.round(np.abs(cfull[-4:]),9))

    print("\n"+"="*78); print("DIFFERENTIAL-ALGEBRAIC ORDER TEST  (s_min/s_max; small => ODE exists)")
    print("="*78)
    degfit=22   # truncate below the noise floor seen in the coeff tail
    # CONTROLS
    run_block("CONTROL exp  g=exp(-1.3(t-0.3)^2)+0.15  [Liouvillian: order-1 d2 expected hit]",
              tnodes, np.exp(-1.3*(tnodes-0.3)**2)+0.15, degfit)
    run_block("CONTROL J0   g=J0(9 t)  [holonomic: order-1 miss, order-2 hit]",
              tnodes, j0(9*tnodes), degfit)
    # SIGMA
    run_block("SIGMA  sigma(a_mid)  [claim: order-1 MISS at all d => non-Liouvillian]",
              tnodes, snodes, degfit)

    print("\nREADING: compare each row to the exp control's order-1 d2 entry (the Liouvillian")
    print("floor). sigma's order-1 residuals staying ORDERS ABOVE that floor = sigma is not")
    print("order-1 differentially algebraic on this slice = non-Liouvillian-classical, the")
    print("conclusion the non-rigidity->transcendence derivation targets.")
