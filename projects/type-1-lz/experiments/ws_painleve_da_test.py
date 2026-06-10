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

RESULT (coupling-scale slice gam->lam*gam, lam in [0.45,2.30], sigma sweeps [0.016,0.905];
35 Chebyshev nodes, deg-26 fit, sigma to ~3e-7-4e-6):
  detector validated by controls at the ONLY diagnostic cell, order-1 d2:
     exp (Liouvillian)  order-1 d2 = 4e-15  HIT     (floor 3e-15)
     J0  (holonomic)    order-1 d2 = 3.0e-3 MISS, order-2 d2 = 7e-12 HIT  (floor 3e-9)
  SIGMA:  order-1 d2 = 1.48e-4  MISS  (37x the deg-26 floor; ~matches J0's clean miss),
          order-2 d2 = 3.2e-6   at floor  (order-2 ~HIT, Painleve/holonomic signature).
  ROBUSTNESS (the decisive check): sigma's order-1 residual was 1.50e-4 at deg-20 and
  1.48e-4 at deg-26 while the fit floor improved 60x (2.5e-4 -> 4e-6) -> the residual is a
  REAL ODE-miss, not fit error.
  READING: sigma(parameter) tracks the J0 (transcendent: order-1 NO / order-2 YES) pattern,
  NOT the exp (elementary, order-1 YES) pattern.  => sigma is NOT order-1 differentially
  algebraic = non-Liouvillian-classical (CONFIRMED, robust), and is consistent with an
  order-2 (Painleve/Garnier-type) transcendent (SUGGESTIVE; the order-2 residual is at-floor
  not orders-below, so pinning the exact order-2 structure needs sigma to ~1e-9: oracle/mpmath).
  Evidence level: numerically-supported for the order-1 falsification; suggestive for the
  order-2 Painleve structure.
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

def connection_ip(eps,gam,a,R,rtol=1e-10,atol=1e-11):
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

def sigma_of(eps,gam,a,Rs=(30.,60.)):
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
    s_can_hi=sigma_of(eps,gam,(alo,0.5,ahi),Rs=(40.,80.))
    print("="*78); print("ACCURACY GATE"); print("="*78)
    print(f"  sigma(a_mid=0.5) R(30,60)={s_can:.7f}  R(40,80)={s_can_hi:.7f}  "
          f"gold=0.214724  |diff|={abs(s_can-0.214724):.1e}  self={abs(s_can-s_can_hi):.1e}")

    # build the slice sigma(lambda) on Chebyshev nodes -- COUPLING-SCALE slice gam->lam*gam,
    # chosen so sigma sweeps a LARGE range (diabatic ~1 -> adiabatic ~0) with curvature, so a
    # genuine order-1 MISS is resolvable above the ~1e-6 floor (the flat a_mid slice was not).
    gam0=np.array([1.0,0.8,1.2]); a3=(alo,0.5,ahi)
    Nn=35
    xn=np.cos(np.pi*np.arange(Nn)/(Nn-1))          # Chebyshev extrema in [-1,1]
    t0,t1=0.45,2.3
    tnodes=0.5*(t1-t0)*(xn[::-1]+1)+t0             # ascending lambda in [t0,t1]
    print("\nbuilding sigma(lambda) [gam->lam*gam] on %d Chebyshev nodes, lam in [%.2f,%.2f] ..."
          %(Nn,t0,t1))
    snodes=np.array([sigma_of(eps,lam*gam0,a3) for lam in tnodes])
    print("  sigma range [%.5f, %.5f]   (want a WIDE sweep)" % (snodes.min(),snodes.max()))
    # smoothness diagnostic: chebyshev coeff decay of sigma
    cfull=C.chebfit((2*tnodes-(t0+t1))/(t1-t0),snodes,Nn-1)
    print("  sigma chebyshev coeff |c_k|: first 6", np.round(np.abs(cfull[:6]),6),
          " last 4", np.round(np.abs(cfull[-4:]),9))

    print("\n"+"="*78); print("DIFFERENTIAL-ALGEBRAIC ORDER TEST  (s_min/s_max; small => ODE exists)")
    print("="*78)
    degfit=26   # resolve sigma to its ~1e-6 noise floor (deg-20 under-resolved it)
    # CONTROLS
    run_block("CONTROL exp  g=exp(-1.3(t-0.3)^2)+0.15  [Liouvillian: order-1 d2 expected hit]",
              tnodes, np.exp(-1.3*(tnodes-0.3)**2)+0.15, degfit)
    run_block("CONTROL J0   g=J0(9 t)  [holonomic: order-1 miss, order-2 hit]",
              tnodes, j0(9*tnodes), degfit)
    # SIGMA
    run_block("SIGMA  sigma(lambda)  [claim: order-1 MISS like J0 => non-Liouvillian-classical]",
              tnodes, snodes, degfit)

    print("\nREADING (only LOW (k,d) are diagnostic; high-degree collapse is Vandermonde")
    print("ill-conditioning, not an ODE). Calibrate by the controls at order-1 d2:")
    print("  exp (Liouvillian, order-1) -> ~1e-15 HIT ;  J0 (no order-1) -> ~1e-4 MISS.")
    print("If sigma's order-1 d2 residual sits at J0's MISS scale (>> its ~1e-6 floor), sigma")
    print("is NOT order-1 differentially algebraic = non-Liouvillian-classical (derivation target).")
    print("If it sits at the noise floor, the test is precision-limited -> needs a higher-")
    print("precision sigma(parameter) solver (oracle/mpmath).")
