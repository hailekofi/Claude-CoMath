"""ws_anti_verdict.py -- final antiderivative discriminator on the high-accuracy slice
(anti_slice_hi.npz, sigma to ~1e-7). Two tests, controls noise-matched at 1e-7:
 [A] DA order-1 battery on f and f' (as before, floors now ~10x lower);
 [B] RATIO test: h := f''/f' rational in t of low degree?  (f = int A e^{int r}, r rational
     <=> h rational; erf: h linear). Better conditioned: linear-in-h monomial system
     {t^k, h t^k}: smin small <=> h rational(deg).
Reproduce after chunks: python3 ws_anti_verdict.py"""
import numpy as np
from numpy.polynomial import chebyshev as C
from itertools import combinations_with_replacement
from scipy.special import j0,erf
rng=np.random.default_rng(7)

d=np.load("anti_slice_hi.npz"); tn,sn=d['tn'],d['vals']
assert not np.any(np.isnan(sn)), "slice incomplete"
print(f"slice: {len(tn)} nodes, sigma in [{sn.min():.6f},{sn.max():.6f}]")
NOISE=1e-7
yerf=0.5*erf(1.2*(tn-1.2))+0.4+rng.standard_normal(len(tn))*NOISE
yj0=j0(9*tn)+rng.standard_normal(len(tn))*NOISE

def monomials(X,deg):
    n,dd=X.shape; cols=[]
    for o in range(0,deg+1):
        for idx in combinations_with_replacement(range(dd),o):
            c=np.ones(n)
            for j in idx: c=c*X[:,j]
            cols.append(c)
    return np.array(cols).T
def smin(jet,deg):
    M=monomials(jet,deg); nm=np.linalg.norm(M,axis=0); nm[nm==0]=1
    return np.linalg.svd(M/nm,compute_uv=False)[-1]/np.linalg.svd(M/nm,compute_uv=False)[0]

def jets(tn,fn,degfit=26,kmax=3):
    a,b=tn.min(),tn.max(); xn=(2*tn-(a+b))/(b-a)
    coef=C.chebfit(xn,fn,degfit)
    te=np.linspace(a+0.06*(b-a),b-0.06*(b-a),120)
    xe=(2*te-(a+b))/(b-a); sc=2.0/(b-a)
    out=[te,C.chebval(xe,coef)]
    ck=coef.copy()
    for _ in range(kmax):
        ck=C.chebder(ck)*sc; out.append(C.chebval(xe,ck))
    return out  # te,f,f1,f2,f3

print("\n[A] DA order-1 battery (d2/d3):")
print("-"*74)
for name,fn in (("erf-type f",yerf),("J0 f",yj0),("SIGMA f",sn)):
    te,f,f1,f2,f3=jets(tn,fn)
    r0=[smin(np.column_stack([te,f,f1]),dd) for dd in (2,3)]
    r1=[smin(np.column_stack([te,f1,f2]),dd) for dd in (2,3)]
    print(f"  {name:12s}: f  d2 {r0[0]:.2e} d3 {r0[1]:.2e}   |  f' d2 {r1[0]:.2e} d3 {r1[1]:.2e}")
print("-"*74)
print("\n[B] RATIO test: h=f''/f' rational of degree m?  (smin of {t^k, h t^k}, k<=m)")
print("-"*74)
for name,fn in (("erf-type",yerf),("J0",yj0),("SIGMA",sn)):
    te,f,f1,f2,f3=jets(tn,fn)
    h=f2/f1
    ok=np.abs(f1)>0.02*np.max(np.abs(f1))   # avoid f'~0 points
    row=[]
    for m in (1,2,3,4):
        cols=[te[ok]**k for k in range(m+1)]+[h[ok]*te[ok]**k for k in range(m+1)]
        M=np.array(cols).T; nm=np.linalg.norm(M,axis=0); nm[nm==0]=1
        s=np.linalg.svd(M/nm,compute_uv=False); row.append(s[-1]/s[0])
    print(f"  {name:10s}: deg1 {row[0]:.2e}  deg2 {row[1]:.2e}  deg3 {row[2]:.2e}  deg4 {row[3]:.2e}")
print("-"*74)
print("READING: erf row = HIT calibration (its h is exactly linear). J0 row = MISS calibration.")
print("SIGMA at the J0 level at all degrees => sigma is NOT an antiderivative of an")
print("elementary (exp-closure) function of the parameter on this slice.")
