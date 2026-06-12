"""ws_anti_chunk.py -- high-accuracy sigma slice for the antiderivative discriminator.
Usage: python3 ws_anti_chunk.py I0 I1   (computes Chebyshev nodes I0..I1-1, appends npz)"""
import sys,os,numpy as np
from scipy.integrate import solve_ivp
I0,I1=int(sys.argv[1]),int(sys.argv[2])
def build(eps,gam,a):
    eps=np.asarray(eps,float);gam=np.asarray(gam,float);a=np.asarray(a,float)
    H0=np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            if i!=j: H0[i,j]=gam[i]*gam[j]*(a[i]-a[j])/(eps[i]-eps[j])
        H0[i,i]=-sum(gam[k]**2*(a[i]-a[k])/(eps[i]-eps[k]) for k in range(3) if k!=i)
    return H0,a
def Pmid(eps,gam,a,Rs=(40.,80.),rtol=1e-12,atol=1e-13):
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
eps=(-2,0,3); gam0=np.array([1.,0.8,1.2]); a3=(-1.,0.5,2.)
Nn=35; xn=np.cos(np.pi*np.arange(Nn)/(Nn-1)); t0,t1=0.45,2.3
tn=0.5*(t1-t0)*(xn[::-1]+1)+t0
fn="anti_slice_hi.npz"
vals=np.full(Nn,np.nan); 
if os.path.exists(fn): vals=np.load(fn)['vals']
for k in range(I0,I1):
    vals[k]=Pmid(eps,tn[k]*gam0,a3)
    np.savez(fn,tn=tn,vals=vals)
    print(f"node {k}: lam={tn[k]:.4f} sigma={vals[k]:.8f}",flush=True)
print("chunk done:",I0,I1)
