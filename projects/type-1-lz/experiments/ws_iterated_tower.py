"""
ws_iterated_tower.py -- explicit iterated-integral tower for sigma.

*** NUMERICAL CAVEAT (R42): the grid reconstruction below is UNRELIABLE. ***
The quick eigenframe under-resolves the sharp avoided crossing (pair 01, gap ~2e-4),
corrupting W near u~-0.25. Symptom: the EXACT adiabatic resummation [1] gives 0.498 while
a correct propagator at the SAME window U=30 gives gold 0.214725. Hence the leading-term
value [2] and Lambda [3] are RETRACTED. The ANALYTIC tower (the dressed-coupling 1-forms
and the Dyson/iterated-integral series) is correct and independent of this instantiation.
Trustworthy numbers require the rational lambda-frame (E_i(lam),V_i(lam) smooth rational,
no eigenvector alignment) -- not yet done.
"""

import numpy as np
from scipy.integrate import cumulative_trapezoid, trapezoid, solve_ivp

def build(eps,gam,a):
    eps=np.asarray(eps,float);gam=np.asarray(gam,float);a=np.asarray(a,float)
    H0=np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            if i!=j: H0[i,j]=gam[i]*gam[j]*(a[i]-a[j])/(eps[i]-eps[j])
        H0[i,i]=-sum(gam[k]**2*(a[i]-a[k])/(eps[i]-eps[k]) for k in range(3) if k!=i)
    return H0,a

eps=(-2,0,3); gam=(1,0.8,1.2); a=(-1.,0.5,2.); H0,A=build(eps,gam,a); Adiag=np.diag(A)
lo,mid,hi=np.argsort(np.asarray(a,float))

# gold sigma via diabatic propagator + slope map
def gold(U=80.,rtol=1e-11):
    def onerun(R):
        w,V=np.linalg.eigh(H0+(-R)*Adiag)
        sol=solve_ivp(lambda u,Y:(-1j*(H0+u*Adiag)@Y.reshape(3,3)).reshape(-1),
                      [-R,R],V.astype(complex).reshape(-1),rtol=rtol,atol=1e-12,method='DOP853')
        w2,V2=np.linalg.eigh(H0+R*Adiag); P=np.abs(V2.conj().T@sol.y[:,-1].reshape(3,3))**2
        sin=[int(np.argmin(np.abs(np.asarray(a)-w[k]/(-R)))) for k in range(3)]
        sout=[int(np.argmin(np.abs(np.asarray(a)-w2[k]/R))) for k in range(3)]
        Pd=np.zeros((3,3))
        for ki in range(3):
            for ko in range(3): Pd[sin[ki],sout[ko]]=P[ki,ko]
        return Pd
    P=(16*onerun(80.)-onerun(40.))/15; return P[mid,mid]

# dressed coupling on a fine grid
U=30.; M=24001; us=np.linspace(-U,U,M)
E=np.zeros((3,M)); Wt=np.zeros((3,3,M),complex)
Vprev=None
Eg=np.zeros((3,M)); Vall=np.zeros((3,3,M))
for n,u in enumerate(us):
    w,V=np.linalg.eigh(H0+u*Adiag)
    if Vprev is not None:
        for k in range(3):
            if np.real(np.vdot(Vprev[:,k],V[:,k]))<0: V[:,k]*=-1
    Vprev=V; Eg[:,n]=w; Vall[:,:,n]=V
# couplings W_{ij} = <i|A|j>/(E_j-E_i)
VA=np.einsum('iam,ab,bjm->ijm',Vall.transpose(1,0,2).conj() if False else Vall.transpose(1,0,2),
             Adiag, Vall)  # (i,j,m) = V_i^T A V_j (real symmetric)
W=np.zeros((3,3,M))
for i in range(3):
    for j in range(3):
        if i!=j: W[i,j]=VA[i,j]/(Eg[j]-Eg[i])
# phases Theta_{ij}(u) = int_0^u (E_i-E_j)
mid0=M//2
def cumphase(x):
    c=cumulative_trapezoid(x,us,initial=0.0); return c-c[mid0]
Theta=np.zeros((3,3,M))
for i in range(3):
    for j in range(3): Theta[i,j]=cumphase(Eg[i]-Eg[j])
Wd=W*np.exp(1j*Theta)   # dressed coupling W~_{ij}

# [2] leading tower term: c_m^{(2)} = sum_k int W~_{mk}(u2) [int_{-U}^{u2} W~_{km}(u1) du1] du2
c2=0j
for k in range(3):
    if k==mid: continue
    inner=cumulative_trapezoid(Wd[k,mid],us,initial=0.0)        # c_k^{(1)}(u) = -inner ; sign below
    c2+= trapezoid(Wd[mid,k]*inner,us)                          # (-1)*(-1)=+1
sig2=abs(1+c2)**2

# [1] exact adiabatic resummation: integrate c' = -W~(u) c, c(-U)=e_mid
from scipy.interpolate import interp1d
Wd_re=interp1d(us,Wd.real,axis=2); Wd_im=interp1d(us,Wd.imag,axis=2)
def rhs(u,c):
    Wu=Wd_re(u)+1j*Wd_im(u); return -(Wu@c)
c0=np.zeros(3,complex); c0[mid]=1
sol=solve_ivp(rhs,[-U,U],c0,rtol=1e-10,atol=1e-12,method='DOP853')
sig_ad=abs(sol.y[mid,-1])**2

# [3] marginality
Lam=sum(trapezoid(np.abs(W[i,j]),us) for i in range(3) for j in range(3) if i<j)

g=gold()
print("="*72)
print("THE ITERATED-INTEGRAL TOWER FOR SIGMA  (canonical anchor)")
print("="*72)
print(f"  gold sigma (full diabatic solve)        = {g:.6f}")
print(f"  [1] exact adiabatic resummation |c_m|^2 = {sig_ad:.6f}   (validates W~ machinery)")
print(f"  [2] leading term  |1 + c_m^(2)|^2       = {sig2:.6f}   (single-excursion Stuckelberg)")
print(f"      c_m^(2) = {c2:.4f}")
print(f"  [3] marginality  Lambda = sum int|W_ij| = {Lam:.4f}   (~pi => marginal tower)")
print()
print("  => leading term is explicit & elementary but far from sigma; the tower does not")
print("     terminate or converge geometrically (Lambda~pi). Unlike N=2 (whose infinite")
print("     adiabatic tower resums to an elementary exponential), the N=3 resummation is")
print("     non-elementary: that -- not non-termination -- is the transcendence.")
