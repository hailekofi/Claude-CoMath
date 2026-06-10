"""
ws_stokes_slots2.py -- convergence-hardened test of the SHARP hypothesis from
ws_stokes_slots.py (blind probe v1):

  H*: The physical connection matrix in the slope-ordered formal frame factorizes as
        S = U . Delta . L          (upper-lateral sectorial product)
      with U unit-upper, L unit-lower (slope order), and Delta PURELY FORMAL:
        |Delta_j| = e^{-pi |b1_j|},   b1_j = sum_k H0_jk^2/(a_j-a_k).
      All transcendence then sits in the 6 unipotent multipliers, and
        P_ll = |Delta_l|^2,  P_hh = |Delta_h|^2          (BE, multiplier-free)
        b    = P[hi,lo] = |L20|^2 |Delta_h|^2            (ONE multiplier)
        sigma= |Delta_m + U12 Delta_h L21|^2             (TWO multipliers + rel. phase)

Hardening vs v1: include the O(1/u) formal prefactor T1 (T1_kj = H0_kj/(a_j-a_k), k!=j)
so frame errors drop 1/R -> 1/R^2; normalize S = C diag(e^{-pi b1}); test R-ladder.
Also: unitarity probe on S; gold gate on P_mm.

Context: wildness fixes WHERE the data lives (sectors/slots/BE-pinning); whether the
slot VALUES are transcendental is a separate (rigidity) question -- the two axes are
deliberately disaggregated in the memo this script feeds.
"""
import numpy as np
from scipy.integrate import solve_ivp

def build_H0(eps,gam,a):
    eps=np.asarray(eps,float); gam=np.asarray(gam,float); a=np.asarray(a,float)
    H0=np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            if i!=j: H0[i,j]=gam[i]*gam[j]*(a[i]-a[j])/(eps[i]-eps[j])
        H0[i,i]=-sum(gam[k]**2*(a[i]-a[k])/(eps[i]-eps[k]) for k in range(3) if k!=i)
    return H0,a

def formal_frame(H0,a,u):
    """Yf = (I + T1/u) diag(e^{Phi_j}), Phi_j = -i[a_j u^2/2 + H0_jj u + b1_j log u]."""
    a=np.asarray(a,float)
    b1=np.array([sum(H0[j,k]**2/(a[j]-a[k]) for k in range(3) if k!=j) for j in range(3)])
    T1=np.zeros((3,3),dtype=complex)
    for k in range(3):
        for j in range(3):
            if k!=j: T1[k,j]=H0[k,j]/(a[j]-a[k])
    phase=np.array([-1j*(a[j]*u**2/2 + H0[j,j]*u + b1[j]*np.log(complex(u))) for j in range(3)])
    return (np.eye(3)+T1/u)@np.diag(np.exp(phase)), b1

def connection_S(eps,gam,a,R,rtol=1e-12,atol=1e-13):
    """Slope-ordered, physically normalized S; |S|^2 should be doubly stochastic."""
    H0,a=build_H0(eps,gam,a); a=np.asarray(a,float)
    Y0,b1=formal_frame(H0,a,-R)
    sol=solve_ivp(lambda u,Y:(-1j*(H0+u*np.diag(a))@Y.reshape(3,3)).reshape(-1),
                  [-R,R],Y0.reshape(-1).astype(complex),rtol=rtol,atol=atol,method='DOP853')
    YR=sol.y[:,-1].reshape(3,3)
    Yf,_=formal_frame(H0,a,R)
    C=np.linalg.solve(Yf,YR)
    lo,mid,hi=np.argsort(a); perm=[lo,mid,hi]
    C=C[np.ix_(perm,perm)]; b1o=b1[perm]
    S=C@np.diag(np.exp(-np.pi*b1o))      # strip incoming half-monodromy moduli
    return S,b1o

def gauss_LDU(C):
    L=np.eye(3,dtype=complex); U=np.eye(3,dtype=complex); D=np.zeros(3,dtype=complex)
    D[0]=C[0,0]; L[1,0]=C[1,0]/D[0]; L[2,0]=C[2,0]/D[0]
    U[0,1]=C[0,1]/D[0]; U[0,2]=C[0,2]/D[0]
    D[1]=C[1,1]-L[1,0]*D[0]*U[0,1]
    L[2,1]=(C[2,1]-L[2,0]*D[0]*U[0,1])/D[1]
    U[1,2]=(C[1,2]-L[1,0]*D[0]*U[0,2])/D[1]
    D[2]=C[2,2]-L[2,0]*D[0]*U[0,2]-L[2,1]*D[1]*U[1,2]
    return L,D,U

def gauss_UDL(C):
    J=np.eye(3)[::-1]
    Lr,Dr,Ur=gauss_LDU(J@C@J)
    return J@Lr@J,(J@np.diag(Dr)@J).diagonal().copy(),J@Ur@J   # U, D, L

def run(nm,eps,gam,a,gold=None):
    print("="*88); print(f"CASE {nm}: eps={eps} gam={gam} a={a}   gold={gold}"); print("="*88)
    out={}
    for R in (60.0,120.0):
        S,b1o=connection_S(eps,gam,a,R)
        BE=np.exp(-np.pi*np.abs(b1o))
        P=np.abs(S)**2
        U,D,L=gauss_UDL(S)
        ratios=np.abs(D)/BE
        out[R]=dict(P=P,D=np.abs(D),ratios=ratios,U=U,L=L,S=S,BE=BE)
        print(f" R={R:5.0f}  rowsum {np.round(P.sum(1),7)} colsum {np.round(P.sum(0),7)}")
        print(f"          P_mm={P[1,1]:.7f}  P_ll={P[0,0]:.7f} (BE^2 {BE[0]**2:.7f})  "
              f"P_hh={P[2,2]:.7f} (BE^2 {BE[2]**2:.7f})")
        print(f"          UDL |D|/BE = {np.round(ratios,7)}")
    # Richardson (error ~ 1/R^2 now): Q_ext = (4 Q(2R) - Q(R))/3
    r_ext=(4*out[120.]['ratios']-out[60.]['ratios'])/3
    P_ext=(4*out[120.]['P']-out[60.]['P'])/3
    print(f"\n  EXTRAP  |D|/BE        = {np.round(r_ext,7)}   <-- H* says (1,1,1)")
    print(f"  EXTRAP  P_mm          = {P_ext[1,1]:.7f}   gold {gold}")
    print(f"  EXTRAP  P_ll/BE^2     = {P_ext[0,0]/out[120.]['BE'][0]**2:.7f},  "
          f"P_hh/BE^2 = {P_ext[2,2]/out[120.]['BE'][2]**2:.7f}   <-- BE exactness")
    # slot equations at R=120
    U,L,S=out[120.]['U'],out[120.]['L'],out[120.]['S']
    Dc=gauss_UDL(S)[1]
    sig=np.abs(Dc[1]+U[1,2]*Dc[2]*L[2,1])**2
    print(f"  sigma slot |D_m + U12 D_h L21|^2 = {sig:.7f}  (=P_mm raw {np.abs(S[1,1])**2:.7f})")
    print(f"  b slot     |L20|^2 |D_h|^2       = {np.abs(L[2,0])**2*np.abs(Dc[2])**2:.7f}"
          f"  (=P_hl raw {np.abs(S[2,0])**2:.7f})")
    print(f"  multipliers: |L10|={abs(L[1,0]):.5f} |L20|={abs(L[2,0]):.5f} |L21|={abs(L[2,1]):.5f}"
          f" |U01|={abs(U[0,1]):.5f} |U02|={abs(U[0,2]):.5f} |U12|={abs(U[1,2]):.5f}")
    # unitarity probe
    G=S@S.conj().T
    print(f"  unitarity:  max|S S^dag - I| = {np.max(np.abs(G-np.eye(3))):.3e}")
    return out

if __name__=="__main__":
    run("canonical",(-2,0,3),(1,0.8,1.2),(-1,0.5,2.0),gold=0.214724)
    run("sampleB",(-1,0,1.5),(0.9,1.1,0.8),(-0.7,0.4,1.3),gold=0.021018)
    run("unsorted-a",(-1.6,0.3,2.1),(1.05,-0.75,0.95),(1.4,-0.9,0.6))
