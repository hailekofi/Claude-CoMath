"""
ws_construction_generic.py -- is the S=U.Delta.L construction (elementary Delta + slot
equations) a TYPE-1 fact, or generic N=3 multistate-Landau-Zener structure?

Reflection driver: a reviewer asks whether S=U.Delta.L is "just the classic LU result"
and whether U,L are anchored in {gamma,eps,a}. To answer the DEEP version we test whether
the content-bearing parts -- (i) Delta purely elementary = e^{-/+ pi |b_j|} with
b_j = sum_k H0_jk^2/(a_j-a_k), and (ii) the slot equations
sigma=|Delta_m+U12 Delta_h L21|^2, b=|L20|^2|Delta_h|^2 -- hold for an ARBITRARY real
symmetric H0 (same slopes), not only the Type-1 Cauchy H0.

If they hold for generic H0, the construction is GENERIC N=3 MLZ sectorial-Stokes
structure; Type-1 (integrability/Cauchy H0) is NOT in the construction -- it lives in the
*value* of the transcendental multipliers and in the integrable spectral geometry.

Reproduce: python3 ws_construction_generic.py   (numpy/scipy)
"""
import numpy as np
from scipy.integrate import solve_ivp

rng=np.random.default_rng(3)

def type1_H0(eps,gam,a):
    eps=np.asarray(eps,float);gam=np.asarray(gam,float);a=np.asarray(a,float)
    H0=np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            if i!=j: H0[i,j]=gam[i]*gam[j]*(a[i]-a[j])/(eps[i]-eps[j])
        H0[i,i]=-sum(gam[k]**2*(a[i]-a[k])/(eps[i]-eps[k]) for k in range(3) if k!=i)
    return H0

def formal_frame(H0,a,u):
    a=np.asarray(a,float)
    b1=np.array([sum(H0[j,k]**2/(a[j]-a[k]) for k in range(3) if k!=j) for j in range(3)])
    T1=np.zeros((3,3),complex)
    for k in range(3):
        for j in range(3):
            if k!=j: T1[k,j]=H0[k,j]/(a[j]-a[k])
    ph=np.array([-1j*(a[j]*u**2/2+H0[j,j]*u+b1[j]*np.log(complex(u))) for j in range(3)])
    return (np.eye(3)+T1/u)@np.diag(np.exp(ph)), b1

def connection_S(H0,a,R,rtol=1e-12,atol=1e-13):
    a=np.asarray(a,float)
    Y0,b1=formal_frame(H0,a,-R)
    sol=solve_ivp(lambda u,Y:(-1j*(H0+u*np.diag(a))@Y.reshape(3,3)).reshape(-1),
                  [-R,R],Y0.reshape(-1).astype(complex),rtol=rtol,atol=atol,method='DOP853')
    YR=sol.y[:,-1].reshape(3,3)
    Yf,_=formal_frame(H0,a,R)
    C=np.linalg.solve(Yf,YR)
    lo,mid,hi=np.argsort(a); perm=[lo,mid,hi]
    C=C[np.ix_(perm,perm)]; b1o=b1[perm]
    S=C@np.diag(np.exp(-np.pi*b1o))
    return S,b1o

def gauss_LDU(C):
    L=np.eye(3,dtype=complex);U=np.eye(3,dtype=complex);D=np.zeros(3,complex)
    D[0]=C[0,0];L[1,0]=C[1,0]/D[0];L[2,0]=C[2,0]/D[0]
    U[0,1]=C[0,1]/D[0];U[0,2]=C[0,2]/D[0]
    D[1]=C[1,1]-L[1,0]*D[0]*U[0,1]
    L[2,1]=(C[2,1]-L[2,0]*D[0]*U[0,1])/D[1]
    U[1,2]=(C[1,2]-L[1,0]*D[0]*U[0,2])/D[1]
    D[2]=C[2,2]-L[2,0]*D[0]*U[0,2]-L[2,1]*D[1]*U[1,2]
    return L,D,U

def gauss_UDL(C):
    J=np.eye(3)[::-1]; L,D,U=gauss_LDU(J@C@J)
    return J@L@J,(J@np.diag(D)@J).diagonal().copy(),J@U@J  # U,D,L

def assess(name,H0,a):
    a=np.asarray(a,float); lo,mid,hi=np.argsort(a)
    S1,b1=connection_S(H0,a,60.); S2,_=connection_S(H0,a,120.)
    S=S2
    BE=np.exp(-np.pi*np.abs(b1))
    U,D,L=gauss_UDL(S)
    rat=np.abs(D)/BE                      # H* says (1,1,1)
    P=np.abs(S)**2
    sig=np.abs(D[1]+U[1,2]*D[2]*L[2,1])**2
    bb =np.abs(L[2,0])**2*np.abs(D[2])**2
    print(f"\n=== {name} ===")
    print(f"  symmetric H0? {np.allclose(H0,H0.T)};  |det S|={abs(np.linalg.det(S)):.4f}; "
          f"unitarity |SS^H-I|={np.max(np.abs(S@S.conj().T-np.eye(3))):.1e}")
    print(f"  Delta-elementarity |D_j|/e^(-/+pi|b_j|) = {np.round(rat,5)}   (want 1,1,1)")
    print(f"  slot sigma=|D_m+U12 D_h L21|^2 = {sig:.6f}  vs P_mm={P[1,1]:.6f}  (d={abs(sig-P[1,1]):.1e})")
    print(f"  slot b=|L20|^2|D_h|^2          = {bb:.6f}  vs P_hl={P[2,0]:.6f}  (d={abs(bb-P[2,0]):.1e})")
    print(f"  P_ll/BE_lo^2={P[0,0]/BE[0]**2:.5f}  P_hh/BE_hi^2={P[2,2]/BE[2]**2:.5f}  (BE pinning)")
    print(f"  the 6 multiplier moduli (UNANCHORED Stokes data): "
          f"|U01|={abs(U[0,1]):.3f} |U02|={abs(U[0,2]):.3f} |U12|={abs(U[1,2]):.3f} "
          f"|L10|={abs(L[1,0]):.3f} |L20|={abs(L[2,0]):.3f} |L21|={abs(L[2,1]):.3f}")
    return P[1,1]

if __name__=="__main__":
    a=[-1.0,0.5,2.0]
    print("Does the construction (elementary Delta + slot eqns) need Type-1, or is it generic?")
    # Type-1 Cauchy H0
    H0_t1=type1_H0([-2,0,3],[1,0.8,1.2],a)
    assess("TYPE-1 (Cauchy H0)  -- canonical anchor",H0_t1,a)
    # generic real symmetric H0, same slopes, NOT Type-1
    for tag in ("generic-1","generic-2"):
        M=rng.standard_normal((3,3)); H0g=0.6*(M+M.T)/2     # random real symmetric, NOT Cauchy
        assess(f"GENERIC symmetric H0 ({tag}) -- NOT Type-1",H0g,a)
    print("\nIf GENERIC H0 also gives Delta-elementary + exact slot eqns => the construction")
    print("is general N=3 MLZ sectorial structure, NOT a Type-1 fact. Type-1 content (integrability,")
    print("genus-0 spectral curve) is NOT in the skeleton; it lives in the multiplier VALUES.")
