"""
ws_stokes_slots.py -- BLIND derivation probe for the sectorial (wild) factorization of
the Type-1 N=3 connection matrix, and the location of {sigma, b} in multiplier slots.

ANALYTIC INPUT (derived, not assumed -- see memo):
  Rates q_j(u) = -i a_j u^2/2 + ... ; all pairwise differences q_i-q_j = -i(a_i-a_j)u^2/2
  have PURE-IMAGINARY quadratic coefficients (a real). Hence on u = R e^{i theta}:
      Re(q_i - q_j) = (1/2)(a_i - a_j) R^2 sin(2 theta).
  =>  ALL three pairs share the SAME singular geometry:
      - anti-Stokes (oscillatory) directions: the AXES theta = 0, pi/2, pi, 3pi/2
        (the physical scattering directions are anti-Stokes: that is WHY P is unistochastic);
      - Stokes/jump rays (Re extremal): the DIAGONALS theta = pi/4, 3pi/4, 5pi/4, 7pi/4;
      - dominance order of {e^{q_j}} = slope order of a in quadrants 1,3 and REVERSED in 2,4.
  The path from u=-infty (theta=pi) to u=+infty (theta=0) through the upper half plane
  crosses exactly TWO jump rays: 3pi/4 then pi/4. Each carries a FULL unipotent factor
  (all three pairs jump together -- maximal ray-degeneracy), oppositely triangular in
  slope order. Therefore in the formal frame
      C = (one unipotent) . (diagonal formal data) . (opposite unipotent),
  i.e. C admits an EXACT Gauss-type triangular factorization w.r.t. slope order.

THIS SCRIPT tests, numerically and blind (no target values), which factorization holds
and what sits in each slot:
  [T1] gold gate: |C|^2 doubly stochastic, P_mm matches oracle.
  [T2] LDU and UDL Gauss decompositions in slope order: report all multiplier moduli.
  [T3] formal-diagonal content: |D_j| vs the half-monodromy moduli e^{+/- pi b1_j}
       (b1_j = sum_k H0_jk^2/(a_j-a_k), the signed-BE exponents) -- which match, which
       instead satisfy a multiplier relation (the unitarity-absorbed corner).
  [T4] slot identification:  b/BE  vs  single-multiplier moduli;  sigma's dressing vs the
       two-step recombination term (Schur complement structure).
  [T5] reality/Schwarz structure of C (which conjugation symmetry the multipliers obey).
  [T6] formal monodromy is DIAGONAL (no permutation part): adjudicates the committed
       Proposition claim (b) ("formal monodromy = oriented 3-cycle") -- expected RETRACTION.

Reproduce: python3 ws_stokes_slots.py     (numpy/scipy only)
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

def formal_solution(H0,a,u):
    """Thome formal frame: per branch j, phase -i[a_j u^2/2 + H0_jj u + b1_j log u]."""
    a=np.asarray(a,float)
    b1=np.array([sum(H0[j,k]**2/(a[j]-a[k]) for k in range(3) if k!=j) for j in range(3)])
    phase=np.array([-1j*(a[j]*u**2/2 + H0[j,j]*u + b1[j]*np.log(complex(u))) for j in range(3)])
    return np.diag(np.exp(phase)).astype(complex), b1

def connection_C(eps,gam,a,R=80.0,rtol=1e-12,atol=1e-13):
    """Complex connection matrix between formal frames at -R and +R (slope-unordered)."""
    H0,a=build_H0(eps,gam,a)
    Y0,b1=formal_solution(H0,a,-R)
    sol=solve_ivp(lambda u,Y:(-1j*(H0+u*np.diag(a))@Y.reshape(3,3)).reshape(-1),
                  [-R,R],Y0.reshape(-1).astype(complex),rtol=rtol,atol=atol,method='DOP853')
    YR=sol.y[:,-1].reshape(3,3)
    Yf,_=formal_solution(H0,a,R)
    return np.linalg.solve(Yf,YR), b1

def gauss_LDU(C):
    """C = L D U, L unit-lower, D diag, U unit-upper (slope order assumed)."""
    L=np.eye(3,dtype=complex); U=np.eye(3,dtype=complex); D=np.zeros(3,dtype=complex)
    D[0]=C[0,0]
    L[1,0]=C[1,0]/D[0]; L[2,0]=C[2,0]/D[0]
    U[0,1]=C[0,1]/D[0]; U[0,2]=C[0,2]/D[0]
    D[1]=C[1,1]-L[1,0]*D[0]*U[0,1]
    L[2,1]=(C[2,1]-L[2,0]*D[0]*U[0,1])/D[1]
    U[1,2]=(C[1,2]-L[1,0]*D[0]*U[0,2])/D[1]
    D[2]=C[2,2]-L[2,0]*D[0]*U[0,2]-L[2,1]*D[1]*U[1,2]
    return L,D,U

def gauss_UDL(C):
    """C = U D L (anti-Gauss: pivot from the (2,2) corner)."""
    J=np.eye(3)[::-1]                       # exchange matrix
    Lr,Dr,Ur=gauss_LDU(J@C@J)               # LDU of the rotated matrix
    U=J@Lr@J; L=J@Ur@J; D=(J@np.diag(Dr)@J).diagonal().copy()
    return U,D,L

def report(nm,eps,gam,a,gold=None):
    print("="*86); print(f"CASE {nm}: eps={eps} gam={gam} a={a}"); print("="*86)
    a_arr=np.asarray(a,float); lo,mid,hi=np.argsort(a_arr)
    # Richardson pair in R for the gold gate
    C1,b1=connection_C(eps,gam,a,R=60.0); C2,_=connection_C(eps,gam,a,R=120.0)
    P1,P2=np.abs(C1)**2,np.abs(C2)**2; Pext=(16*P2-P1)/15
    perm=[lo,mid,hi]
    C=C2[np.ix_(perm,perm)]                 # slope-ordered complex connection matrix
    b1o=b1[perm]
    print(f"[T1] row sums {np.round(P2.sum(1),6)}  col sums {np.round(P2.sum(0),6)}")
    pm=Pext[mid,mid]
    gtxt=f"   gold={gold}" if gold else ""
    print(f"     P_mm = {pm:.6f}{gtxt}   (raw R=120: {P2[mid,mid]:.6f})")
    # BE references (half-monodromy moduli)
    BEmod=np.exp(-np.pi*np.abs(b1o))        # |formal channel factor| expected ~ e^{-pi|b1|}
    print(f"     b1 (slope-ordered) = {np.round(b1o,6)};  e^(-pi|b1|) = {np.round(BEmod,6)}")
    print(f"     |C_ll|={abs(C[0,0]):.6f} |C_mm|={abs(C[1,1]):.6f} |C_hh|={abs(C[2,2]):.6f}")
    # [T2/T3] decompositions
    for tag,fac in (("LDU",gauss_LDU),("UDL",gauss_UDL)):
        if tag=="LDU": L,D,U=fac(C)
        else: U,D,L=fac(C)
        rec=(L@np.diag(D)@U) if tag=="LDU" else (U@np.diag(D)@L)
        err=np.max(np.abs(rec-C))
        print(f"\n[{tag}]  reconstruction err {err:.1e}")
        print(f"  |D| = {np.round(np.abs(D),6)}   vs  e^(-pi|b1|) = {np.round(BEmod,6)}")
        print(f"  ratio |D_j| / e^(-pi|b1_j|) = {np.round(np.abs(D)/BEmod,6)}")
        mods={ "L10":abs(L[1,0]),"L20":abs(L[2,0]),"L21":abs(L[2,1]),
               "U01":abs(U[0,1]),"U02":abs(U[0,2]),"U12":abs(U[1,2]) }
        print("  multiplier moduli:", {k:round(v,5) for k,v in mods.items()})
        # [T4] slots
        if tag=="LDU":
            sig_direct=abs(D[1])**2
            sig_full  =abs(D[1]+L[1,0]*D[0]*U[0,1])**2
            print(f"  sigma slot:  |D_m|^2={sig_direct:.6f}  |D_m + L10 D_l U01|^2={sig_full:.6f}"
                  f"  (C: {abs(C[1,1])**2:.6f})")
            print(f"  b slot:      |C_hl|^2={abs(C[2,0])**2:.6f}  = |L20 D_l|^2="
                  f"{abs(L[2,0]*D[0])**2:.6f}")
            print(f"  hi-corner relation: |D_h + L21 D_m U12 + L20 D_l U02| = {abs(D[2]+L[2,1]*D[1]*U[1,2]+L[2,0]*D[0]*U[0,2]):.6f}"
                  f"  vs |C_hh|={abs(C[2,2]):.6f}  vs e^(-pi|b1_h|)={BEmod[2]:.6f}")
        else:
            sig_full=abs(D[1]+U[1,2]*D[2]*L[2,1])**2
            print(f"  sigma slot (via hi): |D_m + U12 D_h L21|^2={sig_full:.6f}"
                  f"  (C: {abs(C[1,1])**2:.6f})")
            print(f"  b slot:      |C_hl|^2={abs(C[2,0])**2:.6f}  = |L20 D_h?|... raw |L20|^2|D_h|^2="
                  f"{abs(L[2,0]*D[2])**2:.6f}")
    # [T5] conjugation structure
    print("\n[T5] symmetry probes (max |.| of off-diagonal of each product):")
    for lbl,Mtest in (("C C^bar",C@C.conj()),("C C^bar^T",C@C.conj().T),
                      ("C^bar C",C.conj()@C),("C C^T",C@C.T)):
        off=Mtest-np.diag(Mtest.diagonal())
        print(f"   {lbl:10s} off-diag max {np.max(np.abs(off)):.3e}  diag {np.round(np.abs(Mtest.diagonal()),4)}")
    return C,b1o

if __name__=="__main__":
    report("canonical",(-2,0,3),(1,0.8,1.2),(-1,0.5,2.0),gold=0.214724)
    report("sampleB",(-1,0,1.5),(0.9,1.1,0.8),(-0.7,0.4,1.3),gold=0.021018)
    report("unsorted-a",(-1.6,0.3,2.1),(1.05,-0.75,0.95),(1.4,-0.9,0.6))
    print("\n[T6] formal monodromy: by formal_solution construction the only multivalued factor")
    print("     is u^{-i b1_j} per channel (DIAGONAL). No permutation part exists at an")
    print("     unramified point. The committed 'formal monodromy = oriented 3-cycle' claim")
    print("     conflated the incoherent-skeleton cycle (R16) with the formal monodromy.")
