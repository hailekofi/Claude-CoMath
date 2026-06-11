"""
ws_puncture_stokes.py -- (a) extraction of the puncture Stokes data of C (the exact
"Weber replacements"), from real-axis data + branch bookkeeping.

Each puncture of the lambda-system (R36) inherits the wild data of u=infinity: four jump
rays. The UPPER ray pair is the UDL factorization of S with the upper log branch
(ln(-R)=lnR+i pi; R25). The LOWER ray pair is the LDU factorization under the LOWER
branch normalization: S_low = S_up . diag(e^{-2 pi b_j})  (since F_-^{low} = F_-^{up}
diag(e^{-2 pi b})).  Direct semicircle transport is impossible (exponents ~ e^{aR^2/2});
the real-axis encoding is the correct extraction.

Battery (canonical + sampleB):
 [E1] upper factors: S_up = U.Du.L; |Du| vs elementary (e^{+pi|b_lo|}, e^{-pi|b_mid|},
      e^{-pi|b_hi|}); the 6 upper multipliers (modulus, phase) = upper puncture data.
 [E2] lower factors: S_low normalized; LDU: is |Dl| elementary with the MIRRORED dual
      weight? (resolves R25 open item (c): the LDU/lower-branch pairing). Lower 6
      multipliers extracted.
 [E3] conjugation probe (blind): are lower multipliers the complex conjugates of the
      upper ones (Schwarz reflection off the real axis)? cos/ratio table.
 [E4] torus-invariant chirality: Z = (S01 S12 S20)/(S10 S21 S02): modulus + phase
      (the directed-cycle asymmetry of C; modulus is fixed by P alone).
 [E5] rank-1 probe: conj(S) ?= D1 S^{-1} D2 (time-reversal-type identity): singular
      values of the ratio matrix; PASS would tie lower data to upper exactly.
 [E6] Weber-shadow table: extracted |multiplier| vs the 2-level model sqrt(e^{2 pi
      be_ij}-1) for the pair (i,j) the slot couples: the exact-vs-shadow gap per slot.
Reproduce: python3 ws_puncture_stokes.py
"""
import numpy as np
from scipy.integrate import solve_ivp

def build(eps,gam,a):
    eps=np.asarray(eps,float);gam=np.asarray(gam,float);a=np.asarray(a,float)
    H0=np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            if i!=j: H0[i,j]=gam[i]*gam[j]*(a[i]-a[j])/(eps[i]-eps[j])
        H0[i,i]=-sum(gam[k]**2*(a[i]-a[k])/(eps[i]-eps[k]) for k in range(3) if k!=i)
    return H0,a

def connection_S(eps,gam,a,R,rtol=1e-12,atol=1e-13):
    H0,a=build(eps,gam,a)
    b1=np.array([sum(H0[j,k]**2/(a[j]-a[k]) for k in range(3) if k!=j) for j in range(3)])
    T1=np.zeros((3,3),complex)
    for k in range(3):
        for j in range(3):
            if k!=j: T1[k,j]=H0[k,j]/(a[j]-a[k])
    def frame(u):
        ph=np.array([-1j*(a[j]*u**2/2+H0[j,j]*u+b1[j]*np.log(complex(u))) for j in range(3)])
        return (np.eye(3)+T1/u)@np.diag(np.exp(ph))
    Y0=frame(-R)
    sol=solve_ivp(lambda u,Y:(-1j*(H0+u*np.diag(a))@Y.reshape(3,3)).reshape(-1),
                  [-R,R],Y0.reshape(-1).astype(complex),rtol=rtol,atol=atol,method='DOP853')
    YR=sol.y[:,-1].reshape(3,3)
    C=np.linalg.solve(frame(R),YR)
    lo,mid,hi=np.argsort(a); perm=[lo,mid,hi]
    C=C[np.ix_(perm,perm)]; b=b1[perm]
    return C@np.diag(np.exp(-np.pi*b)),b      # upper-branch physical S, slope order

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
    return J@L@J,(J@np.diag(D)@J).diagonal().copy(),J@U@J

def mults(U,L):
    return dict(U01=U[0,1],U02=U[0,2],U12=U[1,2],L10=L[1,0],L20=L[2,0],L21=L[2,1])

def run(nm,eps,gam,a):
    print("="*88); print(f"CASE {nm}"); print("="*88)
    S,b=connection_S(eps,gam,a,120.)
    BEel=np.array([np.exp(+np.pi*abs(b[0])),np.exp(-np.pi*abs(b[1])),np.exp(-np.pi*abs(b[2]))])
    # E1 upper
    Uu,Du,Lu=gauss_UDL(S)
    print(f"[E1] UPPER (UDL, branch +i pi): |Du|/elementary = {np.round(np.abs(Du)/BEel,5)}")
    mu=mults(Uu,Lu)
    print("     upper multipliers:", {k:f"{abs(v):.4f} @ {np.angle(v):+.3f}" for k,v in mu.items()})
    # E2 lower: S_low = S diag(e^{-2 pi b}) renormalized to unit-modulus rows? --
    # physical lower normalization: strip e^{-pi b} with the LOWER branch instead:
    # F_-^{low} = F_-^{up} diag(e^{-2pi b})  =>  S_low = S_up diag(e^{+2pi b})*e^{-2pi b}...
    # operationally: S_low = S @ diag(e^{-2*pi*b}) then re-normalize by e^{+pi b} overall:
    Slow=S@np.diag(np.exp(-np.pi*b))@np.diag(np.exp(-np.pi*b))   # = S e^{-2pi b}
    Slow=Slow@np.diag(np.exp(2*np.pi*b))**0  # keep as is; test both raw and renorm
    # the clean object: physical lower S with its own modulus normalization:
    Slow=S@np.diag(np.exp(-2*np.pi*b))
    Slow=Slow/np.exp(np.log(np.abs(np.linalg.det(Slow)))/3)      # unit |det|
    Ll,Dl,Ul=gauss_LDU(Slow)
    BEel_mirror=np.array([np.exp(-np.pi*abs(b[0])),np.exp(-np.pi*abs(b[1])),np.exp(+np.pi*abs(b[2]))])
    print(f"[E2] LOWER (LDU on branch -i pi normalization): |Dl|/mirror-elementary = "
          f"{np.round(np.abs(Dl)/BEel_mirror,5)}")
    ml=mults(Ul,Ll)
    print("     lower multipliers:", {k:f"{abs(v):.4f} @ {np.angle(v):+.3f}" for k,v in ml.items()})
    # E3 conjugation probe
    print("[E3] conjugation probe |lower| vs |conj-upper-partner| (U<->L transpose pairing):")
    pairs=[("U01","L10"),("U02","L20"),("U12","L21"),("L10","U01"),("L20","U02"),("L21","U12")]
    for kl,ku in pairs[:3]:
        print(f"     |{kl}_low|={abs(ml[kl]):.4f}  vs |{ku}_up|={abs(mu[ku]):.4f}  "
              f"ratio={abs(ml[kl])/max(abs(mu[ku]),1e-12):.4f}")
    # E4 chirality
    Z=(S[0,1]*S[1,2]*S[2,0])/(S[1,0]*S[2,1]*S[0,2])
    print(f"[E4] cyclic chirality Z = (S01 S12 S20)/(S10 S21 S02): |Z|={abs(Z):.4f}  "
          f"arg Z={np.angle(Z):+.4f}   (torus-invariant)")
    # E5 rank-1 probe conj(S) ?= D1 S^{-1} D2
    M=np.conj(S)/np.linalg.inv(S)
    sv=np.linalg.svd(np.log(np.abs(M)),compute_uv=False)
    # rank-1 in log-modulus <=> log|M_ij| = xi_i + eta_j: test via double-centering residual
    LM=np.log(np.abs(M))
    res=LM-LM.mean(1,keepdims=True)-LM.mean(0,keepdims=True)+LM.mean()
    print(f"[E5] conj(S) = D1 S^-1 D2 probe: double-centered residual of log|conjS/S^-1| "
          f"max={np.max(np.abs(res)):.3e}  ({'PASS' if np.max(np.abs(res))<1e-3 else 'FAIL'})")
    # E6 Weber shadow
    eps_=np.asarray(eps,float);gam_=np.asarray(gam,float);a_=np.asarray(a,float)
    lo,mid,hi=np.argsort(a_)
    def be(i,j): return gam_[i]**2*gam_[j]**2*abs(a_[i]-a_[j])/(eps_[i]-eps_[j])**2
    pair={"U01":(lo,mid),"L10":(lo,mid),"U12":(mid,hi),"L21":(mid,hi),"U02":(lo,hi),"L20":(lo,hi)}
    print("[E6] Weber-shadow table (slot couples pair (i,j); model sqrt(e^{2pi be_ij}-1)):")
    for k in ("U01","L10","U12","L21","U02","L20"):
        i,j=pair[k]; model=np.sqrt(np.exp(2*np.pi*be(i,j))-1)
        print(f"     {k}: extracted {abs(mu[k]):.4f}   2-level model {model:.4f}   "
              f"ratio {abs(mu[k])/model:.3f}")
    return S,b,mu,ml

if __name__=="__main__":
    run("canonical",(-2,0,3),(1,.8,1.2),(-1,.5,2.))
    run("sampleB",(-1,0,1.5),(.9,1.1,.8),(-.7,.4,1.3))
