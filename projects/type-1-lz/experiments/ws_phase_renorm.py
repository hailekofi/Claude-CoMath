"""
ws_phase_renorm.py -- the phase-renormalized definition of S (user's point: cutoff S(R)
converges only up to accelerating phases; the definition must reflect it).

S := lim_R F+(R)^{-1} U(R,-R) F-(-R), F_pm = Thome frames with THREE divergent phase tiers
   per channel: quadratic a_j u^2/2, linear H0_jj u, LOGARITHMIC b_j ln u (b_j = signed BE).
Claims tested (canonical anchor):
 [T1] tier-1 stripping only (quad+lin): |S_ij| converge BUT entry phases DRIFT ~ c ln R;
      measured drift slope per entry vs predicted (b_row + b_col-ish) combination.
 [T2] tier-2 stripping (+log): entry phases converge (rate ~1/R).
 [T3] the residual ambiguity is the constant torus S -> D+ S D-: the gauge-invariant
      phases arg(S_ij S_kl S*_il S*_kj) converge IDENTICALLY under both tier-1 and tier-2
      stripping (the drift cancels in cross-ratios) -- the physical phase content.
 [T4] sigma-slot consistency: the relative phase in sigma = |Delta_m + U12 Delta_h L21|^2
      is torus-invariant (numerically: sigma from the factorization is identical under
      arbitrary random torus twists of S).
Reproduce: python3 ws_phase_renorm.py
"""
import numpy as np
from scipy.integrate import solve_ivp

EPS=np.array([-2.,0.,3.]); GAM=np.array([1.,.8,1.2]); A=np.array([-1.,.5,2.])
H0=np.zeros((3,3))
for i in range(3):
    for j in range(3):
        if i!=j: H0[i,j]=GAM[i]*GAM[j]*(A[i]-A[j])/(EPS[i]-EPS[j])
    H0[i,i]=-sum(GAM[k]**2*(A[i]-A[k])/(EPS[i]-EPS[k]) for k in range(3) if k!=i)
b1=np.array([sum(H0[j,k]**2/(A[j]-A[k]) for k in range(3) if k!=j) for j in range(3)])
T1=np.zeros((3,3),complex)
for k in range(3):
    for j in range(3):
        if k!=j: T1[k,j]=H0[k,j]/(A[j]-A[k])

def frame(u,withlog):
    ph=-1j*(A*u**2/2+np.diag(H0)*u+(b1*np.log(complex(u)) if withlog else 0))
    return (np.eye(3)+T1/u)@np.diag(np.exp(ph))

def S_of(R,withlog,rtol=1e-11,atol=1e-12):
    Y0=frame(-R,withlog)
    sol=solve_ivp(lambda u,Y:(-1j*(H0+u*np.diag(A))@Y.reshape(3,3)).reshape(-1),
                  [-R,R],Y0.reshape(-1).astype(complex),rtol=rtol,atol=atol,method='DOP853')
    return np.linalg.solve(frame(R,withlog),sol.y[:,-1].reshape(3,3))

Rs=[40.,80.,160.]
S1={R:S_of(R,False) for R in Rs}   # tier-1: no log stripped
S2={R:S_of(R,True)  for R in Rs}   # tier-2: log stripped

print("="*86)
print("[T1] tier-1 (quad+lin only): moduli converge, phases DRIFT logarithmically")
print("="*86)
print(f"  | |S(160)|-|S(80)| |_max  tier-1: {np.max(np.abs(np.abs(S1[160.])-np.abs(S1[80.]))):.2e}"
      f"   tier-2: {np.max(np.abs(np.abs(S2[160.])-np.abs(S2[80.]))):.2e}   (both converge)")
print("\n  phase drift per e-fold of R, entry (i,j): measured [arg S(160)-arg S(80)]/ln2 ; predicted -(b_i+b_j)")
drift=(np.angle(S1[160.])-np.angle(S1[80.]))/np.log(2)
pred=-(b1[:,None]+b1[None,:])
for i in range(3):
    print("   ",["%+.4f/%+.4f"%(drift[i,j],pred[i,j]) for j in range(3)])

print("\n"+"="*86)
print("[T2] tier-2 (+log): entry phases converge")
print("="*86)
d12=np.max(np.abs(np.angle(S2[160.]*np.conj(S2[80.]))))
d01=np.max(np.abs(np.angle(S2[80.]*np.conj(S2[40.]))))
print(f"  max|arg S(160)-arg S(80)| = {d12:.2e}   (vs {d01:.2e} for 80-40: ratio {d01/max(d12,1e-300):.1f} ~ rate 1/R^p)")

print("\n"+"="*86)
print("[T3] torus-invariant phase cross-ratios: converge under BOTH strippings, equal values")
print("="*86)
def cross(S,i,j,k,l): return np.angle(S[i,j]*S[k,l]*np.conj(S[i,l])*np.conj(S[k,j]))
combos=[(0,0,1,1),(0,1,1,2),(1,0,2,1)]
for c in combos:
    v1=[cross(S1[R],*c) for R in Rs]; v2=[cross(S2[R],*c) for R in Rs]
    print(f"  arg-cross{c}: tier-1 {['%.6f'%x for x in v1]} | tier-2 {['%.6f'%x for x in v2]}"
          f"  agree: {abs(v1[-1]-v2[-1])<1e-4}")

print("\n"+"="*86)
print("[T4] sigma is torus-invariant through the factorization")
print("="*86)
S=S2[160.]@np.diag(np.exp(-np.pi*b1))     # physical normalization (moduli)
rng=np.random.default_rng(5)
def gauss_UDL(C):
    J=np.eye(3)[::-1]
    def LDU(C):
        L=np.eye(3,dtype=complex);U=np.eye(3,dtype=complex);D=np.zeros(3,complex)
        D[0]=C[0,0];L[1,0]=C[1,0]/D[0];L[2,0]=C[2,0]/D[0]
        U[0,1]=C[0,1]/D[0];U[0,2]=C[0,2]/D[0]
        D[1]=C[1,1]-L[1,0]*D[0]*U[0,1]
        L[2,1]=(C[2,1]-L[2,0]*D[0]*U[0,1])/D[1]
        U[1,2]=(C[1,2]-L[1,0]*D[0]*U[0,2])/D[1]
        D[2]=C[2,2]-L[2,0]*D[0]*U[0,2]-L[2,1]*D[1]*U[1,2]
        return L,D,U
    Lr,Dr,Ur=LDU(J@C@J)
    return J@Lr@J,(J@np.diag(Dr)@J).diagonal().copy(),J@Ur@J
vals=[]
for t in range(4):
    Dp=np.diag(np.exp(1j*rng.uniform(0,2*np.pi,3))); Dm=np.diag(np.exp(1j*rng.uniform(0,2*np.pi,3)))
    St=Dp@S@Dm
    U,D,L=gauss_UDL(St)
    vals.append(abs(D[1]+U[1,2]*D[2]*L[2,1])**2)
print(f"  sigma from slot equation under 4 random torus twists: {['%.7f'%v for v in vals]}")
print(f"  spread {max(vals)-min(vals):.1e}   (torus-invariant; gold 0.214724)")
