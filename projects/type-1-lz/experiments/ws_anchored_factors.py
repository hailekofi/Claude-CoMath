"""
ws_anchored_factors.py -- constructive verification of the two REDUCIBLE ingredients of the
anchored factorization  S = W+ D+ . C . D- W-^{-1}  (canonical Type-1 anchor):

 (iii) RATIONAL FRAMES: eigenvectors = adjugate columns of (H0 + u(lam)A - E(lam) I),
       rational in lam, coefficients rational in {gam,eps,a}.  Verified vs eigh.
 (ii)  PHASE DIAGONALS: int E du  pulls back to  int E(lam) u'(lam) dlam, a rational
       differential with poles ONLY at lam = a_i  =>  antiderivative = rational + sum_i
       c_i log(lam - a_i)  (algebro-logarithmic, logs AT THE SLOPES).  Symbolic closed
       form vs numerical quadrature of eigh differences.
 (i)   the central factor C is NOT reducible (a-priori: else sigma Liouvillian, contra
       R26/R27; in-house: WS-E Stuckelberg-band exclusion, M7 negative, R16 2.5-100x).
       Witness number printed.

Reproduce: python3 ws_anchored_factors.py   (sympy/numpy; exact rationals)
"""
import numpy as np, sympy as sp
from scipy.integrate import quad

# ----- canonical anchor, exact rationals -----
EPS=[-2,0,3]; GAM=[1,sp.Rational(4,5),sp.Rational(6,5)]; A=[-1,sp.Rational(1,2),2]
u,E,lam=sp.symbols('u E lam')

H0=sp.zeros(3,3)
for i in range(3):
    for j in range(3):
        if i!=j: H0[i,j]=GAM[i]*GAM[j]*(A[i]-A[j])/sp.Rational(EPS[i]-EPS[j])
    H0[i,i]=-sum(GAM[k]**2*(A[i]-A[k])/sp.Rational(EPS[i]-EPS[k]) for k in range(3) if k!=i)
F=sp.expand((sp.eye(3)*E-H0-u*sp.diag(*A)).det())

# node (exact): solve F=F_u=F_E=0
sols=sp.solve([F,sp.diff(F,u),sp.diff(F,E)],[u,E],dict=True)
node=[s for s in sols if s[u].is_real][0]
u0,E0=node[u],node[E]
print(f"node (u*,E*) = ({u0}, {E0}) = ({float(u0):.6f},{float(E0):.6f})  [exact rationals]")

# uniformization
x,y=sp.symbols('x y')
Fs=sp.expand(F.subs({u:u0+x,E:E0+y}))
P=sp.Poly(Fs,x,y)
Q2l=sum(P.coeff_monomial(x**i*y**(2-i))*lam**(2-i) for i in range(3))
C3l=sp.prod(lam-sp.Rational(a) if isinstance(a,int) else lam-a for a in A)
ul=sp.simplify(u0-Q2l/C3l); El=sp.simplify(E0+lam*(ul-u0))
print(f"u(lam) = {u0} - ({sp.nsimplify(Q2l)})/((lam+1)(lam-1/2)(lam-2))")

# ---------- (iii) rational frames ----------
M=H0+ul*sp.diag(*A)-El*sp.eye(3)
adj=M.adjugate()
vcol=sp.simplify(adj[:,0])           # rational eigenvector (unnormalized) as function of lam
degs=[(sp.degree(sp.numer(sp.together(vcol[k])),lam),sp.degree(sp.denom(sp.together(vcol[k])),lam)) for k in range(3)]
print(f"\n(iii) rational eigenvector v(lam) = adjugate column; entry degrees (num,den): {degs}")
uu=sp.Rational(7,10)
lams=sorted([r for r in sp.solve(sp.Eq(ul,uu),lam) if r.is_real],key=lambda r:float(El.subs(lam,r)))
Hn=np.array((H0+uu*sp.diag(*A)).evalf(),dtype=float)
w,V=np.linalg.eigh(Hn)
errs=[]
for k,l in enumerate(lams):
    vk=np.array([complex(vcol[m].subs(lam,l)) for m in range(3)],dtype=complex).real
    vk/=np.linalg.norm(vk)
    errs.append(1-abs(vk@V[:,k]))
print(f"      vs eigh at u=0.7: 1-|cos| per sheet = {[f'{e:.1e}' for e in errs]}")

# ---------- (ii) phase integral in closed form ----------
g=sp.together(sp.expand(El*sp.diff(ul,lam)))     # E(lam) u'(lam): rational, poles at a_i
Fanti=sp.integrate(sp.apart(g,lam),lam)          # rational + logs at the slopes
print(f"\n(ii) antiderivative of E dlam-pullback (structure): "
      f"{sp.count_ops(Fanti)} ops; log arguments: "
      f"{sorted([str(t.args[0]) for t in Fanti.atoms(sp.log)])}")
uA,uB=sp.Rational(7,10),sp.Rational(3,2)
def sheet_lams(uv):
    return sorted([r for r in sp.solve(sp.Eq(ul,uv),lam) if r.is_real],
                  key=lambda r: float(El.subs(lam,r)))
lA,lB=sheet_lams(uA),sheet_lams(uB)
i,j=2,1                                          # top two sheets
Phi_closed=float(sp.re((Fanti.subs(lam,lB[i])-Fanti.subs(lam,lA[i])
            -(Fanti.subs(lam,lB[j])-Fanti.subs(lam,lA[j]))).evalf(30)))
def Ediff(uv):
    Hn=np.array((H0+uv*sp.diag(*A)).evalf(),dtype=float)
    w=np.linalg.eigvalsh(Hn); return w[i]-w[j]
Phi_num,err=quad(Ediff,float(uA),float(uB),epsabs=1e-13,epsrel=1e-13)
print(f"      Phi_{i}{j} = int_({float(uA)})^({float(uB)}) (E_{i}-E_{j}) du:")
print(f"        closed form (rational + logs at slopes) = {Phi_closed:.12f}")
print(f"        numerical quadrature (eigh)             = {Phi_num:.12f}   |diff|={abs(Phi_closed-Phi_num):.1e}")

# ---------- (i) the witness that C is irreducible ----------
be=lambda i,j: float(GAM[i]**2*GAM[j]**2*abs(A[i]-A[j])/sp.Rational(EPS[i]-EPS[j])**2)
P_inc=np.exp(-2*np.pi*(be(0,1)+be(1,2)))
print(f"\n(i)  central factor C: Weber/incoherent model P_mm = {P_inc:.4f} vs true 0.2147 "
      f"(x{0.2147/P_inc:.1f} off) -> C irreducible (R16/WS-E/M7; a-priori via R26/R27).")
