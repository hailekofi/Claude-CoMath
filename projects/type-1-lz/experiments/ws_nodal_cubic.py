"""
ws_nodal_cubic.py -- the Type-1 matrix-structure anchor: the spectral curve
F(u,E) = det(E I - H0 - u A) = 0 is a PLANE CUBIC, and

   generic N=3 sweep  -> SMOOTH cubic (genus 1, elliptic): frames/periods theta/elliptic
   Type-1             -> NODAL  cubic (genus 0, rational): node = the protected crossing;
                          explicit rational uniformization by lines through the node.

Tests:
 [1] singularity discriminator: solve F = F_u = F_E = 0.  Type-1: a real node (matches
     R3/R8); generic symmetric H0 (same slopes): NO singular point (smooth cubic).
 [2] discriminant root structure: Disc_E F(u) degree 6 in u: Type-1 = 4 complex simple
     (two conj pairs, the Q4 windows) + 1 REAL DOUBLE root (the node);
     generic = 6 simple complex roots (3 conj pairs).  Genus: (branch pts)/2 - 2:
     Type-1 4/2-2 = 0, generic 6/2-2 = 1.
 [3] EXPLICIT rational uniformization (Type-1): shift to the node (u0,E0); then
     F = Q2(x,y) + C3(x,y) (double point: no constant/linear part); on lines y = lam x:
       x(lam) = -Q2(1,lam)/C3(1,lam),  u = u0 + x,  E = E0 + lam x
     -- verify F(u(lam),E(lam)) == 0 identically (sympy) and that the three sheets/
     branches over a given u correspond to the 3 roots lam of u(lam)=u.
All data entering the uniformization -- (u0,E0), Q2, C3 -- are algebraic in {gam,eps,a}:
the anchored seed of the factorization S = W+(params) . T_rational . W-(params)^{-1}.

Reproduce: python3 ws_nodal_cubic.py    (sympy/numpy only; no ODE solves)
"""
import numpy as np, sympy as sp

rng=np.random.default_rng(11)

def H0_type1(eps,gam,a):
    H0=sp.zeros(3,3)
    for i in range(3):
        for j in range(3):
            if i!=j: H0[i,j]=sp.nsimplify(gam[i]*gam[j]*(a[i]-a[j]))/sp.nsimplify(eps[i]-eps[j])
        H0[i,i]=-sum(sp.nsimplify(gam[k]**2*(a[i]-a[k]))/sp.nsimplify(eps[i]-eps[k]) for k in range(3) if k!=i)
    return H0

def curve(H0,a):
    u,E=sp.symbols('u E')
    M=sp.eye(3)*E - H0 - u*sp.diag(*[sp.nsimplify(x) for x in a])
    return sp.expand(M.det()),u,E

def sing_points(F,u,E):
    sols=sp.solve([F,sp.diff(F,u),sp.diff(F,E)],[u,E],dict=True)
    return sols

def disc_structure(F,u,E):
    D=sp.Poly(sp.discriminant(F,E),u)
    rts=sp.Poly(D).nroots(n=30,maxsteps=200)
    # cluster roots to find multiplicities
    rts=[complex(r) for r in rts]
    used=[False]*len(rts); clusters=[]
    for i,r in enumerate(rts):
        if used[i]:continue
        c=[r];used[i]=True
        for j in range(i+1,len(rts)):
            if not used[j] and abs(rts[j]-r)<1e-8: c.append(rts[j]);used[j]=True
        clusters.append((np.mean(c),len(c)))
    return clusters

print("="*84)
print("[1]+[2] SMOOTH (generic) vs NODAL (Type-1) spectral cubic")
print("="*84)
cases=[]
cases.append(("TYPE-1 canonical",H0_type1([-2,0,3],[1,sp.Rational(4,5),sp.Rational(6,5)],[-1,sp.Rational(1,2),2]),[-1,sp.Rational(1,2),2]))
cases.append(("TYPE-1 second",H0_type1([-1,0,sp.Rational(3,2)],[sp.Rational(9,10),sp.Rational(11,10),sp.Rational(4,5)],[sp.Rational(-7,10),sp.Rational(2,5),sp.Rational(13,10)]),[sp.Rational(-7,10),sp.Rational(2,5),sp.Rational(13,10)]))
for tag in ("generic-1","generic-2"):
    M=rng.standard_normal((3,3)); Hf=np.round(0.6*(M+M.T)/2,3)
    Hg=sp.Matrix(3,3,lambda i,j: sp.Rational(int(round(Hf[i,j]*1000)),1000))
    cases.append((f"GENERIC {tag}",Hg,[-1,sp.Rational(1,2),2]))

node_data={}
for name,H0,a in cases:
    F,u,E=curve(H0,a)
    sols=sing_points(F,u,E)
    real_sings=[s for s in sols if all(sp.im(sp.nsimplify(v)).is_zero or abs(complex(v).imag)<1e-10 for v in s.values())]
    clus=disc_structure(F,u,E)
    mults=sorted([m for _,m in clus],reverse=True)
    nreal_double=sum(1 for r,m in clus if m>=2 and abs(r.imag)<1e-7)
    print(f"\n  {name}:")
    print(f"    singular points of the cubic: {len(sols)} total, {len(real_sings)} real "
          f"{'-> NODAL (genus 0)' if sols else '-> SMOOTH (genus 1)'}")
    if real_sings:
        s=real_sings[0]; u0=complex(s[u]).real; E0=complex(s[E]).real
        print(f"    node (u*,E*) = ({u0:.6f}, {E0:.6f})")
        node_data[name]=(F,u,E,s[u],s[E])
    print(f"    Disc_E root multiplicities in u: {mults}  (real double roots: {nreal_double})")
    nb=sum(m for r,m in clus if m==1)
    print(f"    simple branch points: {nb}  =>  genus = {nb//2-2}")

print()
print("="*84)
print("[3] EXPLICIT RATIONAL UNIFORMIZATION of the Type-1 nodal cubic (canonical)")
print("="*84)
F,u,E,u0,E0=node_data["TYPE-1 canonical"]
x,y,lam=sp.symbols('x y lam')
Fs=sp.expand(F.subs({u:u0+x,E:E0+y}))
poly=sp.Poly(Fs,x,y)
const=poly.coeff_monomial(1); linx=poly.coeff_monomial(x); liny=poly.coeff_monomial(y)
print(f"  at the node: |const|={abs(complex(const)):.2e} |F_x|={abs(complex(linx)):.2e} "
      f"|F_y|={abs(complex(liny)):.2e}   (all 0 = double point)")
Q2=sum(poly.coeff_monomial(x**i*y**(2-i))*x**i*y**(2-i) for i in range(3))
C3=sum(poly.coeff_monomial(x**i*y**(3-i))*x**i*y**(3-i) for i in range(4))
Q2l=sp.expand(Q2.subs({x:1,y:lam})); C3l=sp.expand(C3.subs({x:1,y:lam}))
xl=sp.simplify(-Q2l/C3l)
ul=u0+xl; El=E0+lam*xl
resid=sp.simplify(sp.expand(F.subs({u:ul,E:El})*sp.denom(sp.together(xl))**3))
print(f"  uniformization: u(lam) = {sp.nsimplify(u0,rational=False)} - Q2(1,lam)/C3(1,lam),  "
      f"E(lam) = E* + lam*(u(lam)-u*)")
print(f"  Q2(1,lam) = {sp.N(Q2l,6)}")
print(f"  C3(1,lam) = {sp.N(C3l,6)}")
print(f"  F(u(lam),E(lam)) == 0 identically: {sp.simplify(resid)==0 or sp.N(sp.expand(resid),30)==0}")
# sheet check at a sample u
uu=0.7
lams=sp.solve(sp.Eq(ul,uu),lam)
Es=sorted([complex(El.subs(lam,l)).real for l in lams])
H0c=H0_type1([-2,0,3],[1,sp.Rational(4,5),sp.Rational(6,5)],[-1,sp.Rational(1,2),2])
Hnum=np.array(H0c.evalf(),dtype=float)+uu*np.diag([-1,0.5,2])
Etrue=sorted(np.linalg.eigvalsh(Hnum))
print(f"  sheets over u={uu}: E(lam roots) = {np.round(Es,8)}")
print(f"               eigh = {np.round(Etrue,8)}   max|diff|={max(abs(a-b) for a,b in zip(Es,Etrue)):.1e}")
print()
print("  => the eigenvalue sheets, hence eigenprojectors and all curve periods, are RATIONAL")
print("     in lam with coefficients algebraic in {gam,eps,a}: the anchored seed for the")
print("     factorization S = W+(params) . T_rational . W-(params)^{-1}.")
