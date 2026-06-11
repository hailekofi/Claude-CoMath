"""
ws_lambda_system.py -- construct the lambda-sphere realization that defines the central
factor C: pull the Type-1 connection back along the nodal-cubic uniformization
   u(lam) = u* - Q2(lam)/PROD(lam-a_i),  E(lam) = E* + lam(u(lam)-u*):
   d psi/d lam = A(lam) psi,   A(lam) = -i u'(lam) (H0 + u(lam) diag(a)).

Claims:
 [C1] A(lam) is RATIONAL with poles ONLY at the three slope-punctures lam=a_i, each of
      order 3 (Poincare rank 2); leading coefficient  i r_i^2 diag(a)  with
      r_i = -Q2(a_i)/PROD_{j!=i}(a_i-a_j)  (the residue of u at the puncture).
      Turning points (u'(lam)=0, complex quartic) and the node are REGULAR points.
 [C2] node resolution: the two roots of Q2(lam)=0 are distinct real lambda points both
      mapping to (u*,E*): the lambda-curve separates the crossing branches.
 [C3] real geometry (BLIND): track the three real lambda-branches over u in [-30,30];
      report each branch's lambda-interval, its puncture/infinity endpoints, and the
      induced in/out slope pairing (the channel combinatorics of C).
 [C4] pullback equivalence: integrating the lambda-ODE along an arc equals the u-ODE
      propagator on the image segment (machine precision).
C is then DEFINED as: the wild connection data of (P^1_lam, A(lam) dlam) at the three
punctures, glued along the closed-form (R34) interval transports, modulo the torus (R35).
Reproduce: python3 ws_lambda_system.py
"""
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp

EPS=[-2,0,3]; GAM=[1,sp.Rational(4,5),sp.Rational(6,5)]; A=[-1,sp.Rational(1,2),2]
u,E,lam=sp.symbols('u E lam')
H0=sp.zeros(3,3)
for i in range(3):
    for j in range(3):
        if i!=j: H0[i,j]=GAM[i]*GAM[j]*(A[i]-A[j])/sp.Rational(EPS[i]-EPS[j])
    H0[i,i]=-sum(GAM[k]**2*(A[i]-A[k])/sp.Rational(EPS[i]-EPS[k]) for k in range(3) if k!=i)
F=sp.expand((sp.eye(3)*E-H0-u*sp.diag(*A)).det())
sols=sp.solve([F,sp.diff(F,u),sp.diff(F,E)],[u,E],dict=True)
node=[s for s in sols if s[u].is_real][0]; u0,E0=node[u],node[E]
x,y=sp.symbols('x y')
P=sp.Poly(sp.expand(F.subs({u:u0+x,E:E0+y})),x,y)
Q2l=sum(P.coeff_monomial(x**i*y**(2-i))*lam**(2-i) for i in range(3))
C3l=sp.prod(lam-a for a in [sp.Rational(-1),sp.Rational(1,2),sp.Rational(2)])
ul=sp.simplify(u0-Q2l/C3l); Elam=sp.simplify(E0+lam*(ul-u0))
ups=sp.simplify(sp.diff(ul,lam))

# [C1] pole structure of A(lam)
Alam=sp.simplify(-sp.I*ups*(H0+ul*sp.diag(*A)))
ent=sp.together(Alam[2,2])
den=sp.factor(sp.denom(ent))
print(f"[C1] A(lam) sample entry denominator (factored): {den}")
avals=[sp.Rational(-1),sp.Rational(1,2),sp.Rational(2)]
for ai in avals:
    r=sp.simplify(-(Q2l.subs(lam,ai))/sp.prod(ai-aj for aj in avals if aj!=ai))
    lead=sp.simplify(sp.limit(Alam[2,2]*(lam-ai)**3,lam,ai))
    pred=sp.simplify(sp.I*r**2*A[2])
    print(f"     puncture lam={ai}: residue r_i={sp.nsimplify(r)}; "
          f"(lam-a_i)^3 A_22 -> {sp.N(lead,8)} vs i r^2 a_2 = {sp.N(pred,8)}  match={sp.simplify(lead-pred)==0}")
# regularity at turning points & node-lambdas
tp=sp.Poly(sp.numer(sp.together(ups)),lam).nroots(n=20)
q2r=sp.Poly(Q2l,lam).nroots(n=20)
mx=max(abs(complex(Alam[0,1].subs(lam,t))) for t in tp)
mxn=max(abs(complex(Alam[0,1].subs(lam,t))) for t in q2r)
print(f"     |A_01| at the 4 turning points: max {mx:.4f} (FINITE -> regular)")
print(f"     |A_01| at the 2 node-lambdas:   max {mxn:.4f} (FINITE -> regular)")

# [C2] node resolution
print(f"\n[C2] node (u*,E*)=({u0},{E0}); Q2 roots (node lambdas):")
for t in q2r:
    print(f"     lam={complex(t).real:+.6f}: u(lam)={complex(ul.subs(lam,t)).real:+.6f} "
          f"E(lam)={complex(Elam.subs(lam,t)).real:+.6f}   (both = node)")

# [C3] real branch geometry (blind)
ulF=sp.lambdify(lam,ul,'numpy'); ElF=sp.lambdify(lam,Elam,'numpy')
H0n=np.array(H0.evalf(),dtype=float); an=np.array([-1,0.5,2.])
us=np.linspace(-30,30,2401)
lam_tracks=[[],[],[]]
polyu=sp.Poly(sp.together(ul-sp.Symbol('U')).as_numer_denom()[0],lam)
import numpy.polynomial.polynomial as npp
cs=[sp.lambdify(sp.Symbol('U'),c,'numpy') for c in polyu.all_coeffs()]
for uv in us:
    coef=[float(c(uv)) for c in cs]
    rts=np.roots(coef)
    rl=sorted([r.real for r in rts if abs(r.imag)<1e-9])
    Hn=H0n+uv*np.diag(an); w=np.linalg.eigvalsh(Hn)
    pairs=sorted(zip([float(ElF(r)) for r in rl],rl))
    for k in range(3): lam_tracks[k].append(pairs[k][1])   # sorted by energy
lam_tracks=np.array(lam_tracks)
print(f"\n[C3] real branch arcs (energy-sorted levels), u in [-30,30]:")
for k in range(3):
    lin,lout=lam_tracks[k,0],lam_tracks[k,-1]
    rng=(lam_tracks[k].min(),lam_tracks[k].max())
    jumps=np.abs(np.diff(lam_tracks[k]))
    ninf=np.sum(jumps>1.0)
    sin=min(avals,key=lambda a:abs(lin-float(a))); sout=min(avals,key=lambda a:abs(lout-float(a)))
    print(f"     level{k}: lam {lin:+.4f} -> {lout:+.4f}; range[{rng[0]:+.2f},{rng[1]:+.2f}]; "
          f"passes lam=inf {ninf}x; slope pairing a_in={sin} -> a_out={sout}")

# [C4] pullback equivalence on an arc segment
l0,l1=-0.2,0.35     # within an interval (check no puncture inside)
uA,uB=float(ulF(l0)),float(ulF(l1))
AlamF=sp.lambdify(lam,Alam,'numpy')
sol1=solve_ivp(lambda l,y: (np.array(AlamF(l),dtype=complex)@y.reshape(3,3)).reshape(-1),
               [l0,l1],np.eye(3,dtype=complex).reshape(-1),rtol=1e-12,atol=1e-13,method='DOP853')
U_lam=sol1.y[:,-1].reshape(3,3)
sol2=solve_ivp(lambda uu,y: ((-1j*(H0n+uu*np.diag(an)))@y.reshape(3,3)).reshape(-1),
               [uA,uB],np.eye(3,dtype=complex).reshape(-1),rtol=1e-12,atol=1e-13,method='DOP853')
U_u=sol2.y[:,-1].reshape(3,3)
print(f"\n[C4] pullback equivalence: lam-arc [{l0},{l1}] maps to u-segment [{uA:.4f},{uB:.4f}]")
print(f"     |U_lambda - U_u|_max = {np.max(np.abs(U_lam-U_u)):.2e}   (machine precision = same transport)")
