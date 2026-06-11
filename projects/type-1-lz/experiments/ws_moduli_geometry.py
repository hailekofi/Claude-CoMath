"""
ws_moduli_geometry.py -- topology/geometry of the 4-dim Type-1 moduli space with
coordinates (be_lm, be_mh, be_lh, chi)  [R31/R32].

Pure geometry (no ODE): sample valid Type-1 params, compute the 4 invariants, and probe
the IMAGE region and its structure.
 [G1] chi range: is chi in (0,1)?  (analytic: cross-ratio of two conjugate pairs
      {z,zbar,w,wbar} = |z-w|^2/|z-wbar|^2 in (0,1); confirm numerically)
 [G2] be-octant: do (be_lm,be_mh,be_lh) fill R^3_{>0}? marginals of log10(be); coverage.
 [G3] PRODUCT vs COUPLED: holding (be) ~fixed in a bin, what chi-range is reachable?
      and holding chi ~fixed, do the be's still fill the octant?  (product <=> both free)
 [G4] hidden constraint: is there an inequality coupling chi to the be's? scan the
      (sum-window-action, chi) and (overlap be_lh, chi) envelopes.
 [G5] boundary-strata map: chi->1 vs window separation; be_ij->0 vs decoupling.
Reproduce: python3 ws_moduli_geometry.py
"""
import os,sys,numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
import num_S12 as NS
from assay import Params, Geometry
rng=np.random.default_rng(0)

def sample(n):
    out=[]
    while len(out)<n:
        eps=np.sort(rng.uniform(-4,4,3))
        if np.min(np.diff(eps))<0.15: continue
        gam=rng.uniform(0.2,2.0,3)*rng.choice([-1,1],3)
        a=np.sort(rng.uniform(-3,3,3))
        if np.min(np.diff(a))<0.1: continue
        lo,mid,hi=0,1,2
        try:
            geo=Geometry(Params(eps=tuple(eps),gam=tuple(gam),a=tuple(a)))
            chi=float(np.real(NS.q4_cross_ratio(geo)))
        except Exception: continue
        if not np.isfinite(chi): continue
        be=lambda i,j: gam[i]**2*gam[j]**2*abs(a[i]-a[j])/(eps[i]-eps[j])**2
        out.append((be(lo,mid),be(mid,hi),be(lo,hi),chi,eps[1]-eps[0],eps[2]-eps[1]))
    return np.array(out)

D=sample(20000)
be=D[:,:3]; chi=D[:,3]; lbe=np.log10(be)
print("="*82); print(f"MODULI IMAGE from {len(D)} samples"); print("="*82)
print(f"[G1] chi range: [{chi.min():.6f}, {chi.max():.6f}]  in (0,1): {np.all((chi>0)&(chi<1))}")
print(f"     chi quantiles 1/25/50/75/99%: {np.round(np.quantile(chi,[.01,.25,.5,.75,.99]),4)}")
print(f"[G2] log10(be) marginals (min..max):")
for k,nm in enumerate(['be_lm','be_mh','be_lh']):
    print(f"     {nm}: {lbe[:,k].min():+.2f} .. {lbe[:,k].max():+.2f}  "
          f"(quart {np.round(np.quantile(lbe[:,k],[.25,.5,.75]),2)})")
# octant coverage: correlation matrix of log-be (independence => ~0 off-diagonal)
C=np.corrcoef(lbe.T)
print(f"     corr(log be) off-diag: lm-mh {C[0,1]:+.2f}  lm-lh {C[0,2]:+.2f}  mh-lh {C[1,2]:+.2f}")
print(f"[G3] PRODUCT test:")
# chi-spread within a tight be-box (near medians)
med=np.median(lbe,0)
box=np.all(np.abs(lbe-med)<0.25,1)
print(f"     in a tight log-be box ({box.sum()} pts): chi spans [{chi[box].min():.3f},{chi[box].max():.3f}]"
      f"  (vs global [{chi.min():.3f},{chi.max():.3f}]) -> chi {'free' if chi[box].max()-chi[box].min()>0.5 else 'constrained'} at fixed be")
# be-spread within a tight chi-band
band=np.abs(chi-0.5)<0.05
print(f"     in a tight chi band ~0.5 ({band.sum()} pts): log-be_lm spans "
      f"[{lbe[band,0].min():.2f},{lbe[band,0].max():.2f}] -> be {'free' if lbe[band,0].max()-lbe[band,0].min()>2 else 'constrained'} at fixed chi")
print(f"[G4] coupling/envelope test: corr(chi, log be) = "
      f"{[f'{np.corrcoef(chi,lbe[:,k])[0,1]:+.2f}' for k in range(3)]}")
# is chi bounded by a function of be? check min/max chi vs overlap exponent be_lh
qs=np.quantile(lbe[:,2],np.linspace(0,1,9))
print(f"     chi-envelope vs log be_lh (overlap): [lo_bin: (chi_min,chi_max)]")
for i in range(8):
    m=(lbe[:,2]>=qs[i])&(lbe[:,2]<qs[i+1])
    if m.sum()>20: print(f"       be_lh~1e{(qs[i]+qs[i+1])/2:+.1f}: ({chi[m].min():.3f},{chi[m].max():.3f})")
print(f"[G5] boundary map: corr(chi, log eps-gap sum) = "
      f"{np.corrcoef(chi,np.log10(D[:,4]+D[:,5]))[0,1]:+.2f}  "
      f"(chi->1 as windows separate: {np.corrcoef(chi,np.log10(D[:,4]+D[:,5]))[0,1]>0})")
