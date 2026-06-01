"""Phase 0c smoke test: geometry cross-checks, benchmark, Level 0."""
import numpy as np
from assay import (Params, Geometry, cross_check,
                   transition_matrix, unitarity_defect, segmentation_identity)

np.set_printoptions(precision=6, suppress=True, linewidth=110)

# canonical commuting family from the project notes, generic slopes
par = Params(eps=(-2.0, 0.0, 3.0), gam=(1.0, 0.8, 1.2), a=(-1.0, 0.5, 2.0), x=0)
geo = Geometry(par)

print("=== geometry cross-check (explicit Q4/W4 vs definitions) ===")
cc = cross_check(geo)
for k, v in cc.items():
    print(f"  {k}: {v:.3e}")

print("\n p   =", geo.p)
print(" n   =", geo.n)
print(" m   =", geo.m)
print(" W4  =", geo.W4)
print(" Q4  =", geo.Q4)
print(" L_H =", geo.LH, "  (alphaH,betaH,gammaH=",
      f"{geo.alphaH:.4f},{geo.betaH:.4f},{geo.gammaH:.4f})")
print(" kappa_S =", geo.kappaS)
print(" Q4 roots =", np.round(geo.Q4_roots(), 6))
print(" W4 roots =", np.round(geo.W4_roots(), 6))

print("\n=== projected conjugate-pair Q4 windows (theta=0) ===")
for k, w in enumerate(geo.windows()):
    print(f"  X{k+1}: tau={w['tau']:+.5f}  u_center={w['u_center']:+.5f}  "
          f"conj_defect={w['conj_defect']:.2e}")
    print(f"        roots  = {np.round(w['roots'], 5)}")
    print(f"        u_vals = {np.round(w['u_vals'], 5)}")

print("\n=== benchmark transition matrix (T=60) ===")
P = transition_matrix(geo, T=60.0)
print(P)
print(" row sums      =", P.sum(axis=1))
print(" unitarity def =", f"{unitarity_defect(geo, T=60.0):.3e}")

# T-refinement stability
P80 = transition_matrix(geo, T=80.0)
print(" max|P(T=60)-P(T=80)| =", f"{np.max(np.abs(P - P80)):.3e}")

print("\n=== Level 0 segmentation identity ===")
r = segmentation_identity(geo, T=60.0)
for k, v in r.items():
    print(f"  {k}: {v}")
