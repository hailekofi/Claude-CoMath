"""Verify IP propagators against lab-frame, and time them."""
import time
import numpy as np
from assay import Params, Geometry
from assay.benchmark import fundamental_matrix
from assay.level1 import propagate as propagate_lab
from assay.ip import benchmark_ip_U, propagate_ad_ip

np.set_printoptions(precision=4, suppress=True, linewidth=110)
geo = Geometry(Params(eps=(-2., 0., 3.), gam=(1., 0.8, 1.2), a=(-1., 0.5, 2.), x=0))

print("=== correctness: IP vs lab-frame on short intervals ===")
# diabatic benchmark, short span
for T in (3.0, 8.0):
    Ulab = fundamental_matrix(geo, -T, T).U
    Uip = benchmark_ip_U(geo, T=T)
    print(f"  benchmark  T={T}:  max|U_ip - U_lab| = {np.max(np.abs(Uip-Ulab)):.2e}")

# adiabatic constructor, short segments
for (a, b) in [(-8., -2.), (-3., 4.), (2., 9.)]:
    Ulab = propagate_lab(geo, a, b, 0)
    Uip = propagate_ad_ip(geo, a, b, 0)
    print(f"  adiabatic  [{a},{b}]:  max|U_ip - U_lab| = {np.max(np.abs(Uip-Ulab)):.2e}")

print("\n=== speed: tail segment that took 30.5 s lab-frame ===")
t = time.time()
Uip = propagate_ad_ip(geo, -50., -40., 0)
print(f"  adiabatic IP [-50,-40]: {time.time()-t:.3f} s")

t = time.time()
Uip = benchmark_ip_U(geo, T=60.0)
print(f"  benchmark IP [-60,60] full pass: {time.time()-t:.3f} s")
