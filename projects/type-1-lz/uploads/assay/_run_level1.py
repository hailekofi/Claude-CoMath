"""Level 1 confirmation batch (interaction-picture engine)."""
import time
import numpy as np
from assay import Params, Geometry
from assay.level1 import level1_assay, adiabatic_scattering
from assay.windows import segmentation

np.set_printoptions(precision=6, suppress=True, linewidth=110)

SAMPLES = [
    dict(eps=(-2., 0., 3.), gam=(1.0, 0.8, 1.2), a=(-1.0, 0.5, 2.0)),
    dict(eps=(-1.5, 0.4, 2.2), gam=(0.9, 1.1, 0.7), a=(0.3, -0.8, 1.4)),
    dict(eps=(-3.0, -0.5, 1.8), gam=(1.3, 0.6, 1.0), a=(-0.5, 1.2, -1.7)),
]

t0 = time.time()
results = []
for s, sp in enumerate(SAMPLES):
    for x in (0, 1, 2):
        geo = Geometry(Params(x=x, **sp))
        xc = (x == 0)
        r = level1_assay(geo, T=50.0, x=x, cross_check=xc, verbose=False)
        results.append((s, x, r))
        msg = (f"sample {s} x={x}: seg={r['seg_consistency']:.2e}  "
               f"bench={r['benchmark_agreement']:.2e}  PASS={r['passed']}")
        if xc:
            msg += (f"  | diabatic xcheck={r['diabatic_xcheck_err']:.2e} "
                    f"perm={r['diabatic_xcheck_perm']}")
        print(msg)

print("\nT-refinement (sample 0, x=0):")
geo = Geometry(Params(x=0, **SAMPLES[0]))
for T in (40.0, 60.0, 90.0):
    r = level1_assay(geo, T=T, x=0, verbose=False)
    print(f"  T={T}: seg={r['seg_consistency']:.2e}  bench={r['benchmark_agreement']:.2e}")

print("\nwindow half-width stability (sample 0, x=0, T=50):")
geo = Geometry(Params(x=0, **SAMPLES[0]))
base = None
for hw_frac in (0.15, 0.30, 0.45):
    seg0 = segmentation(geo, T=50.0)
    gap = seg0["centers"][1] - seg0["centers"][0]
    seg = segmentation(geo, T=50.0, half_width=hw_frac * gap)
    _, M_seg, _ = adiabatic_scattering(geo, 50.0, 0, seg=seg)
    P = np.abs(M_seg.T) ** 2
    if base is None:
        base = P
    print(f"  half_width={hw_frac:.2f}*gap: max|P-P(0.15)| = {np.max(np.abs(P-base)):.2e}")

allpass = all(r["passed"] for _, _, r in results)
print(f"\nALL PASS: {allpass}   elapsed: {time.time()-t0:.1f}s")
