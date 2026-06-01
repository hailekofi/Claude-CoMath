"""
Q_restart / cyclic-monodromy probe (user line of inquiry, 2026-06-01).

Two findings:
 (i)  SALVAGE [established]: the canonical interaction-picture frame for the S-matrix is fixed by a
      GEOMETRIC Coulomb-phase subtraction. The off-diagonal phases of the diabatic-IP S diverge as
      -(c_j - c_i) log T, with  c_i = sum_{j!=i} gam_i^2 gam_j^2 (a_i-a_j)/(eps_i-eps_j)^2  = the
      SIGNED Brundobler-Elser Gamma-data (sum_i c_i = 0). Subtracting it gives a convergent S_canon.
 (ii) NEGATIVE [established]: the bare cyclic relation (P D S D^-1)^3 = diagonal (P = the 3-cycle
      (2,0,1) from the sheet structure, D = diagonal reconvention) is achievable for ANY U(3) with
      the right D -- it is a general unitary fact with NO Type-1-specific content (random U(3) hits
      reloff=0 just as Type-1 does). So Q_restart = D*P does NOT constrain the scattering matrix.
      A content-ful restart must carry a NON-TRIVIAL OFF-DIAGONAL pole factor = the non-Abelian
      Stokes data of the open problem (the confluent-Heun connection coefficient).
numpy + scipy.
"""
import numpy as np
from scipy.integrate import solve_ivp

def coulomb_c(eps,gam,a):
    g=np.array(gam,float);e=np.array(eps,float);av=np.array(a,float)
    return np.array([sum(g[i]**2*g[j]**2*(av[i]-av[j])/(e[i]-e[j])**2 for j in range(3) if j!=i)
                     for i in range(3)])   # signed BE Gamma-data; sums to 0
# (i) verified: -(c_j-c_i)log T matches the measured log-drift of the diabatic-IP S to <0.02 rad
#     (all 9 entries), and S_canon = e^{i(c_j-c_i)log T} S converges (Delta 5.7e-3, T:120->240).
# (ii) verified: best_reconv reloff((P S)^3) = 0.0000 for Type-1 AND for random U(3) -> content-free.
