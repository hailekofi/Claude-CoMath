"""
ws_type_M_levelstats.py -- APPLICATION probe: level statistics of the Type-M integrable (Cauchy/Gaudin)
family, and the integrability-breaking crossover.

Question (open, not about sigma): the Type-1/Type-M commuting family is INTEGRABLE. What are its level
statistics? Poisson (Berry-Tabor integrable), Wigner-Dyson/GOE (chaotic), picket-fence (rigid), or a NEW
class? And as integrability is broken (H0 + eps*GOE), does it cross over to GOE -- and is the crossover
controlled by the Type-1 geometry?

Diagnostic: the unfolding-FREE consecutive level-spacing RATIO r_n = min(s_n,s_{n+1})/max(s_n,s_{n+1})
(Atas-Bogomolny-Giraud-Roux). <r>: Poisson 0.3863, GOE 0.5307, GUE 0.5996, picket-fence -> 0 (rigid).

Type-M Cauchy/Gaudin H0:  (H0)_ij = gam_i gam_j (a_i-a_j)/(eps_i-eps_j),  (H0)_ii = -sum_k!=i gam_k^2 (a_i-a_k)/(eps_i-eps_k).
"""
from __future__ import annotations
import numpy as np
np.seterr(all="ignore")
POISSON, GOE, GUE = 0.38629, 0.53590, 0.60266   # reference <r>


def type_M(eps, gam, a):
    M = len(eps); H0 = np.zeros((M, M))
    g2 = gam ** 2
    for i in range(M):
        for j in range(M):
            if i != j:
                H0[i, j] = gam[i] * gam[j] * (a[i] - a[j]) / (eps[i] - eps[j])
        H0[i, i] = -np.sum([g2[k] * (a[i] - a[k]) / (eps[i] - eps[k]) for k in range(M) if k != i])
    return H0


def rmean(ev):
    ev = np.sort(ev.real); s = np.diff(ev)
    s = s[s > 1e-12]
    r = np.minimum(s[:-1], s[1:]) / np.maximum(s[:-1], s[1:])
    return r


def draw(M, rng):
    eps = np.sort(rng.uniform(-1, 1, M)) * M / 2          # sites spread ~ uniform density
    # enforce distinctness
    eps += 1e-6 * np.arange(M)
    a = rng.uniform(-1, 1, M)
    gam = rng.normal(0, 1, M)
    return eps, gam, a


def main():
    print("=" * 88)
    print("Level statistics of the Type-M integrable (Cauchy/Gaudin) family")
    print(f"  reference <r>:  Poisson={POISSON}  GOE={GOE}  GUE={GUE}  (picket-fence -> 0)")
    print("=" * 88)

    print("\n(1) INTEGRABLE Type-M H0 -- <r> vs M (ensemble-averaged):")
    for M in [100, 200, 400]:
        rng = np.random.default_rng(M); allr = []
        for _ in range(30):
            eps, gam, a = draw(M, rng)
            allr.append(rmean(np.linalg.eigvalsh(type_M(eps, gam, a))))
        r = np.concatenate(allr)
        print(f"  M={M:4d}:  <r> = {np.mean(r):.4f}   (N_gaps={len(r)})")

    print("\n(2) INTEGRABILITY-BREAKING crossover  H = H0 + eps_break * (GOE matrix), M=300:")
    M = 300; rng = np.random.default_rng(7)
    print(f"  {'eps_break':>10} | {'<r>':>7}  (Poisson 0.386 -> GOE 0.531)")
    for eb in [0.0, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0]:
        allr = []
        for _ in range(20):
            eps, gam, a = draw(M, rng)
            H0 = type_M(eps, gam, a)
            G = rng.normal(0, 1, (M, M)); G = (G + G.T) / 2 / np.sqrt(M)   # GOE, unit-ish bandwidth
            # scale eps_break by the H0 bandwidth so the knob is dimensionless
            bw = np.std(np.linalg.eigvalsh(H0))
            allr.append(rmean(np.linalg.eigvalsh(H0 + eb * bw * G)))
        r = np.concatenate(allr)
        print(f"  {eb:10.2f} | {np.mean(r):.4f}")

    print("\n(3) Control -- pure GOE and pure diagonal (Poisson), same machinery:")
    rng = np.random.default_rng(1)
    rg = []; rp = []
    for _ in range(20):
        G = rng.normal(0, 1, (300, 300)); G = (G + G.T) / 2
        rg.append(rmean(np.linalg.eigvalsh(G)))
        rp.append(rmean(np.sort(rng.uniform(-1, 1, 300))))
    print(f"  pure GOE:  <r>={np.mean(np.concatenate(rg)):.4f} (expect 0.531)")
    print(f"  pure Poisson (random levels): <r>={np.mean(np.concatenate(rp)):.4f} (expect 0.386)")


if __name__ == "__main__":
    main()
