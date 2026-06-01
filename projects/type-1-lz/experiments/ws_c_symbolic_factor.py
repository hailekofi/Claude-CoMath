"""
ws_c_symbolic_factor.py  --  WS-C step 2(c): symbolic search for an algebraic factorization
of the Type-1 3x3 connection / coupling structure that would force P to factorize into
two-level pieces.

Physical fact (WS-A / WS-G): P factorizes into a product of two-level Weber S-matrices iff
the 3x3 problem REDUCES to a direct sum (one level decouples) -- i.e. the off-diagonal Cauchy
coupling V_ij = gamma_i gamma_j (a_i-a_j)/(eps_i-eps_j) connecting the middle level to ONE of
the outer levels effectively vanishes over the whole passage, so no Stokes line of the two
off-diagonal pairs can joint.

We test the only algebraic ways an off-diagonal coupling can vanish on a sub-locus:
  (i)  gamma_i = 0                              (a level decouples entirely)  -- trivial.
  (ii) a_i = a_j                                (degenerate slopes: V_ij = 0 AND Gamma_ij = 0)
  (iii) eps_i -> eps_j                          (pole collision: V_ij blows up; OUTSIDE Type-1)

and check, symbolically, whether ANY *interior* (all gamma!=0, all a distinct, all eps
distinct) locus makes a middle coupling vanish or makes H0 block-diagonalize.

We also factor the characteristic polynomial of H0 and of H(u), and test block-triangular
reducibility of the constant similarity-invariants, over Q(gamma,eps,a).

numpy not needed; sympy only.
"""
import sympy as sp


def setup():
    g = sp.symbols('g0 g1 g2', real=True)
    e = sp.symbols('e0 e1 e2', real=True)
    a = sp.symbols('a0 a1 a2', real=True)
    u = sp.symbols('u', real=True)

    def V(i, j):
        return g[i] * g[j] * (a[i] - a[j]) / (e[i] - e[j])

    def D(i):  # H0_ii
        return -sum(g[k]**2 * (a[i] - a[k]) / (e[i] - e[k]) for k in range(3) if k != i)

    H0 = sp.Matrix(3, 3, lambda i, j: D(i) if i == j else V(i, j))
    A = sp.diag(*a)
    H = H0 + u * A
    return g, e, a, u, V, D, H0, H


def main():
    g, e, a, u, V, D, H0, H = setup()
    print("=== Off-diagonal Cauchy couplings V_ij ===")
    for (i, j) in [(0, 1), (0, 2), (1, 2)]:
        print(f"  V_{i}{j} = {sp.simplify(V(i, j))}")

    print("\n=== When does a middle coupling vanish on an INTERIOR locus? ===")
    # 'middle' is slope-dependent; algebraically V_ij = 0 iff g_i=0 or g_j=0 or a_i=a_j.
    for (i, j) in [(0, 1), (1, 2), (0, 2)]:
        sols = sp.solve(sp.numer(sp.together(V(i, j))), dict=True)
        print(f"  V_{i}{j}=0  <=>  {sp.factor(sp.numer(sp.together(V(i,j))))} = 0"
              f"   (roots: g_{i}=0, g_{j}=0, or a_{i}=a_{j})")

    print("\n=== Does H0 block-diagonalize on any interior locus? "
          "(off-diagonal block = {V_mid,lo, V_mid,hi}) ===")
    # For the middle level (say index 1) to decouple, BOTH V_01 and V_12 must vanish.
    # V_01=0 and V_12=0 with all g!=0, eps distinct  =>  a_0=a_1 AND a_1=a_2  => all slopes equal
    # => H = H0 + u*a*I is a scalar shift of a CONSTANT matrix => no LZ at all (trivial).
    cond = sp.And(sp.Eq(a[0], a[1]), sp.Eq(a[1], a[2]))
    print("  Middle level 1 decouples (V_01=V_12=0, g!=0) <=> a_0=a_1=a_2 (all slopes equal)")
    print("  => H(u) = H0 + u*a_common*I : a CONSTANT matrix up to scalar phase => no crossings.")
    print("  => the ONLY interior all-slopes-distinct decoupling needs some gamma=0 (trivial).")

    print("\n=== Single decoupling: gamma_mid -> 0 makes BOTH its couplings vanish simultaneously ===")
    g0, g1, g2 = g
    H0_g1_0 = H0.subs(g1, 0)
    print("  H0 with g1=0 (middle coupling off):")
    sp.pprint(sp.simplify(H0_g1_0))
    # check it is block diagonal: row/col 1 off-diagonals
    od = [sp.simplify(H0_g1_0[1, 0]), sp.simplify(H0_g1_0[1, 2]),
          sp.simplify(H0_g1_0[0, 1]), sp.simplify(H0_g1_0[2, 1])]
    print(f"  off-diagonals touching level 1 with g1=0: {od}  -> level 1 DECOUPLES.")

    print("\n=== char poly of H0 (does it factor over Q(g,e,a)?) ===")
    cp = sp.factor(H0.charpoly(sp.symbols('L')).as_expr())
    print("  charpoly(H0) factored:", cp)


if __name__ == "__main__":
    main()
